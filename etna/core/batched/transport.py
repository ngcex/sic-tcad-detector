"""Batched steady-state electron/hole continuity (Scharfetter-Gummel, SRH, optical
generation), float64 CPU, Gummel-coupled to the Poisson solver in ``poisson.py``.

Carrier densities in SiC span 1e-35 .. 1e19 cm^-3. Discretisation and contact
values follow devsim's simple_dd / simple_physics so results can be compared
point by point.

Edge i (nodes i -> i+1), Delta_i = (psi_i - psi_{i+1})/V_t, B(x) = x/(e^x - 1):

    Jn_i = q mu_n V_t/h_i [ n_{i+1} B(-Delta_i) - n_i B(Delta_i) ]
    Jp_i = q mu_p V_t/h_i [ p_i B(-Delta_i) - p_{i+1} B(Delta_i) ]
    Jn_i - Jn_{i-1} =  q v_i (U_i - G_i),   Jp_i - Jp_{i-1} = -q v_i (U_i - G_i)
    U = (n p - n_i^2) / (tau_p (n + n_i) + tau_n (p + n_i))        [SRH, midgap]

Each Gummel iteration solves n (p frozen in U), then p (n frozen), then the
Poisson equation with nodal quasi-Fermi levels implied by (n, p).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit, prange

from etna.core.batched.material import Q


@njit(cache=True, parallel=True)
def thomas_batch(a, b, c, d):
    """Tridiagonal solve of B independent systems stored as (N, B) float64 arrays."""
    N, B = b.shape
    x = np.empty_like(d)
    for j in prange(B):
        cp = np.empty(N)
        dp = np.empty(N)
        inv = 1.0 / b[0, j]
        cp[0] = c[0, j] * inv
        dp[0] = d[0, j] * inv
        for i in range(1, N):
            inv = 1.0 / (b[i, j] - a[i, j] * cp[i - 1])
            cp[i] = c[i, j] * inv
            dp[i] = (d[i, j] - a[i, j] * dp[i - 1]) * inv
        x[N - 1, j] = dp[N - 1]
        for i in range(N - 2, -1, -1):
            x[i, j] = dp[i] - cp[i] * x[i + 1, j]
    return x


def bernoulli(x):
    """B(x) = x/(exp(x)-1), stable for all float64 x."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-8
    big = x > 700.0
    neg = x < -700.0
    mid = ~(small | big | neg)
    out[small] = 1.0 - 0.5 * x[small]
    out[mid] = x[mid] / np.expm1(x[mid])
    out[big] = 0.0
    out[neg] = -x[neg]
    return out


@dataclass
class TransportResult:
    psi: np.ndarray  # (B, N) V
    n: np.ndarray  # (B, N)
    p: np.ndarray  # (B, N)
    J_cathode: np.ndarray  # (B,) A/cm^2, cathode Jn + Jp (Jn in the balanced form)
    gummel_iterations: int
    converged: np.ndarray  # (B,) bool


class ContinuityBatch:
    """Fixed per-batch transport data.

    ``mu_n`` / ``mu_p``: EDGE mobilities (N-1, B) (``material.edge_mobilities``, the
    devsim builders' local model), or a scalar / (B,) uniform mobility (kept for unit
    tests and callers that want one value per device). ``tau`` scalar or (B,).
    """

    def __init__(self, pb, mu_n, mu_p, tau_n, tau_p):
        self.pb = pb  # PoissonBatch -- supplies mesh, doping, n_i, V_t
        B, N = pb.B, pb.N
        bc = lambda a: np.broadcast_to(np.asarray(a, dtype=float), (B,)).copy()

        def edges(a):
            a = np.asarray(a, dtype=float)
            if a.ndim == 2:
                if a.shape != (N - 1, B):
                    raise ValueError(f"edge mobility shape {a.shape} != {(N - 1, B)}")
                return a.copy()
            return np.broadcast_to(bc(a)[None, :], (N - 1, B)).copy()

        self.mu_n, self.mu_p = edges(mu_n), edges(mu_p)  # (N-1, B)
        self.tau_n, self.tau_p = bc(tau_n), bc(tau_p)
        self.h = np.diff(pb.x)[:, None]  # (N-1, 1)
        self.v = pb.v[:, None]  # (N, 1), ends unused
        self.kn = (Q * self.mu_n * pb.V_t[None, :]) / self.h  # (N-1, B)
        self.kp = (Q * self.mu_p * pb.V_t[None, :]) / self.h
        D = pb.netdop
        n_i = pb.n_i
        celec = lambda d: 1e-10 + 0.5 * np.abs(d + np.sqrt(d**2 + 4 * n_i**2))
        chole = lambda d: 1e-10 + 0.5 * np.abs(-d + np.sqrt(d**2 + 4 * n_i**2))
        self.n_bc = np.stack(
            [
                np.where(D[:, 0] > 0, celec(D[:, 0]), n_i**2 / chole(D[:, 0])),
                np.where(D[:, -1] > 0, celec(D[:, -1]), n_i**2 / chole(D[:, -1])),
            ]
        )  # (2, B)
        self.p_bc = np.stack(
            [
                np.where(D[:, 0] < 0, chole(D[:, 0]), n_i**2 / celec(D[:, 0])),
                np.where(D[:, -1] < 0, chole(D[:, -1]), n_i**2 / celec(D[:, -1])),
            ]
        )

    def _solve_carrier(self, psi_T, n_other, which, G_T, denom, bc):
        """One linear SG continuity solve for electrons ('n') or holes ('p'); (N, B) arrays."""
        pb = self.pb
        N, B = psi_T.shape
        delta = (psi_T[:-1] - psi_T[1:]) / pb.V_t[None, :]  # (N-1, B)
        Bp, Bm = bernoulli(delta), bernoulli(-delta)
        k = self.kn if which == "n" else self.kp
        a = np.zeros((N, B))
        b = np.ones((N, B))
        c = np.zeros((N, B))
        rhs = np.zeros((N, B))
        n_i2 = (pb.n_i**2)[None, :]
        qv = Q * self.v[1:-1]
        # SRH linearised in the solved carrier: U = (x * other - n_i^2)/denom
        srh_lin = qv * n_other[1:-1] / denom[1:-1]
        srh_const = qv * n_i2 / denom[1:-1]
        if which == "n":
            # Jn_i - Jn_{i-1} - q v (U - G) = 0
            c[1:-1] = k[1:] * Bm[1:]
            a[1:-1] = k[:-1] * Bp[:-1]
            b[1:-1] = -(k[1:] * Bp[1:] + k[:-1] * Bm[:-1]) - srh_lin
            rhs[1:-1] = -qv * G_T[1:-1] - srh_const
        else:
            # Jp_i - Jp_{i-1} + q v (U - G) = 0
            c[1:-1] = -k[1:] * Bp[1:]
            a[1:-1] = -k[:-1] * Bm[:-1]
            b[1:-1] = (k[1:] * Bm[1:] + k[:-1] * Bp[:-1]) + srh_lin
            rhs[1:-1] = qv * G_T[1:-1] + srh_const
        rhs[0], rhs[-1] = bc[0], bc[1]
        return thomas_batch(a, b, c, rhs)

    def currents(self, psi_T, n_T, p_T):
        """Jn, Jp on every edge, (N-1, B), A/cm^2."""
        delta = (psi_T[:-1] - psi_T[1:]) / self.pb.V_t[None, :]
        Bp, Bm = bernoulli(delta), bernoulli(-delta)
        Jn = self.kn * (n_T[1:] * Bm - n_T[:-1] * Bp)
        Jp = self.kp * (p_T[:-1] * Bm - p_T[1:] * Bp)
        return Jn, Jp

    def solve(
        self,
        psi,
        V_anode=0.0,
        V_cathode=0.0,
        G=None,
        n0=None,
        p0=None,
        tol=1e-6,
        max_iter=50,
        poisson_tol=1e-8,
    ):
        """Gummel iteration from a converged Poisson state psi (B, N).

        G: optical generation (B, N) or (N,) cm^-3 s^-1 (None = dark).
        tol: relative change of the cathode current between successive Gummel
        iterations (and of psi [V_t], n, p, floored at 2*poisson_tol).
        """
        pb = self.pb
        B, N = pb.B, pb.N
        V_a = np.broadcast_to(np.asarray(V_anode, dtype=float), (B,))
        V_c = np.broadcast_to(np.asarray(V_cathode, dtype=float), (B,))
        G_T = (
            np.zeros((N, B))
            if G is None
            else np.broadcast_to(np.asarray(G, dtype=float), (B, N)).T.copy()
        )
        psi_T = np.asarray(psi, dtype=float).T.copy()
        n_i = pb.n_i[None, :]
        # start from flat-quasi-Fermi carriers (exact in the dark for reverse bias)
        u = psi_T / pb.V_t[None, :]
        n_T = (
            n_i * np.exp(u - (V_c / pb.V_t)[None, :])
            if n0 is None
            else np.asarray(n0, float).T.copy()
        )
        p_T = (
            n_i * np.exp((V_a / pb.V_t)[None, :] - u)
            if p0 is None
            else np.asarray(p0, float).T.copy()
        )
        tau_n, tau_p = self.tau_n[None, :], self.tau_p[None, :]
        converged = np.zeros(B, dtype=bool)
        J_prev = None
        it = 0
        for it in range(1, max_iter + 1):
            denom = tau_p * (n_T + n_i) + tau_n * (p_T + n_i)
            n_new = self._solve_carrier(psi_T, p_T, "n", G_T, denom, self.n_bc)
            n_new = np.maximum(n_new, 0.0)
            denom = tau_p * (n_new + n_i) + tau_n * (p_T + n_i)
            p_new = self._solve_carrier(psi_T, n_new, "p", G_T, denom, self.p_bc)
            p_new = np.maximum(p_new, 0.0)
            # cathode current of the continuity-consistent state (psi_T, n_new, p_new): it is a
            # tiny difference of two large SG terms, so it is NOT re-evaluated after the Poisson
            # update below
            Jn, Jp = self.currents(psi_T, n_new, p_new)
            # Balanced form (as devsim's extract_steady_state_current): summing the electron
            # equation over the interior nodes gives Jn(cathode edge) = Jn(anode edge) +
            # sum q v (U - G); both terms are well conditioned (electrons are minority at the
            # p+ anode), unlike Jn[-1] itself, which is a difference of huge SG terms when n
            # is large at the cathode (n+ layer). Holes are minority at the cathode: direct.
            denom = tau_p * (n_new + n_i) + tau_n * (p_new + n_i)
            U = (n_new * p_new - n_i**2) / denom
            S = (Q * self.v[1:-1] * (U[1:-1] - G_T[1:-1])).sum(axis=0)
            J = Jn[0] + S + Jp[-1]
            dJ = (
                np.abs(J - J_prev) / (np.abs(J) + 1e-300)
                if J_prev is not None
                else np.full(B, np.inf)
            )
            J_prev = J
            # Poisson with nodal quasi-Fermi levels implied by (n, p)
            phi_n = psi_T - pb.V_t[None, :] * np.log(np.maximum(n_new, 1e-300) / n_i)
            phi_p = psi_T + pb.V_t[None, :] * np.log(np.maximum(p_new, 1e-300) / n_i)
            res = pb.solve(
                psi_T.T,
                V_anode=V_a,
                V_cathode=V_c,
                tol=poisson_tol,
                phi_n=phi_n.T,
                phi_p=phi_p.T,
            )
            psi_new = res.psi.T
            dpsi = np.abs(psi_new - psi_T).max(axis=0) / pb.V_t
            # carriers move with the potential update (Slotboom-consistent)
            shift = (psi_new - psi_T) / pb.V_t[None, :]
            n_T2 = n_new * np.exp(shift)
            p_T2 = p_new * np.exp(-shift)
            # relative change, floored at n_i so swept-out minority tails do not count
            dn = (np.abs(n_T2 - n_T) / (np.abs(n_T) + n_i)).max(axis=0)
            dp = (np.abs(p_T2 - p_T) / (np.abs(p_T) + n_i)).max(axis=0)
            psi_T, n_T, p_T = psi_new, n_T2, p_T2
            # psi, n, p are floored at 2*poisson_tol (the Poisson solver's own tolerance);
            # the cathode current -- the quantity CCE needs -- must still settle to tol.
            floor = max(tol, 2 * poisson_tol)
            converged = (dpsi < floor) & (dn < floor) & (dp < floor) & (dJ < tol)
            if converged.all():
                break
        return TransportResult(
            psi=psi_T.T,
            n=n_T.T,
            p=p_T.T,
            J_cathode=J,
            gummel_iterations=it,
            converged=converged,
        )
