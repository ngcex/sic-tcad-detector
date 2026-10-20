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
from etna.core.batched.poisson import LINEAR_PAR_NODES, use_prange

_NO_AUG = np.zeros((1, 1))  # placeholder ``aug`` when Auger is off


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


@njit(cache=True, inline="always")
def _bern(x):
    """Scalar ``bernoulli`` (same branches and expressions)."""
    if abs(x) < 1e-8:
        return 1.0 - 0.5 * x
    if x > 700.0:
        return 0.0
    if x < -700.0:
        return -x
    return x / np.expm1(x)


@njit(cache=True)
def _end_currents(psi_T, n_T, p_T, kn, kp, V_t):
    """Jn on the anode edge and Jp on the cathode edge, (B,) each: the two entries of
    ``ContinuityBatch.currents`` the Gummel loop uses, same expressions."""
    N, B = psi_T.shape
    Jn0 = np.empty(B)
    JpL = np.empty(B)
    for j in range(B):
        d = (psi_T[0, j] - psi_T[1, j]) / V_t[j]
        Jn0[j] = kn[0, j] * (n_T[1, j] * _bern(-d) - n_T[0, j] * _bern(d))
        d = (psi_T[N - 2, j] - psi_T[N - 1, j]) / V_t[j]
        JpL[j] = kp[N - 2, j] * (p_T[N - 2, j] * _bern(-d) - p_T[N - 1, j] * _bern(d))
    return Jn0, JpL


@njit(cache=True)
def _carrier_column(
    j, x, psi_T, n_other, G_T, denom, aug, has_aug, k, v, V_t, n_i, bc0, bc1, electrons
):
    """Assembly + Thomas solve of the linear SG continuity system of
    ``ContinuityBatch._solve_carrier`` for device j in one pass, with the numpy
    assembly's expressions and operation order; writes x[:, j]."""
    N = psi_T.shape[0]
    cp = np.empty(N)
    dp = np.empty(N)
    vt = V_t[j]
    ni2 = n_i[j] ** 2
    # row 0: identity (a = c = 0, b = 1), rhs = bc0
    cp[0] = 0.0
    dp[0] = bc0[j]
    dl = (psi_T[0, j] - psi_T[1, j]) / vt  # edge i-1 of row i = 1
    Bp_l = _bern(dl)
    Bm_l = _bern(-dl)
    for i in range(1, N - 1):
        dr = (psi_T[i, j] - psi_T[i + 1, j]) / vt  # edge i
        Bp_r = _bern(dr)
        Bm_r = _bern(-dr)
        kl = k[i - 1, j]
        kr = k[i, j]
        qv = Q * v[i]
        srh_lin = qv * n_other[i, j] / denom[i, j]
        srh_const = qv * ni2 / denom[i, j]
        if has_aug:
            srh_lin = srh_lin + qv * aug[i, j] * n_other[i, j]
            srh_const = srh_const + qv * aug[i, j] * ni2
        if electrons:
            ci = kr * Bm_r
            ai = kl * Bp_l
            bi = -(kr * Bp_r + kl * Bm_l) - srh_lin
            ri = -qv * G_T[i, j] - srh_const
        else:
            ci = -kr * Bp_r
            ai = -kl * Bm_l
            bi = (kr * Bm_r + kl * Bp_l) + srh_lin
            ri = qv * G_T[i, j] + srh_const
        inv = 1.0 / (bi - ai * cp[i - 1])
        cp[i] = ci * inv
        dp[i] = (ri - ai * dp[i - 1]) * inv
        Bp_l = Bp_r
        Bm_l = Bm_r
    # row N-1: identity, rhs = bc1 (a = 0: dp = bc1 / (1 - 0 * cp) as thomas_batch)
    inv = 1.0 / (1.0 - 0.0 * cp[N - 2])
    dp[N - 1] = (bc1[j] - 0.0 * dp[N - 2]) * inv
    x[N - 1, j] = dp[N - 1]
    for i in range(N - 2, -1, -1):
        x[i, j] = dp[i] - cp[i] * x[i + 1, j]


@njit(cache=True)
def _carrier_solve_seq(
    psi_T, n_other, G_T, denom, aug, has_aug, k, v, V_t, n_i, bc0, bc1, electrons
):
    """``_carrier_column`` over the B devices, serially (small B: a parallel region's
    launch costs more than the solves)."""
    N, B = psi_T.shape
    x = np.empty((N, B))
    for j in range(B):
        _carrier_column(
            j, x, psi_T, n_other, G_T, denom, aug, has_aug, k, v, V_t, n_i, bc0, bc1,
            electrons,
        )  # fmt: skip
    return x


@njit(cache=True, parallel=True)
def _carrier_solve_par(
    psi_T, n_other, G_T, denom, aug, has_aug, k, v, V_t, n_i, bc0, bc1, electrons
):
    """``_carrier_column`` over the B devices, one prange iteration each."""
    N, B = psi_T.shape
    x = np.empty((N, B))
    for j in prange(B):
        _carrier_column(
            j, x, psi_T, n_other, G_T, denom, aug, has_aug, k, v, V_t, n_i, bc0, bc1,
            electrons,
        )  # fmt: skip
    return x


def _carrier_solve_batch(psi_T, *args):
    """Serial or prange over the devices (``poisson.use_prange``), identical arithmetic."""
    N, B = psi_T.shape
    par = use_prange(B, N, LINEAR_PAR_NODES)
    return (_carrier_solve_par if par else _carrier_solve_seq)(psi_T, *args)


@dataclass
class TransportResult:
    psi: np.ndarray  # (B, N) V
    n: np.ndarray  # (B, N)
    p: np.ndarray  # (B, N)
    J_cathode: np.ndarray  # (B,) A/cm^2, cathode Jn + Jp (Jn in the balanced form)
    gummel_iterations: int
    converged: np.ndarray  # (B,) bool
    G: "np.ndarray | None" = None  # (B, N) generation of the last iteration (G + G_fn)


class ContinuityBatch:
    """Fixed per-batch transport data.

    ``mu_n`` / ``mu_p``: EDGE mobilities (N-1, B) (``material.edge_mobilities``, the
    devsim builders' local model), or a scalar / (B,) uniform mobility (kept for unit
    tests and callers that want one value per device). ``tau`` scalar or (B,).
    """

    def __init__(
        self, pb, mu_n, mu_p, tau_n, tau_p, n1=None, p1=None, C_n=0.0, C_p=0.0
    ):
        """``n1`` / ``p1`` (scalar or (B,), default n_i: midgap SRH) are the SRH trap
        densities of U = (n p - n_i^2) / (tau_p (n + n1) + tau_n (p + p1)): the dark-current
        model passes the Z1/2 level (etna.core.dark_current.tat_constants). ``C_n`` / ``C_p``
        (cm^6/s, default 0 = off) add bulk Auger U_A = (C_n n + C_p p)(n p - n_i^2), as
        etna.core.flash_recombination.add_auger_recombination. With the defaults every
        expression is the pre-2026-10-18 one (bit-identical CCE)."""
        self.pb = pb  # PoissonBatch -- supplies mesh, doping, n_i, V_t
        B, N = pb.B, pb.N
        bc = lambda a: np.broadcast_to(np.asarray(a, dtype=float), (B,)).copy()
        self.n1 = pb.n_i.copy() if n1 is None else bc(n1)
        self.p1 = pb.n_i.copy() if p1 is None else bc(p1)
        self.C_n, self.C_p = bc(C_n), bc(C_p)
        self.auger = bool(np.any(self.C_n != 0.0) or np.any(self.C_p != 0.0))

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
        self._v = np.ascontiguousarray(pb.v, dtype=float)  # (N,)
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

    def _solve_carrier(self, psi_T, n_other, which, G_T, denom, bc, aug=None):
        """One linear SG continuity solve for electrons ('n') or holes ('p'); (N, B) arrays.

        Interior row i (edges i-1, i; k = q mu V_t / h, Delta = (psi_i - psi_{i+1})/V_t):
        SRH linearised in the solved carrier, U = (x * other - n_i^2)/denom; with Auger
        (C_n n + C_p p)(n p - n_i^2) its prefactor ``aug`` is lagged (evaluated at the
        previous iterate): linear in the solved carrier, same fixed point as devsim's
        fully coupled Newton. Electrons: Jn_i - Jn_{i-1} - q v (U - G) = 0; holes:
        Jp_i - Jp_{i-1} + q v (U - G) = 0. Assembly and Thomas sweep are one compiled
        pass (``_carrier_solve_batch``); ``_solve_carrier_numpy`` is the reference.
        """
        pb = self.pb
        k = self.kn if which == "n" else self.kp
        has_aug = aug is not None
        return _carrier_solve_batch(
            np.ascontiguousarray(psi_T, dtype=float),
            np.ascontiguousarray(n_other, dtype=float),
            np.ascontiguousarray(G_T, dtype=float),
            np.ascontiguousarray(denom, dtype=float),
            np.ascontiguousarray(aug, dtype=float) if has_aug else _NO_AUG,
            has_aug,
            k,
            self._v,
            pb.V_t,
            pb.n_i,
            np.ascontiguousarray(bc[0], dtype=float),
            np.ascontiguousarray(bc[1], dtype=float),
            which == "n",
        )

    def _solve_carrier_numpy(self, psi_T, n_other, which, G_T, denom, bc, aug=None):
        """numpy assembly + ``thomas_batch`` (the reference of ``_solve_carrier``)."""
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
        srh_lin = qv * n_other[1:-1] / denom[1:-1]
        srh_const = qv * n_i2 / denom[1:-1]
        if aug is not None:
            srh_lin = srh_lin + qv * aug[1:-1] * n_other[1:-1]
            srh_const = srh_const + qv * aug[1:-1] * n_i2
        if which == "n":
            c[1:-1] = k[1:] * Bm[1:]
            a[1:-1] = k[:-1] * Bp[:-1]
            b[1:-1] = -(k[1:] * Bp[1:] + k[:-1] * Bm[:-1]) - srh_lin
            rhs[1:-1] = -qv * G_T[1:-1] - srh_const
        else:
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
        G_fn=None,
        srh_scale_fn=None,
    ):
        """Gummel iteration from a converged Poisson state psi (B, N).

        G: optical generation (B, N) or (N,) cm^-3 s^-1 (None = dark).
        G_fn: optional field-dependent generation, ``G_fn(psi_T) -> (N, B)`` cm^-3 s^-1,
        re-evaluated from the present potential at the start of every Gummel iteration
        and added to G (the dark-current model's f_depl generation,
        etna.core.batched.dark). At convergence the generation is that of the converged
        field: the frozen-field fixed point of etna.core.dark_current.
        _resolve_self_consistent. The generation of the last iteration is returned in
        ``TransportResult.G``.
        srh_scale_fn: optional field-dependent SRH enhancement, ``srh_scale_fn(psi_T)
        -> (N, B)``, re-evaluated like G_fn: U_SRH is multiplied by it (tau ->
        tau / scale; the Hurkx (1 + Gamma) of the dark-current model). None = off,
        every expression unchanged (FLASH path bit-identical).
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
        n1, p1 = self.n1[None, :], self.p1[None, :]
        C_n, C_p = self.C_n[None, :], self.C_p[None, :]
        aug = (
            (lambda nn, pp: C_n * nn + C_p * pp)
            if self.auger
            else (lambda nn, pp: None)
        )
        G_opt = G_T
        converged = np.zeros(B, dtype=bool)
        J_prev = None
        it = 0
        for it in range(1, max_iter + 1):
            if G_fn is not None:
                G_T = G_opt + G_fn(psi_T)
            if srh_scale_fn is not None:
                inv_scale = 1.0 / srh_scale_fn(psi_T)
                tau_n, tau_p = (
                    self.tau_n[None, :] * inv_scale,
                    self.tau_p[None, :] * inv_scale,
                )
            denom = tau_p * (n_T + n1) + tau_n * (p_T + p1)
            n_new = self._solve_carrier(
                psi_T, p_T, "n", G_T, denom, self.n_bc, aug(n_T, p_T)
            )
            n_new = np.maximum(n_new, 0.0)
            denom = tau_p * (n_new + n1) + tau_n * (p_T + p1)
            p_new = self._solve_carrier(
                psi_T, n_new, "p", G_T, denom, self.p_bc, aug(n_new, p_T)
            )
            p_new = np.maximum(p_new, 0.0)
            # cathode current of the continuity-consistent state (psi_T, n_new, p_new): it is a
            # tiny difference of two large SG terms, so it is NOT re-evaluated after the Poisson
            # update below
            Jn0, JpL = _end_currents(
                np.ascontiguousarray(psi_T), n_new, p_new, self.kn, self.kp, pb.V_t
            )  # psi_T is F-ordered after a Poisson update with B > 1: one signature
            # Balanced form (as devsim's extract_steady_state_current): summing the electron
            # equation over the interior nodes gives Jn(cathode edge) = Jn(anode edge) +
            # sum q v (U - G); both terms are well conditioned (electrons are minority at the
            # p+ anode), unlike Jn[-1] itself, which is a difference of huge SG terms when n
            # is large at the cathode (n+ layer). Holes are minority at the cathode: direct.
            denom = tau_p * (n_new + n1) + tau_n * (p_new + p1)
            U = (n_new * p_new - n_i**2) / denom
            if self.auger:
                U = U + (C_n * n_new + C_p * p_new) * (n_new * p_new - n_i**2)
            S = (Q * self.v[1:-1] * (U[1:-1] - G_T[1:-1])).sum(axis=0)
            J = Jn0 + S + JpL
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
                carriers=False,
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
            G=G_T.T.copy(),
        )
