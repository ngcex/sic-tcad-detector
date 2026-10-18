"""Batched 1D quasi-Fermi Poisson for p+/n- diodes (numba float64, ``prange`` over devices).

Solves, for B independent diodes sharing one mesh x (N nodes),

    d/dx( eps dpsi/dx ) = -q ( p - n + N_D - N_A )
    n = n_i exp((psi - phi_n)/V_t),  p = n_i exp((phi_p - psi)/V_t)

with devsim's ohmic-contact Dirichlet values at both ends. With flat
quasi-Fermi levels phi_n = V_cathode, phi_p = V_anode this is ETNA's
Poisson-only equilibrium at 0 V and, under reverse bias, the electrostatics of
the drift-diffusion solution (currents are negligible), which is all that
C-V / depletion width / field extraction need. Nodal phi_n/phi_p (Gummel
coupling) give the transport stage in ``transport.py``.

Scaled unknown u = psi/V_t; densities scaled by N0 = max |NetDoping| per
device. Residual at interior node i (h_i = x_{i+1}-x_i, v_i = (h_i+h_{i-1})/2):

    F_i = lam2 [ (u_{i+1}-u_i)/h_i - (u_i-u_{i-1})/h_{i-1} ] + v_i (p^ - n^ + D^)
    lam2 = eps V_t / (q N0)

Each device runs its whole Newton loop (assembly, Thomas sweep, clamped
update, convergence test) in one compiled function and stops as soon as it
has converged.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit, prange

from etna.core.batched.material import EPS, Q

EXP_CLAMP = 60.0  # exp(60) ~ 1e26 : above any physical n^/p^


@dataclass
class PoissonResult:
    psi: np.ndarray  # (B, N) V
    n: np.ndarray  # (B, N) cm^-3
    p: np.ndarray  # (B, N) cm^-3
    iterations: int
    converged: np.ndarray  # (B,) bool


def contact_values(netdop_0, netdop_N, n_i, V_t, V_anode, V_cathode):
    """devsim ohmic contact potentials (V) at node 0 (anode, p+) and node N-1."""
    celec = lambda D: 1e-10 + 0.5 * np.abs(D + np.sqrt(D**2 + 4 * n_i**2))
    chole = lambda D: 1e-10 + 0.5 * np.abs(-D + np.sqrt(D**2 + 4 * n_i**2))
    psi_a = np.where(
        netdop_0 > 0,
        V_anode + V_t * np.log(celec(netdop_0) / n_i),
        V_anode - V_t * np.log(chole(netdop_0) / n_i),
    )
    psi_c = np.where(
        netdop_N > 0,
        V_cathode + V_t * np.log(celec(netdop_N) / n_i),
        V_cathode - V_t * np.log(chole(netdop_N) / n_i),
    )
    return psi_a, psi_c


@njit(cache=True, parallel=True)
def _newton_batch(
    u, u_n, u_p, D, lnni, lam2, inv_h, v, tol, max_iter, clamp, iters, conv
):
    B, N = u.shape
    for bi in prange(B):
        a = np.empty(N)
        b = np.empty(N)
        c = np.empty(N)
        F = np.empty(N)
        cp = np.empty(N)
        dp = np.empty(N)
        l2 = lam2[bi]
        ln = lnni[bi]
        converged = False
        it = 0
        for it in range(1, max_iter + 1):
            # assembly (interior rows); boundary rows are identity with zero rhs
            for i in range(1, N - 1):
                ui = u[bi, i]
                en = min(ui - u_n[bi, i] + ln, EXP_CLAMP)
                ep = min(u_p[bi, i] - ui + ln, EXP_CLAMP)
                n_hat = np.exp(en)
                p_hat = np.exp(ep)
                fl = (u[bi, i + 1] - ui) * inv_h[i]
                fr = (ui - u[bi, i - 1]) * inv_h[i - 1]
                F[i] = l2 * (fl - fr) + v[i] * (p_hat - n_hat + D[bi, i])
                a[i] = l2 * inv_h[i - 1]
                c[i] = l2 * inv_h[i]
                b[i] = -l2 * (inv_h[i - 1] + inv_h[i]) - v[i] * (p_hat + n_hat)
            F[0] = 0.0
            F[N - 1] = 0.0
            # Thomas on J delta = -F
            cp[0] = 0.0
            dp[0] = 0.0
            for i in range(1, N - 1):
                m = 1.0 / (b[i] - a[i] * cp[i - 1])
                cp[i] = c[i] * m
                dp[i] = (-F[i] - a[i] * dp[i - 1]) * m
            upd = 0.0
            dnext = 0.0  # delta at node N-1 is zero (Dirichlet)
            for i in range(N - 2, 0, -1):
                d = dp[i] - cp[i] * dnext
                if d > clamp:
                    d = clamp
                elif d < -clamp:
                    d = -clamp
                dnext = dp[i] - cp[i] * dnext  # unclamped for the recursion
                u[bi, i] += d
                if abs(d) > upd:
                    upd = abs(d)
            if upd < tol:
                converged = True
                break
        iters[bi] = it
        conv[bi] = converged


@njit(cache=True, parallel=True)
def _dpsi_dvc_batch(n_hat, p_hat, lam2, inv_h, v, out):
    """Small-signal response du/du_c of a converged flat-quasi-Fermi Poisson state.

    Linearises F(u; u_n = u_c) = 0 in the cathode bias u_c = V_cathode/V_t (phi_n
    follows the cathode, phi_p the anode): J du = -dF/du_c, dF_i/du_c = v_i n^_i, with
    du_0 = 0 (anode) and du_{N-1} = 1 (cathode Dirichlet). J is the Newton Jacobian.
    Then out[b] = sum_i v_i p^_i du_i, the (scaled) hole charge removed per unit u_c.
    """
    B, N = n_hat.shape
    for bi in prange(B):
        l2 = lam2[bi]
        cp = np.empty(N)
        dp = np.empty(N)
        du = np.empty(N)
        cp[0] = 0.0
        dp[0] = 0.0  # du_0 = 0
        for i in range(1, N - 1):
            a = l2 * inv_h[i - 1]
            c = l2 * inv_h[i]
            b = -l2 * (inv_h[i - 1] + inv_h[i]) - v[i] * (p_hat[bi, i] + n_hat[bi, i])
            r = -v[i] * n_hat[bi, i]
            if i == N - 2:
                r -= c  # du_{N-1} = 1 moved to the rhs
                c = 0.0
            m = 1.0 / (b - a * cp[i - 1])
            cp[i] = c * m
            dp[i] = (r - a * dp[i - 1]) * m
        du[N - 1] = 1.0
        du[0] = 0.0
        nxt = 0.0
        for i in range(N - 2, 0, -1):
            du[i] = dp[i] - cp[i] * nxt
            nxt = du[i]
        s = 0.0
        for i in range(1, N - 1):
            s += v[i] * p_hat[bi, i] * du[i]
        out[bi] = s


class PoissonBatch:
    """Fixed per-batch data; call solve() repeatedly along a bias ramp."""

    def __init__(self, x, donors, acceptors, n_i, V_t):
        self.x = np.asarray(x, dtype=float)
        self.N = len(self.x)
        D = np.asarray(donors, dtype=float) - np.asarray(acceptors, dtype=float)
        self.B = D.shape[0]
        self.n_i = np.broadcast_to(np.asarray(n_i, dtype=float), (self.B,)).copy()
        self.V_t = np.broadcast_to(np.asarray(V_t, dtype=float), (self.B,)).copy()
        self.N0 = np.abs(D).max(axis=1)  # (B,)
        self.donors = np.asarray(donors, dtype=float)
        # ionized acceptors (B, N), as devsim's Acceptors node model (local mobility)
        self.acceptors = np.broadcast_to(
            np.asarray(acceptors, dtype=float), D.shape
        )
        self.netdop = D
        self.lam2 = EPS * self.V_t / (Q * self.N0)  # (B,) cm^2
        self.log_ni_hat = np.log(self.n_i / self.N0)  # (B,)
        h = np.diff(self.x)
        v = np.empty(self.N)
        v[1:-1] = 0.5 * (h[1:] + h[:-1])
        v[0] = v[-1] = 1.0
        self.v = v  # (N,) control volumes, ends unused
        self._D_hat = np.ascontiguousarray(D / self.N0[:, None])  # (B, N)
        self._inv_h = 1.0 / h  # (N-1,)

    def neutral_guess(self, V_anode=0.0, V_cathode=0.0):
        """Charge-neutral potential (V), (B, N) -- the standard Newton starting point."""
        u = np.arcsinh(self.netdop / (2 * self.n_i[:, None]))
        V_a = np.broadcast_to(np.asarray(V_anode, float), (self.B,))[:, None]
        V_c = np.broadcast_to(np.asarray(V_cathode, float), (self.B,))[:, None]
        frac = ((self.x - self.x[0]) / (self.x[-1] - self.x[0]))[None, :]
        return u * self.V_t[:, None] + V_a + (V_c - V_a) * frac

    def solve(
        self,
        psi0,
        V_anode=0.0,
        V_cathode=0.0,
        tol=1e-8,
        max_iter=200,
        clamp=20.0,
        phi_n=None,
        phi_p=None,
    ):
        """Newton-solve all B devices at biases V_anode/V_cathode ((B,) or scalar, V).

        phi_n/phi_p: optional nodal quasi-Fermi levels (B, N) in V (Gummel
        coupling); default flat at the contact biases. tol is on max |delta u|
        in V_t units.
        """
        B, N = self.B, self.N
        V_a = np.broadcast_to(np.asarray(V_anode, dtype=float), (B,))
        V_c = np.broadcast_to(np.asarray(V_cathode, dtype=float), (B,))
        Vt = self.V_t[:, None]
        u_n = np.ascontiguousarray(
            np.broadcast_to((V_c[:, None] if phi_n is None else phi_n) / Vt, (B, N))
        )
        u_p = np.ascontiguousarray(
            np.broadcast_to((V_a[:, None] if phi_p is None else phi_p) / Vt, (B, N))
        )
        psi_a, psi_c = contact_values(
            self.netdop[:, 0], self.netdop[:, -1], self.n_i, self.V_t, V_a, V_c
        )
        u = np.ascontiguousarray(np.asarray(psi0, dtype=float) / Vt)
        u[:, 0] = psi_a / self.V_t
        u[:, -1] = psi_c / self.V_t
        iters = np.zeros(B, dtype=np.int64)
        conv = np.zeros(B, dtype=np.bool_)
        _newton_batch(
            u,
            u_n,
            u_p,
            self._D_hat,
            self.log_ni_hat,
            self.lam2,
            self._inv_h,
            self.v,
            tol,
            max_iter,
            clamp,
            iters,
            conv,
        )
        lnni = self.log_ni_hat[:, None]
        n = self.N0[:, None] * np.exp(np.minimum(u - u_n + lnni, EXP_CLAMP))
        p = self.N0[:, None] * np.exp(np.minimum(u_p - u + lnni, EXP_CLAMP))
        return PoissonResult(
            psi=u * Vt, n=n, p=p, iterations=int(iters.max()), converged=conv
        )

    def small_signal_capacitance(self, res):
        """Small-signal C = dQ/dV (F/cm^2, area 1) of each device at a converged state.

        ``res`` is a flat-quasi-Fermi PoissonResult (dark reverse bias). Q is the
        terminal charge: holes leave only through the anode and electrons only through
        the cathode (generation/recombination are negligible in the quasi-static
        limit), so dQ/dV = q d(sum_i v_i p_i)/dV_cathode, evaluated exactly (no finite
        difference) from one linear solve with the Newton Jacobian. This is the C of a
        measured (low-frequency) C-V, i.e. what W = eps*A/C refers to. It agrees with
        devsim's small-signal AC solve to ~1e-8 (lab/microdosimeter/bench/cv_gate.py,
        tests/test_cv_small_signal.py); it is NOT eps/W_edge of the n/N_D = 0.5 edge.
        """
        n_hat = np.ascontiguousarray(res.n / self.N0[:, None])
        p_hat = np.ascontiguousarray(res.p / self.N0[:, None])
        out = np.empty(self.B)
        _dpsi_dvc_batch(n_hat, p_hat, self.lam2, self._inv_h, self.v, out)
        return Q * self.N0 * out / self.V_t
