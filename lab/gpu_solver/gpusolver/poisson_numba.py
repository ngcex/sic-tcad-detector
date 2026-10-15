"""Batched 1D quasi-Fermi Poisson -- numba float64 CPU backend (prange over devices).

Same equations and interface as poisson_np.PoissonBatch. Each device runs its
complete Newton loop (assembly, Thomas sweep, clamped update, convergence
test) inside one compiled function; devices are independent, so the batch is
split over CPU cores with ``prange``. Unlike the vectorised numpy/MLX
backends there is no masked waste: a device stops iterating as soon as it has
converged. Arrays are (B, N) so one device's data is contiguous.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange

from .poisson_np import EXP_CLAMP, PoissonBatch, PoissonResult, contact_values


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


class PoissonBatchNumba(PoissonBatch):
    def __init__(self, x, donors, acceptors, n_i, V_t):
        super().__init__(x, donors, acceptors, n_i, V_t)
        self._D_BN = np.ascontiguousarray(self.D_hat.T)
        self._inv_h1 = np.ascontiguousarray(self.inv_h[:, 0])
        self._v1 = np.ascontiguousarray(self.v[:, 0])

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
        B, N = self.B, self.N
        V_a = np.broadcast_to(np.asarray(V_anode, dtype=float), (B,))
        V_c = np.broadcast_to(np.asarray(V_cathode, dtype=float), (B,))
        u_n_T, u_p_T = self.quasi_fermi(V_a, V_c, phi_n, phi_p)
        u_n = np.ascontiguousarray(np.broadcast_to(u_n_T.T, (B, N)))
        u_p = np.ascontiguousarray(np.broadcast_to(u_p_T.T, (B, N)))
        psi_a, psi_c = contact_values(
            self.netdop[:, 0], self.netdop[:, -1], self.n_i, self.V_t, V_a, V_c
        )
        u = np.ascontiguousarray(np.asarray(psi0, dtype=float) / self.V_t[:, None])
        u[:, 0] = psi_a / self.V_t
        u[:, -1] = psi_c / self.V_t
        iters = np.zeros(B, dtype=np.int64)
        conv = np.zeros(B, dtype=np.bool_)
        _newton_batch(
            u,
            u_n,
            u_p,
            self._D_BN,
            self.log_ni_hat,
            self.lam2,
            self._inv_h1,
            self._v1,
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
            psi=u * self.V_t[:, None],
            n=n,
            p=p,
            iterations=int(iters.max()),
            converged=conv,
            max_update=np.nan,
        )
