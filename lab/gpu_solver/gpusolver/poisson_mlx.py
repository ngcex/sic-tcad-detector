"""Batched 1D quasi-Fermi Poisson -- MLX / Metal float32 backend.

Same discretisation and interface as poisson_np.PoissonBatch (see there for
the equations). Chosen by measurement (notes/spike_mlx_tridiag.md): the only
MLX path that beats numpy float64 is a custom Metal kernel running one
sequential Thomas sweep per GPU thread on an (N, B) layout, with the Newton
assembly fused by mx.compile. Everything on the GPU is float32: the Metal
backend has no float64.

Float32 potential: |u| = |psi/V_t| reaches ~2400 at -60 V, so one ulp is
~2.4e-4 V_t (6 uV). The potential is therefore carried as a float32 pair
(u_hi + u_lo): u_hi is frozen at the start of each solve and the Newton
updates accumulate in u_lo, so neighbour differences u_{i+1}-u_i (the flux)
keep ~1e-7 relative precision instead of ~1e-7 * 2400. Measured
(notes/ablate_f32_split.py, default device, 0..-60 V): without the split the
Newton update stalls between 1e-4 and 1e-3 V_t (tol=1e-4 never converges) and
psi differs from float64 by 4.2e-6 V; with the split tol=1e-4 converges and
psi agrees to 3.2e-7 V, W to 6e-7 relative, E to 1e-5 of E_max.
"""

from __future__ import annotations

import mlx.core as mx
import numpy as np

from .poisson_np import EXP_CLAMP, PoissonBatch, PoissonResult, contact_values

_THOMAS_NB_SRC = r"""
    uint sys = thread_position_in_grid.x;
    uint Bn = (uint)B[0];
    uint Nn = (uint)N[0];
    if (sys >= Bn) return;
    T bprev = b[sys];
    T cp_prev = c[sys] / bprev;
    T dp_prev = d[sys] / bprev;
    cp[sys] = cp_prev;
    x[sys] = dp_prev;
    for (uint i = 1; i < Nn; ++i) {
        uint k = i * Bn + sys;
        T m = b[k] - a[k] * cp_prev;
        T cpi = c[k] / m;
        T dpi = (d[k] - a[k] * dp_prev) / m;
        cp[k] = cpi;
        x[k] = dpi;
        cp_prev = cpi;
        dp_prev = dpi;
    }
    for (int i = (int)Nn - 2; i >= 0; --i) {
        uint k = (uint)i * Bn + sys;
        x[k] = x[k] - cp[k] * x[k + Bn];
    }
"""

_thomas_kernel = mx.fast.metal_kernel(
    name="thomas_NB",
    input_names=["a", "b", "c", "d", "B", "N"],
    output_names=["x", "cp"],
    source=_THOMAS_NB_SRC,
)


def thomas_mlx(a, b, c, d):
    """Solve B tridiagonal systems given as (N, B) float32 arrays; one GPU thread each."""
    N, B = b.shape
    x, _cp = _thomas_kernel(
        inputs=[
            a,
            b,
            c,
            d,
            mx.array([B], dtype=mx.int32),
            mx.array([N], dtype=mx.int32),
        ],
        template=[("T", mx.float32)],
        grid=(B, 1, 1),
        threadgroup=(min(256, B), 1, 1),
        output_shapes=[b.shape, b.shape],
        output_dtypes=[mx.float32, mx.float32],
    )
    return x


class PoissonBatchMLX(PoissonBatch):
    def __init__(self, x, donors, acceptors, n_i, V_t):
        super().__init__(x, donors, acceptors, n_i, V_t)
        f32 = lambda arr: mx.array(np.asarray(arr, dtype=np.float32))
        N, B = self.N, self.B
        self._lam2 = f32(self.lam2[None, :])  # (1, B)
        self._lnni = f32(self.log_ni_hat[None, :])
        self._D = f32(self.D_hat)  # (N, B)
        self._inv_h_l = f32(self.inv_h[:-1])  # rows 1..N-2 : 1/h_{i-1}
        self._inv_h_r = f32(self.inv_h[1:])  # rows 1..N-2 : 1/h_i
        self._inv_h = f32(self.inv_h)  # (N-1, 1)
        self._v = f32(self.v[1:-1])  # (N-2, 1)
        self._zero = mx.zeros((1, B), dtype=mx.float32)
        self._one = mx.ones((1, B), dtype=mx.float32)
        self._Vt = f32(self.V_t[None, :])
        self._step = mx.compile(self._newton_step)

    def _newton_step(self, u_hi, u_lo, u_n, u_p, active, clamp):
        u = u_hi + u_lo
        n_hat = mx.exp(mx.minimum(u - u_n + self._lnni, EXP_CLAMP))
        p_hat = mx.exp(mx.minimum(u_p - u + self._lnni, EXP_CLAMP))
        # flux from the split representation: hi differences are exact in f32
        du = (u_hi[1:] - u_hi[:-1]) + (u_lo[1:] - u_lo[:-1])
        flux = du * self._inv_h
        F = self._lam2 * (flux[1:] - flux[:-1]) + self._v * (
            p_hat[1:-1] - n_hat[1:-1] + self._D[1:-1]
        )
        a = mx.concatenate([self._zero, self._lam2 * self._inv_h_l, self._zero])
        c = mx.concatenate([self._zero, self._lam2 * self._inv_h_r, self._zero])
        b = mx.concatenate(
            [
                self._one,
                -self._lam2 * (self._inv_h_l + self._inv_h_r)
                - self._v * (p_hat[1:-1] + n_hat[1:-1]),
                self._one,
            ]
        )
        rhs = mx.concatenate([self._zero, -F, self._zero])
        delta = mx.clip(thomas_mlx(a, b, c, rhs), -clamp, clamp) * active
        upd = mx.max(mx.abs(delta), axis=0)
        return u_lo + delta, upd

    def solve(
        self,
        psi0,
        V_anode=0.0,
        V_cathode=0.0,
        tol=1e-4,
        max_iter=200,
        clamp=20.0,
        phi_n=None,
        phi_p=None,
    ):
        B, N = self.B, self.N
        V_a = np.broadcast_to(np.asarray(V_anode, dtype=float), (B,))
        V_c = np.broadcast_to(np.asarray(V_cathode, dtype=float), (B,))
        psi_a, psi_c = contact_values(
            self.netdop[:, 0], self.netdop[:, -1], self.n_i, self.V_t, V_a, V_c
        )
        u64 = (np.asarray(psi0, dtype=float) / self.V_t[:, None]).T.copy()
        u64[0], u64[-1] = psi_a / self.V_t, psi_c / self.V_t
        u_hi = mx.array(u64.astype(np.float32))
        u_lo = mx.array((u64 - np.asarray(u_hi)).astype(np.float32))
        u_n64, u_p64 = self.quasi_fermi(V_a, V_c, phi_n, phi_p)
        u_n = mx.array(np.ascontiguousarray(u_n64, dtype=np.float32))
        u_p = mx.array(np.ascontiguousarray(u_p64, dtype=np.float32))
        active = mx.ones((1, B), dtype=mx.float32)
        clamp_a = mx.array(np.float32(clamp))
        converged = np.zeros(B, dtype=bool)
        max_upd = np.inf
        it = 0
        for it in range(1, max_iter + 1):
            u_lo, upd = self._step(u_hi, u_lo, u_n, u_p, active, clamp_a)
            upd_np = np.asarray(upd)  # one GPU sync per Newton iteration
            converged |= upd_np < tol
            max_upd = float(upd_np[~converged].max()) if not converged.all() else 0.0
            if converged.all():
                break
            active = mx.array((~converged).astype(np.float32)[None, :])
        # final densities / potential back in float64 numpy
        u = np.asarray(u_hi).astype(np.float64) + np.asarray(u_lo).astype(np.float64)
        lnni = self.log_ni_hat[None, :]
        n = (self.N0[None, :] * np.exp(np.minimum(u - u_n64 + lnni, EXP_CLAMP))).T
        p = (self.N0[None, :] * np.exp(np.minimum(u_p64 - u + lnni, EXP_CLAMP))).T
        psi = (u * self.V_t[None, :]).T
        return PoissonResult(
            psi=psi, n=n, p=p, iterations=it, converged=converged, max_update=max_upd
        )
