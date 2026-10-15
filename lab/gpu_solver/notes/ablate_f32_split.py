"""Ablation: MLX float32 with and without the u_hi/u_lo potential split, vs numpy float64.

Run from lab/gpu_solver: uv run python notes/ablate_f32_split.py (~20 s). Result is
quoted in gpusolver/poisson_mlx.py and README.md.
"""

import numpy as np, mlx.core as mx, time
from gpusolver.sweep import make_batch, reverse_sweep
from gpusolver import poisson_mlx
from gpusolver.poisson_np import EXP_CLAMP

import os
X = np.load(os.path.join(os.path.dirname(__file__), "..", "bench", "mesh_devsim_default.npy"))
XJ, EPI = 1e-4, 10e-4
V = np.arange(0.0, -61.0, -5.0)
ref = reverse_sweep(make_batch(X, XJ, 1e19, 2.93e15, 8.82e13, 0.987e-4), V, XJ, EPI)


def report(label, sw):
    dW = np.abs(sw.W - ref.W).max() / ref.W.min()
    dpsi = np.abs(sw.psi - ref.psi).max()
    E = -(sw.psi[..., 1:] - sw.psi[..., :-1]) / np.diff(X)
    Er = -(ref.psi[..., 1:] - ref.psi[..., :-1]) / np.diff(X)
    dE = np.abs(E - Er).max() / np.abs(Er).max()
    print(
        f"{label:28s}: its {sw.newton_iterations:4d} conv {sw.converged.all()!s:5s} relW {dW:.2e} |dpsi| {dpsi:.2e} V relE {dE:.2e}"
    )


for tol in (1e-3, 1e-4, 1e-5):
    pb = make_batch(X, XJ, 1e19, 2.93e15, 8.82e13, 0.987e-4, backend="mlx")
    report(f"split      tol={tol:g}", reverse_sweep(pb, V, XJ, EPI, tol=tol))


# ablation: single float32 u (no split): fold u_lo into u_hi each step
def _newton_step_nosplit(self, u_hi, u_lo, u_n, u_p, active, clamp):
    u = u_hi + u_lo
    n_hat = mx.exp(mx.minimum(u - u_n + self._lnni, EXP_CLAMP))
    p_hat = mx.exp(mx.minimum(u_p - u + self._lnni, EXP_CLAMP))
    flux = (u[1:] - u[:-1]) * self._inv_h
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
    delta = mx.clip(poisson_mlx.thomas_mlx(a, b, c, rhs), -clamp, clamp) * active
    return (u + delta) - u_hi, mx.max(
        mx.abs(delta), axis=0
    )  # u_lo absorbs everything -> u = u_hi+u_lo rounds as single f32


poisson_mlx.PoissonBatchMLX._newton_step = _newton_step_nosplit
for tol in (1e-3, 1e-4, 1e-5):
    pb = make_batch(X, XJ, 1e19, 2.93e15, 8.82e13, 0.987e-4, backend="mlx")
    report(
        f"no-split   tol={tol:g}",
        reverse_sweep(
            pb,
            V,
            XJ,
            EPI,
            tol=tol,
        ),
    )
