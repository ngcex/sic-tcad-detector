"""GPU (Apple Metal via MLX, float32) red-black Gauss-Seidel smoother vs the numba f64 CPU
smoother of etna.core.fv3d.mg, and the accuracy cost of an f32 fine-level smoother inside
the f64 MG-preconditioned CG.

MLX is NOT an ETNA dependency; run with
    uv run --with mlx python lab/fv3d/bench/gpu_smoother.py

Measured: (a) time of one full RB-GS sweep (both colours) on n^3 grids, MLX f32 (Metal)
vs numba f64 (10 CPU threads, parallel prange); (b) PCG iterations to 1e-10 and final
true residual with the finest-level smoother on the GPU in f32 (casts per call) vs all-f64.
"""

import time

import numpy as np

from etna.core.fv3d.mg import Multigrid, pcg, rbgs

import mlx.core as mx


def problem(n):
    xs = np.linspace(0, 1, n)
    Z, Y, X = np.meshgrid(xs, xs, xs, indexing="ij")
    free = ~((np.hypot(X, Y) <= 0.2) & (Z <= 0.6) | (Z >= 1.0 - 1e-12))
    return xs, free.astype(np.uint8)


class MLXSmoother:
    """RB-GS on the fine level in float32 on the GPU (masked full-array updates)."""

    def __init__(self, L):
        f = mx.float32
        wx, wy, wz, ihx, ihy, ihz = L.geo
        nz, ny, nx = L.shape
        # couplings to the -/+ neighbour along each axis, 0 outside the box
        ihx_m = np.r_[0.0, ihx]
        ihx_p = np.r_[ihx, 0.0]
        ihy_m = np.r_[0.0, ihy]
        ihy_p = np.r_[ihy, 0.0]
        ihz_m = np.r_[0.0, ihz]
        ihz_p = np.r_[ihz, 0.0]
        ayz = wz[:, None, None] * wy[None, :, None]
        axz = wz[:, None, None] * wx[None, None, :]
        axy = wy[None, :, None] * wx[None, None, :]
        self.c = [mx.array(a.astype(np.float32)) for a in (
            ayz * ihx_m[None, None, :], ayz * ihx_p[None, None, :],
            axz * ihy_m[None, :, None], axz * ihy_p[None, :, None],
            axy * ihz_m[:, None, None], axy * ihz_p[:, None, None])]  # fmt: skip
        k, j, i = np.meshgrid(
            np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij"
        )
        par = (i + j + k) % 2
        fr = L.free != 0
        self.mask = [mx.array(fr & (par == col)) for col in (0, 1)]
        self.L = L

    def sweep(self, u, b, idiag, order=(0, 1)):
        gxm, gxp, gym, gyp, gzm, gzp = self.c
        for col in order:
            up = mx.pad(u, [(1, 1), (1, 1), (1, 1)])
            s = (b + gxm * up[1:-1, 1:-1, :-2] + gxp * up[1:-1, 1:-1, 2:]
                 + gym * up[1:-1, :-2, 1:-1] + gyp * up[1:-1, 2:, 1:-1]
                 + gzm * up[:-2, 1:-1, 1:-1] + gzp * up[2:, 1:-1, 1:-1])  # fmt: skip
            u = mx.where(self.mask[col], s * idiag, u)
        return u


def time_sweeps(n, reps=10):
    xs, free = problem(n)
    mg = Multigrid(xs, xs, xs, free)
    mg.set_c(np.zeros(free.shape))
    L = mg.fine
    rng = np.random.default_rng(0)
    b = rng.standard_normal(free.shape) * (free != 0)
    u = np.zeros(free.shape)
    rbgs(u, b, *L.geo, L.diag, L.free, 0)  # warm-up
    t0 = time.perf_counter()
    for _ in range(reps):
        rbgs(u, b, *L.geo, L.diag, L.free, 0)
        rbgs(u, b, *L.geo, L.diag, L.free, 1)
    tc = (time.perf_counter() - t0) / reps
    S = MLXSmoother(L)
    ub, bb = mx.zeros(free.shape, dtype=mx.float32), mx.array(b.astype(np.float32))
    idg = mx.array((1.0 / L.diag).astype(np.float32))
    mx.eval(S.sweep(ub, bb, idg))
    t0 = time.perf_counter()
    for _ in range(reps):
        ub = S.sweep(ub, bb, idg)
    mx.eval(ub)
    tg = (time.perf_counter() - t0) / reps
    # agreement of one sweep from 0 (f32 vs f64)
    u64 = np.zeros(free.shape)
    rbgs(u64, b, *L.geo, L.diag, L.free, 0)
    rbgs(u64, b, *L.geo, L.diag, L.free, 1)
    u32 = np.array(S.sweep(mx.zeros(free.shape, dtype=mx.float32), bb, idg))
    rel = np.abs(u32 - u64).max() / np.abs(u64).max()
    print(f"n {n}^3 ({n**3:,} nodes): RB-GS sweep numba f64 {tc * 1e3:8.1f} ms, MLX f32 {tg * 1e3:8.1f} ms "
          f"(speed-up {tc / tg:5.2f}x), one-sweep f32 rel. diff {rel:.1e}", flush=True)  # fmt: skip
    return mg, S, b


class MixedMG:
    """The fv3d V-cycle with the finest-level pre/post smoothing done by MLX in f32."""

    def __init__(self, mg, S):
        self.mg, self.S = mg, S
        L = mg.fine
        self.idg = mx.array((1.0 / L.diag).astype(np.float32))
        self.fine = L

    def vcycle(self, b, lev=0):
        mg = self.mg
        if lev > 0:
            return mg.vcycle(b, lev)
        L = mg.levels[0]
        bb = mx.array(b.astype(np.float32))
        u = mx.zeros(L.shape, dtype=mx.float32)
        for _ in range(mg.nu[0]):
            u = self.S.sweep(u, bb, self.idg, (0, 1))
        u = np.array(u).astype(np.float64)
        r = b - L.apply(u)
        ec = mg.vcycle(L.restrict(r) * (L.coarse.free != 0), 1)
        from etna.core.fv3d.mg import prolong_add

        (kx0, kx1, wx0, wx1), (ky0, ky1, wy0, wy1), (kz0, kz1, wz0, wz1) = L.maps
        prolong_add(
            ec, u, kx0, kx1, wx0, wx1, ky0, ky1, wy0, wy1, kz0, kz1, wz0, wz1, L.free
        )
        uu = mx.array(u.astype(np.float32))
        for _ in range(mg.nu[1]):
            uu = self.S.sweep(uu, bb, self.idg, (1, 0))
        return np.array(uu).astype(np.float64) * (L.free != 0)


def pcg_compare(mg, S, b):
    t0 = time.perf_counter()
    x, it, rr = pcg(mg, b, rtol=1e-10, maxiter=300)
    t64 = time.perf_counter() - t0
    mm = MixedMG(mg, S)
    t0 = time.perf_counter()
    x2, it2, rr2 = pcg(mm, b, rtol=1e-10, maxiter=300)
    t32 = time.perf_counter() - t0
    tr = np.linalg.norm(b - mg.fine.apply(x2)) / np.linalg.norm(b)
    print(f"   PCG to 1e-10: all-f64 {it} it {t64:.2f} s | f32 GPU fine smoother {it2} it {t32:.2f} s, "
          f"true rel. residual {tr:.1e}, max|x32 - x64|/max|x64| {np.abs(x2 - x).max() / np.abs(x).max():.1e}")  # fmt: skip


if __name__ == "__main__":
    print("MLX", mx.__version__, mx.default_device())
    for n in (64, 128, 192, 256):
        mg, S, b = time_sweeps(n)
        if n in (128, 192):
            pcg_compare(mg, S, b)
