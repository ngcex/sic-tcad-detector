"""Spike: can MLX/Metal float32 beat numpy float64 (Accelerate) for the inner loop
of a BATCHED 1D semiconductor Newton solver?

B independent devices, each a tridiagonal system with N unknowns on a nonuniform
(geometrically graded) 1D mesh.  Methods compared:

  1. numpy float64 Thomas, vectorized over B, layouts (B,N) and (N,B)
  2. MLX float32 Thomas, Python loop over N (eager and mx.compile)
  3. MLX float32 parallel cyclic reduction (PCR), log2(N) sweeps (eager and compiled)
  4. MLX custom Metal kernel: one GPU thread per system doing sequential Thomas
     (a) (B,N) layout as specified, (b) transposed (N,B) layout (coalesced reads)
  5. Newton assembly (Poisson-like residual + Jacobian diagonal)
  6. Full Newton iteration = assembly + tridiagonal solve + update

Reproduce:
  cd lab/gpu_solver && uv run python notes/spike_mlx_tridiag.py
Quick mode (smaller batch list, fewer repeats):
  uv run python notes/spike_mlx_tridiag.py --quick
"""

from __future__ import annotations

import argparse
import math
import statistics
import sys
import time

import numpy as np
import mlx.core as mx

# ----------------------------------------------------------------------------
# timing helper
# ----------------------------------------------------------------------------
N_WARM = 2
N_RUNS = 5


def timeit(fn, *, warm=None, runs=None) -> float:
    """Median wall time (ms) of `fn()` after warmups. fn must sync itself."""
    warm = N_WARM if warm is None else warm
    runs = N_RUNS if runs is None else runs
    for _ in range(warm):
        fn()
    ts = []
    for _ in range(runs):
        t0 = time.perf_counter()
        fn()
        ts.append((time.perf_counter() - t0) * 1e3)
    return statistics.median(ts)


# ----------------------------------------------------------------------------
# realistic test matrices
# ----------------------------------------------------------------------------
L_CM = 11e-4
H_MIN = 3e-8
H_MAX = 1e-5


def graded_half(n_int: int, length: float) -> np.ndarray:
    """n_int intervals, geometrically graded from H_MIN up to (capped at) H_MAX,
    with ratio r chosen by bisection so that the intervals sum to `length`."""

    def hs(r):
        return np.minimum(H_MIN * r ** np.arange(n_int), H_MAX)

    lo, hi = 1.0, 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if hs(mid).sum() < length:
            lo = mid
        else:
            hi = mid
    h = hs(0.5 * (lo + hi))
    return h * (length / h.sum())  # exact length


def make_mesh(N: int) -> np.ndarray:
    """N points over L_CM, junction in the middle; returns N-1 spacings h."""
    n_int = N - 1
    nl = n_int // 2
    nr = n_int - nl
    hl = graded_half(nl, L_CM * nl / n_int)[
        ::-1
    ]  # coarse at left edge -> fine at junction
    hr = graded_half(nr, L_CM * nr / n_int)  # fine at junction -> coarse at right edge
    h = np.concatenate([hl, hr])
    assert h.shape == (N - 1,)
    return h


def make_systems(B: int, N: int, seed: int = 0, dmode: str = "spec"):
    """Poisson-like tridiagonal systems, layout (B, N), float64.
    sub[i] = 1/h[i-1], sup[i] = 1/h[i], diag[i] = -(1/h[i-1]+1/h[i]) - d_i,
    dmode="spec":      d_i ~ logspace(-2, 8) random per row (as specified)
    dmode="laplacian": d_i = 1e-2 everywhere (depleted region, ~pure Laplacian)
    dmode="dominant":  d_i = 1e8 everywhere (strongly diagonally dominant)
    Dirichlet rows at both ends."""
    rng = np.random.default_rng(seed)
    h = make_mesh(N)
    inv_h = 1.0 / h
    a = np.zeros((B, N))
    c = np.zeros((B, N))
    b = np.zeros((B, N))
    a[:, 1:-1] = inv_h[:-1]
    c[:, 1:-1] = inv_h[1:]
    if dmode == "spec":
        d_i = 10.0 ** rng.uniform(-2, 8, size=(B, N))
    elif dmode == "laplacian":
        d_i = np.full((B, N), 1e-2)
    elif dmode == "dominant":
        d_i = np.full((B, N), 1e8)
    else:
        raise ValueError(dmode)
    b[:, 1:-1] = -(inv_h[:-1] + inv_h[1:]) - d_i[:, 1:-1]
    b[:, 0] = 1.0
    b[:, -1] = 1.0
    rhs = rng.standard_normal((B, N))
    return a, b, c, rhs, h


def row_scale(a, b, c, d):
    """Divide each row by its diagonal (any array lib with broadcasting)."""
    return a / b, b / b, c / b, d / b


# ----------------------------------------------------------------------------
# 1. numpy float64 Thomas, vectorized over B
# ----------------------------------------------------------------------------
def thomas_np_BN(a, b, c, d):
    """Layout (B, N): loop over N, each step touches a strided column."""
    B, N = b.shape
    cp = np.empty_like(b)
    dp = np.empty_like(b)
    cp[:, 0] = c[:, 0] / b[:, 0]
    dp[:, 0] = d[:, 0] / b[:, 0]
    for i in range(1, N):
        m = b[:, i] - a[:, i] * cp[:, i - 1]
        cp[:, i] = c[:, i] / m
        dp[:, i] = (d[:, i] - a[:, i] * dp[:, i - 1]) / m
    x = np.empty_like(b)
    x[:, -1] = dp[:, -1]
    for i in range(N - 2, -1, -1):
        x[:, i] = dp[:, i] - cp[:, i] * x[:, i + 1]
    return x


def thomas_np_NB(a, b, c, d):
    """Layout (N, B): loop over N, each step touches a contiguous row."""
    N, B = b.shape
    cp = np.empty_like(b)
    dp = np.empty_like(b)
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, N):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = np.empty_like(b)
    x[-1] = dp[-1]
    for i in range(N - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


# ----------------------------------------------------------------------------
# 2. MLX float32 Thomas with Python loop over N (layout (N, B))
# ----------------------------------------------------------------------------
def thomas_mx_loop(a, b, c, d):
    N = b.shape[0]
    cp = [None] * N
    dp = [None] * N
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, N):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = [None] * N
    x[-1] = dp[-1]
    for i in range(N - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return mx.stack(x, axis=0)


# ----------------------------------------------------------------------------
# 3. MLX float32 parallel cyclic reduction, layout (B, N), N padded to 2^k
# ----------------------------------------------------------------------------
def pad_pow2(a, b, c, d):
    B, N = b.shape
    N2 = 1 << (N - 1).bit_length()
    if N2 == N:
        return a, b, c, d, N
    p = N2 - N
    z = mx.zeros((B, p), dtype=b.dtype)
    o = mx.ones((B, p), dtype=b.dtype)
    return (
        mx.concatenate([a, z], 1),
        mx.concatenate([b, o], 1),
        mx.concatenate([c, z], 1),
        mx.concatenate([d, z], 1),
        N,
    )


def shift_right(v, s):  # out[i] = v[i-s], zero for i<s
    B = v.shape[0]
    return mx.concatenate([mx.zeros((B, s), dtype=v.dtype), v[:, :-s]], 1)


def shift_left(v, s):  # out[i] = v[i+s], zero for i>=N-s
    B = v.shape[0]
    return mx.concatenate([v[:, s:], mx.zeros((B, s), dtype=v.dtype)], 1)


def pcr_mx(a, b, c, d):
    """Parallel cyclic reduction. Each sweep eliminates neighbours at distance s.
    Padded rows are identity rows (b=1, a=c=d=0) which do not couple."""
    a, b, c, d, N = pad_pow2(a, b, c, d)
    N2 = b.shape[1]
    s = 1
    while s < N2:
        a_m, b_m, c_m, d_m = (shift_right(v, s) for v in (a, b, c, d))
        a_p, b_p, c_p, d_p = (shift_left(v, s) for v in (a, b, c, d))
        # neighbours outside the domain: treat as identity rows (b=1) -> alpha=gamma=0 since a_m=c_p=0
        b_m = mx.where(b_m == 0, mx.ones_like(b_m), b_m)
        b_p = mx.where(b_p == 0, mx.ones_like(b_p), b_p)
        alpha = -a / b_m
        gamma = -c / b_p
        b = b + alpha * c_m + gamma * a_p
        d = d + alpha * d_m + gamma * d_p
        a = alpha * a_m
        c = gamma * c_p
        s *= 2
    x = d / b
    return x[:, :N]


# ----------------------------------------------------------------------------
# 4. MLX custom Metal kernels: one thread per system, sequential Thomas
# ----------------------------------------------------------------------------
_SRC_BN = r"""
    uint sys = thread_position_in_grid.x;
    uint Bn = (uint)B[0];
    uint Nn = (uint)N[0];
    if (sys >= Bn) return;
    uint off = sys * Nn;
    T bprev = b[off];
    T cp_prev = c[off] / bprev;
    T dp_prev = d[off] / bprev;
    cp[off] = cp_prev;
    x[off] = dp_prev;
    for (uint i = 1; i < Nn; ++i) {
        uint k = off + i;
        T m = b[k] - a[k] * cp_prev;
        T cpi = c[k] / m;
        T dpi = (d[k] - a[k] * dp_prev) / m;
        cp[k] = cpi;
        x[k] = dpi;
        cp_prev = cpi;
        dp_prev = dpi;
    }
    for (int i = (int)Nn - 2; i >= 0; --i) {
        uint k = off + (uint)i;
        x[k] = x[k] - cp[k] * x[k + 1];
    }
"""

# transposed layout (N, B): element (i, sys) at i*B + sys -> coalesced across threads
_SRC_NB = r"""
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

_kern_BN = mx.fast.metal_kernel(
    name="thomas_BN",
    input_names=["a", "b", "c", "d", "B", "N"],
    output_names=["x", "cp"],
    source=_SRC_BN,
)
_kern_NB = mx.fast.metal_kernel(
    name="thomas_NB",
    input_names=["a", "b", "c", "d", "B", "N"],
    output_names=["x", "cp"],
    source=_SRC_NB,
)


def thomas_metal_BN(a, b, c, d):
    B, N = b.shape
    outs = _kern_BN(
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
    return outs[0]


def thomas_metal_NB(a, b, c, d):
    N, B = b.shape
    outs = _kern_NB(
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
    return outs[0]


# ----------------------------------------------------------------------------
# 5. Newton assembly: Poisson-like residual F and Jacobian diagonal J, layout (B, N)
# ----------------------------------------------------------------------------
def assemble_np(psi, h, s, dop):
    """F_i = (psi[i+1]-psi[i])/h[i] - (psi[i]-psi[i-1])/h[i-1] + s*(exp(-psi)-exp(psi)+dop)
    J_i = -(1/h[i] + 1/h[i-1]) - s*(exp(-psi)+exp(psi))   (interior points)"""
    xp = np.exp(psi[:, 1:-1])
    xm = np.exp(-psi[:, 1:-1])
    inv_r = 1.0 / h[1:]
    inv_l = 1.0 / h[:-1]
    F = (
        (psi[:, 2:] - psi[:, 1:-1]) * inv_r
        - (psi[:, 1:-1] - psi[:, :-2]) * inv_l
        + s * (xm - xp + dop[:, 1:-1])
    )
    J = -(inv_r + inv_l) - s * (xm + xp)
    return F, J


def assemble_mx(psi, h, s, dop):
    xp = mx.exp(psi[:, 1:-1])
    xm = mx.exp(-psi[:, 1:-1])
    inv_r = 1.0 / h[1:]
    inv_l = 1.0 / h[:-1]
    F = (
        (psi[:, 2:] - psi[:, 1:-1]) * inv_r
        - (psi[:, 1:-1] - psi[:, :-2]) * inv_l
        + s * (xm - xp + dop[:, 1:-1])
    )
    J = -(inv_r + inv_l) - s * (xm + xp)
    return F, J


# ----------------------------------------------------------------------------
# 6. Full Newton iteration: assembly -> tridiagonal solve -> update
# ----------------------------------------------------------------------------
def newton_iter_np(psi_NB, h, s, dop_NB, inv_h):
    """numpy float64, layout (N, B). Interior rows from assembly, Dirichlet ends."""
    N, B = psi_NB.shape
    psi = psi_NB
    xp = np.exp(psi[1:-1])
    xm = np.exp(-psi[1:-1])
    inv_r = inv_h[1:, None]
    inv_l = inv_h[:-1, None]
    F = np.zeros_like(psi)
    bdiag = np.ones_like(psi)
    F[1:-1] = (
        (psi[2:] - psi[1:-1]) * inv_r
        - (psi[1:-1] - psi[:-2]) * inv_l
        + s * (xm - xp + dop_NB[1:-1])
    )
    bdiag[1:-1] = -(inv_r + inv_l) - s * (xm + xp)
    a = np.zeros_like(psi)
    c = np.zeros_like(psi)
    a[1:-1] = inv_l
    c[1:-1] = inv_r
    dx = thomas_np_NB(a, bdiag, c, -F)
    return psi + dx


def _newton_iter_mx_core(psi, h, s, dop, a_fixed, c_fixed):
    """MLX float32, layout (N, B); a_fixed/c_fixed precomputed (N,B) off-diagonals."""
    inv_h = 1.0 / h
    inv_r = inv_h[1:, None]
    inv_l = inv_h[:-1, None]
    xp = mx.exp(psi[1:-1])
    xm = mx.exp(-psi[1:-1])
    F_int = (
        (psi[2:] - psi[1:-1]) * inv_r
        - (psi[1:-1] - psi[:-2]) * inv_l
        + s * (xm - xp + dop[1:-1])
    )
    J_int = -(inv_r + inv_l) - s * (xm + xp)
    B = psi.shape[1]
    zero = mx.zeros((1, B), dtype=psi.dtype)
    one = mx.ones((1, B), dtype=psi.dtype)
    F = mx.concatenate([zero, F_int, zero], 0)
    bdiag = mx.concatenate([one, J_int, one], 0)
    dx = thomas_metal_NB(a_fixed, bdiag, c_fixed, -F)
    return psi + dx


# ----------------------------------------------------------------------------
# accuracy helper
# ----------------------------------------------------------------------------
def rel_err(x32, x64):
    """max over systems of ||x32-x64||_inf / ||x64||_inf  (both (B,N))."""
    num = np.max(np.abs(np.asarray(x32, dtype=np.float64) - x64), axis=1)
    den = np.max(np.abs(x64), axis=1)
    return float(np.max(num / den))


def reference_solution(a, b, c, d):
    return thomas_np_NB(a.T.copy(), b.T.copy(), c.T.copy(), d.T.copy()).T.copy()


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def fmt(v):
    if v is None:
        return "n/a"
    if isinstance(v, str):
        return v
    return f"{v:.3f}" if v < 10 else f"{v:.1f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument(
        "--skip-mlx-loop-compile-large",
        action="store_true",
        help="skip mx.compile of the N-loop Thomas for N=1000 (slow to trace)",
    )
    args = ap.parse_args()
    global N_RUNS, N_WARM
    if args.quick:
        N_RUNS, N_WARM = 3, 1
        Bs = [1, 100, 1000]
    else:
        Bs = [1, 100, 1000, 10000]
    Ns = [327, 1000]

    print(f"mlx {mx.__version__}  numpy {np.__version__}  device {mx.default_device()}")
    print(
        f"median of {N_RUNS} runs after {N_WARM} warmups; wall time, ms; MLX synced with mx.eval\n"
    )

    # --------------------------------------------------------------- solvers
    methods = [
        "np64 Thomas (B,N)",
        "np64 Thomas (N,B)",
        "mx32 Thomas loop eager",
        "mx32 Thomas loop compiled",
        "mx32 PCR eager",
        "mx32 PCR compiled",
        "mx32 Metal kernel (B,N)",
        "mx32 Metal kernel (N,B)",
    ]
    solve_rows = {}  # (N, B) -> dict method -> ms
    acc_rows = {}  # (N, B) -> dict method -> (raw, scaled32, scaled64)

    for N in Ns:
        for B in Bs:
            a, b, c, d, h = make_systems(B, N)
            x64 = reference_solution(a, b, c, d)
            res = {}
            acc = {}

            # numpy (B,N)
            res[methods[0]] = timeit(lambda: thomas_np_BN(a, b, c, d))
            # numpy (N,B)
            aT, bT, cT, dT = (v.T.copy() for v in (a, b, c, d))
            res[methods[1]] = timeit(lambda: thomas_np_NB(aT, bT, cT, dT))
            assert np.allclose(thomas_np_BN(a, b, c, d), x64)

            # float32 copies: raw, scaled in float32, scaled in float64 then cast
            A32 = [mx.array(v.astype(np.float32)) for v in (a, b, c, d)]
            A32s = list(row_scale(*A32))
            A32s64 = [mx.array(v.astype(np.float32)) for v in row_scale(a, b, c, d)]
            mx.eval(*A32, *A32s, *A32s64)
            A32T = [mx.array(v.T.copy().astype(np.float32)) for v in (a, b, c, d)]
            A32sT = list(row_scale(*A32T))
            A32s64T = [
                mx.array(v.T.copy().astype(np.float32)) for v in row_scale(a, b, c, d)
            ]
            mx.eval(*A32T, *A32sT, *A32s64T)

            def acc3(fn, raw, s32, s64, transpose=False):
                out = []
                for inp in (raw, s32, s64):
                    x = fn(*inp)
                    mx.eval(x)
                    xn = np.asarray(x, dtype=np.float64)
                    if transpose:
                        xn = xn.T
                    out.append(rel_err(xn, x64))
                return tuple(out)

            # MLX loop eager
            def run_loop():
                x = thomas_mx_loop(*A32T)
                mx.eval(x)
                return x

            res[methods[2]] = timeit(run_loop)
            acc[methods[2]] = acc3(thomas_mx_loop, A32T, A32sT, A32s64T, transpose=True)

            # MLX loop compiled
            if N >= 1000 and args.skip_mlx_loop_compile_large:
                res[methods[3]] = None
            else:
                loop_c = mx.compile(thomas_mx_loop)
                t0 = time.perf_counter()
                xc = loop_c(*A32T)
                mx.eval(xc)
                print(
                    f"  [N={N} B={B}] mx.compile(Thomas loop) first call: {(time.perf_counter()-t0)*1e3:.0f} ms"
                )

                def run_loop_c():
                    x = loop_c(*A32T)
                    mx.eval(x)
                    return x

                res[methods[3]] = timeit(run_loop_c)

            # PCR eager
            def run_pcr():
                x = pcr_mx(*A32)
                mx.eval(x)
                return x

            res[methods[4]] = timeit(run_pcr)
            acc[methods[4]] = acc3(pcr_mx, A32, A32s, A32s64)

            # PCR compiled
            pcr_c = mx.compile(pcr_mx)
            t0 = time.perf_counter()
            xc = pcr_c(*A32)
            mx.eval(xc)
            print(
                f"  [N={N} B={B}] mx.compile(PCR) first call: {(time.perf_counter()-t0)*1e3:.0f} ms"
            )

            def run_pcr_c():
                x = pcr_c(*A32)
                mx.eval(x)
                return x

            res[methods[5]] = timeit(run_pcr_c)

            # Metal kernel (B,N)
            def run_mk():
                x = thomas_metal_BN(*A32)
                mx.eval(x)
                return x

            res[methods[6]] = timeit(run_mk)
            acc[methods[6]] = acc3(thomas_metal_BN, A32, A32s, A32s64)

            # Metal kernel (N,B)
            def run_mkT():
                x = thomas_metal_NB(*A32T)
                mx.eval(x)
                return x

            res[methods[7]] = timeit(run_mkT)
            acc[methods[7]] = acc3(
                thomas_metal_NB, A32T, A32sT, A32s64T, transpose=True
            )

            solve_rows[(N, B)] = res
            acc_rows[(N, B)] = acc
            print(
                f"done N={N} B={B}: "
                + ", ".join(f"{k}={fmt(v)}" for k, v in res.items())
            )
            sys.stdout.flush()

    # --------------------------------------------------------------- table: solvers
    print("\n## Tridiagonal solve wall time (ms)\n")
    hdr = "| N | B | " + " | ".join(methods) + " |"
    print(hdr)
    print("|" + "---|" * (2 + len(methods)))
    for N in Ns:
        for B in Bs:
            r = solve_rows[(N, B)]
            print(f"| {N} | {B} | " + " | ".join(fmt(r[m]) for m in methods) + " |")

    # --------------------------------------------------------------- table: accuracy
    print("\n## float32 accuracy: max_sys ||x32-x64||_inf / ||x64||_inf\n")
    print(
        "| N | B | method | raw float32 | row-scaled (in f32) | row-scaled (in f64, then cast) |"
    )
    print("|---|---|---|---|---|---|")
    for N in Ns:
        for B in Bs:
            for m, (r0, r1, r2) in acc_rows[(N, B)].items():
                print(f"| {N} | {B} | {m} | {r0:.2e} | {r1:.2e} | {r2:.2e} |")

    # --------------------------------------------------------------- accuracy vs conditioning
    print(
        "\n## float32 accuracy vs conditioning (N=327, B=100; Metal kernel (B,N); numpy f32 Thomas as cross-check)\n"
    )
    print(
        "| d_i regime | cond_2 (system 0, dense) | mx32 Metal raw | mx32 Metal row-scaled | np32 Thomas raw | np32 Thomas row-scaled |"
    )
    print("|---|---|---|---|---|---|")
    for dmode in ("laplacian", "spec", "dominant"):
        a, b, c, d, h = make_systems(100, 327, seed=3, dmode=dmode)
        x64 = reference_solution(a, b, c, d)
        A = np.diag(b[0]) + np.diag(a[0, 1:], -1) + np.diag(c[0, :-1], 1)
        cond = float(np.linalg.cond(A))
        A32 = [mx.array(v.astype(np.float32)) for v in (a, b, c, d)]
        A32s = list(row_scale(*A32))
        e_raw = rel_err(np.asarray(thomas_metal_BN(*A32), dtype=np.float64), x64)
        e_sc = rel_err(np.asarray(thomas_metal_BN(*A32s), dtype=np.float64), x64)
        N32 = [v.T.copy().astype(np.float32) for v in (a, b, c, d)]
        N32s = [v.T.copy().astype(np.float32) for v in row_scale(a, b, c, d)]
        n_raw = rel_err(thomas_np_NB(*N32).T.astype(np.float64), x64)
        n_sc = rel_err(thomas_np_NB(*N32s).T.astype(np.float64), x64)
        print(
            f"| {dmode} | {cond:.1e} | {e_raw:.2e} | {e_sc:.2e} | {n_raw:.2e} | {n_sc:.2e} |"
        )
    sys.stdout.flush()

    # --------------------------------------------------------------- 5. assembly
    print("\n## Newton assembly (F and J, ~10 elementwise ops), wall time (ms)\n")
    print("| N | B | np64 | mx32 eager | mx32 compiled |")
    print("|---|---|---|---|---|")
    rng = np.random.default_rng(1)
    for N in Ns:
        h = make_mesh(N)
        for B in Bs:
            psi = rng.uniform(-1, 1, size=(B, N))
            dop = rng.uniform(-1, 1, size=(B, N))
            s = 1e-3
            t_np = timeit(lambda: assemble_np(psi, h, s, dop))
            psi32 = mx.array(psi.astype(np.float32))
            dop32 = mx.array(dop.astype(np.float32))
            h32 = mx.array(h.astype(np.float32))
            mx.eval(psi32, dop32, h32)

            def run_asm():
                F, J = assemble_mx(psi32, h32, s, dop32)
                mx.eval(F, J)

            t_mx = timeit(run_asm)
            asm_c = mx.compile(lambda p, hh, dd: assemble_mx(p, hh, s, dd))
            F, J = asm_c(psi32, h32, dop32)
            mx.eval(F, J)

            def run_asm_c():
                F, J = asm_c(psi32, h32, dop32)
                mx.eval(F, J)

            t_mxc = timeit(run_asm_c)
            Fn, Jn = assemble_np(psi, h, s, dop)
            F, J = assemble_mx(psi32, h32, s, dop32)
            errF = float(
                np.max(np.abs(np.asarray(F, dtype=np.float64) - Fn))
                / np.max(np.abs(Fn))
            )
            print(
                f"| {N} | {B} | {fmt(t_np)} | {fmt(t_mx)} | {fmt(t_mxc)} |   (F rel err f32 vs f64: {errF:.1e})"
            )
            sys.stdout.flush()

    # --------------------------------------------------------------- 6. full Newton iteration
    print(
        "\n## One full Newton iteration (assembly + tridiagonal solve + update), wall time (ms)\n"
    )
    print(
        "| N | B | np64 (N,B) Thomas | mx32 Metal (N,B) eager | mx32 Metal (N,B) compiled | rel diff psi (f32 vs f64) |"
    )
    print("|---|---|---|---|---|---|")
    for N in [327]:
        h = make_mesh(N)
        inv_h = 1.0 / h
        for B in [b_ for b_ in (1000, 10000) if b_ in Bs or not args.quick]:
            psi = rng.uniform(-0.5, 0.5, size=(N, B))
            dop = rng.uniform(-1, 1, size=(N, B))
            s = 1e-3
            t_np = timeit(lambda: newton_iter_np(psi, h, s, dop, inv_h))
            psi_ref = newton_iter_np(psi, h, s, dop, inv_h)

            psi32 = mx.array(psi.astype(np.float32))
            dop32 = mx.array(dop.astype(np.float32))
            h32 = mx.array(h.astype(np.float32))
            a_fixed = np.zeros((N, B), dtype=np.float32)
            c_fixed = np.zeros((N, B), dtype=np.float32)
            a_fixed[1:-1] = inv_h[:-1, None]
            c_fixed[1:-1] = inv_h[1:, None]
            a32 = mx.array(a_fixed)
            c32 = mx.array(c_fixed)
            mx.eval(psi32, dop32, h32, a32, c32)

            def run_newton():
                out = _newton_iter_mx_core(psi32, h32, s, dop32, a32, c32)
                mx.eval(out)
                return out

            t_mx = timeit(run_newton)
            out = run_newton()
            rel = float(
                np.max(np.abs(np.asarray(out, dtype=np.float64) - psi_ref))
                / np.max(np.abs(psi_ref))
            )

            t_mxc = None
            try:
                core_c = mx.compile(
                    lambda p, hh, dd, aa, cc: _newton_iter_mx_core(p, hh, s, dd, aa, cc)
                )
                o = core_c(psi32, h32, dop32, a32, c32)
                mx.eval(o)

                def run_newton_c():
                    o = core_c(psi32, h32, dop32, a32, c32)
                    mx.eval(o)
                    return o

                t_mxc = timeit(run_newton_c)
            except Exception as e:  # noqa: BLE001
                t_mxc = f"failed: {type(e).__name__}"
            print(
                f"| {N} | {B} | {fmt(t_np)} | {fmt(t_mx)} | {fmt(t_mxc)} | {rel:.1e} |"
            )
            sys.stdout.flush()


if __name__ == "__main__":
    main()
