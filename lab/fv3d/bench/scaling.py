"""Time / memory scaling of the fv3d MG-PCG solver vs the direct solver (SuperLU, as fv2d),
and one full nonlinear solve of the issue #4 named point (a ~0.9 M-node 3D column cell).

Part 1: n^3 grids, operator = the fv3d box Laplacian with a Dirichlet column and a
Dirichlet back plane (c = 0: the hardest case for the preconditioner), rtol 1e-10.
Memory = tracemalloc peak of the numpy allocations of setup + solve (numba kernels
allocate nothing). SuperLU (MMD on A^T + A) on the same matrices up to n = 40.
Part 2: column_cell (30 um columns + p+ back plane, 50 um) refined to ~0.9 M nodes:
equilibrium + 50 V + 100 V, wall time and peak memory.

Usage: uv run python lab/fv3d/bench/scaling.py [--big]
"""

import sys
import time
import tracemalloc

import numpy as np
import scipy.sparse.linalg as spla

from etna.core.fv3d import Poisson3D, column_cell
from etna.core.fv3d.mg import Multigrid, pcg


def cube(n):
    xs = np.linspace(0, 1, n)
    Z, Y, X = np.meshgrid(xs, xs, xs, indexing="ij")
    free = ~((np.hypot(X, Y) <= 0.2) & (Z <= 0.6) | (Z >= 1.0 - 1e-12))
    return xs, free.astype(np.uint8)


def mg_point(n):
    xs, free = cube(n)
    rng = np.random.default_rng(0)
    tracemalloc.start()
    t0 = time.perf_counter()
    mg = Multigrid(xs, xs, xs, free)
    mg.set_c(np.zeros(free.shape))
    t1 = time.perf_counter()
    b = rng.standard_normal(free.shape) * (free != 0)
    x, it, rr = pcg(mg, b, rtol=1e-10, maxiter=300)
    t2 = time.perf_counter()
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    return mg, b, x, t1 - t0, t2 - t1, it, rr, peak


def part1():
    mg_point(16)  # JIT warm-up
    rows = []
    print(
        "n^3 grid       N   MG levels  setup s  PCG s  it  rel.res   peak MB  B/node | SuperLU factor s  fill MB  max|x_mg - x_lu|/max|x|"
    )
    for n in (32, 40, 64, 96, 128, 160, 200):
        mg, b, x, ts, tsol, it, rr, peak = mg_point(n)
        N = n**3
        line = (f"{n:4d}^3 {N:9d} {len(mg.levels):5d} {ts:9.2f} {tsol:7.2f} {it:3d} {rr:9.1e} "
                f"{peak / 1e6:8.1f} {peak / N:6.0f}")  # fmt: skip
        if n <= 40:
            A = mg.fine.sparse().tocsc()
            t0 = time.perf_counter()
            lu = spla.splu(A, permc_spec="MMD_AT_PLUS_A")
            tf = time.perf_counter() - t0
            xl = lu.solve(b.ravel()).reshape(b.shape)
            fill = (lu.L.nnz + lu.U.nnz) * 12 / 1e6
            line += f" | {tf:8.2f} {fill:9.1f} {np.abs(x - xl).max() / np.abs(xl).max():.1e}"
            rows.append((N, tsol, peak, tf, fill))
        else:
            rows.append((N, tsol, peak, np.nan, np.nan))
        print(line, flush=True)
    a = np.array(rows)
    big = a[:, 0] >= 2.5e5  # small grids are overhead-dominated
    pt = np.polyfit(np.log(a[big, 0]), np.log(a[big, 1]), 1)
    pm = np.polyfit(np.log(a[big, 0]), np.log(a[big, 2]), 1)
    print(f"MG-PCG fit (N >= 2.5e5): time ~ N^{pt[0]:.2f}, memory ~ N^{pm[0]:.2f}")
    Nt = 114 * 114 * 69
    for name, NN in (("issue #4 cell 114x114x69", Nt), ("8 GB budget", None)):
        if NN is None:
            # largest N whose extrapolated peak fits 8 GB
            NN = float(np.exp((np.log(8e9) - pm[1]) / pm[0]))
        print(f"  {name}: N = {NN:.3g}: PCG solve {np.exp(np.polyval(pt, np.log(NN))):.1f} s, "
              f"peak {np.exp(np.polyval(pm, np.log(NN))) / 1e9:.2f} GB")  # fmt: skip


def part2():
    st = column_cell(depth_n=30e-4, depth_p=30e-4, back="p+", h_fine=0.15e-4,
                     h_bulk=0.8e-4)  # fmt: skip
    print(f"column cell {st.shape} = {st.N} nodes")
    tracemalloc.start()
    t0 = time.perf_counter()
    p = Poisson3D(st)
    p.equilibrium()
    t1 = time.perf_counter()
    for V in (50.0, 100.0):
        s = p.solve({"cathode": V})
        print(f"  {V:5.0f} V: Newton {s.newton_iterations}, CG {s.cg_iterations}, {s.wall_s:.1f} s, "
              f"depleted {p.depleted_fraction(s):.5f}", flush=True)  # fmt: skip
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    print(f"  setup + equilibrium {t1 - t0:.1f} s, total {time.perf_counter() - t0:.1f} s, "
          f"peak {peak / 1e9:.2f} GB ({peak / st.N:.0f} B/node)")  # fmt: skip


if __name__ == "__main__":
    part1()
    if "--big" in sys.argv:
        part2()
