"""Is a true 3D DD solver tractable with ETNA's direct-solver approach? Measures SuperLU
(MMD on A^T + A, as fv2d) on 7-point 3D Laplacians n^3 (an optimistic stand-in for the
Poisson / SG Jacobians: same sparsity, better conditioned), and Jacobi-preconditioned CG
for comparison, then extrapolates to a 3D unit cell meshed like the 2D structures
(about 114 x 114 x 69 lines = 0.9 M nodes, the native 3D-electrode lines in x, y and z).

Usage: uv run python lab/geometry2d/bench/tractability_3d.py
"""

import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def lap(n, dim):
    e = np.ones(n)
    T = sp.diags([-e[:-1], 2 * e, -e[:-1]], [-1, 0, 1])
    I = sp.identity(n)
    if dim == 2:
        return (sp.kron(T, I) + sp.kron(I, T)).tocsc()
    return (
        sp.kron(sp.kron(T, I), I)
        + sp.kron(sp.kron(I, T), I)
        + sp.kron(sp.kron(I, I), T)
    ).tocsc()


rows = []
for dim, ns in ((2, (100, 200, 300)), (3, (20, 30, 40, 50))):
    for n in ns:
        A = lap(n, dim)
        N = A.shape[0]
        b = np.ones(N)
        t0 = time.perf_counter()
        lu = spla.splu(A, permc_spec="MMD_AT_PLUS_A")
        tf = time.perf_counter() - t0
        t0 = time.perf_counter()
        lu.solve(b)
        ts = time.perf_counter() - t0
        nnz = lu.L.nnz + lu.U.nnz
        it = [0]
        t0 = time.perf_counter()
        Minv = sp.diags(1.0 / A.diagonal())
        spla.cg(
            A,
            b,
            M=Minv,
            rtol=1e-10,
            maxiter=20000,
            callback=lambda xk: it.__setitem__(0, it[0] + 1),
        )
        tc = time.perf_counter() - t0
        rows.append((dim, N, tf, ts, nnz))
        print(f"{dim}D N={N:8d}: splu {tf:7.3f} s, solve {ts:6.3f} s, nnz(L+U) {nnz:10d} "
              f"({nnz * 12 / 1e6:7.1f} MB), Jacobi-CG {it[0]:5d} it {tc:6.2f} s", flush=True)  # fmt: skip

d3 = np.array([(N, tf, nnz) for dim, N, tf, ts, nnz in rows if dim == 3], float)
pt = np.polyfit(np.log(d3[:, 0]), np.log(d3[:, 1]), 1)
pm = np.polyfit(np.log(d3[:, 0]), np.log(d3[:, 2]), 1)
Nt = 114 * 114 * 69
print(f"3D fit: factor time ~ N^{pt[0]:.2f}, fill ~ N^{pm[0]:.2f}")
print(
    f"extrapolated to N = {Nt}: factor {np.exp(np.polyval(pt, np.log(Nt))):.0f} s, "
    f"fill {np.exp(np.polyval(pm, np.log(Nt))) * 12 / 1e9:.1f} GB per factorisation "
    f"(x3 for the Gummel chord set, x tens of Newton steps for a bias ramp)"
)
