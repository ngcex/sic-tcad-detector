"""Linear-solver choice of fv2d: numba banded LU vs SuperLU (MMD on A^T + A), for the dark
Newton ramp and the illuminated Gummel chord solves, at the devsim mesh size (12k nodes)
and on refined native meshes (up to ~50k nodes). Same session, one timed run per entry (the JIT is warmed by the first mesh).

Usage: uv run python lab/geometry2d/bench/linear_solvers.py
"""

import time

import numpy as np

from etna.core.fv2d import Solver2D, alpha_stripe_generation_fv2d
from etna.core.fv2d.refine import refine_lines
from etna.core.fv2d.structures import build_mesh, electrode3d_structure, initial_lines


def best(f, k=1):
    t = []
    for _ in range(k):
        t0 = time.perf_counter()
        r = f()
        t.append(time.perf_counter() - t0)
    return min(t), r


st = electrode3d_structure()
xs, ys = initial_lines(st)
for level in range(4):
    mesh = build_mesh(st, xs, ys)
    Gs = np.array(
        [alpha_stripe_generation_fv2d(mesh, r) for r in np.linspace(8e-4, 4.3e-3, 6)]
    )
    row = [f"{mesh.nx}x{mesh.ny}={mesh.N}"]
    ref = None
    for lin in ("superlu", "band") if mesh.N < 40000 else ("superlu",):

        def ramp():
            s = Solver2D(mesh)
            s.linear_solver = lin
            s.equilibrium().set_bias({"cathode": 50.0})
            return s

        t, s = best(ramp)
        row.append(f"ramp[{lin}] {t:.2f}s ({s.newton_iterations} Newton)")
    for g in ("band", "superlu") if mesh.N < 40000 else ("superlu",):
        s.gummel_solver = g
        s._lu_cont = None
        t, (states, Gf) = best(lambda: s.solve_generation(Gs))
        cce = np.array(
            [
                s.contact_currents(x, G)["cathode"] / s.injected_generation(G)
                for x, G in zip(states, Gf)
            ]
        )
        if ref is None:
            ref = cce
        row.append(f"scan6[{g}] {t:.2f}s (dCCE {np.max(np.abs(cce - ref)):.1e})")
    print(" | ".join(row), flush=True)
    if level % 2 == 0:
        ys = refine_lines(ys, np.ones(len(ys) - 1, bool))
    else:
        xs = refine_lines(xs, np.ones(len(xs) - 1, bool))
