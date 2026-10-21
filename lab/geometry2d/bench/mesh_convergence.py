"""Uniform-bisection mesh convergence of the centre CCE (fv2d, native mesa = planar
half-device, ALT_STACK) starting from devsim's own mesh lines of create_mesa_device,
at 10 V (partially depleted) and 50 V. Separates the discretisation error of devsim's
mesh from the solver: level 0 IS devsim's mesh (parity 1e-16), each level bisects every
depth interval (lateral lines kept: the problem is 1D at the centre).

Usage: uv run python lab/geometry2d/bench/mesh_convergence.py
"""

import numpy as np

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.alternative_structures import create_mesa_device  # noqa: E402
from etna.core.fv2d import (
    cce_lateral_scan_fv2d,
    mesh_from_builder,
    prepared_solver,
)  # noqa: E402
from etna.core.fv2d.refine import refine_lines  # noqa: E402
from etna.core.fv2d.structures import build_mesh, mesa_structure  # noqa: E402

dm, _ = mesh_from_builder(create_mesa_device)
st = mesa_structure()
for V in (10.0, 50.0):
    ys = dm.ys
    prev = None
    for lev in range(5):
        mesh = build_mesh(st, dm.xs, ys)
        s = prepared_solver(mesh, {"cathode": V}, dark_check=False)
        c = cce_lateral_scan_fv2d(s, x_positions_cm=[0.0])["cce_values"][0]
        d = "" if prev is None else f"  dCCE {c - prev:+.3e}"
        print(f"V {V:4.0f}  level {lev}: ny {len(ys):5d}  CCE(0) {c:.8f}{d}")
        prev = c
        ys = refine_lines(ys, np.ones(len(ys) - 1, bool))
