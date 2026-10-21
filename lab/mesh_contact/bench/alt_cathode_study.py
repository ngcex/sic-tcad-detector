"""Convergence study of the depth mesh of the contact-on-epi builders (ALT_STACK mesa =
planar half-device, centre CCE at 10 V (partial depletion) and 50 V).

Meshes are devsim's own (create_mesa_device(depth_mesh=spec), mesh-only build) solved by
fv2d (parity with devsim on the same mesh: 1e-16, lab/geometry2d), on a 3-column strip
at the centre (the centre problem is 1D: the strip differs from the full mesa by
< 2e-7). The REFERENCE is uniform bisection of every depth interval, 7-8 levels, from
two different starting meshes: the increments become asymptotically first order (they
halve), so the limit is L_k + d_k (Richardson, p = 1) with the error estimate
|d_k - d_{k-1} / 2|; the two starts must agree within it.

Usage: uv run python lab/mesh_contact/bench/alt_cathode_study.py [quick]
"""

import sys
import time

import numpy as np

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.alternative_structures import (  # noqa: E402
    create_guard_ring_device,
    create_mesa_device,
)
from etna.core.fv2d import (  # noqa: E402
    cce_lateral_scan_fv2d,
    mesh_from_builder,
    prepared_solver,
)
from etna.core.fv2d.refine import refine_lines  # noqa: E402
from etna.core.fv2d.structures import build_mesh, mesa_structure  # noqa: E402

VOLTS = (10.0, 50.0)
ST = mesa_structure()


def centre_cce(xs, ys, V):
    mesh = build_mesh(ST, xs, ys)
    s = prepared_solver(mesh, {"cathode": V}, dark_check=False)
    return cce_lateral_scan_fv2d(s, x_positions_cm=[0.0])["cce_values"][0]


def lines(spec):
    dm, _ = mesh_from_builder(create_mesa_device, depth_mesh=spec)
    return dm


def reference(spec, V, levels):
    """Uniform-bisection sequence; returns (values, limit, error estimate)."""
    dm = lines(spec)
    xs, y = dm.xs[:3], dm.ys
    vals = []
    for _ in range(levels):
        vals.append(centre_cce(xs, y, V))
        y = refine_lines(y, np.ones(len(y) - 1, bool))
    d = np.diff(vals)
    return vals, vals[-1] + d[-1], abs(d[-1] - 0.5 * d[-2])


SPECS = {
    "legacy (100 nm / 10 nm)": "legacy",
    "cathode 1 nm, junction 10 nm": {"first_cell_cm": 1e-7, "junction_ps_cm": 1e-6},
    "cathode 1 nm, junction 1 nm": {"first_cell_cm": 1e-7},
    "graded default (0.1 nm / 1 nm)": "graded",
    "cathode 0.01 nm, junction 1 nm": {"first_cell_cm": 1e-9},
    "default, near 100 nm": {"near_ps_cm": 1e-5},
    "default, near 25 nm": {"near_ps_cm": 2.5e-6},
    "default, ratio 1.25": {"ratio": 1.25},
    "default, near depth 1 um": {"near_depth_cm": 1e-4},
}

if __name__ == "__main__":
    quick = "quick" in sys.argv
    ref = {}
    for V in VOLTS:
        lims = []
        for start in ({"first_cell_cm": 1e-7, "junction_ps_cm": 1e-6}, "graded"):
            t = time.perf_counter()
            vals, lim, err = reference(start, V, 6 if quick else 8)
            lims.append((lim, err))
            print(f"REF {V:4.0f} V start {start}: last {vals[-1]:.9f} increments "
                  f"{np.array2string(np.diff(vals)[-3:], precision=3)} -> limit {lim:.9f}"
                  f" +- {err:.1e} ({time.perf_counter() - t:.1f}s)", flush=True)  # fmt: skip
        ref[V] = lims[-1][0]
        print(f"REF {V:4.0f} V: limits agree to {abs(lims[0][0] - lims[1][0]):.1e}")
    print()
    for name, spec in SPECS.items():
        dm = lines(spec)
        cell = (dm.ys[-1] - dm.ys[-2]) * 1e7
        msg = f"{name:32s} ny {len(dm.ys):4d} cathode cell {cell:7.3f} nm"
        for V in VOLTS:
            c = centre_cce(dm.xs[:3], dm.ys, V)
            msg += f" | {V:2.0f} V {c:.8f} err {c - ref[V]:+.2e}"
        print(msg, flush=True)
    for spec in ("legacy", "graded"):
        gr, _ = mesh_from_builder(create_guard_ring_device, depth_mesh=spec)
        me = lines(spec)
        print(f"nodes {spec:7s}: mesa {me.N} ({me.nx} x {me.ny}), guard ring {gr.N} "
              f"({gr.nx} x {gr.ny})")  # fmt: skip
