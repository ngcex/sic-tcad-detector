"""Default (n+) 2D diode: does the cathode cell matter? fv2d on devsim's own
create_sic_2d_device mesh (== engine2d to < 1e-9 on the legacy mesh, tests/test_fv2d.py),
legacy cathode cell (13.4 nm on the n+ stack, 100 nm without n+) vs a 0.1 nm graded cell, centre
and edge CCE of the alpha stripe (hw 50 um). Also the legacy-stack planar diode (nplus
None) for contrast, and engine2d itself on the production mesh.

Usage: uv run python lab/mesh_contact/bench/default_2d.py
"""

import numpy as np

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.device2d import create_sic_2d_device  # noqa: E402
from etna.core.engine2d import lateral_cce_scan_2d  # noqa: E402
from etna.core.fv2d import (
    cce_lateral_scan_fv2d,
    mesh_from_builder,
    prepared_solver,
)  # noqa: E402

HW = 50.0
VOLTS = (2.0, 5.0, 10.0, 50.0)
LEGACY_STACK = dict(
    substrate_thickness_cm=1e-4,
    nplus_thickness_cm=None,
    buffer_thickness_cm=None,
    buffer_end_doping=None,
    N_D_junction=2.958237e15,
    N_D_bulk=8.252732e13,
    L_transition=0.987e-4,
)


def scan(depth_mesh, **kw):
    m, _ = mesh_from_builder(
        create_sic_2d_device, half_width_um=HW, depth_mesh=depth_mesh, **kw
    )
    out = {}
    for V in VOLTS:
        s = prepared_solver(m, {"cathode": V}, dark_check=False)
        out[V] = cce_lateral_scan_fv2d(s, n_points=3)["cce_values"]
    return m, out


if __name__ == "__main__":
    for label, kw in (
        ("DEFAULT (n+)", {}),
        ("legacy stack (cathode on epi)", LEGACY_STACK),
    ):
        mL, a = scan("legacy", **kw)
        # n+ stack: grade only the last 20 nm (inside the 1 um n+): the cathode cell alone
        # ("graded" there also re-meshes the transition and the epi bottom, a general
        # resolution change of 1e-5 at 2 V, see the 1D uniform-bisection study)
        mG, b = scan({"near_depth_cm": 2e-6} if not kw else "graded", **kw)
        print(f"{label}: ny {mL.ny} -> {mG.ny}, cathode cell "
              f"{(mL.ys[-1] - mL.ys[-2]) * 1e7:.1f} -> {(mG.ys[-1] - mG.ys[-2]) * 1e7:.2f} nm")  # fmt: skip
        for V in VOLTS:
            print(f"   {V:4.0f} V: centre {a[V][0]:.9f} -> {b[V][0]:.9f} "
                  f"(d {b[V][0] - a[V][0]:+.2e}); edge d {b[V][-1] - a[V][-1]:+.2e}")  # fmt: skip
        if not kw:
            for V in VOLTS[1:]:
                e2 = lateral_cce_scan_2d(HW, V, n_lateral_points=3)["cce_values"]
                print(f"   engine2d {V:4.0f} V: {np.array2string(e2, precision=9)} "
                      f"(fv2d legacy mesh - engine2d {np.max(np.abs(a[V] - e2)):.1e})")  # fmt: skip
