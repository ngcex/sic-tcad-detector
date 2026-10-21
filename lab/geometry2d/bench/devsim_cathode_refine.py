"""Independent (devsim) check of the contact-cell finding of mesh_convergence.py: the mesa
(ALT_STACK) built by devsim with EXTRA depth lines graded toward the cathode, ramped to
10 V and 50 V, centre CCE by devsim's cce_lateral_scan, and fv2d on the same imported
mesh. If devsim's CCE moves like fv2d's native refinement, the discretisation error of
devsim's default mesh is confirmed by devsim itself.

Usage: uv run python lab/geometry2d/bench/devsim_cathode_refine.py
"""

import devsim
import numpy as np

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)

import etna.core.alternative_structures as A  # noqa: E402
from etna.core.charge_collection_2d import cce_lateral_scan  # noqa: E402
from etna.core.drift_diffusion import ramp_bias, setup_sic_drift_diffusion  # noqa: E402
from etna.core.fv2d import (
    cce_lateral_scan_fv2d,
    mesh_from_devsim,
    prepared_solver,
)  # noqa: E402
from etna.core.poisson import setup_poisson, solve_equilibrium  # noqa: E402

_orig = A._add_depth_mesh_lines


def run(extra, V):
    def patched(mesh_name, junction_pos, total_depth, air_buffer, **kw):
        _orig(mesh_name, junction_pos, total_depth, air_buffer, **kw)
        for d in extra:
            devsim.add_2d_mesh_line(
                mesh=mesh_name, dir="y", pos=total_depth - d, ps=d / 2
            )

    A._add_depth_mesh_lines = patched
    try:
        info = A.create_mesa_device(
            device_name="cathref", pillar_half_width_um=10.0, depth_mesh="legacy"
        )
    finally:
        A._add_depth_mesh_lines = _orig
    try:
        setup_poisson(info)
        solve_equilibrium(info)
        setup_sic_drift_diffusion(info)
        ramp_bias(info, V_target=V, contact="cathode", V_step=0.5, V_step_max=2.5)
        c_dev = cce_lateral_scan(info, n_points=1)["cce_values"][0]
        mesh, _ = mesh_from_devsim(info)
        s = prepared_solver(mesh, {"cathode": V}, dark_check=False)
        c_fv = cce_lateral_scan_fv2d(s, x_positions_cm=[0.0])["cce_values"][0]
        y = np.unique(
            np.array(
                devsim.get_node_model_values(device="cathref", region="sic", name="y")
            )
        )
        return len(y), y[-1] - y[-2], c_dev, c_fv
    finally:
        devsim.delete_device(device="cathref")


if __name__ == "__main__":
    grades = {
        "default": [],
        "to 10 nm": [5e-6, 2e-6, 1e-6],
        "to 1 nm": [5e-6, 2e-6, 1e-6, 5e-7, 2e-7, 1e-7],
    }
    for V in (10.0, 50.0):
        for name, extra in grades.items():
            ny, h, cd, cf = run(extra, V)
            print(f"V {V:4.0f}  {name:9s} ny {ny:4d}  cathode cell {h * 1e7:7.2f} nm  "
                  f"devsim CCE(0) {cd:.8f}  fv2d same mesh {cf:.8f}  diff {cf - cd:+.1e}", flush=True)  # fmt: skip
