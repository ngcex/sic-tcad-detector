"""devsim convergence of the guard ring on the graded depth mesh up to 50 V (proxy with a
10 um SV half-width to stay inside the time budget; the full 50 um device is the
notebook-19 run). Reports the ramp tiers and the centre anode CCE vs fv2d on the same
mesh.

Usage: uv run python lab/mesh_contact/bench/gr_devsim_50v.py [sv_half_width_um]
"""

import sys
import time

import devsim

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.alternative_structures import (  # noqa: E402
    create_guard_ring_device,
    guard_ring_node_sets,
)
from etna.core.charge_collection_2d import cce_lateral_scan  # noqa: E402
from etna.core.drift_diffusion import ramp_bias, setup_sic_drift_diffusion  # noqa: E402
from etna.core.fv2d import (  # noqa: E402
    cce_lateral_scan_fv2d,
    cut_node_sets,
    mesh_from_devsim,
    prepared_solver,
)
from etna.core.poisson import setup_poisson, solve_equilibrium  # noqa: E402

sv = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
t = time.perf_counter()
info = create_guard_ring_device(
    device_name="gr50", sv_half_width_um=sv, guard_ring_width_um=5.0,
    guard_ring_gap_um=3.0, guard_ring_depth_um=1.0,
)  # fmt: skip
try:
    setup_poisson(info)
    solve_equilibrium(info)
    setup_sic_drift_diffusion(info)
    ramp_bias(info, V_target=50.0, contact="cathode", V_step=0.5, V_step_max=2.5)
    ns = guard_ring_node_sets(info)
    sc = cce_lateral_scan(info, n_points=1, contact="anode", node_set=ns["anode"])
    mesh, _ = mesh_from_devsim(info)
    s = prepared_solver(mesh, {"cathode": 50.0}, dark_check=False)
    fv = cce_lateral_scan_fv2d(
        s, x_positions_cm=[0.0], contact="anode",
        node_sets=cut_node_sets(mesh, info["sv_half_width_cm"] + 0.5 * info["gap_cm"]),
    )  # fmt: skip
    print(f"RESULT sv {sv:g} um, {mesh.N} nodes: devsim reached 50 V in "
          f"{time.perf_counter() - t:.0f} s (relaxed {info.get('ramp_relaxed_retries', 0)}, "
          f"aggressive {info.get('ramp_aggressive_retries', 0)}); CCE(0) devsim "
          f"{sc['cce_values'][0]:.7f}, fv2d {fv['cce_values'][0]:.7f}")  # fmt: skip
finally:
    devsim.delete_device(device="gr50")
