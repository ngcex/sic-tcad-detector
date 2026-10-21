"""Does the guard-ring structure give a lateral CCE map with edge physics? (B-11 check)

devsim guard-ring device (SV half-width 50 um, ring 5 um, gap 3 um, depth 1 um), cathode
ramped to +10 V, anode read-out with its node set, alpha stripe sigma 2 um; then fv2d on
the same imported mesh (parity check).

Usage: uv run python lab/edge_cce/bench/gr_scan.py [legacy]
  legacy: the depth mesh before 2026-10-21 (100 nm cathode cell; CCE +6e-3 too high)
"""

import sys
import time

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

depth_mesh = "legacy" if "legacy" in sys.argv else "graded"
t0 = time.perf_counter()
info = create_guard_ring_device(
    sv_half_width_um=50.0,
    guard_ring_width_um=5.0,
    guard_ring_gap_um=3.0,
    guard_ring_depth_um=1.0,
    depth_mesh=depth_mesh,
)
setup_poisson(info)
solve_equilibrium(info)
setup_sic_drift_diffusion(info)
ramp_bias(info, V_target=10.0, contact="cathode", V_step=0.5, V_step_max=2.5)
ns = guard_ring_node_sets(info)
x_max = info["sv_half_width_cm"] + 1.5e-4  # middle of the 3 um gap
scan = cce_lateral_scan(
    info, n_points=11, contact="anode", node_set=ns["anode"], x_max_cm=x_max
)
t_dev = time.perf_counter() - t0
mesh, _ = mesh_from_devsim(info)
t1 = time.perf_counter()
s = prepared_solver(mesh, {"cathode": 10.0}, dark_check=False)
fv = cce_lateral_scan_fv2d(
    s,
    n_points=11,
    contact="anode",
    node_sets=cut_node_sets(mesh, info["sv_half_width_cm"] + 0.5 * info["gap_cm"]),
    x_max_cm=x_max,
)
t_fv = time.perf_counter() - t1
print(f"depth_mesh {depth_mesh}: {mesh.nx} x {mesh.ny} = {mesh.N} nodes; "
      f"devsim {t_dev:.1f} s (build + ramp + scan), fv2d {t_fv:.1f} s (ramp + scan)")  # fmt: skip
c0 = scan["cce_values"][0]
for x, c, cf in zip(scan["x_positions_um"], scan["cce_values"], fv["cce_values"]):
    print(
        f"RESULT x {x:6.2f} um  CCE_anode {c:.6f}  f {c / c0:.6f}  fv2d-devsim {cf - c:+.1e}"
    )
print("RESULT edge_to_center", scan["edge_to_center_ratio"])
