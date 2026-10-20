"""Does the guard-ring structure give a lateral CCE map with edge physics? (B-11 check)"""

import numpy as np

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.alternative_structures import (
    create_guard_ring_device,
    guard_ring_node_sets,
)
from etna.core.charge_collection_2d import cce_lateral_scan
from etna.core.drift_diffusion import ramp_bias, setup_sic_drift_diffusion
from etna.core.poisson import setup_poisson, solve_equilibrium

info = create_guard_ring_device(
    sv_half_width_um=50.0,
    guard_ring_width_um=5.0,
    guard_ring_gap_um=3.0,
    guard_ring_depth_um=1.0,
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
for x, c in zip(scan["x_positions_um"], scan["cce_values"]):
    print(f"RESULT x {x:6.2f} um  CCE_anode {c:.6f}  f {c / scan['cce_values'][0]:.6f}")
print("RESULT edge_to_center", scan["edge_to_center_ratio"])
