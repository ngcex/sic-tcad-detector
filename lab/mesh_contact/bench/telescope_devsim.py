"""Delta-E/E telescope (notebook-19 settings: 2 um delta-E stage at 10 V, 50 um E stage
at 50 V) on the legacy and the graded depth mesh, devsim (the two-region device has no
fv2d twin). Centre CCE of each stage, read at its cathode.

Both cathodes are ohmic contacts on n- epi. The E stage is far from full depletion at
50 V (V_fd ~ 180 V for 49 um of 8e13) and the alpha (range ~18 um) stops long before its
cathode, so only the delta-E stage (1 um epi, fully depleted, the alpha goes through)
can see the cathode cell.

Usage: uv run python lab/mesh_contact/bench/telescope_devsim.py [legacy|graded|fine]
  fine: a 10x finer contact cell, 1.1x growth and 25 nm near-cathode cells (check)
"""

import sys
import time

import devsim

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.alternative_structures import create_delta_e_e_device  # noqa: E402
from etna.core.charge_collection_2d import cce_lateral_scan  # noqa: E402
from etna.core.drift_diffusion import ramp_bias, setup_sic_drift_diffusion  # noqa: E402
from etna.core.poisson import setup_poisson, solve_equilibrium  # noqa: E402

FINE = {"first_cell_cm": 1e-9, "ratio": 1.1, "near_ps_cm": 2.5e-6}
meshes = [a for a in sys.argv[1:] if a in ("legacy", "graded", "fine")] or [
    "legacy",
    "graded",
]
for dm in meshes:
    spec = FINE if dm == "fine" else dm
    t = time.perf_counter()
    info = create_delta_e_e_device(
        device_name=f"tel_{dm}",
        half_width_um=50.0,
        delta_e_thickness_um=2.0,
        e_stop_thickness_um=50.0,
        depth_mesh=spec,
    )
    L = info["layers"]
    try:
        for lay in L.values():
            setup_poisson(lay)
        solve_equilibrium(info)
        for lay in L.values():
            setup_sic_drift_diffusion(lay)
        ramp_bias(L["delta_e"], V_target=10.0, contact="de_cathode", V_step=0.5,
                  V_step_max=2.5)  # fmt: skip
        ramp_bias(L["e_stop"], V_target=50.0, contact="e_cathode", V_step=0.5,
                  V_step_max=2.5)  # fmt: skip
        de = cce_lateral_scan(L["delta_e"], n_points=1, contact="de_cathode")
        e = cce_lateral_scan(L["e_stop"], n_points=1, contact="e_cathode")
        print(f"RESULT {dm}: nodes {info['num_nodes']}, delta-E CCE(0) "
              f"{de['cce_values'][0]:.6f}, E CCE(0) {e['cce_values'][0]:.6f} "
              f"({time.perf_counter() - t:.0f} s)", flush=True)  # fmt: skip
    finally:
        devsim.delete_device(device=info["device_name"])
