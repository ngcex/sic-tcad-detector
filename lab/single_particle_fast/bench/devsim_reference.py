"""devsim single-particle references for the fast-engine parity (run ONCE, ~12 min).

    uv run python lab/single_particle_fast/bench/devsim_reference.py OUT_DIR

For each (half-width, LET) of CASES: create_2d_dd_device(hw, 50 V), ion_track_generation_2d,
simulate_single_particle (defaults: 1 ps injection, 10 % BDF1 grid, rtol_converged
1e-7), analyze_current_pulse; the series (times, conduction, displacement, total),
I_dark, Q_generated, Q_collected and t_collection are saved to OUT_DIR/hw{hw}_let{LET}.npz
(lab/single_particle_fast/data holds the copy the parity bench reads).
"""

import sys
import time

import numpy as np

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)

from etna.core.charge_collection_2d import create_2d_dd_device  # noqa: E402
from etna.core.devsim_reset import delete_device_and_mesh  # noqa: E402
from etna.core.single_particle import (  # noqa: E402
    _restore_dark_state,
    _snapshot_dark_state,
    analyze_current_pulse,
    ion_track_generation_2d,
    simulate_single_particle,
)

CASES = {
    2.0: [0.5, 5.0, 15.811388300841896, 500.0],
    50.0: [0.5, 500.0],
    150.0: [0.5],
}

if __name__ == "__main__":
    out = sys.argv[1]
    for hw, lets in CASES.items():
        t0 = time.time()
        info = create_2d_dd_device(half_width_um=hw, V_bias=50.0)
        build_s = time.time() - t0
        dark = _snapshot_dark_state(info)
        for LET in lets:
            t1 = time.time()
            _restore_dark_state(info, dark)
            gen, qg = ion_track_generation_2d(info, LET)
            res = simulate_single_particle(info, gen)
            pulse = analyze_current_pulse(
                res["times"], res["currents"], res["I_dark"],
                currents_total=res["currents_total"],
            )  # fmt: skip
            np.savez(
                f"{out}/hw{hw:g}_let{LET:.6g}.npz",
                Q_gen=qg, times=res["times"], currents=res["currents"],
                disp=res["currents_displacement"], total=res["currents_total"],
                I_dark=res["I_dark"], Q_collected=res["Q_collected"],
                converged=res["converged"], t_collection=pulse["t_collection"],
                build_s=build_s, run_s=time.time() - t1,
            )  # fmt: skip
            print(f"hw {hw:g} LET {LET:.6g}: CCE {res['Q_collected'] / qg:.10f} "
                  f"{time.time() - t1:.1f} s", flush=True)  # fmt: skip
        delete_device_and_mesh(info["device_name"])
