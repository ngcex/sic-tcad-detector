"""Reference CCE(LET) points with the build-per-LET loop (pre-optimisation code).

Runs the original build_cce_let_table inner loop (fresh 2D device per LET) and saves
the per-step transient series so later variants (build once, linear scaling,
looser stop) can be compared point by point.

    uv run python lab/engine_box/bench/let_table_reference.py HW OUT.npz LET [LET ...]
"""

import sys
import time

import numpy as np

import etna
from etna.core import single_particle as sp
from etna.core.charge_collection_2d import create_2d_dd_device
from etna.core.devsim_reset import delete_device_and_mesh

etna.set_devsim_quiet(True)

hw = float(sys.argv[1])
out = sys.argv[2]
lets = [float(v) for v in sys.argv[3:]]
save = {}
for L in lets:
    t0 = time.perf_counter()
    info = create_2d_dd_device(half_width_um=hw, V_bias=50.0)
    tb = time.perf_counter() - t0
    gen, q_gen = sp.ion_track_generation_2d(info, L)
    r = sp.simulate_single_particle(info, gen)
    tt = time.perf_counter() - t0 - tb
    delete_device_and_mesh(info["device_name"])
    cce = r["Q_collected"] / q_gen
    print(
        f"LET {L:g}: CCE_raw={cce:.10f} steps={len(r['times'])} "
        f"build={tb:.1f}s transient={tt:.1f}s",
        flush=True,
    )
    key = f"{L:g}"
    save[f"times_{key}"] = r["times"]
    save[f"currents_{key}"] = r["currents"]
    save[f"idark_{key}"] = r["I_dark"]
    save[f"qgen_{key}"] = q_gen
    save[f"cce_{key}"] = cce
    np.savez(out, **save)
