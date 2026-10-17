"""Serial warm timing, 2D engine vs devsim WITH n+ (lateral CCE scan, 10 points).

usage: uv run python lab/solver2d/bench/nplus2d_bench.py
"""

import os
import time
import uuid

import devsim
import numpy as np

from etna.core.charge_collection_2d import cce_lateral_scan, create_2d_dd_device
from etna.core.engine2d import lateral_cce_scan_2d

KW = dict(
    epi_thickness_cm=10e-4, N_D_bulk=8.5e13, nplus_thickness_cm=1e-4, N_D_nplus=1e19
)
V = 60.0

lateral_cce_scan_2d(25.0, V, n_lateral_points=2, **KW)  # JIT warm-up
print("load avg", os.getloadavg())
for hw in (50.0, 150.0):
    te = []
    for _ in range(3):
        t = time.perf_counter()
        e = lateral_cce_scan_2d(hw, V, n_lateral_points=10, **KW)
        te.append(time.perf_counter() - t)
    t = time.perf_counter()
    name = f"bench_{uuid.uuid4().hex[:8]}"
    info = create_2d_dd_device(half_width_um=hw, V_bias=V, device_name=name, **KW)
    t_build = time.perf_counter() - t
    d = cce_lateral_scan(info, n_points=10)
    td = time.perf_counter() - t
    devsim.delete_device(device=name)
    dC = np.max(np.abs(e["cce_values"] - d["cce_values"]))
    print(
        f"hw={hw:g}: devsim build+ramp {t_build:.2f} s, total {td:.2f} s; engine "
        f"{min(te):.3f} s (runs {', '.join(f'{x:.3f}' for x in te)}); "
        f"speedup {td / min(te):.0f}x; max dCCE {dC:.1e}",
        flush=True,
    )
print("load avg", os.getloadavg())
