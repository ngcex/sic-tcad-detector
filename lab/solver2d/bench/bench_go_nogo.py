"""GO/NO-GO: single-process wall time of (equilibrium + ramp to 50 V + 10-point lateral
CCE scan), devsim vs the 2D engine, back to back in one process, at hw = 50 and 150.

usage: uv run python lab/solver2d/bench/bench_go_nogo.py [hw ...]
Engine runs single-threaded (no prange). JIT is warmed on hw=25 first (compile time
reported separately). Engine time includes building devsim's mesh (no equations).
"""

import os
import sys
import time
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import devsim  # noqa: E402
import numpy as np  # noqa: E402

from etna.core.charge_collection_2d import (
    cce_lateral_scan,
    create_2d_dd_device,
)  # noqa: E402
from solver2d.engine import Device2D  # noqa: E402
from solver2d.grid import grid_from_devsim  # noqa: E402
from solver2d.kernels import cathode_current  # noqa: E402
from etna.core.drift_diffusion import extract_contact_current  # noqa: E402

V = 50.0
NPTS = 10


def engine(hw, full2d=False):
    t0 = time.perf_counter()
    g = grid_from_devsim.__wrapped__(hw)
    t1 = time.perf_counter()
    dev = Device2D(g)
    dev.dark(V, full2d=full2d)
    Ie, Ih = cathode_current(dev.u, dev.n, dev.p, dev.kn, dev.kp, dev.g2.hy, dev.g2.dxd, dev.g2.nx, dev.g2.ny)
    t2 = time.perf_counter()
    dev.factor_continuity()
    t3 = time.perf_counter()
    _, cce = dev.cce_lateral_scan(NPTS)
    t4 = time.perf_counter()
    return cce, dict(
        mesh=t1 - t0,
        dark=t2 - t1,
        factor=t3 - t2,
        scan=t4 - t3,
        total=t4 - t0,
        gummel=dev.stats.gummel,
        newton_2d=dev.stats.newton_2d,
        I_dark=1.602e-19 * (Ie + Ih),
    )


def ref(hw):
    t0 = time.perf_counter()
    info = create_2d_dd_device(half_width_um=hw, V_bias=V)
    t1 = time.perf_counter()
    I_dark = extract_contact_current(info, "cathode")
    res = cce_lateral_scan(info, n_points=NPTS)
    t2 = time.perf_counter()
    devsim.delete_device(device=info["device_name"])
    return res["cce_values"], dict(build=t1 - t0, scan=t2 - t1, total=t2 - t0, I_dark=I_dark)


tc = time.perf_counter()
engine(25.0)
print(
    f"JIT warm-up (compile from cache or fresh) + hw=25 run: {time.perf_counter() - tc:.2f}s"
)
saved = {}
for hw in [float(a) for a in sys.argv[1:]] or [50.0, 150.0]:
    print(f"\n=== hw={hw:g} um, load {os.getloadavg()[0]:.1f}")
    ce, te = min((engine(hw) for _ in range(3)), key=lambda r: r[1]["total"])
    cd, td = ref(hw)
    ce2, te2 = engine(hw)
    ce3, te3 = engine(hw, full2d=True)
    saved[f"hw{hw:g}"] = dict(cce_engine=ce, cce_devsim=np.asarray(cd), cce_engine_full2d=ce3)
    print("engine full-2D ramp:", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in te3.items()}, f"speedup {td['total'] / te3['total']:.1f}x, max|dCCE| {np.abs(ce3 - cd).max():.2e}")
    print(
        "engine (best of 3):",
        {k: (round(v, 4) if isinstance(v, float) else v) for k, v in te.items()},
    )
    print("engine (after devsim):", round(te2["total"], 4))
    print("devsim:", {k: round(v, 3) for k, v in td.items()})
    print("cce engine", np.array2string(ce, precision=8))
    print("cce devsim", np.array2string(np.asarray(cd), precision=8))
    print(
        f"max |dCCE| = {np.abs(ce - cd).max():.2e}   speedup = {td['total'] / te['total']:.1f}x (vs post-devsim engine run {td['total'] / te2['total']:.1f}x)"
    )

np.savez("lab/solver2d/notes/cce_ref_lateral_scan.npz", **{f"{h}_{k}": v for h, d in saved.items() for k, v in d.items()})
