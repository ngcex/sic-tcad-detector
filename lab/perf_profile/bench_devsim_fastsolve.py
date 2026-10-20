"""A/B: devsim stock umfshim callback vs etna.core.devsim_fastsolve (cached init).

Interleaved in one process, min of N; parity = bit-identical node arrays.
usage: uv run python lab/perf_profile/bench_devsim_fastsolve.py [repeats]
"""

from __future__ import annotations

import sys
import time

import numpy as np

import devsim
from etna import DeviceConfig, run_dark_current, run_field, set_devsim_quiet
from etna.core import devsim_fastsolve as fs

set_devsim_quiet(True)
reps = int(sys.argv[1]) if len(sys.argv) > 1 else 5
from devsim.umfpack import umfshim  # noqa: E402

print(
    "callback identity preserved by get_parameter:",
    devsim.get_parameter(name="solver_callback")
    in (umfshim.local_solver_callback, fs.solver_callback),
)

CASES = {
    "run_field 1D (18 solves)": lambda: run_field(DeviceConfig(), -30.0).metadata[
        "potential"
    ],
    "run_dark_current devsim (47 solves)": lambda: run_dark_current(
        DeviceConfig(), 0.0, -60.0, 7, solver="devsim"
    ).y,
}
for name, fn in CASES.items():
    res = {"stock": [], "fast": []}
    out = {}
    for _ in range(reps):
        for mode in ("stock", "fast"):
            fs.uninstall() if mode == "stock" else fs.install()
            t, c = time.perf_counter(), time.process_time()
            out[mode] = np.asarray(fn())
            res[mode].append((time.perf_counter() - t, time.process_time() - c))
    fs.install()
    ws = min(r[0] for r in res["stock"])
    wf = min(r[0] for r in res["fast"])
    cs = min(r[1] for r in res["stock"])
    cf = min(r[1] for r in res["fast"])
    print(
        f"{name}: wall {ws:.3f} -> {wf:.3f} s (x{ws / wf:.3f}), cpu {cs:.3f} -> {cf:.3f} s "
        f"(x{cs / cf:.3f}); bit-identical {np.array_equal(out['stock'], out['fast'])}"
    )
