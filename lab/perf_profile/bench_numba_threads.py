"""Single-device engine calls: NUMBA_NUM_THREADS=1 vs the default (all cores).

The batched kernels prange over devices; a facade call on one config has B = 1, so
the parallel launch is pure overhead. Fresh interpreter per setting, interleaved,
min of N in-process repeats after a warm-up.

usage: uv run python lab/perf_profile/bench_numba_threads.py [rounds]
"""

from __future__ import annotations

import os
import subprocess
import sys

CHILD = r"""
import time
from etna import DeviceConfig, microdosimeter_config, run_cce, run_cv, run_cce_tilted_alpha, run_dark_current, run_cv_batch
cases = {
  "run_cv": lambda: run_cv(DeviceConfig()),
  "run_cce": lambda: run_cce(DeviceConfig()),
  "run_dark_current": lambda: run_dark_current(DeviceConfig()),
  "run_cce_tilted_alpha": lambda: run_cce_tilted_alpha(microdosimeter_config()),
  "run_cv_batch x32": lambda: run_cv_batch([DeviceConfig(N_D_bulk=3e13 + 1e12 * i) for i in range(32)], 0.0, -50.0, 40),
}
for name, fn in cases.items():
    fn()
    best = 1e9
    for _ in range(5):
        t = time.perf_counter(); fn(); best = min(best, time.perf_counter() - t)
    print(f"RES {name} {best:.4f}")
"""


def run(threads):
    env = dict(os.environ)
    if threads:
        env["NUMBA_NUM_THREADS"] = str(threads)
    else:
        env.pop("NUMBA_NUM_THREADS", None)
    out = subprocess.run(
        [sys.executable, "-c", CHILD], capture_output=True, text=True, env=env
    )
    return {
        ln.split()[1]: float(ln.split()[-1])
        for ln in out.stdout.splitlines()
        if ln.startswith("RES")
    } | ({"_err": out.stderr[-500:]} if out.returncode else {})


if __name__ == "__main__":
    rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    agg = {"default": {}, "1": {}}
    for _ in range(rounds):
        for key, th in (("default", None), ("1", 1)):
            for k, v in run(th).items():
                if k != "_err":
                    agg[key][k] = min(agg[key].get(k, 1e9), v)
                else:
                    print(v)
    for k in agg["default"]:
        d, s = agg["default"][k], agg["1"].get(k, float("nan"))
        print(f"{k:24s} default {d:.4f} s | 1 thread {s:.4f} s | x{d / s:.2f}")
