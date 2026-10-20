"""Share of devsim.solve spent in devsim's Python UMFPACK shim (ctypes), per path.

devsim (macOS wheel) solves its linear systems through a Python solver_callback
(devsim/umfpack/umfshim.py): every solve 'init' builds a dsobject that dlopens
libumfpack_lgpl + the BLAS again; factor/solve go through ctypes.

usage: uv run python lab/perf_profile/devsim_shim_share.py [field1d|transient|dark]
"""

from __future__ import annotations

import cProfile
import pstats
import sys
import time

from etna import (
    DeviceConfig,
    run_dark_current,
    run_field,
    run_transient,
    set_devsim_quiet,
)

set_devsim_quiet(True)
CALLS = {
    "field1d": lambda: run_field(DeviceConfig(), -30.0),
    "transient": lambda: run_transient(DeviceConfig()),
    "dark": lambda: run_dark_current(DeviceConfig(), 0.0, -60.0, 7, solver="devsim"),
}

for name in sys.argv[1:] or ["field1d"]:
    fn = CALLS[name]
    fn()
    pr = cProfile.Profile()
    t = time.perf_counter()
    pr.enable()
    fn()
    pr.disable()
    wall = time.perf_counter() - t
    st = pstats.Stats(pr)
    agg = {
        "solve_total": 0.0,
        "callback": 0.0,
        "init(dlopen)": 0.0,
        "factor": 0.0,
        "lin_solve": 0.0,
    }
    n_solve = 0
    for (f, line, fname), (cc, nc, tt, ct, callers) in st.stats.items():
        if fname == "<built-in method devsim.devsim_py3.solve>":
            agg["solve_total"] += ct
            n_solve += nc
        if f.endswith("umfshim.py"):
            if fname == "local_solver_callback":
                agg["callback"] += ct
            elif fname == "initialize_umfpack":
                agg["init(dlopen)"] += ct
            elif fname == "factor":
                agg["factor"] += ct
            elif fname == "solve":
                agg["lin_solve"] += ct
    print(
        f"{name}: wall {wall:.3f}s, {n_solve} devsim.solve, "
        + ", ".join(f"{k} {v:.3f}s" for k, v in agg.items())
    )
