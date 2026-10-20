"""Who calls ctypes dlopen during devsim solves? (run_field 1D, profiled)

usage: uv run python lab/perf_profile/who_dlopens.py
"""

from __future__ import annotations

import cProfile
import pstats

from etna import DeviceConfig, run_field, set_devsim_quiet

set_devsim_quiet(True)
run_field(DeviceConfig(), -30.0)
pr = cProfile.Profile()
pr.enable()
run_field(DeviceConfig(), -30.0)
pr.disable()
st = pstats.Stats(pr)
st.print_callers("dlopen")
st.print_callers("ctypes/__init__.py:343")
st.print_callers("LoadLibrary")
