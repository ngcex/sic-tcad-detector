"""Import-time breakdown: fresh interpreter per module, wall + cpu.

usage: uv run python lab/perf_profile/import_time.py [module ...]
"""

from __future__ import annotations

import subprocess
import sys

MODS = sys.argv[1:] or [
    "numpy",
    "scipy.sparse",
    "scipy.stats",
    "scipy.optimize",
    "pandas",
    "matplotlib",
    "matplotlib.pyplot",
    "numba",
    "devsim",
    "etna.core.batched",
    "etna.api.device",
    "etna",
]

CODE = (
    "import time,resource;t=time.perf_counter();import {m};"
    "w=time.perf_counter()-t;r=resource.getrusage(resource.RUSAGE_SELF);"
    "import sys;heavy=[k for k in ('matplotlib.pyplot','pandas','scipy.stats','devsim','numba','plotly','streamlit') if k in sys.modules];"
    "print(f'{{w:.3f}} {{r.ru_utime+r.ru_stime:.3f}} '+','.join(heavy))"
)

if __name__ == "__main__":
    for m in MODS:
        best = None
        for _ in range(3):
            out = subprocess.run(
                [sys.executable, "-c", CODE.format(m=m)], capture_output=True, text=True
            )
            line = (
                out.stdout.strip().splitlines()[-1]
                if out.stdout.strip()
                else "ERR " + out.stderr[-200:]
            )
            w = float(line.split()[0]) if line[0].isdigit() else 1e9
            if best is None or w < best[0]:
                best = (w, line)
        print(f"{m:28s} wall_s cpu_s loaded: {best[1]}")
