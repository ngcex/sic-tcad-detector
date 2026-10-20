"""Static scan: heavy calls per notebook code cell (no execution).

usage: uv run python lab/perf_profile/nb_calls.py [notebook-stem ...]
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

NB = Path(__file__).resolve().parents[2] / "notebooks"
PAT = re.compile(r"\b([A-Za-z_0-9]+)\(")
HEAVY = re.compile(
    r"(_2d|sweep|scan|table|device|_vs_|cce|run_|transient|dark|simulate|ramp|solve|calibrat|structure|mesa|guard|electrode)"
)

stems = sys.argv[1:] or [p.stem for p in sorted(NB.glob("*.ipynb"))]
for stem in stems:
    d = json.loads((NB / f"{stem}.ipynb").read_text())
    print("=====", stem)
    for k, c in enumerate(d["cells"]):
        if c["cell_type"] != "code":
            continue
        s = "".join(c["source"])
        calls = sorted(
            {m for m in PAT.findall(s) if HEAVY.search(m) and not m[0].isupper()}
        )
        loops = len(re.findall(r"^\s*for .* in ", s, re.M))
        if calls:
            print(
                f"  cell {k:2d} loops={loops} {'n_jobs ' if 'n_jobs' in s else ''}{calls}"
            )
