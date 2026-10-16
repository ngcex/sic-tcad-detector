"""Box-widening parity grid: _sweep_one(solver="devsim") vs ("auto") on a chunk of points.

usage: uv run python lab/solver2d/bench/parity_box.py START STOP   (indices into POINTS)
appends rows to lab/solver2d/notes/parity_box.csv
"""

import csv
import sys
import time
from pathlib import Path

from etna.core.optimization import _sweep_one

D = 8.823059e13  # the default N_D_bulk, passed explicitly (_sweep_one needs a number)
POINTS = [  # (hw um, epi um, N_D_bulk, V)
    (25, 5, 5e13, 5),
    (25, 5, D, 20),
    (25, 5, 5e14, 35),
    (25, 5, 5e14, 37),
    (50, 5, 5e13, 10),
    (25, 5, 8.5e13, 50),
    (25, 7.5, 5e13, 20),
    (25, 7.5, 5e14, 65),
    (25, 7.5, D, 35),
    (50, 7.5, 5e14, 50),
    (25, 10, 5e13, 5),
    (25, 10, D, 10),
    (25, 10, 5e14, 35),
    (50, 10, 8.5e13, 65),
    (25, 10, 5e14, 75),
    (25, 15, 5e13, 80),
    (25, 15, D, 50),
    (25, 15, 5e14, 5),
    (50, 15, 5e13, 35),
    (25, 20, 5e13, 80),
    (25, 20, D, 20),
    (25, 20, 5e14, 10),
    (50, 20, 5e14, 65),
    (25, 20, 5e13, 5),
]
FIELDS = [
    "hw",
    "epi_um",
    "nd",
    "V",
    "solver_a",
    "center_d",
    "center_a",
    "edge_d",
    "edge_a",
    "ratio_d",
    "ratio_a",
    "valid_d",
    "valid_a",
    "t_d",
    "t_a",
]
out = Path("lab/solver2d/notes/parity_box.csv")
new = not out.exists()
with out.open("a", newline="") as f:
    w = csv.writer(f)
    if new:
        w.writerow(FIELDS)
    for hw, epi, nd, V in POINTS[int(sys.argv[1]) : int(sys.argv[2])]:
        task = (0, 1, float(hw), epi * 1e-4, nd, float(V), 10)
        t = time.perf_counter()
        a = _sweep_one(task + ("auto",))
        ta = time.perf_counter() - t
        t = time.perf_counter()
        d = _sweep_one(task + ("devsim",))
        td = time.perf_counter() - t
        row = [
            hw,
            epi,
            nd,
            V,
            a["solver"],
            d["center_cce"],
            a["center_cce"],
            d["edge_cce"],
            a["edge_cce"],
            d["edge_center_ratio"],
            a["edge_center_ratio"],
            d["is_valid"],
            a["is_valid"],
            round(td, 2),
            round(ta, 3),
        ]
        w.writerow(row)
        f.flush()
        print("ROW", row, f"dC={a['center_cce'] - d['center_cce']:.2e}", flush=True)
