"""build_cce_let_table(method="affine") vs method="full", row by row, plus the stored table.

    uv run python lab/engine_box/bench/let_table_affine.py HW N STORED.json OUTDIR [skip_full]

"full" simulates every LET on one device (bit-identical to the former build-per-LET
loop, let_table_build_once.py); "affine" simulates let_min / let_max / the middle LET
and fills the rest from Q_col = a * Q_gen + b. Grid: np.logspace(0.5, 500, N), the
notebook-16 grid. Saves both tables to OUTDIR (scratch, never data/).
"""

import json
import os
import sys
import time

import numpy as np
import pandas as pd

import etna
from etna.core.single_particle import build_cce_let_table

etna.set_devsim_quiet(True)

hw, n, stored_path, outdir = (
    float(sys.argv[1]),
    int(sys.argv[2]),
    sys.argv[3],
    sys.argv[4],
)
skip_full = len(sys.argv) > 5
os.makedirs(outdir, exist_ok=True)
kw = dict(half_width_um=hw, V_bias=50.0, n_let_points=n, let_min=0.5, let_max=500.0)
tables, walls = {}, {}
for method in ("affine",) if skip_full else ("affine", "full"):
    path = os.path.join(outdir, f"table_hw{hw:g}_{method}.csv")
    t0 = time.perf_counter()
    tables[method] = build_cce_let_table(method=method, **kw)
    walls[method] = time.perf_counter() - t0
    tables[method].to_csv(path, index=False)
    print(f"{method}: {walls[method]:.1f} s ({walls[method] / 60:.1f} min)", flush=True)
if skip_full:
    path = os.path.join(outdir, f"table_hw{hw:g}_full.csv")
    tables["full"] = pd.read_csv(path)

with open(stored_path) as f:
    stored = json.load(f)
s_cce = np.array(stored["CCE_raw"])
a, fu = tables["affine"], tables["full"]
print(
    f"{'LET':>9} {'full':>14} {'affine':>14} {'aff-full':>9} {'stored':>14} "
    f"{'full-st':>9} {'aff-st':>9} sim t95_f t95_a"
)
for k in range(n):
    cf, ca = fu["CCE_raw"][k], a["CCE_raw"][k]
    print(
        f"{a['LET_keV_um'][k]:9.4g} {cf:14.10f} {ca:14.10f} {ca - cf:9.1e} "
        f"{s_cce[k]:14.10f} {cf - s_cce[k]:9.1e} {ca - s_cce[k]:9.1e} "
        f"{'Y' if a['simulated'][k] else 'n'} {fu['t_collection_ns'][k]:.4f} "
        f"{a['t_collection_ns'][k]:.4f}"
    )
d = np.abs(a["CCE_raw"] - fu["CCE_raw"])
print(
    f"max |affine - full| = {d.max():.2e}; max |full - stored| = "
    f"{np.abs(fu['CCE_raw'] - s_cce).max():.2e}; max |affine - stored| = "
    f"{np.abs(a['CCE_raw'] - s_cce).max():.2e}"
)
print(
    f"max |t95 affine - full| = {np.abs(a['t_collection_ns'] - fu['t_collection_ns']).max():.2e} ns"
)
if "full" in walls:
    print(f"wall full/affine = {walls['full'] / walls['affine']:.2f}x")
