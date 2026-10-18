"""Parity + timing of build_cce_let_table against the build-per-LET reference.

    uv run python lab/engine_box/bench/let_table_build_once.py HW REF.npz N LET_MIN LET_MAX [method]

REF.npz from let_table_reference.py at the same LETs (np.logspace(LET_MIN, LET_MAX, N)
must reproduce its keys). Prints per-point |dCCE| and the wall time.
"""

import sys
import time

import numpy as np

import etna
from etna.core.single_particle import build_cce_let_table

etna.set_devsim_quiet(True)

hw, ref_path = float(sys.argv[1]), sys.argv[2]
n, lo, hi = int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
kw = {"method": sys.argv[6]} if len(sys.argv) > 6 else {}
ref = np.load(ref_path)
t0 = time.perf_counter()
table = build_cce_let_table(
    half_width_um=hw, V_bias=50.0, n_let_points=n, let_min=lo, let_max=hi, **kw
)
wall = time.perf_counter() - t0
print(f"build_cce_let_table {kw}: {wall:.1f} s for {n} points")
worst = 0.0
for _, row in table.iterrows():
    key = f"{row['LET_keV_um']:g}"
    line = (
        f"LET {key}: CCE_raw={row['CCE_raw']:.10f} t95={row['t_collection_ns']:.4f} ns"
    )
    if f"cce_{key}" in ref:
        d = abs(row["CCE_raw"] - float(ref[f"cce_{key}"]))
        worst = max(worst, d)
        line += f" ref={float(ref[f'cce_{key}']):.10f} |d|={d:.2e}"
    print(line)
print(f"max |dCCE| vs reference = {worst:.2e}")
