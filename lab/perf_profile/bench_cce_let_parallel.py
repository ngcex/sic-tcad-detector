"""Timing + parity: build_cce_let_table serial vs n_jobs, on a reduced size.

usage: uv run python lab/perf_profile/bench_cce_let_parallel.py [hw_um] [n_let] [n_jobs]
  default: hw 2 um, 8 LETs, n_jobs 4 (the full notebook-16 tables, hw 50/150 um with
  20/10 LETs, take ~15/26 min serial: not run here).
"""

from __future__ import annotations

import sys
import time

import numpy as np

from etna import set_devsim_quiet
from etna.core.single_particle import build_cce_let_table

set_devsim_quiet(True)
hw = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
n = int(sys.argv[2]) if len(sys.argv) > 2 else 8
nj = int(sys.argv[3]) if len(sys.argv) > 3 else 4
kw = dict(half_width_um=hw, V_bias=50.0, n_let_points=n, let_min=0.5, let_max=500.0)

t = time.perf_counter()
c = time.process_time()
ser = build_cce_let_table(**kw)
t_ser, c_ser = time.perf_counter() - t, time.process_time() - c
t = time.perf_counter()
par = build_cce_let_table(n_jobs=nj, **kw)
t_par = time.perf_counter() - t
d = np.nanmax(np.abs(par["CCE_raw"].to_numpy() - ser["CCE_raw"].to_numpy()))
print(
    f"hw {hw} um, {n} LETs: serial {t_ser:.1f} s (parent cpu {c_ser:.1f} s), "
    f"n_jobs={nj} {t_par:.1f} s -> x{t_ser / t_par:.2f}; simulated {int(ser['simulated'].sum())}/{n}; "
    f"same flags {bool((par['simulated'] == ser['simulated']).all())}; max |dCCE_raw| {d:.2e}"
)
