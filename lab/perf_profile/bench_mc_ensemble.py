"""process_mc_ensemble: per-event CCE lookup loop vs the vectorized cce_interp.

usage: uv run python lab/perf_profile/bench_mc_ensemble.py [n_events ...]
Old path = the same cce_interp with ``vectorized`` removed (falls back to the loop).
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd

from etna.core.mc_coupling import process_mc_ensemble
from etna.core.single_particle import load_cce_let_table

interp, _ = load_cce_let_table("data/cce_let_table_100um.json")


def scalar_only(let):  # the pre-change behaviour: one np.interp per event
    return interp(let)


rng = np.random.default_rng(0)
for n in [int(a) for a in sys.argv[1:]] or [10_000, 100_000, 1_000_000]:
    ev = pd.DataFrame(
        {
            "event_id": np.repeat(np.arange(n), 3),
            "edep_keV": rng.lognormal(3.0, 1.2, 3 * n),
        }
    )
    best = {}
    for name, fn in (("loop", scalar_only), ("vectorized", interp)):
        ts = []
        for _ in range(3):
            t = time.process_time()
            r = process_mc_ensemble(ev, fn)
            ts.append(time.process_time() - t)
        best[name] = (min(ts), r["event_CCE"])
    same = np.array_equal(best["loop"][1], best["vectorized"][1])
    print(
        f"{n:>8d} events: loop {best['loop'][0]:.3f} s -> vectorized {best['vectorized'][0]:.3f} s cpu "
        f"(x{best['loop'][0] / best['vectorized'][0]:.1f}); bit-identical {same}"
    )
