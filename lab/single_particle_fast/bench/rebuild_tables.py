"""Rebuild the committed CCE(LET) tables with build_cce_let_table(engine="fast") into an
output directory (never over data/) and compare them with data/cce_let_table_*.json.

    uv run python lab/single_particle_fast/bench/rebuild_tables.py OUT_DIR [n_jobs] [scheme] [width]

scheme: "adaptive" (default, TR-BDF2, converged in time) or "devsim_grid" (devsim's own
10 % BDF1 grid and estimator: the twin of the committed tables). width: lateral cut in
um beyond the track (default 10; "full" = the whole half-width). The build keywords are
those of notebook 16 (lab/mobility_caxis/bench/rebuild_let_tables.py TABLES).
"""

import json
import os
import sys
import time

import numpy as np

from etna.core.devsim_lock import set_devsim_quiet
from etna.core.single_particle import build_cce_let_table, save_cce_let_table

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "mobility_caxis", "bench"))
from rebuild_let_tables import TABLES  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "..", "data")

if __name__ == "__main__":
    set_devsim_quiet(True)
    out_dir = sys.argv[1]
    n_jobs = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    scheme = sys.argv[3] if len(sys.argv) > 3 else "adaptive"
    width = sys.argv[4] if len(sys.argv) > 4 else "10"
    width = None if width == "full" else float(width)
    os.makedirs(out_dir, exist_ok=True)
    for name, (kw, geo) in TABLES.items():
        t0 = time.time()
        table = build_cce_let_table(
            n_jobs=n_jobs, engine="fast",
            fast_options={"scheme": scheme, "lateral_width_um": width}, **kw,
        )  # fmt: skip
        wall = time.time() - t0
        path = os.path.join(out_dir, name)
        save_cce_let_table(table, path, geometry_info=geo)
        with open(os.path.join(DATA, name)) as f:
            ref = json.load(f)
        dc = np.asarray(table["CCE_raw"]) - np.asarray(ref["CCE_raw"])
        dt = (
            np.asarray(table["t_collection_ns"]) / np.asarray(ref["t_collection_ns"])
            - 1
        )
        dq = np.asarray(table["Q_generated_fC"]) / np.asarray(ref["Q_generated_fC"]) - 1
        print(f"{name} [{scheme}, cut {width}, n_jobs {n_jobs}]: {wall:.1f} s "
              f"(build {table.attrs['build_s']:.1f} s, {table.attrs['nodes']} nodes)")  # fmt: skip
        print(f"  Q_gen ratio-1 max {np.abs(dq).max():.1e}; all converged "
              f"{bool(table['converged'].all())}")  # fmt: skip
        for L, c, r, d, t in zip(
            table["LET_keV_um"], table["CCE_raw"], ref["CCE_raw"], dc, dt
        ):
            print(f"  LET {L:9.4f}: CCE fast {c:.10f} table {r:.10f} diff {d:+.2e} "
                  f"t95 rel diff {t:+.2e}")  # fmt: skip
        print(f"  max |dCCE| {np.abs(dc).max():.2e}, max |dt95|/t95 {np.abs(dt).max():.2e}",
              flush=True)  # fmt: skip
