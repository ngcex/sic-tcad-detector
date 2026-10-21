"""Rebuild the committed CCE(LET) tables (data/cce_let_table_{100,300}um.json) with a
given low-field mobility set, with the exact build keywords of notebook 16
(scripts/create_notebook_16.py), into an output directory (never over data/).

    uv run python lab/mobility_caxis/bench/rebuild_let_tables.py OUT_DIR [mobility] [n_jobs] [100|300]

mobility: "c_axis" (default) or "legacy_perp". The 2D device builder takes the set via
device_kwargs={"mobility": ...}.
"""

import json
import os
import sys
import time

from etna.core.devsim_lock import set_devsim_quiet
from etna.core.single_particle import build_cce_let_table, save_cce_let_table

TABLES = {
    "cce_let_table_100um.json": (
        dict(
            half_width_um=50.0,
            V_bias=50.0,
            n_let_points=20,
            let_min=0.5,
            let_max=500.0,
            x_ion_cm=0.0,
        ),  # fmt: skip
        {"half_width_um": 50.0, "epi_um": 10.0},
    ),
    "cce_let_table_300um.json": (
        dict(
            half_width_um=150.0,
            V_bias=50.0,
            n_let_points=10,
            let_min=0.5,
            let_max=500.0,
            x_ion_cm=0.0,
        ),  # fmt: skip
        {"half_width_um": 150.0, "epi_um": 10.0},
    ),
}

if __name__ == "__main__":
    set_devsim_quiet(True)
    out_dir = sys.argv[1]
    mobility = sys.argv[2] if len(sys.argv) > 2 else "c_axis"
    n_jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    only = sys.argv[4] if len(sys.argv) > 4 else None
    os.makedirs(out_dir, exist_ok=True)
    for name, (kw, geo) in TABLES.items():
        if only is not None and f"_{only}um" not in name:
            continue
        t0 = time.time()
        dk = None if mobility == "c_axis" else {"mobility": mobility}
        table = build_cce_let_table(n_jobs=n_jobs, device_kwargs=dk, **kw)
        path = os.path.join(out_dir, name)
        save_cce_let_table(table, path, geometry_info=geo)
        with open(path) as f:
            d = json.load(f)
        print(f"{name} [{mobility}]: {time.time() - t0:.0f} s, CCE "
              f"{min(d['CCE']):.7f}..{max(d['CCE']):.7f}", flush=True)  # fmt: skip
