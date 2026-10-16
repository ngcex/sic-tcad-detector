"""GATE: notebook 20's 36-config grid, microdosimetric_sweep(solver="devsim") vs ("auto").

usage: uv run python lab/solver2d/bench/gate_nb20.py [--jobs N]
Both runs n_jobs=1 (serial), same process, devsim first. Per-config wall time via a
timing wrapper around etna.core.optimization._sweep_one (the serial path calls it
through the module global). Optional --jobs N adds an auto run with n_jobs=N.
Writes lab/solver2d/notes/gate_nb20.npz-compatible CSVs: gate_devsim.csv, gate_auto.csv.
"""

import argparse
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

import etna.core.optimization as opt

GRID = dict(
    half_widths_um=[25, 50, 100, 150],
    epi_thicknesses_cm=[10e-4],
    N_D_bulks=[5e13, 8.5e13, 5e14],
    V_biases=[20, 50, 80],
    n_lateral_points=10,
)
OUT = Path("lab/solver2d/notes")
KEY = ["half_width_um", "epi_thickness_cm", "N_D_bulk", "V_bias"]

_orig = opt._sweep_one
_times = []


def _timed(task):
    t = time.perf_counter()
    r = _orig(task)
    r["t_config_s"] = time.perf_counter() - t
    return r




def main():
    opt._sweep_one = _timed


    def run(solver, n_jobs=1):
        t = time.perf_counter()
        df = opt.microdosimetric_sweep(**GRID, n_jobs=n_jobs, solver=solver)
        return df, time.perf_counter() - t


    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=0)
    args = ap.parse_args()

    print(f"load {os.getloadavg()}", flush=True)
    dev, t_dev = run("devsim")
    dev.to_csv(OUT / "gate_devsim.csv", index=False)
    print(f"devsim total {t_dev:.1f}s", flush=True)
    auto, t_auto = run("auto")
    auto.to_csv(OUT / "gate_auto.csv", index=False)
    print(f"auto total {t_auto:.2f}s  speedup {t_dev / t_auto:.1f}x", flush=True)

    m = dev.merge(auto, on=KEY, suffixes=("_d", "_a"))
    for c in ("center_cce", "edge_cce", "edge_center_ratio", "cce_std"):
        d = (m[f"{c}_a"] - m[f"{c}_d"]).abs()
        print(f"max |d {c}| = {d.max():.2e}")
    nan_d, nan_a = m["center_cce_d"].isna(), m["center_cce_a"].isna()
    print("NaN rows identical:", bool((nan_d == nan_a).all()), " count", int(nan_d.sum()))
    print("is_valid identical:", bool((m["is_valid_d"] == m["is_valid_a"]).all()))
    print(
        "ranking identical:",
        dev[KEY].reset_index(drop=True).equals(auto[KEY].reset_index(drop=True)),
    )
    print("solver used (auto):", auto["solver"].value_counts().to_dict())
    print(
        m[
            KEY
            + [
                "solver_a",
                "t_config_s_d",
                "t_config_s_a",
                "center_cce_d",
                "center_cce_a",
                "edge_center_ratio_d",
                "edge_center_ratio_a",
                "is_valid_d",
            ]
        ].to_string()
    )
    if args.jobs > 1:
        auto2, t2 = run("auto", n_jobs=args.jobs)
        same = np.allclose(
            auto2["center_cce"].to_numpy(float),
            auto["center_cce"].to_numpy(float),
            equal_nan=True,
            rtol=0,
            atol=0,
        )
        print(
            f"auto n_jobs={args.jobs}: {t2:.2f}s (serial {t_auto:.2f}s), identical {same}"
        )


if __name__ == "__main__":  # spawn workers re-import this module
    main()
