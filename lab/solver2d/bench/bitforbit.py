"""Dump default-model outputs (no n+) to an .npz, for bit-for-bit comparison across commits.

usage: PYTHONPATH=<exported etna root> uv run python lab/solver2d/bench/bitforbit.py OUT.npz
       uv run python lab/solver2d/bench/bitforbit.py --compare A.npz B.npz
"""

from __future__ import annotations

import sys

import numpy as np


def dump(path):
    import devsim

    import etna
    from etna.api.device import DeviceConfig
    from etna.api.simulation import run_cce, run_cv, run_dark_current
    from etna.core.batched import node_positions
    from etna.core.charge_collection import cce_vs_bias
    from etna.core.cv_analysis import cv_sweep
    from etna.core.dark_current import create_dark_current_device, dark_current_sweep
    from etna.core.drift_diffusion import create_dd_device

    print("etna from", etna.__file__, flush=True)
    out = {}
    out["x_nodes"] = np.asarray(node_positions(10e-4, 1e-4))
    info = create_dd_device(device_name="bfb_cv")
    r = cv_sweep(info, -np.array([0.0, 5.0, 10.0, 30.0, 50.0]), area=0.05)
    out["W"] = np.asarray(r["depletion_widths"])
    out["C"] = np.asarray(r["capacitance"])
    devsim.delete_device(device="bfb_cv")
    c = cce_vs_bias(-np.arange(10.0, 55.0, 10.0))
    out["cce"] = np.asarray(c["cce_values"])
    out["I_coll"] = np.asarray(c["I_collected"])
    info = create_dark_current_device(device_name="bfb_dark")
    d = dark_current_sweep(info, np.array([-10.0, -30.0, -60.0]), area=0.04)
    out["dark"] = np.asarray(d["I_total"])
    out["dark_tat"] = np.asarray(d["I_TAT"])
    devsim.delete_device(device="bfb_dark")
    cfg = DeviceConfig()
    out["api_cv"] = np.asarray(run_cv(cfg).y)
    out["api_cce"] = np.asarray(run_cce(cfg).y)
    out["api_dark"] = np.asarray(run_dark_current(cfg).y)
    np.savez(path, **out)


def compare(a, b):
    A, B = np.load(a), np.load(b)
    worst = 0.0
    for k in A.files:
        same = A[k].shape == B[k].shape and np.array_equal(A[k], B[k])
        rel = (
            np.max(np.abs(A[k] - B[k]) / np.maximum(np.abs(A[k]), 1e-300))
            if A[k].shape == B[k].shape
            else np.inf
        )
        worst = max(worst, rel)
        print(f"{k:10s} identical={same} max_rel={rel:.3e}")
    print("ALL IDENTICAL" if worst == 0.0 else f"WORST {worst:.3e}")


if __name__ == "__main__":
    if sys.argv[1] == "--compare":
        compare(sys.argv[2], sys.argv[3])
    else:
        dump(sys.argv[1])
