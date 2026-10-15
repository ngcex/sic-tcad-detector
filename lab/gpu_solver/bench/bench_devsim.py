"""Baseline: time the existing ETNA/devsim path on P1 and P2.

Run from lab/gpu_solver:
    uv run python bench/bench_devsim.py            # P1 + P2 + P3 (~4 min)
    uv run python bench/bench_devsim.py --only P3   # one problem, merged into results/devsim.json

Writes bench/results/devsim.json with wall times and the reference W(V) for
every device (used to cross-check the batched solvers at the same points).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from problems import P1, P2_BIASES, P3_BIASES, p2_devices  # noqa: E402

from gpusolver.ref_devsim import (
    build_ref_device,
    delete_ref_device,
    ref_cv_sweep,
)  # noqa: E402


def run_p1():
    kw = {
        k: P1[k]
        for k in (
            "epi_thickness_um",
            "N_A",
            "N_D_junction",
            "N_D_bulk",
            "L_transition_um",
            "T",
        )
    }
    t0 = time.perf_counter()
    ref = build_ref_device(**kw)
    t_build = time.perf_counter() - t0
    V = np.arange(0.0, P1["V_stop"] + P1["V_step"] / 2, P1["V_step"])
    t0 = time.perf_counter()
    r = ref_cv_sweep(ref, V)
    t_sweep = time.perf_counter() - t0
    delete_ref_device(ref)
    n_solves = len(r["voltages"]) - 1
    out = dict(
        nodes=int(len(ref.x)),
        t_build_s=t_build,
        t_sweep_s=t_sweep,
        n_bias_points_reached=int(len(r["voltages"])),
        V_min_reached=float(r["voltages"].min()),
        truncated=bool(r["truncated"]),
        ms_per_solve=1e3 * t_sweep / max(n_solves, 1),
        voltages=r["voltages"].tolist(),
        W_cm=r["depletion_widths"].tolist(),
    )
    print(
        f"P1: {out['nodes']} nodes, build {t_build:.3f}s, ramp 0->{out['V_min_reached']:.1f}V "
        f"in {t_sweep:.2f}s ({out['ms_per_solve']:.1f} ms/solve), truncated={out['truncated']}"
    )
    return out


def run_p2():
    devices = []
    t_total0 = time.perf_counter()
    for i, (nb, nj) in enumerate(p2_devices()):
        t0 = time.perf_counter()
        ref = build_ref_device(N_D_bulk=nb, N_D_junction=nj)
        r = ref_cv_sweep(ref, P2_BIASES)
        delete_ref_device(ref)
        dt = time.perf_counter() - t0
        devices.append(
            dict(
                N_D_bulk=nb,
                N_D_junction=nj,
                t_s=dt,
                truncated=bool(r["truncated"]),
                voltages=r["voltages"].tolist(),
                W_cm=r["depletion_widths"].tolist(),
            )
        )
        if i % 10 == 0:
            print(f"  P2 device {i:3d}/100  {dt:.2f}s  truncated={r['truncated']}")
    t_total = time.perf_counter() - t_total0
    n_points = sum(len(d["voltages"]) for d in devices)
    print(f"P2: 100 devices x 10 biases = {n_points} points reached in {t_total:.1f}s")
    return dict(t_total_s=t_total, n_points_reached=n_points, devices=devices)


def run_p3():
    from etna.core.charge_collection import cce_vs_bias

    from gpusolver.ref_devsim import quiet

    devices = []
    t_total0 = time.perf_counter()
    for i, (nb, nj) in enumerate(p2_devices()):
        t0 = time.perf_counter()
        with quiet():
            r = cce_vs_bias(P3_BIASES, device_kwargs=dict(N_D_junction=nj, N_D_bulk=nb, L_transition=0.987e-4))
        dt = time.perf_counter() - t0
        devices.append(
            dict(N_D_bulk=nb, N_D_junction=nj, t_s=dt, truncated=bool(r["truncated"]),
                 voltages=r["voltages"].tolist(), cce=r["cce_values"].tolist(),
                 I_collected=r["I_collected"].tolist(), I_generated=float(r["I_generated"]))
        )
        if i % 10 == 0:
            print(f"  P3 device {i:3d}/100  {dt:.2f}s  truncated={r['truncated']}")
    t_total = time.perf_counter() - t_total0
    n_points = sum(len(d["voltages"]) for d in devices)
    print(f"P3: 100 devices x 5 biases CCE = {n_points} points reached in {t_total:.1f}s")
    return dict(t_total_s=t_total, n_points_reached=n_points, devices=devices)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["P1", "P2", "P3"], default=None)
    args = ap.parse_args()
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results", "devsim.json")
    results = json.load(open(path)) if os.path.exists(path) else {}
    runners = {"P1": run_p1, "P2": run_p2, "P3": run_p3}
    for name, fn in runners.items():
        if args.only is None or args.only == name:
            results[name] = fn()
    with open(path, "w") as f:
        json.dump(results, f)
    print("wrote", path)
