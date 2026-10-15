"""Batched solver on P1 and P2 (same mesh and bias points as the devsim baseline).

Run from lab/gpu_solver (after bench_devsim.py has written results/devsim.json):
    uv run python bench/bench_batched.py --backend numpy
    uv run python bench/bench_batched.py --backend mlx

Reports wall time and max relative W error vs devsim's DD results, and writes
results/batched_<backend>.json.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from problems import P1, P2_BIASES, p2_devices  # noqa: E402

from gpusolver.sweep import make_batch, reverse_sweep  # noqa: E402

HERE = os.path.dirname(__file__)
X = np.load(os.path.join(HERE, "mesh_devsim_default.npy"))
XJ = 1e-4
EPI = 10e-4


def _time(fn, repeats=3):
    best = np.inf
    out = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        out = fn()
        best = min(best, time.perf_counter() - t0)
    return best, out


def run_p1(backend, devsim_res, tol):
    V = np.array(devsim_res["P1"]["voltages"])
    W_ref = np.array(devsim_res["P1"]["W_cm"])
    kw = dict(
        N_A=P1["N_A"],
        N_D_junction=P1["N_D_junction"],
        N_D_bulk=P1["N_D_bulk"],
        L_transition_cm=P1["L_transition_um"] * 1e-4,
    )
    out = {}
    for label, step in (("0.5V steps (as devsim)", 0.5), ("direct jumps", None)):
        pb = make_batch(X, XJ, backend=backend, **kw)
        t, sw = _time(
            lambda: reverse_sweep(pb, V, XJ, EPI, V_step=step, tol=tol, keep_psi=False)
        )
        err = float(np.abs(sw.W[:, 0] - W_ref).max() / W_ref.min())
        out[label] = dict(
            t_s=t,
            newton_its=sw.newton_iterations,
            max_rel_W_err=err,
            converged=bool(sw.converged.all()),
        )
        print(
            f"P1 [{backend}] {label}: {t*1e3:.0f} ms, {sw.newton_iterations} Newton its, max rel W err {err:.2e}, converged={sw.converged.all()}"
        )
    return out


def run_p2(backend, devsim_res, tol):
    devs = devsim_res["P2"]["devices"]
    nb = np.array([d["N_D_bulk"] for d in devs])
    nj = np.array([d["N_D_junction"] for d in devs])
    assert np.allclose(list(zip(nb, nj)), list(p2_devices()))
    pb = make_batch(X, XJ, 1e19, nj, nb, 0.987e-4, backend=backend)
    t, sw = _time(
        lambda: reverse_sweep(pb, P2_BIASES, XJ, EPI, tol=tol, keep_psi=False)
    )
    W_ref = np.array([d["W_cm"] for d in devs]).T  # (K, B)
    err = float((np.abs(sw.W - W_ref) / W_ref).max())
    print(
        f"P2 [{backend}] 100 devices x 10 biases: {t*1e3:.0f} ms, {sw.newton_iterations} Newton its, max rel W err {err:.2e}, converged={sw.converged.all()}"
    )
    return dict(
        t_s=t,
        newton_its=sw.newton_iterations,
        max_rel_W_err=err,
        converged=bool(sw.converged.all()),
    )


def run_scaling(backend, tol, sizes=(1, 10, 100, 1000, 10000)):
    """Wall time of one 0 -> -45 V sweep (10 points) vs batch size, random dopings."""
    rng = np.random.default_rng(0)
    out = {}
    for B in sizes:
        nb = 10 ** rng.uniform(13.5, 14.5, B)
        nj = 10 ** rng.uniform(15.0, 15.7, B)
        pb = make_batch(X, XJ, 1e19, nj, nb, 0.987e-4, backend=backend)
        t, sw = _time(
            lambda: reverse_sweep(pb, P2_BIASES, XJ, EPI, tol=tol, keep_psi=False),
            repeats=2,
        )
        out[B] = dict(
            t_s=t, newton_its=sw.newton_iterations, converged=bool(sw.converged.all())
        )
        print(
            f"scaling [{backend}] B={B:6d}: {t*1e3:8.0f} ms  ({t/B*1e3:.3f} ms/device, {sw.newton_iterations} its, converged={sw.converged.all()})"
        )
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", default="numpy", choices=["numpy", "mlx"])
    ap.add_argument("--no-scaling", action="store_true")
    args = ap.parse_args()
    tol = 1e-8 if args.backend == "numpy" else 1e-4
    with open(os.path.join(HERE, "results", "devsim.json")) as f:
        devsim_res = json.load(f)
    results = {
        "P1": run_p1(args.backend, devsim_res, tol),
        "P2": run_p2(args.backend, devsim_res, tol),
    }
    if not args.no_scaling:
        results["scaling"] = run_scaling(args.backend, tol)
    path = os.path.join(HERE, "results", f"batched_{args.backend}.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=1)
    print("wrote", path)
