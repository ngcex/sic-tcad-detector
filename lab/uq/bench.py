"""etna.uq bench: throughput (samples/s), batching gain, Ishigami convergence, and the
bands / Sobol indices of the default dosimeter (numbers in notes/LOG.md).

    uv run python lab/uq/bench.py            # ~2-3 min on a 10-core M4
    uv run python lab/uq/bench.py --quick    # throughput + Ishigami only (~30 s)

Throughput = model evaluations / wall time of the evaluator calls (warm numba cache,
kernels warmed once per observable first), excluding the nominal / mesh-error extras.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from etna import DeviceConfig, run_cv
from etna.uq import (
    ParameterSpace,
    SaltelliDesign,
    UncertainParameter,
    indices_table,
    ishigami,
    ishigami_indices,
    run_uq,
    sobol_indices,
)

CFG = DeviceConfig()


def ishigami_table():
    ref = ishigami_indices()
    sp = ParameterSpace([UncertainParameter(f"x{i}", 0, "uniform", -np.pi, np.pi, "model:x")
                         for i in range(3)])  # fmt: skip
    print("Ishigami (a=7, b=0.1): max |error| of S1 / ST vs closed form, 95 % CI width")
    for N in (128, 512, 2048, 8192):
        d = SaltelliDesign(sp, seed=1)
        P = d.new_points(N)
        s = sobol_indices(*SaltelliDesign.split(ishigami(P), N, 3), n_boot=200)
        print(
            f"  N={N:5d}  S1 err {np.abs(s.S1 - ref['S1']).max():.4f}  "
            f"ST err {np.abs(s.ST - ref['ST']).max():.4f}  CI {s.max_ci_width:.3f}"
        )


def throughput():
    print("\nThroughput (default registry per observable, Saltelli, warm):")
    runs = [("v_fd", 64), ("cv", 64), ("cce", 32), ("alpha_charge", 32),
            ("dark_current", 32), ("transit_time", 64), ("breakdown", 16),
            ("compensation_fluence", 1024), ("microdosimetry", 256), ("cce_fluence", 2)]  # fmt: skip
    for obs, n in runs:
        kw = {"fluences": np.geomspace(1e11, 1e14, 4)} if obs == "cce_fluence" else {}
        run_uq(CFG, obs, n=2, numerical_error=False, model_form=False, n_boot=0, **kw)
        r = run_uq(
            CFG, obs, n=n, numerical_error=False, model_form=False, n_boot=0, **kw
        )
        md = r.metadata
        print(
            f"  {obs:21s} n={n:5d} evals {md['n_evaluations']:6d}  "
            f"{md['samples_per_s']:9.1f} /s  groups {md['info'].get('n_groups', '-')}  "
            f"failed {md['n_failed']}"
        )


def batching_gain():
    print(
        "\nBatching: C-V (31 biases), doping-only inputs (one mesh group) vs a run_cv loop"
    )
    params = ["N_D_bulk", "N_D_junction", "L_transition_um"]
    run_uq(
        CFG, "cv", params=params, n=4, numerical_error=False, model_form=False, n_boot=0
    )
    for n in (64, 256, 1024):
        r = run_uq(CFG, "cv", params=params, n=n, numerical_error=False,
                   model_form=False, n_boot=0)  # fmt: skip
        print(
            f"  batched  n={n:5d} evals {r.metadata['n_evaluations']:6d}  "
            f"{r.metadata['samples_per_s']:8.1f} /s"
        )
    import dataclasses

    X = r.metadata["samples"][:64]
    t = time.perf_counter()
    for row in X:
        run_cv(dataclasses.replace(CFG, N_D_bulk=row[0], N_D_junction=row[1],
                                   L_transition_um=row[2]), 0.0, -30.0, 31)  # fmt: skip
    dt = time.perf_counter() - t
    print(
        f"  run_cv loop (engine, one device per call, with reach-through refinement): "
        f"{64 / dt:.1f} /s"
    )


def science():
    print("\nBands and total indices, default dosimeter (Saltelli):")
    for obs, n, kw, col in [
        ("v_fd", 256, {}, -1),
        ("cv", 128, {"voltages": [0.0, -2.0, -10.0, -30.0]}, -1),
        ("cce", 128, {"voltages": [-1.0, -3.0, -10.0, -30.0]}, 0),
        ("cce", 128, {"voltages": [-1.0, -3.0, -10.0, -30.0]}, -1),
        ("dark_current", 128, {"voltages": [-10.0, -60.0]}, -1),
        ("transit_time", 128, {"voltages": [-10.0, -50.0]}, -1),
        ("breakdown", 64, {}, -1),
        ("compensation_fluence", 2048, {}, -1),
        ("microdosimetry", 256, {}, 1),
    ]:
        t = time.perf_counter()
        r = run_uq(CFG, obs, n=n, **kw)
        md = r.metadata
        x = r.x[col]
        print(
            f"  {obs} at x={x:g}: nominal {md['nominal'][col]:.5g}, 5-95 % "
            f"[{md['band_low'][col]:.5g}, {md['band_high'][col]:.5g}], "
            f"mesh err {None if md['numerical_error'] is None else f'{md['numerical_error'][col]:.2g}'}"
            f", {md['n_evaluations']} evals in {time.perf_counter() - t:.1f} s"
        )
        rows = indices_table(r, col)
        print("     ST " + ", ".join(f"{q['factor']} {q['ST']:.3f} [{q['ST_lo']:.3f}, {q['ST_hi']:.3f}]"
                                     for q in rows))  # fmt: skip
        if md["model_form"]:
            print(f"     model form: {md['model_form']}")


if __name__ == "__main__":
    t0 = time.perf_counter()
    ishigami_table()
    throughput()
    batching_gain()
    if "--quick" not in sys.argv:
        science()
    print(f"\ntotal {time.perf_counter() - t0:.1f} s")
