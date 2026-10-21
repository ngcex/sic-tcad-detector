"""Bench of etna.design (issue #10): evaluator throughput, BO vs random, robust mode,
alternative physics. Run: uv run python lab/design_opt/bench/design_bench.py [section]
Sections: throughput, analytic, dosimeter, robust, micro (default: all but micro).
Numbers go to lab/design_opt/notes/LOG.md by hand."""

from __future__ import annotations

import sys
import time

import numpy as np

import etna
from etna import DeviceConfig
from etna.design import DesignSettings, evaluate_configs, run_design
from etna.design import testproblems as tp
from etna.design.bo import optimize


def throughput():
    rng = np.random.default_rng(0)
    nb = np.geomspace(3e13, 1e15, 64)
    # warm the numba caches
    evaluate_configs([DeviceConfig()], metrics=("cce_fluence", "V_bd"))
    print("== cheap metrics, per device (warm), shared mesh (epi 10 um)")
    for B in (1, 8, 64):
        cfgs = [DeviceConfig(N_D_bulk=n) for n in nb[:B]]
        t = time.perf_counter()
        evaluate_configs(cfgs)
        dt = time.perf_counter() - t
        print(f"  B={B:3d}: {dt * 1e3 / B:7.2f} ms/device  ({B / dt:6.1f} devices/s)")
    print("== B=64 with continuous epi (one mesh group per device) vs 0.25 um lattice")
    epi_c = rng.uniform(8.0, 12.0, 64)
    epi_l = np.round(epi_c / 0.25) * 0.25
    for label, epi in (("continuous", epi_c), ("lattice", epi_l)):
        cfgs = [
            DeviceConfig(N_D_bulk=n, epi_thickness_um=float(e)) for n, e in zip(nb, epi)
        ]
        t = time.perf_counter()
        m = evaluate_configs(cfgs)
        dt = time.perf_counter() - t
        print(
            f"  {label:10s}: {len(set(m['group']))} groups, {dt * 1e3 / 64:7.2f} ms/device"
        )
    print("== expensive metrics, per device")
    cfgs = [DeviceConfig(N_D_bulk=n) for n in nb[:8]]
    for met in ("cce_fluence", "V_bd"):
        t0 = time.perf_counter()
        evaluate_configs(cfgs)
        t1 = time.perf_counter()
        evaluate_configs(cfgs, metrics=(met,))
        t2 = time.perf_counter()
        print(f"  {met}: +{((t2 - t1) - (t1 - t0)) * 1e3 / 8:.1f} ms/device (B=8)")


def analytic():
    print("== analytic problems: BO vs Sobol (same budget), 5 seeds")
    for name, pr, kw, score in (
        ("branin", tp.branin(), dict(n_init=8, n_iter=14, batch_size=2),
         lambda h: h.F[h.best_index(), 0] - tp.BRANIN_OPTIMUM),
        ("constrained", tp.constrained_quadratic(), dict(n_init=8, n_iter=8, batch_size=2),
         lambda h: h.F[h.best_index(), 0] - tp.CONSTRAINED_OPTIMUM[2]),
        ("zdt1 HV ratio", tp.zdt1(3), dict(n_init=10, n_iter=10, batch_size=4),
         lambda h: h.hypervolume((1.1, 1.1)) / tp.zdt1_true_hypervolume()),
    ):  # fmt: skip
        for m in ("bo", "random"):
            vals, times = [], []
            for seed in range(5):
                t = time.perf_counter()
                h = optimize(pr, seed=seed, method=m, **kw)
                times.append(time.perf_counter() - t)
                vals.append(score(h))
            print(f"  {name:14s} {m:6s} n={h.n}: median {np.median(vals):.4g} "
                  f"[{min(vals):.4g}, {max(vals):.4g}]  {np.mean(times):.2f} s/run")  # fmt: skip


def _front(r, label):
    md = r.metadata
    print(f"-- {label}: {md['n_evaluations']} designs, {md['n_feasible']} feasible, "
          f"{md['n_failed']} failed, eval {md['eval_time_s']:.1f} s, model "
          f"{md['model_time_s']:.1f} s, {md['throughput_designs_per_s'] or 0:.2f} designs/s, "
          f"HV trace {np.round(md['hypervolume_trace'], 3)}")  # fmt: skip
    for p in md["pareto"]:
        v, m = p["variables"], p["metrics"]
        print(
            f"   epi {v['epi_thickness_um']:6.2f} Nb {v['N_D_bulk']:.3g} Nj "
            f"{v['N_D_junction']:.3g} L {v['L_transition_um']:.2f} buf "
            f"{v['buffer_thickness_um']:.2f} | S {m['sensitivity_nC_per_Gy']:6.2f} nC/Gy "
            f"CCE(phi) {m['cce_fluence']:.4f} gate {m['vfd_gate']:.2f} "
            f"I {m['dark_current_A'] * 1e12:.2f} pA"
        )
    for k, a in md["alt_physics"].items():
        if "metrics" in a:
            cf = np.asarray(a["metrics"].get("cce_fluence", []))
            print(f"   alt [{k}]: feasible {a['feasible_fraction']:.2f}, non-dominated "
                  f"{a['non_dominated_fraction']:.2f}, cce_fluence {np.round(cf, 4)}")  # fmt: skip
    if md["breakdown"] is not None:
        print(f"   V_bd {np.round(md['breakdown']['V_bd'])}")


def dosimeter():
    print("== default dosimeter problem (sensitivity vs CCE at 1e14 p/cm2, -10 V)")
    for m in ("bo", "random"):
        hv = []
        for seed in range(3):
            r = run_design(n_iter=6, batch_size=8, seed=seed, method=m, alt_physics=False,
                           breakdown_check=False, cache=False)  # fmt: skip
            md = r.metadata
            # common reference: fixed nadir in the minimised space
            from etna.design.pareto import hypervolume

            F = np.column_stack([-np.asarray(r.x), -np.asarray(r.y)])
            hv.append(hypervolume(F, [0.0, -0.8]) if F.size else 0.0)
        print(
            f"  {m:6s}: n={md['n_evaluations']} HV(ref S=0, CCE=0.8) per seed {np.round(hv, 3)}"
        )
    r = run_design(n_iter=6, batch_size=8, seed=0)
    _front(r, "BO seed 0 with alt physics + breakdown")
    r = run_design(n_iter=6, batch_size=8, seed=0, V_work=-50.0)
    _front(r, "BO seed 0 at -50 V")


def robust():
    print("== robust (worst case over the DESIGN-1 tolerances) vs nominal")
    for rb in (None, "worst", "p5"):
        t = time.perf_counter()
        r = run_design(n_iter=5, batch_size=8, seed=0, robust=rb, alt_physics=False)
        print(f"   [{rb}] wall {time.perf_counter() - t:.1f} s")
        _front(r, f"robust={rb}")


def micro():
    print("== microdosimeter cell (fv2d), -30 V")
    t = time.perf_counter()
    r = etna.run_design(
        problem="microdosimeter", n_init=8, n_iter=4, batch_size=4, seed=0
    )
    md = r.metadata
    print(f"   {time.perf_counter() - t:.1f} s, {md['n_evaluations']} designs, "
          f"{md['throughput_designs_per_s'] or 0:.2f} designs/s")  # fmt: skip
    for p in md["pareto"]:
        print("  ", p["variables"], {k: round(p["metrics"][k], 4) for k in
              ("efficient_fraction", "cce_min", "cell_area_um2", "scan_spacing_um")})  # fmt: skip


if __name__ == "__main__":
    which = sys.argv[1:] or ["throughput", "analytic", "dosimeter", "robust"]
    for w in which:
        globals()[w]()
        sys.stdout.flush()
    _ = DesignSettings
