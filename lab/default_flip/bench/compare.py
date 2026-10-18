"""Reproducibility check and old -> new impact table of the default flip.

usage: uv run python lab/default_flip/bench/compare.py PRE.json POST.json
Checks: post.default == pre.measured, post.nplus == pre.default, post.legacy == pre.legacy,
post.measured == post.default (alias); core defaults (direct core calls) == the new config.
"""

from __future__ import annotations

import json
import sys

import numpy as np

SKIP = {"config", "run_all_rows", "validated_domain_reason"}
TIMING = ("t_", "run_all_wall_s", "core2d_build_s", "engine2d_point_s")


def flat(d, prefix=""):
    out = {}
    for k, v in d.items():
        if k in SKIP or k.startswith(TIMING):
            continue
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(flat(v, key + "."))
        else:
            out[key] = v
    return out


def worst(a, b):
    fa, fb = flat(a), flat(b)
    keys = sorted(set(fa) & set(fb))
    w, wk, bad = 0.0, None, []
    for k in keys:
        x, y = fa[k], fb[k]
        if isinstance(x, list) and x and isinstance(x[0], (int, float)):
            x, y = np.asarray(x, float), np.asarray(y, float)
            if x.shape != y.shape:
                bad.append(f"{k}: shape {x.shape} vs {y.shape}")
                continue
            den = np.maximum(np.abs(x), 1e-300)
            r = float(np.max(np.abs(x - y) / den)) if x.size else 0.0
        elif (
            isinstance(x, (int, float))
            and isinstance(y, (int, float))
            and x is not None
        ):
            r = abs(x - y) / max(abs(x), 1e-300)
        else:
            if x != y:
                bad.append(f"{k}: {x!r} vs {y!r}")
            continue
        if r > w:
            w, wk = r, k
    missing = sorted(set(fa) ^ set(fb))
    return w, wk, bad, missing


def main():
    pre = json.load(open(sys.argv[1]))
    post = json.load(open(sys.argv[2]))
    for new, old in (
        ("default", "measured"),
        ("nplus", "default"),
        ("legacy", "legacy"),
        ("measured", "measured"),
    ):
        if new in post and old in pre:
            w, wk, bad, missing = worst(pre[old], post[new])
            print(
                f"CHECK post.{new} vs pre.{old}: worst rel {w:.3e} ({wk}); "
                f"non-numeric diffs {bad[:5]}; keys only on one side {missing[:6]}"
            )
    c = post.get("core", {})
    d = post["default"]
    print(
        "CHECK core1d donors == DeviceConfig() donors:",
        (
            np.max(np.abs(np.asarray(c["core1d_donors"]) - np.asarray(d["donors"])) / np.maximum(np.abs(np.asarray(d["donors"])), 1.0))
            if len(c["core1d_donors"]) == len(d["donors"])
            else "node count differs"
        ),
    )
    cc = np.asarray(c["core_cce_vs_bias"])
    dd = np.array([d["cce"][k] for k in ("0", "-2", "-10", "-50")])
    print(
        "CHECK core cce_vs_bias() vs run_cce(DeviceConfig()):", np.max(np.abs(cc - dd))
    )

    o, n = pre["default"], post["default"]
    print("\nIMPACT old default -> new default (4 mm2, 300 K)")
    print(
        f"V_fd analytic {o['vfd_analytic']:.2f} -> {n['vfd_analytic']:.2f} V; gate "
        f"{o['vfd_gate']:.2f} -> {n['vfd_gate']:.2f} V; reach-through {o['reach_through_V']} -> {n['reach_through_V']} V"
    )
    print(
        "eps*A/C um: "
        + "; ".join(
            f"{k} V {o['W_um'][k]:.3f} -> {n['W_um'][k]:.3f}" for k in o["W_um"]
        )
    )
    print(
        "CCE: "
        + "; ".join(
            f"{k} V {o['cce'][k]:.4f} -> {n['cce'][k]:.4f}"
            for k in ("0", "-1", "-2", "-3", "-5", "-10")
        )
    )
    print(f"CCE -50 V {o['cce']['-50']:.6f} -> {n['cce']['-50']:.6f}")
    print(
        "dark pA: "
        + "; ".join(
            f"{k} V {o['dark_pA'][k]:.2f} -> {n['dark_pA'][k]:.2f}"
            for k in o["dark_pA"]
        )
    )
    print(
        "peak E (excl. spikes): "
        + "; ".join(
            f"{k} V {o['peakE'][k][0]:.3g} -> {n['peakE'][k][0]:.3g} V/cm (x {n['peakE'][k][1]:.2f} um)"
            for k in o["peakE"]
        )
    )
    print(
        f"compensation fluence 5.6 MeV p: {o['comp_fluence_5p6MeV']:.3e} -> {n['comp_fluence_5p6MeV']:.3e} p/cm2"
    )
    print(
        f"run_all sweep grid: {o['sweep_valid']} {o['sweep_solvers']} -> {n['sweep_valid']} {n['sweep_solvers']}"
    )
    print(
        f"gate by bias hw 25: {o['gate_by_bias']} -> {n['gate_by_bias']} ({n['gate_solvers']})"
    )
    print(
        f"timings: run_all (all sections) {o['run_all_wall_s']:.1f} -> {n['run_all_wall_s']:.1f} s; "
        f"sweep grid {o['t_sweep_grid_s']:.2f} -> {n['t_sweep_grid_s']:.2f} s; run_cv 241 pts "
        f"{o['t_cv_s']:.3f} -> {n['t_cv_s']:.3f} s ({n['cv_solver']}); CCE solver {n['cce_solver']}"
    )
    co, cn = pre["core"], post["core"]
    print(
        f"2D default mesh (hw 50): {co['core2d_nodes']} -> {cn['core2d_nodes']} nodes, build "
        f"{co['core2d_build_s']:.2f} -> {cn['core2d_build_s']:.2f} s; engine2d point (hw 50, 40 V) "
        f"{co['engine2d_point_s']:.2f} -> {cn['engine2d_point_s']:.2f} s {cn['engine2d_point_solver']}, "
        f"CCE {co['engine2d_point_cce']:.6f} -> {cn['engine2d_point_cce']:.6f}"
    )


if __name__ == "__main__":
    main()
