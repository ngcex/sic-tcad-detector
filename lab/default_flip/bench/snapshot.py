"""Snapshot of every number the default flip moves, per named config, to JSON.

Run at the pre-flip HEAD (labels default / measured / legacy) and after the flip
(default / nplus / measured / legacy); compare.py checks new default == old measured,
nplus == old default, legacy == legacy, and prints the old -> new impact table.

usage: uv run python lab/default_flip/bench/snapshot.py OUT.json [label ...]
"""
# NOTE 2026-10-21: the defaults changed (NIEL convention sic-damage, damage model
# "auto" = calibrated kinetic traps). The numbers logged from this bench predate the flip
# (si-neq); pass niel_convention="si-neq" / niel_table=hardness_table("proton", "si-neq")
# and radiation_traps=False or True explicitly to reproduce them.


from __future__ import annotations

import json
import sys
import time

import numpy as np

import etna
from etna.api.results import peak_field_1d
from etna.api.run_all import _SWEEP_BIASES, _SWEEP_HALF_WIDTHS, SECTION_NAMES
from etna.api.simulation import (
    compensation_fluence,
    run_cce,
    run_cv,
    run_dark_current,
    run_field,
    validated_domain_reason,
)
from etna.core.batched.material import EPS
from etna.core.optimization import V_FD_MARGIN, microdosimetric_sweep
from etna.export.labels import analytic_vfd


def configs():
    out = {
        "default": etna.DeviceConfig(area_cm2=0.04),
        "legacy": etna.legacy_config(area_cm2=0.04),
    }
    if hasattr(etna, "measured_cv_config"):
        out["measured"] = etna.measured_cv_config(area_cm2=0.04)
    if hasattr(etna, "nplus_calibrated_config"):
        out["nplus"] = etna.nplus_calibrated_config(area_cm2=0.04)
    return out


def donors_1d(cfg):
    import devsim

    from etna.api.device import build_device
    from etna.core.devsim_reset import reset_devsim_fully

    reset_devsim_fully()
    info = build_device(cfg, device_name="snap1d")
    d = np.array(
        devsim.get_node_model_values(
            device="snap1d", region=info["region_name"], name="Donors"
        )
    )
    x = np.array(
        devsim.get_node_model_values(
            device="snap1d", region=info["region_name"], name="x"
        )
    )
    reset_devsim_fully()
    return x, d


def sweep_kw(cfg):
    kw = dict(
        nplus_thickness_cm=(
            None if cfg.nplus_thickness_um is None else cfg.nplus_thickness_um * 1e-4
        ),
        substrate_thickness_cm=cfg.substrate_thickness_um * 1e-4,
        N_D_junction=cfg.N_D_junction,
        L_transition=cfg.L_transition_um * 1e-4,
        buffer_thickness_cm=(
            None if cfg.buffer_thickness_um is None else cfg.buffer_thickness_um * 1e-4
        ),
        buffer_end_doping=cfg.buffer_end_doping,
    )
    if cfg.nplus_thickness_um is not None:
        kw["N_D_nplus"] = cfg.N_D_nplus
    return kw


def one(cfg, heavy=True):
    r = {"config": {k: v for k, v in cfg.__dict__.items()}}
    r["validated_domain_reason"] = validated_domain_reason(cfg)
    x, d = donors_1d(cfg)
    r["donors_x"] = x.tolist()
    r["donors"] = d.tolist()
    r["vfd_analytic"] = analytic_vfd(cfg)
    r["vfd_gate"] = V_FD_MARGIN * r["vfd_analytic"]
    t = time.perf_counter()
    cv = run_cv(cfg, 0.0, -60.0, 241)
    r["t_cv_s"] = time.perf_counter() - t
    W = EPS * cfg.area_cm2 / cv.y * 1e4
    r["cv_solver"] = cv.metadata.get("solver")
    r["reach_through_V"] = cv.metadata["reach_through_V"]
    r["W_um"] = {
        f"{v:g}": float(W[np.argmin(np.abs(cv.x - v))])
        for v in (0, -2, -5, -10, -30, -60)
    }
    cce = run_cce(cfg, 0.0, -10.0, 11)
    c50 = run_cce(cfg, -50.0, -50.0, 1)
    r["cce_solver"] = cce.metadata.get("solver")
    r["cce"] = {f"{v:g}": float(c) for v, c in zip(cce.x, cce.y)}
    r["cce"]["-50"] = float(c50.y[0])
    dk = run_dark_current(cfg, 0.0, -60.0, 7)
    r["dark_pA"] = {
        f"{v:g}": float(abs(dk.y[np.argmin(np.abs(dk.x - v))]) * 1e12)
        for v in (-10, -30, -60)
    }
    r["peakE"] = {}
    for vb in (-30.0, -60.0):
        f = run_field(cfg, vb)
        pk = peak_field_1d(f.x, f.y, cfg)
        r["peakE"][f"{vb:g}"] = [float(pk[0]), float(pk[1])]
    r["comp_fluence_5p6MeV"] = compensation_fluence(cfg, 5.6)
    if cfg.nplus_thickness_um is None:
        return r
    kw = sweep_kw(cfg)
    t = time.perf_counter()
    df = microdosimetric_sweep(
        half_widths_um=_SWEEP_HALF_WIDTHS,
        epi_thicknesses_cm=(10e-4,),
        N_D_bulks=(cfg.N_D_bulk,),
        V_biases=_SWEEP_BIASES,
        solver="auto",
        **kw,
    )
    r["t_sweep_grid_s"] = time.perf_counter() - t
    r["sweep_valid"] = f"{int(df['is_valid'].sum())}/{len(df)}"
    r["sweep_solvers"] = sorted(set(df["solver"]))
    r["sweep_center_cce"] = [
        float(df["center_cce"].min()),
        float(df["center_cce"].max()),
    ]
    df2 = microdosimetric_sweep(
        half_widths_um=(25.0,),
        epi_thicknesses_cm=(10e-4,),
        N_D_bulks=(cfg.N_D_bulk,),
        V_biases=(5, 10, 15, 20, 30, 40, 60),
        solver="auto",
        **kw,
    )
    r["gate_by_bias"] = {
        int(q.V_bias): (
            "valid"
            if q.is_valid
            else ("shallow" if not q.is_fully_depleted else "lowCCE")
        )
        for q in df2.itertuples()
    }
    r["gate_solvers"] = sorted(set(df2["solver"]))
    if heavy:
        t = time.perf_counter()
        rep = etna.run_all(cfg, sections=SECTION_NAMES)
        r["run_all_wall_s"] = time.perf_counter() - t
        r["run_all_ok"] = bool(rep.ok)
        r["run_all_rows"] = [str(x) for x in rep.rows()]
    return r


def core_defaults():
    """Direct core calls with no device keywords (the core defaults)."""
    import devsim

    from etna.core.charge_collection import cce_vs_bias
    from etna.core.device import create_sic_device
    from etna.core.device2d import create_sic_2d_device
    from etna.core.devsim_reset import reset_devsim_fully

    out = {}
    reset_devsim_fully()
    info = create_sic_device(device_name="core1d", doping_profile="graded")
    out["core1d_donors"] = list(
        devsim.get_node_model_values(
            device="core1d", region=info["region_name"], name="Donors"
        )
    )
    reset_devsim_fully()
    t = time.perf_counter()
    info = create_sic_2d_device(device_name="core2d", half_width_um=50.0)
    out["core2d_build_s"] = time.perf_counter() - t
    out["core2d_nodes"] = len(
        devsim.get_node_model_values(
            device="core2d", region=info["region_name"], name="x"
        )
    )
    reset_devsim_fully()
    res = cce_vs_bias(np.array([0.0, -2.0, -10.0, -50.0]))
    out["core_cce_vs_bias"] = [float(c) for c in res["cce_values"]]
    reset_devsim_fully()
    # 2D engine time on the default (one lateral-scan point, half-width 50 um, 40 V)
    t = time.perf_counter()
    df = microdosimetric_sweep(
        half_widths_um=(50.0,),
        epi_thicknesses_cm=(10e-4,),
        N_D_bulks=(None,) if False else (etna.DeviceConfig().N_D_bulk,),
        V_biases=(40,),
        solver="auto",
        N_D_junction=etna.DeviceConfig().N_D_junction,
        L_transition=etna.DeviceConfig().L_transition_um * 1e-4,
        **(
            {}
            if etna.DeviceConfig().buffer_thickness_um is None
            else dict(
                buffer_thickness_cm=etna.DeviceConfig().buffer_thickness_um * 1e-4,
                buffer_end_doping=etna.DeviceConfig().buffer_end_doping,
            )
        ),
    )
    out["engine2d_point_s"] = time.perf_counter() - t
    out["engine2d_point_solver"] = list(df["solver"])
    out["engine2d_point_cce"] = float(df["center_cce"].iloc[0])
    return out


def main():
    out_path = sys.argv[1]
    labels = sys.argv[2:]
    cfgs = configs()
    res = {}
    for lab, cfg in cfgs.items():
        if labels and lab not in labels and lab != "core":
            continue
        print("==", lab, flush=True)
        res[lab] = one(cfg, heavy=lab == "default")
    if not labels or "core" in labels:
        print("== core", flush=True)
        res["core"] = core_defaults()
    with open(out_path, "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    print("wrote", out_path)


if __name__ == "__main__":
    main()
