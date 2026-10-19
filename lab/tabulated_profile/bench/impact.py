"""Impact table: parametric default stack vs the opt-in tabulated measured-C-V preset
(same 4 mm2 area, 300 K).

V_fd (analytic, reach-through), the 1.2 V_fd gate, eps*A/C, CCE 0..-10 V and -50 V, dark
current at -10/-30/-60 V, peak field (excluding the step spikes), compensation fluence, the
run_all sweep grid validity and the microdosimetric gate by bias. Prints RES lines.

usage: uv run python lab/tabulated_profile/bench/impact.py
"""

from __future__ import annotations

import time

import numpy as np

import etna
from etna.api.results import peak_field_1d
from etna.api.run_all import _SWEEP_BIASES, _SWEEP_HALF_WIDTHS
from etna.api.simulation import (
    compensation_fluence,
    run_cce,
    run_cv,
    run_dark_current,
    run_field,
)
from etna.core.batched.material import EPS
from etna.core.optimization import V_FD_MARGIN, microdosimetric_sweep
from etna.export.labels import analytic_vfd


def one(label, cfg):
    print(f"RES ===== {label}")
    vfd = analytic_vfd(cfg)
    cv = run_cv(cfg, 0.0, -60.0, 241)
    W = EPS * cfg.area_cm2 / cv.y * 1e4
    rt = cv.metadata["reach_through_V"]
    print(
        f"RES V_fd analytic {vfd:.3f} V; 1.2 V_fd gate {V_FD_MARGIN * vfd:.3f} V; "
        f"reach-through {rt} V ({cv.metadata['reach_through_method']}); solver {cv.metadata['solver']}"
    )
    at = lambda v: W[np.argmin(np.abs(cv.x - v))]
    print(
        "RES eps*A/C um: "
        + " ".join(f"{v:g}V={at(v):.3f}" for v in (0, -2, -5, -10, -30, -60))
    )
    cce = run_cce(cfg, 0.0, -10.0, 11)
    c50 = run_cce(cfg, -50.0, -50.0, 1)
    print(
        "RES CCE: "
        + " ".join(f"{v:g}V={c:.4f}" for v, c in zip(cce.x, cce.y))
        + f" | -50V={c50.y[0]:.6f} ({cce.metadata['solver']})"
    )
    t = time.time()
    dk = run_dark_current(cfg, 0.0, -60.0, 7)
    atd = lambda v: abs(dk.y[np.argmin(np.abs(dk.x - v))]) * 1e12
    print(
        f"RES dark current (N_t {cfg.N_t:.7g}): -10V {atd(-10):.3f} pA, -30V {atd(-30):.3f} pA, "
        f"-60V {atd(-60):.3f} pA on {cfg.area_cm2 * 100:g} mm2 ({time.time() - t:.1f} s)"
    )
    for vb in (-30.0, -60.0):
        f = run_field(cfg, vb)
        pk = peak_field_1d(f.x, f.y, cfg)
        print(
            f"RES peak |E| at {vb:g} V: {pk[0]:.4g} V/cm at x = {pk[1]:.3f} um (step spikes "
            f"excluded: {pk[2]}); |E| at the junction node {np.max(f.y):.3g}"
        )
    nd_min = (
        min(n for _, n in cfg.donor_table)
        if cfg.donor_table is not None
        else min(cfg.N_D_junction, cfg.N_D_bulk)
    )
    print(
        f"RES compensation fluence (5.6 MeV p): {compensation_fluence(cfg, 5.6):.3e} p/cm2 "
        f"(min N_D {nd_min:.3g})"
    )
    kw = dict(
        nplus_thickness_cm=cfg.nplus_thickness_um * 1e-4,
        N_D_nplus=cfg.N_D_nplus,
        substrate_thickness_cm=cfg.substrate_thickness_um * 1e-4,
        N_D_junction=cfg.N_D_junction,
        L_transition=cfg.L_transition_um * 1e-4,
        buffer_thickness_cm=(
            None if cfg.buffer_thickness_um is None else cfg.buffer_thickness_um * 1e-4
        ),
        buffer_end_doping=cfg.buffer_end_doping,
    )
    if cfg.donor_table is not None:
        kw["donor_table"] = cfg.donor_table
    t = time.time()
    df = microdosimetric_sweep(
        half_widths_um=_SWEEP_HALF_WIDTHS,
        epi_thicknesses_cm=(10e-4,),
        N_D_bulks=(cfg.N_D_bulk,),
        V_biases=_SWEEP_BIASES,
        solver="auto",
        **kw,
    )
    print(
        f"RES run_all sweep grid (hw {_SWEEP_HALF_WIDTHS}, V {_SWEEP_BIASES}): valid "
        f"{int(df['is_valid'].sum())}/{len(df)}, fully depleted {int(df['is_fully_depleted'].sum())}, "
        f"V_fd column {df['V_fd'].iloc[0]:.3f} V, solvers {sorted(set(df['solver']))}, "
        f"center CCE {df['center_cce'].min():.5f}..{df['center_cce'].max():.5f} ({time.time() - t:.1f} s)"
    )
    df2 = microdosimetric_sweep(
        half_widths_um=(25.0,),
        epi_thicknesses_cm=(10e-4,),
        N_D_bulks=(cfg.N_D_bulk,),
        V_biases=(5, 10, 15, 20, 30, 40, 60),
        solver="auto",
        **kw,
    )
    print(
        "RES gate by bias (hw 25): "
        + " ".join(
            f"{int(r.V_bias)}V:{'valid' if r.is_valid else ('shallow' if not r.is_fully_depleted else 'lowCCE')}"
            for r in df2.sort_values("V_bias").itertuples()
        )
        + f" solvers {sorted(set(df2['solver']))}"
    )


def main():
    one("PARAMETRIC DEFAULT (DeviceConfig)", etna.DeviceConfig(area_cm2=0.04))
    one(
        "TABULATED PRESET (measured_cv_tabulated_config)",
        etna.measured_cv_tabulated_config(area_cm2=0.04),
    )


if __name__ == "__main__":
    main()
