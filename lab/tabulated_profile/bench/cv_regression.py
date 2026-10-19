"""Tabulated measured-C-V preset vs the digitised Photons Fig. 4b curve, through run_cv
(batched engine), eps*A/C = eps*A/C(dQ/dV), A = 4 mm2.

Gates: knee markers (-3.25..-0.25 V) rms <= 0.05 um; plateau (-5..-59 V, 250 mV grid) rms
<= 0.02 um; labels 0 / -10 / -30 V. Also: the parametric default on the same grid, the
inversion's robustness (the lam 0.3 and lam 3 knot sets of lab/cv_gap/notes/knots.json run
through the same path), and the sensitive-epi-end depth of each knot set.

usage: uv run python lab/tabulated_profile/bench/cv_regression.py
"""

from __future__ import annotations

import json

import numpy as np

import etna
from etna.api.simulation import run_cv
from etna.core.batched.material import EPS
from etna.core.device import (
    _MEASURED_CV_KNOT_DEPTHS_UM,
    donor_table_spec,
    table_sensitive_end_cm,
    table_sensitive_end_doping,
)

ROOT = "lab/cv_gap/notes/"


def digitised():
    band = np.genfromtxt(ROOT + "photons_w.csv", delimiter=",", names=True)
    mk = np.genfromtxt(ROOT + "photons_w_markers.csv", delimiter=",", names=True)
    plateau = [
        (float(v), float(w))
        for v, w, wt in zip(band["V"], band["W_um"], band["weight"])
        if np.isfinite(w)
        and wt >= 10
        and v <= -5.0
        and v not in (-10.0, -30.0, -9.75, -10.25, -10.5, -9.5, -29.75, -30.25)
    ]
    knee = [
        (round(float(v) * 4) / 4, float(w))
        for v, w in zip(mk["V"], mk["W_um"])
        if -3.3 < v < -0.1
    ]
    labels = [(0.0, 1.68), (-10.0, 9.52), (-30.0, 9.73)]
    return plateau, knee, labels


def model_W(config, V):
    r = run_cv(config, v_start=float(V[0]), v_stop=float(V[-1]), n_points=len(V))
    assert r.metadata["solver"] == "batched", r.metadata
    assert np.allclose(r.x, V)
    return EPS * config.area_cm2 / r.y * 1e4, r


def report(tag, cfg, plateau, knee, labels):
    V = np.round(np.arange(0.0, -60.001, -0.25), 2)
    W, r = model_W(cfg, V)
    at = lambda v: W[np.argmin(np.abs(V - v))]
    res_p = np.array([w - at(v) for v, w in plateau])
    res_k = np.array([w - at(v) for v, w in knee])
    slope = np.gradient(W, V)
    sig_k = np.sqrt(
        0.03**2
        + (0.04 * np.array([slope[np.argmin(np.abs(V - v))] for v, _ in knee])) ** 2
    )
    rms_p, rms_k = np.sqrt(np.mean(res_p**2)), np.sqrt(np.mean(res_k**2))
    chi2 = float(
        np.sum((res_p / 0.02) ** 2)
        + np.sum((res_k / sig_k) ** 2)
        + sum(((w - at(v)) / (0.02 if v == 0 else 0.01)) ** 2 for v, w in labels)
    )
    print(
        f"RES [{tag}] plateau N={len(plateau)} rms {rms_p:.4f} um max|res| {np.abs(res_p).max():.4f} "
        f"(gate <= 0.02) | knee N={len(knee)} rms {rms_k:.4f} um (gate <= 0.05) max|res|/sigma "
        f"{np.max(np.abs(res_k) / sig_k):.2f}, within 1 sigma {int(np.sum(np.abs(res_k) <= sig_k))}/{len(knee)}"
    )
    print(
        f"RES [{tag}] labels 0/-10/-30: model {at(0):.4f} / {at(-10):.4f} / {at(-30):.4f} "
        f"(res {1.68 - at(0):+.4f} / {9.52 - at(-10):+.4f} / {9.73 - at(-30):+.4f}); "
        f"W(-2) {at(-2):.3f} (7.18) W(-5) {at(-5):.3f} (9.27) W(-58) {at(-58):.3f} (9.89); "
        f"chi2 {chi2:.1f} / {len(plateau) + len(knee) + 3}; reach-through {r.metadata['reach_through_V']} V"
    )
    return W, V


def main():
    plateau, knee, labels = digitised()
    cfg = etna.measured_cv_tabulated_config(area_cm2=0.04)
    W_tab, V = report("tabulated preset", cfg, plateau, knee, labels)
    for v, w in knee:
        print(
            f"    knee V {v:6.2f} data {w:.3f} model {W_tab[np.argmin(np.abs(V - v))]:.3f}"
        )
    report(
        "parametric default", etna.DeviceConfig(area_cm2=0.04), plateau, knee, labels
    )
    # robustness of the inversion: the other smoothings of knots.json through the same path
    with open(ROOT + "knots.json") as f:
        knots = json.load(f)
    for key in (
        "knots lam 0.3 Qf None",
        "knots lam 3.0 Qf None",
        "knots lam 1.0 Qf 0.0",
    ):
        table = tuple(
            (d, 10.0**l)
            for d, l in zip(_MEASURED_CV_KNOT_DEPTHS_UM, knots[key]["logN_knots"])
        )
        c = etna.measured_cv_tabulated_config(area_cm2=0.04, donor_table=table)
        report(key, c, plateau, knee, labels)
        spec = donor_table_spec(10e-4, table)
        print(
            f"RES [{key}] sensitive-epi end {table_sensitive_end_cm(spec) * 1e4:.3f} um, N_D there "
            f"{table_sensitive_end_doping(spec):.3g}, min {10 ** min(knots[key]['logN_knots']):.3g}"
        )
    spec = donor_table_spec(10e-4, cfg.donor_table)
    print(
        f"RES [tabulated preset] sensitive-epi end {table_sensitive_end_cm(spec) * 1e4:.3f} um, "
        f"N_D there {table_sensitive_end_doping(spec):.3g}"
    )


if __name__ == "__main__":
    main()
