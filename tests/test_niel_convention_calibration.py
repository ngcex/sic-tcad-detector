"""NIEL-convention consistency of the calibrated removal (issues #5, #8; 2026-10-21).

The Radmanovac slopes are per PROTON: the calibrated preset must reproduce the same
per-proton C-V removal slope whichever NIEL convention it is built (and run) in, while
its per-neq d differs (re-fitted, not rescaled). Full refit:
lab/radiation_traps/bench/convention_refit.py (~6 min); here one device, one evaluation
per convention.
"""

from __future__ import annotations

import dataclasses

import pytest

from etna.core.radiation_damage import get_hardness_factor, hardness_table
from etna.core.trap_calibration import RemovalObservation, removal_observable
from etna.core.trap_presets import (
    CALIBRATED_SLOW_REMOVAL_BY_CONVENTION,
    radmanovac_calibrated,
)

ONSEMI = RemovalObservation(
    "onsemi g", "slope", 5.60, 0.81, 4e13, 50.0, 252.7,
    (0.0, 1.4e12, 2.8e12), (10.0, 40.0), 150.0, 5.0, 1e6,
)  # fmt: skip  (3 of the 5 published fluences: the test stays < 10 s)


def test_per_proton_slope_is_convention_independent():
    g = {}
    for conv in ("si-neq", "sic-damage"):
        obs = dataclasses.replace(ONSEMI, niel_convention=conv)
        g[conv] = removal_observable(obs, radmanovac_calibrated(niel_convention=conv))
    # refit values on the 5 fluences (convention_refit.py): 4.948 / 4.858 per p/cm2
    # (the slope is linear in fluence here: 3 fluences give the same to 1e-3)
    assert g["si-neq"] == pytest.approx(4.948, abs=0.01)
    assert g["sic-damage"] == pytest.approx(4.858, abs=0.01)
    assert abs(g["si-neq"] - g["sic-damage"]) < 0.15  # << data sigma 0.81


def test_refit_is_not_a_kappa_rescaling():
    k = {c: get_hardness_factor(252.7, hardness_table("proton", c)) for c in
         ("si-neq", "sic-damage", "sic-1mev")}  # fmt: skip
    d = {c: v[0] for c, v in CALIBRATED_SLOW_REMOVAL_BY_CONVENTION.items()}
    assert k["sic-damage"] == pytest.approx(0.639, abs=1e-3)
    naive = d["si-neq"] * k["si-neq"] / k["sic-damage"]  # 3.83: ignores the Burin share
    assert d["sic-damage"] - naive > 0.5
    # the per-proton direct removal d * kappa still lands in the same decade
    assert d["sic-damage"] * k["sic-damage"] == pytest.approx(3.0, abs=0.1)
    assert d["sic-1mev"] > d["sic-damage"] > d["si-neq"]
