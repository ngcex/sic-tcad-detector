"""Batched engine path of the fluence sweeps vs their devsim path (same mesh, physics).

cce_vs_fluence / cce_uncertainty_envelope / cv_at_fluence solve the damaged devices
with the numba engine by default. Measured on notebook 12's grids (2026-10-17, 11 + 9*30
fluences, 1e10-5e13 p/cm^2, -20..-100 V, three- and single-defect models): max |dCCE|
3.3e-15, max rel dW / dC 8e-16 (200-point C-V to -100 V).
"""

import numpy as np
import pytest

import etna.core.charge_collection as cc
from etna.core.charge_collection import cce_vs_fluence
from etna.core.cv_analysis import cv_at_fluence
from etna.core.radiation_damage import (
    RadiationDamageParams,
    cce_uncertainty_envelope,
    make_single_defect_params,
)

# default stack (measured-C-V profile since 2026-10-18, N_D_bulk 3.38e13): compensation at
# 1.93e13 p/cm2 (62 MeV); 1.8e13: near compensation; 3e13: NaN. (Previous default,
# N_D_bulk 8.25e13: 0 / 1e11 / 3e12 / 2e13 / 3.5e13 / 5e13.)
FLUENCES = [0.0, 1e11, 3e12, 1e13, 1.8e13, 3e13]


def _same_cce(a, b, tol=1e-9):
    assert np.array_equal(np.isnan(a["cce_values"]), np.isnan(b["cce_values"]))
    assert np.array_equal(a["valid"], b["valid"])
    m = ~np.isnan(a["cce_values"])
    assert np.max(np.abs(a["cce_values"][m] - b["cce_values"][m])) < tol


@pytest.mark.parametrize("V_bias", [-40.0, -100.0])
@pytest.mark.parametrize("single", [False, True])
def test_cce_vs_fluence_engine_matches_devsim(V_bias, single):
    dp = make_single_defect_params() if single else None
    eng = cce_vs_fluence(FLUENCES, V_bias=V_bias, damage_params=dp, solver="engine")
    dev = cce_vs_fluence(FLUENCES, V_bias=V_bias, damage_params=dp, solver="devsim")
    _same_cce(eng, dev)
    assert eng["solver"][:-1] == ["batched"] * (len(FLUENCES) - 1)
    assert dev["solver"][:-1] == ["devsim"] * (len(FLUENCES) - 1)
    assert eng["solver"][-1] == "" and "full donor compensation" in eng["reason"][-1]


def test_cce_vs_fluence_legacy_geometry_matches_devsim():
    kw = dict(nplus_thickness_cm=None, substrate_thickness_cm=1e-4, buffer_thickness_cm=None, buffer_end_doping=None)
    eng = cce_vs_fluence(
        [0.0, 1e12, 1e13], V_bias=-40.0, device_kwargs=kw, solver="engine"
    )
    dev = cce_vs_fluence(
        [0.0, 1e12, 1e13], V_bias=-40.0, device_kwargs=kw, solver="devsim"
    )
    _same_cce(eng, dev)  # asserts NaN pattern, validity and |dCCE| < 1e-9
    assert eng["solver"] == ["batched"] * 3 and dev["solver"] == ["devsim"] * 3
    assert all(eng["valid"])


def test_legacy_geometry_beyond_cross_check_stays_on_devsim():
    kw = dict(nplus_thickness_cm=None, substrate_thickness_cm=1e-4, buffer_thickness_cm=None, buffer_end_doping=None)
    why = cc._fluence_engine_reason(kw, -60.0)
    assert isinstance(why, str) and why
    assert cc._fluence_engine_reason(kw, -40.0) is None
    assert cc._fluence_engine_reason({}, -100.0) is None  # default stack has n+
    with pytest.raises(ValueError, match="legacy geometry"):
        cce_vs_fluence([1e12], V_bias=-60.0, device_kwargs=kw, solver="engine")


def test_out_of_scope_uses_devsim():
    r = cce_vs_fluence([1e12], V_bias=-40.0, device_kwargs={"T": 450.0})
    assert r["solver"] == ["devsim"] and r["valid"][0]
    with pytest.raises(ValueError, match="engine unavailable"):
        cce_vs_fluence([1e12], device_kwargs={"T": 450.0}, solver="engine")


def test_unconverged_engine_points_fall_back_to_devsim(monkeypatch):
    real = cc._cce_fluence_engine

    def half_converged(*args, **kwargs):
        cce, conv = real(*args, **kwargs)
        conv = conv.copy()
        conv[::2] = False
        return cce, conv

    monkeypatch.setattr(cc, "_cce_fluence_engine", half_converged)
    fl = [0.0, 1e12, 1e13]
    r = cce_vs_fluence(fl, V_bias=-40.0)
    assert r["solver"] == ["devsim", "batched", "devsim"]
    _same_cce(r, cce_vs_fluence(fl, V_bias=-40.0, solver="devsim"))


def test_envelope_equals_per_set_sweeps():
    fl = np.array([1e11, 5e12, 4e13])  # default N_D_bulk 3.38e13 since 2026-10-18 (was 1e13)
    env = cce_uncertainty_envelope(fl, V_bias=-40.0)
    nominal = cce_vs_fluence(fl, V_bias=-40.0, damage_params=RadiationDamageParams())
    # batch composition moves the Gummel stop by a few ulp
    np.testing.assert_allclose(env["cce_nominal"], nominal["cce_values"], rtol=0, atol=1e-12)
    ok = np.isfinite(env["cce_min"]) & np.isfinite(env["cce_max"])
    assert ok[:2].all()  # 4e13 is fully compensated for the nominal/2x-removal sets
    assert np.all(env["cce_min"][ok] <= env["cce_nominal"][ok])
    assert np.all(env["cce_nominal"][ok] <= env["cce_max"][ok])


@pytest.mark.parametrize(
    "fluence, V_range",
    [
        (5e12, np.linspace(-100, 0, 41)),  # starts deep, as notebook 12
        (1.2e13, [0.0, -1.0, -5.0, -20.0, -40.0]),  # was 2e13 (compensated on the new default)
    ],
)
def test_cv_at_fluence_engine_matches_devsim(fluence, V_range):
    eng = cv_at_fluence(fluence, V_range, solver="engine")
    dev = cv_at_fluence(fluence, V_range, solver="devsim")
    assert eng["solver"] == "batched" and dev["solver"] == "devsim"
    assert not dev["truncated"]
    np.testing.assert_array_equal(eng["voltages"], dev["voltages"])
    np.testing.assert_allclose(
        eng["depletion_widths"], dev["depletion_widths"], rtol=1e-9
    )
    # C = dQ/dV since 2026-10-17: engine exact linearisation vs devsim +-1 mV central
    # difference of the hole charge (independent; measured <= 3e-8). Was rtol 1e-9
    # when both sides were eps*A/W_edge of the same W.
    np.testing.assert_allclose(eng["capacitance"], dev["capacitance"], rtol=1e-6)


def test_cv_at_fluence_beyond_phi_crit_still_none():
    assert cv_at_fluence(2e14, [0.0, -10.0]) is None
