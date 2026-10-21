"""Regression tests of the minor findings of the v8 damage-chain audit
(.planning/AUDIT_v8_MODELS_DAMAGE.md, N1-N4, N9)."""

from __future__ import annotations

import dataclasses
import inspect

import numpy as np
import pytest

from etna.api.device import DeviceConfig


# ------------------------------------------------------------------ N1
def test_semi_insulating_fluence_uses_the_mobility_set():
    from etna.core.radiation_traps import electron_mobility, trap_diagnostics
    from etna.core.trap_calibration import semi_insulating_fluence
    from etna.core.trap_presets import radmanovac_calibrated

    m = radmanovac_calibrated(niel_convention="si-neq")
    mu = electron_mobility(4e13, 300.0)  # c_axis default
    assert mu == pytest.approx(1173.0, abs=1.0)
    assert electron_mobility(4e13, 300.0, "legacy_perp") == pytest.approx(
        947.0, abs=5.0
    )
    F = semi_insulating_fluence(m, 4e13, freq_Hz=1e6, kappa=0.8514)
    # before: the fixed legacy 950 (4.94e12 p/cm2, onsemi 1 MHz, si-neq preset)
    assert F == pytest.approx(
        semi_insulating_fluence(m, 4e13, freq_Hz=1e6, mu_n=mu, kappa=0.8514), rel=1e-12
    )
    assert F == pytest.approx(5.11e12, rel=5e-3)
    d = trap_diagnostics(4e13, 1e12, model=m)
    assert d["resistivity_ohm_cm"] == pytest.approx(
        1.0 / (1.602e-19 * mu * d["n_neutral"]), rel=1e-12
    )


# ------------------------------------------------------------------ N2
def test_calibrate_removal_reports_nan_and_warns(monkeypatch):
    import etna.core.trap_calibration as tc
    from etna.core.trap_presets import literature_kinetic

    def fake(o, m):  # linear in donor_removal; observation "b" unmeasurable
        return np.nan if o.name == "b" else 2.0 + m.donor_removal

    monkeypatch.setattr(tc, "removal_observable", fake)
    obs = [
        tc.RemovalObservation(n, "slope", v, 0.5, 4e13, fluences=(0, 1e12))
        for n, v in (("a", 5.0), ("b", 4.0), ("c", 4.5))
    ]
    cal = tc.calibrate_removal(obs, ("donor_removal",), base_model=literature_kinetic())
    mv = cal.settings["model_values"]
    # before: [4.75, 504.0, 4.75] (the 1e3 penalty reported as a model value), no warning
    assert np.isnan(mv[1]) and mv[0] == pytest.approx(4.75, abs=1e-6)
    assert any("not finite" in w and "'b'" in w for w in cal.warnings)


# ------------------------------------------------------------------ N3
def test_out_of_table_energy_raises_in_core_cce_vs_fluence():
    from etna.api.device import device_kwargs
    from etna.core.charge_collection import cce_vs_fluence

    kw = device_kwargs(DeviceConfig())
    eth = kw.pop("epi_thickness_cm")
    # before: per-point NaN labelled "solver failure: ValueError ..."
    with pytest.raises(ValueError, match="outside|range|table"):
        cce_vs_fluence(
            np.array([0.0, 1e12]), V_bias=-40.0, epi_thickness_cm=eth,
            energy_MeV=2000.0, device_kwargs=kw,
        )  # fmt: skip


# ------------------------------------------------------------------ N4
def test_auto_helpers_default_to_the_run_bias():
    from etna.api import simulation as s

    v_run = inspect.signature(s.run_radiation_damage).parameters["V_bias"].default
    for f in (s.damage_model_of, s.compensation_fluence, s.default_fluence_grid):
        # before: 0 V, so "auto" was decided at another bias than the run's
        assert inspect.signature(f).parameters["V_bias"].default == v_run


def test_run_all_summary_shows_the_fallback_note():
    from etna.api import run_all as ra
    from etna.api.results import SimResult

    res = SimResult(
        config=DeviceConfig(), sim_type="damage", x=np.array([1e12]),
        y=np.array([0.99]),
        metadata={"compensation_fluence": 4.36e12, "damage_model": "carrier_removal",
                  "damage_model_note": "radiation_traps='auto': legacy model (test)"},
    )  # fmt: skip
    assert "legacy model (test)" in ra._compensation_note(res)


# ------------------------------------------------------------------ N9
def test_a_rescaled_d_breaks_the_per_proton_invariance():
    """The re-fitted sic-damage d reproduces the si-neq per-proton slope (tested in
    test_niel_convention_calibration); a kappa-RESCALED d (3.83 /neq) does not: the
    slope then misses by ~0.6 /p (USAGE "NIEL convention and the calibrated
    removal"). A wrong d in the preset table fails here."""
    from etna.core.radiation_damage import get_hardness_factor, hardness_table
    from etna.core.trap_calibration import removal_observable
    from etna.core.trap_presets import (
        CALIBRATED_SLOW_REMOVAL_BY_CONVENTION,
        radmanovac_calibrated,
    )
    from tests.test_niel_convention_calibration import ONSEMI

    k = {c: get_hardness_factor(252.7, hardness_table("proton", c))
         for c in ("si-neq", "sic-damage")}  # fmt: skip
    d_si = CALIBRATED_SLOW_REMOVAL_BY_CONVENTION["si-neq"][0]
    naive = d_si * k["si-neq"] / k["sic-damage"]
    # 2 fluences, 0-70 V: the shipped d gives the full-setup 4.948 / 4.857 to 1e-3
    obs = dataclasses.replace(ONSEMI, niel_convention="sic-damage",
                              fluences=(0.0, 2.8e12), V_max=70.0)  # fmt: skip
    shipped = radmanovac_calibrated(niel_convention="sic-damage")
    g_naive = removal_observable(obs, dataclasses.replace(shipped, donor_removal=naive))
    # 4.304 /p (rescaled d 3.83) vs 4.857 /p (shipped d 4.69) and 4.948 /p (si-neq)
    assert g_naive == pytest.approx(4.304, abs=0.02)
    assert 4.948 - g_naive > 0.4  # vs |4.948 - 4.857| = 0.09 with the shipped d


def test_auto_fallback_end_to_end_in_run_metadata():
    """N9: the "auto" fallback seen through run_radiation_damage, not damage_model_of:
    legacy geometry beyond 40 V runs the legacy removal and says why; the helper at
    the same bias agrees."""
    from etna import run_radiation_damage
    from etna.api.device import legacy_config
    from etna.api.simulation import damage_model_of

    cfg = legacy_config()
    r = run_radiation_damage(cfg, fluences=np.array([1e11]), V_bias=-60.0,
                             proton_energy_MeV=62.0)  # fmt: skip
    md = r.metadata
    assert md["damage_model"] == "carrier_removal"
    assert "legacy geometry" in md["damage_model_note"]
    assert np.all(md["valid"])
    model, note = damage_model_of(cfg, 62.0, V_bias=-60.0)
    assert model is None and note == md["damage_model_note"]
