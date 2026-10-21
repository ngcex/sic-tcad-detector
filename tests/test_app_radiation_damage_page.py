"""Radiation Damage page AppTest coverage: persistent banner + Run->cache->
render->download + NaN tolerance.

Each test monkeypatches etna.run_radiation_damage as a MODULE ATTRIBUTE
(the seam proven in tests/test_app_run_mockability.py / 39-01) BEFORE
calling at.run(), so the real (expensive) devsim solve never executes.
Mirrors tests/test_app_field_page.py's structure: a `_run_radiation_damage_page()`
wrapper importing and calling app.workflows.radiation_damage.render, and
assertions limited to at.exception, at.session_state, at.button, at.warning,
at.info, and at.error (AppTest 1.55 has no plotly_chart / download_button
accessor, per 39/40-RESEARCH).
"""

from __future__ import annotations

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

import etna
from etna import DeviceConfig, SimResult


def _fake_run_radiation_damage(cfg, **kwargs):
    """Clean fake: 6-element result, no NaN."""
    fluences = np.geomspace(1e13, 1e16, 6)
    cce = np.array([0.9928, 0.9820, 0.9706, 0.9208, 0.7564, 0.4191])
    return SimResult(
        config=cfg,
        sim_type="damage",
        x=fluences,
        y=cce,
        metadata={"V_bias": -20.0, "energy_MeV": 5.6},
        mesh=None,
    )


def _fake_run_radiation_damage_with_nan(cfg, **kwargs):
    """NaN fake: one NaN mid-array, per 41-01-SPIKE-NOTES.md Check B."""
    fluences = np.geomspace(1e13, 1e16, 6)
    cce = np.array([0.9928, np.nan, 0.9706, 0.9208, 0.7564, 0.4191])
    return SimResult(
        config=cfg,
        sim_type="damage",
        x=fluences,
        y=cce,
        metadata={"V_bias": -20.0, "energy_MeV": 5.6},
        mesh=None,
    )


def _run_radiation_damage_page():
    from app.workflows.radiation_damage import render

    render()


def test_kappa_banner_persistent(monkeypatch):
    monkeypatch.setattr(etna, "run_radiation_damage", _fake_run_radiation_damage)

    at = AppTest.from_function(_run_radiation_damage_page)
    # session_state is empty by default (do not pre-seed device_config).
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    warning_texts = [el.value for el in at.warning]
    # SR-NIEL data since 2026-10-20 (B-07); since 2026-10-21 the SiC neutron NIEL is
    # computed (NJOY): the banner states r_n and the opt-in sic-damage convention
    assert any(
        "r_n = 1.069" in (w or "").lower()
        and "sr-niel" in (w or "").lower()
        and "sic-damage" in (w or "").lower()
        for w in warning_texts
    ), f"kappa NIEL-normalisation banner not found in warnings: {warning_texts}"
    assert not any("data-blocked" in (w or "").lower() for w in warning_texts)


def test_empty_state_guard(monkeypatch):
    monkeypatch.setattr(etna, "run_radiation_damage", _fake_run_radiation_damage)

    at = AppTest.from_function(_run_radiation_damage_page)
    at.run()

    assert at.exception == []
    info_texts = [el.value for el in at.info]
    assert "Configure a device in the sidebar to begin." in info_texts


def test_run_caches_damage_result(monkeypatch):
    monkeypatch.setattr(etna, "run_radiation_damage", _fake_run_radiation_damage)

    at = AppTest.from_function(_run_radiation_damage_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on Run click: {at.exception}"
    result = at.session_state["damage_result"]
    assert result.sim_type == "damage"
    assert len(result.x) == 6  # proves the fake ran, not a real devsim solve


def test_nan_in_result_does_not_crash(monkeypatch):
    monkeypatch.setattr(
        etna, "run_radiation_damage", _fake_run_radiation_damage_with_nan
    )

    at = AppTest.from_function(_run_radiation_damage_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on NaN result: {at.exception}"
    warning_texts = [el.value or "" for el in at.warning]
    assert any(
        "1 of 6 fluence point is not shown" in w and "not converged" in w
        for w in warning_texts
    ), f"excluded-point warning not found: {warning_texts}"


def test_compensated_points_listed_and_defaults_below_compensation(monkeypatch):
    from etna.api.simulation import compensation_fluence

    seen = {}

    def _fake(cfg, fluences, **kwargs):
        seen["fluences"] = fluences
        valid = np.array([True, True, False])
        reason = ["", "", "full donor compensation: ...; not solved"]
        return SimResult(
            config=cfg,
            sim_type="damage",
            x=np.array([1e13, 3e13, 1e14]),
            y=np.array([0.99, 0.98, np.nan]),
            metadata={
                "V_bias": -20.0,
                "energy_MeV": 62,
                "valid": valid,
                "reason": reason,
            },
        )

    monkeypatch.setattr(etna, "run_radiation_damage", _fake)
    at = AppTest.from_function(_run_radiation_damage_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    # the legacy carrier-removal model (the default until 2026-10-21) refuses points
    at.selectbox(key="rad_damage_model").set_value("Legacy carrier removal").run()
    phi_c = compensation_fluence(DeviceConfig(radiation_traps=False), 62)
    assert phi_c == pytest.approx(4.358e12, rel=1e-3)  # sic-damage default
    assert any(f"{phi_c:.2e}" in (i.value or "") for i in at.info)

    at.button[0].click()
    at.run()
    assert at.exception == []
    # Default range: 3 decades ending at 0.7x the compensation fluence.
    assert seen["fluences"].max() < phi_c
    assert seen["fluences"].max() == pytest.approx(0.7 * phi_c, rel=1e-2)
    assert seen["fluences"].min() == pytest.approx(0.7e-3 * phi_c, rel=1e-2)
    assert any(
        "epi fully compensated: the model solves an ill-posed problem there"
        in (w.value or "")
        for w in at.warning
    )


def test_solver_convergence_failure_shows_error_not_crash(monkeypatch):
    def _raise_run_radiation_damage(cfg, **kwargs):
        raise RuntimeError(
            "ramp_bias: failed to converge at V=-20.000V: Convergence failure!"
        )

    monkeypatch.setattr(etna, "run_radiation_damage", _raise_run_radiation_damage)

    at = AppTest.from_function(_run_radiation_damage_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert (
        at.exception == []
    ), f"page crashed instead of showing st.error: {at.exception}"
    assert len(at.error) >= 1
    assert "damage_result" not in at.session_state


def test_2d_config_shows_1d_only_warning(monkeypatch):
    monkeypatch.setattr(etna, "run_radiation_damage", _fake_run_radiation_damage)

    at = AppTest.from_function(_run_radiation_damage_page)
    at.session_state["device_config"] = DeviceConfig(half_width_um=50.0)
    at.run()

    assert at.exception == [], f"page crashed on 2D config: {at.exception}"
    warning_texts = [el.value for el in at.warning]
    assert any(
        "1D-only" in (w or "") for w in warning_texts
    ), f"1D-only warning not found: {warning_texts}"


def test_invalid_fluence_range_is_rejected_before_the_facade(monkeypatch):
    """audit C-08: min 0, negative min or min > max -> st.error, facade not called."""
    calls = []
    monkeypatch.setattr(
        etna, "run_radiation_damage", lambda cfg, **k: calls.append(k) or None
    )
    for fmin, fmax in ((0.0, 1e14), (-1e13, 1e14), (1e15, 1e13)):
        at = AppTest.from_function(_run_radiation_damage_page)
        at.session_state["device_config"] = DeviceConfig()
        at.run()
        by_label = {n.label: n for n in at.number_input}
        by_label["Min fluence (p/cm²)"].set_value(fmin)
        by_label["Max fluence (p/cm²)"].set_value(fmax)
        at.button[0].click()
        at.run()
        assert at.exception == [], (fmin, fmax, at.exception)
        assert any("min ≤ max" in (e.value or "") for e in at.error), (fmin, fmax)
        assert "damage_result" not in at.session_state
    assert calls == []


def test_damage_model_selector(monkeypatch):
    """Re-pinned 2026-10-21 (was the opt-in charged-trap checkbox): the default model is
    the calibrated kinetic preset ("auto", the sidebar config unchanged) with a range
    from 0.01x its compensation fluence to 1e16 neq; "Legacy carrier removal" passes
    radiation_traps="legacy", "Burin Table I traps" passes True."""
    from etna.api.simulation import compensation_fluence
    from etna.core.radiation_damage import get_hardness_factor

    seen = {}

    def _fake(cfg, fluences, **kwargs):
        seen["cfg"], seen["fluences"] = cfg, fluences
        return _fake_run_radiation_damage(cfg)

    monkeypatch.setattr(etna, "run_radiation_damage", _fake)
    at = AppTest.from_function(_run_radiation_damage_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    assert at.selectbox(key="rad_damage_model").value.startswith("Calibrated")
    at.button[0].click()
    at.run()
    assert at.exception == []
    assert seen["cfg"].radiation_traps == "auto"
    phi = compensation_fluence(DeviceConfig(), 62)
    assert phi == pytest.approx(1.802e12, rel=2e-3)
    assert seen["fluences"].min() == pytest.approx(1e-2 * phi, rel=1e-2)
    assert seen["fluences"].max() == pytest.approx(
        1e16 / get_hardness_factor(62), rel=1e-2
    )
    assert any("semi-insulating" in (i.value or "") for i in at.info)
    assert any(f"{phi:.2e}" in (i.value or "") for i in at.info)
    for label, spec in (("Legacy carrier removal", "legacy"),
                        ("Burin Table I traps (steady state)", True)):  # fmt: skip
        at.selectbox(key="rad_damage_model").set_value(label).run()
        at.button[0].click()
        at.run()
        assert at.exception == []
        assert seen["cfg"].radiation_traps == spec


def test_neutron_sic_damage_choices_reach_the_facade(monkeypatch):
    """The opt-in particle / NIEL-convention selectors are passed to
    run_radiation_damage (default runs pass neither: byte-identical calls)."""
    seen = {}

    def fake(cfg, **kwargs):
        seen.update(kwargs)
        r = _fake_run_radiation_damage(cfg, **kwargs)
        r.metadata.update(particle=kwargs.get("particle", "proton"))
        return r

    monkeypatch.setattr(etna, "run_radiation_damage", fake)
    at = AppTest.from_function(_run_radiation_damage_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    at.radio(key="rad_particle").set_value("neutron").run()
    assert at.selectbox(key="rad_niel_convention").value == "sic-damage"  # default
    at.selectbox(key="rad_niel_convention").set_value("si-neq").run()
    assert at.exception == []
    assert any("neutron" in (c.value or "") and "si-neq" in (c.value or "")
               for c in at.caption)
    at.button[0].click()
    at.run()
    assert at.exception == []
    assert seen["particle"] == "neutron"
    assert seen["niel_convention"] == "si-neq"
    assert seen["proton_energy_MeV"] == pytest.approx(1.0)
    assert "damage_result" in at.session_state


def test_default_run_passes_no_niel_keywords(monkeypatch):
    seen = {}

    def fake(cfg, **kwargs):
        seen.update(kwargs)
        return _fake_run_radiation_damage(cfg, **kwargs)

    monkeypatch.setattr(etna, "run_radiation_damage", fake)
    at = AppTest.from_function(_run_radiation_damage_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    at.button[0].click()
    at.run()
    assert at.exception == []
    assert "particle" not in seen and "niel_convention" not in seen
