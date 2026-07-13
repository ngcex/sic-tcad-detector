"""Radiation Damage page AppTest coverage: persistent banner + Run->cache->
render->download + NaN tolerance.

Each test monkeypatches petringa.run_radiation_damage as a MODULE ATTRIBUTE
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
from streamlit.testing.v1 import AppTest

import petringa
from petringa import DeviceConfig, SimResult


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
    monkeypatch.setattr(petringa, "run_radiation_damage", _fake_run_radiation_damage)

    at = AppTest.from_function(_run_radiation_damage_page)
    # session_state is empty by default (do not pre-seed device_config).
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    warning_texts = [el.value for el in at.warning]
    assert any(
        "data-blocked" in (w or "").lower() for w in warning_texts
    ), f"kappa data-blocked banner not found in warnings: {warning_texts}"


def test_empty_state_guard(monkeypatch):
    monkeypatch.setattr(petringa, "run_radiation_damage", _fake_run_radiation_damage)

    at = AppTest.from_function(_run_radiation_damage_page)
    at.run()

    assert at.exception == []
    info_texts = [el.value for el in at.info]
    assert "Configure a device in the sidebar to begin." in info_texts


def test_run_caches_damage_result(monkeypatch):
    monkeypatch.setattr(petringa, "run_radiation_damage", _fake_run_radiation_damage)

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
        petringa, "run_radiation_damage", _fake_run_radiation_damage_with_nan
    )

    at = AppTest.from_function(_run_radiation_damage_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on NaN result: {at.exception}"
    info_texts = [el.value for el in at.info]
    assert any(
        "did not converge" in (i or "") for i in info_texts
    ), f"partial-failure NaN info message not found: {info_texts}"


def test_solver_convergence_failure_shows_error_not_crash(monkeypatch):
    def _raise_run_radiation_damage(cfg, **kwargs):
        raise RuntimeError(
            "ramp_bias: failed to converge at V=-20.000V: Convergence failure!"
        )

    monkeypatch.setattr(petringa, "run_radiation_damage", _raise_run_radiation_damage)

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
    monkeypatch.setattr(petringa, "run_radiation_damage", _fake_run_radiation_damage)

    at = AppTest.from_function(_run_radiation_damage_page)
    at.session_state["device_config"] = DeviceConfig(half_width_um=50.0)
    at.run()

    assert at.exception == [], f"page crashed on 2D config: {at.exception}"
    warning_texts = [el.value for el in at.warning]
    assert any(
        "1D-only" in (w or "") for w in warning_texts
    ), f"1D-only warning not found: {warning_texts}"
