"""Field map page AppTest coverage: 1D pre-check, Run->cache->render->download.

Each test monkeypatches petringa.run_field as a MODULE ATTRIBUTE (the seam
proven in tests/test_app_run_mockability.py / 39-01) BEFORE calling at.run(),
so the real (expensive) devsim solve never executes. AppTest 1.55 has no
plotly_chart / download_button accessor (39-RESEARCH.md Validation), so
assertions are limited to at.exception, at.session_state, at.button,
at.warning, and at.info.

The 2D guard here is CRITICAL and subtly different from the C-V/CCE pages:
run_field does NOT raise for a 2D config — it silently returns empty x/y
arrays — so the page's guard MUST be a pre-check (before calling run_field),
never a try/except around the call. test_2d_config_warns_and_skips guards
against that trap directly by asserting "field_result" never lands in
session_state for a 2D config.
"""

from __future__ import annotations

import numpy as np
from streamlit.testing.v1 import AppTest

import petringa
from petringa import DeviceConfig, SimResult


def _fake_run_field(cfg, **kwargs):
    return SimResult(
        config=cfg,
        sim_type="field",
        x=np.array([0.0, 1.0, 2.0]),
        y=np.array([1e5, 8e4, 5e4]),
        metadata={
            "bias_V": -100.0,
            "potential": np.array([0.0, -20.0, -45.0]),
            "net_doping": np.array([3e15, 2e15, 9e13]),
        },
        mesh=None,
    )


def _run_field_page():
    from app.workflows.field_map import render

    render()


def test_run_caches_field_result(monkeypatch):
    monkeypatch.setattr(petringa, "run_field", _fake_run_field)

    at = AppTest.from_function(_run_field_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on Run click: {at.exception}"
    result = at.session_state["field_result"]
    assert result.sim_type == "field"
    assert len(result.x) == 3  # proves the fake ran, not a real devsim solve


def test_2d_config_warns_and_skips(monkeypatch):
    monkeypatch.setattr(petringa, "run_field", _fake_run_field)

    at = AppTest.from_function(_run_field_page)
    at.session_state["device_config"] = DeviceConfig(half_width_um=50.0)
    at.run()

    assert at.exception == []
    assert len(at.warning) >= 1
    assert "field_result" not in at.session_state


def test_empty_state_guard(monkeypatch):
    monkeypatch.setattr(petringa, "run_field", _fake_run_field)

    at = AppTest.from_function(_run_field_page)
    # session_state is empty by default (do not pre-seed device_config).
    at.run()

    assert at.exception == []
    info_texts = [el.value for el in at.info]
    assert "Configure a device in the sidebar to begin." in info_texts


def test_solver_convergence_failure_shows_error_not_crash(monkeypatch):
    def _raise_run_field(cfg, **kwargs):
        raise RuntimeError(
            "ramp_bias: failed to converge at V=66.000V: Convergence failure!"
        )

    monkeypatch.setattr(petringa, "run_field", _raise_run_field)

    at = AppTest.from_function(_run_field_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert (
        at.exception == []
    ), f"page crashed instead of showing st.error: {at.exception}"
    assert len(at.error) >= 1
    assert "field_result" not in at.session_state
