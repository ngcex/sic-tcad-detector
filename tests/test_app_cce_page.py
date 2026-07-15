"""CCE page AppTest coverage: 1D guard, Run->cache->render->download.

Mirrors tests/test_app_cv_page.py, substituting cce for cv. Each test
monkeypatches petringa.run_cce as a MODULE ATTRIBUTE BEFORE calling
at.run() (the seam proven in tests/test_app_run_mockability.py / 39-01),
so the real devsim solve never executes. AppTest 1.55 has no plotly_chart
/ download_button accessor, so assertions are limited to at.exception,
at.session_state, at.button, at.warning, and at.info.
"""

from __future__ import annotations

import numpy as np
from streamlit.testing.v1 import AppTest

import petringa
from petringa import DeviceConfig, SimResult


def _fake_run_cce(cfg, **kwargs):
    return SimResult(
        config=cfg,
        sim_type="cce",
        x=np.array([-10.0, -100.0]),
        y=np.array([0.85, 0.98]),
        metadata={
            "I_collected": np.array([1e-6, 2e-6]),
            "I_generated": 5e-6,
        },
        mesh=None,
    )


def _run_cce_page():
    from app.workflows.cce import render

    render()


def test_run_caches_result(monkeypatch):
    monkeypatch.setattr(petringa, "run_cce", _fake_run_cce)

    at = AppTest.from_function(_run_cce_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on Run click: {at.exception}"
    result = at.session_state["cce_result"]
    assert result.sim_type == "cce"
    assert len(result.x) == 2  # proves the fake ran, not a real devsim solve


def test_configurable_bias_range_passed_to_facade(monkeypatch):
    captured_kwargs = {}

    def _capturing_run_cce(cfg, **kwargs):
        captured_kwargs.update(kwargs)
        return _fake_run_cce(cfg, **kwargs)

    monkeypatch.setattr(petringa, "run_cce", _capturing_run_cce)

    at = AppTest.from_function(_run_cce_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.number_input(key="cce_v_start").set_value(-1.0)
    at.number_input(key="cce_v_stop").set_value(-20.0)
    at.run()

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on Run click: {at.exception}"
    assert captured_kwargs["v_start"] == -1.0
    assert captured_kwargs["v_stop"] == -20.0


def test_bias_inputs_default_to_facade_defaults(monkeypatch):
    captured_kwargs = {}

    def _capturing_run_cce(cfg, **kwargs):
        captured_kwargs.update(kwargs)
        return _fake_run_cce(cfg, **kwargs)

    monkeypatch.setattr(petringa, "run_cce", _capturing_run_cce)

    at = AppTest.from_function(_run_cce_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert at.exception == []
    assert captured_kwargs["v_start"] == -10.0
    assert captured_kwargs["v_stop"] == -40.0


def test_2d_config_warns_and_skips(monkeypatch):
    monkeypatch.setattr(petringa, "run_cce", _fake_run_cce)

    at = AppTest.from_function(_run_cce_page)
    at.session_state["device_config"] = DeviceConfig(half_width_um=50.0)
    at.run()

    assert at.exception == []
    assert len(at.warning) >= 1
    assert "cce_result" not in at.session_state


def test_empty_state_guard(monkeypatch):
    monkeypatch.setattr(petringa, "run_cce", _fake_run_cce)

    at = AppTest.from_function(_run_cce_page)
    # session_state is empty by default (do not pre-seed device_config).
    at.run()

    assert at.exception == []
    info_texts = [el.value for el in at.info]
    assert "Configure a device in the sidebar to begin." in info_texts


def test_solver_convergence_failure_shows_error_not_crash(monkeypatch):
    def _raise_run_cce(cfg, **kwargs):
        raise RuntimeError(
            "ramp_bias: failed to converge at V=60.542V: Convergence failure!"
        )

    monkeypatch.setattr(petringa, "run_cce", _raise_run_cce)

    at = AppTest.from_function(_run_cce_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert (
        at.exception == []
    ), f"page crashed instead of showing st.error: {at.exception}"
    assert len(at.error) >= 1
    assert "cce_result" not in at.session_state
