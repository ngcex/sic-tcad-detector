"""C-V page AppTest coverage: 1D guard, Run->cache->render->download.

Each test monkeypatches etna.run_cv as a MODULE ATTRIBUTE (the seam
proven in tests/test_app_run_mockability.py / 39-01) BEFORE calling at.run(),
so the real (expensive) devsim solve never executes. AppTest 1.55 has no
plotly_chart / download_button accessor (39-RESEARCH.md Validation), so
assertions are limited to at.exception, at.session_state, at.button,
at.warning, and at.info.
"""

from __future__ import annotations

import numpy as np
from streamlit.testing.v1 import AppTest

import etna
from etna import DeviceConfig, SimResult


def _fake_run_cv(cfg, **kwargs):
    return SimResult(
        config=cfg,
        sim_type="cv",
        x=np.array([0.0, -10.0]),
        y=np.array([1e-9, 2e-9]),
        metadata={
            "one_over_C_squared": np.array([1e18, 2.5e17]),
            "depletion_widths": np.array([1e-4, 2e-4]),
            "area_cm2": 1e-4,
        },
        mesh=None,
    )


def _run_cv_page():
    from app.workflows.cv import render

    render()


def test_run_caches_result(monkeypatch):
    monkeypatch.setattr(etna, "run_cv", _fake_run_cv)

    at = AppTest.from_function(_run_cv_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on Run click: {at.exception}"
    result = at.session_state["cv_result"]
    assert result.sim_type == "cv"
    assert len(result.x) == 2  # proves the fake ran, not a real devsim solve


def test_2d_config_warns_and_skips(monkeypatch):
    monkeypatch.setattr(etna, "run_cv", _fake_run_cv)

    at = AppTest.from_function(_run_cv_page)
    at.session_state["device_config"] = DeviceConfig(half_width_um=50.0)
    at.run()

    assert at.exception == []
    assert len(at.warning) >= 1
    assert "cv_result" not in at.session_state


def test_empty_state_guard(monkeypatch):
    monkeypatch.setattr(etna, "run_cv", _fake_run_cv)

    at = AppTest.from_function(_run_cv_page)
    # session_state is empty by default (do not pre-seed device_config).
    at.run()

    assert at.exception == []
    info_texts = [el.value for el in at.info]
    assert "Configure a device in the sidebar to begin." in info_texts


def test_solver_convergence_failure_shows_error_not_crash(monkeypatch):
    def _raise_run_cv(cfg, **kwargs):
        raise RuntimeError(
            "ramp_bias: failed to converge at V=60.000V: Convergence failure!"
        )

    monkeypatch.setattr(etna, "run_cv", _raise_run_cv)

    at = AppTest.from_function(_run_cv_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert (
        at.exception == []
    ), f"page crashed instead of showing st.error: {at.exception}"
    assert len(at.error) >= 1
    assert "cv_result" not in at.session_state
