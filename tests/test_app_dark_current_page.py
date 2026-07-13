"""Dark current page AppTest coverage: temperature sweep via the sweep
orchestration utility from petringa/api/sweep.py.

ARCHITECTURE UNDER TEST: this page sweeps TEMPERATURE (not bias) at a fixed
operating bias, via `petringa.<sweep-utility>(param="T",
sim_fn=petringa.run_dark_current, ...)`, per 41-RESEARCH.md's Decision
Addendum and 41-UI-SPEC.md. Each test monkeypatches `petringa.run_dark_current`
as a MODULE ATTRIBUTE (the seam proven in tests/test_app_run_mockability.py /
39-01) BEFORE calling `at.run()`, so the real (expensive) devsim solve never
executes -- but the sweep utility class itself is NEVER monkeypatched: its
real `.run()` orchestration logic executes against the faked
`run_dark_current`, proving the page genuinely wires the sweep utility
end-to-end rather than a hand-rolled loop.

AppTest 1.55 has no `plotly_chart` / `download_button` accessor (39/40-
RESEARCH precedent, see tests/test_app_field_page.py), so assertions are
limited to at.exception, at.session_state, at.button, at.warning, at.info,
and at.error.
"""

from __future__ import annotations

import numpy as np
from streamlit.testing.v1 import AppTest

import petringa
from petringa import DeviceConfig, SimResult


def _fake_run_dark_current(cfg, v_start, v_stop, n_points=1, **kwargs):
    """Clean single-point fake: returns one operating-point value at v_stop."""
    return SimResult(
        config=cfg,
        sim_type="dark_current",
        x=np.array([v_stop]),
        y=np.array([1e-13]),
        metadata={
            "I_SRH": np.array([1e-13]),
            "I_TAT": np.array([-1e-14]),
            "I_SRV": np.array([0.0]),
            "area_cm2": 1e-4,
        },
    )


def _run_dark_current_page():
    from app.workflows.dark_current import render

    render()


def test_empty_state_guard(monkeypatch):
    monkeypatch.setattr(petringa, "run_dark_current", _fake_run_dark_current)

    at = AppTest.from_function(_run_dark_current_page)
    # session_state is empty by default (do not pre-seed device_config).
    at.run()

    assert at.exception == []
    info_texts = [el.value for el in at.info]
    assert "Configure a device in the sidebar to begin." in info_texts


def test_2d_config_shows_1d_only_warning(monkeypatch):
    monkeypatch.setattr(petringa, "run_dark_current", _fake_run_dark_current)

    at = AppTest.from_function(_run_dark_current_page)
    at.session_state["device_config"] = DeviceConfig(half_width_um=50.0)
    at.run()

    assert at.exception == []
    assert any("1D-only" in (w.value or "") for w in at.warning)


def test_run_uses_parametric_sweep_and_caches_result(monkeypatch):
    monkeypatch.setattr(petringa, "run_dark_current", _fake_run_dark_current)

    at = AppTest.from_function(_run_dark_current_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on Run click: {at.exception}"

    result = at.session_state["dark_current_result"]
    assert result.sim_type == "dark_current"
    n_requested = at.session_state["dark_current_n_requested"]
    assert len(result.x) == n_requested
    assert n_requested == 6  # default n_temperatures

    # Proves x holds TEMPERATURES (e.g. 250-400 K), not bias (0 to -100 V):
    # monotonically increasing and within the default T_min/T_max range.
    assert np.all(np.diff(result.x) > 0)
    assert result.x.min() >= 250.0
    assert result.x.max() <= 400.0


def test_partial_temperature_failure_shows_warning_not_crash(monkeypatch):
    calls = {"n": 0}

    def _fake_with_one_failure(cfg, v_start, v_stop, n_points=1, **kwargs):
        calls["n"] += 1
        if calls["n"] == 3:
            return SimResult(
                config=cfg,
                sim_type="dark_current",
                x=np.array([]),
                y=np.array([]),
                metadata={
                    "I_SRH": np.array([]),
                    "I_TAT": np.array([]),
                    "I_SRV": np.array([]),
                },
            )
        return _fake_run_dark_current(cfg, v_start, v_stop, n_points, **kwargs)

    monkeypatch.setattr(petringa, "run_dark_current", _fake_with_one_failure)

    at = AppTest.from_function(_run_dark_current_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on partial failure: {at.exception}"
    assert len(at.warning) >= 1
    warning_texts = [w.value or "" for w in at.warning]
    assert any(
        "temperature" in text and "5" in text and "6" in text for text in warning_texts
    ), f"expected a partial-failure warning mentioning counts, got: {warning_texts}"

    result = at.session_state["dark_current_result"]
    n_requested = at.session_state["dark_current_n_requested"]
    assert len(result.x) == n_requested - 1


def test_first_call_runtime_error_shows_error_not_crash(monkeypatch):
    def _raise_run_dark_current(cfg, v_start, v_stop, n_points=1, **kwargs):
        raise RuntimeError("dark_current_sweep: failed to converge")

    monkeypatch.setattr(petringa, "run_dark_current", _raise_run_dark_current)

    at = AppTest.from_function(_run_dark_current_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert (
        at.exception == []
    ), f"page crashed instead of showing st.error: {at.exception}"
    assert len(at.error) >= 1
    assert "dark_current_result" not in at.session_state
