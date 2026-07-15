"""Batch sweep page AppTest coverage: a general-case parametric sweep via the
real ``etna.ParametricSweep(...).run()`` orchestration.

Each test monkeypatches the FACADE (``etna.run_cce``) as a MODULE ATTRIBUTE
(the seam proven in ``tests/test_app_run_mockability.py`` / 39-01) BEFORE calling
``at.run()``, so the expensive real devsim solve never executes -- but
``ParametricSweep`` itself is NEVER monkeypatched: its real ``.run()`` logic runs
against the faked facade, proving the page genuinely wires ParametricSweep
end-to-end rather than a hand-rolled loop.

AppTest has no ``plotly_chart`` / ``download_button`` accessor (39/40/41-RESEARCH
precedent), so assertions are limited to ``at.exception``, ``at.session_state``,
``at.button``, ``at.warning``, ``at.info``, and ``at.error``.
"""

from __future__ import annotations

import numpy as np
from streamlit.testing.v1 import AppTest

import etna
from etna import DeviceConfig, SimResult


def _fake_run_cce(cfg, **kwargs):
    """Clean overlayable CCE-vs-bias curve fake (no devsim solve)."""
    return SimResult(
        config=cfg,
        sim_type="cce",
        x=np.linspace(0, -40, 5),
        y=np.linspace(1.0, 0.6, 5),
        metadata={"I_collected": np.zeros(5), "I_generated": 1.0},
    )


def _run_batch_sweep_page():
    from app.workflows.batch_sweep import render

    render()


def test_empty_state_guard(monkeypatch):
    monkeypatch.setattr(etna, "run_cce", _fake_run_cce)

    at = AppTest.from_function(_run_batch_sweep_page)
    # session_state is empty by default (do not pre-seed device_config).
    at.run()

    assert at.exception == []
    info_texts = [el.value for el in at.info]
    assert "Configure a device in the sidebar to begin." in info_texts


def test_2d_config_shows_1d_only_warning(monkeypatch):
    monkeypatch.setattr(etna, "run_cce", _fake_run_cce)

    at = AppTest.from_function(_run_batch_sweep_page)
    at.session_state["device_config"] = DeviceConfig(half_width_um=50.0)
    at.run()

    assert at.exception == []
    assert any("1D-only" in (w.value or "") for w in at.warning)


def test_run_uses_parametric_sweep_and_caches_results(monkeypatch):
    # Monkeypatch the FACADE only; the real ParametricSweep.run() must execute.
    monkeypatch.setattr(etna, "run_cce", _fake_run_cce)

    at = AppTest.from_function(_run_batch_sweep_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on Run click: {at.exception}"

    # Default value list "10, 15, 20" -> 3 swept values -> >= 3 cached results.
    assert len(at.session_state["sweep_results"]) >= 3
    # Read the RENAMED run-snapshot key (NOT the live sweep_param widget value).
    assert at.session_state["sweep_run_param"] == "epi_thickness_um"


def test_bad_value_list_shows_error(monkeypatch):
    monkeypatch.setattr(etna, "run_cce", _fake_run_cce)

    at = AppTest.from_function(_run_batch_sweep_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.text_input(key="sweep_values").set_value("a, b, c")
    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on bad value list: {at.exception}"
    assert len(at.error) >= 1
    assert any("comma-separated list of numbers" in (e.value or "") for e in at.error)
    assert "sweep_results" not in at.session_state


def test_partial_value_failure_shows_warning_not_crash(monkeypatch):
    calls = {"n": 0}

    def _fake_with_one_failure(cfg, **kwargs):
        calls["n"] += 1
        if calls["n"] == 2:
            return SimResult(
                config=cfg,
                sim_type="cce",
                x=np.array([]),
                y=np.array([]),
                metadata={"I_collected": np.array([]), "I_generated": 1.0},
            )
        return _fake_run_cce(cfg, **kwargs)

    monkeypatch.setattr(etna, "run_cce", _fake_with_one_failure)

    at = AppTest.from_function(_run_batch_sweep_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on partial failure: {at.exception}"
    assert len(at.warning) >= 1
    warning_texts = [w.value or "" for w in at.warning]
    assert any(
        "2" in text and "3" in text and "completed successfully" in text
        for text in warning_texts
    ), f"expected a partial-failure warning mentioning counts, got: {warning_texts}"

    n_requested = at.session_state["sweep_n_requested"]
    assert len(at.session_state["sweep_results"]) == n_requested - 1


def test_facade_runtime_error_shows_error_not_crash(monkeypatch):
    def _raise_run_cce(cfg, **kwargs):
        raise RuntimeError("failed to converge")

    monkeypatch.setattr(etna, "run_cce", _raise_run_cce)

    at = AppTest.from_function(_run_batch_sweep_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed instead of st.error: {at.exception}"
    assert len(at.error) >= 1
    assert "sweep_results" not in at.session_state
