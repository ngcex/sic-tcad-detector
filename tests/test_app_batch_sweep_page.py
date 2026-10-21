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


def test_non_finite_or_non_positive_values_are_rejected(monkeypatch):
    """audit C-10: nan / inf / negative values never reach the facade."""
    calls = []
    monkeypatch.setattr(etna, "run_cce", lambda cfg, **k: calls.append(cfg))
    for raw in ("10, nan", "10, inf", "10, -5", "0, 10"):
        at = AppTest.from_function(_run_batch_sweep_page)
        at.session_state["device_config"] = DeviceConfig()
        at.run()
        at.text_input(key="sweep_values").set_value(raw)
        at.button[0].click()
        at.run()
        assert at.exception == [], raw
        assert any("finite, positive" in (e.value or "") for e in at.error), raw
        assert "sweep_results" not in at.session_state
    assert calls == []


def test_valid_values_still_run(monkeypatch):
    monkeypatch.setattr(etna, "run_cce", _fake_run_cce)
    at = AppTest.from_function(_run_batch_sweep_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    at.text_input(key="sweep_values").set_value("10, 15")
    at.button[0].click()
    at.run()
    assert at.exception == []
    assert len(at.session_state["sweep_results"]) == 2


def test_area_sweep_with_cce_facade_is_disabled_with_warning(monkeypatch):
    """audit C-02: area cannot change a dimensionless CCE -> no identical overlay."""
    monkeypatch.setattr(etna, "run_cce", _fake_run_cce)
    from app.workflows.batch_sweep import SIM_FACADES

    for label, attr in SIM_FACADES.items():
        at = AppTest.from_function(_run_batch_sweep_page)
        at.session_state["device_config"] = DeviceConfig()
        at.run()
        at.selectbox(key="sweep_param").set_value("area_cm2")
        at.selectbox(key="sweep_sim").set_value(label)
        at.run()
        assert at.exception == []
        warned = any(
            "does not depend on the device area" in (w.value or "") for w in at.warning
        )
        run_button = [b for b in at.button if b.label == "Run simulation"][0]
        if attr == "run_cv":  # C = eps*A/W: the area matters
            assert not warned and not run_button.disabled
        else:
            assert warned and run_button.disabled, label


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


def _tabulated_page_run(monkeypatch, raw):
    from etna.api.device import measured_cv_tabulated_config

    monkeypatch.setattr(etna, "run_cce", _fake_run_cce)
    at = AppTest.from_function(_run_batch_sweep_page)
    at.session_state["device_config"] = measured_cv_tabulated_config()
    at.run()
    at.selectbox(key="sweep_param").set_value("epi_thickness_um")
    at.run()
    at.text_input(key="sweep_values").set_value(raw)
    at.button[0].click()
    at.run()
    assert at.exception == [], at.exception
    return at


def test_invalid_value_reason_shown_inline_per_point(monkeypatch):
    """An epi thicker than the tabulated profile: the reason is on the page, not
    only inside the 'Error details' expander (which keeps the full message)."""
    at = _tabulated_page_run(monkeypatch, "6, 15")
    assert len(at.session_state["sweep_results"]) == 1
    texts = [w.value or "" for w in at.warning]
    hit = [t for t in texts if "1 of 2 values completed successfully" in t]
    assert hit, texts
    assert "= 15: the tabulated epi profile covers only 10 µm" in hit[0]
    assert "epi thickness ≤ 10 µm" in hit[0]
    assert "Invalid DeviceConfig" not in hit[0]  # raw text stays in the expander
    assert any(e.label == "Error details" for e in at.expander)
    assert any("must cover the epi" in (c.value or "") for c in at.code)


def test_all_values_invalid_reason_shown_inline(monkeypatch):
    at = _tabulated_page_run(monkeypatch, "15, 20")
    assert "sweep_results" not in at.session_state
    texts = [e.value or "" for e in at.error]
    hit = [t for t in texts if "All swept values failed" in t]
    assert hit, texts
    assert "= 15: the tabulated epi profile covers only 10 µm" in hit[0]
    assert "= 20: the tabulated epi profile covers only 10 µm" in hit[0]
    assert any("must cover the epi" in (c.value or "") for c in at.code)


def test_failure_reason_wording():
    from app.workflows.batch_sweep import NO_DATA_REASON, failure_reason

    assert failure_reason(None) == NO_DATA_REASON
    assert "did not converge" in failure_reason("RuntimeError: Newton diverged")
    assert failure_reason(
        "ValueError: Invalid DeviceConfig: T must be between 100 and 600 K (got 900.0)"
    ) == "invalid for this device: T must be between 100 and 600 K (got 900.0)"
    assert len(failure_reason("ValueError: " + "x" * 1000)) == 200
