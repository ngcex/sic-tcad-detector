"""Track engine page (app/workflows/track.py): AppTest coverage with fake facades
(a real track event is 2-15 s; tests/test_track_engine.py covers the physics)."""

from __future__ import annotations

import numpy as np

import etna
from etna import DeviceConfig, SimResult
from streamlit.testing.v1 import AppTest

TIMEOUT = 30
CALLS = []


def _fake_event(cfg, LET_keV_um=100.0, V_bias=50.0, **kw):
    CALLS.append(("event", LET_keV_um, V_bias, kw))
    t = np.geomspace(1e-15, 1e-8, 50)
    md = dict(cce=0.9999, recombined_fraction=1e-4, wrong_contact=1e-8, N_gen=1e5,
              ledger_residual=1e-9, ledger_floor=1e-13, converged=True, n_steps=50,
              wall_s=0.1, excess=True, peak_excess_cm3=1e19, LET_keV_um=LET_keV_um,
              V_bias=V_bias, mesh=kw.get("mesh"), solver="fake", label="predicted (fake)")  # fmt: skip
    return SimResult(cfg, "track_event", t, 1e9 * np.exp(-t / 1e-10), md)


def _fake_let(cfg, LET_values=(1.0,), V_bias=50.0, **kw):
    CALLS.append(("let", tuple(LET_values), V_bias, kw))
    x = np.asarray(LET_values, float)
    y = 1 - 1e-5 * x
    md = dict(PHD=1 - y, recombined_fraction=1e-5 * x, ledger_residual=np.zeros_like(x),
              converged=np.ones_like(x, bool), wall_s=np.zeros_like(x), V_bias=V_bias,
              mesh=kw.get("mesh"), solver="fake", label="predicted (fake)")  # fmt: skip
    return SimResult(cfg, "track_cce_let", x, y, md)


def _page():
    from app.workflows.track import render

    render()


def _start(monkeypatch, cfg=DeviceConfig()):
    CALLS.clear()
    monkeypatch.setattr(etna, "run_track_event", _fake_event)
    monkeypatch.setattr(etna, "run_track_cce_let", _fake_let)
    at = AppTest.from_function(_page, default_timeout=TIMEOUT)
    if cfg is not None:
        at.session_state["device_config"] = cfg
    at.run()
    assert at.exception == [], at.exception
    return at


def test_guards(monkeypatch):
    at = _start(monkeypatch, cfg=None)
    assert "Configure a device in the sidebar to begin." in [i.value for i in at.info]
    at = _start(monkeypatch, cfg=DeviceConfig(half_width_um=20.0))
    assert any("1D" in w.value for w in at.warning)


def test_single_event_inputs_reach_the_facade_and_results_render(monkeypatch):
    at = _start(monkeypatch)
    at.number_input(key="tr_let").set_value(1000.0)
    at.selectbox(key="tr_radial").set_value("gaussian")
    at.selectbox(key="tr_mesh").set_value("tiny")
    at.run()
    at.button(key="tr_run").click()
    at.run()
    assert at.exception == [], at.exception
    kind, let, V, kw = CALLS[-1]
    assert (kind, let, V) == ("event", 1000.0, 50.0)
    assert kw["radial"] == "gaussian" and kw["mesh"] == "tiny"
    assert any(m.value == "0.999900" for m in at.metric)
    assert any("PREDICTED (FAKE)" in c.value for c in at.caption)
    # same inputs: served from the session cache, no new facade call
    n = len(CALLS)
    at.button(key="tr_run").click()
    at.run()
    assert len(CALLS) == n
    # changed input -> stale warning until rerun
    at.number_input(key="tr_V").set_value(10.0)
    at.run()
    assert any(
        "different device or different page inputs" in w.value for w in at.warning
    )


def test_cce_vs_let_parses_and_validates(monkeypatch):
    at = _start(monkeypatch)
    at.text_input(key="tl_lets").set_value("1000, 1, 10")
    at.run()
    at.button(key="tl_run").click()
    at.run()
    assert CALLS[-1][1] == (1.0, 10.0, 1000.0)
    assert at.session_state["tl_result"].sim_type == "track_cce_let"
    at.text_input(key="tl_lets").set_value("1, abc")
    at.run()
    assert at.exception == []
    at.text_input(key="tl_lets").set_value("1, -2")
    at.run()
    assert any("positive" in e.value for e in at.error)


def test_tilted_alpha_is_a_script_with_the_calibration_caveat(monkeypatch):
    at = _start(monkeypatch)
    assert any("run_tilted_alpha_track" in c.value for c in at.code)
    assert any("does NOT reproduce" in w.value for w in at.warning)


def test_failure_is_reported(monkeypatch):
    at = _start(monkeypatch)

    def boom(*a, **k):
        raise RuntimeError("newton failed")

    monkeypatch.setattr(etna, "run_track_event", boom)
    at.button(key="tr_run").click()
    at.run()
    assert at.exception == []
    assert any("did not converge" in e.value for e in at.error)
