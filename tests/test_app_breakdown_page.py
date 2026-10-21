"""Breakdown & Avalanche page (app/workflows/breakdown.py): AppTest coverage.

The three facades are the 1D batched engine (milliseconds without the mesh check), so
the default runs are real; the error / guard paths use the module-attribute seam
(monkeypatch etna.run_breakdown).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import etna
from etna import DeviceConfig
from streamlit.testing.v1 import AppTest

TIMEOUT = 60


def _page():
    from app.workflows.breakdown import render

    render()


def _start(cfg=None):
    at = AppTest.from_function(_page, default_timeout=TIMEOUT)
    if cfg is not None:
        at.session_state["device_config"] = cfg
    at.run()
    assert at.exception == [], at.exception
    return at


def test_empty_state_guard():
    at = _start()
    assert "Configure a device in the sidebar to begin." in [i.value for i in at.info]


def test_2d_guard():
    at = _start(DeviceConfig(half_width_um=20.0))
    assert any("1D-only" in w.value for w in at.warning)
    assert len(at.tabs) == 0


def test_breakdown_real_run_and_labels():
    at = _start(DeviceConfig())
    at.checkbox(key="bd_mesh").uncheck()
    at.run()
    at.button(key="bd_run").click()
    at.run()
    assert at.exception == [], at.exception
    r = at.session_state["bd_result"]
    assert r.sim_type == "breakdown"
    assert 1000 < r.metadata["V_BD"] < 5000
    assert any("PREDICTED, not validated" in w.value for w in at.warning)
    assert any("V_BD" in m.label for m in at.metric)
    # cached: a second click with the same inputs reuses the same object
    at.button(key="bd_run").click()
    at.run()
    assert at.session_state["bd_result"] is r


def test_breakdown_needs_impact_ionization():
    at = _start(DeviceConfig())
    at.selectbox(key="bd_ii").set_value("(none)")
    at.run()
    at.button(key="bd_run").click()
    at.run()
    assert any("impact-ionization set" in e.value for e in at.error)


def test_avalanche_dark_and_transit_real_runs():
    at = _start(DeviceConfig())
    at.radio(key="av_src").set_value("dark")
    at.selectbox(key="bd_hf").set_value("lades2000")
    at.run()
    at.button(key="av_run").click()
    at.run()
    assert at.exception == [], at.exception
    r = at.session_state["av_result"]
    assert r.sim_type == "avalanche_dark"
    assert np.all(np.asarray(r.metadata["gain"]) >= 0.99)
    assert any("Dark model" in i.value for i in at.info)
    at.button(key="tt_run").click()
    at.run()
    t = at.session_state["tt_result"]
    assert t.sim_type == "transit_time"
    # velocity saturation slows the carriers: t(models) >= t(low field)
    assert np.all(
        np.asarray(t.metadata["t_e"]) >= np.asarray(t.metadata["t_e_low_field"]) * 0.999
    )


def test_temperature_out_of_range_is_refused_then_allowed():
    at = _start(DeviceConfig(T=350.0))
    at.checkbox(key="bd_mesh").uncheck()
    at.run()
    at.button(key="bd_run").click()
    at.run()
    assert at.exception == []
    assert len(at.error) >= 1  # hatakeyama2004 is a 300 K set
    at.checkbox(key="bd_allowT").check()
    at.run()
    at.button(key="bd_run").click()
    at.run()
    assert at.session_state["bd_result"].metadata["T_extrapolated"] is True
    assert any("extrapolated" in w.value for w in at.warning)


def test_solver_failure_is_reported_not_raised(monkeypatch):
    def boom(*a, **k):
        raise RuntimeError("no convergence")

    monkeypatch.setattr(etna, "run_breakdown", boom)
    at = _start(DeviceConfig())
    at.button(key="bd_run").click()
    at.run()
    assert at.exception == []
    assert any("did not converge" in e.value for e in at.error)
    assert "bd_result" not in at.session_state


def test_publication_figure_preview():
    at = _start(DeviceConfig())
    at.checkbox(key="bd_mesh").uncheck()
    at.run()
    at.button(key="bd_run").click()
    at.run()
    at.toggle(key="bd_fig_prev").set_value(True)
    at.run()
    assert at.exception == []
    assert len(at.get("imgs")) == 1


@pytest.mark.parametrize(
    "url", ["breakdown", "mc-deposition", "track", "trap-kinetics", "structures-2d"]
)
def test_new_pages_in_main_navigation(url):
    src = Path("app/main.py").read_text()
    assert f'url_path="{url}"' in src
