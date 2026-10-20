"""Measurements page AppTest: upload (SYNTHETIC files generated in-test from the model),
validation messages, model overlay, results cached for the downloads.

The C-V inversion is wrapped to a few iterations (max_nfev=4) through the
``etna.measurements`` module attribute the page calls, so the page path runs end to end
without the full fit; the other analyses run for real (batched engine, < 1 s each).
"""

from __future__ import annotations

import numpy as np
from streamlit.testing.v1 import AppTest

import etna
import etna.measurements as M
from etna import DeviceConfig
from tests._measurement_synth import cv_csv, it_csv, sims_csv

TIMEOUT = 180


def _page():
    from app.workflows.measurements import render

    render()


def _app(cfg=None):
    at = AppTest.from_function(_page, default_timeout=TIMEOUT)
    if cfg is not None:
        at.session_state["device_config"] = cfg
    return at


def _uploader(at, key):
    return next(u for u in at.get("file_uploader") if u.key == key)


def _button(at, key):
    return next(b for b in at.button if b.key == key)


def _cv_bytes(area=True):
    r = etna.run_cv(etna.measured_cv_tabulated_config(), 0.0, -30.0, 61)
    return cv_csv(r.x, r.y, area_cm2=0.04 if area else None).encode()


def test_empty_state_and_2d_guard():
    at = _app()
    at.run()
    assert at.exception == []
    assert "Configure a device in the sidebar to begin." in [i.value for i in at.info]
    at = _app(DeviceConfig(half_width_um=50.0))
    at.run()
    assert at.exception == []
    assert any("1D-only" in (w.value or "") for w in at.warning)


def test_tabs_and_upload_prompts():
    at = _app(DeviceConfig())
    at.run()
    assert at.exception == []
    assert [t.label for t in at.tabs] == ["C-V", "SIMS / SRP", "I(T)"]
    infos = [i.value for i in at.info]
    assert "Upload a C-V file to begin." in infos
    assert "Upload an I(T) file to begin." in infos


def test_bad_cv_file_shows_line_and_column():
    at = _app(DeviceConfig())
    at.run()
    _uploader(at, "meas_cv_file").set_value(
        ("bad.csv", b"# area_cm2: 0.04\nV_V,C_pF\n0,200\n-1,abc\n", "text/csv")
    )
    at.run()
    assert at.exception == []
    errs = [e.value for e in at.error]
    assert any("1 problem(s)" in e for e in errs)
    assert any("line 4, column 'C_pF'" in e for e in errs)


def test_cv_compare_and_invert(monkeypatch):
    real = M.invert_cv
    monkeypatch.setattr(M, "invert_cv", lambda d, c, **kw: real(d, c, max_nfev=4, **kw))
    at = _app(DeviceConfig())
    at.run()
    _uploader(at, "meas_cv_file").set_value(
        ("cv.csv", _cv_bytes(area=False), "text/csv")
    )
    at.run()
    assert at.exception == []
    assert any("61 points loaded" in s.value for s in at.success)
    assert any("sidebar device area" in i.value for i in at.info)  # file without area
    _button(at, "meas_cv_compare").click()
    at.run()
    assert at.exception == []
    cmp = at.session_state["meas_cv_cmp"]
    assert cmp.n_points == 61 and 0.0 < cmp.rms_W_um < 0.5
    _button(at, "meas_cv_invert").click()
    at.run()
    assert at.exception == []
    inv = at.session_state["meas_cv_inv"]
    assert inv.config.donor_table == inv.donor_table and len(inv.donor_table) == 19


def test_sims_resample():
    cfg = DeviceConfig()
    depth = np.arange(0.0, 11.0, 0.02)
    N = M.model_donor_profile(cfg, depth - 0.3)
    N = np.where(depth < 0.3, N.max(), N)
    at = _app(cfg)
    at.run()
    _uploader(at, "meas_sims_file").set_value(
        (
            "sims.csv",
            sims_csv(depth, N, meta=["# junction_depth_um: 0.3"]).encode(),
            "text/csv",
        )
    )
    at.run()
    assert at.exception == []
    _button(at, "meas_sims_run").click()
    at.run()
    assert at.exception == []
    tab = at.session_state["meas_sims_tab"]
    assert tab.junction_source == "file metadata" and len(tab.knots_um) == 41


def test_it_fit_and_model_comparison():
    cfg = DeviceConfig()
    T = np.array([290.0, 300.0, 310.0, 320.0])
    Im, _ = M.model_current_vs_T(cfg, T, -30.0)
    at = _app(cfg)
    at.run()
    _uploader(at, "meas_it_file").set_value(
        ("it.csv", it_csv(T, Im, meta=["# bias_V: -30"]).encode(), "text/csv")
    )
    at.run()
    assert at.exception == []
    _button(at, "meas_it_run").click()
    at.run()
    assert at.exception == []
    res = at.session_state["meas_it_res"]
    assert abs(res.comparisons[0].delta_eV) < 1e-6
