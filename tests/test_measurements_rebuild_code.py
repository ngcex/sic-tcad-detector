"""The donor-table rebuild snippet reproduces the device the C-V was fitted on (issue #20).

Before the fix the Measurements page showed ``cfg = etna.DeviceConfig(donor_table=table)``
("same geometry as the sidebar device"): with an 8 um sidebar epi it raised (the table
stops at 8 um < the default 10 um epi), with the legacy preset it silently changed the
substrate, n+, mobility, N_t and E_t.
"""

from __future__ import annotations

import json

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

import etna
import etna.measurements as M
from etna import DeviceConfig
from etna.measurements.cv import default_knots_um, table_config
from tests._measurement_synth import cv_csv


def _fitted(base):
    k = default_knots_um(base.epi_thickness_um)
    logn = np.linspace(15.5, 13.8, k.size)
    return table_config(base, k, logn)


def _exec_snippet(code, cfg, tmp_path, monkeypatch, name="donor_table_cv.json"):
    rec = {
        "schema": "etna.donor_table/1",
        "donor_table": [list(p) for p in cfg.donor_table],
    }
    (tmp_path / name).write_text(json.dumps(rec))
    monkeypatch.chdir(tmp_path)
    ns: dict = {}
    exec(code, ns)
    return ns["cfg"]


@pytest.mark.parametrize(
    "base",
    [
        DeviceConfig(epi_thickness_um=8.0),
        etna.legacy_config(),
        etna.nplus_calibrated_config(),
        DeviceConfig(),
    ],
    ids=["epi8", "legacy", "nplus_calibrated", "default"],
)
def test_rebuild_code_round_trips(base, tmp_path, monkeypatch):
    fitted = _fitted(base)
    code = M.device_rebuild_code(fitted)
    assert _exec_snippet(code, fitted, tmp_path, monkeypatch) == fitted


def test_old_snippet_was_wrong(tmp_path, monkeypatch):
    """Documents the defect: the defaults-based snippet fails or changes the device."""
    old = (
        "import etna, etna.measurements as M\n"
        "table = M.load_donor_table_json('donor_table_cv.json')\n"
        "cfg = etna.DeviceConfig(donor_table=table)"
    )
    with pytest.raises(ValueError, match="cover the epi"):
        _exec_snippet(
            old, _fitted(DeviceConfig(epi_thickness_um=8.0)), tmp_path, monkeypatch
        )
    leg = _fitted(etna.legacy_config())
    assert _exec_snippet(old, leg, tmp_path, monkeypatch) != leg


def _page():
    from app.workflows.measurements import render

    render()


def test_page_snippet_rebuilds_the_sidebar_device(tmp_path, monkeypatch):
    real = M.invert_cv
    monkeypatch.setattr(M, "invert_cv", lambda d, c, **kw: real(d, c, max_nfev=2, **kw))
    at = AppTest.from_function(_page, default_timeout=180)
    at.session_state["device_config"] = DeviceConfig(epi_thickness_um=8.0)
    at.run()
    r = etna.run_cv(etna.measured_cv_tabulated_config(), 0.0, -30.0, 31)
    up = next(u for u in at.get("file_uploader") if u.key == "meas_cv_file")
    up.set_value(("cv.csv", cv_csv(r.x, r.y, area_cm2=0.04).encode(), "text/csv"))
    at.run()
    next(b for b in at.button if b.key == "meas_cv_invert").click()
    at.run()
    assert at.exception == []
    inv = at.session_state["meas_cv_inv"]
    code = next(c.value for c in at.code if "load_donor_table_json" in c.value)
    (tmp_path / "donor_table_cv.json").write_text(inv.to_json())
    monkeypatch.chdir(tmp_path)
    ns: dict = {}
    exec(code, ns)
    assert ns["cfg"] == inv.config
    assert ns["cfg"].epi_thickness_um == 8.0
