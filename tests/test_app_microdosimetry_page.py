"""Microdosimetry page AppTest coverage: empty-state guard + no-file warning +
upload->Run->cache happy path + malformed-CSV error-not-crash.

Each test monkeypatches etna.run_microdosimetry as a MODULE ATTRIBUTE
(the seam proven in tests/test_app_run_mockability.py / 39-01) BEFORE
calling at.run(), so the real data pipeline never executes and the happy
path stays deterministic. Mirrors tests/test_app_radiation_damage_page.py's
structure: a `_run_microdosimetry_page()` wrapper importing and calling
app.workflows.microdosimetry.render, and assertions limited to at.exception,
at.session_state, at.button, at.warning, at.info, and at.error (no
plotly_chart / download_button accessor).

The file_uploader is driven via at.file_uploader[0].upload(filename,
content, mime_type) -- a streamlit 1.58 capability with no in-repo
precedent. `.upload` takes THREE POSITIONAL args (filename, content bytes,
mime_type), NOT a single (name, bytes, mime) tuple -- confirmed against the
installed streamlit.testing.v1.element_tree.FileUploader.upload signature.
No 1D-guard test and no NaN-tolerance test are present: those mechanisms are
intentionally absent from this page.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import etna
from etna import DeviceConfig, SimResult
from streamlit.testing.v1 import AppTest


def _fake_run_microdosimetry(cfg, mc_csv_path, sv_thickness_um=10.0, sv_width_um=150.0):
    """Clean fake: 300-point log-x spectrum with a valid y_F/y_D readout."""
    return SimResult(
        config=cfg,
        sim_type="microdosimetry",
        x=np.geomspace(0.01, 9000, 300),
        y=np.random.rand(300),
        metadata={"y_F": 17.23, "y_D": 53.22, "l_bar_um": 20.0},
    )


def _run_microdosimetry_page():
    from app.workflows.microdosimetry import render

    render()


def test_empty_state_guard(monkeypatch):
    monkeypatch.setattr(etna, "run_microdosimetry", _fake_run_microdosimetry)

    at = AppTest.from_function(_run_microdosimetry_page)
    # session_state is empty by default (do not pre-seed device_config).
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"
    info_texts = [el.value for el in at.info]
    assert "Configure a device in the sidebar to begin." in info_texts


def test_no_file_on_run_warns(monkeypatch):
    monkeypatch.setattr(etna, "run_microdosimetry", _fake_run_microdosimetry)

    at = AppTest.from_function(_run_microdosimetry_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on Run click: {at.exception}"
    warning_texts = [el.value for el in at.warning]
    assert any(
        "Upload an MC events CSV" in (w or "") for w in warning_texts
    ), f"no-file warning not found: {warning_texts}"
    assert "microdosimetry_result" not in at.session_state


def test_upload_run_caches_spectrum(monkeypatch):
    monkeypatch.setattr(etna, "run_microdosimetry", _fake_run_microdosimetry)

    at = AppTest.from_function(_run_microdosimetry_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    data = Path("data/synthetic_mc_events.csv").read_bytes()
    at.file_uploader[0].upload("events.csv", data, "text/csv")
    at.run()

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on upload+Run: {at.exception}"
    result = at.session_state["microdosimetry_result"]
    assert result.sim_type == "microdosimetry"


def test_malformed_csv_shows_error_not_crash(monkeypatch):
    def _raise_run_microdosimetry(cfg, mc_csv_path, **kwargs):
        raise ValueError("bad columns")

    monkeypatch.setattr(etna, "run_microdosimetry", _raise_run_microdosimetry)

    at = AppTest.from_function(_run_microdosimetry_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.file_uploader[0].upload("junk.csv", b"not,a,valid\n1,2", "text/csv")
    at.run()

    at.button[0].click()
    at.run()

    assert (
        at.exception == []
    ), f"page crashed instead of showing st.error: {at.exception}"
    assert len(at.error) >= 1
    assert "microdosimetry_result" not in at.session_state
