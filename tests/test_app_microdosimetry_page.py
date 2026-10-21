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


def _fake_run_microdosimetry(
    cfg, mc_csv_path, sv_thickness_um=10.0, sv_width_um=100.0, **kwargs
):
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

    data = (
        Path(__file__).resolve().parents[1] / "data" / "synthetic_mc_events.csv"
    ).read_bytes()
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


_FIXTURE = Path(__file__).parent.parent / "data" / "synthetic_mc_events.csv"


def _uploaded_page():
    at = AppTest.from_function(_run_microdosimetry_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    at.file_uploader[0].upload("events.csv", _FIXTURE.read_bytes(), "text/csv")
    at.run()
    return at


def test_real_pipeline_reads_the_upload_in_memory(monkeypatch):
    """audit C-17: the REAL run_microdosimetry gets an in-memory stream (no temp file)."""
    seen = []
    real = etna.run_microdosimetry

    def spy(cfg, mc_csv_path, **kw):
        seen.append(type(mc_csv_path).__name__)
        return real(cfg, mc_csv_path=mc_csv_path, **kw)

    monkeypatch.setattr(etna, "run_microdosimetry", spy)
    at = _uploaded_page()
    at.button[0].click()
    at.run()
    assert at.exception == []
    assert seen == ["BytesIO"]
    # the page passes its own width widget value (100 um, DESIGN-2's SV, = the facade
    # default); since 2026-10-20 the width enters l_bar (finite-SV Cauchy 4V/S), so
    # compare at the same geometry
    expected = real(DeviceConfig(), mc_csv_path=str(_FIXTURE), sv_width_um=100.0)
    assert at.session_state["microdosimetry_result"].metadata["y_F"] == (
        expected.metadata["y_F"]
    )


def test_changed_sensitive_volume_flags_the_cached_spectrum(monkeypatch):
    """audit C-11: the spectrum depends on page inputs, not the sidebar."""
    monkeypatch.setattr(etna, "run_microdosimetry", _fake_run_microdosimetry)
    at = _uploaded_page()
    at.button[0].click()
    at.run()
    assert not any("different page inputs" in (w.value or "") for w in at.warning)
    at.number_input(key="micro_sv_t").set_value(20.0)
    at.run()
    stale = [w.value for w in at.warning if "different page inputs" in (w.value or "")]
    assert stale and "sensitive-volume thickness" in stale[0]


def test_failed_rerun_clears_the_previous_spectrum(monkeypatch):
    """audit C-11 (b): a failed run does not leave the old spectrum on screen."""
    monkeypatch.setattr(etna, "run_microdosimetry", _fake_run_microdosimetry)
    at = _uploaded_page()
    at.button[0].click()
    at.run()
    assert "microdosimetry_result" in at.session_state

    def _bad(cfg, mc_csv_path, **kw):
        raise ValueError("the MC events CSV contains no events")

    monkeypatch.setattr(etna, "run_microdosimetry", _bad)
    at.button[0].click()
    at.run()
    assert at.exception == []
    assert any("no events" in (e.value or "") for e in at.error)
    assert "microdosimetry_result" not in at.session_state
