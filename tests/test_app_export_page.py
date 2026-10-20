"""Export page AppTest coverage (etna.run_all / etna.run_field mocked: no solver).

AppTest has no download_button accessor; the file bytes are covered by the pure
builders in tests/test_export.py. Here: boot, the Compute buttons, report reuse,
session caching, and the main-script navigation.
"""

from __future__ import annotations

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

import etna
from etna import DeviceConfig, MeshData, SimResult
from etna.api.run_all import Report, SectionResult

LABEL_START = "simulation output, validated against devsim reference"

# the default layout parameter file has PLACEHOLDER values; the page shows them as
# st.warning (asserted below), the Python warning is silenced here
pytestmark = pytest.mark.filterwarnings(
    "ignore::etna.export.layout_params.PlaceholderDimensionWarning"
)


def _page():
    from app.workflows.export import render

    render()


def _fake_report(cfg, sections=("cv", "cce", "field", "dark")):
    bias = np.linspace(0.0, -50.0, 11)
    W = np.linspace(1.7, 9.9, 11) * 1e-4
    md = dict(
        depletion_widths=W,
        solver="batched",
        validated_domain=True,
        domain_note=None,
        reach_through_V=10.3,
        reach_through_known=True,
    )
    res = {
        "cv": SimResult(cfg, "cv", bias, 1e-12 / W, md),
        "cce": SimResult(cfg, "cce", bias[2:], np.full(9, 0.999), dict(md)),
        "field": SimResult(
            cfg,
            "field",
            np.array([0.0, 0.3, 11.3]),
            np.array([0, 1e5, 0.0]),
            dict(bias_V=-50.0, solver="devsim"),
        ),
        "dark": SimResult(
            cfg,
            "dark_current",
            np.linspace(0, -100, 11),
            np.full(11, -4e-14),
            dict(solver="devsim"),
        ),
    }
    return Report(
        cfg,
        [SectionResult(n, n, "ok", solver="x", result=res[n]) for n in sections],
        0.1,
    )


def _fake_field(cfg, bias_V=-50.0):
    x = np.array([0.0, 0.3, 11.3]) * 1e-4
    nv = {
        "NetDoping": np.array([-1e18, 3e15, 1e19]),
        "Potential": np.zeros(3),
        "ElectricField": np.array([0.0, 1e5, 0.0]),
    }
    return SimResult(
        cfg,
        "field",
        x * 1e4,
        nv["ElectricField"],
        dict(bias_V=bias_V, solver="devsim", validated_domain=True),
        MeshData(x, None, nv, [], []),
    )


def _start(monkeypatch, calls):
    def fake_run_all(cfg, sections=None):
        calls.append(("run_all", tuple(sections)))
        return _fake_report(cfg)

    def fake_run_field(cfg, bias_V=-50.0):
        calls.append(("run_field", bias_V))
        return _fake_field(cfg, bias_V)

    monkeypatch.setattr(etna, "run_all", fake_run_all)
    monkeypatch.setattr(etna, "run_field", fake_run_field)
    at = AppTest.from_function(_page, default_timeout=60)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    assert at.exception == [], at.exception
    return at


def _button(at, key):
    return next(b for b in at.button if b.key == key)


def test_boot_shows_label_and_cross_section(monkeypatch):
    calls = []
    at = _start(monkeypatch, calls)
    assert any(i.value.startswith(LABEL_START) for i in at.info)
    assert calls == []  # nothing solved on boot
    assert "Depth to scale" in " ".join(c.value for c in at.caption)


def test_spec_sheet_button_runs_spec_sections(monkeypatch):
    calls = []
    at = _start(monkeypatch, calls)
    _button(at, "ex_spec_btn").click()
    at.run()
    assert at.exception == [], at.exception
    assert calls == [("run_all", ("cv", "cce", "field", "dark"))]
    assert at.session_state["export_report"].config == DeviceConfig()
    assert len(at.error) == 0


def test_spec_sheet_reuses_matching_run_all_report(monkeypatch):
    calls = []
    at = _start(monkeypatch, calls)
    at.session_state["run_all_report"] = _fake_report(DeviceConfig())
    _button(at, "ex_spec_btn").click()
    at.run()
    assert at.exception == [] and calls == []


def test_profile_button_runs_field_at_bias(monkeypatch):
    calls = []
    at = _start(monkeypatch, calls)
    at.number_input(key="ex_bias").set_value(-20.0)
    _button(at, "ex_field_btn").click()
    at.run()
    assert at.exception == [], at.exception
    assert calls == [("run_field", -20.0)]
    assert at.session_state["export_field"].metadata["bias_V"] == -20.0
    assert "3 nodes at -20 V" in " ".join(c.value for c in at.caption)


def test_page_without_config():
    at = AppTest.from_function(_page)
    at.run()
    assert at.exception == []


def test_export_page_in_main_navigation():
    at = AppTest.from_file("app/main.py", default_timeout=60)
    at.run()
    assert at.exception == []
    # callable st.Page entries cannot be switched to by path in AppTest; check the wiring
    src = open("app/main.py").read()
    assert 'url_path="export"' in src and "render_export" in src


def test_animation_unavailable_shows_reason(monkeypatch):
    import etna.viz as viz

    monkeypatch.setattr(
        viz, "manim_available", lambda: (False, "manim is not importable (test)")
    )
    at = _start(monkeypatch, [])
    assert any("manim is not importable (test)" in i.value for i in at.info)
    assert not [b for b in at.button if b.key == "ex_anim_btn"]


def test_animation_button_renders_videos(monkeypatch, tmp_path):
    import etna.viz as viz

    rendered = []

    def fake_render(data, scene, out_dir=None, **k):
        p = tmp_path / f"{scene}.mp4"
        p.write_bytes(b"\x00\x00\x00\x18ftypmp42fake")
        rendered.append(scene)
        return p, 0.1

    monkeypatch.setattr(viz, "manim_available", lambda: (True, ""))
    monkeypatch.setattr(viz, "animation_data", lambda cfg: {"frames": []})
    monkeypatch.setattr(viz, "render_animation", fake_render)
    at = _start(monkeypatch, [])
    _button(at, "ex_anim_btn").click()
    at.run()
    assert at.exception == [], at.exception
    assert rendered == ["StackScene", "DepletionScene"]
    cfg, vids = at.session_state["export_animation"]
    assert cfg == DeviceConfig() and set(vids) == {"StackScene", "DepletionScene"}


def test_layout_block_shows_disclaimer_and_degrades_without_gdstk(monkeypatch):
    from etna import export as ex

    at = _start(monkeypatch, [])
    assert ex.LAYOUT_DISCLAIMER in [w.value for w in at.warning]
    assert "REFERENCE-NOT-FAB" in " ".join(h.value for h in at.subheader)

    monkeypatch.setattr(ex, "gds_available", lambda: (False, "gdstk not installed"))
    at = _start(monkeypatch, [])
    assert any("GDSII download unavailable" in i.value for i in at.info)
    assert ex.LAYOUT_DISCLAIMER in [w.value for w in at.warning]
