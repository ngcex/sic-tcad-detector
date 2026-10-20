"""etna.viz: animation data (real-solver seam mocked), availability, one real render.

The render test skips cleanly when manim (optional extra ``viz``) is not installed.
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

import etna
from etna import DeviceConfig, SimResult
from etna.export import VALIDATION_LABEL
from etna.viz import SCENES, animation_data, manim_available, render_animation


def _fake_cv(cfg, v_start=0.0, v_stop=-60.0, n_points=40):
    x = np.linspace(v_start, v_stop, n_points)
    W = np.minimum(1.7 + 0.8 * np.abs(x), 9.9) * 1e-4
    return SimResult(
        cfg,
        "cv",
        x,
        1e-12 / W,
        dict(depletion_widths=W, solver="batched", reach_through_V=10.3),
    )


def _fake_field(cfg, bias_V=-50.0):
    x = np.array([0.0, 0.3, 0.5, 5.0, 10.25, 10.3, 10.35, 11.3])
    E = np.array([0.0, 9e5, 1e5 + 100 * abs(bias_V), 3e4, 2e4, 9e5, 1e4, 0.0])
    return SimResult(cfg, "field", x, E, dict(bias_V=bias_V, solver="devsim"))


def test_animation_data_uses_solvers(monkeypatch):
    calls = []
    monkeypatch.setattr(
        etna, "run_cv", lambda *a, **k: (calls.append("cv"), _fake_cv(*a, **k))[1]
    )
    monkeypatch.setattr(
        etna,
        "run_field",
        lambda *a, **k: (calls.append("field"), _fake_field(*a, **k))[1],
    )
    d = animation_data(DeviceConfig(), biases=(0.0, -5.0, -20.0))
    assert calls == ["cv", "field", "field", "field"]
    assert d["label"] == VALIDATION_LABEL
    assert [f["bias_V"] for f in d["frames"]] == [0.0, -5.0, -20.0]
    assert d["frames"][1]["W_um"] == pytest.approx(1.7 + 0.8 * 5)
    # the 9e5 spikes at both abrupt steps (0.3 and 10.3 um) are excluded from the peak
    assert d["frames"][2]["peak_E_V_per_cm"] == pytest.approx(1e5 + 2000)
    assert all(f["peak_excludes_interface"] for f in d["frames"])
    assert d["E_axis_max"] == pytest.approx(1.15 * (1e5 + 2000))
    assert [L["name"] for L in d["layers"]] == [
        "p+ anode",
        "n- epi",
        "n-/n+ transition",  # default stack since 2026-10-18
        "n+ substrate",
    ]


def test_animation_data_rejects_2d_and_forward_bias():
    with pytest.raises(ValueError, match="1D"):
        animation_data(DeviceConfig(half_width_um=25.0))
    with pytest.raises(ValueError, match="reverse"):
        animation_data(DeviceConfig(), biases=(0.0, 5.0))


def test_manim_unavailable_gives_reason(monkeypatch):
    import importlib.util

    real_find = importlib.util.find_spec
    monkeypatch.setattr(
        importlib.util,
        "find_spec",
        lambda name, *a: None if name == "manim" else real_find(name, *a),
    )
    ok, why = manim_available()
    assert not ok and "uv sync --extra dev --extra viz" in why and "pkgconf" in why
    with pytest.raises(RuntimeError, match="extra viz"):
        render_animation({}, "StackScene")
    # an import-blocked manim (sys.modules entry None) is unavailable too
    monkeypatch.setattr(importlib.util, "find_spec", real_find)
    monkeypatch.setitem(sys.modules, "manim", None)
    ok, why = manim_available()
    assert not ok and "extra viz" in why


def test_manim_available_does_not_import_manim():
    code = (
        "import sys; from etna.viz import manim_available; manim_available(); "
        "print('LOADED=' + str('manim' in sys.modules))"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    (line,) = [ln for ln in out.stdout.splitlines() if ln.startswith("LOADED=")]
    assert line == "LOADED=False"


def test_animate_auto_falls_back_when_manim_render_fails(monkeypatch, tmp_path):
    import etna.viz.animations as A

    monkeypatch.setattr(A, "manim_available", lambda: (True, ""))

    def broken(*a, **k):
        raise RuntimeError("manim render of X failed (exit 1): cairo")

    monkeypatch.setattr(A, "render_animation", broken)
    seen = {}

    def frames(data, scene, fmt, out_dir=None):
        seen["fmt"] = fmt
        return Path(tmp_path) / "x.mp4", 0.1

    monkeypatch.setattr("etna.viz.frames.render_frames", frames)
    _, _, engine = A.animate({}, "StackScene", "mp4", out_dir=tmp_path)
    assert engine == "matplotlib" and seen["fmt"] == "mp4"
    with pytest.raises(RuntimeError, match="cairo"):
        A.animate({}, "StackScene", "mp4", engine="manim", out_dir=tmp_path)


def test_render_rejects_unknown_scene():
    with pytest.raises(ValueError, match="scene"):
        render_animation({}, "Nope")


def test_real_render_both_scenes(tmp_path):
    pytest.importorskip("manim")
    data = animation_data(DeviceConfig(), biases=(0.0, -10.0, -30.0))
    assert {f["solver"] for f in data["frames"]} == {"batched"}
    t0 = time.perf_counter()
    for scene in SCENES:
        path, wall = render_animation(data, scene, out_dir=tmp_path / scene)
        assert path.suffix == ".mp4" and path.stat().st_size > 10_000, scene
        assert wall < 120
    assert time.perf_counter() - t0 < 180
