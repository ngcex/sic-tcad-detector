"""etna.viz figure layer: every kind builds, exports are non-empty and exactly sized,
provenance footer, style does not leak; carrier-map parity vs devsim; frame animations;
the two validation variants (real solves, pinned numbers).
"""

from __future__ import annotations

import re
import struct

import matplotlib as mpl
import numpy as np
import pytest

import etna
from etna import DeviceConfig, MeshData, SimResult
from etna.viz import (
    KINDS,
    config_hash,
    export_figure,
    figure_for,
    kinds_for,
    save_figure,
)
from etna.viz.style import SIZES_MM, size_inches

CFG = DeviceConfig()
V = np.linspace(0.0, -60.0, 13)


def _cv():
    W = np.minimum(1.7 + 0.8 * np.abs(V), 9.9) * 1e-4
    return SimResult(
        CFG,
        "cv",
        V,
        1e-12 / W,
        dict(
            one_over_C_squared=(W / 1e-12) ** 2,
            depletion_widths=W,
            W_from_C=W * 1.01,
            solver="batched",
            reach_through_V=3.1,
            cross_checked_to_V=40.0,
            validated_domain=True,
            area_cm2=0.04,
        ),
    )


def _field():
    x = np.linspace(0, 11.3, 200)
    E = 9e4 * np.exp(-x / 3)
    return SimResult(
        CFG,
        "field",
        x,
        E,
        dict(
            bias_V=-30.0,
            potential=np.linspace(0, 30, 200),
            net_doping=np.ones(200),
            solver="devsim",
        ),
    )


def _map(carriers=True):
    xs, ys = np.linspace(0, 50e-4, 6), np.linspace(0, 11.3e-4, 30)
    X, Y = np.meshgrid(xs, ys)
    x, y = X.ravel(), Y.ravel()
    nv = {
        "NetDoping": np.where(y < 0.3e-4, -1e19, 1e15),
        "Potential": y * 3e4,
        "ElectricField": 3e4 + 0 * y,
    }
    if carriers:
        nv["Electrons"] = np.exp(-y * 1e4) * 1e15
        nv["Holes"] = np.exp(y * 1e4 - 11) * 1e15
    mesh = MeshData(x, y, nv, [], [{"name": "anode", "position": 0.0}])
    return SimResult(
        DeviceConfig(half_width_um=50.0),
        "field",
        np.array([]),
        np.array([]),
        dict(bias_V=-30.0, solver="engine2d"),
        mesh,
    )


def _simple(sim_type, **md):
    x = np.array([20.0, 50.0, 100.0, 230.0])
    return SimResult(
        CFG,
        sim_type,
        x,
        np.array([0.99, 0.98, 0.97, 0.96]),
        dict(solver="batched", **md),
    )


INPUTS = {
    "cv": _cv,
    "field": _field,
    "map": _map,
    "cce": lambda: SimResult(
        CFG,
        "cce",
        V,
        np.ones_like(V),
        dict(solver="batched", reach_through_V=3.1, I_collected=V, I_generated=1.0),
    ),
    "tilted": lambda: SimResult(
        CFG,
        "cce",
        V,
        np.ones_like(V),
        dict(
            solver="batched",
            theta_deg=65.0,
            cce_vs_E0=0.97 * np.ones_like(V),
            generation_shape="bragg",
        ),
    ),
    "dark": lambda: SimResult(
        CFG,
        "dark_current",
        V,
        1e-11 * np.ones_like(V),
        dict(
            solver="batched",
            past_reach_through=np.abs(V) > 3.1,
            cross_checked_to_V=40.0,
        ),
    ),
    "dark_T": lambda: SimResult(
        CFG,
        "dark_current_vs_T",
        np.array([250.0, 300.0]),
        np.array([1e-13, 1e-11]),
        dict(V_bias=-20),
    ),
    "flash": lambda: _simple("flash", V_bias=-30, E_MeV=62),
    "transient": lambda: _simple("transient", V_bias=-30, transient_cce_raw=np.ones(4)),
    "temperature": lambda: SimResult(
        CFG,
        "temperature",
        np.array([250.0, 300.0, 350.0]),
        np.ones(3),
        dict(voltage_V=-30, method="hecht"),
    ),
    "damage": lambda: SimResult(
        CFG,
        "damage",
        np.array([1e11, 1e12, 1e14]),
        np.array([1.0, 0.99, np.nan]),
        dict(
            valid=[True, True, False],
            reason=["", "", "full donor compensation"],
            compensation_fluence=1.35e13,
            V_bias=-40,
            energy_MeV=5.6,
        ),
    ),
    "micro": lambda: SimResult(
        CFG,
        "microdosimetry",
        np.logspace(-2, 3, 50),
        np.ones(50),
        dict(y_F=1.0, y_D=2.0),
    ),
    "let": lambda: {
        "LET_keV_um": [0.5, 5.0, 50.0],
        "CCE": [0.99997, 0.99997, 0.99996],
        "converged": [True, True, False],
        "geometry": {"half_width_um": 50.0, "epi_um": 10.0},
        "bias_V": 50.0,
    },
    "config": lambda: CFG,
}


def _all_cases():
    for name, make in INPUTS.items():
        for kind in kinds_for(make()):
            yield name, kind


@pytest.mark.parametrize("name,kind", list(_all_cases()))
def test_every_kind_builds_and_exports(name, kind):
    fig = figure_for(INPUTS[name](), kind)
    assert fig.etna_kind == kind
    png = export_figure(fig, "png", dpi=72)
    assert png[:8] == b"\x89PNG\r\n\x1a\n" and len(png) > 2000


def test_kinds_cover_the_requested_observables():
    seen = {k for _, k in _all_cases()}
    for k in (
        "cv",
        "field",
        "cce",
        "cce_let",
        "dark_current",
        "flash",
        "transient",
        "map_potential",
        "map_field",
        "map_electrons",
        "stack",
    ):
        assert k in seen
    assert seen <= set(KINDS)


def _png_size(b):
    return struct.unpack(">II", b[16:24])


@pytest.mark.parametrize("size", ["single", "double", (120.0, 80.0)])
@pytest.mark.parametrize("dpi", [150, 300])
def test_png_pixel_size_is_exact(size, dpi):
    w, h = size_inches(size)
    png = export_figure(figure_for(_cv()), "png", dpi=dpi, size=size)
    assert _png_size(png) == (int(w * dpi), int(h * dpi))  # Agg truncates


def test_svg_and_pdf_page_size_and_editable_text(tmp_path):
    fig = figure_for(_cv(), "cv", size="double")
    w_pt, h_pt = (v * 72 / 25.4 for v in SIZES_MM["double"])
    svg = export_figure(fig, "svg").decode()
    m = re.search(r'<svg[^>]*width="([\d.]+)pt" height="([\d.]+)pt"', svg)
    assert (
        m
        and float(m.group(1)) == pytest.approx(w_pt, abs=0.1)
        and float(m.group(2)) == pytest.approx(h_pt, abs=0.1)
    )
    assert (
        "<text" in svg and config_hash(CFG) in svg
    )  # svg.fonttype none: real text, footer
    pdf = export_figure(fig, "pdf")
    assert pdf[:5] == b"%PDF-"
    box = re.search(rb"/MediaBox \[ *0 0 ([\d.]+) ([\d.]+) *\]", pdf)
    assert box and float(box.group(1)) == pytest.approx(w_pt, abs=0.1)
    assert float(box.group(2)) == pytest.approx(h_pt, abs=0.1)
    p = save_figure(fig, tmp_path / "cv.pdf")
    assert p.stat().st_size > 1000


def test_footer_carries_provenance():
    fig = figure_for(_cv())
    text = "\n".join(t.get_text() for t in fig.texts)
    assert f"config {config_hash(CFG)}" in text and "solver: batched" in text
    assert (
        f"ETNA {etna.__version__}" in text
        and "not validated against measurement" in text
    )
    assert re.search(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2} UTC", text)


def test_validity_bands_drawn():
    fig = figure_for(_cv(), "cv")
    labels = fig.axes[0].get_legend_handles_labels()[1]
    assert any("not cross-checked" in s for s in labels)
    assert any("full depletion" in s for s in labels)
    out = SimResult(
        CFG,
        "cv",
        V,
        np.ones_like(V),
        dict(
            _cv().metadata, validated_domain=False, domain_note="N_D_bulk outside box"
        ),
    )
    text = "\n".join(t.get_text() for t in figure_for(out).texts)
    assert "OUTSIDE the engine/devsim cross-check box" in text


def test_honesty_notes():
    for name, phrase in [
        ("dark", "budget estimate"),
        ("flash", "sensitivity bound"),
        ("damage", "absolute SiC fluence scale unvalidated"),  # SR-NIEL since B-07
        ("let", "regeneration pending"),
    ]:
        fig = figure_for(INPUTS[name]())
        assert phrase in " ".join(
            " ".join(t.get_text() for t in fig.texts).split()
        ), name


def test_config_hash_is_stable_and_sensitive():
    assert config_hash(DeviceConfig()) == config_hash(DeviceConfig())
    assert config_hash(DeviceConfig()) != config_hash(DeviceConfig(T=310.0))
    assert len(config_hash(etna.measured_cv_tabulated_config())) == 12  # tuple knots


def test_bad_inputs_rejected():
    with pytest.raises(ValueError, match="unknown kind"):
        figure_for(_cv(), "nope")
    with pytest.raises(ValueError, match="does not apply"):
        figure_for(_cv(), "map_field")
    with pytest.raises(ValueError, match="fmt"):
        export_figure(figure_for(_cv()), "jpg")
    with pytest.raises(ValueError, match="size"):
        figure_for(_cv(), size="huge")
    assert "map_electrons" not in kinds_for(_map(carriers=False))


def test_style_does_not_leak_into_global_rcparams():
    before = dict(mpl.rcParams)
    export_figure(figure_for(_cv()), "svg")
    assert mpl.rcParams["svg.fonttype"] == before["svg.fonttype"]
    assert mpl.rcParams["xtick.direction"] == before["xtick.direction"]


def test_let_table_loader_records_source():
    from etna.viz import load_let_table

    t = load_let_table("data/cce_let_table_100um.json")
    assert t["_source"].startswith("cce_let_table_100um.json sha ")
    text = "\n".join(x.get_text() for x in figure_for(t).texts)
    assert "cce_let_table_100um.json" in text


# ----------------------------------------------------------------- real-solver checks


def test_carrier_map_matches_devsim():
    """engine2d n, p vs devsim Electrons / Holes on the same 2D device and bias."""
    import devsim

    from etna.api.device import build_device
    from etna.core.devsim_reset import delete_device_and_mesh, reset_devsim_fully
    from etna.core.drift_diffusion import ramp_bias
    from etna.viz.maps import CARRIER_FLOOR_CM3, carrier_state_2d

    cfg = DeviceConfig(half_width_um=25.0)
    r = carrier_state_2d(cfg, -10.0)
    assert r.metadata["solver"] == "engine2d"
    reset_devsim_fully()
    di = build_device(cfg)
    try:
        ramp_bias(di, V_target=10.0, contact="cathode", V_step=0.5, V_step_max=2.0)
        get = lambda n: np.array(
            devsim.get_node_model_values(
                device=di["device_name"], region=di["region_name"], name=n
            )
        )  # noqa: E731
        for ours, theirs in (("Electrons", "Electrons"), ("Holes", "Holes")):
            a, b = r.mesh.node_values[ours], get(theirs)
            m = b >= CARRIER_FLOOR_CM3
            assert m.sum() > 500
            assert np.max(np.abs(a[m] - b[m]) / b[m]) < 1e-10
        assert (
            np.max(np.abs(r.mesh.node_values["Potential"] - get("Potential"))) < 1e-10
        )
    finally:
        delete_device_and_mesh(di["device_name"])
    assert "map_electrons" in kinds_for(r)


def test_carrier_map_unavailable_is_explicit():
    from etna.viz.maps import CarrierMapUnavailable, carrier_state_2d

    with pytest.raises(CarrierMapUnavailable, match="2D device"):
        carrier_state_2d(DeviceConfig(), -10.0)
    with pytest.raises(CarrierMapUnavailable, match="carrier map unavailable"):
        carrier_state_2d(DeviceConfig(half_width_um=25.0, T=600.0), -10.0)


def test_w_creep_validation_pinned():
    """1D graded epi + n-/n+ transition W creep vs the Photons Fig. 4b labels.

    Pins the 2026-10-20 numbers of the default stack (a FIT to that C-V curve, so this is a
    regression of the fit, not an independent validation): eps*A/C 1.655 / 9.519 / 9.750
    um at 0 / -10 / -30 V (labels 1.68 / 9.52 / 9.73), creep -10 -> -30 V +0.232 um
    (labels +0.21). Tolerance 0.005 um = the labels' rounding.
    """
    from etna.viz.validation import w_creep

    r = w_creep()
    assert r.metadata["solver"] == "batched"
    np.testing.assert_allclose(
        r.metadata["W_C_at_labels_um"], [1.655, 9.519, 9.750], atol=5e-3
    )
    assert np.max(np.abs(r.metadata["residual_C_um"])) < 0.03
    assert r.metadata["creep_C_um"] == pytest.approx(0.232, abs=5e-3)
    assert r.metadata["creep_label_um"] == pytest.approx(0.21)
    assert "not an independent validation" in r.metadata["verdict"]
    assert figure_for(r).etna_kind == "w_creep"


def test_tilted_alpha_shape_not_reproduced():
    """Microdosimeter preset, run_cce_tilted_alpha defaults, vs digitised Micro Fig. 4b.

    Shape only. Documents the known non-reproduction (USAGE section 8): the model is flat
    (CCE(0 V)/CCE(-50 V) ~0.994, rise -20 -> -50 V ~1e-4) where the data drop to 0.908 at
    0 V and rise ~1 % from -20 to -50 V.
    """
    from etna.viz.validation import tilted_alpha_shape

    r = tilted_alpha_shape()
    md = r.metadata
    assert md["solver"] == "batched" and md["reach_through_V"] == pytest.approx(
        9.0, abs=0.1
    )
    assert list(md["data_V"]) == [0, -5, -10, -20, -30, -40, -50]
    # 0.9928 (c-axis electron default of 2026-10-21 and legacy_perp alike)
    assert md["model_at_data"][0] == pytest.approx(0.9928, abs=1e-3)
    assert abs(md["rise_20_50_model"]) < 2e-3 < md["rise_20_50_data"]
    assert md["rms"] > 0.02  # not reproduced
    assert "NOT reproduced" in md["verdict"]
    assert figure_for(r).etna_kind == "tilted_alpha_shape"


# ----------------------------------------------------------------- animations


def _fake_data(monkeypatch):
    from tests.test_viz import _fake_cv, _fake_field

    monkeypatch.setattr(etna, "run_cv", _fake_cv)
    monkeypatch.setattr(etna, "run_field", _fake_field)
    from etna.viz import animation_data

    return animation_data(CFG, biases=(0.0, -10.0, -30.0))


@pytest.mark.parametrize("scene", ["StackScene", "DepletionScene"])
def test_frame_animation_gif_and_mp4(monkeypatch, tmp_path, scene):
    from PIL import Image

    from etna.viz import ffmpeg_available, render_frames

    d = _fake_data(monkeypatch)
    assert d["config_hash"] == config_hash(CFG) and "ETNA" in d["provenance"]
    gif, _ = render_frames(
        d, scene, "gif", out_dir=tmp_path, size_px=(320, 180), dpi=60
    )
    with Image.open(gif) as im:
        assert im.size == (320, 180) and im.n_frames >= 3
    if not ffmpeg_available()[0]:
        pytest.skip("ffmpeg not installed")
    mp4, _ = render_frames(
        d, scene, "mp4", out_dir=tmp_path, size_px=(320, 180), dpi=60
    )
    blob = mp4.read_bytes()
    assert blob[4:8] == b"ftyp" and len(blob) > 2000


def test_animate_falls_back_without_manim(monkeypatch, tmp_path):
    import etna.viz.animations as A

    d = _fake_data(monkeypatch)
    monkeypatch.setattr(A, "manim_available", lambda: (False, "no manim"))
    monkeypatch.setattr(
        A, "render_animation", lambda *a, **k: pytest.fail("manim called")
    )
    path, _, engine = A.animate(d, "StackScene", "gif", out_dir=tmp_path)
    assert engine == "matplotlib" and path.suffix == ".gif"
    with pytest.raises(ValueError, match="mp4 only"):
        A.animate(d, "StackScene", "gif", engine="manim")
