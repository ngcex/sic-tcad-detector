"""etna.export: spec sheet, cross-section, mesh/profile exports (pure, synthetic data).

Synthetic Reports / SimResults keep these tests solver-free; one test runs the real
2D engine (no devsim equations) for the tensor-grid VTK path.
"""

from __future__ import annotations

import csv
import io
import json

import numpy as np
import pytest

import etna
from etna import DeviceConfig, MeshData, SimResult
from etna.api.results import peak_field_1d
from etna.api.run_all import Report, SectionResult
from etna.export import (
    VALIDATION_LABEL,
    build_spec,
    cross_section_png,
    cross_section_svg,
    profile_to_csv,
    profile_to_vtk,
    spec_to_csv,
    spec_to_json,
    spec_to_markdown,
    spec_to_pdf,
    stack_layers,
)
from etna.export.labels import analytic_vfd
from etna.export.layout_params import PlaceholderDimensionWarning, load_layout_params


LABEL = (
    "simulation output; fast engines cross-checked against the devsim reference solver; "
    "not validated against measurement except where stated"
)


def _sec(name, res, solver="batched"):
    return SectionResult(name, name, "ok", solver=solver, result=res)


def _fake_report(cfg, w30_um=9.9788):  # eps*A/C(-30 V) of the default stack
    bias = np.linspace(0.0, -50.0, 51)
    W = np.interp(-bias, [0, 10, 30, 50], [1.7, 9.49, w30_um, 9.93]) * 1e-4
    md = dict(
        depletion_widths=W,
        solver="batched",
        validated_domain=True,
        domain_note=None,
        reach_through_V=10.3,
        reach_through_known=True,
        reach_through_method="engine dark sweep, 0.05 V resolution",
    )
    cv = SimResult(cfg, "cv", bias, 9.7 * 8.854e-14 * cfg.area_cm2 / W, md)
    cb = np.linspace(-10, -40, 4)
    cce = SimResult(cfg, "cce", cb, np.array([0.999, 0.9995, 0.9997, 0.9998]), dict(md))
    db = np.linspace(0, -100, 11)
    dark = SimResult(
        cfg,
        "dark_current",
        db,
        # J(-60 V) = 4.5e-10 A/cm2 at any area (18 pA on the 4 mm2 default since
        # 2026-10-18; the fake was a fixed 4.5e-14 A = the old 1e-4 cm2 default)
        -4.5e-10 * cfg.area_cm2 * np.abs(db) / 60.0,
        dict(solver="devsim", validated_domain=True, domain_note=None),
    )
    x = np.array([0.0, 0.3, 0.5, 5.0, 10.29, 10.3, 10.31, 11.3])
    E = np.array([0.0, 8e5, 1.14e5, 3e4, 2.2e4, 9e5, 2.1e4, 0.0])
    field = SimResult(
        cfg,
        "field",
        x,
        E,
        dict(bias_V=-50.0, solver="devsim", validated_domain=True, domain_note=None),
    )
    secs = [
        _sec("cv", cv),
        _sec("cce", cce),
        _sec("field", field, "devsim"),
        _sec("dark", dark, "devsim"),
    ]
    return Report(cfg, secs, 1.0)


def test_label_is_verbatim():
    assert VALIDATION_LABEL == LABEL


def test_analytic_vfd_pinned_default_stack():
    # same V_bi inputs as the 2D sweep gate (etna.core.optimization): 2.74 V on the
    # measured-C-V default (to the transition start) since 2026-10-18; 9.93 V on the
    # previous default (8.21 V before the 2026-10-17 small-signal C recalibration)
    assert analytic_vfd(DeviceConfig()) == pytest.approx(2.736, abs=0.005)
    assert analytic_vfd(etna.nplus_calibrated_config()) == pytest.approx(
        9.93, abs=0.005
    )


def test_spec_default_reference_device():
    cfg = DeviceConfig()
    spec = build_spec(cfg, _fake_report(cfg))
    assert spec["label"] == LABEL and spec["reference_device"]
    fd = spec["full_depletion"]
    assert (
        fd["reach_through_V"] == 10.3
        and "Debye length" in fd["reach_through_definition"]
    )
    assert fd["analytic_V_fd_V"] == pytest.approx(2.736, abs=0.005)
    assert "transition" in fd["reach_through_definition"]
    assert [L["name"] for L in spec["layers"]] == [
        "p+ anode",
        "n- epi",
        "n-/n+ transition",
        "n+ substrate",
    ]
    assert [L["thickness_um"] for L in spec["layers"]] == [0.3, 8.863366, 1.136634, 1.0]
    assert "not a fabrication thickness" in spec["layers"][3]["note"]
    s = spec["epi_doping"]["samples"]
    assert s[0]["N_D_cm3"] == pytest.approx(cfg.N_D_junction)
    assert s[-1]["depth_below_junction_um"] == 10.0
    # the default is fitted on the whole digitised curve (2026-10-18): no out-of-fit point
    assert (
        "No out-of-fit C-V comparison" in spec["cv"]["systematic"]
        and "in-fit residual" in spec["cv"]["systematic"]
    )
    assert "calibration input" in spec["cv"]["calibration"]
    assert "not by SIMS" in spec["cv"]["calibration"]
    assert spec["cv"]["W_um"]["-30 V"] == pytest.approx(9.9788, abs=1e-4)
    assert "dQ/dV" in spec["cv"]["W_definition"]
    assert (
        "65 deg" in spec["cce"]["generation"]
        and "never reach the n+" in spec["cce"]["generation"]
    )
    # 4 mm2 default area since 2026-10-18: 1.8e-11 A (was 4.5e-14 A on 1e-4 cm2, x400)
    assert spec["dark_current"]["I_at_minus60V_A"] == pytest.approx(1.8e-11, rel=1e-12)
    assert spec["dark_current"]["J_at_minus60V_A_per_cm2"] == pytest.approx(
        4.5e-10, rel=1e-12
    )
    assert "calibration" in spec["dark_current"]["note"]
    # the spikes at both abrupt steps (8e5 at 0.3 um, 9e5 at 10.3 um) are excluded from the peak
    assert spec["field"]["peak_E_V_per_cm"] == pytest.approx(1.14e5)
    assert [p["solver"] for p in spec["provenance"]] == [
        "batched",
        "batched",
        "devsim",
        "devsim",
    ]


def test_spec_other_config_does_not_claim_measurement():
    cfg = DeviceConfig(epi_thickness_um=12.0)
    spec = build_spec(cfg, _fake_report(cfg))
    assert not spec["reference_device"]
    assert "not re-evaluated" in spec["cv"]["systematic"]
    assert "measured_W_minus30V_um" not in spec["cv"]


def test_spec_dark_calibration_only_for_calibrated_n_t():
    cfg = DeviceConfig(N_t=1e13)
    dk = build_spec(cfg, _fake_report(cfg))["dark_current"]
    assert dk["N_t_used"] == 1e13 and "does not apply" in dk["note"]
    md = spec_to_markdown(build_spec(cfg, _fake_report(cfg)))
    assert "cross-check box" in md and "not a validation against measurement" in md


def test_spec_legacy_definitions():
    cfg = etna.legacy_config()
    spec = build_spec(cfg, None)
    assert "last mesh interval" in spec["full_depletion"]["reach_through_definition"]
    assert "65 deg" not in spec["cce"]["generation"]
    assert [L["name"] for L in spec["layers"]] == ["p+ anode", "n- epi"]
    assert all(p["status"] == "not run" for p in spec["provenance"])


def test_spec_partial_2d_report():
    cfg = DeviceConfig(half_width_um=25.0)
    secs = [
        SectionResult("cv", "C-V", "skipped", reason="needs a 1D device"),
        SectionResult("dark", "Dark", "error", reason="RuntimeError: x"),
    ]
    spec = build_spec(cfg, Report(cfg, secs, 0.0))
    assert spec["dimension"] == "2D" and not spec["cv"]["available"]
    prov = {p["section"]: p for p in spec["provenance"]}
    assert prov["cv"]["status"] == "skipped" and "1D" in prov["cv"]["reason"]
    assert prov["dark"]["status"] == "error"
    for fmt in (spec_to_json, spec_to_csv, spec_to_markdown):
        assert LABEL in fmt(spec)


def test_spec_rejects_mismatched_report():
    with pytest.raises(ValueError, match="different DeviceConfig"):
        build_spec(DeviceConfig(), _fake_report(DeviceConfig(epi_thickness_um=12.0)))


def test_spec_formats_carry_label():
    cfg = DeviceConfig()
    spec = build_spec(cfg, _fake_report(cfg))
    assert json.loads(spec_to_json(spec))["label"] == LABEL
    rows = list(csv.reader(io.StringIO(spec_to_csv(spec))))
    assert rows[0] == ["# " + LABEL] and rows[1] == ["key", "value"]
    assert ["full_depletion.reach_through_V", "10.3"] in rows
    md = spec_to_markdown(spec)
    assert LABEL in md and "Reach-through: 10.30 V" in md and "9.93 V" in md
    pdf = spec_to_pdf(spec, cfg)
    assert pdf.startswith(b"%PDF") and len(pdf) > 5000


def _svg_text(svg: str) -> str:
    """All <text> contents joined (long notes are wrapped over several elements)."""
    import re

    return " ".join(re.findall(r">([^<]*)</text>", svg))


def test_cross_section_svg_to_scale_and_labelled():
    cfg = DeviceConfig()
    svg = cross_section_svg(cfg)
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert LABEL in svg
    txt = _svg_text(svg)
    assert "not to scale" in txt and "layout dimensions (dosimeter" in txt
    assert "die 2.7 +- 0.15 mm [MEASURED_FROM_FIGURE]" in txt
    assert "PLACEHOLDER dimensions (dosimeter)" in txt
    assert "0.3 µm" in svg and "10 µm" in svg and "1 µm" in svg
    import re
    import xml.etree.ElementTree as ET

    ET.fromstring(svg)  # well-formed
    hs = [
        float(h)
        for h in re.findall(
            r'<rect x="110.0" y="[^"]+" width="300.0" height="([^"]+)"', svg
        )
    ]
    assert (
        len(hs) == 4
    )  # p+ / n- epi / graded n-/n+ transition / n+ (default since 2026-10-18)
    np.testing.assert_allclose(
        np.array(hs) / sum(hs), np.array([0.3, 8.863366, 1.136634, 1]) / 11.3, rtol=1e-3
    )


def test_cross_section_png():
    png = cross_section_png(DeviceConfig(), dpi=60)
    assert png[:8] == b"\x89PNG\r\n\x1a\n"


def test_stack_layers_uniform_and_legacy():
    L = stack_layers(DeviceConfig(doping_profile="uniform", N_D=1e15))
    assert "uniform" in L[1]["doping"]
    assert len(stack_layers(etna.legacy_config())) == 2


def test_peak_field_helper_excludes_the_abrupt_steps():
    x = np.array([0.3, 1.0, 5.0, 10.3, 11.0])
    E = np.array([9.0, 7.0, 2.0, 8.0, 0.5])
    # default stack: p+ 0.3 um step and epi/n+ step (10.3) excluded
    assert peak_field_1d(x, E, DeviceConfig())[:2] == (7.0, 1.0)
    # legacy (p+ 1 um, no n+): only the step at 1.0 um is excluded
    assert peak_field_1d(x, E, etna.legacy_config())[:2] == (9.0, 0.3)


# ------------------------------------------------------------------ mesh / profile


def _fake_field_1d():
    x_cm = np.array([0.0, 0.3, 5.0, 11.3]) * 1e-4
    return SimResult(
        DeviceConfig(),
        "field",
        x_cm * 1e4,
        np.array([0.0, 1e5, np.nan, 0.0]),
        dict(bias_V=0.0, solver="devsim", validated_domain=True, domain_note=None),
        mesh=MeshData(
            x_cm,
            None,
            {
                "NetDoping": np.array([-1e18, 3e15, 8e13, 1e19]),
                "Potential": np.array([-1.5, -1.4, 18.0, 31.6]),
                "ElectricField": np.array([0.0, 1e5, np.nan, 0.0]),
            },
            [],
            [],
        ),
    )


def _parse_vtk(text):
    lines = text.splitlines()
    assert lines[0] == "# vtk DataFile Version 3.0" and lines[2] == "ASCII"
    out = {"title": lines[1]}
    i = 3
    while i < len(lines):
        tok = lines[i].split()
        if tok and tok[0] in ("X_COORDINATES", "Y_COORDINATES"):
            out[tok[0]] = np.array(lines[i + 1].split(), dtype=float)
            i += 2
        elif tok and tok[0] == "DIMENSIONS":
            out["dims"] = tuple(int(t) for t in tok[1:])
            i += 1
        elif tok and tok[0] == "SCALARS":
            n = int(np.prod(out["dims"]))
            vals, i = [], i + 2
            while len(vals) < n:
                vals += lines[i].split()
                i += 1
            out[tok[1]] = np.array(vals, dtype=float)
        else:
            i += 1
    return out


def test_profile_1d_csv_and_vtk():
    res = _fake_field_1d()
    text = profile_to_csv(res)
    assert text.splitlines()[0] == "# " + LABEL
    body = [r for r in csv.reader(io.StringIO(text)) if r and not r[0].startswith("#")]
    assert body[0] == [
        "x_um",
        "NetDoping_cm-3",
        "Potential_V",
        "ElectricField_V_per_cm",
    ]
    assert body[2][0] == "0.3" and body[3][3] == "nan"
    vtk = _parse_vtk(profile_to_vtk(res))
    assert vtk["title"].startswith(LABEL) and len(vtk["title"]) <= 256
    assert vtk["dims"] == (4, 1, 1)
    np.testing.assert_allclose(vtk["X_COORDINATES"], [0, 0.3, 5, 11.3])
    np.testing.assert_allclose(vtk["Potential_V"], [-1.5, -1.4, 18.0, 31.6])
    assert np.isnan(vtk["ElectricField_V_per_cm"][2])
    assert "bias_V 1 1 double\n0.0" in profile_to_vtk(res)  # 0 V is not lost


def test_profile_2d_engine_vtk_reorders_consistently():
    res = etna.run_field(DeviceConfig(half_width_um=25.0), bias_V=-30.0)
    assert res.metadata["solver"] == "engine2d"
    m = res.mesh
    vtk = _parse_vtk(profile_to_vtk(res))
    nx, ny, _ = vtk["dims"]
    assert nx * ny == m.x_coords.size
    X, Y = np.meshgrid(vtk["X_COORDINATES"], vtk["Y_COORDINATES"])  # x fastest
    # every exported value must sit at the node it belongs to
    lookup = {
        (round(a * 1e4, 9), round(b * 1e4, 9)): i
        for i, (a, b) in enumerate(zip(m.x_coords, m.y_coords))
    }
    idx = np.array(
        [lookup[(round(a, 9), round(b, 9))] for a, b in zip(X.ravel(), Y.ravel())]
    )
    np.testing.assert_array_equal(vtk["Potential_V"], m.node_values["Potential"][idx])
    np.testing.assert_array_equal(
        vtk["NetDoping_cm-3"], m.node_values["NetDoping"][idx]
    )
    text = profile_to_csv(res)
    assert "x = lateral, y = depth" in text and LABEL in text


def test_profile_2d_non_grid_falls_back_to_polydata():
    x = np.array([0.0, 1.0, 0.5]) * 1e-4
    y = np.array([0.0, 0.0, 1.0]) * 1e-4
    nv = {
        k: np.array([1.0, 2.0, 3.0])
        for k in ("NetDoping", "Potential", "ElectricField")
    }
    res = SimResult(
        DeviceConfig(half_width_um=1.0),
        "field",
        np.array([]),
        np.array([]),
        dict(bias_V=-1.0, solver="devsim"),
        MeshData(x, y, nv, [], []),
    )
    text = profile_to_vtk(res)
    assert "DATASET POLYDATA" in text and "POINTS 3 double" in text


def test_profile_requires_mesh():
    with pytest.raises(ValueError, match="no mesh"):
        profile_to_csv(SimResult(DeviceConfig(), "cv", np.array([]), np.array([])))


# ------------------------------------------------------------ layout parameter file


def test_spec_and_cross_section_read_the_layout_file(tmp_path):
    cfg = DeviceConfig()
    with pytest.warns(PlaceholderDimensionWarning, match="spec sheet"):
        spec = build_spec(cfg)
    lay = spec["layout"]
    assert lay["device"] == "dosimeter"
    assert lay["params_file"] == "data/foundry/layout_params.toml"
    assert "dosimeter.guard_ring_gap" in lay["placeholders"]
    names = {d["name"]: d for d in lay["dimensions"]}
    assert names["die_side"]["value_um"] == 2700.0
    assert names["guard_ring_gap"]["cls"] == "PLACEHOLDER"
    md = spec_to_markdown(spec)
    assert "## Layout dimensions" in md and "PLACEHOLDER dimensions (dosimeter)" in md
    json.loads(spec_to_json(spec))  # valid JSON
    rows = list(csv.reader(io.StringIO(spec_to_csv(spec))))
    assert ["layout.device", "dosimeter"] in rows
    with pytest.warns(PlaceholderDimensionWarning, match="cross-section"):
        cross_section_svg(cfg)
    # an alternative file changes what is quoted
    alt = tmp_path / "alt.toml"
    text = load_layout_params().path
    alt.write_text(
        open(text)
        .read()
        .replace(
            "[dosimeter.guard_ring_gap]\nvalue = 150.0",
            "[dosimeter.guard_ring_gap]\nvalue = 120.0",
        )
    )
    spec2 = build_spec(cfg, layout_params=alt)
    gap = {d["name"]: d for d in spec2["layout"]["dimensions"]}["guard_ring_gap"]
    assert gap["value"] == 120.0
    assert "ring gap from p+ edge 120 um" in _svg_text(
        cross_section_svg(cfg, layout_params=alt)
    )
    pdf = spec_to_pdf(spec2, cfg, layout_params=alt)
    assert pdf.startswith(b"%PDF")


def test_layout_section_without_matching_device():
    cfg = DeviceConfig(area_cm2=0.1)
    spec = build_spec(cfg)
    assert spec["layout"]["device"] is None and spec["layout"]["placeholders"] == []
    assert "no layout dimensions quoted" in spec["layout"]["note"]
    assert "lateral dimensions not quoted" in _svg_text(cross_section_svg(cfg))
    assert etna.microdosimeter_config().area_cm2 == 0.25
    assert build_spec(etna.microdosimeter_config())["layout"]["device"] == (
        "microdosimeter"
    )


def test_spec_pdf_refuses_to_guess_the_layout_file():
    import tomllib

    raw = tomllib.loads(open(load_layout_params().path).read())
    cfg = DeviceConfig()
    spec = build_spec(cfg, layout_params=raw)  # built from a mapping
    assert spec["layout"]["params_file"] == "<mapping>"
    with pytest.raises(ValueError, match="pass layout_params"):
        spec_to_pdf(spec, cfg)
    assert spec_to_pdf(spec, cfg, layout_params=raw).startswith(b"%PDF")
