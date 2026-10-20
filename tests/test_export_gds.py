"""Reference layout (etna.export.gds) and its parameter file (etna.export.layout_params).

No devsim; gdstk optional (readback tests skip without it). The default parameter file
is checked against deliverables/DIMENSIONS-provenance.md (non-PLACEHOLDER rows must be
that file's rows), and no number of the simulation model's guard-ring mesh defaults may
appear in the layout code, the parameter file or any label.
"""

from __future__ import annotations

import copy
import re
import tomllib
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from etna.export import gds
from etna.export import layout_params as lp

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "deliverables" / "DIMENSIONS-provenance.md"



def _polygon_area_um2(points) -> float:
    """Shoelace area (um^2) of a simple / keyhole polygon (test oracle for the layout)."""
    s = 0.0
    n = len(points)
    for i in range(n):
        x0, y0 = points[i]
        x1, y1 = points[(i + 1) % n]
        s += x0 * y1 - x1 * y0
    return abs(s) / 2.0


# The simulation model's guard-ring mesh defaults (written here only as regexes).
FORBIDDEN = [
    re.compile(r"(?<![\d.])(50|5|3|1)(\.0+)?\s*(um|µm|micron)", re.I),
    re.compile(r"(?<![\d.])5(\.0*)?\s*[eE]\+?18"),
    re.compile(r"5\s*x\s*10\^?\s*18"),
    re.compile(r"50\s*/\s*5\s*/\s*3\s*/\s*1"),
]


def _doc_rows() -> dict[str, dict[str, str]]:
    rows = {}
    for line in DOC.read_text().splitlines():
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 4 or cells[0] == "item" or set(cells[0]) <= {"-"}:
            continue
        rows[cells[0]] = dict(value=cells[1], cls=cells[2], source=cells[3])
    return rows


def _no_forbidden(text: str, where: str) -> None:
    for rx in FORBIDDEN:
        m = rx.search(text)
        assert m is None, f"mesh-default number {m.group(0)!r} in {where}"


def _raw() -> dict:
    return tomllib.loads(lp.DEFAULT_PARAMS_PATH.read_text())


def _all_real(raw: dict, cls="DOCUMENTED") -> dict:
    """Copy of the default mapping with every PLACEHOLDER re-classed (test only)."""
    raw = copy.deepcopy(raw)
    for dev in lp.DEVICES:
        for t in raw[dev].values():
            if t["class"] == "PLACEHOLDER":
                t["class"] = cls
                t["source"] = "test mask set (fictitious, test only)"
    return raw


# --------------------------------------------------------------------------- file


def test_default_file_valid_and_seeded_from_provenance_doc():
    params = lp.load_layout_params()
    rows = _doc_rows()
    assert len(rows) >= 15
    for dev in lp.DEVICES:
        for name, d in params.devices[dev].items():
            if d.is_placeholder:
                assert "PLACEHOLDER" in d.source, (dev, name)
                continue
            # non-placeholder: must be a row of the provenance doc, same class / source
            assert d.doc_item in rows, (dev, name)
            r = rows[d.doc_item]
            assert r["cls"].replace(" ", "_").startswith(d.cls), (dev, name)
            assert d.source == r["source"], (dev, name)
            assert f"{d.value:g}" in r["value"], (dev, name)
    dos = params.devices["dosimeter"]
    assert dos["active_side"].value == 2.0 and dos["active_side"].cls == "DOCUMENTED"
    assert dos["die_side"].value == 2.7 and dos["die_side"].uncertainty == 0.15
    assert dos["die_side"].cls == "MEASURED_FROM_FIGURE"
    mic = params.devices["microdosimeter"]
    assert mic["active_side"].value == 5.0
    # unknown per the doc -> PLACEHOLDER
    for name in ("guard_ring_gap", "guard_ring_width", "pad_size", "dicing_street"):
        assert dos[name].is_placeholder and mic[name].is_placeholder
    assert mic["die_side"].is_placeholder and dos["oxide_thickness"].is_placeholder
    # measured periphery 0.2-0.4 mm per side holds for the dosimeter
    band = (dos["die_side"].um - dos["active_side"].um) / 2e3
    assert 0.2 <= band <= 0.4


def test_sources_free_of_mesh_defaults():
    for f in (
        Path(gds.__file__),
        Path(lp.__file__),
        lp.DEFAULT_PARAMS_PATH,
    ):
        _no_forbidden(f.read_text(), str(f.relative_to(ROOT)))


def test_env_var_and_argument_select_the_file(tmp_path, monkeypatch):
    raw = _raw()
    raw_text = lp.DEFAULT_PARAMS_PATH.read_text()
    alt = tmp_path / "alt.toml"
    alt.write_text(
        raw_text.replace('value = 100.0\nunit = "um"', 'value = 80.0\nunit = "um"')
    )
    monkeypatch.setenv(lp.ENV_VAR, str(alt))
    assert lp.load_layout_params().path == str(alt)
    assert lp.load_layout_params().dim("dosimeter", "pad_size").value == 80.0
    # an explicit argument wins over the env var
    assert (
        lp.load_layout_params(lp.DEFAULT_PARAMS_PATH).dim("dosimeter", "pad_size").value
        == raw["dosimeter"]["pad_size"]["value"]
    )


def _bad(mutate, match):
    raw = _raw()
    mutate(raw)
    with pytest.raises(lp.LayoutParamsError, match=match):
        lp.parse_layout_params(raw)


def test_invalid_files_rejected_with_field_names(tmp_path):
    _bad(
        lambda r: r["dosimeter"]["guard_ring_width"].update(value=-2.0),
        r"dosimeter\.guard_ring_width\.value: must be > 0",
    )
    _bad(
        lambda r: r["dosimeter"]["pad_size"].update(unit="mil"),
        r"dosimeter\.pad_size\.unit: must be one of",
    )
    _bad(
        lambda r: r["microdosimeter"]["die_side"].update({"class": "GUESS"}),
        r"microdosimeter\.die_side\.class",
    )
    _bad(
        lambda r: r["dosimeter"]["die_side"].pop("source"),
        r"dosimeter\.die_side: missing key 'source'",
    )
    _bad(
        lambda r: r["dosimeter"].pop("dicing_street"),
        r"dosimeter\.dicing_street: missing dimension",
    )
    _bad(
        lambda r: r["dosimeter"].update(guard_ring_widht=r["dosimeter"]["pad_size"]),
        r"dosimeter\.guard_ring_widht: unknown dimension",
    )
    _bad(
        lambda r: r["dosimeter"]["guard_ring_count"].update(value=1.5),
        r"guard_ring_count\.value: must be an integer",
    )
    _bad(
        lambda r: r["dosimeter"]["guard_ring_gap"].update(uncertainty=-1),
        r"guard_ring_gap\.uncertainty",
    )
    _bad(lambda r: r.update(schema_version=2), "schema_version")
    # geometry: rings must fit inside the die, the pad between ring and die edge
    _bad(
        lambda r: r["dosimeter"]["guard_ring_gap"].update(value=0.4, unit="mm"),
        r"dosimeter\.guard_ring_\*: guard ring\(s\) do not fit inside the die",
    )
    _bad(
        lambda r: r["dosimeter"]["pad_size"].update(value=0.3, unit="mm"),
        r"dosimeter\.pad_size: pad .* does not fit",
    )
    _bad(
        lambda r: r["microdosimeter"]["die_side"].update(value=4.0),
        r"microdosimeter\.die_side: die .* must be larger",
    )
    bad = tmp_path / "broken.toml"
    bad.write_text("schema_version = 1\n[dosimeter.active_side\n")
    with pytest.raises(lp.LayoutParamsError, match="not valid TOML"):
        lp.load_layout_params(bad)
    with pytest.raises(lp.LayoutParamsError, match="not found"):
        lp.load_layout_params(tmp_path / "missing.toml")


def test_unused_ring_spacing_is_not_a_blocking_placeholder():
    params = lp.load_layout_params()
    assert "dosimeter.guard_ring_spacing" not in params.placeholders("dosimeter")
    raw = _raw()
    raw["dosimeter"]["guard_ring_count"]["value"] = 2
    two = lp.parse_layout_params(raw)
    assert "dosimeter.guard_ring_spacing" in two.placeholders("dosimeter")


# --------------------------------------------------------------------------- geometry


def test_placeholder_warning_and_label():
    with pytest.warns(lp.PlaceholderDimensionWarning, match="guard_ring_gap"):
        lay = gds.build_layout("dosimeter")
    assert lay.cell_name.endswith("REFERENCE-NOT-FAB")
    assert "dosimeter.guard_ring_gap" in lay.placeholders
    note = [n.text for n in lay.notes if "PLACEHOLDER dimensions" in n.text]
    assert note and "guard_ring_gap" in note[0] and "pad_size" in note[0]
    assert any(
        n.layer == "DISCLAIMER" and gds.REFERENCE_LABEL in n.text for n in lay.notes
    )
    with pytest.warns(lp.PlaceholderDimensionWarning):
        svg = gds.layout_svg("microdosimeter")
    assert "die_side" in svg and "PLACEHOLDER dimensions" in svg


def test_no_placeholder_left_drops_not_fab_and_warning():
    params = lp.parse_layout_params(_all_real(_raw()))
    assert params.placeholders() == []
    with warnings.catch_warnings():
        warnings.simplefilter("error", lp.PlaceholderDimensionWarning)
        lay = gds.build_layout("dosimeter", params)
    assert lay.cell_name == "ETNA_DOSIMETER_2x2mm_REFERENCE"
    assert "REFERENCE-NOT-FAB" not in lay.cell_name
    assert any(gds.REFERENCE_LABEL in n.text for n in lay.notes)  # always labelled


def test_layout_geometry_follows_parameters():
    params = lp.load_layout_params()
    lay = gds.build_layout("dosimeter", params)
    by = {}
    for s in lay.shapes:
        by.setdefault(s.layer, []).append(s)
        assert s.label.split(" | ")[0] in gds.CLASSES and "src: " in s.label
    d = params.devices["dosimeter"]
    a, b = d["active_side"].um, d["die_side"].um
    assert _polygon_area_um2(by["ACTIVE"][0].points_um) == pytest.approx(a * a)
    assert _polygon_area_um2(by["ACTIVE"][0].points_um) / 1e6 == pytest.approx(4.0)
    assert _polygon_area_um2(by["DIE_OUTLINE"][0].points_um) == pytest.approx(b * b)
    assert len(by["GUARD_RING"]) == int(d["guard_ring_count"].value)
    r_in = a / 2 + d["guard_ring_gap"].um
    r_out = r_in + d["guard_ring_width"].um
    assert _polygon_area_um2(by["GUARD_RING"][0].points_um) == pytest.approx(
        (2 * r_out) ** 2 - (2 * r_in) ** 2
    )
    assert _polygon_area_um2(by["PAD"][0].points_um) == pytest.approx(
        d["pad_size"].um ** 2
    )
    s = d["dicing_street"].um / 2
    assert _polygon_area_um2(by["DICING_STREET"][0].points_um) == pytest.approx(
        (b + 2 * s) ** 2 - b * b
    )
    mic = gds.build_layout("microdosimeter", params)
    act = next(x for x in mic.shapes if x.layer == "ACTIVE")
    assert _polygon_area_um2(act.points_um) / 1e6 == pytest.approx(25.0)
    with pytest.raises(ValueError):
        gds.build_layout("pixel", params)


def test_three_rings_and_mapping_input():
    raw = _raw()
    raw["dosimeter"]["guard_ring_count"]["value"] = 3
    lay = gds.build_layout("dosimeter", raw)  # a mapping is accepted too
    rings = [s for s in lay.shapes if s.layer == "GUARD_RING"]
    assert len(rings) == 3 and "spacing" in rings[1].label


@pytest.mark.parametrize("device", gds.DEVICES)
def test_svg_preview_well_formed(device):
    svg = gds.layout_svg(device)
    root = ET.fromstring(svg)
    assert root.tag.endswith("svg")
    assert "REFERENCE-NOT-FAB" in svg
    assert gds.REFERENCE_LABEL in svg and gds.VALIDATION_LABEL in svg
    lay = gds.build_layout(device)
    assert svg.count("<polygon") == len(lay.shapes)
    text = " ".join(re.findall(r">([^<]*)</text>", svg))
    _no_forbidden(text, f"{device} SVG text")


# --------------------------------------------------------------------------- GDSII


@pytest.fixture(scope="module")
def readback(tmp_path_factory):
    gdstk = pytest.importorskip("gdstk")
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", lp.PlaceholderDimensionWarning)
        for device in gds.DEVICES:
            path = tmp_path_factory.mktemp("gds") / f"{device}.gds"
            data = gds.layout_gds(device)
            assert data[:4] == b"\x00\x06\x00\x02"  # HEADER record
            path.write_bytes(data)
            out[device] = gdstk.read_gds(str(path))
    return gdstk, out


@pytest.mark.parametrize("device", gds.DEVICES)
def test_gds_top_cell_and_disclaimer(readback, device):
    _, libs = readback
    lib = libs[device]
    tops = lib.top_level()
    assert len(tops) == 1
    top = tops[0]
    assert top.name.endswith("REFERENCE-NOT-FAB")
    assert top.name == gds.build_layout(device).cell_name
    dl, dt = gds.LAYERS["DISCLAIMER"]
    assert any(
        lb.text == gds.DISCLAIMER and lb.layer == dl and lb.texttype == dt
        for lb in top.labels
    )
    assert gds.DISCLAIMER.startswith(gds.REFERENCE_LABEL)
    assert any(gds.VALIDATION_LABEL in lb.text for lb in top.labels)
    assert any("PLACEHOLDER dimensions" in lb.text for lb in top.labels)
    assert lib.unit == pytest.approx(1e-6, rel=1e-12)


@pytest.mark.parametrize("device", gds.DEVICES)
def test_gds_every_polygon_labelled(readback, device):
    gdstk, libs = readback
    top = libs[device].top_level()[0]
    assert top.polygons and not top.paths and not top.references
    for poly in top.polygons:
        own = [
            lb
            for lb in top.labels
            if lb.layer == poly.layer
            and lb.texttype == poly.datatype
            and gdstk.inside([lb.origin], poly)[0]
        ]
        assert own, f"unlabelled polygon on layer {poly.layer}"
        assert any(
            lb.text.split(" | ")[0] in gds.CLASSES and "src: " in lb.text for lb in own
        )


@pytest.mark.parametrize("device", gds.DEVICES)
def test_gds_round_trip_geometry(readback, device):
    _, libs = readback
    top = libs[device].top_level()[0]
    lay = gds.build_layout(device)
    assert len(top.polygons) == len(lay.shapes)
    got = sorted(((p.layer, p.datatype), round(p.area(), 3)) for p in top.polygons)
    want = sorted(
        (gds.LAYERS[s.layer], round(_polygon_area_um2(s.points_um), 3))
        for s in lay.shapes
    )
    assert got == want
    for lb in top.labels:
        _no_forbidden(lb.text, f"{device} GDS label")


def test_changing_the_file_changes_the_gds(tmp_path):
    gdstk = pytest.importorskip("gdstk")
    raw_text = lp.DEFAULT_PARAMS_PATH.read_text()
    edited = tmp_path / "edited.toml"
    # dosimeter ring width 20 um -> 35 um (first occurrence is the dosimeter's)
    old = "[dosimeter.guard_ring_width]\nvalue = 20.0"
    assert old in raw_text
    edited.write_text(
        raw_text.replace(old, "[dosimeter.guard_ring_width]\nvalue = 35.0")
    )

    def ring_area(path):
        f = tmp_path / (Path(path).stem + ".gds")
        f.write_bytes(gds.layout_gds("dosimeter", path))
        top = gdstk.read_gds(str(f)).top_level()[0]
        rings = [
            p for p in top.polygons if (p.layer, p.datatype) == gds.LAYERS["GUARD_RING"]
        ]
        assert len(rings) == 1
        lab = [lb.text for lb in top.labels if lb.layer == gds.LAYERS["GUARD_RING"][0]]
        return rings[0].area(), lab[0]

    a0, l0 = ring_area(lp.DEFAULT_PARAMS_PATH)
    a1, l1 = ring_area(edited)
    act = 2000.0 / 2 + 150.0  # ring inner half-width (active/2 + gap), um
    assert a0 == pytest.approx((2 * (act + 20.0)) ** 2 - (2 * act) ** 2)
    assert a1 == pytest.approx((2 * (act + 35.0)) ** 2 - (2 * act) ** 2)
    assert "guard_ring_width 20 um" in l0 and "guard_ring_width 35 um" in l1
