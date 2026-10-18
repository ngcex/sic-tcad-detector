"""Schematic reference layout (etna.export.gds): no devsim, gdstk optional.

The DIMENSIONS table is checked against deliverables/DIMENSIONS-provenance.md (the only
allowed source of dimensions) so the layout and the doc cannot drift.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from etna.export import gds

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


# --------------------------------------------------------------------------- table


def test_table_matches_provenance_doc():
    rows = _doc_rows()
    assert len(rows) >= 15
    for d in gds.DIMENSIONS.values():
        assert d.doc_item in rows, d.doc_item
        r = rows[d.doc_item]
        assert d.value_text == r["value"], d.key
        assert d.doc_class == r["cls"], d.key
        assert d.source == r["source"], d.key
        if d.layout_class is None:  # never drawn: only MISSING rows
            assert r["cls"].startswith("MISSING"), d.key
            assert d.value_mm is None
            continue
        assert d.layout_class in gds.CLASSES
        if d.key == "bond_pad":
            assert r["cls"] == "QUALITATIVE" and d.layout_class == "PLACEHOLDER"
        else:
            assert d.layout_class == r["cls"].replace(" ", "-"), d.key


def test_table_numbers_parse_from_doc_values():
    t = gds.DIMENSIONS
    assert t["active_dosimeter"].value_text.startswith("2 x 2 mm")
    assert t["active_dosimeter"].value_mm == 2.0
    assert t["active_microdosimeter"].value_text.startswith("5 x 5 mm")
    assert t["active_microdosimeter"].value_mm == 5.0
    die = t["die_dosimeter"]
    assert die.value_text.startswith(f"{die.value_mm:g} +- {die.uncertainty_mm:g} mm")
    lo, hi = gds.PERIPHERY_RANGE_MM
    assert t["periphery"].value_text.startswith(f"{lo:g}-{hi:g} mm")


def test_termination_zone_width_within_documented_range():
    t = gds.DIMENSIONS
    band = (t["die_dosimeter"].value_mm - t["active_dosimeter"].value_mm) / 2
    lo, hi = gds.PERIPHERY_RANGE_MM
    assert lo <= band <= hi


# --------------------------------------------------------------------------- geometry


def test_layouts_pure_python():
    dos = gds.build_layout("dosimeter")
    mic = gds.build_layout("microdosimeter")
    for lay in (dos, mic):
        assert "REFERENCE-NOT-FAB" in lay.cell_name
        assert any(
            n.layer == "DISCLAIMER" and n.text == gds.DISCLAIMER for n in lay.notes
        )
        for s in lay.shapes:
            assert s.label.split(" | ")[0] in gds.CLASSES
            assert "src: " in s.label
    layers = [s.layer for s in dos.shapes]
    assert layers.count("TERMINATION_ZONE") == 1  # one zone, no discrete rings
    assert layers.count("PAD_PLACEHOLDER") == 1
    pad = next(s for s in dos.shapes if s.layer == "PAD_PLACEHOLDER")
    assert "SIZE UNKNOWN" in pad.label and pad.label.startswith("PLACEHOLDER")
    assert [s.layer for s in mic.shapes] == ["ACTIVE"]  # no die outline
    assert any("die size MISSING" in n.text for n in mic.notes)
    assert not any("DICING" in name for name in gds.LAYERS)
    with pytest.raises(ValueError):
        gds.build_layout("pixel")


def test_active_areas_4_and_25_mm2():
    act = {
        d: next(s for s in gds.build_layout(d).shapes if s.layer == "ACTIVE")
        for d in gds.DEVICES
    }
    assert _polygon_area_um2(act["dosimeter"].points_um) / 1e6 == pytest.approx(4.0)
    assert _polygon_area_um2(act["microdosimeter"].points_um) / 1e6 == pytest.approx(
        25.0
    )
    zone = next(
        s for s in gds.build_layout("dosimeter").shapes if s.layer == "TERMINATION_ZONE"
    )
    assert _polygon_area_um2(zone.points_um) / 1e6 == pytest.approx(2.7**2 - 4.0)


@pytest.mark.parametrize("device", gds.DEVICES)
def test_svg_preview(device):
    svg = gds.layout_svg(device)
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert "REFERENCE-NOT-FAB" in svg
    assert gds.DISCLAIMER in svg
    assert gds.VALIDATION_LABEL in svg
    lay = gds.build_layout(device)
    assert svg.count("<polygon") == len(lay.shapes)
    text = " ".join(re.findall(r">([^<]*)</text>", svg))
    _no_forbidden(text, f"{device} SVG text")
    if device == "dosimeter":
        assert 'fill="url(#hatch)"' in svg and "SIZE UNKNOWN" in svg
    else:
        assert "die size MISSING" in svg


def test_module_source_free_of_mesh_defaults():
    _no_forbidden(Path(gds.__file__).read_text(), "etna/export/gds.py")


# --------------------------------------------------------------------------- GDSII


@pytest.fixture(scope="module")
def readback(tmp_path_factory):
    gdstk = pytest.importorskip("gdstk")
    out = {}
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
    assert "REFERENCE-NOT-FAB" in top.name
    assert top.name == gds.build_layout(device).cell_name
    dl, dt = gds.LAYERS["DISCLAIMER"]
    assert any(
        lb.text == gds.DISCLAIMER and lb.layer == dl and lb.texttype == dt
        for lb in top.labels
    )
    assert any(gds.VALIDATION_LABEL in lb.text for lb in top.labels)
    assert lib.unit == pytest.approx(1e-6)


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
    by_layer = {(p.layer, p.datatype): p for p in top.polygons}
    for s in lay.shapes:
        p = by_layer[gds.LAYERS[s.layer]]
        assert p.area() == pytest.approx(_polygon_area_um2(s.points_um), rel=1e-9)
    act = by_layer[gds.LAYERS["ACTIVE"]]
    assert act.area() / 1e6 == pytest.approx(4.0 if device == "dosimeter" else 25.0)
    assert set(by_layer) <= set(gds.LAYERS.values())
    # No label text carries the model's mesh defaults; no polygon side equals one.
    for lb in top.labels:
        _no_forbidden(lb.text, f"{device} GDS label")
    for p in top.polygons:
        (x0, y0), (x1, y1) = p.bounding_box()
        for side in (x1 - x0, y1 - y0):
            assert min(abs(side - v) for v in (50.0, 5.0, 3.0, 1.0)) > 1e-6
