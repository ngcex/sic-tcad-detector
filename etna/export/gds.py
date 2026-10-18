"""Labelled SCHEMATIC reference layouts (GDSII + SVG preview): NOT fabrication layouts.

Two cells, one per device, each a top cell whose name contains ``REFERENCE-NOT-FAB``:

- dosimeter (2 x 2 mm active): active square, die outline (2.7 +- 0.15 mm, measured
  from a figure), the periphery drawn as ONE 'termination zone' polygon (hatched in
  the SVG; GDSII has no fill patterns) instead of discrete rings, and ONE corner pad
  placeholder marked SIZE UNKNOWN. No dicing layer, no ring widths.
- microdosimeter (5 x 5 mm active): the active square only, plus a note that the die
  size is MISSING. DESIGN-2's sensitive-volume array has no documented counterpart and
  is not drawn.

Every dimension comes from ``DIMENSIONS`` (generated from
``deliverables/DIMENSIONS-provenance.md``; ``tests/test_export_gds.py`` checks the two
cannot drift). Every polygon carries a text label (same layer, origin inside the
polygon) with its class (DOCUMENTED / MEASURED-FROM-FIGURE / PLACEHOLDER) and source.
The guard-ring numbers of the simulation model are mesh defaults, not device
dimensions, and appear nowhere here.

The SVG is pure Python; the GDSII needs the optional ``gdstk`` (``viz`` extra):
``gds_available()`` tells whether it can be written. Units: text in mm, GDSII
coordinates in um (unit 1e-6 m, precision 1e-9 m).
"""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from html import escape

from etna.export.labels import VALIDATION_LABEL

DISCLAIMER = (
    "schematic, not a fabrication layout; missing dimensions must come from the mask "
    "data (SiCILIA STM/IMM-CNR via G. Petringa)"
)
NO_SIM_NOTE = (
    "This layout contains no simulation output (geometry from the literature only)."
)

CLASSES = ("DOCUMENTED", "MEASURED-FROM-FIGURE", "PLACEHOLDER")
DEVICES = ("dosimeter", "microdosimeter")


# --------------------------------------------------------------------------- table


@dataclass(frozen=True)
class Dim:
    """One row of deliverables/DIMENSIONS-provenance.md used by the layout.

    ``doc_item`` / ``value_text`` / ``source`` are verbatim from the doc row;
    ``doc_class`` is the doc's class, ``layout_class`` the class the drawn polygon
    carries (the bond pad is QUALITATIVE in the doc and a PLACEHOLDER in the layout).
    ``value_mm`` / ``uncertainty_mm`` are None when the doc gives no number.
    """

    key: str
    doc_item: str
    value_text: str
    value_mm: "float | None"
    uncertainty_mm: "float | None"
    doc_class: str
    layout_class: "str | None"
    source: str


DIMENSIONS: dict[str, Dim] = {
    d.key: d
    for d in (
        Dim(
            "active_dosimeter",
            "active p+ area, dosimeter",
            "2 x 2 mm (4 mm2)",
            2.0,
            None,
            "DOCUMENTED",
            "DOCUMENTED",
            "JINST 20 C08019 p.3; Photons II.A p.3",
        ),
        Dim(
            "active_microdosimeter",
            "active p+ area, microdosimeter",
            "5 x 5 mm (25 mm2)",
            5.0,
            None,
            "DOCUMENTED",
            "DOCUMENTED",
            "Phys. Med. 140 (2025) 105675 s2.1; Altana, Sensors 23:6522",
        ),
        Dim(
            "die_dosimeter",
            "die size, 2 mm device",
            "2.7 +- 0.15 mm square",
            2.7,
            0.15,
            "MEASURED FROM FIGURE",
            "MEASURED-FROM-FIGURE",
            "Photons Fig.1b (81/79 px = 2 mm, ~25 um/px)",
        ),
        Dim(
            "periphery",
            "periphery active edge -> die edge",
            "0.2-0.4 mm per side (+-0.1)",
            None,  # a range: (0.2, 0.4) in PERIPHERY_RANGE_MM
            0.1,
            "MEASURED FROM FIGURE",
            "MEASURED-FROM-FIGURE",
            "Photons Fig.1b",
        ),
        Dim(
            "bond_pad",
            "bond pad",
            "1 wire bond at a die corner; size unknown",
            None,
            None,
            "QUALITATIVE",
            "PLACEHOLDER",
            "Photons Fig.1b",
        ),
        Dim(
            "dicing_street",
            "dicing street",
            "unknown",
            None,
            None,
            "MISSING",
            None,
            "-",
        ),
        Dim(
            "missing_mask_data",
            "ring widths/gaps, oxide thickness, pad size, 5x5 die size",
            "unknown",
            None,
            None,
            "MISSING: owner/fab (SiCILIA mask set: STM / IMM-CNR via G. Petringa)",
            None,
            "-",
        ),
    )
}
PERIPHERY_RANGE_MM = (0.2, 0.4)

# GDSII layers have no names: (layer, datatype) per name; the map is also written as a
# NOTES label in every cell. No dicing layer.
LAYERS: dict[str, tuple[int, int]] = {
    "ACTIVE": (1, 0),
    "DIE_OUTLINE": (2, 0),
    "TERMINATION_ZONE": (3, 0),
    "PAD_PLACEHOLDER": (4, 0),
    "NOTES": (10, 0),
    "DISCLAIMER": (63, 0),
}
LAYER_MAP_NOTE = "LAYER MAP: " + ", ".join(
    f"{n}={l}/{d}" for n, (l, d) in LAYERS.items()
)


# --------------------------------------------------------------------------- geometry


@dataclass(frozen=True)
class Shape:
    layer: str
    points_um: tuple[tuple[float, float], ...]
    label: str
    label_xy_um: tuple[float, float]
    hatched: bool = False


@dataclass(frozen=True)
class Note:
    layer: str
    text: str
    xy_um: tuple[float, float]
    size: float  # magnification (GDSII) / relative font size (SVG)


@dataclass(frozen=True)
class Layout:
    device: str
    cell_name: str
    shapes: tuple[Shape, ...]
    notes: tuple[Note, ...]
    extent_um: float  # half-width of the drawing area (for the SVG view box)


def _square(half: float) -> tuple[tuple[float, float], ...]:
    return ((-half, -half), (half, -half), (half, half), (-half, half))


def _frame(outer: float, inner: float) -> tuple[tuple[float, float], ...]:
    """One keyhole polygon: outer square minus inner square (a single boundary)."""
    return (
        (outer, outer),
        (-outer, outer),
        (-outer, -inner),
        (-inner, -inner),
        (-inner, inner),
        (inner, inner),
        (inner, -inner),
        (-inner, -inner),
        (-outer, -inner),
        (-outer, -outer),
        (outer, -outer),
    )


def _label(d: Dim, what: str, drawn: str = "") -> str:
    tail = f"; {drawn}" if drawn else ""
    return f"{d.layout_class} | {what}: {d.value_text}{tail} | src: {d.source}"


def _common_notes(extent: float, extra: tuple[str, ...]) -> tuple[Note, ...]:
    top = extent * 0.92
    lines = (
        ("DISCLAIMER", DISCLAIMER, (-extent * 0.95, top), 1.8),
        ("NOTES", LAYER_MAP_NOTE, (-extent * 0.95, -extent * 0.80), 1.0),
        (
            "NOTES",
            "ETNA export label: " + VALIDATION_LABEL,
            (-extent * 0.95, -extent * 0.86),
            1.0,
        ),
        ("NOTES", NO_SIM_NOTE, (-extent * 0.95, -extent * 0.92), 1.0),
    )
    out = [Note(*ln) for ln in lines]
    for i, t in enumerate(extra):
        out.append(
            Note("NOTES", t, (-extent * 0.95, top - extent * (0.10 + 0.06 * i)), 1.2)
        )
    return tuple(out)


def build_layout(device: str) -> Layout:
    """Pure-Python geometry + labels of one device (no gdstk needed)."""
    if device == "dosimeter":
        act = DIMENSIONS["active_dosimeter"]
        die = DIMENSIONS["die_dosimeter"]
        per = DIMENSIONS["periphery"]
        pad = DIMENSIONS["bond_pad"]
        a = act.value_mm * 1e3 / 2.0
        b = die.value_mm * 1e3 / 2.0
        band = b - a  # per-side periphery implied by the two values (0.35 mm)
        # Pad glyph: symbolic extent (a fraction of the corner band), NOT a dimension.
        g = 0.3 * band
        cx = cy = a + band / 2.0  # corner choice arbitrary (paper: "a die corner")
        shapes = (
            Shape(
                "ACTIVE",
                _square(a),
                _label(act, "active p+ area, dosimeter"),
                (-a * 0.9, a * 0.5),
            ),
            Shape(
                "DIE_OUTLINE",
                _square(b),
                _label(
                    die,
                    "die outline (reference only, not a dicing line)",
                ),
                (0.0, -(a + band * 0.3)),
            ),
            Shape(
                "TERMINATION_ZONE",
                _frame(b, a),
                _label(
                    per,
                    "termination zone",
                    "drawn = (die - active)/2; rings/JTE/field plate inside unknown",
                ),
                (0.0, -(a + band * 0.7)),
                hatched=True,
            ),
            Shape(
                "PAD_PLACEHOLDER",
                (
                    (cx - g, cy - g),
                    (cx + g, cy - g),
                    (cx + g, cy + g),
                    (cx - g, cy + g),
                ),
                _label(
                    pad,
                    "bond pad SIZE UNKNOWN",
                    "drawn extent symbolic, corner arbitrary",
                ),
                (cx, cy),
            ),
        )
        extent = b * 1.6
        extra = (
            "DOSIMETER 2 x 2 mm: periphery drawn as ONE termination zone, not as rings",
            "MISSING (not drawn): ring widths/gaps, pad size, dicing street, oxide thickness",
        )
        return Layout(
            device,
            "ETNA_DOSIMETER_2x2_REFERENCE-NOT-FAB",
            shapes,
            _common_notes(extent, extra),
            extent,
        )
    if device == "microdosimeter":
        act = DIMENSIONS["active_microdosimeter"]
        a = act.value_mm * 1e3 / 2.0
        shapes = (
            Shape(
                "ACTIVE",
                _square(a),
                _label(act, "active p+ area, microdosimeter"),
                (-a * 0.9, a * 0.5),
            ),
        )
        extent = a * 1.6
        extra = (
            "MICRODOSIMETER 5 x 5 mm: die size MISSING (only the active square is drawn)",
            "MISSING (not drawn): die size, ring widths/gaps, pad size, dicing street",
            "One 5 x 5 mm diode; DESIGN-2's sensitive-volume array has no documented "
            "counterpart and is not drawn",
        )
        return Layout(
            device,
            "ETNA_MICRODOSIMETER_5x5_REFERENCE-NOT-FAB",
            shapes,
            _common_notes(extent, extra),
            extent,
        )
    raise ValueError(f"device must be one of {DEVICES}, got {device!r}")


# --------------------------------------------------------------------------- GDSII


def gds_available() -> "tuple[bool, str]":
    try:
        import gdstk  # noqa: F401
    except Exception as e:  # pragma: no cover - depends on the environment
        return False, f"gdstk not installed ({e}); install the 'viz' extra"
    return True, ""


def layout_library(device: str):
    """gdstk.Library with one top cell (name contains REFERENCE-NOT-FAB)."""
    import gdstk

    lay = build_layout(device)
    lib = gdstk.Library("ETNA_REFERENCE-NOT-FAB", unit=1e-6, precision=1e-9)
    cell = lib.new_cell(lay.cell_name)
    for s in lay.shapes:
        layer, dt = LAYERS[s.layer]
        cell.add(gdstk.Polygon(s.points_um, layer=layer, datatype=dt))
        cell.add(gdstk.Label(s.label, s.label_xy_um, layer=layer, texttype=dt))
    for n in lay.notes:
        layer, dt = LAYERS[n.layer]
        cell.add(
            gdstk.Label(
                n.text,
                n.xy_um,
                magnification=n.size * (100.0 if n.layer == "DISCLAIMER" else 30.0),
                layer=layer,
                texttype=dt,
            )
        )
    return lib


def layout_gds(device: str) -> bytes:
    """GDSII bytes of one device's schematic (needs gdstk)."""
    lib = layout_library(device)
    fd, path = tempfile.mkstemp(suffix=".gds", prefix="etna_ref_")
    os.close(fd)
    try:
        lib.write_gds(path)
        with open(path, "rb") as f:
            return f.read()
    finally:
        os.unlink(path)


# --------------------------------------------------------------------------- SVG

_SVG_STYLE = {
    "ACTIVE": 'fill="#4c8bf5" fill-opacity="0.35" stroke="#1a56c4"',
    "DIE_OUTLINE": 'fill="none" stroke="#444" stroke-dasharray="40 20"',
    "TERMINATION_ZONE": 'fill="url(#hatch)" stroke="#b5651d"',
    "PAD_PLACEHOLDER": 'fill="#f5c04c" fill-opacity="0.6" stroke="#a07000" stroke-dasharray="10 10"',
}


_SVG_ANCHOR = {
    "DIE_OUTLINE": "middle",
    "TERMINATION_ZONE": "middle",
    "PAD_PLACEHOLDER": "end",
}


def layout_svg(device: str) -> str:
    """SVG preview (pure Python): same geometry and labels as the GDSII, y up."""
    lay = build_layout(device)
    e = lay.extent_um
    fs = e / 60.0  # base font size in drawing units

    def pt(x, y):
        return f"{x:.1f},{-y:.1f}"

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{-e:.1f} {-e:.1f} {2 * e:.1f} '
        f'{2 * e:.1f}" width="900" height="900" font-family="sans-serif">',
        f"<title>{escape(lay.cell_name)}</title>",
        "<defs>"
        '<pattern id="hatch" patternUnits="userSpaceOnUse" width="40" height="40" '
        'patternTransform="rotate(45)">'
        '<line x1="0" y1="0" x2="0" y2="40" stroke="#b5651d" stroke-width="8"/>'
        "</pattern></defs>",
        f'<rect x="{-e:.1f}" y="{-e:.1f}" width="{2 * e:.1f}" height="{2 * e:.1f}" fill="white"/>',
    ]
    for s in lay.shapes:
        pts = " ".join(pt(x, y) for x, y in s.points_um)
        out.append(
            f'<polygon data-layer="{s.layer}" points="{pts}" {_SVG_STYLE[s.layer]} '
            f'stroke-width="{fs / 4:.1f}" fill-rule="evenodd"/>'
        )
    for s in lay.shapes:
        x, y = s.label_xy_um
        anchor = _SVG_ANCHOR.get(s.layer, "start")
        out.append(
            f'<text data-layer="{s.layer}" x="{x:.1f}" y="{-y:.1f}" font-size="{fs * 0.7:.1f}" '
            f'text-anchor="{anchor}" paint-order="stroke" stroke="white" '
            f'stroke-width="{fs * 0.25:.1f}">{escape(s.label)}</text>'
        )
    for n in lay.notes:
        x, y = n.xy_um
        color = "#c00000" if n.layer == "DISCLAIMER" else "#222"
        weight = ' font-weight="bold"' if n.layer == "DISCLAIMER" else ""
        out.append(
            f'<text data-layer="{n.layer}" x="{x:.1f}" y="{-y:.1f}" '
            f'font-size="{fs * 0.8 * n.size:.1f}" fill="{color}"{weight}>{escape(n.text)}</text>'
        )
    out.append("</svg>")
    return "\n".join(out) + "\n"
