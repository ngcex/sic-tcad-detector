"""Labelled reference layouts (GDSII + SVG preview) from the layout parameter file.

Reference layout, not TCAD-validated, and not a fabrication layout while any dimension
is a PLACEHOLDER. Every dimension comes from ``data/foundry/layout_params.toml`` (or an
alternative file / mapping, see ``etna.export.layout_params``); nothing is hard-coded
here. One top cell per device, drawn centred on the active area:

- ACTIVE: the active p+ square (its boundary is the p+ edge);
- GUARD_RING: ``guard_ring_count`` square rings, ``guard_ring_gap`` from the p+ edge,
  ``guard_ring_width`` wide, ``guard_ring_spacing`` apart;
- PAD: one square bond pad (``pad_size``) in a die corner, between the outermost ring and
  the die edge (corner choice arbitrary: the source says "a die corner");
- DIE_OUTLINE: the die square; DICING_STREET: half the street around it.

Every polygon carries a text label (same layer, origin inside the polygon)
``CLASS | what: value | src: source``. The top cell name ends in ``REFERENCE-NOT-FAB``
while the device has any PLACEHOLDER dimension (``REFERENCE`` once none is left); the
placeholder names are written as a NOTES label and announced with a
``PlaceholderDimensionWarning``. The simulation model's guard-ring numbers are mesh
defaults, not device dimensions, and appear nowhere here.

The SVG is pure Python; the GDSII needs the optional ``gdstk`` (``viz`` extra):
``gds_available()`` tells whether it can be written. GDSII coordinates in um (unit
1e-6 m, precision 1e-9 m).
"""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from html import escape

from etna.export.labels import VALIDATION_LABEL
from etna.export.layout_params import (
    CLASSES,
    DEVICES,
    LayoutParams,
    as_layout_params,
    warn_placeholders,
)

REFERENCE_LABEL = "reference layout, not TCAD-validated"
DISCLAIMER = (
    f"{REFERENCE_LABEL}; not a fabrication layout while any dimension is a "
    "PLACEHOLDER: missing dimensions must come from the mask data (SiCILIA STM/IMM-CNR "
    "via G. Petringa)"
)
NO_SIM_NOTE = (
    "This layout contains no simulation output (geometry from the parameter file only)."
)

__all__ = [
    "CLASSES",
    "DEVICES",
    "DISCLAIMER",
    "REFERENCE_LABEL",
    "LAYERS",
    "build_layout",
    "gds_available",
    "layout_library",
    "layout_gds",
    "layout_svg",
]

# GDSII layers have no names: (layer, datatype) per name; the map is also written as a
# NOTES label in every cell.
LAYERS: dict[str, tuple[int, int]] = {
    "ACTIVE": (1, 0),
    "DIE_OUTLINE": (2, 0),
    "GUARD_RING": (3, 0),
    "PAD": (4, 0),
    "DICING_STREET": (5, 0),
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
    dims: tuple[str, ...] = ()  # parameter names the shape is drawn from


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
    placeholders: tuple[str, ...]  # 'device.name' still PLACEHOLDER
    params_path: str


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


def _worst(*dims) -> str:
    """Class of a shape drawn from several dimensions: the weakest one."""
    return max((d.cls for d in dims), key=CLASSES.index)


def _label(what: str, *dims) -> str:
    vals = "; ".join(f"{d.name} {d.text()}" for d in dims)
    srcs = "; ".join(dict.fromkeys(d.source for d in dims))
    return f"{_worst(*dims)} | {what}: {vals} | src: {srcs}"


def _notes(extent: float, lines: "tuple[str, ...]") -> tuple[Note, ...]:
    top = extent * 0.92
    base = (
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
    out = [Note(*ln) for ln in base]
    for i, t in enumerate(lines):
        out.append(
            Note("NOTES", t, (-extent * 0.95, top - extent * (0.10 + 0.06 * i)), 1.2)
        )
    return tuple(out)


def _check_device(device: str) -> None:
    if device not in DEVICES:
        raise ValueError(f"device must be one of {DEVICES}, got {device!r}")


def _build(device: str, params: LayoutParams) -> Layout:
    _check_device(device)
    d = params.devices[device]
    act, die, street = d["active_side"], d["die_side"], d["dicing_street"]
    n_ring = d["guard_ring_count"]
    gap, width, spacing = (
        d["guard_ring_gap"],
        d["guard_ring_width"],
        d["guard_ring_spacing"],
    )
    pad = d["pad_size"]
    a, b = act.um / 2.0, die.um / 2.0
    s = street.um / 2.0
    shapes = [
        Shape(
            "ACTIVE",
            _square(a),
            _label("active p+ area (boundary = p+ edge)", act),
            (-a * 0.9, a * 0.5),
            ("active_side",),
        )
    ]
    n = int(n_ring.value)
    r_out = a
    for k in range(n):
        r_in = a + gap.um + k * (width.um + spacing.um)
        r_out = r_in + width.um
        used = (n_ring, gap, width) + ((spacing,) if k else ())
        shapes.append(
            Shape(
                "GUARD_RING",
                _frame(r_out, r_in),
                _label(f"guard ring {k + 1} of {n}", *used),
                (0.0, -(r_in + r_out) / 2.0),
                tuple(x.name for x in used),
            )
        )
    c = (r_out + b) / 2.0  # pad centre: corner, between outer ring and die edge
    g = pad.um / 2.0
    shapes.append(
        Shape(
            "PAD",
            ((c - g, c - g), (c + g, c - g), (c + g, c + g), (c - g, c + g)),
            _label("bond pad (corner position arbitrary)", pad),
            (c, c),
            ("pad_size",),
        )
    )
    shapes.append(
        Shape(
            "DIE_OUTLINE",
            _square(b),
            _label("die outline (reference only, not a dicing line)", die),
            (0.0, -(b - (b - r_out) * 0.25)),
            ("die_side",),
        )
    )
    shapes.append(
        Shape(
            "DICING_STREET",
            _frame(b + s, b),
            _label("half dicing street", street),
            (0.0, -(b + s / 2.0)),
            ("dicing_street",),
        )
    )
    ph = tuple(params.placeholders(device))
    side = f"{act.value:g}x{act.value:g}{act.unit}".replace(".", "p")
    cell = f"ETNA_{device.upper()}_{side}_" + (
        "REFERENCE-NOT-FAB" if ph else "REFERENCE"
    )
    extent = (b + s) * 1.6
    lines = (
        f"{device.upper()}: {params.placeholder_note(device)}",
        f"dimensions from {params.display_path}",
    )
    return Layout(
        device,
        cell,
        tuple(shapes),
        _notes(extent, lines),
        extent,
        ph,
        params.display_path,
    )


def build_layout(device: str, params=None) -> Layout:
    """Pure-Python geometry + labels of one device (no gdstk needed).

    ``params``: None (default file / $ETNA_LAYOUT_PARAMS), a path, a mapping or a
    ``LayoutParams``. Warns ``PlaceholderDimensionWarning`` if placeholders remain.
    """
    _check_device(device)
    p = as_layout_params(params)
    warn_placeholders(p, device, f"layout {device}")
    return _build(device, p)


# --------------------------------------------------------------------------- GDSII


def gds_available() -> "tuple[bool, str]":
    try:
        import gdstk  # noqa: F401
    except Exception as e:  # pragma: no cover - depends on the environment
        return False, f"gdstk not installed ({e}); install the 'viz' extra"
    return True, ""


def layout_library(device: str, params=None):
    """gdstk.Library with one top cell (name ends in REFERENCE-NOT-FAB while any
    PLACEHOLDER remains)."""
    import gdstk

    _check_device(device)
    p = as_layout_params(params)
    warn_placeholders(p, device, f"GDSII layout {device}")
    lay = _build(device, p)
    lib = gdstk.Library(
        "ETNA_" + lay.cell_name.split("_")[-1], unit=1e-6, precision=1e-9
    )
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


def layout_gds(device: str, params=None) -> bytes:
    """GDSII bytes of one device's reference layout (needs gdstk)."""
    lib = layout_library(device, params)  # warns on placeholders
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
    "GUARD_RING": 'fill="#b5651d" fill-opacity="0.55" stroke="#7a3f0c"',
    "PAD": 'fill="#f5c04c" fill-opacity="0.6" stroke="#a07000"',
    "DICING_STREET": 'fill="#999" fill-opacity="0.25" stroke="#777"',
}
_SVG_ANCHOR = {
    "DIE_OUTLINE": "middle",
    "GUARD_RING": "middle",
    "DICING_STREET": "middle",
    "PAD": "end",
}


def layout_svg(device: str, params=None) -> str:
    """SVG preview (pure Python): same geometry and labels as the GDSII, y up."""
    lay = build_layout(device, params)
    e = lay.extent_um
    fs = e / 60.0  # base font size in drawing units

    def pt(x, y):
        return f"{x:.1f},{-y:.1f}"

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{-e:.1f} {-e:.1f} {2 * e:.1f} '
        f'{2 * e:.1f}" width="900" height="900" font-family="sans-serif">',
        f"<title>{escape(lay.cell_name)}</title>",
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
