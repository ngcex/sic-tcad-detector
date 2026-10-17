"""Dimensioned cross-section of the simulated stack (SVG hand-written, PNG via matplotlib).

The vertical (depth) axis is to scale. The lateral extent is NOT: only the active area
``area_cm2`` is known, not the die shape or the guard ring / pad layout, so the width is
schematic and labelled as such (no dimension is invented).
"""

from __future__ import annotations

import io
import textwrap
from xml.sax.saxutils import escape

from etna.export.labels import (
    NPLUS_THICKNESS_NOTE,
    VALIDATION_LABEL,
    epi_donor_profile,
    has_nplus,
)

_COLORS = {"p+": "#d6604d", "epi": "#92c5de", "n+": "#2166ac"}


def _fmt_um(v: float) -> str:
    return f"{v:g} µm"


def stack_layers(config) -> "list[dict]":
    """Layers top (anode) to bottom (cathode): name, depth range (um), doping text, note."""
    t_p = config.substrate_thickness_um
    t_e = config.epi_thickness_um
    layers = [
        dict(
            key="p+",
            name="p+ anode",
            top_um=0.0,
            bottom_um=t_p,
            thickness_um=t_p,
            doping=f"N_A = {config.N_A:.3g} cm^-3 (acceptors, total)",
            note="",
        )
    ]
    if config.doping_profile == "uniform":
        epi_dop = f"N_D = {float(epi_donor_profile(config, 0.0)):.3g} cm^-3 (uniform)"
    else:
        epi_dop = (
            f"graded N_D: {config.N_D_junction:.3g} at junction -> "
            f"{float(epi_donor_profile(config, t_e)):.3g} cm^-3 at epi end "
            f"(bulk {config.N_D_bulk:.3g}, L = {config.L_transition_um:g} µm)"
        )
    layers.append(
        dict(
            key="epi",
            name="n- epi",
            top_um=t_p,
            bottom_um=t_p + t_e,
            thickness_um=t_e,
            doping=epi_dop,
            note="",
        )
    )
    if has_nplus(config):
        t_n = config.nplus_thickness_um
        layers.append(
            dict(
                key="n+",
                name="n+ substrate",
                top_um=t_p + t_e,
                bottom_um=t_p + t_e + t_n,
                thickness_um=t_n,
                doping=f"N_D = {config.N_D_nplus:.3g} cm^-3",
                note=NPLUS_THICKNESS_NOTE.format(t=t_n),
            )
        )
    return layers


def _lateral_note(config) -> str:
    return (
        f"lateral: schematic, not to scale; active area {config.area_cm2:g} cm^2 "
        f"({config.area_cm2 * 100:g} mm^2), shape not specified"
    )


def cross_section_svg(config, width_px: int = 900) -> str:
    """SVG string of the stack, depth to scale, with thickness quotes and doping."""
    layers = stack_layers(config)
    total = layers[-1]["bottom_um"]
    top, left, block_w = 70.0, 110.0, 300.0
    body_h = 460.0
    s = body_h / total  # px per um, vertical
    height = top + body_h + 130
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width_px}" height="{height:.0f}" '
        f'viewBox="0 0 {width_px} {height:.0f}" font-family="Helvetica, Arial, sans-serif">',
        '<defs><marker id="arr" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="6" '
        'markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" '
        'fill="#222"/></marker></defs>',
        f'<rect x="0" y="0" width="{width_px}" height="{height:.0f}" fill="white"/>',
        f'<text x="{left}" y="28" font-size="18" font-weight="bold">'
        "Device cross-section (depth to scale)</text>",
        f'<text x="{left}" y="50" font-size="12">{escape(_lateral_note(config))}</text>',
    ]
    # anode / cathode contacts (drawn as lines; metal thickness not modelled)
    out.append(
        f'<line x1="{left}" y1="{top}" x2="{left + block_w}" y2="{top}" stroke="#555" '
        'stroke-width="4"/>'
    )
    for L in layers:
        y0 = top + L["top_um"] * s
        h = L["thickness_um"] * s
        out.append(
            f'<rect x="{left}" y="{y0:.2f}" width="{block_w}" height="{h:.2f}" '
            f'fill="{_COLORS[L["key"]]}" stroke="#222" stroke-width="1"/>'
        )
        # thickness quote (dimension line) on the right
        xq = left + block_w + 30
        out.append(
            f'<line x1="{xq}" y1="{y0:.2f}" x2="{xq}" y2="{y0 + h:.2f}" stroke="#222" '
            'stroke-width="1" marker-start="url(#arr)" marker-end="url(#arr)"/>'
        )
        for yy in (y0, y0 + h):
            out.append(
                f'<line x1="{left + block_w}" y1="{yy:.2f}" x2="{xq + 8}" y2="{yy:.2f}" '
                'stroke="#888" stroke-width="0.6"/>'
            )
        ym = y0 + h / 2 + 4
        out.append(
            f'<text x="{xq + 14}" y="{ym:.2f}" font-size="13" font-weight="bold">'
            f'{escape(L["name"])}: {escape(_fmt_um(L["thickness_um"]))}</text>'
        )
        for k, part in enumerate(textwrap.wrap(L["doping"], 62)):
            out.append(
                f'<text x="{xq + 14}" y="{ym + 15 + 14 * k:.2f}" font-size="11">'
                f"{escape(part)}</text>"
            )
    bottom = top + total * s
    out.append(
        f'<line x1="{left}" y1="{bottom:.2f}" x2="{left + block_w}" y2="{bottom:.2f}" '
        'stroke="#555" stroke-width="4"/>'
    )
    out.append(
        f'<text x="{left + 6}" y="{top - 6}" font-size="11">anode contact</text>'
    )
    out.append(
        f'<text x="{left + 6}" y="{bottom + 16:.2f}" font-size="11">cathode contact</text>'
    )
    # depth axis on the left
    step = 1.0 if total <= 15 else (2.0 if total <= 30 else 5.0)
    out.append(
        f'<line x1="{left - 20}" y1="{top}" x2="{left - 20}" y2="{bottom:.2f}" '
        'stroke="#222" stroke-width="1"/>'
    )
    d = 0.0
    while d <= total + 1e-9:
        yy = top + d * s
        out.append(
            f'<line x1="{left - 25}" y1="{yy:.2f}" x2="{left - 20}" y2="{yy:.2f}" '
            'stroke="#222"/>'
        )
        out.append(
            f'<text x="{left - 30}" y="{yy + 4:.2f}" font-size="10" text-anchor="end">'
            f"{d:g}</text>"
        )
        d += step
    out.append(
        f'<text x="20" y="{top + body_h / 2:.0f}" font-size="12" '
        f'transform="rotate(-90 20 {top + body_h / 2:.0f})" text-anchor="middle">'
        "depth (µm)</text>"
    )
    yb = bottom + 45
    notes = [n["note"] for n in layers if n["note"]]
    for n in notes:
        for part in textwrap.wrap("n+: " + n, 110):
            out.append(
                f'<text x="{left}" y="{yb:.0f}" font-size="11">{escape(part)}</text>'
            )
            yb += 16
    out.append(
        f'<text x="{left}" y="{yb:.0f}" font-size="11">total modelled stack '
        f"{escape(_fmt_um(total))}; metal contacts drawn as lines (not modelled)</text>"
    )
    yb += 22
    out.append(
        f'<text x="{left}" y="{yb:.0f}" font-size="12" font-style="italic" '
        f'fill="#b2182b">{escape(VALIDATION_LABEL)}</text>'
    )
    out.append("</svg>")
    return "\n".join(out)


def _draw_cross_section(fig, config) -> None:
    """Draw the stack on a matplotlib Figure (shared by PNG and the PDF spec sheet)."""
    from matplotlib.patches import Rectangle

    layers = stack_layers(config)
    total = layers[-1]["bottom_um"]
    ax = fig.add_axes([0.12, 0.2, 0.36, 0.68])
    for L in layers:
        ax.add_patch(
            Rectangle(
                (0, L["top_um"]),
                1.0,
                L["thickness_um"],
                facecolor=_COLORS[L["key"]],
                edgecolor="#222",
                lw=0.8,
            )
        )
        ax.annotate(
            "",
            xy=(1.08, L["top_um"]),
            xytext=(1.08, L["bottom_um"]),
            arrowprops=dict(arrowstyle="<->", lw=0.8),
            annotation_clip=False,
        )
        ym = (L["top_um"] + L["bottom_um"]) / 2
        ax.text(
            1.14,
            ym,
            f"{L['name']}: {L['thickness_um']:g} µm\n"
            + "\n".join(textwrap.wrap(L["doping"], 55)),
            va="center",
            fontsize=7,
            clip_on=False,
        )
    ax.plot([0, 1], [0, 0], color="#555", lw=3)
    ax.plot([0, 1], [total, total], color="#555", lw=3)
    ax.set_xlim(0, 1)
    ax.set_ylim(total, 0)
    ax.set_xticks([])
    ax.set_ylabel("depth (µm)  [to scale]")
    ax.set_xlabel(_lateral_note(config), fontsize=7)
    fig.suptitle("Device cross-section (depth to scale)", fontsize=11)
    lines = [f"n+: {L['note']}" for L in layers if L["note"]]
    lines = [p for ln in lines for p in textwrap.wrap(ln, 140)]
    lines.append(VALIDATION_LABEL)
    fig.text(
        0.02, 0.03, "\n".join(lines), fontsize=6.5, color="#b2182b", style="italic"
    )


def cross_section_png(config, dpi: int = 150) -> bytes:
    """PNG bytes of the same drawing (matplotlib Agg canvas; no pyplot global state)."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(figsize=(9, 6))
    FigureCanvasAgg(fig)
    _draw_cross_section(fig, config)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    return buf.getvalue()
