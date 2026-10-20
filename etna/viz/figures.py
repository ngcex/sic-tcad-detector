"""Publication-quality figures of etna results: ``figure_for(result, kind)`` + export.

One entry point turns a finished result into a matplotlib ``Figure`` in the house style
(``etna.viz.style``), with physical units on the axes, the validity information that
comes with the result drawn on it, and a provenance footer (etna version, solver,
config hash, date, plus ``etna.export.VALIDATION_LABEL``); ``export_figure`` writes PNG /
SVG / PDF bytes at a set DPI and size. Nothing is simulated here and nothing is invented:
every drawn number comes from the result (or a committed data file) passed in.

Validity bands, not error bands: no etna solver returns an uncertainty, and the
engine-vs-devsim parity (<= 3e-6 relative, usually 1e-9..1e-12) is far too small to see,
so simulated curves have no error band. What IS drawn, when the result carries it:

- |V| beyond ``cross_checked_to_V`` hatched grey ("not cross-checked vs devsim");
- the full-depletion (reach-through) voltage as a dotted line;
- results outside the engine/devsim cross-check box, truncated sweeps, dark-current
  points past reach-through, unsolved (fully compensated) damage points with their reason;
- per-observable honesty notes (dark current = calibrated budget, FLASH = sensitivity
  bound, kappa placeholder, tilted alpha does not reproduce the data, ...);
- error bars on digitised measured points (their read-off uncertainty).

Accepted inputs: a ``SimResult`` (any ``run_*`` facade, ``etna.viz.maps.carrier_state_2d``,
``etna.viz.validation``), a ``DeviceConfig`` (kind ``"stack"``) or a CCE(LET) table dict
from ``load_let_table`` (kind ``"cce_let"``). ``kinds_for(obj)`` lists the kinds that apply.
"""

from __future__ import annotations

import io
import textwrap
from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib.figure import Figure
from matplotlib.layout_engine import ConstrainedLayoutEngine

from etna.viz import style as S
from etna.viz.provenance import bytes_hash, footer_text, provenance

KINDS = {
    "cv": "C-V: capacitance vs reverse bias",
    "mott_schottky": "Mott-Schottky: 1/C^2 vs reverse bias",
    "depletion": "Depletion width vs reverse bias (edge and eps*A/C)",
    "field": "E(x) and potential vs depth, with the layer stack",
    "cce": "CCE vs reverse bias",
    "cce_let": "CCE vs LET (committed single-particle table)",
    "dark_current": "Dark current vs reverse bias (or vs T)",
    "flash": "FLASH steady-state CCE vs dose rate",
    "transient": "FLASH transient CCE vs dose rate",
    "temperature": "CCE vs temperature",
    "damage": "CCE vs proton fluence",
    "microdosimetry": "Microdosimetric spectrum y d(y)",
    "map_potential": "2D potential map with contours",
    "map_field": "2D |E| map with contours",
    "map_doping": "2D net doping map",
    "map_electrons": "2D electron density map (engine2d only)",
    "map_holes": "2D hole density map (engine2d only)",
    "stack": "Device stack cross-section (depth to scale)",
    "tilted_alpha_shape": "Validation: tilted-alpha CCE shape vs Micro Fig. 4b",
    "w_creep": "Validation: W(V) creep vs Photons Fig. 4b labels",
}
MAP_KINDS = ("map_potential", "map_field", "map_doping", "map_electrons", "map_holes")

NOTES = {
    "dark_current": (
        "Dark current: calibration / budget estimate (N_t fitted to < 18 pA at -60 V, "
        "4 mm2, 300 K, used as equality); bias and T dependence heuristic."
    ),
    "flash": (
        "FLASH: numerical sensitivity bound, not a validated plasma-recombination "
        "prediction (high injection not modelled; SRH only, Auger inactive)."
    ),
    "transient": (
        "FLASH transient: low-injection DD pulse response; high-injection plasma "
        "physics not modelled."
    ),
    "damage": (
        "kappa (NIEL hardness) is a data-blocked placeholder: absolute fluence "
        "scale unvalidated."
    ),
    "temperature": "CCE vs T: model trend, not validated against a T measurement.",
    "cce_let": (
        "Committed table (built with the device model of its build date, not the "
        "current default); differs ~6e-6 from fresh builds, regeneration pending."
    ),
    "microdosimetry": "Spectrum from the supplied MC events; depends only on them and the SV.",
}


# --------------------------------------------------------------------------- public


def load_let_table(path) -> dict:
    """A committed CCE(LET) table (``data/cce_let_table_*.json``) as a figure input."""
    from etna.core.single_particle import load_cce_let_table

    p = Path(path)
    _, data = load_cce_let_table(p)
    data = dict(data)
    data["_source"] = f"{p.name} sha {bytes_hash(p.read_bytes())}"
    return data


def _is_config(obj) -> bool:
    return type(obj).__name__ == "DeviceConfig"


def _is_let_table(obj) -> bool:
    return isinstance(obj, dict) and "LET_keV_um" in obj and "CCE" in obj


def kinds_for(obj) -> "list[str]":
    """Figure kinds that apply to ``obj`` (first = default)."""
    if _is_config(obj):
        return ["stack"]
    if _is_let_table(obj):
        return ["cce_let"]
    st = getattr(obj, "sim_type", None)
    md = getattr(obj, "metadata", {}) or {}
    if st == "cv":
        return ["cv", "depletion", "mott_schottky", "stack"]
    if st == "field":
        if obj.mesh is not None and obj.mesh.y_coords is not None:
            ks = ["map_field", "map_potential", "map_doping"]
            if "Electrons" in obj.mesh.node_values:
                ks += ["map_electrons", "map_holes"]
            return ks + ["stack"]
        return ["field", "stack"]
    if st == "cce":
        return ["cce", "stack"]
    if st in ("dark_current", "dark_current_vs_T"):
        return ["dark_current", "stack"]
    if st in ("flash", "transient", "temperature", "damage", "microdosimetry"):
        return [st, "stack"]
    if st == "validation_tilted_alpha":
        return ["tilted_alpha_shape"]
    if st == "validation_w_creep":
        return ["w_creep"]
    raise ValueError(
        f"no figure kind for {type(obj).__name__} sim_type={st!r} ({md.get('solver')})"
    )


def figure_for(
    obj, kind: "str | None" = None, size=S.DEFAULT_SIZE, footer: bool = True, **opts
) -> Figure:
    """Publication figure of ``obj``; ``kind`` defaults to ``kinds_for(obj)[0]``.

    ``size``: "single" (89 x 70 mm), "double", "slide" or (w_mm, h_mm). ``footer``: draw
    the provenance / validation footer (on by default; turn off only for composites that
    carry it elsewhere). ``opts``: per-kind options (``title``, ``bias_label``).
    """
    allowed = kinds_for(obj)
    kind = kind or allowed[0]
    if kind not in KINDS:
        raise ValueError(f"unknown kind {kind!r}; known: {sorted(KINDS)}")
    if kind not in allowed:
        raise ValueError(
            f"kind {kind!r} does not apply to this input; use one of {allowed}"
        )
    with mpl.rc_context(S.RC):
        fig = Figure(figsize=S.size_inches(size), layout=ConstrainedLayoutEngine())
        notes: list[str] = []
        if kind == "stack":
            cfg = obj if _is_config(obj) else obj.config
            _stack(fig, cfg, notes)
            prov = provenance(config=cfg, solver="geometry (no solve)")
        elif kind == "cce_let":
            _cce_let(fig, obj, notes)
            g = obj.get("geometry", {})
            prov = provenance(
                solver="devsim 2D transient (table)", source=obj.get("_source", "table")
            )
            notes.insert(
                0,
                f"Table geometry: half-width {g.get('half_width_um')} um, epi "
                f"{g.get('epi_um')} um, reverse bias {abs(float(obj.get('bias_V', float('nan')))):g} V.",
            )
        else:
            _BUILDERS[kind](fig, obj, notes, **opts)
            prov = provenance(obj, source=obj.metadata.get("data_source"))
        if opts.get("title"):
            fig.suptitle(opts["title"])
        fig.etna_provenance = prov  # read by the app / tests
        fig.etna_kind = kind
        if footer:
            _footer(fig, prov, notes)
    return fig


def export_figure(
    fig: Figure, fmt: str = "png", dpi: int = S.DEFAULT_DPI, size=None
) -> bytes:
    """``fig`` as PNG / SVG / PDF bytes at ``dpi`` (raster) and ``size`` (None = as built).

    The canvas size is exact (no tight bbox): a PNG is int(w_in*dpi) x int(h_in*dpi)
    pixels, SVG/PDF are w_in x h_in at 72 pt/in.
    """
    fmt = fmt.lower()
    if fmt not in S.FORMATS:
        raise ValueError(f"fmt must be one of {S.FORMATS}")
    if not (30 <= int(dpi) <= 2400):
        raise ValueError("dpi must lie in 30..2400")
    with mpl.rc_context(S.RC):
        if size is not None:
            fig.set_size_inches(*S.size_inches(size))
            if getattr(fig, "_etna_footer_h_pt", None):
                _set_rect(fig)
        buf = io.BytesIO()
        meta = None
        prov = getattr(fig, "etna_provenance", None)
        if prov and fmt in ("pdf", "svg"):
            desc = footer_text(prov)
            meta = {
                "Title": f"ETNA {getattr(fig, 'etna_kind', '')}",
                "Creator": f"ETNA {prov['etna_version']}",
            }
            meta["Subject" if fmt == "pdf" else "Description"] = desc
        fig.savefig(buf, format=fmt, dpi=int(dpi), metadata=meta)
    return buf.getvalue()


def save_figure(fig: Figure, path, dpi: int = S.DEFAULT_DPI, size=None) -> Path:
    """Write ``fig`` to ``path`` (format from the suffix: .png / .svg / .pdf)."""
    p = Path(path)
    p.write_bytes(export_figure(fig, p.suffix.lstrip(".") or "png", dpi=dpi, size=size))
    return p


# --------------------------------------------------------------------------- helpers


def _footer(fig, prov, notes):
    from etna.export.labels import VALIDATION_LABEL

    w_in = fig.get_figwidth()
    chars = max(40, int(w_in * 72 / (0.62 * S.FOOTER_FONT_PT)))
    lines = []
    for n in notes:
        lines += textwrap.wrap(n, chars)
    lines += textwrap.wrap(footer_text(prov), chars)
    lines += textwrap.wrap(VALIDATION_LABEL, chars)
    fig.text(
        0.008,
        0.006,
        "\n".join(lines),
        fontsize=S.FOOTER_FONT_PT,
        color="#444444",
        ha="left",
        va="bottom",
        linespacing=1.15,
    )
    fig._etna_footer_h_pt = len(lines) * S.FOOTER_FONT_PT * 1.2 + 4.0
    _set_rect(fig)


def _set_rect(fig):
    frac = min(0.6, fig._etna_footer_h_pt / 72.0 / fig.get_figheight())
    fig.get_layout_engine().set(rect=(0.0, frac, 1.0, 1.0 - frac))


def _arr(v):
    return np.asarray(v, dtype=float)


def _common_notes(res, notes):
    md = res.metadata
    if md.get("validated_domain") is False:
        notes.append(
            "OUTSIDE the engine/devsim cross-check box"
            + (f" ({md['domain_note']})" if md.get("domain_note") else "")
            + "; not validated against measurement."
        )
    if md.get("fallback_reason"):
        notes.append(f"Solver fallback: {md['fallback_reason']}.")
    if md.get("truncated") and len(res.x):
        notes.append(
            f"Sweep truncated at {np.min(_arr(res.x)):.1f} V (requested "
            f"{md.get('requested_v_stop', float('nan')):.1f} V): devsim did not converge further."
        )


def _validity_bands(ax, res, xs_rev):
    """Reach-through line and the not-cross-checked band on a -V axis."""
    md = res.metadata
    v_rt = md.get("reach_through_V")
    if v_rt is not None and np.isfinite(v_rt):
        ax.axvline(
            v_rt, ls=":", color="#555555", lw=0.9, label=f"full depletion {v_rt:.2f} V"
        )
    v_cc = md.get("cross_checked_to_V")
    if v_cc is not None and xs_rev.size and np.max(xs_rev) > v_cc + 1e-9:
        ax.axvspan(
            v_cc,
            float(np.max(xs_rev)) * 1.02,
            facecolor="none",
            edgecolor=S.INVALID_COLOR,
            hatch="////",
            lw=0,
            label=f"|V| > {v_cc:.0f} V: not cross-checked vs devsim",
        )


def _rev(res):
    x = _arr(res.x)
    keep = x <= 0
    return keep, -x[keep]


def _legend(ax, **kw):
    h, _ = ax.get_legend_handles_labels()
    if h:
        ax.legend(**kw)


# --------------------------------------------------------------------------- builders


def _cv(fig, res, notes, **_):
    ax = fig.add_subplot()
    keep, v = _rev(res)
    f, lab = S.UNITS["capacitance"]
    ax.plot(v, _arr(res.y)[keep] * f, "o-", ms=2.0, label="small-signal C = dQ/dV")
    ax.set_xlabel(S.UNITS["reverse_bias"][1])
    ax.set_ylabel(lab)
    ax.set_title("C-V")
    area = res.metadata.get("area_cm2", getattr(res.config, "area_cm2", None))
    if area:
        notes.insert(0, f"Area {area * 100:g} mm2.")
    _validity_bands(ax, res, v)
    _legend(ax)
    _common_notes(res, notes)


def _mott(fig, res, notes, **_):
    ax = fig.add_subplot()
    keep, v = _rev(res)
    ax.plot(v, _arr(res.metadata["one_over_C_squared"])[keep], "o-", ms=2.0)
    ax.set_xlabel(S.UNITS["reverse_bias"][1])
    ax.set_ylabel(S.UNITS["inv_c2"][1])
    ax.set_title("Mott-Schottky")
    _validity_bands(ax, res, v)
    _legend(ax)
    _common_notes(res, notes)


def _epi_marks(ax, cfg):
    from etna.api.device import has_buffer

    ax.axhline(
        cfg.epi_thickness_um, color="#777777", lw=0.7, ls="--", label="epi thickness"
    )
    if has_buffer(cfg):
        ax.axhline(
            cfg.epi_thickness_um - cfg.buffer_thickness_um,
            color="#777777",
            lw=0.7,
            ls="-.",
            label="n-/n+ transition start",
        )


def _depletion(fig, res, notes, **_):
    ax = fig.add_subplot()
    keep, v = _rev(res)
    f, lab = S.UNITS["width_um"]
    ax.plot(
        v,
        _arr(res.metadata["depletion_widths"])[keep] * f,
        "-",
        label="depletion edge (n/N_D = 0.5)",
    )
    if "W_from_C" in res.metadata:
        ax.plot(
            v,
            _arr(res.metadata["W_from_C"])[keep] * f,
            "--",
            label=r"$\varepsilon A/C$ (small-signal)",
        )
    _epi_marks(ax, res.config)
    ax.set_xlabel(S.UNITS["reverse_bias"][1])
    ax.set_ylabel(lab)
    ax.set_title("Depletion width")
    _validity_bands(ax, res, v)
    _legend(ax, fontsize=S.BASE_FONT_PT - 2)
    _common_notes(res, notes)


def _layer_bands(ax, cfg, label=True):
    from etna.export.cross_section import stack_layers

    for L in stack_layers(cfg):
        ax.axvspan(
            L["top_um"],
            L["bottom_um"],
            color=S.LAYER_COLORS.get(L["key"], "#cccccc"),
            alpha=0.18,
            lw=0,
        )
        if label:
            ax.text(
                0.5 * (L["top_um"] + L["bottom_um"]),
                1.0,
                L["key"],
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=S.BASE_FONT_PT - 2.5,
                color="#333333",
            )


def _field(fig, res, notes, **_):
    from etna.api.results import peak_field_1d

    ax1, ax2 = fig.subplots(2, 1, sharex=True)
    x = _arr(res.x)
    E = _arr(res.y)
    f, lab = S.UNITS["field"]
    _layer_bands(ax1, res.config)
    _layer_bands(ax2, res.config, label=False)
    ax1.plot(x, E * f, color=S.PALETTE[0])
    pk = peak_field_1d(x, E, res.config)
    if pk is not None:
        ax1.plot(
            [pk[1]],
            [pk[0] * f],
            "v",
            color=S.PALETTE[1],
            ms=3.5,
            label=f"peak {pk[0] * f:.3g} kV/cm",
        )
        if pk[2]:
            ax1.set_ylim(0, 1.2 * pk[0] * f)
            notes.append(
                "E axis clipped: the nm-scale spikes at the abrupt doping steps "
                "(+-0.1 um, mesh dependent) are excluded from the peak and not shown."
            )
    ax1.set_ylabel(lab)
    bias = res.metadata.get("bias_V")
    ax1.set_title(f"Field at {bias:g} V" if bias is not None else "Field")
    _legend(ax1, loc="center right")
    ax2.plot(x, _arr(res.metadata["potential"]), color=S.PALETTE[2])
    ax2.set_ylabel(S.UNITS["potential"][1])
    ax2.set_xlabel(S.UNITS["depth"][1])
    ax2.set_xlim(x.min() if x.size else 0, x.max() if x.size else 1)
    _common_notes(res, notes)


def _cce(fig, res, notes, **_):
    ax = fig.add_subplot()
    keep, v = _rev(res)
    md = res.metadata
    ax.plot(v, _arr(res.y)[keep], "o-", ms=2.0, label="CCE (pairs in the epi)")
    if "cce_vs_E0" in md:
        ax.plot(
            v,
            _arr(md["cce_vs_E0"])[keep],
            "s--",
            ms=1.8,
            label=r"collected / ($E_0/E_{pair}$)",
        )
    if "theta_deg" in md:
        ax.set_title(f"Tilted alpha ({md['theta_deg']:g} deg), CCE vs bias")
        notes.append(
            f"Alpha generation shape '{md.get('generation_shape', '')}' is qualitative "
            "(not ASTAR/SRIM). Steady-state 1D low-injection model."
        )
    else:
        ax.set_title("CCE vs reverse bias")
    ax.axhline(1.0, color="#999999", lw=0.6, ls="--")
    ax.set_xlabel(S.UNITS["reverse_bias"][1])
    ax.set_ylabel(S.UNITS["cce"][1])
    _validity_bands(ax, res, v)
    _legend(ax, loc="lower right")
    _common_notes(res, notes)


def _cce_let(fig, table, notes, **_):
    ax = fig.add_subplot()
    let = _arr(table["LET_keV_um"])
    cce = _arr(table["CCE"])
    conv = np.asarray(table.get("converged", [True] * let.size), dtype=bool)
    ax.semilogx(let[conv], cce[conv], "o-", ms=2.5, label="converged")
    if (~conv).any():
        ax.semilogx(
            let[~conv], cce[~conv], "x", color=S.INVALID_COLOR, label="not converged"
        )
    ax.set_xlabel(S.UNITS["let"][1])
    ax.set_ylabel(S.UNITS["cce"][1])
    ax.set_title("CCE vs LET (single particle)")
    ax.ticklabel_format(axis="y", useOffset=False)
    _legend(ax)
    notes.append(NOTES["cce_let"])


def _dark(fig, res, notes, **_):
    ax = fig.add_subplot()
    md = res.metadata
    y = np.abs(_arr(res.y))
    if res.sim_type == "dark_current_vs_T":
        ax.semilogy(_arr(res.x), y, "o-", ms=2.5)
        ax.set_xlabel(S.UNITS["temperature"][1])
        ax.set_title(f"Dark current vs T at {md.get('V_bias', '')} V")
    else:
        keep, v = _rev(res)
        yk = y[keep]
        past = np.asarray(
            md.get("past_reach_through", np.zeros(len(res.x), bool)), dtype=bool
        )
        past = past[keep] if past.size == keep.size else np.zeros(v.size, bool)
        ax.semilogy(v, yk, "-", color=S.PALETTE[0])
        ax.semilogy(
            v[~past], yk[~past], "o", ms=2.5, color=S.PALETTE[0], label="I_total"
        )
        if past.any():
            ax.semilogy(
                v[past],
                yk[past],
                "o",
                ms=2.5,
                mfc="none",
                color=S.PALETTE[0],
                label="past reach-through",
            )
        ax.set_xlabel(S.UNITS["reverse_bias"][1])
        ax.set_title(f"Dark current at {getattr(res.config, 'T', 300):g} K")
        _validity_bands(ax, res, v)
        if md.get("validity_note"):
            notes.append(f"Validity: {md['validity_note']}")
    ax.set_ylabel(S.UNITS["current"][1])
    _legend(ax, loc="lower right")
    notes.insert(0, NOTES["dark_current"])
    _common_notes(res, notes)


def _flash(fig, res, notes, **_):
    ax = fig.add_subplot()
    md = res.metadata
    ax.plot(_arr(res.x), _arr(res.y), "o-", ms=2.5, label="steady state")
    ax.set_xlabel(S.UNITS["dose_rate"][1])
    ax.set_ylabel(S.UNITS["cce"][1])
    ax.ticklabel_format(axis="y", useOffset=False)
    ax.set_title(f"FLASH CCE at {md.get('V_bias', '')} V, {md.get('E_MeV', '')} MeV")
    notes.insert(0, NOTES["flash"])
    _common_notes(res, notes)


def _transient(fig, res, notes, **_):
    ax = fig.add_subplot()
    md = res.metadata
    ax.plot(_arr(res.x), _arr(res.y), "o-", ms=2.5, label="CCE (clipped to [0, 1])")
    raw = md.get("transient_cce_raw")
    if raw is not None and np.size(raw) == len(res.x):
        ax.plot(_arr(res.x), _arr(raw), "s--", ms=1.8, label="CCE raw")
    ax.set_xlabel(S.UNITS["dose_rate"][1])
    ax.set_ylabel(S.UNITS["cce"][1])
    ax.ticklabel_format(axis="y", useOffset=False)
    ax.set_title(f"FLASH transient CCE at {md.get('V_bias', '')} V")
    _legend(ax)
    notes.insert(0, NOTES["transient"])
    _common_notes(res, notes)


def _temperature(fig, res, notes, **_):
    ax = fig.add_subplot()
    md = res.metadata
    ax.plot(_arr(res.x), _arr(res.y), "o-", ms=2.5)
    ax.set_xlabel(S.UNITS["temperature"][1])
    ax.set_ylabel(S.UNITS["cce"][1])
    ax.ticklabel_format(axis="y", useOffset=False)
    ax.set_title(f"CCE vs T at {md.get('voltage_V', '')} V ({md.get('method', '')})")
    notes.insert(0, NOTES["temperature"])
    _common_notes(res, notes)


def _damage(fig, res, notes, **_):
    ax = fig.add_subplot()
    md = res.metadata
    x, y = _arr(res.x), _arr(res.y)
    valid = np.asarray(md.get("valid", np.isfinite(y)), dtype=bool)
    ax.semilogx(x[valid], y[valid], "o-", ms=2.5, label="solved")
    if (~valid).any():
        from matplotlib.transforms import blended_transform_factory

        tr = blended_transform_factory(ax.transData, ax.transAxes)
        ax.plot(
            x[~valid],
            np.full((~valid).sum(), 0.04),
            "x",
            color=S.INVALID_COLOR,
            transform=tr,
            label="not solved (no CCE)",
        )
        reasons = [r for r, ok in zip(md.get("reason") or [], valid) if not ok and r]
        if reasons:
            notes.append(f"{len(reasons)} point(s) not solved, e.g.: {reasons[0]}")
    phi_c = md.get("compensation_fluence")
    if phi_c:
        ax.axvline(
            phi_c,
            ls=":",
            color="#555555",
            lw=0.9,
            label=f"compensation {phi_c:.2g} cm$^{{-2}}$",
        )
    ax.set_ylim(0, 1.08)
    ax.set_xlabel(S.UNITS["fluence"][1])
    ax.set_ylabel(S.UNITS["cce"][1])
    ax.set_title(
        f"CCE vs fluence ({md.get('energy_MeV', '')} MeV p, {md.get('V_bias', '')} V)"
    )
    _legend(ax, loc="lower left")
    notes.insert(0, NOTES["damage"])
    _common_notes(res, notes)


def _microdosimetry(fig, res, notes, **_):
    ax = fig.add_subplot()
    md = res.metadata
    ax.semilogx(_arr(res.x), _arr(res.y), "-")
    ax.set_xlabel(S.UNITS["lineal"][1])
    ax.set_ylabel(r"$y\,d(y)$")
    ax.set_title(
        rf"Microdosimetric spectrum ($y_F$ {md.get('y_F', float('nan')):.3g}, "
        rf"$y_D$ {md.get('y_D', float('nan')):.3g} keV/$\mu$m)"
    )
    notes.append(NOTES["microdosimetry"])


_MAP_TITLES = {
    "map_potential": "Potential",
    "map_field": "|E|",
    "map_doping": "Net doping",
    "map_electrons": "Electron density",
    "map_holes": "Hole density",
}
_MAP_SPEC = {
    "map_potential": ("Potential", "Potential (V)", "RdBu_r", "lin"),
    "map_field": ("ElectricField", "|E| (kV/cm)", "viridis", "field"),
    "map_doping": (
        "NetDoping",
        r"log$_{10}$|N$_D$ - N$_A$| (cm$^{-3}$)",
        "plasma",
        "log",
    ),
    "map_electrons": ("Electrons", r"log$_{10}$ n (cm$^{-3}$)", "magma", "carrier"),
    "map_holes": ("Holes", r"log$_{10}$ p (cm$^{-3}$)", "magma", "carrier"),
}


def _map2d(fig, res, notes, kind="map_field", **_):
    from matplotlib.tri import Triangulation

    from etna.viz.maps import CARRIER_FLOOR_CM3

    key, clabel, cmap, scale = _MAP_SPEC[kind]
    nv = res.mesh.node_values
    if key not in nv:
        raise ValueError(
            f"{key} not available from solver {res.metadata.get('solver')}"
        )
    x = _arr(res.mesh.x_coords) * 1e4
    y = _arr(res.mesh.y_coords) * 1e4
    z = _arr(nv[key])
    if scale == "field":
        z = z * 1e-3
    elif scale == "log":
        z = np.log10(np.maximum(np.abs(z), 1.0))
    elif scale == "carrier":
        z = np.log10(np.maximum(z, CARRIER_FLOOR_CM3))
        notes.append(
            f"Carrier densities below {CARRIER_FLOOR_CM3:g} cm^-3 clipped (numerical floor; "
            "engine2d vs devsim parity <= 1.5e-13 above it)."
        )
    from matplotlib.ticker import MaxNLocator

    ax = fig.add_subplot()
    ax.grid(False)
    tri = Triangulation(x, y)
    if scale == "field":  # the abrupt-step spikes would swamp the colour scale
        from etna.api.results import INTERFACE_SPIKE_WINDOW_UM

        steps = [res.config.substrate_thickness_um]
        if res.config.nplus_thickness_um:
            steps.append(res.config.substrate_thickness_um + res.config.epi_thickness_um)
        win = np.ones_like(y, bool)
        for s in steps:
            win &= np.abs(y - s) > INTERFACE_SPIKE_WINDOW_UM
        if win.any() and not win.all():
            z = np.minimum(z, float(np.nanmax(z[win])))
            notes.append("Colour scale clipped: nm-scale spikes at the abrupt doping steps not shown.")
    zmin, zmax = float(np.nanmin(z)), float(np.nanmax(z))
    if zmax <= zmin:
        zmax = zmin + 1.0
    levels = MaxNLocator(24).tick_values(zmin, zmax)
    cs = ax.tricontourf(tri, z, levels=levels, cmap=cmap, extend="neither")
    lines = ax.tricontour(tri, z, levels=MaxNLocator(8).tick_values(zmin, zmax), colors="k", linewidths=0.35)
    ax.clabel(lines, fontsize=S.BASE_FONT_PT - 3, fmt="%.3g")
    for c in res.mesh.contacts or []:
        ax.axhline(c["position"] * 1e4, color="#222222", lw=1.2)
    ax.set_ylim(y.max(), y.min())  # depth downwards
    ax.set_xlabel(S.UNITS["lateral"][1])
    ax.set_ylabel(S.UNITS["depth"][1])
    bias = res.metadata.get("bias_V")
    name = _MAP_TITLES[kind]
    ax.set_title(f"{name} at {bias:g} V" if bias is not None else name)
    cb = fig.colorbar(cs, ax=ax)
    cb.set_label(clabel)
    _common_notes(res, notes)


def _stack(fig, cfg, notes):
    from matplotlib.patches import Rectangle

    from etna.export.cross_section import stack_layers

    layers = stack_layers(cfg)
    total = layers[-1]["bottom_um"]
    ax, axt = fig.subplots(1, 2, gridspec_kw={"width_ratios": [1, 3.2], "wspace": 0}, sharey=True)
    for L in layers:
        ax.add_patch(
            Rectangle(
                (0, L["top_um"]),
                1.0,
                L["thickness_um"],
                facecolor=S.LAYER_COLORS.get(L["key"], "#cccccc"),
                edgecolor="k",
                lw=0.5,
            )
        )
    # label positions: layer centres, pushed apart so thin layers (p+) stay readable
    centres = [0.5 * (L["top_um"] + L["bottom_um"]) for L in layers]
    gap = 0.18 * total
    ys = list(centres)
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + gap)
    ys[-1] = min(ys[-1], 0.97 * total)
    for i in range(len(ys) - 2, -1, -1):
        ys[i] = min(ys[i], ys[i + 1] - gap)
    ys = [max(y, 0.03 * total) for y in ys]
    for L, yc, yl in zip(layers, centres, ys):
        axt.plot([0.0, 0.06], [yc, yl], color="#333333", lw=0.4)
        txt = f"{L['name']}: {L['thickness_um']:.3g} " + r"$\mu$m" + "\n" + "\n".join(
            textwrap.wrap(L["doping"], 52)
        )
        axt.text(0.08, yl, txt, va="center", fontsize=S.BASE_FONT_PT - 2.5, linespacing=1.1)
    ax.set_xlim(0, 1)
    axt.set_xlim(0, 1)
    ax.set_ylim(total, 0)
    for a in (ax, axt):
        a.set_xticks([])
        a.grid(False)
        a.tick_params(axis="x", which="both", top=False, bottom=False)
        for s in ("top", "right", "bottom"):
            a.spines[s].set_visible(False)
    ax.tick_params(axis="y", which="both", right=False)
    axt.axis("off")
    ax.set_ylabel(S.UNITS["depth"][1])
    fig.suptitle(f"Device stack (depth to scale; area {cfg.area_cm2 * 100:g} mm2)")
    notes.append(
        "Lateral width schematic, not to scale (only the active area is known: die, "
        "guard ring and pads are not modelled)."
    )


def _tilted_alpha_shape(fig, res, notes, **_):
    ax = fig.add_subplot()
    md = res.metadata
    keep, v = _rev(res)
    ax.plot(
        v,
        _arr(res.y)[keep],
        "-",
        color=S.PALETTE[0],
        label="model, normalised at -50 V",
    )
    dv, dc, de = _arr(md["data_V"]), _arr(md["data_cce"]), _arr(md["data_err"])
    ax.errorbar(
        -dv,
        dc,
        yerr=de,
        fmt="o",
        ms=3,
        color=S.DATA_COLOR,
        capsize=1.5,
        lw=0.7,
        label="Micro Fig. 4b (digitised)",
    )
    v_rt = md.get("reach_through_V")
    if v_rt is not None:
        ax.axvline(
            v_rt, ls=":", color="#555555", lw=0.9, label=f"full depletion {v_rt:.2f} V"
        )
    ax.set_xlabel(S.UNITS["reverse_bias"][1])
    ax.set_ylabel("CCE / CCE(-50 V)")
    ax.ticklabel_format(axis="y", useOffset=False)
    ax.set_title("Tilted alpha: shape only (not reproduced)")
    _legend(ax, loc="lower right")
    notes.append(md["verdict"])
    notes.append(
        "Data: Petringa et al. 2025, Phys. Med. 140, 105675, Fig. 4b, digitised in this repo "
        "(read-off estimates, +-0.0005; lab/microdosimeter/notes/LOG.md section 1). Absolute "
        "CCE not compared (calibration ambiguity ~8 %)."
    )


def _w_creep(fig, res, notes, **_):
    ax = fig.add_subplot()
    md = res.metadata
    keep, v = _rev(res)
    f = 1e4
    ax.plot(
        v,
        _arr(md["depletion_widths"])[keep] * f,
        "-",
        color=S.PALETTE[0],
        label="depletion edge (n/N_D = 0.5)",
    )
    ax.plot(
        v,
        _arr(md["W_from_C"])[keep] * f,
        "--",
        color=S.PALETTE[1],
        label=r"$\varepsilon A/C$ (small-signal)",
    )
    ax.errorbar(
        -_arr(md["label_V"]),
        _arr(md["label_W_um"]),
        yerr=_arr(md["label_err_um"]),
        fmt="o",
        ms=3,
        color=S.DATA_COLOR,
        capsize=1.5,
        lw=0.7,
        label="Photons Fig. 4b labels",
    )
    _epi_marks(ax, res.config)
    ax.set_xlabel(S.UNITS["reverse_bias"][1])
    ax.set_ylabel(S.UNITS["width_um"][1])
    ax.set_title(r"W creep -10 $\to$ -30 V vs Photons Fig. 4b")
    _legend(ax, loc="lower right", fontsize=S.BASE_FONT_PT - 2)
    notes.append(md["verdict"])


_BUILDERS = {
    "cv": _cv,
    "mott_schottky": _mott,
    "depletion": _depletion,
    "field": _field,
    "cce": _cce,
    "dark_current": _dark,
    "flash": _flash,
    "transient": _transient,
    "temperature": _temperature,
    "damage": _damage,
    "microdosimetry": _microdosimetry,
    "tilted_alpha_shape": _tilted_alpha_shape,
    "w_creep": _w_creep,
}
for _k in MAP_KINDS:
    _BUILDERS[_k] = (
        lambda k: lambda fig, res, notes, **o: _map2d(fig, res, notes, kind=k, **o)
    )(_k)
