"""Shared helpers of the opt-in model pages (breakdown, MC deposition, track engine,
trap kinetics, 2D structures): publication figures in the ``etna.viz`` house style,
a provenance-headed CSV writer, and a small per-session result cache.

``etna.viz.figure_for`` knows only the facade kinds (C-V, CCE, field, ...). The pages
here draw their own matplotlib figures under the same rc set (``etna.viz.style.RC``),
sizes, palette and provenance / validation footer, and export them through
``etna.viz.export_figure`` (PNG / SVG / PDF at a chosen DPI). Nothing here simulates.
"""

from __future__ import annotations

import textwrap
from datetime import datetime, timezone

import matplotlib as mpl
import numpy as np
import pandas as pd
import streamlit as st
from matplotlib.figure import Figure
from matplotlib.layout_engine import ConstrainedLayoutEngine

import etna
from etna.viz import style as S
from etna.viz.provenance import config_hash, footer_text, provenance

MIME = {"png": "image/png", "svg": "image/svg+xml", "pdf": "application/pdf"}
CACHE_SIZE = 8
STALE_TEXT = (
    "These results were computed for a different device or different page inputs; "
    "press the run button to update."
)


# ----------------------------------------------------------------- figures (pure)
def publication_figure(draw, prov: dict, notes=(), size=S.DEFAULT_SIZE) -> Figure:
    """A house-style matplotlib Figure: ``draw(fig)`` adds the axes; the footer carries
    ``notes``, the provenance and the validation label (as ``etna.viz.figure_for``)."""
    from etna.export.labels import VALIDATION_LABEL

    with mpl.rc_context(S.RC):
        fig = Figure(figsize=S.size_inches(size), layout=ConstrainedLayoutEngine())
        draw(fig)
        fig.etna_provenance = prov
        w_in = fig.get_figwidth()
        chars = max(40, int(w_in * 72 / (0.62 * S.FOOTER_FONT_PT)))
        lines = []
        for n in notes:
            lines += textwrap.wrap(str(n), chars)
        lines += textwrap.wrap(footer_text(prov), chars)
        clauses = VALIDATION_LABEL.split("; ")
        for i, c in enumerate(clauses):
            lines += textwrap.wrap(c + (";" if i < len(clauses) - 1 else ""), chars)
        fig.text(
            0.008, 0.006, "\n".join(lines), fontsize=S.FOOTER_FONT_PT,
            color="#444444", ha="left", va="bottom", linespacing=1.15,
        )  # fmt: skip
        h_pt = len(lines) * S.FOOTER_FONT_PT * 1.2 + 4.0
        fig._etna_footer_h_pt = h_pt  # read by etna.viz.export_figure on resize
        frac = min(0.6, h_pt / 72.0 / fig.get_figheight())
        fig.get_layout_engine().set(rect=(0.0, frac, 1.0, 1.0 - frac))
    return fig


def figure_bytes(build, fmt="png", dpi=S.DEFAULT_DPI, size=S.DEFAULT_SIZE) -> bytes:
    """``build(size)`` -> Figure, exported by ``etna.viz.export_figure``."""
    return etna.viz.export_figure(build(size), fmt, dpi=dpi)


# ----------------------------------------------------------------- CSV (pure)
def table_csv_bytes(title: str, columns: dict, header=(), config=None) -> bytes:
    """CSV with a '#' header (title, version, date, device fields, ``header`` lines)."""
    from dataclasses import asdict, is_dataclass

    lines = [
        f"# ETNA SiC TCAD Simulator — {title}",
        f"# software_version: {etna.__version__}",
        f"# generated: {datetime.now(timezone.utc).isoformat()}",
    ]
    if config is not None and is_dataclass(config):
        lines.append(
            "# device: " + ", ".join(f"{k}={v}" for k, v in asdict(config).items())
        )
    lines += [f"# {h}" for h in header]
    df = pd.DataFrame({k: np.asarray(v) for k, v in columns.items()})
    return ("\n".join(lines) + "\n" + df.to_csv(index=False)).encode("utf-8")


# ----------------------------------------------------------------- session cache
def token(cfg, inputs: dict) -> tuple:
    """Hashable key of (device, page inputs)."""
    h = config_hash(cfg) if cfg is not None else ""
    return (h, tuple(sorted((k, repr(v)) for k, v in inputs.items())))


def cached_run(page: str, tok: tuple, compute):
    """Result of ``compute()`` for ``tok``, reused within the session (last
    ``CACHE_SIZE`` keys). Exceptions are not cached."""
    cache = st.session_state.setdefault(f"{page}_cache", {})
    if tok in cache:
        cache[tok] = cache.pop(tok)  # most recent last
        return cache[tok]
    out = compute()
    cache[tok] = out
    while len(cache) > CACHE_SIZE:
        cache.pop(next(iter(cache)))
    return out


# ----------------------------------------------------------------- figure widget
def figure_block(key: str, stem: str, build, tok=None, prov_hint: str = "") -> None:
    """Publication preview + download of one figure. ``build(size)`` -> Figure; it is
    called only when the preview is on or a download is requested. ``tok`` identifies
    the plotted result (figure bytes are reused while it is unchanged)."""
    with st.expander("Publication figure (etna.viz style): preview and download"):
        c = st.columns(3)
        fmt = c[0].selectbox("Format", ["png", "svg", "pdf"], key=f"{key}_fmt")
        dpi = c[1].number_input(
            "Resolution (dpi)", 72, 1200, S.DEFAULT_DPI, 50, key=f"{key}_dpi",
            help="Raster resolution of the PNG; SVG / PDF are vector.",
        )  # fmt: skip
        size = c[2].selectbox(
            "Size", list(S.SIZES_MM), key=f"{key}_size",
            help="single column 89 x 70 mm, double column 183 x 105 mm, slide 254 x 143 mm",
        )  # fmt: skip
        if prov_hint:
            st.caption(prov_hint)
        if st.toggle("Show preview (150 dpi PNG)", value=False, key=f"{key}_prev"):
            st.image(figure_bytes(build, "png", 150, size))
        store = st.session_state.setdefault("pubfig_bytes", {})

        def _data(build=build, fmt=fmt, dpi=int(dpi), size=size):
            k = (key, tok, fmt, dpi, size)
            if k not in store:
                store.clear()  # keep one figure's bytes per session
                store[k] = figure_bytes(build, fmt, dpi, size)
            return store[k]

        st.download_button(
            f"Download figure ({fmt.upper()})",
            data=_data,
            file_name=f"etna_{stem}.{fmt}",
            mime=MIME[fmt],
            key=f"{key}_dl",
            on_click="ignore",
        )


def prov_for(result=None, config=None, solver=None, source=None) -> dict:
    return provenance(result, config=config, solver=solver, source=source)
