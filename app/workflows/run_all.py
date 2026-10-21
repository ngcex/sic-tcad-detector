"""Run All page: every simulation for the sidebar device, in one click.

Calls ``etna.run_all`` once and shows a summary table (status, solver, key
figures, reason, wall time), one expander per section with its primary curve,
and a combined JSON / CSV download. A section that fails or does not apply to
the configuration is a labelled row; the others still run. The report is cached
against the device configuration (the stale-result warning applies).
"""

from __future__ import annotations

import io

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import etna
from etna.api.run_all import DEFAULT_SECTIONS, SECTION_NAMES, _sections
from app.components.page_utils import (
    about,
    require_config,
    show_error,
    show_stale_warning,
)
from app.components.results import build_result_figure


def _curve(sec) -> "go.Figure | None":
    """The section's primary curve with labelled axes (None: no 1D curve)."""
    return build_result_figure(sec.result, sec.title)


def render() -> None:
    st.title("Run All")
    st.write(
        "Runs every simulation for the device in the sidebar and reports one "
        "table. Sections that fail or do not apply (the 1D-only ones on a 2D "
        "device) are listed with the reason; the rest still run."
    )
    about(
        "Each section runs its facade with the library defaults (the single-run pages "
        "let you change them):\n"
        "- **C-V** 0 → −50 V; **CCE** −10 → −40 V (Am-241 alpha); **Field** at −50 V.\n"
        "- **Radiation damage**: 5.6 MeV protons, 8 fluences from 0.01× to 3× the "
        "compensation fluence; the points beyond it are refused (NaN with a reason) "
        "unless charged deep traps are on.\n"
        "- **Dark current** 0 → −100 V: calibrated budget estimate (N_t fitted to "
        "< 18 pA at −60 V).\n"
        "- **Temperature sweep**: Hecht CCE at −30 V, 250-350 K.\n"
        "- **FLASH**: steady-state CCE vs dose rate at −30 V, 62 MeV, SRH + Auger; a "
        "numerical sensitivity bound, not a validated plasma-recombination "
        "prediction.\n"
        "- **Microdosimetry sweep** (opt-in): 2D engine only, half-widths 25/50 µm at "
        "20/40/60 V.\n"
        "Status 'warning' = the section ran but returned no valid point."
    )
    cfg = require_config()

    chosen = st.multiselect(
        "Sections",
        options=list(SECTION_NAMES),
        default=list(DEFAULT_SECTIONS),
        format_func=lambda n: _sections()[n][0],
        key="ra_sections",
    )
    if st.button("Run all", disabled=not chosen):
        try:
            with st.spinner("Running all sections..."):
                st.session_state["run_all_report"] = etna.run_all(cfg, tuple(chosen))
        except Exception as e:
            show_error(e)

    report = st.session_state.get("run_all_report")
    if report is None:
        return
    show_stale_warning(report.config, cfg, None, None)
    n_bad = sum(s.status == "error" for s in report.sections)
    st.caption(f"{len(report.sections)} sections, {report.wall_s:.1f} s total")
    n_warn = sum(s.status == "warning" for s in report.sections)
    if n_bad:
        st.warning(f"{n_bad} section(s) failed; see the reason column.")
    if n_warn:
        st.warning(f"{n_warn} section(s) ran but returned no valid point; see the reason column.")
    st.dataframe(pd.DataFrame(report.rows()), hide_index=True)
    for sec in report.sections:
        if sec.result is not None and sec.status == "ok":
            fig = _curve(sec)
            if fig is None:  # 2D field: no 1D curve (see the Field Map page)
                continue
            with st.expander(f"{sec.title} ({sec.solver})"):
                st.plotly_chart(fig, key=f"ra_fig_{sec.name}")
    buf = io.StringIO()
    pd.DataFrame(report.rows()).to_csv(buf, index=False)
    c1, c2 = st.columns(2)
    c1.download_button("Download summary CSV", buf.getvalue(), "run_all_summary.csv", "text/csv")
    c2.download_button("Download full JSON", report.to_json(), "run_all.json", "application/json")
