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
from etna.api.run_all import SECTION_NAMES, _sections
from app.components.page_utils import require_config, show_error, show_stale_warning


def _curve(sec) -> go.Figure:
    r = sec.result
    fig = go.Figure(go.Scatter(x=r.x, y=r.y, mode="lines+markers"))
    fig.update_layout(title=sec.title, xaxis_title="x", yaxis_title="y")
    return fig


def render() -> None:
    st.title("Run All")
    st.write(
        "Runs every simulation for the device in the sidebar and reports one "
        "table. Sections that fail or do not apply (the 1D-only ones on a 2D "
        "device) are listed with the reason; the rest still run."
    )
    cfg = require_config()

    chosen = st.multiselect(
        "Sections",
        options=list(SECTION_NAMES),
        default=list(SECTION_NAMES),
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
    if n_bad:
        st.warning(f"{n_bad} section(s) failed; see the reason column.")
    st.dataframe(pd.DataFrame(report.rows()), hide_index=True)
    for sec in report.sections:
        if sec.result is not None and sec.status == "ok":
            with st.expander(f"{sec.title} ({sec.solver})"):
                st.plotly_chart(_curve(sec))
    buf = io.StringIO()
    pd.DataFrame(report.rows()).to_csv(buf, index=False)
    c1, c2 = st.columns(2)
    c1.download_button("Download summary CSV", buf.getvalue(), "run_all_summary.csv", "text/csv")
    c2.download_button("Download full JSON", report.to_json(), "run_all.json", "application/json")
