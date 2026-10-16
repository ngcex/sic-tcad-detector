"""Charge Collection (CCE) page: 1D guard, Run -> cache -> render -> download.

All Streamlit calls live inside render() -- no module-level side effects --
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests
(tests/test_app_cce_page.py). etna.run_cce is referenced as a MODULE
ATTRIBUTE (not `from etna import run_cce`) so tests can intercept it via
monkeypatch.setattr(etna, "run_cce", fake).
"""

from __future__ import annotations

import streamlit as st

import etna
from app.components.page_utils import (
    show_error,
    show_result_notes,
    show_stale_warning,
)
from app.components.results import build_cce_figure, to_csv_bytes

FORWARD_BIAS_ERROR = (
    "CCE is defined here for reverse bias only: both V_start and V_stop must be "
    "≤ 0 V (negative = reverse bias)."
)


def render() -> None:
    st.title("Charge Collection (CCE)")

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    if cfg.half_width_um is not None:
        st.warning(
            "This workflow is 1D-only (2D devices are shown on the Field Map "
            "page). Set Dimensionality to 1D in the sidebar."
        )
        st.stop()

    col1, col2, col3 = st.columns(3)
    with col1:
        v_start = st.number_input(
            "Start bias V_start (V)",
            value=-10.0,
            key="cce_v_start",
        )
    with col2:
        v_stop = st.number_input(
            "End bias V_stop (V)",
            value=-40.0,
            key="cce_v_stop",
        )
    with col3:
        n_points = st.number_input(
            "Number of bias points",
            value=30,
            min_value=2,
            max_value=400,
            step=1,
            key="cce_n_points",
        )
    st.caption("Reverse bias is negative (e.g. −40 V means 40 V reverse bias).")

    if st.button("Run simulation"):
        if v_start > 0 or v_stop > 0:
            st.error(FORWARD_BIAS_ERROR)
        else:
            try:
                with st.spinner("Running CCE sweep..."):
                    st.session_state["cce_result"] = etna.run_cce(
                        cfg,
                        v_start=float(v_start),
                        v_stop=float(v_stop),
                        n_points=int(n_points),
                    )
            except Exception as e:  # incl. devsim.error (not a RuntimeError)
                show_error(
                    e,
                    "Try a smaller reverse bias or adjust epi thickness / doping "
                    "in the sidebar.",
                )

    result = st.session_state.get("cce_result")
    if result is not None:
        show_stale_warning(result.config, cfg)
        show_result_notes(result)
        st.plotly_chart(build_cce_figure(result))
        st.download_button(
            "Download CSV",
            data=to_csv_bytes(result),
            file_name="cce_result.csv",
            mime="text/csv",
        )
        st.caption(
            "Generated current (alpha-like profile, 15 µm range): "
            f"{float(result.metadata.get('I_generated', float('nan'))):.3e} A/cm²"
        )
