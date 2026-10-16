"""C-V Analysis page: 1D guard, Run -> cache -> render -> download.

All Streamlit calls live inside render() -- no module-level side effects --
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests
(tests/test_app_cv_page.py). etna.run_cv is referenced as a MODULE
ATTRIBUTE (not `from etna import run_cv`) so tests can intercept it via
monkeypatch.setattr(etna, "run_cv", fake).
"""

from __future__ import annotations

import streamlit as st

import etna
from app.components.page_utils import (
    require_config,
    show_error,
    show_result_notes,
    show_stale_warning,
)
from app.components.results import (
    build_cv_figure,
    build_mott_schottky_figure,
    to_csv_bytes,
)


def render() -> None:
    st.title("C-V Analysis")

    cfg = require_config()

    if cfg.half_width_um is not None:
        st.warning(
            "This workflow is 1D-only (2D devices are shown on the Field Map "
            "page). Set Dimensionality to 1D in the sidebar."
        )
        st.stop()

    col1, col2 = st.columns(2)
    with col1:
        v_stop = st.number_input(
            "Max reverse bias V_stop (V)",
            value=-50.0,
            min_value=-500.0,
            max_value=-0.5,
            step=5.0,
            key="cv_v_stop",
            help="The sweep runs from 0 V to V_stop (negative = reverse bias).",
        )
    with col2:
        n_points = st.number_input(
            "Number of bias points",
            value=40,
            min_value=2,
            max_value=400,
            step=1,
            key="cv_n_points",
        )

    inputs = {"V_stop": float(v_stop), "number of points": int(n_points)}
    if st.button("Run simulation"):
        try:
            with st.spinner("Running C-V sweep..."):
                st.session_state["cv_result"] = etna.run_cv(
                    cfg, v_start=0.0, v_stop=float(v_stop), n_points=int(n_points)
                )
            st.session_state["cv_inputs"] = inputs
        except Exception as e:  # incl. devsim.error (not a RuntimeError)
            show_error(
                e,
                "Try a smaller reverse bias or adjust epi thickness / doping in "
                "the sidebar.",
            )

    result = st.session_state.get("cv_result")
    if result is not None:
        show_stale_warning(
            result.config, cfg, st.session_state.get("cv_inputs"), inputs
        )
        show_result_notes(result)
        st.plotly_chart(build_cv_figure(result))
        st.plotly_chart(build_mott_schottky_figure(result))
        st.download_button(
            "Download CSV",
            data=to_csv_bytes(result),
            file_name="cv_result.csv",
            mime="text/csv",
        )

        with st.expander("Depletion width vs V"):
            import plotly.graph_objects as go

            depletion_um = result.metadata["depletion_widths"] * 1e4
            fig = go.Figure(
                data=go.Scatter(x=result.x, y=depletion_um, mode="lines+markers")
            )
            fig.update_layout(
                title="Depletion Width vs Bias",
                xaxis_title="Voltage (V)",
                yaxis_title="Depletion Width (µm)",
            )
            st.plotly_chart(fig)
