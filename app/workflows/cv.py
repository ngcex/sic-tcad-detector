"""C-V Analysis page: 1D guard, Run -> cache -> render -> download.

All Streamlit calls live inside render() -- no module-level side effects --
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests
(tests/test_app_cv_page.py). petringa.run_cv is referenced as a MODULE
ATTRIBUTE (not `from petringa import run_cv`) so tests can intercept it via
monkeypatch.setattr(petringa, "run_cv", fake) -- the seam proven in
tests/test_app_run_mockability.py (39-01).
"""

from __future__ import annotations

import streamlit as st

import petringa
from app.components.results import (
    build_cv_figure,
    build_mott_schottky_figure,
    to_csv_bytes,
)


def render() -> None:
    st.title("C-V Analysis")

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    if cfg.half_width_um is not None:
        st.warning(
            "These workflows are 1D-only. 2D field visualization arrives in "
            "Phase 40 (geometry viewer). Set Dimensionality to 1D in the "
            "sidebar."
        )
        st.stop()

    if st.button("Run simulation"):
        try:
            st.session_state["cv_result"] = petringa.run_cv(cfg)
        except RuntimeError as e:
            st.error(
                f"Simulation failed to converge: {e}\n\n"
                "Try adjusting device parameters (e.g. epi thickness, doping) "
                "in the sidebar and running again."
            )

    result = st.session_state.get("cv_result")
    if result is not None:
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
