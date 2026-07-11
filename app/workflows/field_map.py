"""Field Map page: run_field -> cache -> render (E-field + potential) -> download.

All Streamlit calls live inside render() — no module-level side effects —
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests.

The 2D guard here is CRITICAL and subtly different from the C-V/CCE pages:
`petringa.run_field` does NOT raise for a 2D config — it silently returns
empty x/y arrays (2D x/y are lateral coordinates, not a well-defined depth
profile) — so the guard MUST be a pre-check before calling run_field, never
a try/except around the call (an exception handler would catch nothing and
silently plot empty arrays). The 2D geometry viewer for this page arrives in
Phase 40.
"""

from __future__ import annotations

import streamlit as st

import petringa
from app.components.results import build_field_figures, to_csv_bytes


def render() -> None:
    st.title("Field Map")

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
            st.session_state["field_result"] = petringa.run_field(cfg)
        except RuntimeError as e:
            st.error(
                f"Simulation failed to converge: {e}\n\n"
                "Try adjusting device parameters (e.g. epi thickness, doping) "
                "in the sidebar and running again."
            )

    result = st.session_state.get("field_result")
    if result is not None:
        efield_fig, potential_fig = build_field_figures(result)
        st.plotly_chart(efield_fig)
        st.plotly_chart(potential_fig)

        st.download_button(
            "Download CSV",
            data=to_csv_bytes(result),
            file_name="field_result.csv",
            mime="text/csv",
        )

        with st.expander("Net doping vs depth"):
            import plotly.graph_objects as go

            net_doping = result.metadata["net_doping"]
            doping_fig = go.Figure(
                data=go.Scatter(x=result.x, y=net_doping, mode="lines")
            )
            doping_fig.update_layout(
                title="Net Doping vs Depth",
                xaxis_title="Depth (µm)",
                yaxis_title="Net Doping (cm⁻³)",
                yaxis_type="log",
            )
            st.plotly_chart(doping_fig)
