"""Field Map page: run_field -> cache -> dimensionality-branched render.

All Streamlit calls live inside render() — no module-level side effects —
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests.

This page routes BOTH 1D and 2D configs THROUGH `etna.run_field` (there is
no 2D pre-check that short-circuits before the solve). The render then branches
on `result.mesh.y_coords is None`:

- 1D (y_coords is None): the existing E-field/potential line charts + CSV
  download + net-doping expander (which read `result.x`/`result.y`/`metadata`,
  populated only for 1D) render as before, AND a geometry-viewer bar is added
  below them (a supplement, not a replacement).
- 2D (y_coords is not None): `result.x`/`result.y` are empty and `to_csv_bytes`
  has no 2D branch, so the line charts and CSV download are skipped; only the
  geometry-viewer heatmap renders.

The `try/except RuntimeError` around `run_field` is kept because a 2D
`ramp_bias` may not converge (documented upstream blocker); the geometry viewer
renders only when `result.mesh is not None`. A quantity `st.selectbox` (fixed
`QUANTITIES` options, persistent key) feeds `build_geometry_figure`; changing it
reruns the page, reads the cached `session_state["field_result"]`, and re-runs
only the griddata interpolation — never the devsim solve.
"""

from __future__ import annotations

import streamlit as st

import etna
from app.components.geometry_viewer import build_geometry_figure, QUANTITIES
from app.components.results import build_field_figures, to_csv_bytes


def render() -> None:
    st.title("Field Map")

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    if st.button("Run simulation"):
        try:
            st.session_state["field_result"] = etna.run_field(cfg)
        except RuntimeError as e:
            st.error(
                f"Simulation failed to converge: {e}\n\n"
                "Try adjusting device parameters (e.g. epi thickness, doping) "
                "in the sidebar and running again."
            )

    result = st.session_state.get("field_result")
    if result is not None and result.mesh is not None:
        if result.mesh.y_coords is None:
            # 1D branch: existing line charts + CSV + expander (read result.x/y),
            # THEN the supplemental geometry-viewer bar below them.
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

        # Geometry viewer (both branches): quantity dropdown feeds the pure
        # builder. On a selectbox change Streamlit reruns, the Run button
        # returns False, the cached result is read, and only griddata re-runs.
        quantity = st.selectbox(
            "Quantity", list(QUANTITIES.keys()), index=0, key="geo_quantity"
        )
        st.plotly_chart(build_geometry_figure(result.mesh, quantity))
