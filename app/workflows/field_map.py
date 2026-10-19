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

The `try/except Exception` around `run_field` catches solver non-convergence
(RuntimeError or devsim.error) and shows a friendly error; the geometry viewer
renders only when `result.mesh is not None`. A quantity `st.selectbox` (fixed
`QUANTITIES` options, persistent key) feeds `build_geometry_figure`; changing it
reruns the page, reads the cached `session_state["field_result"]`, and re-runs
only the griddata interpolation — never the devsim solve.
"""

from __future__ import annotations

import streamlit as st

import etna
from app.components.geometry_viewer import build_geometry_figure, QUANTITIES
from app.components.page_utils import (
    require_config,
    show_error,
    show_result_notes,
    show_stale_warning,
)
from app.components.results import build_field_figures, to_csv_bytes


def render() -> None:
    st.title("Field Map")

    cfg = require_config()

    bias_V = st.number_input(
        "Reverse bias (V)",
        value=-50.0,
        min_value=-200.0,
        max_value=0.0,
        step=5.0,
        key="field_bias_V",
        help=(
            "Negative = reverse bias. 2D devices use the fast 2D engine for "
            "|V| between 5 and 80 V; other cases use devsim, which may not "
            "converge at deep reverse bias."
        ),
    )

    inputs = {"reverse bias": float(bias_V)}
    if st.button("Run simulation"):
        try:
            with st.spinner("Solving the field..."):
                st.session_state["field_result"] = etna.run_field(
                    cfg, bias_V=float(bias_V)
                )
            st.session_state["field_inputs"] = inputs
        except Exception as e:  # incl. devsim.error (not a RuntimeError)
            show_error(
                e,
                "Try a smaller reverse bias or adjust epi thickness / doping in "
                "the sidebar.",
            )

    result = st.session_state.get("field_result")
    if result is not None and result.mesh is not None:
        show_stale_warning(
            result.config, cfg, st.session_state.get("field_inputs"), inputs
        )
        st.caption(f"Bias: {result.metadata.get('bias_V', float('nan')):.1f} V")
        show_result_notes(result)
        if result.mesh.y_coords is None:
            st.caption(
                "Sharp nm-scale spikes at the abrupt doping steps (p+/n- junction, epi/n+ "
                "interface) depend on the mesh and on the field estimator (the 2D result "
                "shows none); they are not breakdown sites and are excluded from the "
                "reported peak field."
            )
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

        if result.mesh.y_coords is not None:
            # 2D: export the mesh node table (x lateral, y depth, fields).
            st.download_button(
                "Download CSV (2D mesh nodes)",
                data=to_csv_bytes(result),
                file_name="field_2d_result.csv",
                mime="text/csv",
            )
