"""Charge Collection (CCE) page: 1D guard, Run -> cache -> render -> download.

All Streamlit calls live inside render() -- no module-level side effects --
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests
(tests/test_app_cce_page.py). petringa.run_cce is referenced as a MODULE
ATTRIBUTE (not `from petringa import run_cce`) so tests can intercept it via
monkeypatch.setattr(petringa, "run_cce", fake) -- the seam proven in
tests/test_app_run_mockability.py (39-01).
"""

from __future__ import annotations

import streamlit as st

import petringa
from app.components.results import build_cce_figure, to_csv_bytes


def render() -> None:
    st.title("Charge Collection (CCE)")

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
            st.session_state["cce_result"] = petringa.run_cce(cfg)
        except RuntimeError as e:
            st.error(
                f"Simulation failed to converge: {e}\n\n"
                "Try adjusting device parameters (e.g. epi thickness, doping) "
                "in the sidebar and running again."
            )

    result = st.session_state.get("cce_result")
    if result is not None:
        st.plotly_chart(build_cce_figure(result))
        st.download_button(
            "Download CSV",
            data=to_csv_bytes(result),
            file_name="cce_result.csv",
            mime="text/csv",
        )

        st.caption(
            f"I_collected: {result.metadata['I_collected']} A/cm² | "
            f"I_generated: {result.metadata['I_generated']} A/cm²"
        )
