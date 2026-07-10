"""Home landing page: orientation line + current device config summary.

All Streamlit calls live inside render() — no module-level side effects —
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests.
"""

from __future__ import annotations

import streamlit as st


def render() -> None:
    st.title("petringa — SiC TCAD Simulator")
    st.write(
        "Configure a device in the sidebar, then choose a simulation workflow "
        "from the navigation to explore its results."
    )

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    st.json({k: getattr(cfg, k) for k in cfg.__dataclass_fields__})
