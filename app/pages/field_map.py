"""Field Map placeholder page (behavior arrives in Phase 39; geometry viewer
component arrives in Phase 40 as part of this page, not a separate nav entry).
"""

from __future__ import annotations

import streamlit as st


def render() -> None:
    st.title("Field Map")

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    st.json({k: getattr(cfg, k) for k in cfg.__dataclass_fields__})
    st.caption("Running this simulation is implemented in Phase 39.")
