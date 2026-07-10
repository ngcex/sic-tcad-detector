"""Microdosimetry placeholder page (behavior arrives in Phase 42)."""

from __future__ import annotations

import streamlit as st


def render() -> None:
    st.title("Microdosimetry")

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    st.json({k: getattr(cfg, k) for k in cfg.__dataclass_fields__})
    st.caption("Running this simulation is implemented in Phase 42.")
