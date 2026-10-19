"""Home landing page: orientation line + current device config summary.

All Streamlit calls live inside render() — no module-level side effects —
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests.
"""

from __future__ import annotations

import streamlit as st

from app.components.page_utils import require_config


def render() -> None:
    st.title("ETNA — SiC TCAD Simulator")
    st.write(
        "Configure a device in the sidebar, then choose a simulation workflow "
        "from the navigation to explore its results."
    )
    st.caption(
        "C-V and CCE run on a fast batched engine inside the validated parameter "
        "domain (where it was cross-checked against devsim) and on devsim "
        "otherwise; 2D field maps use a 2D engine where eligible. Every result "
        "page shows the solver used."
    )

    cfg = require_config()

    st.json({k: getattr(cfg, k) for k in cfg.__dataclass_fields__})
