"""st.navigation entry script: petringa — SiC TCAD Simulator.

Run with: streamlit run app/main.py

Wires the device-config sidebar (renders on every page, before pg.run(),
so device_config persists across navigation per UI-07) and registers
Home + all 7 workflow placeholder pages via callable-based st.Page.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import streamlit as st

from app.components.device_sidebar import render_device_sidebar
from app.workflows.home import render as render_home
from app.workflows.cv import render as render_cv
from app.workflows.cce import render as render_cce
from app.workflows.field_map import render as render_field_map
from app.workflows.radiation_damage import render as render_radiation_damage
from app.workflows.dark_current import render as render_dark_current
from app.workflows.microdosimetry import render as render_microdosimetry
from app.workflows.batch_sweep import render as render_batch_sweep

st.set_page_config(
    page_title="petringa — SiC TCAD Simulator",
    page_icon=":material/memory:",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Runs on EVERY page (called before pg.run()) — writes st.session_state["device_config"]
# so the config survives navigation (UI-02 "on any page", UI-07 persistence).
render_device_sidebar()

# Keep this list easy to append to — Phases 39-42 add behavior to these same pages.
# Explicit url_path= is required: every page module exposes a function named
# `render`, so Streamlit's filename/callable-name URL-pathname inference would
# otherwise collide across all 8 pages.
pages = [
    st.Page(
        render_home, title="Home", icon=":material/home:", url_path="home", default=True
    ),
    st.Page(
        render_cv, title="C-V Analysis", icon=":material/show_chart:", url_path="cv"
    ),
    st.Page(
        render_cce,
        title="Charge Collection (CCE)",
        icon=":material/bolt:",
        url_path="cce",
    ),
    st.Page(
        render_field_map, title="Field Map", icon=":material/map:", url_path="field-map"
    ),
    st.Page(
        render_radiation_damage,
        title="Radiation Damage",
        icon=":material/warning:",
        url_path="radiation-damage",
    ),
    st.Page(
        render_dark_current,
        title="Dark Current",
        icon=":material/dark_mode:",
        url_path="dark-current",
    ),
    st.Page(
        render_microdosimetry,
        title="Microdosimetry",
        icon=":material/scatter_plot:",
        url_path="microdosimetry",
    ),
    st.Page(
        render_batch_sweep,
        title="Batch Sweep",
        icon=":material/stacked_line_chart:",
        url_path="batch-sweep",
    ),
]

pg = st.navigation(pages)
pg.run()
