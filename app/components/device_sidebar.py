"""Device configuration sidebar: pure config-assembly seam + Streamlit renderer.

`assemble_config(values: dict) -> DeviceConfig` is a pure function with ZERO
Streamlit dependency — it is imported directly by unit tests
(tests/test_app_device_sidebar.py, tests/test_app_session.py) so the
config-assembly and mode-mapping logic (UI-02, UI-07) is testable without a
running Streamlit server.

`render_device_sidebar()` is the only function in this module that touches
`st.*` — it renders the reactive (non-form) sidebar widgets (mode selectors
render outside any submit-gated wrapper so conditional fields toggle
immediately), collects their values into a dict, calls `assemble_config()`,
and stores the result under the single non-widget key
`st.session_state["device_config"]`.
"""

from __future__ import annotations

import streamlit as st

from petringa import DeviceConfig


def assemble_config(values: dict) -> DeviceConfig:
    """Turn sidebar form values into a DeviceConfig.

    Pure function, no Streamlit dependency. Constructs a DeviceConfig from
    all 11 supplied fields and applies the single mode-consistency rule that
    belongs to config assembly (not to the renderer): doping_profile
    "graded" forces N_D=None, "uniform" keeps N_D as supplied. All other
    fields (including the graded triplet and half_width_um) pass through
    unchanged — resetting hidden/gated fields to their defaults is
    render_device_sidebar's responsibility, not this pure seam's.
    """
    doping_profile = values["doping_profile"]
    N_D = None if doping_profile == "graded" else values["N_D"]

    return DeviceConfig(
        epi_thickness_um=values["epi_thickness_um"],
        substrate_thickness_um=values["substrate_thickness_um"],
        half_width_um=values["half_width_um"],
        N_A=values["N_A"],
        doping_profile=doping_profile,
        N_D=N_D,
        N_D_junction=values["N_D_junction"],
        N_D_bulk=values["N_D_bulk"],
        L_transition_um=values["L_transition_um"],
        T=values["T"],
        area_cm2=values["area_cm2"],
    )


def render_device_sidebar() -> None:
    """Render the reactive device-config sidebar and store the assembled
    DeviceConfig under st.session_state["device_config"].

    Not wrapped in a submit-gated container — the dimensionality and
    doping-profile mode selectors must be reactive so conditional fields
    appear/disappear immediately (38-UI-SPEC Group 1). Gated fields are set
    programmatically to DeviceConfig defaults when hidden, so the dict
    passed to assemble_config always has all 11 keys populated.
    """
    st.sidebar.header("Device configuration")

    # --- Group 1: reactive mode selectors (rendered plainly, no submit gate) ---
    dim = st.sidebar.radio("Dimensionality", ["1D", "2D"], index=0)
    profile = st.sidebar.selectbox("Doping profile", ["graded", "uniform"], index=0)

    st.sidebar.divider()

    # --- Group 2: geometry ---
    st.sidebar.header("Geometry")
    epi_thickness_um = st.sidebar.number_input(
        "Epi thickness (µm)", value=10.0, min_value=0.1
    )
    substrate_thickness_um = st.sidebar.number_input(
        "Substrate thickness (µm)", value=1.0, min_value=0.0
    )
    if dim == "2D":
        half_width_um = st.sidebar.number_input(
            "Half-width (µm)", value=50.0, min_value=0.1
        )
    else:
        half_width_um = None

    st.sidebar.divider()

    # --- Group 3: doping ---
    st.sidebar.header("Doping")
    N_A = st.sidebar.number_input("N_A substrate (cm⁻³)", value=1e19, format="%.3e")
    if profile == "uniform":
        N_D = st.sidebar.number_input("N_D uniform (cm⁻³)", value=1e15, format="%.3e")
        # Graded triplet hidden this mode — fall back to DeviceConfig defaults.
        N_D_junction = DeviceConfig.N_D_junction
        N_D_bulk = DeviceConfig.N_D_bulk
        L_transition_um = DeviceConfig.L_transition_um
    else:  # graded
        N_D_junction = st.sidebar.number_input(
            "N_D junction (cm⁻³)", value=2.93e15, format="%.3e"
        )
        N_D_bulk = st.sidebar.number_input(
            "N_D bulk (cm⁻³)", value=8.82e13, format="%.3e"
        )
        L_transition_um = st.sidebar.number_input(
            "Transition length (µm)", value=0.987, min_value=0.0
        )
        # N_D hidden this mode — programmatically set to default; assemble_config
        # also forces N_D=None for "graded", this default is a harmless input.
        N_D = DeviceConfig.N_D

    st.sidebar.divider()

    # --- Group 4: operating conditions ---
    st.sidebar.header("Operating conditions")
    T = st.sidebar.number_input("Temperature (K)", value=300.0, min_value=1.0)
    area_cm2 = st.sidebar.number_input("Area (cm²)", value=1e-4, format="%.3e")

    values = {
        "epi_thickness_um": epi_thickness_um,
        "substrate_thickness_um": substrate_thickness_um,
        "half_width_um": half_width_um,
        "N_A": N_A,
        "doping_profile": profile,
        "N_D": N_D,
        "N_D_junction": N_D_junction,
        "N_D_bulk": N_D_bulk,
        "L_transition_um": L_transition_um,
        "T": T,
        "area_cm2": area_cm2,
    }

    st.session_state["device_config"] = assemble_config(values)
