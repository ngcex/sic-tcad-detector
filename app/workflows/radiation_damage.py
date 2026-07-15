"""Radiation Damage page: 1D guard, Run -> cache -> render -> download.

All Streamlit calls live inside render() -- no module-level side effects --
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests
(tests/test_app_radiation_damage_page.py). etna.run_radiation_damage is
referenced as a MODULE ATTRIBUTE (not `from etna import
run_radiation_damage`) so tests can intercept it via
monkeypatch.setattr(etna, "run_radiation_damage", fake) -- the seam
proven in tests/test_app_run_mockability.py (39-01).

The kappa (NIEL hardness factor) data-blocked warning banner is rendered
UNCONDITIONALLY, before the `cfg is None` empty-state guard, so it is
visible on every page load regardless of device-config or run state (per
project memory: stopping power CSV data produces unrealistically flat
kappa -- this must stay visible until real SR-NIEL SiC proton NIEL data
replaces the placeholder).

The `run_radiation_damage` facade can return a `result.y` array containing
a single `np.nan` mid-array even at the spike-confirmed safe V_bias=-20.0
default (see 41-01-SPIKE-NOTES.md Check B) -- this is expected, not a bug,
and is handled here via NaN-tolerant rendering (build_damage_figure does
not drop NaNs; Plotly renders them as a native line gap) plus an
informational message.
"""

from __future__ import annotations

import numpy as np
import streamlit as st

import etna
from app.components.results import build_damage_figure, to_csv_bytes


def render() -> None:
    st.title("Radiation Damage")

    st.warning(
        "**Data-blocked placeholder:** kappa (NIEL hardness factor) values "
        "used here are unvalidated placeholders. The energy TREND is "
        "physically motivated, but no ABSOLUTE Phi_crit or "
        "defect-concentration number is citable until real SR-NIEL SiC "
        "proton NIEL data replaces these placeholders. Treat the "
        "CCE-vs-fluence curve below as a relative sensitivity shape only, "
        "not an absolute prediction."
    )

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    if cfg.half_width_um is not None:
        st.warning(
            "These workflows are 1D-only. Set Dimensionality to 1D in the " "sidebar."
        )
        st.stop()

    col1, col2, col3 = st.columns(3)
    with col1:
        fluence_min = st.number_input(
            "Min fluence (p/cm²)",
            format="%.3e",
            value=1e13,
            key="rad_fluence_min",
        )
    with col2:
        fluence_max = st.number_input(
            "Max fluence (p/cm²)",
            format="%.3e",
            value=1e16,
            key="rad_fluence_max",
        )
    with col3:
        n_points = st.number_input(
            "Number of points",
            min_value=2,
            step=1,
            value=6,
            key="rad_n_points",
        )

    proton_energy_MeV = st.selectbox(
        "Proton energy (MeV)",
        [30, 62, 70, 150],
        index=1,
        key="rad_proton_energy",
    )
    V_bias = st.number_input(
        "Reverse bias V_bias (V)",
        value=-20.0,
        key="rad_v_bias",
    )

    if st.button("Run simulation"):
        try:
            fluences = np.geomspace(fluence_min, fluence_max, int(n_points))
            st.session_state["damage_result"] = etna.run_radiation_damage(
                cfg,
                fluences=fluences,
                V_bias=V_bias,
                proton_energy_MeV=proton_energy_MeV,
            )
        except RuntimeError as e:
            st.error(
                f"Simulation failed to converge: {e}\n\n"
                "Try adjusting device parameters (e.g. epi thickness, "
                "doping) in the sidebar and running again."
            )

    result = st.session_state.get("damage_result")
    if result is not None:
        if np.any(np.isnan(result.y)):
            st.info(
                "One or more fluence points did not converge and are "
                "shown as gaps in the curve below. This does not affect "
                "the other points."
            )
        st.plotly_chart(build_damage_figure(result))
        st.download_button(
            "Download CSV",
            data=to_csv_bytes(result),
            file_name="radiation_damage_result.csv",
            mime="text/csv",
        )
