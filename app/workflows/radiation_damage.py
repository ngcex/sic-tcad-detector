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
from etna.api.simulation import compensation_fluence
from app.components.page_utils import show_error, show_stale_warning
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
            "This workflow is 1D-only. Set Dimensionality to 1D in the sidebar."
        )
        st.stop()

    # Energy first: the compensation fluence (and thus the default range)
    # depends on it and on the sidebar doping.
    proton_energy_MeV = st.selectbox(
        "Proton energy (MeV)",
        [30, 62, 70, 150],
        index=1,
        key="rad_proton_energy",
    )
    phi_c = compensation_fluence(cfg, proton_energy_MeV)
    st.info(
        f"Full donor compensation at about **{phi_c:.2e} p/cm²** for this device "
        f"at {proton_energy_MeV} MeV (carrier removal × κ(E) × Φ reaches the "
        "lowest epi N_D). Beyond it the epi is fully compensated and the model "
        "is ill-posed, so those points are not solved. The default range spans "
        "three decades up to 0.7× that fluence."
    )
    # Defaults are keyed on (energy, compensation fluence) so they refresh when
    # either changes; user edits persist while both stay the same.
    ctx = f"{proton_energy_MeV}_{phi_c:.4g}"
    col1, col2, col3 = st.columns(3)
    with col1:
        fluence_min = st.number_input(
            "Min fluence (p/cm²)",
            format="%.3e",
            value=float(f"{0.7e-3 * phi_c:.3g}"),
            key=f"rad_fluence_min_{ctx}",
        )
    with col2:
        fluence_max = st.number_input(
            "Max fluence (p/cm²)",
            format="%.3e",
            value=float(f"{0.7 * phi_c:.3g}"),
            key=f"rad_fluence_max_{ctx}",
        )
    with col3:
        n_points = st.number_input(
            "Number of points",
            min_value=2,
            max_value=50,
            step=1,
            value=6,
            key="rad_n_points",
        )

    V_bias = st.number_input(
        "Reverse bias V_bias (V)",
        value=-20.0,
        key="rad_v_bias",
    )

    if st.button("Run simulation"):
        if not (0 < fluence_min <= fluence_max):
            st.error("Fluences must be positive, with min ≤ max.")
        elif V_bias > 0:
            st.error("V_bias must be ≤ 0 V (negative = reverse bias).")
        else:
            try:
                with st.spinner("Running fluence sweep..."):
                    fluences = np.geomspace(fluence_min, fluence_max, int(n_points))
                    st.session_state["damage_result"] = etna.run_radiation_damage(
                        cfg,
                        fluences=fluences,
                        V_bias=V_bias,
                        proton_energy_MeV=proton_energy_MeV,
                    )
            except Exception as e:  # incl. devsim.error (not a RuntimeError)
                show_error(
                    e,
                    "Try a smaller reverse bias or adjust epi thickness / doping "
                    "in the sidebar.",
                )

    result = st.session_state.get("damage_result")
    if result is not None:
        show_stale_warning(result.config, cfg)
        if result.metadata.get("solver"):
            st.caption(f"Solver: {result.metadata['solver']}")
        valid = np.asarray(
            result.metadata.get("valid", np.isfinite(result.y)), dtype=bool
        )
        reasons = result.metadata.get("reason") or [""] * len(result.x)
        excluded = [
            (fl, reason or "not converged")
            for fl, ok, reason in zip(result.x, valid, reasons)
            if not ok
        ]
        if excluded:
            st.warning(
                f"{len(excluded)} of {len(result.x)} fluence points are not shown: "
                + "; ".join(
                    f"{fl:.2e} p/cm² ("
                    + (
                        "epi fully compensated: the model solves an ill-posed "
                        "problem there"
                        if "compensation" in reason
                        else reason
                    )
                    + ")"
                    for fl, reason in excluded
                )
            )
        if not valid.any():
            st.error(
                "No valid fluence point: lower the fluence range below the "
                "compensation fluence shown above."
            )
        st.plotly_chart(build_damage_figure(result))
        st.download_button(
            "Download CSV",
            data=to_csv_bytes(result),
            file_name="radiation_damage_result.csv",
            mime="text/csv",
        )
