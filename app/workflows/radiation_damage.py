"""Radiation Damage page: 1D guard, Run -> cache -> render -> download.

All Streamlit calls live inside render() -- no module-level side effects --
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests
(tests/test_app_radiation_damage_page.py). etna.run_radiation_damage is
referenced as a MODULE ATTRIBUTE (not `from etna import
run_radiation_damage`) so tests can intercept it via
monkeypatch.setattr(etna, "run_radiation_damage", fake) -- the seam
proven in tests/test_app_run_mockability.py (39-01).

The kappa (NIEL hardness factor) caveat banner is rendered UNCONDITIONALLY,
before the `cfg is None` empty-state guard, so it is visible on every page
load regardless of device-config or run state. Since 2026-10-20 (audit B-07)
kappa is real SR-NIEL proton NIEL data (data/niel/) in the Si 1-MeV-neq
convention of the defect introduction rates; the SiC neutron normalisation
is still open, which is what the banner says.

The `run_radiation_damage` facade can return a `result.y` array containing
a single `np.nan` mid-array even at the spike-confirmed safe V_bias=-20.0
default (see 41-01-SPIKE-NOTES.md Check B) -- this is expected, not a bug,
and is handled here via NaN-tolerant rendering (build_damage_figure does
not drop NaNs; Plotly renders them as a native line gap) plus an
informational message.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import streamlit as st

import etna
from etna.api.simulation import compensation_fluence, trap_compensation_fluence
from etna.core.radiation_damage import get_hardness_factor
from app.components.page_utils import (
    about,
    require_config,
    show_error,
    show_stale_warning,
    solver_badge_text,
)
from app.components.results import build_damage_figure, to_csv_bytes


def render() -> None:
    st.title("Radiation Damage")
    st.caption(
        "CCE of the alpha test signal (as on the Charge Collection page) at one reverse "
        "bias vs the proton fluence the device has received."
    )
    about(
        "- **Damage scale:** 1 MeV-neutron-equivalent fluence Φ_neq = κ(E)·Φ_p, with "
        "κ(E) the SR-NIEL proton NIEL hardness factor (1-1000 MeV, Si 1-MeV-neq "
        "convention of the introduction rates).\n"
        "- **Default model (carrier removal):** each Φ_neq removes donors at a constant "
        "rate and adds SRH lifetime damage. When the removal reaches the lowest epi "
        "N_D (the *compensation fluence*) the epi is fully compensated: those points "
        "are **refused**, not solved (the model is ill-posed there).\n"
        "- **Opt-in charged deep traps:** Z1/2, EH4, EH6/7 (Burin Table I) as charged "
        "species with dark SRH occupancy; the compensated epi turns semi-insulating and "
        "is solved up to 1e16 neq. Batched engine only; exclusive with incomplete "
        "ionization.\n"
        "- **Validity:** introduction rates are literature values with a wide spread "
        "and the SiC neutron normalisation of κ is open, so absolute fluences "
        "(Φ_crit) are not validated hardness predictions; at low energies (~5 MeV) "
        "NIEL scaling overpredicts measured Si proton damage by tens of percent."
    )

    st.warning(
        "**NIEL normalisation open:** kappa (NIEL hardness factor) is "
        "SR-NIEL proton NIEL data (1-1000 MeV, data/niel/) in the Si "
        "1-MeV-neq convention of the defect introduction rates. The SiC "
        "neutron normalisation is still open (no SiC reactor-neutron NIEL "
        "available), so the absolute SiC fluence scale (Phi_crit, defect "
        "concentrations) is unvalidated; if the SiC/Si NIEL ratio of the "
        "reactor neutrons exceeds the proton one, kappa overstates the SiC "
        "damage per proton."
    )

    cfg = require_config()

    if cfg.half_width_um is not None:
        st.warning(
            "This workflow is 1D-only. Set Dimensionality to 1D in the sidebar."
        )
        st.stop()

    # Energy first: the compensation fluence (and thus the default range)
    # depends on it and on the sidebar doping.
    proton_energy_MeV = float(
        st.number_input(
            "Proton energy (MeV)",
            min_value=1.0,
            max_value=1000.0,
            value=62.0,
            step=1.0,
            format="%g",
            key="rad_proton_energy",
            help="Kinetic energy of the incident protons; sets κ(E) (SR-NIEL table, "
            "1-1000 MeV). 62 MeV: the CATANA / LNS proton line.",
        )
    )
    kappa = get_hardness_factor(proton_energy_MeV)
    phi_c = compensation_fluence(cfg, proton_energy_MeV)
    ii_on = cfg.incomplete_ionization not in (False, None)
    use_traps = st.checkbox(
        "Charged deep traps in Poisson (opt-in, solves beyond compensation)",
        value=False,
        key="rad_traps",
        disabled=ii_on,
        help=(
            "Replaces the donor-removal scalar by the Burin Table I deep levels "
            "(Z1/2, EH4, EH6/7) as charged species with dark steady-state SRH "
            "occupancy: the compensated epi turns semi-insulating and is solved. "
            "Batched engine only. See docs/USAGE.md for what is calibrated and what "
            "is assumed."
        ),
    )
    if ii_on:
        st.caption(
            "Charged deep traps are unavailable while 'Incomplete ionization' is on in "
            "the sidebar (the two models are exclusive)."
        )
        use_traps = False
    st.caption(
        f"κ({proton_energy_MeV:g} MeV) = {kappa:.3g} (1 MeV neq per proton), so "
        f"1e12 p/cm² ≙ {kappa * 1e12:.3g} neq/cm²."
    )
    run_cfg = cfg
    if use_traps:
        try:
            run_cfg = dataclasses.replace(cfg, radiation_traps=True)
        except ValueError as e:
            show_error(e, "Turn off incomplete ionization in the sidebar.")
            st.stop()
        top = 1e16 / kappa
        phi_t = trap_compensation_fluence(run_cfg, proton_energy_MeV)
        st.info(
            f"Charged-trap model: the epi is compensated at about **{phi_t:.2e} "
            f"p/cm²** (emergent acceptor removal; legacy carrier-removal reference "
            f"{phi_c:.2e} p/cm²); beyond it the epi is semi-insulating (Fermi level "
            "pinned by the deep levels) and still solved. The default range reaches "
            f"1e16 neq ({top:.2e} p/cm²)."
        )
        default_min, default_max = 1e-2 * phi_c, top
    else:
        st.info(
            f"Full donor compensation at about **{phi_c:.2e} p/cm²** for this device "
            f"at {proton_energy_MeV:g} MeV (carrier removal × κ(E) × Φ reaches the "
            "lowest epi N_D). Beyond it the epi is fully compensated and the model "
            "is ill-posed, so those points are not solved. The default range spans "
            "three decades up to 0.7× that fluence."
        )
        default_min, default_max = 0.7e-3 * phi_c, 0.7 * phi_c
    # Defaults are keyed on (energy, compensation fluence, model) so they refresh when
    # either changes; user edits persist while all stay the same.
    ctx = f"{proton_energy_MeV:g}_{phi_c:.4g}_{int(use_traps)}"
    col1, col2, col3 = st.columns(3)
    with col1:
        fluence_min = st.number_input(
            "Min fluence (p/cm²)",
            format="%.3e",
            value=float(f"{default_min:.3g}"),
            key=f"rad_fluence_min_{ctx}",
            help="Proton fluence Φ_p (not neq); the default follows the energy and the "
            "model.",
        )
    with col2:
        fluence_max = st.number_input(
            "Max fluence (p/cm²)",
            format="%.3e",
            value=float(f"{default_max:.3g}"),
            key=f"rad_fluence_max_{ctx}",
            help="Without charged deep traps keep it below the compensation fluence "
            "above: points beyond it are refused.",
        )
    with col3:
        n_points = st.number_input(
            "Number of points",
            min_value=2,
            max_value=50,
            step=1,
            value=6,
            key="rad_n_points",
            help="Log-spaced between the min and max fluence.",
        )

    V_bias = st.number_input(
        "Reverse bias V_bias (V)",
        value=-20.0,
        key="rad_v_bias",
        help="Bias at which the CCE is computed for every fluence (negative = reverse).",
    )

    inputs = {
        "min fluence": float(fluence_min),
        "max fluence": float(fluence_max),
        "number of points": int(n_points),
        "proton energy": proton_energy_MeV,
        "V_bias": float(V_bias),
    }
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
                        run_cfg,
                        fluences=fluences,
                        V_bias=V_bias,
                        proton_energy_MeV=proton_energy_MeV,
                    )
                st.session_state["damage_inputs"] = inputs
            except Exception as e:  # incl. devsim.error (not a RuntimeError)
                show_error(
                    e,
                    "Try a smaller reverse bias or adjust epi thickness / doping "
                    "in the sidebar.",
                )

    result = st.session_state.get("damage_result")
    if result is not None:
        show_stale_warning(
            result.config, run_cfg, st.session_state.get("damage_inputs"), inputs
        )
        badge = solver_badge_text(result.metadata)
        if badge:
            st.caption(badge)
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
            n_ex = len(excluded)
            st.warning(
                f"{n_ex} of {len(result.x)} fluence "
                + ("point is" if n_ex == 1 else "points are")
                + " not shown: "
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
                + (
                    ". Enable 'Charged deep traps' to solve beyond the compensation "
                    "fluence."
                    if any("compensation" in r for _, r in excluded) and not use_traps
                    else ""
                )
            )
        if not valid.any():
            st.error(
                "No valid fluence point: lower the fluence range below the "
                "compensation fluence shown above"
                + ("." if use_traps else ", or enable 'Charged deep traps'.")
            )
        st.plotly_chart(build_damage_figure(result))
        st.download_button(
            "Download CSV",
            data=to_csv_bytes(result),
            file_name="radiation_damage_result.csv",
            mime="text/csv",
        )
