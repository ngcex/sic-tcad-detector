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
convention of the defect introduction rates; since 2026-10-21 the SiC neutron
NIEL is computed (NJOY, data/niel/) and the page offers neutrons and the opt-in
sic-damage / sic-1mev conventions; the banner states r_n and what stays open.

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
from etna.core.radiation_damage import (
    NIEL_CONVENTIONS,
    get_hardness_factor,
    hardness_table,
)
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
        "and the default Si-neq κ overstates the SiC damage by r_n/r_p (opt-in "
        "'sic-damage' convention), so absolute fluences "
        "(Φ_crit) are not validated hardness predictions; at low energies (~5 MeV) "
        "NIEL scaling overpredicts measured Si proton damage by tens of percent."
    )

    st.warning(
        "**NIEL normalisation:** kappa (NIEL hardness factor) is SR-NIEL proton "
        "NIEL data (1-1000 MeV, data/niel/); by default in the Si 1-MeV-neq "
        "convention of the defect introduction rates. The SiC neutron NIEL is now "
        "computed (NJOY/ENDF-B-VIII.0: SiC/Si fission-spectrum ratio r_n = 1.069, "
        "> the proton r_p), so the default Si-neq kappa overstates the SiC damage "
        "per proton by r_n/r_p (+25 % at 62 MeV); the opt-in 'sic-damage' "
        "convention removes that. The introduction rates themselves stay "
        "unvalidated literature values, so absolute fluences (Phi_crit) are not "
        "validated hardness predictions."
    )

    cfg = require_config()

    if cfg.half_width_um is not None:
        st.warning(
            "This workflow is 1D-only. Set Dimensionality to 1D in the sidebar."
        )
        st.stop()

    # Particle / NIEL convention (opt-in, defaults = protons in Si-neq), then the
    # energy: the compensation fluence (and thus the default range) depends on them.
    pcol, ccol = st.columns(2)
    with pcol:
        particle = st.radio(
            "Particle",
            ["proton", "neutron"],
            horizontal=True,
            key="rad_particle",
            help="Neutrons: NJOY HEATR / ENDF-B-VIII.0 NIEL (data/niel/), "
            "lethargy-averaged, 1e-5 to 18.9 MeV.",
        )
    with ccol:
        niel_convention = st.selectbox(
            "NIEL convention",
            list(NIEL_CONVENTIONS),
            index=0,
            key="rad_niel_convention",
            help="si-neq (default): silicon NIEL / ASTM 95 MeV mb, the label of the "
            "introduction rates. sic-damage: the Si-labelled reactor-neq fluence that "
            "gives the same SiC damage (divides by the SiC/Si fission ratio r_n). "
            "sic-1mev: SiC-native 1 MeV neq. See data/niel/README.md.",
        )
    niel = dict(particle=particle, niel_convention=niel_convention)
    # only non-default choices are passed on (default calls stay byte-identical)
    niel_kw = {} if niel == {"particle": "proton", "niel_convention": "si-neq"} else niel
    if particle == "proton":
        proton_energy_MeV = float(
            st.number_input(
                "Proton energy (MeV)",
                min_value=1.0,
                max_value=1000.0,
                value=62.0,
                step=1.0,
                format="%g",
                key="rad_proton_energy",
                help="Kinetic energy of the incident protons; sets κ(E) (SR-NIEL "
                "table, 1-1000 MeV). 62 MeV: the CATANA / LNS proton line.",
            )
        )
    else:
        proton_energy_MeV = float(
            st.number_input(
                "Neutron energy (MeV)",
                min_value=0.001,
                max_value=18.8,
                value=1.0,
                format="%g",
                key="rad_neutron_energy",
                help="Kinetic energy of the incident neutrons; sets κ(E) from the "
                "bin-averaged NJOY table (resonances averaged over 1/20 decade).",
            )
        )
    pu = "p/cm²" if particle == "proton" else "n/cm²"
    pname = "proton" if particle == "proton" else "neutron"
    kappa = get_hardness_factor(
        proton_energy_MeV, hardness_table(particle, niel_convention)
    )
    phi_c = compensation_fluence(cfg, proton_energy_MeV, **niel_kw)
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
        f"κ({proton_energy_MeV:g} MeV) = {kappa:.3g} (1 MeV neq per {pname}, "
        f"{niel_convention}), so 1e12 {pu} ≙ {kappa * 1e12:.3g} neq/cm²."
    )
    run_cfg = cfg
    if use_traps:
        try:
            run_cfg = dataclasses.replace(cfg, radiation_traps=True)
        except ValueError as e:
            show_error(e, "Turn off incomplete ionization in the sidebar.")
            st.stop()
        top = 1e16 / kappa
        phi_t = trap_compensation_fluence(run_cfg, proton_energy_MeV, **niel_kw)
        st.info(
            f"Charged-trap model: the epi is compensated at about **{phi_t:.2e} "
            f"{pu}** (emergent acceptor removal; legacy carrier-removal reference "
            f"{phi_c:.2e} {pu}); beyond it the epi is semi-insulating (Fermi level "
            "pinned by the deep levels) and still solved. The default range reaches "
            f"1e16 neq ({top:.2e} {pu})."
        )
        default_min, default_max = 1e-2 * phi_c, top
    else:
        st.info(
            f"Full donor compensation at about **{phi_c:.2e} {pu}** for this device "
            f"at {proton_energy_MeV:g} MeV (carrier removal × κ(E) × Φ reaches the "
            "lowest epi N_D). Beyond it the epi is fully compensated and the model "
            "is ill-posed, so those points are not solved. The default range spans "
            "three decades up to 0.7× that fluence."
        )
        default_min, default_max = 0.7e-3 * phi_c, 0.7 * phi_c
    # Defaults are keyed on (energy, compensation fluence, model) so they refresh when
    # either changes; user edits persist while all stay the same.
    ctx = f"{proton_energy_MeV:g}_{phi_c:.4g}_{int(use_traps)}"
    if niel_kw:
        ctx += f"_{particle}_{niel_convention}"  # default keys unchanged
    col1, col2, col3 = st.columns(3)
    with col1:
        fluence_min = st.number_input(
            f"Min fluence ({pu})",
            format="%.3e",
            value=float(f"{default_min:.3g}"),
            key=f"rad_fluence_min_{ctx}",
            help="Particle fluence Φ (not neq); the default follows the energy and the "
            "model.",
        )
    with col2:
        fluence_max = st.number_input(
            f"Max fluence ({pu})",
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
        **niel_kw,
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
                        **niel_kw,
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
                    f"{fl:.2e} {pu} ("
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
