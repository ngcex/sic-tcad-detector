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
kappa is real SR-NIEL proton NIEL data (data/niel/); since 2026-10-21 the SiC
neutron NIEL is computed (NJOY, data/niel/), the DEFAULT convention is
"sic-damage" (si-neq / sic-1mev selectable) and the DEFAULT damage model is the
calibrated kinetic trap preset (legacy carrier removal and the Burin steady traps
selectable); the banner states what is calibrated and what stays open.

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
from etna.api.simulation import (
    compensation_fluence,
    damage_model_of,
    legacy_compensation_fluence,
)
from etna.core.radiation_damage import (
    DEFAULT_NIEL_CONVENTION,
    NIEL_CONVENTIONS,
    get_hardness_factor,
    hardness_table,
)

# Damage-model choices of the page -> DeviceConfig.radiation_traps
DAMAGE_MODELS = {
    "Calibrated kinetic traps (default)": "auto",
    "Legacy carrier removal": "legacy",
    "Burin Table I traps (steady state)": True,
}
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
        "bias vs the particle (proton or neutron) fluence the device has received."
    )
    about(
        "- **Damage scale:** 1 MeV-neutron-equivalent fluence Φ_neq = κ(E)·Φ_p, with "
        "κ(E) the SR-NIEL proton NIEL hardness factor (1-1000 MeV). Default "
        "convention **sic-damage**: the reactor-neq fluence that does the same SiC "
        "damage (the Burin introduction rates are neutron-anchored); 'si-neq' (the "
        "previous default) overstates the SiC damage per proton by r_n/r_p.\n"
        "- **Default model (calibrated kinetic traps):** Z1/Z2 (negative-U, Capan), "
        "EH4, EH6/7 as charged species with kinetic occupancy (0 V start, ramp, "
        "hold) and occupancy-weighted lifetimes, plus a direct donor removal "
        "**calibrated** to published 252.7 MeV C-V removal slopes (per convention). "
        "The compensated epi turns semi-insulating and is solved up to 1e16 neq. "
        "Batched engine only: with incomplete ionization (or outside the engine's "
        "box) it runs the legacy model and says so.\n"
        "- **Legacy carrier removal:** each Φ_neq removes donors at 5 /neq; points "
        "beyond the *compensation fluence* are **refused** (ill-posed).\n"
        "- **Validity:** the removal is calibrated on other devices' C-V at one "
        "temperature (the mechanism, direct removal vs a slow deep acceptor, is not "
        "identified); CCE beyond compensation is a prediction, unvalidated. At low "
        "energies (~5 MeV) NIEL scaling overpredicts measured Si proton damage by "
        "tens of percent."
    )

    st.warning(
        "**NIEL normalisation:** kappa (NIEL hardness factor) is SR-NIEL proton "
        "NIEL data (1-1000 MeV, data/niel/), by default in the 'sic-damage' "
        "convention: SiC NIEL normalised by the computed SiC/Si fission-neutron "
        "ratio r_n = 1.069 (NJOY/ENDF-B-VIII.0), so a proton fluence maps onto the "
        "SiC damage of the reactor-neutron fluence the introduction rates are "
        "labelled in ('si-neq', the previous default, gives +25 % at 62 MeV). The "
        "literature introduction rates stay unvalidated and the calibrated removal "
        "comes from other devices, so absolute fluences (Phi_crit) are not "
        "validated hardness predictions."
    )

    cfg = require_config()

    if cfg.half_width_um is not None:
        st.warning(
            "This workflow is 1D-only. Set Dimensionality to 1D in the sidebar."
        )
        st.stop()

    # Particle / NIEL convention (defaults = protons in sic-damage), then the
    # energy: the compensation fluence (and thus the default range) depends on them.
    pcol, ccol = st.columns(2)
    with pcol:
        particle = st.radio(
            "Particle",
            ["proton", "neutron"],
            horizontal=True,
            key="rad_particle",
            help="Neutrons: NJOY HEATR / ENDF-B-VIII.0 NIEL (data/niel/), "
            "lethargy-averaged; energies 1e-3 to 18.8 MeV on this page.",
        )
    with ccol:
        niel_convention = st.selectbox(
            "NIEL convention",
            list(NIEL_CONVENTIONS),
            index=list(NIEL_CONVENTIONS).index(DEFAULT_NIEL_CONVENTION),
            key="rad_niel_convention",
            help="sic-damage (default): the Si-labelled reactor-neq fluence that gives "
            "the same SiC damage (SiC NIEL divided by the SiC/Si fission ratio r_n). "
            "si-neq (the previous default): silicon NIEL / ASTM 95 MeV mb. sic-1mev: "
            "SiC-native 1 MeV neq. The calibrated removal is re-fitted per convention, "
            "so its per-proton calibration holds in each. See docs/USAGE.md.",
        )
    niel = dict(particle=particle, niel_convention=niel_convention)
    # only non-default choices are passed on (default calls stay minimal)
    default_niel = {"particle": "proton", "niel_convention": DEFAULT_NIEL_CONVENTION}
    niel_kw = {} if niel == default_niel else niel
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
    ii_on = cfg.incomplete_ionization not in (False, None)
    model_name = st.selectbox(
        "Damage model",
        list(DAMAGE_MODELS),
        index=0,
        key="rad_damage_model",
        help=(
            "Calibrated kinetic traps (default): Z1/Z2, EH4, EH6/7 as charged species "
            "with kinetic occupancy plus a direct donor removal CALIBRATED to "
            "published 252.7 MeV C-V slopes; solves beyond compensation. Legacy "
            "carrier removal: 5 /neq donor removal, refuses points beyond "
            "compensation (the default until 2026-10-21). Burin Table I traps: the "
            "three literature levels, steady-state occupancy (opt-in since "
            "2026-10-18). See docs/USAGE.md for what is calibrated and what is "
            "predicted."
        ),
    )
    spec = DAMAGE_MODELS[model_name]
    if ii_on and spec is True:
        st.caption(
            "The Burin trap model is unavailable while 'Incomplete ionization' is on "
            "in the sidebar (the two models are exclusive): running the default."
        )
        spec = "auto"
    try:
        run_cfg = (
            cfg
            if cfg.radiation_traps == spec
            else dataclasses.replace(cfg, radiation_traps=spec)
        )
    except ValueError as e:
        show_error(e, "Turn off incomplete ionization in the sidebar.")
        st.stop()
    model, note = damage_model_of(run_cfg, proton_energy_MeV, **niel_kw)
    use_traps = model is not None
    phi_c = compensation_fluence(run_cfg, proton_energy_MeV, **niel_kw)
    phi_legacy = legacy_compensation_fluence(cfg, proton_energy_MeV, **niel_kw)
    st.caption(
        f"κ({proton_energy_MeV:g} MeV) = {kappa:.3g} (1 MeV neq per {pname}, "
        f"{niel_convention}), so 1e12 {pu} ≙ {kappa * 1e12:.3g} neq/cm²."
    )
    if note:
        st.caption(note)
    if use_traps:
        top = 1e16 / kappa
        st.info(
            f"Charged-trap model ({model.label.split(':')[0] or 'custom'}): the epi is "
            f"compensated at about **{phi_c:.2e} {pu}** (neutral-bulk removal "
            f"{model.neutral_removal_rate:.3g} /neq; legacy carrier-removal "
            f"reference {phi_legacy:.2e} {pu}); beyond it the epi is "
            "semi-insulating (Fermi level pinned by the deep levels) and still "
            f"solved. The default range reaches 1e16 neq ({top:.2e} {pu})."
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
        # always recorded: a proton -> neutron switch must flag the result stale
        **niel,
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
                    ". Choose a charged-trap damage model to solve beyond the "
                    "compensation fluence."
                    if any("compensation" in r for _, r in excluded) and not use_traps
                    else ""
                )
            )
        if not valid.any():
            st.error(
                "No valid fluence point: lower the fluence range below the "
                "compensation fluence shown above"
                + ("." if use_traps else ", or choose a charged-trap damage model.")
            )
        st.plotly_chart(build_damage_figure(result))
        st.download_button(
            "Download CSV",
            data=to_csv_bytes(result),
            file_name="radiation_damage_result.csv",
            mime="text/csv",
        )
