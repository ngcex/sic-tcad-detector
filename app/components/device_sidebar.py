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

from etna import DeviceConfig


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
        # n+ layer and dark-current fit keys are optional (older callers/tests pass the
        # original 11 fields; the defaults are the default stack)
        nplus_thickness_um=values.get(
            "nplus_thickness_um", DeviceConfig.nplus_thickness_um
        ),
        N_D_nplus=values.get("N_D_nplus", DeviceConfig.N_D_nplus),
        N_t=values.get("N_t", DeviceConfig.N_t),
        E_t=values.get("E_t", DeviceConfig.E_t),
        # graded n-/n+ transition: optional keys, default = the default stack's (both
        # None = abrupt epi end, e.g. the legacy / previous calibrated presets)
        buffer_thickness_um=values.get(
            "buffer_thickness_um", DeviceConfig.buffer_thickness_um
        ),
        buffer_end_doping=values.get("buffer_end_doping", DeviceConfig.buffer_end_doping),
    )


# Sidebar labels for DeviceConfig fields, used to phrase validation errors.
FIELD_LABELS = {
    "epi_thickness_um": "Epi thickness (µm)",
    "substrate_thickness_um": "Substrate thickness (µm)",
    "half_width_um": "Half-width (µm)",
    "N_A": "N_A substrate (cm⁻³)",
    "N_D_junction": "N_D junction (cm⁻³)",
    "N_D_bulk": "N_D bulk (cm⁻³)",
    "N_D": "N_D uniform (cm⁻³)",
    "L_transition_um": "Transition length (µm)",
    "T": "Temperature (K)",
    "area_cm2": "Area (cm²)",
    "nplus_thickness_um": "n+ layer thickness (µm)",
    "N_D_nplus": "N_D n+ (cm⁻³)",
    "N_t": "N_t dark-current fit (cm⁻³ s⁻¹)",
    "E_t": "E_t trap level (eV)",
    "buffer_thickness_um": "n-/n+ transition thickness (µm)",
    "buffer_end_doping": "N_D at the epi end (cm⁻³)",
}


def humanize_config_error(msg: str) -> str:
    """Replace DeviceConfig field names in a validation message by sidebar labels."""
    import re

    names = sorted(FIELD_LABELS, key=len, reverse=True)
    pattern = re.compile(r"(?<![\w])(" + "|".join(map(re.escape, names)) + r")(?![\w])")
    return pattern.sub(lambda m: f"'{FIELD_LABELS[m.group(1)]}'", msg)


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
    dim = st.sidebar.radio(
        "Dimensionality", ["1D", "2D"], index=0, key="cfg_dimensionality"
    )
    profile = st.sidebar.selectbox(
        "Doping profile", ["graded", "uniform"], index=0, key="cfg_doping_profile"
    )

    st.sidebar.divider()

    # --- Group 2: geometry ---
    st.sidebar.header("Geometry")
    # No widget min/max bounds: an out-of-range entry would be silently clamped
    # by the browser while the field still shows it. DeviceConfig validation
    # names the rejected value instead (st.sidebar.error below) and the pages
    # refuse to run until it is fixed.
    epi_thickness_um = st.sidebar.number_input(
        "Epi thickness (µm)",
        value=10.0,
        step=0.5,
        format="%.2f",
        key="cfg_epi_thickness_um",
    )
    substrate_thickness_um = st.sidebar.number_input(
        "Substrate thickness (µm)",
        value=DeviceConfig.substrate_thickness_um,
        step=0.1,
        format="%.2f",
        key="cfg_substrate_thickness_um",
        help="p+ anode layer (not the n+ substrate, see 'n+ layer under the epi').",
    )
    legacy_on = st.sidebar.checkbox(
        "historical model (cathode directly on n- epi)",
        value=False,
        key="cfg_legacy_on",
        help="Applies the legacy preset (the default before the n+ stack): no n+ layer, "
        "cathode contact directly on the n- epi, p+ 1 um, N_D_junction 2.93e15, N_D_bulk "
        "8.82e13, L 0.987 um, N_t 2.8e12, E_t 0.65 (overrides the fields above). It "
        "reproduces the old numbers; its non-convergence wall at deep reverse bias is a "
        "contact artefact. Off = the default stack: p+ / graded n- epi with a graded n-/n+ "
        "transition / 1 um n+ (1e19).",
    )
    nplus_cal_on = st.sidebar.checkbox(
        "previous calibrated model (n+, no buffer)",
        value=False,
        key="cfg_nplus_cal_on",
        help="Applies the previous default (etna.nplus_calibrated_config, 2026-10-17): "
        "exponential epi calibrated on two C-V points (N_D_junction 2.958237e15, N_D_bulk "
        "8.252732e13, L 0.987 um held), ABRUPT 1e19 n+ at the epi end (no graded "
        "transition), N_t 2.805336e12 (overrides the doping fields). Its C-V knee is at "
        "-8..-10 V where the measurement has it at -1..-4 V. The default (unchecked) is the "
        "measured-C-V profile: a fit to the digitised Photons Fig. 4b (plateau rms 0.013 um, "
        "knee rms 0.24 um; N_D_bulk 3.4e13 below the documented 0.5-1e14; epi/n+ position "
        "set by the fit, not by SIMS).",
    )
    nplus_thickness_um = DeviceConfig.nplus_thickness_um
    N_D_nplus = DeviceConfig.N_D_nplus
    if dim == "2D":
        half_width_um = st.sidebar.number_input(
            "Half-width (µm)",
            value=50.0,
            step=5.0,
            format="%.1f",
            key="cfg_half_width_um",
        )
    else:
        half_width_um = None

    st.sidebar.divider()

    # --- Group 3: doping ---
    st.sidebar.header("Doping")
    N_A = st.sidebar.number_input(
        "N_A substrate (cm⁻³)",
        value=1e19,
        step=1e18,
        format="%.3e",
        key="cfg_N_A",
    )
    if profile == "uniform":
        N_D = st.sidebar.number_input(
            "N_D uniform (cm⁻³)",
            value=1e15,
            step=1e14,
            format="%.3e",
            key="cfg_N_D",
        )
        # Graded triplet hidden this mode — fall back to DeviceConfig defaults.
        N_D_junction = DeviceConfig.N_D_junction
        N_D_bulk = DeviceConfig.N_D_bulk
        L_transition_um = DeviceConfig.L_transition_um
    else:  # graded
        N_D_junction = st.sidebar.number_input(
            "N_D junction (cm⁻³)",
            value=DeviceConfig.N_D_junction,
            step=1e9,  # the 7-digit exact-solve value must lie on the step grid
            format="%.3e",
            key="cfg_N_D_junction",
        )
        N_D_bulk = st.sidebar.number_input(
            "N_D bulk (cm⁻³)",
            value=DeviceConfig.N_D_bulk,
            step=1e7,  # idem
            format="%.3e",
            key="cfg_N_D_bulk",
        )
        L_transition_um = st.sidebar.number_input(
            "Transition length (µm)",
            value=DeviceConfig.L_transition_um,
            step=1e-7,  # idem (0.4638630)
            format="%.7f",
            key="cfg_L_transition_um",
        )
        # N_D hidden this mode — programmatically set to default; assemble_config
        # also forces N_D=None for "graded", this default is a harmless input.
        N_D = DeviceConfig.N_D
    # graded n-/n+ transition before the epi end (the default stack's C-V fit); off =
    # abrupt epi/n+ step
    buffer_on = st.sidebar.checkbox(
        "graded n-/n+ transition (C-V fit)",
        value=True,
        key="cfg_buffer_on",
        help="N_D rises log-linearly from the epi value to 'N_D at the epi end' over the "
        "last part of the epi (1e19 n+ beyond): the form fitted on the digitised measured "
        "C-V; not a measured profile (SIMS would settle it). Off = abrupt epi/n+ step.",
    )
    if buffer_on:
        buffer_thickness_um = st.sidebar.number_input(
            "n-/n+ transition thickness (µm)",
            value=DeviceConfig.buffer_thickness_um,
            step=1e-6,
            format="%.6f",
            key="cfg_buffer_thickness_um",
        )
        buffer_end_doping = st.sidebar.number_input(
            "N_D at the epi end (cm⁻³)",
            value=DeviceConfig.buffer_end_doping,
            step=1e9,
            format="%.3e",
            key="cfg_buffer_end_doping",
        )
    else:
        buffer_thickness_um = buffer_end_doping = None

    st.sidebar.divider()

    # --- Group 4: operating conditions ---
    st.sidebar.header("Operating conditions")
    T = st.sidebar.number_input(
        "Temperature (K)",
        value=300.0,
        step=5.0,
        format="%.1f",
        key="cfg_T",
    )
    area_cm2 = st.sidebar.number_input(
        "Area (cm²)",
        value=DeviceConfig.area_cm2,  # the 4 mm² dosimeter (0.04 cm²)
        step=1e-3,
        format="%.3e",
        key="cfg_area_cm2",
        help="Default: the 4 mm² dosimeter (4.000e-02 cm²). Capacitance and "
        "leakage scale with area.",
    )

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
        "nplus_thickness_um": nplus_thickness_um,
        "N_D_nplus": N_D_nplus,
        "buffer_thickness_um": buffer_thickness_um,
        "buffer_end_doping": buffer_end_doping,
    }

    if legacy_on:
        from etna import LEGACY_PRESET

        values.update(LEGACY_PRESET)
        st.sidebar.info(
            "Historical model active: cathode directly on the n- epi (no n+ layer), "
            "pre-flip doping and dark-current fit."
        )
        if nplus_cal_on:
            st.sidebar.warning(
                "Both presets checked: the historical model wins; untick it to apply the "
                "previous calibrated model."
            )
    elif nplus_cal_on:
        from etna import NPLUS_CALIBRATED_PRESET

        values.update(NPLUS_CALIBRATED_PRESET)
        st.sidebar.info(
            "Previous calibrated model active: 3-point n+ calibration, abrupt epi/n+ step "
            "(no graded transition), N_t 2.805336e12. Its C-V knee (-8..-10 V) disagrees "
            "with the digitised measurement (-1..-4 V)."
        )

    try:
        st.session_state["device_config"] = assemble_config(values)
        st.session_state.pop("device_config_error", None)
    except ValueError as e:
        # Never keep a stale valid config: drop it, record why, and every page
        # refuses to run (page_utils.require_config) until the value is fixed.
        st.session_state.pop("device_config", None)
        msg = humanize_config_error(str(e).removeprefix("Invalid DeviceConfig: "))
        st.session_state["device_config_error"] = msg
        st.sidebar.error(f"Invalid device configuration: {msg}")
        return

    # Cross-check domain hint (where the fast engine was checked against devsim).
    from etna.api.simulation import validated_domain_reason

    cfg = st.session_state["device_config"]
    if cfg.half_width_um is None:
        reason = validated_domain_reason(cfg)
        if reason is None:
            st.sidebar.caption(
                "Inside the 1D engine/devsim cross-check domain (the two solvers "
                "agree there; not a validation against measurement)."
            )
        else:
            st.sidebar.caption(
                f"Outside the 1D engine/devsim cross-check domain ({reason}): "
                "C-V/CCE use devsim only, with no cross-check."
            )
