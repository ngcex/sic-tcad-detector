"""The validation label every export carries, and the device-model facts it qualifies."""

from __future__ import annotations

from dataclasses import asdict

import numpy as np

# Verbatim on every export (JSON, CSV, Markdown, PDF, SVG, PNG, VTK, animations).
VALIDATION_LABEL = (
    "simulation output; fast engines cross-checked against the devsim reference solver; "
    "not validated against measurement except where stated"
)

# What "except where stated" refers to (the only measurement comparisons).
MEASURED_W_30V_UM = 9.73  # Photons Fig. 4b, W = eps*A/C at -30 V (A = 4 mm2)
CALIBRATION_NOTE = (
    "Default stack doping (N_D_junction, N_D_bulk, L_transition and the graded n-/n+ "
    "transition) is a FIT to the digitised Photons Fig. 4b C-V curve (eps*A/C every 250 mV, "
    "0 to -60 V, A = 4 mm2, C the small-signal dQ/dV; chi2 198 / 224 dof): plateau rms "
    "0.013 um over -5..-59 V, knee rms 0.24 um (-3.25..-0.25 V). The whole curve is a "
    "calibration input, not a validation. N_D_bulk 3.4e13 cm^-3 is below the documented "
    "0.5-1e14 cm^-3 statements for this epi; the epi/n+ position (10 um epi, transition over "
    "its last 1.14 um, 1e19 beyond) is set by the fit, not by SIMS (C-V does not show where "
    "1e19 is reached). Residuals at the labelled points 0 / -10 / -30 V: +0.025 / +0.001 / "
    "-0.020 um (model 1.655 / 9.519 / 9.750 vs 1.68 / 9.52 / 9.73 um)."
)
NPLUS_CALIBRATION_NOTE = (
    "Previous default (nplus_calibrated_config): doping (N_D_junction, N_D_bulk) solved "
    "exactly on the measured C-V points at 0 V and -10 V (W = eps*A/C = 1.7 / 9.5 um, C the "
    "small-signal dQ/dV; L_transition held at 0.987 um): these two points are calibration "
    "inputs, not reproduced or validated points."
)
DARK_CALIBRATION_NOTE = (
    "Dark current is a calibration / budget estimate, not a prediction: N_t is fitted so "
    "that I(300 K, -60 V, 4 mm2) = 18 pA, the measured bound '< 18 pA up to -60 V' used "
    "as equality (N_t is the magnitude of an effective-generation proxy, not a trap density)."
)
CV_SYSTEMATIC_NOTE = (
    "No out-of-fit C-V comparison: the whole digitised curve is a fit input (see "
    "cv.calibration). W = eps*A/C at -30 V simulated {w:.3f} um vs measured 9.73 um "
    "({d:+.2f} %, an in-fit residual). The knee (-3.25..-0.25 V) is reproduced only to "
    "0.24 um rms: the exponential epi form cannot follow the measured dip to ~1.3e13 at 4-6 um."
)
CV_SYSTEMATIC_NPLUS_NOTE = (
    "Known C-V residual: W = eps*A/C at -30 V simulated {w:.2f} um vs measured 9.73 um "
    "({d:+.1f} %), out of fit. The full digitised Fig. 4b curve excludes this abrupt "
    "epi/n+ form: it puts the C-V knee at -8..-10 V where the measurement has it at "
    "-1..-4 V (W(-2 V) 2.99 vs 7.18 um); the default stack's graded n-/n+ transition fits it "
    "(lab/cv_gap/notes/LOG.md)."
)
CV_SYSTEMATIC_OTHER = (
    "C-V comparisons with the measurement (digitised Photons Fig. 4b) are stated for the "
    "reference device (default stack, 300 K) only; not re-evaluated for this config."
)
ALPHA_NOTE = (
    "Epi-only alpha generation matches the measurement geometry (alphas tilted 65 deg, "
    "projected depth ~7.6 um, never reach the n+); CCE is relative to that in-epi generation."
)
NPLUS_THICKNESS_NOTE = (
    "{t:g} um modelled: the conductive n+ substrate truncated in the model (W, CCE and dark "
    "current do not depend on it); not a fabrication thickness"
)
BUFFER_NOTE = (
    "graded n-/n+ transition: N_D rises log-linearly from the epi value to {n_end:.3g} cm^-3 "
    "at the epi end over the last {t:g} um of the epi (the 1e19 substrate is not reached at "
    "the nominal epi end); a FITTED form from the digitised measured C-V, not a measured "
    "profile (SIMS / SRP would settle it)"
)
# The measured-C-V profile is the default since 2026-10-18: same statement.
MEASURED_CV_NOTE = CALIBRATION_NOTE
# The OPT-IN tabulated profile (measured_cv_tabulated_config): what it is and is not.
TABULATED_CALIBRATION_NOTE = (
    "Opt-in TABULATED epi profile: 19 knots (depth below the junction, N_D), log-linear "
    "between them, from the model-independent inversion of the digitised Photons Fig. 4b "
    "C-V (eps*A/C every 250 mV, 0 to -60 V, A = 4 mm2; lab/cv_gap/notes/knots.json, curvature "
    "penalty lam 1): chi2 32 / 211 dof, knee rms 0.014 um (-3.25..-0.25 V), plateau rms "
    "0.007 um (-5..-59 V); the parametric default leaves 0.24 um rms in the knee. The whole "
    "curve is the fit input, not a validation. HONESTY: this is an inversion of ONE C-V "
    "curve (a regularised fit of W(V)), NOT a SIMS / SRP profile. Robustness: the knots are "
    "stable to +-0.2 dex for smoothing lam 0.3-3 EXCEPT the first 0.5 um (inside the "
    "zero-bias depletion width, unresolved by C-V) and the 8.5-9.25 um stretch (degenerate "
    "with a fringe term). The 1e19 substrate is not reached at 10 um (table end 1.3e16); "
    "where it is reached is not known from C-V."
)
CV_SYSTEMATIC_TABULATED_NOTE = (
    "No out-of-fit C-V comparison: the whole digitised curve is the inversion's input "
    "(see cv.calibration). W = eps*A/C at -30 V simulated {w:.3f} um vs measured 9.73 um "
    "({d:+.2f} %, an in-fit residual); knee and plateau reproduced within the digitisation "
    "uncertainty (13/13 knee markers within 1 sigma)."
)
TABLE_NOTE = (
    "tabulated N_D (log-linear between {k} knots): the low-doped epi ends at {t_end:.2f} um "
    "below the junction where N_D has risen one decade above the table minimum "
    "({n_min:.3g} -> {n_end:.3g} cm^-3); beyond it the table rises to {n_last:.3g} cm^-3 at "
    "the epi end (the n+ step follows). An inversion of one measured C-V curve, not a SIMS "
    "profile"
)


def has_nplus(config) -> bool:
    return getattr(config, "nplus_thickness_um", None) is not None


def has_buffer(config) -> bool:
    return getattr(config, "buffer_thickness_um", None) is not None


def has_table(config) -> bool:
    """True when ``config`` carries a tabulated epi donor profile (``donor_table``)."""
    return getattr(config, "donor_table", None) is not None


def table_spec(config):
    """``donor_table_spec`` of ``config`` (None without a table)."""
    from etna.core.device import donor_table_spec

    if not has_table(config):
        return None
    return donor_table_spec(config.epi_thickness_um * 1e-4, config.donor_table)


def table_sensitive_end_um(config) -> "float | None":
    """Depth (um below the junction) where the low-doped epi of a tabulated profile
    ends (one decade above the table minimum; etna.core.device.table_sensitive_end_cm).
    """
    from etna.core.device import table_sensitive_end_cm

    spec = table_spec(config)
    return None if spec is None else table_sensitive_end_cm(spec) * 1e4


def _same_device(config, other) -> bool:
    skip = {"area_cm2", "N_t", "E_t"}
    a = {k: v for k, v in asdict(config).items() if k not in skip}
    b = {k: v for k, v in asdict(other).items() if k not in skip}
    return a == b


def is_reference_device(config) -> bool:
    """The default device (geometry, doping, transition, 300 K); area / N_t / E_t free.

    Only for this device do the stated measurement comparisons apply. Since 2026-10-18
    the default is the measured-C-V profile (DEFAULT_STACK).
    """
    from etna.api.device import DeviceConfig

    return _same_device(config, DeviceConfig())


def is_measured_cv_device(config) -> bool:
    """The measured-C-V profile (== the reference device since the 2026-10-18 flip)."""
    return is_reference_device(config)


def is_nplus_calibrated_device(config) -> bool:
    """The previous default (nplus_calibrated_config: 3-point n+ calibration, abrupt
    epi/n+); area / N_t / E_t free."""
    from etna.api.device import nplus_calibrated_config

    return _same_device(config, nplus_calibrated_config())


def is_measured_cv_tabulated_device(config) -> bool:
    """The opt-in tabulated measured-C-V preset (measured_cv_tabulated_config: the 19-knot
    inversion on the default geometry); area / N_t / E_t free."""
    from etna.api.device import measured_cv_tabulated_config

    return _same_device(config, measured_cv_tabulated_config())


def reach_through_definition(config) -> str:
    if has_table(config):
        t_end = table_sensitive_end_um(config)
        return (
            "smallest |V| at which the depletion edge (n/N_D = 0.5) lies within one "
            "extrinsic Debye length L_D of the END OF THE LOW-DOPED EPI of the tabulated "
            f"profile, the depth ({t_end:.2f} um) where N_D has risen one decade above the "
            "table minimum (L_D at the doping there). A definition, not a measurement, and "
            "not the same point as the parametric transition start (a table has no "
            "analytic start of the rise; the decade rule is stable to 0.04 um against the "
            "inversion's smoothing). 4.85 V on the tabulated preset (parametric default "
            "3.10 V)"
        )
    if has_buffer(config):
        return (
            "smallest |V| at which the depletion edge (n/N_D = 0.5) lies within one "
            "extrinsic Debye length L_D of the START of the graded n-/n+ transition "
            f"(epi end - {config.buffer_thickness_um:g} um; L_D at the epi doping there): "
            "full depletion of the low-doped epi; the edge then creeps into the "
            "transition for tens of volts. A definition, not a measurement (3.10 V on "
            "the default stack)"
        )
    if has_nplus(config):
        return (
            "smallest |V| at which the depletion edge (n/N_D = 0.5) lies within one "
            "extrinsic Debye length L_D = sqrt(eps V_t / (q N_D(epi end))) of the "
            "n-/n+ step. A definition, not a measurement (12.30 V on the previous default, "
            "nplus_calibrated_config)"
        )
    return (
        "smallest |V| at which the depletion edge lies in the last mesh interval before "
        "the cathode (legacy model: ohmic cathode directly on the n- epi; 17.0 V on the "
        "legacy device). A definition, not a measurement"
    )


VFD_DEFINITION = (
    "analytic graded full-depletion voltage V_fd = (q/eps) int_0^t y N_D(y) dy - V_bi "
    "(Poisson double integral of the graded epi profile; V_bi from ionized N_A and "
    "N_D_junction). With a graded n-/n+ transition t is the depth of the TRANSITION START "
    "(the low-doped epi; depleting the transition itself is the -10..-60 V creep, not full "
    "depletion). The looser of the two definitions; used by the 2D sweep's full-depletion "
    "gate (|V| >= 1.2 V_fd). Default stack: 2.74 V (gate 3.28 V); previous default "
    "(abrupt epi/n+): 9.93 V. With a TABULATED profile t is the depth where N_D has "
    "risen one decade above the table minimum (the end of the low-doped epi; "
    "etna.core.device.table_sensitive_end_cm): 3.37 V on the tabulated preset (gate "
    "4.05 V; its engine reach-through is 4.85 V)"
)


def epi_donor_profile(config, depth_um):
    """N_D (cm^-3) at ``depth_um`` below the metallurgical junction (p+/epi interface),
    the epi profile with its graded n-/n+ transition if any (no n+ step; a depth at or
    beyond the epi end reads the epi value just before the n+ step)."""
    from etna.core.device import donor_profile, tabulated_donor_profile

    t = config.epi_thickness_um * 1e-4
    d = np.minimum(np.asarray(depth_um, dtype=float) * 1e-4, t * (1.0 - 1e-12))
    if has_table(config):
        return tabulated_donor_profile(d, 0.0, table_spec(config), n_side=False)
    buffer = None
    if has_buffer(config):
        buffer = (t, config.buffer_thickness_um * 1e-4, config.buffer_end_doping)
    if config.doping_profile == "uniform":
        nd = config.N_D  # required with "uniform" (C-12)
        return donor_profile(d, 0.0, nd, nd, 1e-4, buffer=buffer, n_side=False)
    return donor_profile(
        d,
        0.0,
        config.N_D_junction,
        config.N_D_bulk,
        config.L_transition_um * 1e-4,
        buffer=buffer,
        n_side=False,
    )


def doping_formula(config) -> str:
    if has_table(config):
        return (
            f"tabulated N_D(y): {len(config.donor_table)} knots (depth below the p+/epi "
            "junction, N_D), log10 N_D linear between knots, last value held beyond the "
            "last knot; N_D_junction / N_D_bulk / L_transition are informative only"
        )
    buf = (
        (
            f"; over the last {config.buffer_thickness_um:g} um of the epi N_D rises "
            f"log-linearly to {config.buffer_end_doping:.3g} cm^-3 at the epi end "
            "(graded n-/n+ transition)"
        )
        if has_buffer(config)
        else ""
    )
    if config.doping_profile == "uniform":
        return "uniform N_D" + buf
    return (
        "N_D(y) = N_D_bulk + (N_D_junction - N_D_bulk) * exp(-y / L_transition), "
        "y = depth below the p+/epi junction" + buf
    )


def analytic_vfd(config) -> "float | None":
    """Analytic graded V_fd (V), same inputs as the 2D sweep gate (etna.core.optimization)."""
    from etna.core.analytical import (
        built_in_potential,
        full_depletion_voltage_graded,
        full_depletion_voltage_profile,
    )
    from etna.core.incomplete_ionization import ionized_acceptor_concentration
    from etna.core.sic_material import SiC4H_Parameters, intrinsic_concentration

    try:
        n_i = intrinsic_concentration(config.T, SiC4H_Parameters())[0]
        na = ionized_acceptor_concentration(config.N_A, config.T)
        uniform = config.doping_profile == "uniform"
        nd_j = (
            float(epi_donor_profile(config, 0.0))
            if (uniform or has_table(config))
            else config.N_D_junction
        )
        v_bi = built_in_potential(na, nd_j, n_i, T=config.T)
        if has_table(config):  # to the end of the low-doped epi (see VFD_DEFINITION)
            y_um = np.linspace(0.0, table_sensitive_end_um(config), 20001)
            return float(
                full_depletion_voltage_profile(
                    y_um * 1e-4, epi_donor_profile(config, y_um), v_bi
                )
            )
        if has_buffer(config):  # to the transition start (see VFD_DEFINITION)
            t_low = config.epi_thickness_um - config.buffer_thickness_um
            y_um = np.linspace(0.0, t_low, 20001)
            return float(
                full_depletion_voltage_profile(
                    y_um * 1e-4, epi_donor_profile(config, y_um), v_bi
                )
            )
        return float(
            full_depletion_voltage_graded(
                epi_thickness=config.epi_thickness_um * 1e-4,
                N_D_bulk=config.N_D_bulk,
                N_D_junction=config.N_D_junction,
                L_transition=config.L_transition_um * 1e-4,
                V_bi=v_bi,
                n_d_uniform=nd_j if uniform else None,
            )
        )
    except Exception:  # pragma: no cover - defensive (never block an export)
        return None
