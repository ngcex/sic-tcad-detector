"""The validation label every export carries, and the device-model facts it qualifies."""

from __future__ import annotations

from dataclasses import asdict

import numpy as np

# Verbatim on every export (JSON, CSV, Markdown, PDF, SVG, PNG, VTK, animations).
VALIDATION_LABEL = (
    "simulation output, validated against devsim reference; "
    "not validated against measurement except where stated"
)

# What "except where stated" refers to (the only measurement comparisons).
MEASURED_W_30V_UM = 9.73  # Photons Fig. 4b, W = eps*A/C at -30 V (A = 4 mm2)
CALIBRATION_NOTE = (
    "Doping (N_D_junction, N_D_bulk) was solved exactly on the measured C-V points at 0 V "
    "and -10 V (W = eps*A/C = 1.7 / 9.5 um, with C the small-signal dQ/dV): these two points "
    "are calibration inputs, not reproduced or validated points."
)
DARK_CALIBRATION_NOTE = (
    "Dark current is a calibration / budget estimate, not a prediction: N_t is fitted so "
    "that I(300 K, -60 V, 4 mm2) = 18 pA, the measured bound '< 18 pA up to -60 V' used "
    "as equality (N_t is the magnitude of an effective-generation proxy, not a trap density)."
)
CV_SYSTEMATIC_NOTE = (
    "Known C-V residual: W = eps*A/C at -30 V simulated {w:.2f} um vs measured 9.73 um "
    "({d:+.1f} %), the one out-of-fit measurement comparison. Three C-V points cannot tell "
    "apart its candidate causes (a ~0.5 um graded n-/n+ buffer, an epi ~2.5 % thinner than "
    "nominal, lateral edge/fringe capacitance falling with bias): none is picked; not refitted."
)
CV_SYSTEMATIC_OTHER = (
    "Known C-V residual (+2.6 % in eps*A/C at -30 V, cause not identified, not refitted) is "
    "stated for the reference device (default stack, 300 K); not re-evaluated for this config."
)
ALPHA_NOTE = (
    "Epi-only alpha generation matches the measurement geometry (alphas tilted 65 deg, "
    "projected depth ~7.6 um, never reach the n+); CCE is relative to that in-epi generation."
)
NPLUS_THICKNESS_NOTE = (
    "{t:g} um modelled: the conductive n+ substrate truncated in the model (W, CCE and dark "
    "current do not depend on it); not a fabrication thickness"
)


def has_nplus(config) -> bool:
    return getattr(config, "nplus_thickness_um", None) is not None


def is_reference_device(config) -> bool:
    """The calibrated default device (geometry, doping, 300 K); area / N_t / E_t free.

    Only for this device do the stated measurement comparisons apply.
    """
    from etna.api.device import DeviceConfig

    skip = {"area_cm2", "N_t", "E_t"}
    a = {k: v for k, v in asdict(config).items() if k not in skip}
    b = {k: v for k, v in asdict(DeviceConfig()).items() if k not in skip}
    return a == b


def reach_through_definition(config) -> str:
    if has_nplus(config):
        return (
            "smallest |V| at which the depletion edge (n/N_D = 0.5) lies within one "
            "extrinsic Debye length L_D = sqrt(eps V_t / (q N_D(epi end))) of the "
            "n-/n+ step. A definition, not a measurement (12.30 V on the default stack; "
            "paper: ~10 V)"
        )
    return (
        "smallest |V| at which the depletion edge lies in the last mesh interval before "
        "the cathode (legacy model: ohmic cathode directly on the n- epi; 17.0 V on the "
        "legacy device). A definition, not a measurement"
    )


VFD_DEFINITION = (
    "analytic graded full-depletion voltage V_fd = (q/eps) int_0^t y N_D(y) dy - V_bi "
    "(Poisson double integral of the graded epi profile; V_bi from ionized N_A and "
    "N_D_junction). The looser of the two definitions; used by the 2D sweep's "
    "full-depletion gate (|V| >= 1.2 V_fd)"
)


def epi_donor_profile(config, depth_um):
    """N_D (cm^-3) at ``depth_um`` below the metallurgical junction (p+/epi interface)."""
    d = np.asarray(depth_um, dtype=float)
    if config.doping_profile == "uniform":
        from etna.core.device import DEFAULT_N_D

        nd = DEFAULT_N_D if config.N_D is None else config.N_D
        return np.full_like(d, nd, dtype=float)
    return config.N_D_bulk + (config.N_D_junction - config.N_D_bulk) * np.exp(
        -d / config.L_transition_um
    )


def doping_formula(config) -> str:
    if config.doping_profile == "uniform":
        return "uniform N_D"
    return "N_D(y) = N_D_bulk + (N_D_junction - N_D_bulk) * exp(-y / L_transition), y = depth below the p+/epi junction"


def analytic_vfd(config) -> "float | None":
    """Analytic graded V_fd (V), same inputs as the 2D sweep gate (etna.core.optimization)."""
    from etna.core.analytical import built_in_potential, full_depletion_voltage_graded
    from etna.core.incomplete_ionization import ionized_acceptor_concentration
    from etna.core.sic_material import SiC4H_Parameters, intrinsic_concentration

    try:
        n_i = intrinsic_concentration(config.T, SiC4H_Parameters())[0]
        na = ionized_acceptor_concentration(config.N_A, config.T)
        uniform = config.doping_profile == "uniform"
        nd_j = float(epi_donor_profile(config, 0.0)) if uniform else config.N_D_junction
        v_bi = built_in_potential(na, nd_j, n_i, T=config.T)
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
