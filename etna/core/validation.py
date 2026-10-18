"""Validation framework for comparing simulation against reference experimental data.

Provides experimental data targets for I-V and C-V characterization of the
4H-SiC p+/n- diode, and functions for computing agreement metrics (R-squared,
RMSE, max deviation, relative errors) without sklearn dependency.

Experimental targets from:
    Petringa et al. -- 4H-SiC Schottky/p-n diode characterization:
    - I-V: dark current < 18 pA at -60V, rectification ratio ~1e5 at +/-2V
    - C-V: W(0V)=1.7um, W(-10V)=9.5um, W(-30V)=9.73um at 1 kHz

All units SI/CGS as noted per quantity.
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)


# --- Experimental Data Targets ---

EXPERIMENTAL_IV = {
    "dark_current_60V": 18e-12,  # A (< 18 pA at -60V, absolute value)
    "rectification_ratio_2V": 1e5,  # I(+2V) / I(-2V), dimensionless
    "series_resistance": 3e3,  # Ohm (from high-forward slope)
}

EXPERIMENTAL_CV = {
    "voltages": [0.0, -10.0, -30.0],  # V
    "depletion_widths_cm": [1.7e-4, 9.5e-4, 9.73e-4],  # cm
    "frequency_hz": 1000,  # Hz
}


def compute_agreement_metrics(sim_values, exp_values):
    """Compute agreement metrics between simulated and experimental data.

    All metrics computed with numpy (no sklearn dependency).

    Parameters
    ----------
    sim_values : array_like
        Simulated values.
    exp_values : array_like
        Experimental/reference values.

    Returns
    -------
    metrics : dict
        Dictionary with:
        - "r_squared": coefficient of determination (1.0 = perfect)
        - "max_deviation": maximum absolute difference
        - "max_relative_deviation": maximum |sim-exp|/|exp| (excludes exp==0)
        - "rmse": root mean square error
        - "mean_relative_error": mean |sim-exp|/|exp| (excludes exp==0)
    """
    sim = np.asarray(sim_values, dtype=float)
    exp = np.asarray(exp_values, dtype=float)

    if len(sim) != len(exp):
        raise ValueError(f"Array length mismatch: sim={len(sim)}, exp={len(exp)}")

    # R-squared (coefficient of determination)
    ss_res = np.sum((exp - sim) ** 2)
    ss_tot = np.sum((exp - np.mean(exp)) ** 2)
    if ss_tot == 0:
        r_squared = 1.0 if ss_res == 0 else 0.0
    else:
        r_squared = 1.0 - ss_res / ss_tot

    # Absolute deviations
    abs_dev = np.abs(sim - exp)
    max_deviation = float(np.max(abs_dev))
    rmse = float(np.sqrt(np.mean(abs_dev**2)))

    # Relative deviations (exclude zero-valued experimental points)
    nonzero = np.abs(exp) > 0
    if np.any(nonzero):
        rel_dev = abs_dev[nonzero] / np.abs(exp[nonzero])
        max_relative_deviation = float(np.max(rel_dev))
        mean_relative_error = float(np.mean(rel_dev))
    else:
        max_relative_deviation = float("inf")
        mean_relative_error = float("inf")

    return {
        "r_squared": float(r_squared),
        "max_deviation": max_deviation,
        "max_relative_deviation": max_relative_deviation,
        "rmse": rmse,
        "mean_relative_error": mean_relative_error,
    }


# A validation target counts as reached only if a simulated bias lies this close to
# it (V); otherwise the metric is NaN and the target is listed in "missing_targets"
# (audit-v7 A-12: the nearest point used to be taken at any distance).
TARGET_V_TOL = 0.25


def _pick(V, target, tol=TARGET_V_TOL):
    """Index of the simulated bias nearest ``target``, or None if farther than tol."""
    if len(V) == 0:
        return None
    i = int(np.argmin(np.abs(V - target)))
    return i if abs(V[i] - target) <= tol else None


def validate_iv(iv_data, area=1.0):
    """Validate simulated I-V data against experimental targets.

    Parameters
    ----------
    iv_data : dict
        Output from iv_sweep() with "voltages" and "currents" (A/cm^2).
    area : float
        Device area (cm^2) for converting current density to absolute current.

    Returns
    -------
    result : dict
        Dictionary with:
        - "dark_current_60V": simulated dark current at -60V (A)
        - "dark_current_target": experimental target (A)
        - "dark_current_pass": bool, True if sim <= target (order-of-magnitude)
        - "rectification_ratio": I(+2V)/I(-2V)
        - "rectification_target": experimental target
        - "rectification_pass": bool
        - "series_resistance": estimated from dV/dI at high forward bias (Ohm)
        - "series_resistance_target": experimental target (Ohm)
        - "missing_targets": biases (-60, +2, -2 V) with no simulated point within
          TARGET_V_TOL; the metrics built on them are NaN and their pass flags False
          (e.g. an iv_sweep truncated before -60 V)
    """
    V = np.asarray(iv_data["voltages"], dtype=float)
    I = np.asarray(iv_data["currents"]) * area  # Convert to absolute current (A)

    result = {}
    missing = []

    # Dark current at -60V
    idx_60 = _pick(V, -60.0)
    if idx_60 is None:
        missing.append(-60.0)
        I_dark = float("nan")
    else:
        I_dark = abs(I[idx_60])
    result["dark_current_60V"] = I_dark
    result["dark_current_target"] = EXPERIMENTAL_IV["dark_current_60V"]
    # Pass if within 2 orders of magnitude (simulation vs measurement tolerance)
    result["dark_current_pass"] = bool(
        I_dark < EXPERIMENTAL_IV["dark_current_60V"] * 100
    )

    # Ideal-SRH floor detection
    # For 4H-SiC with n_i ~ 5e-9 cm^-3, the ideal SRH generation current is
    # ~1e-49 A -- far below any measurable dark current. Experimental dark
    # currents are dominated by surface leakage and trap-assisted tunneling,
    # which are not modeled here. If the simulated dark current is more than
    # 10 orders of magnitude below the experimental target, it indicates the
    # simulation is at the ideal-SRH floor, not a physically meaningful result.
    target = EXPERIMENTAL_IV["dark_current_60V"]
    result["ideal_srh_floor"] = bool(I_dark > 0 and I_dark < target * 1e-10)
    result["dark_current_physically_meaningful"] = bool(
        np.isfinite(I_dark) and not result["ideal_srh_floor"]
    )

    # Rectification ratio at +/- 2V
    idx_p2 = _pick(V, 2.0)
    idx_m2 = _pick(V, -2.0)
    if idx_p2 is None:
        missing.append(2.0)
    if idx_m2 is None:
        missing.append(-2.0)
    if idx_p2 is None or idx_m2 is None:
        rectification = float("nan")
    else:
        I_forward = abs(I[idx_p2])
        I_reverse = abs(I[idx_m2])
        rectification = I_forward / I_reverse if I_reverse > 0 else float("inf")
    result["rectification_ratio"] = rectification
    result["rectification_target"] = EXPERIMENTAL_IV["rectification_ratio_2V"]
    # Pass if within 2 orders of magnitude
    result["rectification_pass"] = bool(
        rectification > EXPERIMENTAL_IV["rectification_ratio_2V"] / 100
    )

    # Series resistance from high-forward slope (dV/dI)
    forward_mask = V > 1.5  # Use voltages above turn-on
    if np.sum(forward_mask) >= 2:
        V_fwd = V[forward_mask]
        I_fwd = I[forward_mask]
        # Linear fit: V = R_s * I + V_offset
        if np.max(I_fwd) > np.min(I_fwd):
            coeffs = np.polyfit(I_fwd, V_fwd, 1)
            R_s = coeffs[0]
        else:
            R_s = float("nan")
    else:
        R_s = float("nan")

    result["series_resistance"] = R_s
    result["series_resistance_target"] = EXPERIMENTAL_IV["series_resistance"]
    result["missing_targets"] = missing

    return result


def validate_cv(cv_data):
    """Validate simulated C-V data against experimental depletion width targets.

    Parameters
    ----------
    cv_data : dict
        Dictionary with "voltages" and "W_from_C" (eps*A/C, cm; cv_sweep output) or,
        for manually constructed data, "depletion_widths" (cm). The experimental W are
        eps*A/C of a measured C, so "W_from_C" is used when present (the depletion edge
        "depletion_widths" is not eps*A/C near full depletion).

    Returns
    -------
    result : dict
        Dictionary with:
        - "sim_W": simulated depletion widths at experimental voltages (cm)
        - "exp_W": experimental depletion widths (cm)
        - "exp_voltages": experimental voltages (V)
        - "metrics": output of compute_agreement_metrics()
        - "per_point_error": relative error at each experimental voltage
        - "extrapolated": True if an experimental voltage lies outside the simulated
          bias range (or cv_data["truncated"] is True). Such points are NaN in
          "sim_W" / "per_point_error" and excluded from "metrics" (np.interp would
          silently return the edge value; audit-v7 A-12)
        - "n_extrapolated": number of experimental voltages outside the range

    Raises
    ------
    ValueError
        If no experimental voltage lies inside the simulated range.
    """
    V_sim = np.asarray(cv_data["voltages"], dtype=float)
    W_sim = np.asarray(cv_data.get("W_from_C", cv_data["depletion_widths"]))

    V_exp = np.array(EXPERIMENTAL_CV["voltages"])
    W_exp = np.array(EXPERIMENTAL_CV["depletion_widths_cm"])

    if len(V_sim) == 0:
        raise ValueError("validate_cv: empty C-V data (no simulated voltage)")
    out = (V_exp < V_sim.min() - 1e-9) | (V_exp > V_sim.max() + 1e-9)
    if out.all():
        raise ValueError(
            f"validate_cv: no experimental voltage {V_exp.tolist()} inside the "
            f"simulated range [{V_sim.min():g}, {V_sim.max():g}] V"
        )

    # Interpolate simulated W at the experimental voltages inside the range only
    order = np.argsort(V_sim)
    W_sim_at_exp = np.where(out, np.nan, np.interp(V_exp, V_sim[order], W_sim[order]))

    inside = ~out
    metrics = compute_agreement_metrics(W_sim_at_exp[inside], W_exp[inside])

    per_point_error = np.abs(W_sim_at_exp - W_exp) / W_exp

    return {
        "sim_W": W_sim_at_exp,
        "exp_W": W_exp,
        "exp_voltages": V_exp,
        "metrics": metrics,
        "per_point_error": per_point_error.tolist(),
        "extrapolated": bool(out.any() or cv_data.get("truncated", False)),
        "n_extrapolated": int(out.sum()),
    }
