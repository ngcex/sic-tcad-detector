"""Validation framework for comparing simulation against Petringa experimental data.

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
    """
    V = np.asarray(iv_data["voltages"])
    I = np.asarray(iv_data["currents"]) * area  # Convert to absolute current (A)

    result = {}

    # Dark current at -60V
    idx_60 = np.argmin(np.abs(V - (-60.0)))
    I_dark = abs(I[idx_60])
    result["dark_current_60V"] = I_dark
    result["dark_current_target"] = EXPERIMENTAL_IV["dark_current_60V"]
    # Pass if within 2 orders of magnitude (simulation vs measurement tolerance)
    result["dark_current_pass"] = I_dark < EXPERIMENTAL_IV["dark_current_60V"] * 100

    # Rectification ratio at +/- 2V
    idx_p2 = np.argmin(np.abs(V - 2.0))
    idx_m2 = np.argmin(np.abs(V - (-2.0)))
    I_forward = abs(I[idx_p2])
    I_reverse = abs(I[idx_m2])
    if I_reverse > 0:
        rectification = I_forward / I_reverse
    else:
        rectification = float("inf")
    result["rectification_ratio"] = rectification
    result["rectification_target"] = EXPERIMENTAL_IV["rectification_ratio_2V"]
    # Pass if within 2 orders of magnitude
    result["rectification_pass"] = (
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

    return result


def validate_cv(cv_data):
    """Validate simulated C-V data against experimental depletion width targets.

    Parameters
    ----------
    cv_data : dict
        Dictionary with "voltages" and "depletion_widths" arrays (cm).
        Can be output from cv_sweep() or manually constructed.

    Returns
    -------
    result : dict
        Dictionary with:
        - "sim_W": simulated depletion widths at experimental voltages (cm)
        - "exp_W": experimental depletion widths (cm)
        - "exp_voltages": experimental voltages (V)
        - "metrics": output of compute_agreement_metrics()
        - "per_point_error": relative error at each experimental voltage
    """
    V_sim = np.asarray(cv_data["voltages"])
    W_sim = np.asarray(cv_data["depletion_widths"])

    V_exp = np.array(EXPERIMENTAL_CV["voltages"])
    W_exp = np.array(EXPERIMENTAL_CV["depletion_widths_cm"])

    # Interpolate simulated W at experimental voltages
    W_sim_at_exp = np.interp(V_exp, np.sort(V_sim), W_sim[np.argsort(V_sim)])

    metrics = compute_agreement_metrics(W_sim_at_exp, W_exp)

    per_point_error = np.abs(W_sim_at_exp - W_exp) / W_exp

    return {
        "sim_W": W_sim_at_exp,
        "exp_W": W_exp,
        "exp_voltages": V_exp,
        "metrics": metrics,
        "per_point_error": per_point_error.tolist(),
    }
