"""Charge collection efficiency (CCE) computation utilities.

Provides:
- Hecht equation (two-carrier form) for analytical CCE benchmark
- Extended Hecht for partially-depleted detectors with diffusion
- CCE extraction from DD simulation current ratio

All units in CGS (cm, cm^-3, V, s) per devsim convention.

References:
    - Hecht, Z. Physik 77, 235 (1932)
    - Sze & Ng, "Physics of Semiconductor Devices", 3rd ed.
    - Knoll, "Radiation Detection and Measurement", 4th ed.
"""

import numpy as np


def hecht_cce(V, d, mu_e=950.0, tau_e=1e-9, mu_p=125.0, tau_p=6e-7):
    """Two-carrier Hecht equation for charge collection efficiency.

    Assumes uniform electric field E = |V|/d across the active region.
    Valid for fully-depleted detectors in the low-injection regime.

    Parameters
    ----------
    V : float or array_like
        Applied voltage (V). Uses absolute value for reverse bias.
    d : float
        Active region thickness (cm), typically depletion width.
    mu_e : float
        Electron mobility (cm^2/Vs). Default: 950 (4H-SiC low doping).
    tau_e : float
        Electron SRH lifetime (s). Default: 1e-9 (4H-SiC).
    mu_p : float
        Hole mobility (cm^2/Vs). Default: 125 (4H-SiC low doping).
    tau_p : float
        Hole SRH lifetime (s). Default: 6e-7 (4H-SiC).

    Returns
    -------
    cce : float or ndarray
        Charge collection efficiency, clipped to [0, 1].

    Notes
    -----
    CCE = (lambda_e/d)*(1-exp(-d/lambda_e)) + (lambda_h/d)*(1-exp(-d/lambda_h))

    where lambda_e = mu_e * tau_e * |V| / d (electron drift length),
    and similarly for holes.

    For 4H-SiC at V=40V, d=10um: lambda_e ~ 380 um >> d, so CCE -> 1.0.
    """
    V = np.abs(np.asarray(V, dtype=float))
    d = float(d)

    # Drift lengths: lambda = mu * tau * E = mu * tau * |V| / d
    lambda_e = mu_e * tau_e * V / d
    lambda_h = mu_p * tau_p * V / d

    # Avoid division by zero at V=0 (lambda=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        cce_e = np.where(
            lambda_e > 0, (lambda_e / d) * (1.0 - np.exp(-d / lambda_e)), 0.0
        )
        cce_h = np.where(
            lambda_h > 0, (lambda_h / d) * (1.0 - np.exp(-d / lambda_h)), 0.0
        )

    cce = cce_e + cce_h
    return np.clip(cce, 0.0, 1.0)


def compute_cce_from_current(I_collected, I_generated):
    """Compute CCE as ratio of collected to generated current.

    Used by the DD solver (Plan 02) to extract CCE from simulation.

    Parameters
    ----------
    I_collected : float or array_like
        Collected current (A/cm^2) at the contact.
    I_generated : float or array_like
        Total generated current (A/cm^2) from radiation.

    Returns
    -------
    cce : float or ndarray
        Charge collection efficiency, clipped to [0, 1].
    """
    I_collected = np.abs(np.asarray(I_collected, dtype=float))
    I_generated = np.asarray(I_generated, dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        cce = np.where(I_generated > 0, I_collected / I_generated, 0.0)

    return np.clip(cce, 0.0, 1.0)


def hecht_cce_partial_depletion(
    V, d_epi, W_func, mu_e=950.0, tau_e=1e-9, mu_p=125.0, tau_p=6e-7, L_diff_p=7e-4
):
    """Extended Hecht equation for partially-depleted detector.

    Combines drift collection in the depletion region with diffusion
    collection from the neutral region. This is an approximation for
    comparison with numerical DD results.

    Parameters
    ----------
    V : float or array_like
        Applied voltage (V). Uses absolute value.
    d_epi : float
        Total epitaxial layer thickness (cm).
    W_func : callable
        Function W(V) returning depletion width (cm) at voltage V.
        Can be obtained from C-V data or DD simulation.
    mu_e : float
        Electron mobility (cm^2/Vs). Default: 950.
    tau_e : float
        Electron lifetime (s). Default: 1e-9.
    mu_p : float
        Hole mobility (cm^2/Vs). Default: 125.
    tau_p : float
        Hole lifetime (s). Default: 6e-7.
    L_diff_p : float
        Minority carrier (hole) diffusion length (cm). Default: 7e-4 (7 um).
        From literature CCE fitting (Sciencedirect S0168900205006443).

    Returns
    -------
    cce : float or ndarray
        Charge collection efficiency, clipped to [0, 1].

    Notes
    -----
    The model assumes:
    - Uniform generation within the epitaxial layer
    - Drift collection (standard Hecht) for charge in the depletion region
    - Exponential diffusion collection probability exp(-x/L_diff) for charge
      in the neutral region (from W to d_epi)
    - This is an approximation; the numerical DD solver is more accurate
    """
    V_arr = np.atleast_1d(np.abs(np.asarray(V, dtype=float)))
    cce_vals = np.zeros_like(V_arr)

    for i, v in enumerate(V_arr):
        W = float(W_func(v))
        W = min(W, d_epi)  # can't exceed epi thickness

        if W <= 0 or v <= 0:
            cce_vals[i] = 0.0
            continue

        # Drift component: standard Hecht within depletion region
        cce_drift = float(hecht_cce(v, W, mu_e, tau_e, mu_p, tau_p))

        # Fraction of generation in depletion region
        f_depl = W / d_epi

        if W >= d_epi:
            # Fully depleted: all charge collected by drift
            cce_vals[i] = cce_drift
        else:
            # Diffusion component: charge in neutral region (W to d_epi)
            # Collection probability = exp(-(x-W)/L_diff), integrated over [W, d_epi]
            neutral_thickness = d_epi - W
            if L_diff_p > 0:
                # Average collection probability in neutral region
                # = (1/t) * integral_0^t exp(-x/L) dx = (L/t)*(1-exp(-t/L))
                avg_coll = (L_diff_p / neutral_thickness) * (
                    1.0 - np.exp(-neutral_thickness / L_diff_p)
                )
            else:
                avg_coll = 0.0

            f_neutral = 1.0 - f_depl
            cce_vals[i] = f_depl * cce_drift + f_neutral * avg_coll

    result = np.clip(cce_vals, 0.0, 1.0)

    # Return scalar if input was scalar
    if np.ndim(V) == 0:
        return float(result[0])
    return result
