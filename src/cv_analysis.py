"""Capacitance-voltage (C-V) analysis for 4H-SiC p+/n- diode.

Computes junction capacitance from depletion width using the parallel-plate
approximation C = eps * A / W, and provides a C-V sweep function that
uses the numerical Poisson/DD solver to extract W(V) at each bias point.

The 1/C^2 vs V (Mott-Schottky) analysis enables extraction of doping
concentration and built-in potential from the linear fit slope/intercept.

All units CGS (cm, F/cm^2, V) per devsim convention unless otherwise noted.

References:
    - Petringa et al.: C-V measurements at 1 kHz, W(0V)=1.7um, W(-10V)=9.5um
    - Sze & Ng, "Physics of Semiconductor Devices", Ch. 2: Junction capacitance
"""

import logging

import numpy as np

from src.poisson import ramp_voltage, extract_depletion_width_numerical

logger = logging.getLogger(__name__)

# Physical constants (CGS)
EPS_0 = 8.854e-14  # F/cm (vacuum permittivity)
Q = 1.602e-19  # C (elementary charge)


def junction_capacitance(W, eps_r=9.7, area=1.0):
    """Compute junction capacitance from depletion width.

    Uses the parallel-plate approximation: C = eps_r * eps_0 * area / W.

    Parameters
    ----------
    W : float or array_like
        Depletion width (cm).
    eps_r : float
        Relative permittivity of 4H-SiC. Default 9.7.
    area : float
        Junction area (cm^2). Default 1.0 (gives F/cm^2).

    Returns
    -------
    C : float or ndarray
        Capacitance (F if area given, F/cm^2 if area=1.0).
    """
    W = np.asarray(W, dtype=float)
    return eps_r * EPS_0 * area / W


def depletion_width_from_capacitance(C, eps_r=9.7, area=1.0):
    """Compute depletion width from capacitance (inverse of junction_capacitance).

    Parameters
    ----------
    C : float or array_like
        Capacitance (F or F/cm^2).
    eps_r : float
        Relative permittivity of 4H-SiC. Default 9.7.
    area : float
        Junction area (cm^2). Default 1.0.

    Returns
    -------
    W : float or ndarray
        Depletion width (cm).
    """
    C = np.asarray(C, dtype=float)
    return eps_r * EPS_0 * area / C


def compute_cv_from_depletion(voltages, depletion_widths, eps_r=9.7, area=1.0):
    """Compute C-V curve from voltage and depletion width arrays.

    Parameters
    ----------
    voltages : array_like
        Applied bias voltages (V).
    depletion_widths : array_like
        Depletion widths at each voltage (cm).
    eps_r : float
        Relative permittivity. Default 9.7.
    area : float
        Junction area (cm^2). Default 1.0.

    Returns
    -------
    result : dict
        Dictionary with:
        - "voltages": numpy array of voltages (V)
        - "capacitance": numpy array of C values (F or F/cm^2)
        - "one_over_C_squared": numpy array of 1/C^2 values
    """
    voltages = np.asarray(voltages, dtype=float)
    W = np.asarray(depletion_widths, dtype=float)

    C = junction_capacitance(W, eps_r=eps_r, area=area)

    return {
        "voltages": voltages,
        "capacitance": C,
        "one_over_C_squared": 1.0 / C**2,
    }


def cv_sweep(device_info, V_range, eps_r=9.7, area=1.0):
    """Sweep reverse bias and compute C-V from numerical depletion widths.

    At each voltage, ramps bias using the Poisson/DD solver, extracts the
    depletion width numerically from the E-field profile, and computes
    capacitance via the parallel-plate approximation.

    Parameters
    ----------
    device_info : dict
        Device info dict (Poisson or DD equations must be set up).
    V_range : array_like
        Array of reverse bias voltages (V, should be <= 0).
    eps_r : float
        Relative permittivity. Default 9.7.
    area : float
        Junction area (cm^2). Default 1.0.

    Returns
    -------
    result : dict
        Dictionary with:
        - "voltages": numpy array of voltages (V)
        - "depletion_widths": numpy array of W values (cm)
        - "capacitance": numpy array of C values (F or F/cm^2)
    """
    V_range = np.asarray(V_range, dtype=float)

    depletion_widths = []
    solved_voltages = []

    # Extract W at 0V first (current equilibrium state)
    W0 = extract_depletion_width_numerical(device_info)
    current_V = 0.0

    for V_target in V_range:
        if abs(V_target) < 1e-12:
            # Use equilibrium W
            depletion_widths.append(W0)
            solved_voltages.append(0.0)
            continue

        # Ramp from current state to target
        V_step = -0.5 if V_target < current_V else 0.5
        ramp_results = ramp_voltage(
            device_info,
            contact_name="cathode",
            V_start=current_V,
            V_end=V_target,
            V_step=V_step,
        )

        if ramp_results:
            # Extract W from the final solved state
            W = extract_depletion_width_numerical(device_info)
            depletion_widths.append(W)
            solved_voltages.append(V_target)
            current_V = V_target
        else:
            logger.warning(f"cv_sweep: no results at V={V_target:.2f}V")

    voltages = np.array(solved_voltages)
    W_arr = np.array(depletion_widths)
    C_arr = junction_capacitance(W_arr, eps_r=eps_r, area=area)

    return {
        "voltages": voltages,
        "depletion_widths": W_arr,
        "capacitance": C_arr,
    }
