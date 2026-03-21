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

import devsim
import devsim.python_packages.simple_physics as simple_physics

from src.poisson import extract_depletion_width_numerical

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

    At each voltage, ramps bias using the DD/Poisson solver, extracts the
    depletion width numerically from the carrier profile, and computes
    capacitance via the parallel-plate approximation.

    Reverse bias convention: negative V_range values mean reverse bias on
    the diode. Internally, these are applied as positive cathode bias in
    devsim (positive V on cathode = reverse bias for p+/n- diode).

    Parameters
    ----------
    device_info : dict
        Device info dict (DD equations must be set up via create_dd_device).
    V_range : array_like
        Array of reverse bias voltages (V, should be <= 0 for reverse bias).
    eps_r : float
        Relative permittivity. Default 9.7.
    area : float
        Junction area (cm^2). Default 1.0.

    Returns
    -------
    result : dict
        Dictionary with:
        - "voltages": numpy array of voltages (V), in conventional form
        - "depletion_widths": numpy array of W values (cm)
        - "capacitance": numpy array of C values (F or F/cm^2)
    """
    V_range = np.asarray(V_range, dtype=float)
    device = device_info["device_name"]
    bias_name = simple_physics.GetContactBiasName("cathode")

    depletion_widths = []
    solved_voltages = []

    # Extract W at 0V first (current equilibrium state)
    W0 = extract_depletion_width_numerical(device_info)
    current_V_cathode = 0.0

    for V_target in V_range:
        if abs(V_target) < 1e-12:
            # Use equilibrium W
            depletion_widths.append(W0)
            solved_voltages.append(0.0)
            continue

        # Convert conventional reverse bias to cathode voltage
        # Negative V_target (reverse bias) -> positive cathode voltage
        V_cathode_target = -V_target

        # Ramp cathode in 0.5V steps from current state
        V_step = 0.5
        if V_cathode_target < current_V_cathode:
            V_step = -0.5

        V = current_V_cathode + V_step
        converged = True

        if V_step > 0:
            cond = lambda v: v <= V_cathode_target + 1e-10
        else:
            cond = lambda v: v >= V_cathode_target - 1e-10

        while cond(V):
            devsim.set_parameter(device=device, name=bias_name, value=V)
            try:
                devsim.solve(
                    type="dc",
                    absolute_error=1e10,
                    relative_error=1e-10,
                    maximum_iterations=40,
                )
            except devsim.error:
                try:
                    devsim.solve(
                        type="dc",
                        absolute_error=1e12,
                        relative_error=1e-8,
                        maximum_iterations=100,
                    )
                except devsim.error as e:
                    logger.warning(f"cv_sweep: failed at V_cathode={V:.2f}V: {e}")
                    converged = False
                    break
            V += V_step
            V = round(V, 10)

        if converged:
            W = extract_depletion_width_numerical(device_info)
            depletion_widths.append(W)
            solved_voltages.append(V_target)
            current_V_cathode = V_cathode_target
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
