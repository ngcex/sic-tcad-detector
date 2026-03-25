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
import uuid

import matplotlib.pyplot as plt
import numpy as np

import devsim
import devsim.python_packages.simple_physics as simple_physics

from src.device import apply_damaged_params, create_sic_device
from src.drift_diffusion import setup_sic_drift_diffusion
from src.poisson import (
    extract_depletion_width_numerical,
    setup_poisson,
    solve_equilibrium,
)
from src.radiation_damage import compute_damaged_params, compute_phi_crit
from src.sic_material import srh_lifetime

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


def cv_at_fluence(
    fluence,
    V_range,
    area=1.0,
    epi_thickness_cm=10e-4,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
    phi_crit_threshold=0.90,
):
    """Compute C-V curve for a device at a given proton fluence.

    Creates a fresh DD device with radiation-damaged parameters using the
    staged device creation pattern (fluence-as-temperature), runs a C-V
    sweep, and cleans up the device.

    Parameters
    ----------
    fluence : float
        Proton fluence (protons/cm^2).
    V_range : array_like
        Array of reverse bias voltages (V, should be <= 0 for reverse bias).
    area : float
        Junction area (cm^2). Default 1.0.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default 10 um.
    energy_MeV : float
        Proton energy (MeV). Default 62.0.
    lifetime_model : str
        "linear" or "logarithmic". Default "linear".
    damage_params : RadiationDamageParams or None
        Custom damage parameters. Default: RadiationDamageParams().
    phi_crit_threshold : float
        Fraction of Phi_crit at which to warn. Default 0.90.

    Returns
    -------
    dict or None
        C-V result dict with keys: voltages, depletion_widths, capacitance,
        fluence. Returns None if fluence >= Phi_crit (solver would diverge).
    """
    V_range = np.asarray(V_range, dtype=float)

    # Extract pristine lifetimes
    pristine_tau_n = srh_lifetime(300.0, "electron")
    pristine_tau_p = srh_lifetime(300.0, "hole")

    # Create a reference device to extract pristine N_D profile
    ref_id = uuid.uuid4().hex[:8]
    ref_name = f"cv_ref_{ref_id}"
    ref_info = create_sic_device(
        device_name=ref_name,
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=2.90e15,
        N_D_bulk=8.50e13,
        L_transition=1.0e-4,
    )
    try:
        ref_x = np.array(
            devsim.get_node_model_values(
                device=ref_name, region=ref_info["region_name"], name="x"
            )
        )
        ref_donors = np.array(
            devsim.get_node_model_values(
                device=ref_name, region=ref_info["region_name"], name="Donors"
            )
        )
        epi_mask = ref_x >= ref_info["junction_pos"]
        pristine_N_D_profile = ref_donors[epi_mask]
    finally:
        try:
            devsim.delete_device(device=ref_name)
        except Exception:
            pass

    # Check Phi_crit
    phi_crit_info = compute_phi_crit(
        pristine_N_D_profile,
        energy_MeV=energy_MeV,
    )
    phi_crit = phi_crit_info["phi_crit_proton"]

    if fluence >= phi_crit:
        logger.error(
            "Fluence %.3e >= Phi_crit %.3e protons/cm^2: full compensation, "
            "cannot solve. Returning None.",
            fluence,
            phi_crit,
        )
        return None

    if fluence >= phi_crit_threshold * phi_crit:
        logger.warning(
            "Fluence %.3e is >= %.0f%% of Phi_crit (%.3e): approaching "
            "full carrier compensation.",
            fluence,
            phi_crit_threshold * 100,
            phi_crit,
        )

    # Create device for this fluence point
    dev_id = uuid.uuid4().hex[:8]
    dev_name = f"cv_fluence_{dev_id}"

    try:
        # Compute damaged parameters
        damaged = compute_damaged_params(
            pristine_tau_n=pristine_tau_n,
            pristine_tau_p=pristine_tau_p,
            N_D_profile=pristine_N_D_profile,
            fluence=fluence,
            energy_MeV=energy_MeV,
            damage_params=damage_params,
            lifetime_model=lifetime_model,
        )

        # Staged device creation
        device_info = create_sic_device(
            device_name=dev_name,
            epi_thickness_cm=epi_thickness_cm,
            doping_profile="graded",
            N_D_junction=2.90e15,
            N_D_bulk=8.50e13,
            L_transition=1.0e-4,
        )

        # Apply damage BEFORE Poisson setup
        if fluence > 0:
            apply_damaged_params(device_info, damaged)

        # Continue staged setup
        setup_poisson(device_info)
        solve_equilibrium(device_info)
        setup_sic_drift_diffusion(device_info)

        # Run C-V sweep
        cv_result = cv_sweep(device_info, V_range, area=area)
        cv_result["fluence"] = fluence

        logger.info(
            "cv_at_fluence: fluence=%.3e, %d voltage points solved",
            fluence,
            len(cv_result["voltages"]),
        )

        return cv_result

    except Exception as e:
        logger.error("cv_at_fluence failed at fluence=%.3e: %s", fluence, e)
        raise
    finally:
        try:
            devsim.delete_device(device=dev_name)
        except Exception:
            pass


def plot_cv_evolution(cv_results, fluences, ax=None, title=None):
    """Overlay C-V curves at different fluence levels.

    Uses a viridis colormap gradient to distinguish fluence levels, with
    lower fluence in darker colors and higher fluence in lighter colors.

    Parameters
    ----------
    cv_results : list of dict
        List of C-V result dicts (from cv_at_fluence or cv_sweep).
        None entries are skipped (above Phi_crit).
    fluences : array_like
        Fluence values corresponding to each cv_result (protons/cm^2).
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    title : str, optional
        Plot title. Default: "C-V Evolution with Fluence".

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    fluences = np.asarray(fluences, dtype=float)
    cmap = plt.cm.viridis

    # Normalize fluence range for colormap
    valid_fluences = [f for f, r in zip(fluences, cv_results) if r is not None]
    if len(valid_fluences) == 0:
        logger.warning("plot_cv_evolution: no valid C-V results to plot")
        return ax

    f_min = min(valid_fluences)
    f_max = max(valid_fluences)
    if f_max > f_min:
        norm = plt.Normalize(vmin=f_min, vmax=f_max)
    else:
        norm = plt.Normalize(vmin=0, vmax=max(f_max, 1.0))

    for fluence, cv_result in zip(fluences, cv_results):
        if cv_result is None:
            continue
        color = cmap(norm(fluence))
        if fluence == 0:
            label = "Pristine"
        else:
            label = f"{fluence:.1e} p/cm$^2$"
        ax.plot(
            cv_result["voltages"], cv_result["capacitance"], color=color, label=label
        )

    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("Capacitance (F/cm$^2$)")
    ax.set_title(title or "C-V Evolution with Fluence")
    ax.legend(fontsize=8)

    return ax
