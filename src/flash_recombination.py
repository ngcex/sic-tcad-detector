"""Auger recombination and high-injection continuation solver for FLASH dose rates.

At conventional dose rates, Shockley-Read-Hall (SRH) recombination dominates
carrier loss in 4H-SiC detectors. Under FLASH dose rates (20-230 Gy/s),
carrier densities can become high enough that Auger recombination -- which
scales as n^3 at high injection -- becomes significant. This additional
recombination channel reduces the charge collection efficiency (CCE),
contributing to the observed dose-rate-dependent response in SiC dosimeters.

The Auger recombination rate is:
    R_Auger = (C_n * n + C_p * p) * (n * p - n_i^2)

where C_n and C_p are the Auger coefficients for electron- and hole-initiated
processes, respectively. For 4H-SiC, C_n = 5e-31 cm^6/s and C_p = 2e-31 cm^6/s
(Ioffe NSM Archive).

This module provides:
- add_auger_recombination: extend DD equations with Auger model and Jacobian
- solve_with_continuation: generation-rate ramping for convergence at high injection

All units CGS (cm, cm^-3, A/cm^2, V) per devsim convention.
"""

import logging
import uuid

import devsim
from devsim.python_packages.model_create import (
    CreateNodeModel,
    CreateNodeModelDerivative,
)
import numpy as np

from src.charge_collection import add_generation_to_dd, compute_cce_from_dd
from src.generation_profiles import proton_generation_profile

logger = logging.getLogger(__name__)


def add_auger_recombination(device_info):
    """Add Auger recombination to the drift-diffusion continuity equations.

    Must be called after setup_sic_drift_diffusion(). Can be called before
    or after add_generation_to_dd() -- if RadGenRate does not yet exist,
    a zero-valued model is created automatically.

    Parameters
    ----------
    device_info : dict
        Device info dict from create_dd_device() with dd_initialized=True.

    Notes
    -----
    Creates UAuger node model and updates ElectronGeneration / HoleGeneration
    to include combined SRH + Auger recombination. All Jacobian derivatives
    w.r.t. Electrons and Holes are registered for Newton solver convergence.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    params = device_info["params"]

    # --- Step 1: Set Auger coefficients as region parameters ---
    C_n = params.C_n
    C_p = params.C_p
    devsim.set_parameter(device=device, region=region, name="C_n", value=C_n)
    devsim.set_parameter(device=device, region=region, name="C_p", value=C_p)
    logger.info(f"Set Auger coefficients: C_n={C_n:.2e}, C_p={C_p:.2e} cm^6/s")

    # --- Step 2: Create UAuger node model ---
    # R_Auger = (C_n * n + C_p * p) * (n * p - n_i^2)
    UAuger = "(C_n * Electrons + C_p * Holes) * (Electrons * Holes - n_i^2)"
    CreateNodeModel(device, region, "UAuger", UAuger)

    # Jacobian derivatives (CRITICAL for Newton convergence)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "UAuger", UAuger, var)

    # --- Step 3: Ensure RadGenRate model exists ---
    try:
        devsim.get_node_model_values(device=device, region=region, name="RadGenRate")
    except devsim.error:
        # RadGenRate not yet created; create a zero placeholder
        CreateNodeModel(device, region, "RadGenRate", "0")
        logger.debug("Created zero RadGenRate placeholder for Auger setup")

    # --- Step 4: Update ElectronGeneration to include Auger ---
    # Gn = -q * (USRH + UAuger) + q * RadGenRate
    Gn = "-ElectronCharge * (USRH + UAuger) + ElectronCharge * RadGenRate"
    CreateNodeModel(device, region, "ElectronGeneration", Gn)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn, var)

    # --- Step 5: Update HoleGeneration to include Auger ---
    # Gp = +q * (USRH + UAuger) - q * RadGenRate
    Gp = "+ElectronCharge * (USRH + UAuger) - ElectronCharge * RadGenRate"
    CreateNodeModel(device, region, "HoleGeneration", Gp)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "HoleGeneration", Gp, var)

    # --- Step 6: Re-register continuity equations with updated models ---
    devsim.equation(
        device=device,
        region=region,
        name="ElectronContinuityEquation",
        variable_name="Electrons",
        node_model="ElectronGeneration",
        edge_model="ElectronCurrent",
        variable_update="positive",
    )
    devsim.equation(
        device=device,
        region=region,
        name="HoleContinuityEquation",
        variable_name="Holes",
        node_model="HoleGeneration",
        edge_model="HoleCurrent",
        variable_update="positive",
    )

    device_info["auger_enabled"] = True
    logger.info("Auger recombination added to DD equations")


def solve_with_continuation(device_info, generation_values, n_steps=5):
    """Solve DD equations with generation-rate continuation for convergence.

    Ramps the radiation generation rate from 10% to 100% of the target
    in equal fractional steps. At each step, the solver is invoked with
    relaxed tolerances. If a step fails, the remaining ramp is bisected
    to find a convergent path.

    Parameters
    ----------
    device_info : dict
        Device info dict with DD and (optionally) Auger already set up.
    generation_values : array_like
        Target generation rate profile at mesh nodes (cm^-3 s^-1).
    n_steps : int
        Number of continuation steps (default 5: 20%, 40%, ..., 100%).

    Returns
    -------
    converged : bool
        True if the solver converged at 100% generation rate.
    """
    generation_values = np.asarray(generation_values, dtype=float)

    # Build list of fractional steps from 1/n_steps to 1.0
    fractions = np.linspace(1.0 / n_steps, 1.0, n_steps)

    # Track current fraction achieved
    current_frac = 0.0
    max_bisections = 3  # limit bisection depth

    for i, target_frac in enumerate(fractions):
        # Scale generation to this fraction
        scaled_gen = generation_values * target_frac

        add_generation_to_dd(device_info, scaled_gen)

        try:
            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=1e-10,
                maximum_iterations=60,
            )
            current_frac = target_frac
            logger.info(
                f"Continuation step {i+1}/{n_steps}: "
                f"frac={target_frac:.2f} CONVERGED"
            )
        except devsim.error:
            # Bisection: try intermediate fractions between current and target
            logger.warning(
                f"Continuation step {i+1}/{n_steps}: "
                f"frac={target_frac:.2f} FAILED, attempting bisection"
            )
            converged_at_target = False
            lo = current_frac
            hi = target_frac

            for bisect_iter in range(max_bisections):
                mid = (lo + hi) / 2.0
                scaled_mid = generation_values * mid
                add_generation_to_dd(device_info, scaled_mid)

                try:
                    devsim.solve(
                        type="dc",
                        absolute_error=1e10,
                        relative_error=1e-10,
                        maximum_iterations=60,
                    )
                    lo = mid
                    current_frac = mid
                    logger.info(
                        f"  Bisection {bisect_iter+1}: frac={mid:.4f} CONVERGED"
                    )
                except devsim.error:
                    hi = mid
                    logger.warning(
                        f"  Bisection {bisect_iter+1}: frac={mid:.4f} FAILED"
                    )

            # After bisection, try the original target again
            if current_frac < target_frac:
                scaled_gen = generation_values * target_frac
                add_generation_to_dd(device_info, scaled_gen)
                try:
                    devsim.solve(
                        type="dc",
                        absolute_error=1e10,
                        relative_error=1e-10,
                        maximum_iterations=60,
                    )
                    current_frac = target_frac
                    converged_at_target = True
                    logger.info(
                        f"  Retry at frac={target_frac:.2f} CONVERGED after bisection"
                    )
                except devsim.error:
                    logger.error(
                        f"  Failed to reach frac={target_frac:.2f} "
                        f"even after bisection (stuck at {current_frac:.4f})"
                    )
                    return False
            else:
                converged_at_target = True

            if not converged_at_target:
                return False

    logger.info(f"Continuation complete: converged at 100% generation rate")
    return True


def cce_vs_dose_rate(
    dose_rates_Gy_s,
    V_bias=-30.0,
    epi_thickness_cm=10e-4,
    E_MeV=62,
    n_continuation_steps=5,
):
    """Compute CCE vs dose rate across the FLASH range.

    Sweeps dose rates (typically 20-230 Gy/s) at fixed reference conditions
    to reveal whether Auger recombination degrades charge collection at
    FLASH dose rates.

    Parameters
    ----------
    dose_rates_Gy_s : array_like
        Dose rates to sweep (Gy/s). Sorted ascending internally.
    V_bias : float
        Reverse bias voltage (V, negative). Default: -30.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default: 10 um.
    E_MeV : float
        Proton kinetic energy (MeV). Default: 62.
    n_continuation_steps : int
        Steps for generation-rate continuation ramp. Default: 5.

    Returns
    -------
    result : dict
        Dictionary with:
        - "dose_rates": numpy array of dose rates (Gy/s)
        - "cce_values": numpy array of CCE values
        - "cce_no_auger_ref": float, CCE at lowest dose rate without Auger
        - "V_bias": bias voltage used
        - "epi_thickness_cm": epi thickness used
        - "E_MeV": proton energy used
    """
    from src.drift_diffusion import create_dd_device, ramp_bias

    dose_rates = np.asarray(dose_rates_Gy_s, dtype=float)
    sorted_idx = np.argsort(dose_rates)
    dose_rates_sorted = dose_rates[sorted_idx]

    cce_sorted = np.zeros(len(dose_rates))

    # --- Device with Auger ---
    dev_id = uuid.uuid4().hex[:8]
    device_info = create_dd_device(
        device_name=f"flash_sweep_{dev_id}",
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=2.90e15,
        N_D_bulk=8.50e13,
        L_transition=1.0e-4,
    )
    device = device_info["device_name"]
    region = device_info["region_name"]

    try:
        # Add Auger recombination
        add_auger_recombination(device_info)

        # Ramp to bias (cathode voltage is -V_bias for reverse bias)
        cathode_V = -V_bias
        ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5)

        # Get mesh nodes
        x_nodes = np.array(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        junction_pos = device_info["junction_pos"]

        # Loop over dose rates (ascending for continuation stability)
        for i, dose_rate in enumerate(dose_rates_sorted):
            gen_values = proton_generation_profile(x_nodes, E_MeV, dose_rate)

            # Zero generation in p+ substrate (x < junction_pos)
            gen_values[x_nodes < junction_pos] = 0.0

            converged = solve_with_continuation(
                device_info, gen_values, n_steps=n_continuation_steps
            )

            if converged:
                cce = compute_cce_from_dd(device_info, gen_values, contact="cathode")
            else:
                logger.warning(f"Continuation failed at dose_rate={dose_rate:.1f} Gy/s")
                cce = np.nan

            cce_sorted[i] = cce
            logger.info(f"cce_vs_dose_rate: {dose_rate:.0f} Gy/s -> CCE={cce:.6f}")

            # Reset generation to zero for next iteration
            zero_gen = np.zeros_like(gen_values)
            add_generation_to_dd(device_info, zero_gen)
            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=1e-10,
                maximum_iterations=40,
            )

    finally:
        try:
            devsim.delete_device(device=device)
        except Exception:
            pass

    # --- No-Auger reference at lowest dose rate ---
    dev_id2 = uuid.uuid4().hex[:8]
    device_info2 = create_dd_device(
        device_name=f"flash_noauger_{dev_id2}",
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=2.90e15,
        N_D_bulk=8.50e13,
        L_transition=1.0e-4,
    )
    device2 = device_info2["device_name"]
    region2 = device_info2["region_name"]

    try:
        # No Auger -- just SRH
        cathode_V = -V_bias
        ramp_bias(device_info2, cathode_V, contact="cathode", V_step=0.5)

        x_nodes2 = np.array(
            devsim.get_node_model_values(device=device2, region=region2, name="x")
        )
        junction_pos2 = device_info2["junction_pos"]

        lowest_rate = dose_rates_sorted[0]
        gen_ref = proton_generation_profile(x_nodes2, E_MeV, lowest_rate)
        gen_ref[x_nodes2 < junction_pos2] = 0.0

        add_generation_to_dd(device_info2, gen_ref)
        devsim.solve(
            type="dc",
            absolute_error=1e10,
            relative_error=1e-10,
            maximum_iterations=40,
        )
        cce_no_auger = compute_cce_from_dd(device_info2, gen_ref, contact="cathode")
        logger.info(
            f"No-Auger reference CCE at {lowest_rate:.0f} Gy/s: {cce_no_auger:.6f}"
        )

    finally:
        try:
            devsim.delete_device(device=device2)
        except Exception:
            pass

    # Unsort back to original order
    cce_original_order = np.zeros(len(dose_rates))
    for i, idx in enumerate(sorted_idx):
        cce_original_order[idx] = cce_sorted[i]

    return {
        "dose_rates": dose_rates,
        "cce_values": cce_original_order,
        "cce_no_auger_ref": float(cce_no_auger),
        "V_bias": V_bias,
        "epi_thickness_cm": epi_thickness_cm,
        "E_MeV": E_MeV,
    }
