"""Charge collection efficiency (CCE) computation utilities.

Provides:
- Hecht equation (two-carrier form) for analytical CCE benchmark
- Extended Hecht for partially-depleted detectors with diffusion
- CCE extraction from DD simulation current ratio
- DD-based CCE computation with radiation generation injection
- CCE vs bias sweep with Hecht comparison

All units in CGS (cm, cm^-3, V, s) per devsim convention.

References:
    - Hecht, Z. Physik 77, 235 (1932)
    - Sze & Ng, "Physics of Semiconductor Devices", 3rd ed.
    - Knoll, "Radiation Detection and Measurement", 4th ed.
"""

import logging
import uuid

import devsim
from devsim.python_packages.model_create import (
    CreateNodeModel,
    CreateNodeModelDerivative,
)
import numpy as np

from src.generation_profiles import alpha_generation_profile

logger = logging.getLogger(__name__)


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


# ---------------------------------------------------------------------------
# DD-based CCE computation (Plan 02)
# ---------------------------------------------------------------------------


def add_generation_to_dd(device_info, generation_values):
    """Add spatially-varying radiation generation rate to DD equations.

    Modifies the ElectronGeneration and HoleGeneration node models to
    include a RadGenRate source term alongside the existing SRH
    recombination.

    Parameters
    ----------
    device_info : dict
        Device info dict from create_dd_device (must have dd_initialized=True).
    generation_values : array_like
        Array of G(x) at each mesh node (cm^-3 s^-1). Length must match
        the number of mesh nodes.

    Notes
    -----
    Sign convention for continuity equations in devsim:
    - ElectronGeneration is the node_model in ElectronContinuityEquation.
      Positive ElectronGeneration = electron source (creation).
    - HoleGeneration is the node_model in HoleContinuityEquation.
      Positive HoleGeneration = hole source (creation), but the
      convention uses opposite sign: HoleGeneration = +q*USRH for
      recombination (loss), so generation requires -q*G.

    Verification: after adding generation, carrier densities should
    INCREASE relative to equilibrium values.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]

    # Step 1: Create RadGenRate node model initialized to "0"
    # Then set actual values via set_node_values
    try:
        CreateNodeModel(device, region, "RadGenRate", "0")
    except devsim.error:
        # Model already exists (e.g., updating generation), just overwrite values
        pass

    devsim.set_node_values(
        device=device,
        region=region,
        name="RadGenRate",
        values=list(np.asarray(generation_values, dtype=float)),
    )

    # Step 2: Update ElectronGeneration to include radiation generation
    # Original: Gn = "-ElectronCharge * USRH"
    # Updated: Gn = "-ElectronCharge * USRH + ElectronCharge * RadGenRate"
    # The +ElectronCharge * RadGenRate term creates electrons (positive = source)
    Gn = "-ElectronCharge * USRH + ElectronCharge * RadGenRate"
    CreateNodeModel(device, region, "ElectronGeneration", Gn)
    # RadGenRate has no carrier dependence, so only USRH derivatives matter
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn, var)

    # Step 3: Update HoleGeneration to include radiation generation
    # Original: Gp = "+ElectronCharge * USRH"
    # Updated: Gp = "+ElectronCharge * USRH - ElectronCharge * RadGenRate"
    # The -ElectronCharge * RadGenRate term creates holes
    # (HoleContinuityEquation uses opposite sign convention)
    Gp = "+ElectronCharge * USRH - ElectronCharge * RadGenRate"
    CreateNodeModel(device, region, "HoleGeneration", Gp)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "HoleGeneration", Gp, var)

    # Step 4: Re-register continuity equations with updated node models
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

    device_info["generation_added"] = True
    logger.debug("Added radiation generation to DD equations")


def compute_cce_from_dd(device_info, generation_values, contact="cathode"):
    """Extract CCE from a solved DD device with radiation generation.

    CCE = |I_collected| / I_generated, where I_generated is the total
    generation current integrated over the device volume.

    Parameters
    ----------
    device_info : dict
        Device info dict with DD solved and generation added.
    generation_values : array_like
        Array of G(x) at each mesh node (cm^-3 s^-1).
    contact : str
        Contact at which to extract collected current.

    Returns
    -------
    cce : float
        Charge collection efficiency, clipped to [0, 1].
    """
    from src.drift_diffusion import extract_contact_current

    Q = 1.602e-19  # C

    device = device_info["device_name"]
    region = device_info["region_name"]

    # Get mesh node positions
    x_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    gen_vals = np.asarray(generation_values, dtype=float)

    # Total generated current density: Q * integral(G(x) dx) [A/cm^2 for 1D]
    I_generated = Q * np.trapz(gen_vals, x_nodes)

    # Collected current at contact
    I_collected = abs(extract_contact_current(device_info, contact))

    if I_generated <= 0:
        return 0.0

    cce = I_collected / I_generated
    return float(np.clip(cce, 0.0, 1.0))


def cce_vs_bias(
    V_range,
    epi_thickness_cm=10e-4,
    alpha_range_cm=15e-4,
    generation_rate=1e18,
    device_kwargs=None,
):
    """Compute CCE vs reverse bias voltage for alpha particle irradiation.

    Creates a DD device with calibrated graded doping, sweeps reverse bias,
    and computes CCE at each voltage point.

    Parameters
    ----------
    V_range : array_like
        Array of voltages (V). Negative = reverse bias on diode.
        Convention: V applied to anode; reverse bias is negative.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default: 10 um.
    alpha_range_cm : float
        Alpha particle range in SiC (cm). Default: 15 um.
    generation_rate : float
        Peak generation rate (cm^-3 s^-1). Default: 1e18.
        Low injection since N_D ~ 1e13-1e15.
    device_kwargs : dict or None
        Additional keyword arguments for create_dd_device.

    Returns
    -------
    result : dict
        Dictionary with:
        - "voltages": numpy array of applied voltages (V)
        - "cce_values": numpy array of CCE values
        - "I_collected": numpy array of collected currents (A/cm^2)
        - "I_generated": total generated current (A/cm^2)
    """
    from src.drift_diffusion import create_dd_device, ramp_bias

    V_range = np.asarray(V_range, dtype=float)

    # Unique device name to avoid devsim name conflicts
    dev_id = uuid.uuid4().hex[:8]
    dev_name = f"cce_sweep_{dev_id}"

    # Calibrated graded doping from Phase 2
    kwargs = dict(
        device_name=dev_name,
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=2.90e15,
        N_D_bulk=8.50e13,
        L_transition=1.0e-4,
    )
    if device_kwargs:
        kwargs.update(device_kwargs)

    device_info = create_dd_device(**kwargs)

    device = device_info["device_name"]
    region = device_info["region_name"]

    # Get mesh node positions for generation profile
    x_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )

    # Generate alpha particle profile at mesh nodes
    # x_nodes are absolute positions; generation starts at detector entrance
    # For p+/n- diode, radiation enters from cathode side (n- epi, right side)
    # so generation is relative to epi entrance
    junction_pos = device_info["junction_pos"]
    x_epi = x_nodes - junction_pos  # relative to junction (epi start)

    # Alpha generation profile (pairs/cm per alpha)
    gen_profile = alpha_generation_profile(x_epi, alpha_range_cm=alpha_range_cm)

    # Scale to generation rate (cm^-3 s^-1)
    # gen_profile is in pairs/cm; to get cm^-3 s^-1 we need to scale
    # such that the peak is approximately generation_rate
    max_profile = np.max(gen_profile)
    if max_profile > 0:
        gen_values = gen_profile * (generation_rate / max_profile)
    else:
        gen_values = np.zeros_like(x_nodes)

    # Zero out generation in p+ substrate (x < junction_pos)
    gen_values[x_epi < 0] = 0.0

    Q = 1.602e-19
    I_generated = Q * np.trapz(gen_values, x_nodes)

    cce_values = []
    I_collected_list = []

    # Sort voltages for stable ramping (0V first, then increasingly negative)
    # We need to track the ramp order
    sorted_indices = np.argsort(-V_range)  # descending: 0, -10, -20, ...
    sorted_V = V_range[sorted_indices]

    cce_sorted = np.zeros(len(V_range))
    I_sorted = np.zeros(len(V_range))

    try:
        for idx, V_target in zip(sorted_indices, sorted_V):
            # For anode contact: negative V = reverse bias
            # Ramp to bias point FIRST (without generation)
            # The cathode bias is -V_target for reverse bias convention
            # Actually, for iv_sweep convention: V is applied to anode
            # Reverse bias on diode = negative anode voltage
            # But devsim convention from calibrate_graded_doping:
            # cathode_voltages = [-v for v in voltages] where v is negative
            # So reverse bias -40V on diode = cathode at +40V
            cathode_V = -V_target  # e.g., V_target=-40 -> cathode_V=+40

            ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5)

            # Add generation and re-solve
            add_generation_to_dd(device_info, gen_values)

            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=1e-10,
                maximum_iterations=40,
            )

            # Extract CCE
            cce = compute_cce_from_dd(device_info, gen_values, contact="cathode")
            from src.drift_diffusion import extract_contact_current

            I_coll = abs(extract_contact_current(device_info, "cathode"))

            cce_sorted[idx] = cce
            I_sorted[idx] = I_coll

            logger.info(
                f"cce_vs_bias: V={V_target:.1f}V, CCE={cce:.4f}, "
                f"I_coll={I_coll:.4e} A/cm^2"
            )

            # Remove generation for next bias ramp (reset to SRH-only)
            # by setting generation values to zero
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

    return {
        "voltages": V_range,
        "cce_values": cce_sorted,
        "I_collected": I_sorted,
        "I_generated": I_generated,
    }
