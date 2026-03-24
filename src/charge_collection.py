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
from src.sic_material import (
    SiC4H_Parameters,
    mobility_caughey_thomas_T,
    srh_lifetime,
)

logger = logging.getLogger(__name__)

# Sentinel for detecting when caller did not pass explicit values
_UNSET = object()


def hecht_cce(
    V,
    d,
    mu_e=_UNSET,
    tau_e=_UNSET,
    mu_p=_UNSET,
    tau_p=_UNSET,
    T=300,
    params=None,
):
    """Two-carrier Hecht equation for charge collection efficiency.

    Assumes uniform electric field E = |V|/d across the active region.
    Valid for fully-depleted detectors in the low-injection regime.

    Parameters
    ----------
    V : float or array_like
        Applied voltage (V). Uses absolute value for reverse bias.
    d : float
        Active region thickness (cm), typically depletion width.
    mu_e : float or None
        Electron mobility (cm^2/Vs). If not provided, computed from T.
    tau_e : float or None
        Electron SRH lifetime (s). If not provided, computed from T.
    mu_p : float or None
        Hole mobility (cm^2/Vs). If not provided, computed from T.
    tau_p : float or None
        Hole SRH lifetime (s). If not provided, computed from T.
    T : float
        Temperature (K). Default 300. Used to compute defaults when
        mu_e/tau_e/mu_p/tau_p are not explicitly provided.
    params : SiC4H_Parameters, optional
        Material parameters. Defaults to SiC4H_Parameters().

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
    if params is None:
        params = SiC4H_Parameters()

    # Use T-dependent defaults when caller doesn't pass explicit values
    if mu_e is _UNSET:
        mu_e = mobility_caughey_thomas_T(params.N_ref_n * 0.001, T, "electron", params)
        # Low-doping limit: use very low N to get mu_max(T)
        # More precisely, use the mu_n_max formula directly
        mu_e = (
            params.mu_n_max * (T / 300.0) ** params.gamma_n
            if T != 300
            else params.mu_n_max
        )
    if tau_e is _UNSET:
        tau_e = srh_lifetime(T, "electron", params)
    if mu_p is _UNSET:
        mu_p = (
            params.mu_p_max * (T / 300.0) ** params.gamma_p
            if T != 300
            else params.mu_p_max
        )
    if tau_p is _UNSET:
        tau_p = srh_lifetime(T, "hole", params)

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
    V,
    d_epi,
    W_func,
    mu_e=_UNSET,
    tau_e=_UNSET,
    mu_p=_UNSET,
    tau_p=_UNSET,
    L_diff_p=7e-4,
    T=300,
    params=None,
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
    if params is None:
        params = SiC4H_Parameters()
    if mu_e is _UNSET:
        mu_e = (
            params.mu_n_max * (T / 300.0) ** params.gamma_n
            if T != 300
            else params.mu_n_max
        )
    if tau_e is _UNSET:
        tau_e = srh_lifetime(T, "electron", params)
    if mu_p is _UNSET:
        mu_p = (
            params.mu_p_max * (T / 300.0) ** params.gamma_p
            if T != 300
            else params.mu_p_max
        )
    if tau_p is _UNSET:
        tau_p = srh_lifetime(T, "hole", params)

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
    # CRITICAL: time_node_model must match original registration in
    # setup_sic_drift_diffusion to preserve transient solve capability.
    devsim.equation(
        device=device,
        region=region,
        name="ElectronContinuityEquation",
        variable_name="Electrons",
        time_node_model="NCharge",
        node_model="ElectronGeneration",
        edge_model="ElectronCurrent",
        variable_update="positive",
    )
    devsim.equation(
        device=device,
        region=region,
        name="HoleContinuityEquation",
        variable_name="Holes",
        time_node_model="PCharge",
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
    I_generated = Q * np.trapezoid(gen_vals, x_nodes)

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
    I_generated = Q * np.trapezoid(gen_values, x_nodes)

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


def cce_vs_fluence(
    fluence_range,
    V_bias=-40.0,
    epi_thickness_cm=10e-4,
    alpha_range_cm=15e-4,
    generation_rate=1e18,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
):
    """Compute CCE vs proton fluence at fixed reverse bias.

    Creates a fresh DD device for each fluence point (fluence-as-temperature
    pattern) using staged device creation so that damaged doping is applied
    before Poisson equilibrium.

    For the internal fluence grid, np.geomspace is recommended since damage
    physics spans orders of magnitude (e.g., 1e10 to 1e14 protons/cm^2).

    Parameters
    ----------
    fluence_range : array_like
        Array of proton fluences (protons/cm^2). Zero fluence returns
        pristine CCE (regression safety).
    V_bias : float
        Fixed reverse bias voltage applied to anode (V, negative).
        Default: -40V.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default: 10 um.
    alpha_range_cm : float
        Alpha particle range in SiC (cm). Default: 15 um.
    generation_rate : float
        Peak generation rate (cm^-3 s^-1). Default: 1e18.
    energy_MeV : float
        Proton energy (MeV). Default: 62.0.
    lifetime_model : str
        "linear" or "logarithmic". Default: "linear".
    damage_params : RadiationDamageParams or None
        Custom damage parameters. Default: RadiationDamageParams().

    Returns
    -------
    result : dict
        Dictionary with:
        - "fluences": numpy array of fluences (protons/cm^2)
        - "cce_values": numpy array of CCE values
        - "V_bias": float, the bias voltage used
        - "energy_MeV": float
        - "lifetime_model": str
    """
    from src.device import apply_damaged_params, create_sic_device
    from src.drift_diffusion import ramp_bias, setup_sic_drift_diffusion
    from src.poisson import setup_poisson, solve_equilibrium
    from src.radiation_damage import compute_damaged_params

    fluence_range = np.asarray(fluence_range, dtype=float)

    # --- Extract pristine parameters for compute_damaged_params ---
    pristine_tau_n = srh_lifetime(300.0, "electron")
    pristine_tau_p = srh_lifetime(300.0, "hole")

    # Extract pristine N_D profile from a reference device (epi-only nodes)
    ref_id = uuid.uuid4().hex[:8]
    ref_name = f"fluence_ref_{ref_id}"
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
        epi_mask_ref = ref_x >= ref_info["junction_pos"]
        pristine_N_D_profile = ref_donors[epi_mask_ref]
    finally:
        try:
            devsim.delete_device(device=ref_name)
        except Exception:
            pass

    # --- Sweep fluence ---
    cce_values = np.zeros(len(fluence_range))

    for i, fluence in enumerate(fluence_range):
        dev_id = uuid.uuid4().hex[:8]
        dev_name = f"fluence_sweep_{dev_id}"

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
            apply_damaged_params(device_info, damaged)

            # Continue staged setup
            setup_poisson(device_info)
            solve_equilibrium(device_info)
            setup_sic_drift_diffusion(device_info)
            device_info["dd_initialized"] = True

            device = device_info["device_name"]
            region = device_info["region_name"]

            # Ramp bias (cathode voltage = -V_bias for reverse bias)
            cathode_V = -V_bias
            ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5)

            # Prepare generation profile
            x_nodes = np.array(
                devsim.get_node_model_values(device=device, region=region, name="x")
            )
            junction_pos = device_info["junction_pos"]
            x_epi = x_nodes - junction_pos

            gen_profile = alpha_generation_profile(x_epi, alpha_range_cm=alpha_range_cm)
            max_profile = np.max(gen_profile)
            if max_profile > 0:
                gen_values = gen_profile * (generation_rate / max_profile)
            else:
                gen_values = np.zeros_like(x_nodes)
            gen_values[x_epi < 0] = 0.0

            # Add generation and solve
            add_generation_to_dd(device_info, gen_values)
            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=1e-10,
                maximum_iterations=40,
            )

            # Extract CCE
            cce = compute_cce_from_dd(device_info, gen_values, contact="cathode")
            cce_values[i] = cce

            logger.info(f"cce_vs_fluence: fluence={fluence:.2e} p/cm^2, CCE={cce:.4f}")

        except Exception as e:
            logger.warning(f"cce_vs_fluence: failed at fluence={fluence:.2e}: {e}")
            cce_values[i] = np.nan
        finally:
            try:
                devsim.delete_device(device=dev_name)
            except Exception:
                pass

    return {
        "fluences": fluence_range,
        "cce_values": cce_values,
        "V_bias": V_bias,
        "energy_MeV": energy_MeV,
        "lifetime_model": lifetime_model,
    }


def cce_vs_bias_at_fluence(
    V_range,
    fluence,
    epi_thickness_cm=10e-4,
    alpha_range_cm=15e-4,
    generation_rate=1e18,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
):
    """Compute CCE vs reverse bias voltage at a fixed proton fluence.

    Creates a single damaged DD device using staged creation, then sweeps
    bias voltage to show how higher reverse bias recovers CCE at a given
    damage level.

    Parameters
    ----------
    V_range : array_like
        Array of voltages (V). Negative = reverse bias on diode.
    fluence : float
        Proton fluence (protons/cm^2).
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default: 10 um.
    alpha_range_cm : float
        Alpha particle range in SiC (cm). Default: 15 um.
    generation_rate : float
        Peak generation rate (cm^-3 s^-1). Default: 1e18.
    energy_MeV : float
        Proton energy (MeV). Default: 62.0.
    lifetime_model : str
        "linear" or "logarithmic". Default: "linear".
    damage_params : RadiationDamageParams or None
        Custom damage parameters. Default: RadiationDamageParams().

    Returns
    -------
    result : dict
        Dictionary with:
        - "voltages": numpy array of applied voltages (V)
        - "cce_values": numpy array of CCE values
        - "fluence": float
        - "energy_MeV": float
        - "lifetime_model": str
    """
    from src.device import apply_damaged_params, create_sic_device
    from src.drift_diffusion import ramp_bias, setup_sic_drift_diffusion
    from src.poisson import setup_poisson, solve_equilibrium
    from src.radiation_damage import compute_damaged_params

    V_range = np.asarray(V_range, dtype=float)

    # --- Extract pristine parameters ---
    pristine_tau_n = srh_lifetime(300.0, "electron")
    pristine_tau_p = srh_lifetime(300.0, "hole")

    # Extract pristine N_D profile from reference device
    ref_id = uuid.uuid4().hex[:8]
    ref_name = f"bias_fluence_ref_{ref_id}"
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
        epi_mask_ref = ref_x >= ref_info["junction_pos"]
        pristine_N_D_profile = ref_donors[epi_mask_ref]
    finally:
        try:
            devsim.delete_device(device=ref_name)
        except Exception:
            pass

    # --- Compute damaged parameters ---
    damaged = compute_damaged_params(
        pristine_tau_n=pristine_tau_n,
        pristine_tau_p=pristine_tau_p,
        N_D_profile=pristine_N_D_profile,
        fluence=fluence,
        energy_MeV=energy_MeV,
        damage_params=damage_params,
        lifetime_model=lifetime_model,
    )

    # --- Staged device creation ---
    dev_id = uuid.uuid4().hex[:8]
    dev_name = f"bias_at_fluence_{dev_id}"

    device_info = create_sic_device(
        device_name=dev_name,
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=2.90e15,
        N_D_bulk=8.50e13,
        L_transition=1.0e-4,
    )

    apply_damaged_params(device_info, damaged)

    setup_poisson(device_info)
    solve_equilibrium(device_info)
    setup_sic_drift_diffusion(device_info)
    device_info["dd_initialized"] = True

    device = device_info["device_name"]
    region = device_info["region_name"]

    # --- Prepare generation profile ---
    x_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    junction_pos = device_info["junction_pos"]
    x_epi = x_nodes - junction_pos

    gen_profile = alpha_generation_profile(x_epi, alpha_range_cm=alpha_range_cm)
    max_profile = np.max(gen_profile)
    if max_profile > 0:
        gen_values = gen_profile * (generation_rate / max_profile)
    else:
        gen_values = np.zeros_like(x_nodes)
    gen_values[x_epi < 0] = 0.0

    Q = 1.602e-19
    I_generated = Q * np.trapezoid(gen_values, x_nodes)

    # --- Sweep bias on single damaged device ---
    # Sort voltages for stable ramping (descending: 0, -10, -20, ...)
    sorted_indices = np.argsort(-V_range)
    sorted_V = V_range[sorted_indices]

    cce_sorted = np.zeros(len(V_range))

    try:
        for idx, V_target in zip(sorted_indices, sorted_V):
            cathode_V = -V_target  # reverse bias convention

            ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5)

            # Add generation and solve
            add_generation_to_dd(device_info, gen_values)
            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=1e-10,
                maximum_iterations=40,
            )

            cce = compute_cce_from_dd(device_info, gen_values, contact="cathode")
            cce_sorted[idx] = cce

            logger.info(
                f"cce_vs_bias_at_fluence: V={V_target:.1f}V, fluence={fluence:.2e}, "
                f"CCE={cce:.4f}"
            )

            # Reset generation for next bias ramp
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
            devsim.delete_device(device=dev_name)
        except Exception:
            pass

    return {
        "voltages": V_range,
        "cce_values": cce_sorted,
        "fluence": fluence,
        "energy_MeV": energy_MeV,
        "lifetime_model": lifetime_model,
    }


def cce_vs_epi_thickness(
    epi_range_cm,
    V_bias=-3.0,
    alpha_range_cm=15e-4,
    generation_rate=1e18,
):
    """Compute CCE vs epitaxial layer thickness at fixed reverse bias.

    Sweeps epi thickness (e.g. 5-20 um) at a single bias voltage to study
    how detector thickness affects charge collection.

    Parameters
    ----------
    epi_range_cm : array_like
        Array of epi thicknesses to sweep (cm).
    V_bias : float
        Fixed reverse bias voltage (V, negative). Default: -3V.
        Use a low bias where partial depletion occurs for the larger
        thicknesses to see the CCE-decreasing-with-thickness effect.
    alpha_range_cm : float
        Alpha particle range in SiC (cm). Default: 15 um.
    generation_rate : float
        Peak generation rate (cm^-3 s^-1). Default: 1e18.

    Returns
    -------
    result : dict
        Dictionary with:
        - "epi_thicknesses": numpy array of epi thicknesses (cm)
        - "cce_values": numpy array of CCE values
        - "V_bias": the bias voltage used

    Notes
    -----
    Expected physics: CCE decreases with epi thickness at low bias
    because thicker epi is harder to fully deplete, and charge generated
    in the neutral (undepleted) region has incomplete collection via
    diffusion only. At high bias (>30V for this doping) all thicknesses
    are fully depleted and CCE ~ 1.0.
    """
    from src.drift_diffusion import (
        create_dd_device,
        ramp_bias,
        extract_contact_current,
    )

    epi_range_cm = np.asarray(epi_range_cm, dtype=float)
    cce_values = np.zeros(len(epi_range_cm))

    Q = 1.602e-19

    for i, epi_thick in enumerate(epi_range_cm):
        dev_id = uuid.uuid4().hex[:8]
        dev_name = f"cce_epi_{dev_id}"

        device_info = create_dd_device(
            device_name=dev_name,
            epi_thickness_cm=epi_thick,
            doping_profile="graded",
            N_D_junction=2.90e15,
            N_D_bulk=8.50e13,
            L_transition=1.0e-4,
        )

        device = device_info["device_name"]
        region = device_info["region_name"]

        try:
            # Get mesh nodes
            x_nodes = np.array(
                devsim.get_node_model_values(device=device, region=region, name="x")
            )

            # Alpha generation profile relative to epi entrance
            junction_pos = device_info["junction_pos"]
            x_epi = x_nodes - junction_pos

            gen_profile = alpha_generation_profile(x_epi, alpha_range_cm=alpha_range_cm)
            max_profile = np.max(gen_profile)
            if max_profile > 0:
                gen_values = gen_profile * (generation_rate / max_profile)
            else:
                gen_values = np.zeros_like(x_nodes)

            # Zero generation in p+ substrate
            gen_values[x_epi < 0] = 0.0

            I_generated = Q * np.trapezoid(gen_values, x_nodes)

            # Ramp to bias (cathode voltage is -V_bias)
            cathode_V = -V_bias
            ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5)

            # Add generation and solve
            add_generation_to_dd(device_info, gen_values)
            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=1e-10,
                maximum_iterations=40,
            )

            # Extract CCE
            cce = compute_cce_from_dd(device_info, gen_values, contact="cathode")
            cce_values[i] = cce

            logger.info(
                f"cce_vs_epi: epi={epi_thick*1e4:.1f}um, V={V_bias:.0f}V, CCE={cce:.4f}"
            )

        finally:
            try:
                devsim.delete_device(device=device)
            except Exception:
                pass

    return {
        "epi_thicknesses": epi_range_cm,
        "cce_values": cce_values,
        "V_bias": V_bias,
    }


def compare_cce_hecht_vs_dd(V_range, epi_thickness_cm=10e-4, **kwargs):
    """Compare DD-computed CCE with Hecht equation benchmark.

    Runs both DD simulation and analytical Hecht equation for the same
    voltage range, documenting agreement and divergence regimes.

    Parameters
    ----------
    V_range : array_like
        Array of voltages (V). Negative = reverse bias.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default: 10 um.
    **kwargs
        Additional keyword arguments passed to cce_vs_bias.

    Returns
    -------
    result : dict
        Dictionary with:
        - "voltages": voltage array
        - "cce_dd": DD-computed CCE values
        - "cce_hecht": Hecht equation CCE (fully depleted approx)
        - "cce_hecht_partial": partial-depletion Hecht CCE (if W_func available)
        - "max_deviation": maximum |DD - Hecht| across all voltages
        - "regime_notes": string documenting regime of validity

    Notes
    -----
    Regime of validity:
    - Hecht assumes uniform E-field (E = V/d), valid when fully depleted
      with uniform doping.
    - DD handles non-uniform E-field from graded doping, diffusion
      collection, and partial depletion naturally.
    - Expected agreement at high bias (>-30V) where field is nearly
      uniform across the depleted region.
    - Expected divergence at low bias where non-uniform field and
      diffusion transport dominate.
    """
    V_range = np.asarray(V_range, dtype=float)

    # DD-computed CCE
    dd_result = cce_vs_bias(V_range, epi_thickness_cm=epi_thickness_cm, **kwargs)
    cce_dd = dd_result["cce_values"]

    # Hecht equation CCE (fully depleted approximation, d = epi thickness)
    cce_hecht = hecht_cce(V_range, d=epi_thickness_cm)

    # Partial depletion Hecht using a simple W(V) model
    # W ~ sqrt(2 * eps * V / (q * N_D)) with calibrated N_D
    # For graded doping, use an effective N_D ~ geometric mean
    N_D_eff = np.sqrt(2.90e15 * 8.50e13)  # ~ 1.57e14
    eps_sic = 9.66 * 8.854e-14  # F/cm
    Q_val = 1.602e-19

    def W_func(v):
        v_abs = abs(v)
        if v_abs < 1e-12:
            return 0.0
        # Add built-in potential (~2.7V for SiC p+n)
        V_total = v_abs + 2.7
        W = np.sqrt(2 * eps_sic * V_total / (Q_val * N_D_eff))
        return min(W, epi_thickness_cm)

    cce_hecht_partial = hecht_cce_partial_depletion(V_range, epi_thickness_cm, W_func)

    # Compute maximum deviation
    max_dev = float(np.max(np.abs(cce_dd - np.asarray(cce_hecht, dtype=float))))

    # Compute agreement metrics (R-squared, RMSE, relative errors)
    from src.validation import compute_agreement_metrics

    metrics_vs_hecht = compute_agreement_metrics(cce_dd, cce_hecht)
    metrics_vs_partial = compute_agreement_metrics(cce_dd, cce_hecht_partial)

    # Regime notes
    regime_notes = (
        "Hecht equation assumes uniform E-field (E=V/d), valid for fully "
        "depleted detectors with uniform doping. DD solver handles non-uniform "
        "E-field from graded doping, diffusion collection, and partial depletion. "
        "Agreement expected at high reverse bias (>30V) where field is nearly "
        "uniform. Divergence expected at low bias where diffusion transport "
        "and non-uniform field dominate."
    )

    return {
        "voltages": V_range,
        "cce_dd": cce_dd,
        "cce_hecht": np.asarray(cce_hecht, dtype=float),
        "cce_hecht_partial": np.asarray(cce_hecht_partial, dtype=float),
        "max_deviation": max_dev,
        "agreement_metrics_hecht": metrics_vs_hecht,
        "agreement_metrics_partial": metrics_vs_partial,
        "regime_notes": regime_notes,
    }
