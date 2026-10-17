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

from etna.core.solver_tol import DC_REL_ERROR
import logging
import uuid

import devsim
from devsim.python_packages.model_create import (
    CreateNodeModel,
    CreateNodeModelDerivative,
)
import numpy as np

from etna.core.generation_profiles import (
    alpha_generation_profile,
    confine_generation_to_epi,
)
from etna.core.sic_material import (
    SiC4H_Parameters,
    mobility_caughey_thomas_T,
    srh_lifetime,
)

logger = logging.getLogger(__name__)


# graded-doping defaults of the fluence sweeps = cce_vs_bias = the default stack
from etna.core.device import DEFAULT_STACK as _STACK  # noqa: E402


def full_compensation_reason(damaged):
    """Reason string if carrier removal fully compensates the epi donors, else "".

    ``damaged`` is the dict of ``compute_damaged_params`` / ``compute_annealed_params``.
    Carrier removal there is N_D(x) - eta_removal * Phi_neq floored at 0 (annealed:
    eta_removal * (1 - f_Z12)), so a zero minimum means eta * Phi_neq >= min N_D over
    the epi: the n- epi has no net donors left somewhere. The doping-based device model
    (ionised N_D, ohmic contact on the n- epi) is ill-posed there -- devsim either fails or
    converges to non-physical values -- so callers return NaN with this reason instead
    of solving. Uses the damage model's own constants (kappa(E), eta_removal, anneal).
    """
    nd = np.asarray(damaged.get("N_D_profile", []), dtype=float)
    if nd.size == 0 or float(np.min(nd)) > 0.0:
        return ""
    return (
        "full donor compensation: eta_removal * Phi_neq >= min N_D over the epi "
        f"(Phi_neq = {damaged.get('fluence_neq', float('nan')):.3g} cm^-2); "
        "outside model validity, not solved"
    )


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
    Two-carrier Hecht for charge generated UNIFORMLY through the active region,
    averaged over the generation depth x in [0, d]:

        CCE = (le/d) - (le/d)^2*(1-exp(-d/le))
            + (lh/d) - (lh/d)^2*(1-exp(-d/lh))

    where le = mu_e*tau_e*|V|/d (electron drift length), lh similarly for holes.
    This is the depth-average of the single-pair induced charge
        q(x)/q0 = (le/d)*(1-exp(-(d-x)/le)) + (lh/d)*(1-exp(-x/lh)),
    in which each electron-hole pair contributes ONCE. In the long-drift limit
    (le, lh >> d) each carrier term -> 1/2, so CCE -> 1.0 without clipping.

    The earlier form (le/d)*(1-exp(-d/le)) + (lh/d)*(1-exp(-d/lh)) double-counted
    each pair (each term -> 1, CCE -> 2.0, then clipped). See audit C1.

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

    # Depth-averaged single-pair Hecht: each carrier term is intrinsically <= 1/2,
    # so the e+h sum is bounded by 1 without clipping (audit C1).
    # Avoid division by zero at V=0 (lambda=0).
    with np.errstate(divide="ignore", invalid="ignore"):
        re = lambda_e / d
        rh = lambda_h / d
        cce_e = np.where(lambda_e > 0, re - re**2 * (1.0 - np.exp(-d / lambda_e)), 0.0)
        cce_h = np.where(lambda_h > 0, rh - rh**2 * (1.0 - np.exp(-d / lambda_h)), 0.0)

    cce = cce_e + cce_h
    # Defensive clip only (formula is already bounded to [0, 1]).
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
    from etna.core.drift_diffusion import extract_steady_state_current

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
    # Balanced form: the direct majority-carrier contact current has a round-off floor
    # and is unusable with an n+ cathode layer (see extract_steady_state_current).
    I_collected = abs(extract_steady_state_current(device_info, contact))

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
        - "voltages": numpy array of applied voltages actually reached (V) —
          shorter than the requested V_range if the DD solver hit its
          convergence wall partway through the sweep; see "truncated". The wall
          is a model-validity limit (ohmic cathode directly on the n- epi, no n+
          substrate; devsim stops converging well past reach-through at V_fd),
          not physical punch-through.
        - "cce_values": numpy array of CCE values, same length as "voltages"
        - "I_collected": numpy array of collected currents (A/cm^2), same
          length as "voltages"
        - "I_generated": total generated current (A/cm^2)
        - "truncated": True if the solver failed to converge before reaching
          the end of V_range (a partial result is returned instead of
          raising); False if the full requested sweep converged
        - "requested_v_stop": the most extreme (deepest reverse bias) voltage
          in V_range, for comparison against the last entry of "voltages"
          when "truncated" is True
    """
    from etna.core.drift_diffusion import create_dd_device, ramp_bias

    V_range = np.asarray(V_range, dtype=float)

    # Unique device name to avoid devsim name conflicts
    dev_id = uuid.uuid4().hex[:8]
    dev_name = f"cce_sweep_{dev_id}"

    # graded doping of the default stack (etna.core.device.DEFAULT_STACK); was a
    # fourth set 2.90e15 / 8.50e13 / 1.0e-4 before the n+ flip. Geometry (p+, n+) comes
    # from device_kwargs (the facades pass the full DeviceConfig stack).
    from etna.core.device import DEFAULT_STACK

    kwargs = dict(
        device_name=dev_name,
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=DEFAULT_STACK["N_D_junction"],
        N_D_bulk=DEFAULT_STACK["N_D_bulk"],
        L_transition=DEFAULT_STACK["L_transition"],
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
    gen_values = confine_generation_to_epi(device_info, x_nodes, gen_values)

    Q = 1.602e-19
    I_generated = Q * np.trapezoid(gen_values, x_nodes)

    # Sort voltages for stable ramping (0V first, then increasingly negative)
    # We need to track the ramp order
    sorted_indices = np.argsort(-V_range)  # descending: 0, -10, -20, ...
    sorted_V = V_range[sorted_indices]

    cce_sorted = np.zeros(len(V_range))
    I_sorted = np.zeros(len(V_range))
    reached_mask = np.zeros(len(V_range), dtype=bool)
    truncated = False

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

            try:
                ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5)
            except RuntimeError as e:
                # Outside model validity: the cathode is an ohmic contact on the
                # n- epi (no n+ substrate); after reach-through (V_fd) the field
                # at that contact grows with bias and devsim's Newton stops
                # converging. This is a solver/model limit, NOT a physical
                # endpoint (not punch-through) -- stop here and return the
                # converged portion rather than raising (mirrors cv_sweep).
                logger.warning(
                    f"cce_vs_bias: stopped sweep at V={V_target:.1f}V "
                    f"(requested down to {V_range.min():.1f}V): {e}"
                )
                truncated = True
                break

            # Add generation and re-solve
            add_generation_to_dd(device_info, gen_values)

            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=DC_REL_ERROR,
                maximum_iterations=40,
            )

            # Extract CCE
            cce = compute_cce_from_dd(device_info, gen_values, contact="cathode")
            from etna.core.drift_diffusion import extract_steady_state_current

            I_coll = abs(extract_steady_state_current(device_info, "cathode"))

            cce_sorted[idx] = cce
            I_sorted[idx] = I_coll
            reached_mask[idx] = True

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
                relative_error=DC_REL_ERROR,
                maximum_iterations=40,
            )

    finally:
        try:
            devsim.delete_device(device=device)
        except Exception:
            pass

    # Keep only the bias points that actually converged, in their original
    # V_range order — drops the zero-filled placeholders for points past a
    # mid-sweep convergence failure (an unreached point defaulting to
    # cce_values=0 would be indistinguishable from a genuine CCE=0).
    return {
        "voltages": V_range[reached_mask],
        "cce_values": cce_sorted[reached_mask],
        "I_collected": I_sorted[reached_mask],
        "I_generated": I_generated,
        "truncated": truncated,
        "requested_v_stop": float(V_range.min()),
    }


def cce_vs_fluence(
    fluence_range,
    V_bias=-40.0,
    epi_thickness_cm=10e-4,
    N_D_junction=_STACK["N_D_junction"],
    N_D_bulk=_STACK["N_D_bulk"],
    L_transition=_STACK["L_transition"],
    alpha_range_cm=15e-4,
    generation_rate=1e18,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
    device_kwargs=None,
    solver="auto",
):
    """Compute CCE vs proton fluence at fixed reverse bias.

    Each fluence point is a separate damaged device (damaged lifetimes and donor
    profile applied before the Poisson equilibrium). By default all points are
    solved together by the batched numba engine (same mesh, doping, contacts,
    mobilities, SRH and generation as the devsim path; CCE agrees with devsim to
    ~1e-7, tests/test_fluence_engine.py); a point the engine does not converge
    goes through the devsim path (staged device creation, one device per point).

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
    N_D_junction : float
        Junction-side donor concentration (cm^-3). Default: DEFAULT_STACK (2.974271e15).
    N_D_bulk : float
        Bulk-side donor concentration (cm^-3). Default: DEFAULT_STACK (6.335936e13).
    L_transition : float
        Doping transition length (cm). Default: DEFAULT_STACK (0.987 um).
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
    device_kwargs : dict or None
        Additional keyword arguments for create_sic_device (e.g. N_A, T,
        substrate_thickness_cm, doping_profile, N_D); override the defaults.
    solver : {"auto", "engine", "devsim"}
        "auto" (default): batched engine when the device is in its scope (graded
        or uniform epi, 250-400 K, reverse bias), devsim per point otherwise and
        for unconverged points. "engine": as auto but raise if out of scope.
        "devsim": the devsim path only (reference).

    Returns
    -------
    result : dict
        Dictionary with:
        - "fluences": numpy array of fluences (protons/cm^2)
        - "cce_values": numpy array of CCE values (NaN where not valid)
        - "valid": bool array, False where the point was not solved or failed
        - "reason": list of str, "" when valid; otherwise why (full donor
          compensation -- see full_compensation_reason -- or the solver error)
        - "solver": list of str per point: "batched", "devsim" or "" (not solved)
        - "V_bias": float, the bias voltage used
        - "energy_MeV": float
        - "lifetime_model": str

    Points at or beyond full donor removal (eta_removal * Phi_neq >= min N_D of
    the epi; ~4.86e13 p/cm^2 at 62 MeV for N_D_bulk 8.5e13 with the default
    placeholder kappa 0.35 and eta 5) are NOT solved: CCE is NaN with a reason.
    """
    return _cce_vs_fluence_multi(
        fluence_range,
        [damage_params],
        V_bias=V_bias,
        epi_thickness_cm=epi_thickness_cm,
        N_D_junction=N_D_junction,
        N_D_bulk=N_D_bulk,
        L_transition=L_transition,
        alpha_range_cm=alpha_range_cm,
        generation_rate=generation_rate,
        energy_MeV=energy_MeV,
        lifetime_model=lifetime_model,
        device_kwargs=device_kwargs,
        solver=solver,
    )[0]


# Batched-engine scope for the fluence sweeps: the profiles the engine was cross-checked
# on (graded/uniform epi, with or without n+) and the T box of
# etna.api.simulation._ENGINE_BOX. The damaged donor profile itself goes below that box
# (carrier removal), which is why every engine point is pinned against devsim by
# tests/test_fluence_engine.py and unconverged points fall back to devsim.
_FLUENCE_ENGINE_PROFILES = ("graded", "uniform")
_FLUENCE_ENGINE_T = (250.0, 400.0)
# Largest cathode step of the engine's dark Poisson ramp to V_bias (V).
_FLUENCE_ENGINE_STEP_V = 25.0
# Gummel stop: relative change of the cathode current between iterations. 1e-9 keeps
# the engine-vs-devsim CCE difference at the 1e-7 level (devsim's own DC tolerance).
_FLUENCE_ENGINE_GUMMEL_TOL = 1e-9
_FLUENCE_ENGINE_GUMMEL_MAX_ITER = 200


def _fluence_reference(dev_kw):
    """Mesh, doping and material constants of the pristine device, read from devsim.

    One mesh-and-doping-only devsim device (no equations solved). The batched engine
    takes everything from here (node positions, ionised acceptors, pristine donors,
    n_i, mobilities), so both solvers work on the identical discretisation.
    """
    from etna.core.device import create_sic_device

    ref_name = f"fluence_ref_{uuid.uuid4().hex[:8]}"
    info = create_sic_device(device_name=ref_name, **dev_kw)
    try:
        get = lambda name: np.array(
            devsim.get_node_model_values(
                device=ref_name, region=info["region_name"], name=name
            )
        )
        x, donors, acceptors = get("x"), get("Donors"), get("Acceptors")
    finally:
        try:
            devsim.delete_device(device=ref_name)
        except Exception:
            pass
    epi_mask = x >= info["junction_pos"]
    return dict(
        info=info,
        x=x,
        acceptors=acceptors,
        epi_mask=epi_mask,
        pristine_N_D_profile=donors[epi_mask],
    )


def _alpha_generation_values(device_info, x_nodes, alpha_range_cm, generation_rate):
    """Alpha generation G(x) on the mesh nodes, as the CCE sweeps inject it."""
    x_epi = x_nodes - device_info["junction_pos"]
    gen_profile = alpha_generation_profile(x_epi, alpha_range_cm=alpha_range_cm)
    max_profile = np.max(gen_profile)
    if max_profile > 0:
        gen_values = gen_profile * (generation_rate / max_profile)
    else:
        gen_values = np.zeros_like(x_nodes)
    gen_values[x_epi < 0] = 0.0
    return confine_generation_to_epi(device_info, x_nodes, gen_values)


def _fluence_engine_reason(dev_kw, V_bias):
    """Why the batched engine cannot serve this fluence sweep (None = it can)."""
    profile = dev_kw.get("doping_profile", "graded")
    if profile not in _FLUENCE_ENGINE_PROFILES:
        return f"doping_profile={profile!r}"
    T = dev_kw.get("T", 300.0)
    lo, hi = _FLUENCE_ENGINE_T
    if not lo <= T <= hi:
        return f"T={T:g} outside [{lo:g}, {hi:g}]"
    if not (np.isfinite(V_bias) and V_bias <= 0):
        return f"V_bias={V_bias:g} is not a reverse bias"
    return _legacy_depth_reason(dev_kw, abs(V_bias))


def _legacy_depth_reason(dev_kw, v_max):
    """Legacy geometry (no n+, ``nplus_thickness_cm=None``) beyond the engine's
    devsim cross-check depth: devsim hits its punch-through wall there, so the
    engine would return numbers devsim never produced. None = fine."""
    from etna.api.simulation import ENGINE_CROSS_CHECKED_TO_V

    if "nplus_thickness_cm" in dev_kw and dev_kw["nplus_thickness_cm"] is None:
        if v_max > ENGINE_CROSS_CHECKED_TO_V:
            return (
                f"legacy geometry (no n+) beyond {ENGINE_CROSS_CHECKED_TO_V:g} V "
                "(engine not cross-checked there)"
            )
    return None


def _cce_fluence_engine(ref, damaged_list, V_bias, alpha_range_cm, generation_rate):
    """CCE at ``V_bias`` of every damaged device in ``damaged_list``, in one batch.

    Same equations, mesh, doping, contacts, mobilities and generation as the devsim
    path of cce_vs_fluence (Poisson + SG continuity + SRH with the damaged lifetimes
    and donors), solved by the numba engine: dark Poisson ramp to V_bias, then Gummel
    with generation. Returns (cce, converged), both (B,).
    """
    from etna.core.batched.material import Q as Q_ENGINE
    from etna.core.batched.material import thermal_voltage
    from etna.core.batched.poisson import PoissonBatch
    from etna.core.batched.sweeps import _Ramp
    from etna.core.batched.transport import ContinuityBatch

    info, x = ref["info"], ref["x"]
    B, N = len(damaged_list), x.size
    donors = np.zeros((B, N))
    for b, damaged in enumerate(damaged_list):
        donors[b, ref["epi_mask"]] = damaged["N_D_profile"]
    acceptors = np.broadcast_to(ref["acceptors"], (B, N))
    pb = PoissonBatch(x, donors, acceptors, info["n_i"], thermal_voltage(info["T"]))
    cb = ContinuityBatch(
        pb,
        info["mu_n"],
        info["mu_p"],
        np.array([d["tau_n"] for d in damaged_list], dtype=float),
        np.array([d["tau_p"] for d in damaged_list], dtype=float),
    )
    G = _alpha_generation_values(info, x, alpha_range_cm, generation_rate)
    V_cathode = -float(V_bias)
    ramp = _Ramp(pb, 1e-8)
    if V_cathode > 0:
        n_steps = max(1, int(np.ceil(V_cathode / _FLUENCE_ENGINE_STEP_V)))
        for Vc in np.linspace(0.0, V_cathode, n_steps + 1)[1:]:
            ramp.to(Vc)
    tr = cb.solve(
        ramp.psi,
        V_cathode=V_cathode,
        G=G,
        tol=_FLUENCE_ENGINE_GUMMEL_TOL,
        max_iter=_FLUENCE_ENGINE_GUMMEL_MAX_ITER,
        poisson_tol=1e-8,
    )
    J_gen = Q_ENGINE * np.trapezoid(G, x)
    if J_gen <= 0:
        cce = np.zeros(B)
    else:
        cce = np.clip(np.abs(tr.J_cathode) / J_gen, 0.0, 1.0)
    return cce, ramp.converged & tr.converged & np.isfinite(cce)


def _cce_fluence_devsim(dev_kw, damaged, V_bias, alpha_range_cm, generation_rate):
    """CCE of one damaged device at ``V_bias`` with devsim (staged creation)."""
    from etna.core.device import apply_damaged_params, create_sic_device
    from etna.core.drift_diffusion import ramp_bias, setup_sic_drift_diffusion
    from etna.core.poisson import setup_poisson, solve_equilibrium

    dev_name = f"fluence_sweep_{uuid.uuid4().hex[:8]}"
    try:
        # Staged device creation; damage applied BEFORE Poisson setup
        device_info = create_sic_device(device_name=dev_name, **dev_kw)
        apply_damaged_params(device_info, damaged)
        setup_poisson(device_info)
        solve_equilibrium(device_info)
        setup_sic_drift_diffusion(device_info)
        device_info["dd_initialized"] = True

        # Ramp bias (cathode voltage = -V_bias for reverse bias)
        ramp_bias(device_info, -V_bias, contact="cathode", V_step=0.5, V_step_max=2.0)

        x_nodes = np.array(
            devsim.get_node_model_values(
                device=dev_name, region=device_info["region_name"], name="x"
            )
        )
        gen_values = _alpha_generation_values(
            device_info, x_nodes, alpha_range_cm, generation_rate
        )
        add_generation_to_dd(device_info, gen_values)
        devsim.solve(
            type="dc",
            absolute_error=1e10,
            relative_error=DC_REL_ERROR,
            maximum_iterations=40,
        )
        return compute_cce_from_dd(device_info, gen_values, contact="cathode")
    finally:
        try:
            devsim.delete_device(device=dev_name)
        except Exception:
            pass


def _cce_vs_fluence_multi(
    fluence_range,
    damage_params_list,
    V_bias=-40.0,
    epi_thickness_cm=10e-4,
    N_D_junction=_STACK["N_D_junction"],
    N_D_bulk=_STACK["N_D_bulk"],
    L_transition=_STACK["L_transition"],
    alpha_range_cm=15e-4,
    generation_rate=1e18,
    energy_MeV=62.0,
    lifetime_model="linear",
    device_kwargs=None,
    solver="auto",
):
    """cce_vs_fluence for several damage-parameter sets on one device and bias.

    One reference device and (with the engine) one batched solve for all
    len(damage_params_list) * len(fluence_range) points. Returns one
    cce_vs_fluence result dict per damage-parameter set, in order.
    """
    from etna.core.radiation_damage import compute_damaged_params

    if solver not in ("auto", "engine", "devsim"):
        raise ValueError(
            f"solver must be 'auto', 'engine' or 'devsim' (got {solver!r})"
        )
    fluence_range = np.asarray(fluence_range, dtype=float)
    dev_kw = dict(
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=N_D_junction,
        N_D_bulk=N_D_bulk,
        L_transition=L_transition,
    )
    dev_kw.update(device_kwargs or {})
    # --- Pristine parameters for compute_damaged_params, at the device T ---
    T_dev = dev_kw.get("T", 300.0)
    pristine_tau_n = srh_lifetime(T_dev, "electron")
    pristine_tau_p = srh_lifetime(T_dev, "hole")
    ref = _fluence_reference(dev_kw)

    P, F = len(damage_params_list), len(fluence_range)
    cce_values = np.zeros((P, F))
    valid = np.zeros((P, F), dtype=bool)
    reasons = [[""] * F for _ in range(P)]
    solvers = [[""] * F for _ in range(P)]
    todo = []  # (p, i, damaged) to solve
    for p, damage_params in enumerate(damage_params_list):
        for i, fluence in enumerate(fluence_range):
            try:
                damaged = compute_damaged_params(
                    pristine_tau_n=pristine_tau_n,
                    pristine_tau_p=pristine_tau_p,
                    N_D_profile=ref["pristine_N_D_profile"],
                    fluence=fluence,
                    energy_MeV=energy_MeV,
                    damage_params=damage_params,
                    lifetime_model=lifetime_model,
                    T=T_dev,
                )
            except Exception as e:
                logger.warning(f"cce_vs_fluence: failed at fluence={fluence:.2e}: {e}")
                cce_values[p, i] = np.nan
                reasons[p][i] = f"solver failure: {type(e).__name__}: {e}"
                continue
            reason = full_compensation_reason(damaged)
            if reason:
                logger.warning(f"cce_vs_fluence: fluence={fluence:.2e}: {reason}")
                cce_values[p, i] = np.nan
                reasons[p][i] = reason
                continue
            todo.append((p, i, damaged))

    engine_reason = _fluence_engine_reason(dev_kw, V_bias)
    if solver == "engine" and engine_reason is not None:
        raise ValueError(
            f"cce_vs_fluence: batched engine unavailable ({engine_reason})"
        )
    use_engine = solver != "devsim" and engine_reason is None
    if not use_engine and solver == "auto":
        logger.info("cce_vs_fluence: using devsim (%s)", engine_reason)
    remaining = todo
    if use_engine and todo:
        cce_b, conv_b = _cce_fluence_engine(
            ref, [d for _, _, d in todo], V_bias, alpha_range_cm, generation_rate
        )
        remaining = []
        for (p, i, damaged), cce, ok in zip(todo, cce_b, conv_b):
            if ok:
                cce_values[p, i], valid[p, i], solvers[p][i] = cce, True, "batched"
            else:
                logger.warning(
                    f"cce_vs_fluence: batched engine did not converge at "
                    f"fluence={fluence_range[i]:.2e}, using devsim"
                )
                remaining.append((p, i, damaged))

    for p, i, damaged in remaining:
        fluence = fluence_range[i]
        try:
            cce = _cce_fluence_devsim(
                dev_kw, damaged, V_bias, alpha_range_cm, generation_rate
            )
            cce_values[p, i], valid[p, i], solvers[p][i] = cce, True, "devsim"
            logger.info(f"cce_vs_fluence: fluence={fluence:.2e} p/cm^2, CCE={cce:.4f}")
        except Exception as e:
            logger.warning(f"cce_vs_fluence: failed at fluence={fluence:.2e}: {e}")
            cce_values[p, i] = np.nan
            reasons[p][i] = f"solver failure: {type(e).__name__}: {e}"

    return [
        {
            "fluences": fluence_range,
            "cce_values": cce_values[p],
            "valid": valid[p],
            "reason": reasons[p],
            "solver": solvers[p],
            "V_bias": V_bias,
            "energy_MeV": energy_MeV,
            "lifetime_model": lifetime_model,
        }
        for p in range(P)
    ]


def cce_vs_bias_at_fluence(
    V_range,
    fluence,
    epi_thickness_cm=10e-4,
    N_D_junction=_STACK["N_D_junction"],
    N_D_bulk=_STACK["N_D_bulk"],
    L_transition=_STACK["L_transition"],
    alpha_range_cm=15e-4,
    generation_rate=1e18,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
    device_kwargs=None,
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
    N_D_junction : float
        Junction-side donor concentration (cm^-3). Default: DEFAULT_STACK (2.974271e15).
    N_D_bulk : float
        Bulk-side donor concentration (cm^-3). Default: DEFAULT_STACK (6.335936e13).
    L_transition : float
        Doping transition length (cm). Default: DEFAULT_STACK (0.987 um).
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

    device_kwargs : dict or None
        Extra create_sic_device keywords (geometry: substrate_thickness_cm,
        nplus_thickness_cm, N_D_nplus, ...; or doping), overriding the defaults.
        None = the default stack (DEFAULT_STACK); nplus_thickness_cm=None and
        substrate_thickness_cm=1e-4 (or etna.api.device.device_kwargs(legacy_config()))
        give the legacy geometry.

    Returns
    -------
    result : dict
        Dictionary with:
        - "voltages": numpy array of applied voltages (V)
        - "cce_values": numpy array of CCE values (all NaN when not valid)
        - "valid": bool array (all False at full donor compensation)
        - "reason": list of str, "" when valid (see full_compensation_reason)
        - "fluence": float
        - "energy_MeV": float
        - "lifetime_model": str
    """
    dev_kw = dict(
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=N_D_junction,
        N_D_bulk=N_D_bulk,
        L_transition=L_transition,
    )
    dev_kw.update(device_kwargs or {})
    from etna.core.device import apply_damaged_params, create_sic_device
    from etna.core.drift_diffusion import ramp_bias, setup_sic_drift_diffusion
    from etna.core.poisson import setup_poisson, solve_equilibrium
    from etna.core.radiation_damage import compute_damaged_params

    V_range = np.asarray(V_range, dtype=float)

    # --- Extract pristine parameters ---
    pristine_tau_n = srh_lifetime(300.0, "electron")
    pristine_tau_p = srh_lifetime(300.0, "hole")

    # Extract pristine N_D profile from reference device
    ref_id = uuid.uuid4().hex[:8]
    ref_name = f"bias_fluence_ref_{ref_id}"
    ref_info = create_sic_device(
        device_name=ref_name,
        **dev_kw,
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
    reason = full_compensation_reason(damaged)
    if reason:
        logger.warning(f"cce_vs_bias_at_fluence: fluence={fluence:.2e}: {reason}")
        return {
            "voltages": V_range,
            "cce_values": np.full(len(V_range), np.nan),
            "valid": np.zeros(len(V_range), dtype=bool),
            "reason": [reason] * len(V_range),
            "fluence": fluence,
            "energy_MeV": energy_MeV,
            "lifetime_model": lifetime_model,
        }

    # --- Staged device creation ---
    dev_id = uuid.uuid4().hex[:8]
    dev_name = f"bias_at_fluence_{dev_id}"

    device_info = create_sic_device(
        device_name=dev_name,
        **dev_kw,
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
    gen_values = confine_generation_to_epi(device_info, x_nodes, gen_values)

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

            ramp_bias(
                device_info, cathode_V, contact="cathode", V_step=0.5, V_step_max=2.0
            )

            # Add generation and solve
            add_generation_to_dd(device_info, gen_values)
            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=DC_REL_ERROR,
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
                relative_error=DC_REL_ERROR,
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
        "valid": np.ones(len(V_range), dtype=bool),
        "reason": [""] * len(V_range),
        "fluence": fluence,
        "energy_MeV": energy_MeV,
        "lifetime_model": lifetime_model,
    }


def cce_post_anneal(
    fluence,
    T_anneal,
    t_anneal,
    V_bias=-40.0,
    epi_thickness_cm=10e-4,
    alpha_range_cm=15e-4,
    generation_rate=1e18,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
    anneal_params=None,
    device_kwargs=None,
) -> dict:
    """Compute CCE for an irradiated + annealed device at a single operating point.

    Composes annealing recovery (Plan 17-01) with DD-based CCE extraction.
    Uses the same "fresh device per point" staged creation pattern as
    cce_vs_fluence but with compute_annealed_params instead of
    compute_damaged_params.

    Parameters
    ----------
    fluence : float
        Proton fluence (protons/cm^2). Zero returns pristine CCE.
    T_anneal : float
        Annealing temperature (K).
    t_anneal : float
        Annealing time (s).
    V_bias : float
        Reverse bias voltage applied to anode (V, negative). Default -40V.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default 10 um.
    alpha_range_cm : float
        Alpha particle range in SiC (cm). Default 15 um.
    generation_rate : float
        Peak generation rate (cm^-3 s^-1). Default 1e18.
    energy_MeV : float
        Proton energy (MeV). Default 62.0.
    lifetime_model : str
        "linear" or "logarithmic". Default "linear".
    damage_params : RadiationDamageParams or None
        Custom damage parameters.
    anneal_params : AnnealingParams or None
        Custom annealing parameters.

    device_kwargs : dict or None
        Extra create_sic_device keywords (geometry: substrate_thickness_cm,
        nplus_thickness_cm, N_D_nplus, ...; or doping), overriding the defaults.
        None = the default stack (DEFAULT_STACK); nplus_thickness_cm=None and
        substrate_thickness_cm=1e-4 (or etna.api.device.device_kwargs(legacy_config()))
        give the legacy geometry.

    Returns
    -------
    result : dict
        Dictionary with: cce, fluence, T_anneal, t_anneal, V_bias,
        energy_MeV, f_Z12, f_EH67, f_EH4, tau_n, tau_p, valid (bool) and
        reason (str, "" when valid). cce is NaN (valid False) at full donor
        compensation after annealing (see full_compensation_reason) or when
        the solve fails.
    """
    dev_kw = dict(
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=_STACK["N_D_junction"],
        N_D_bulk=_STACK["N_D_bulk"],
        L_transition=_STACK["L_transition"],
    )
    dev_kw.update(device_kwargs or {})
    from etna.core.device import apply_damaged_params, create_sic_device
    from etna.core.drift_diffusion import ramp_bias, setup_sic_drift_diffusion
    from etna.core.poisson import setup_poisson, solve_equilibrium
    from etna.core.radiation_damage import compute_annealed_params

    fluence = float(fluence)

    # Extract pristine parameters
    pristine_tau_n = srh_lifetime(300.0, "electron")
    pristine_tau_p = srh_lifetime(300.0, "hole")

    # Extract pristine N_D profile from a reference device
    ref_id = uuid.uuid4().hex[:8]
    ref_name = f"cce_anneal_ref_{ref_id}"
    ref_info = create_sic_device(
        device_name=ref_name,
        **dev_kw,
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

    # Compute annealed parameters (damage + recovery)
    annealed = compute_annealed_params(
        pristine_tau_n=pristine_tau_n,
        pristine_tau_p=pristine_tau_p,
        N_D_profile=pristine_N_D_profile,
        fluence=fluence,
        energy_MeV=energy_MeV,
        T_anneal=T_anneal,
        t_anneal=t_anneal,
        damage_params=damage_params,
        anneal_params=anneal_params,
        lifetime_model=lifetime_model,
    )

    def _invalid(why):
        return {
            "cce": np.nan,
            "fluence": fluence,
            "T_anneal": T_anneal,
            "t_anneal": t_anneal,
            "V_bias": V_bias,
            "energy_MeV": energy_MeV,
            "f_Z12": annealed.get("f_Z12", 0.0),
            "f_EH67": annealed.get("f_EH67", 0.0),
            "f_EH4": annealed.get("f_EH4", 0.0),
            "tau_n": annealed.get("tau_n", np.nan),
            "tau_p": annealed.get("tau_p", np.nan),
            "valid": False,
            "reason": why,
        }

    reason = full_compensation_reason(annealed)
    if reason:
        logger.warning(f"cce_post_anneal: fluence={fluence:.2e}: {reason}")
        return _invalid(reason)

    # Create fresh device with staged creation
    dev_id = uuid.uuid4().hex[:8]
    dev_name = f"cce_anneal_{dev_id}"

    try:
        device_info = create_sic_device(
            device_name=dev_name,
            **dev_kw,
        )

        # Apply annealed params before Poisson setup
        apply_damaged_params(device_info, annealed)

        # Continue staged setup
        setup_poisson(device_info)
        solve_equilibrium(device_info)
        setup_sic_drift_diffusion(device_info)
        device_info["dd_initialized"] = True

        device = device_info["device_name"]
        region = device_info["region_name"]

        # Ramp bias
        cathode_V = -V_bias
        ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5, V_step_max=2.0)

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
        gen_values = confine_generation_to_epi(device_info, x_nodes, gen_values)

        # Add generation and solve
        add_generation_to_dd(device_info, gen_values)
        devsim.solve(
            type="dc",
            absolute_error=1e10,
            relative_error=DC_REL_ERROR,
            maximum_iterations=40,
        )

        # Extract CCE
        cce = compute_cce_from_dd(device_info, gen_values, contact="cathode")

        logger.info(
            f"cce_post_anneal: fluence={fluence:.2e}, T_anneal={T_anneal}K, "
            f"t_anneal={t_anneal}s, CCE={cce:.4f}"
        )

        return {
            "cce": cce,
            "fluence": fluence,
            "T_anneal": T_anneal,
            "t_anneal": t_anneal,
            "V_bias": V_bias,
            "energy_MeV": energy_MeV,
            "f_Z12": annealed.get("f_Z12", 0.0),
            "f_EH67": annealed.get("f_EH67", 0.0),
            "f_EH4": annealed.get("f_EH4", 0.0),
            "tau_n": annealed["tau_n"],
            "tau_p": annealed["tau_p"],
            "valid": True,
            "reason": "",
        }

    except Exception as e:
        logger.warning(f"cce_post_anneal: failed: {e}")
        return _invalid(f"solver failure: {type(e).__name__}: {e}")
    finally:
        try:
            devsim.delete_device(device=dev_name)
        except Exception:
            pass


def cce_anneal_vs_temperature(
    fluence,
    T_anneal_range,
    t_anneal,
    V_bias=-40.0,
    epi_thickness_cm=10e-4,
    alpha_range_cm=15e-4,
    generation_rate=1e18,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
    anneal_params=None,
) -> dict:
    """Sweep annealing temperature and compute CCE at each point.

    Uses a fresh device per temperature point. Enables plotting CCE
    recovery vs annealing temperature at fixed fluence and anneal time.

    Parameters
    ----------
    fluence : float
        Proton fluence (protons/cm^2).
    T_anneal_range : array_like
        Array of annealing temperatures (K) to sweep.
    t_anneal : float
        Annealing time (s), fixed across sweep.
    V_bias : float
        Reverse bias voltage (V, negative). Default -40V.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default 10 um.
    alpha_range_cm : float
        Alpha particle range (cm). Default 15 um.
    generation_rate : float
        Peak generation rate (cm^-3 s^-1). Default 1e18.
    energy_MeV : float
        Proton energy (MeV). Default 62.0.
    lifetime_model : str
        "linear" or "logarithmic". Default "linear".
    damage_params : RadiationDamageParams or None
        Custom damage parameters.
    anneal_params : AnnealingParams or None
        Custom annealing parameters.

    Returns
    -------
    result : dict
        Dictionary with:
        - "T_anneal_range": numpy array of temperatures (K)
        - "cce_values": numpy array of CCE at each temperature
        - "valid", "reason": per temperature, as in cce_post_anneal
        - "fluence": float
        - "t_anneal": float
        - "V_bias": float
    """
    T_anneal_range = np.asarray(T_anneal_range, dtype=float)
    cce_values = np.zeros(len(T_anneal_range))
    valid = np.zeros(len(T_anneal_range), dtype=bool)
    reasons = [""] * len(T_anneal_range)

    for i, T_ann in enumerate(T_anneal_range):
        result = cce_post_anneal(
            fluence=fluence,
            T_anneal=T_ann,
            t_anneal=t_anneal,
            V_bias=V_bias,
            epi_thickness_cm=epi_thickness_cm,
            alpha_range_cm=alpha_range_cm,
            generation_rate=generation_rate,
            energy_MeV=energy_MeV,
            lifetime_model=lifetime_model,
            damage_params=damage_params,
            anneal_params=anneal_params,
        )
        cce_values[i] = result["cce"]
        valid[i] = result["valid"]
        reasons[i] = result["reason"]

        logger.info(
            f"cce_anneal_vs_temperature: T={T_ann:.0f}K, CCE={result['cce']:.4f}"
        )

    return {
        "T_anneal_range": T_anneal_range,
        "cce_values": cce_values,
        "valid": valid,
        "reason": reasons,
        "fluence": fluence,
        "t_anneal": t_anneal,
        "V_bias": V_bias,
    }


def cce_vs_epi_thickness(
    epi_range_cm,
    V_bias=-3.0,
    alpha_range_cm=15e-4,
    generation_rate=1e18,
    device_kwargs=None,
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

    device_kwargs : dict or None
        Extra create_dd_device keywords except epi_thickness_cm (geometry or doping),
        overriding the defaults (DEFAULT_STACK); nplus_thickness_cm=None and
        substrate_thickness_cm=1e-4 give the legacy geometry.

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
    from etna.core.drift_diffusion import (
        create_dd_device,
        ramp_bias,
    )

    epi_range_cm = np.asarray(epi_range_cm, dtype=float)
    cce_values = np.zeros(len(epi_range_cm))

    Q = 1.602e-19

    for i, epi_thick in enumerate(epi_range_cm):
        dev_id = uuid.uuid4().hex[:8]
        dev_name = f"cce_epi_{dev_id}"

        dev_kw = dict(
            doping_profile="graded",
            N_D_junction=_STACK["N_D_junction"],
            N_D_bulk=_STACK["N_D_bulk"],
            L_transition=_STACK["L_transition"],
        )
        dev_kw.update(device_kwargs or {})  # geometry / doping overrides (not epi)
        device_info = create_dd_device(
            device_name=dev_name, epi_thickness_cm=epi_thick, **dev_kw
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
            gen_values = confine_generation_to_epi(device_info, x_nodes, gen_values)

            I_generated = Q * np.trapezoid(gen_values, x_nodes)

            # Ramp to bias (cathode voltage is -V_bias)
            cathode_V = -V_bias
            ramp_bias(
                device_info, cathode_V, contact="cathode", V_step=0.5, V_step_max=2.0
            )

            # Add generation and solve
            add_generation_to_dd(device_info, gen_values)
            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=DC_REL_ERROR,
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
    # compare on the biases the DD sweep reached (a truncated sweep is shorter)
    V_range = np.asarray(dd_result["voltages"], dtype=float)

    # Hecht equation CCE (fully depleted approximation, d = epi thickness)
    cce_hecht = hecht_cce(V_range, d=epi_thickness_cm)

    # Partial depletion Hecht using a simple W(V) model
    # W ~ sqrt(2 * eps * V / (q * N_D)) with calibrated N_D
    # For graded doping, use an effective N_D ~ geometric mean
    # (the device doping of the DD sweep: device_kwargs, else DEFAULT_STACK; was a
    # hardcoded 2.90e15 / 8.50e13)
    _dk = kwargs.get("device_kwargs") or {}
    N_D_eff = np.sqrt(
        (_dk.get("N_D_junction") or _STACK["N_D_junction"])
        * (_dk.get("N_D_bulk") or _STACK["N_D_bulk"])
    )  # ~ 4.34e14 for the default stack
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
    from etna.core.validation import compute_agreement_metrics

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
        "truncated": bool(dd_result.get("truncated", False)),
        "cce_dd": cce_dd,
        "cce_hecht": np.asarray(cce_hecht, dtype=float),
        "cce_hecht_partial": np.asarray(cce_hecht_partial, dtype=float),
        "max_deviation": max_dev,
        "agreement_metrics_hecht": metrics_vs_hecht,
        "agreement_metrics_partial": metrics_vs_partial,
        "regime_notes": regime_notes,
    }
