"""Hurkx trap-assisted tunneling and surface recombination dark current models.

Implements effective dark current generation for 4H-SiC detectors, combining:
1. Modified SRH recombination with Z1/2 deep-level n1/p1 (not midgap)
2. Effective depletion-region generation calibrated via N_t parameter
3. Hurkx field-enhancement factor Gamma for voltage dependence
4. Surface recombination velocity at contacts

Physical context:
    Standard midgap SRH in 4H-SiC produces ~1e-49 A dark current because
    n_i ~ 5e-9 cm^-3. Real dark current (~18 pA at -30V) arises from multiple
    mechanisms including perimeter leakage (2D), surface states, and
    trap-assisted tunneling that cannot be fully captured in 1D.

    This module uses an effective generation rate G_eff = N_t (cm^-3 s^-1)
    applied in the depletion region, where N_t is calibrated to match
    experimental dark current. The E-field-dependent Gamma factor provides
    voltage-dependent scaling. Default N_t = 1e12 cm^-3 s^-1 targets ~18 pA.

All units CGS (cm, cm^-3, A/cm^2, V) per devsim convention.

References:
    - Hurkx et al., IEEE TED 39(2), 1992
    - Schenk, SSE 35(11), 1992 (Gamma approximation)
    - Kimoto & Cooper, "Fundamentals of Silicon Carbide Technology"
"""

import logging

import devsim
import numpy as np
from devsim.python_packages.model_create import (
    CreateNodeModel,
    CreateNodeModelDerivative,
)

from src.drift_diffusion import (
    create_dd_device,
    extract_contact_current,
    ramp_bias,
)

logger = logging.getLogger(__name__)

# Physical constants in SI
_HBAR_SI = 1.0546e-34  # J*s
_M0_SI = 9.109e-31  # kg
_Q_SI = 1.602e-19  # C (also J/eV)
_K_B_EV = 8.617e-5  # eV/K


def setup_tat_model(device_info, E_t=None, m_t=None, N_t=None):
    """Add TAT-based dark current generation to drift-diffusion equations.

    Replaces the existing SRH-only ElectronGeneration/HoleGeneration with
    a combined model:
    1. SRH recombination with Z1/2 trap-level n1/p1 (physical)
    2. Effective generation rate G_eff = N_t * Gamma in depletion region

    The N_t parameter (cm^-3 s^-1) is the primary calibration knob. It
    represents the effective volumetric generation rate from all dark
    current mechanisms (TAT, perimeter leakage, surface states).

    Parameters
    ----------
    device_info : dict
        Device info dict from create_dd_device(). Must have DD equations
        already set up (dd_initialized=True).
    E_t : float, optional
        Trap energy below Ec (eV). Default from device_info["params"].E_t.
    m_t : float, optional
        Tunneling effective mass in units of m0. Default from params.m_t.
    N_t : float, optional
        Effective generation rate (cm^-3 s^-1). Default from params.N_t.
        Primary calibration parameter. Default 1e12 targets ~8 pA at -30V.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    params = device_info["params"]

    if E_t is None:
        E_t = params.E_t
    if m_t is None:
        m_t = params.m_t
    if N_t is None:
        N_t = params.N_t

    T = device_info["T"]
    kT_eV = _K_B_EV * T
    n_i = device_info["n_i"]
    E_g = device_info["E_g"]

    # Store TAT parameters
    devsim.set_parameter(device=device, region=region, name="E_t", value=E_t)
    devsim.set_parameter(device=device, region=region, name="m_t", value=m_t)
    devsim.set_parameter(device=device, region=region, name="N_t", value=N_t)

    # --- Step 1: Node-averaged electric field ---
    _compute_node_efield(device, region)

    # --- Step 2: Kt and Gamma from Hurkx model ---
    m_tunnel_kg = m_t * _M0_SI
    E_t_J = E_t * _Q_SI

    Kt_numerator = (4.0 / 3.0) * np.sqrt(2.0 * m_tunnel_kg) * E_t_J**1.5
    Kt_denom_const = _Q_SI * _HBAR_SI * 100.0

    Kt_expr = f"{Kt_numerator} / ({Kt_denom_const} * E_field_node)"
    CreateNodeModel(device, region, "Kt_TAT", Kt_expr)

    _compute_gamma_factors(device, region)

    # --- Step 3: Trap-level n1, p1 for Z1/2 ---
    E_t_from_midgap = E_g / 2.0 - E_t
    n1_tat = n_i * np.exp(E_t_from_midgap / kT_eV)
    p1_tat = n_i * np.exp(-E_t_from_midgap / kT_eV)

    devsim.set_parameter(device=device, region=region, name="n1_tat", value=n1_tat)
    devsim.set_parameter(device=device, region=region, name="p1_tat", value=p1_tat)

    ni2 = n_i**2

    # --- Step 4: Effective dark current generation ---
    # The total recombination/generation rate has two parts:
    #
    # Part A: SRH with Z1/2 trap level (replaces midgap SRH)
    #   U_SRH = (np - ni^2) / (taup*(n+n1) + taun*(p+p1))
    #   where n1, p1 are for Z1/2 trap (not midgap n_i)
    #
    # Part B: Effective generation (always negative = generation)
    #   G_eff = -N_t * Gamma * f_depl
    #   where f_depl identifies the depletion region using E-field:
    #   f_depl = E_field_node / E_ref, clamped to [0, 1]
    #   This gives ~1 in depletion, ~0 in quasi-neutral regions.
    #
    # The contact current from devsim captures the net effect of both.

    # Reference field for depletion identification
    # Built-in field at junction is ~70 kV/cm; use half as reference
    E_ref = 3.5e4  # V/cm

    # Store G0_TAT = N_t as generation rate parameter
    devsim.set_parameter(device=device, region=region, name="G0_TAT", value=N_t)
    devsim.set_parameter(device=device, region=region, name="E_ref_TAT", value=E_ref)

    # Standard SRH with Z1/2 trap level (for component decomposition)
    U_SRH_tat = (
        f"(Electrons * Holes - {ni2}) / "
        f"(taup * (Electrons + {n1_tat}) + taun * (Holes + {p1_tat}))"
    )
    CreateNodeModel(device, region, "U_SRH_only", U_SRH_tat)

    # Full generation model: SRH_trap + effective generation
    # f_depl = min(E_field_node / E_ref, 1.0) — larger in high-field depletion region
    # G_eff = -N_t * Gamma_n * f_depl (negative = net generation)
    # U_TAT = U_SRH_tat + G_eff
    U_TAT = (
        f"(Electrons * Holes - {ni2}) / "
        f"(taup * (Electrons + n1_tat * Gamma_n) + "
        f"taun * (Holes + p1_tat * Gamma_p)) - "
        f"G0_TAT * Gamma_n * min(E_field_node / E_ref_TAT, 1.0)"
    )
    CreateNodeModel(device, region, "U_TAT", U_TAT)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "U_TAT", U_TAT, var)

    # --- Step 5: Replace generation models ---
    Gn_TAT = "-ElectronCharge * U_TAT"
    Gp_TAT = "+ElectronCharge * U_TAT"
    CreateNodeModel(device, region, "ElectronGeneration", Gn_TAT)
    CreateNodeModel(device, region, "HoleGeneration", Gp_TAT)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn_TAT, var)
        CreateNodeModelDerivative(device, region, "HoleGeneration", Gp_TAT, var)

    device_info["tat_initialized"] = True
    logger.info(
        f"TAT model setup: E_t={E_t} eV, m_t={m_t} m0, N_t={N_t:.1e} cm^-3/s, "
        f"n1_tat={n1_tat:.2e}, p1_tat={p1_tat:.2e}"
    )


def _compute_node_efield(device, region):
    """Average absolute edge electric field to nodes.

    Computes the arithmetic mean of adjacent edge |E| values for each node.
    Boundary nodes use the single adjacent edge. Result clamped to minimum
    1e3 V/cm to avoid division by zero in Kt computation.
    """
    E_edge = np.array(
        devsim.get_edge_model_values(device=device, region=region, name="ElectricField")
    )
    node_x = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    n_nodes = len(node_x)
    n_edges = len(E_edge)

    E_node = np.zeros(n_nodes)
    abs_E = np.abs(E_edge)

    for i in range(n_nodes):
        if i == 0:
            E_node[i] = abs_E[0]
        elif i == n_nodes - 1:
            E_node[i] = abs_E[n_edges - 1]
        else:
            E_node[i] = 0.5 * (abs_E[i - 1] + abs_E[i])

    E_node = np.maximum(E_node, 1e3)

    devsim.node_model(device=device, region=region, name="E_field_node", equation="1e3")
    devsim.set_node_values(
        device=device, region=region, name="E_field_node", values=E_node.tolist()
    )


def _compute_gamma_factors(device, region):
    """Compute Gamma_n and Gamma_p field-enhancement factors numerically.

    Uses the Schenk approximation:
        Kt < 4: Gamma = sqrt(pi)/(2*Kt) * exp(1/(4*Kt^2))
        Kt >= 4: Gamma = 1 (low field, no enhancement)

    For simplicity, Gamma_n = Gamma_p (same trap, symmetric tunneling).
    """
    Kt = np.array(
        devsim.get_node_model_values(device=device, region=region, name="Kt_TAT")
    )

    Gamma = np.ones_like(Kt)
    high_field = Kt < 4.0

    if np.any(high_field):
        Kt_hf = Kt[high_field]
        Kt_hf = np.maximum(Kt_hf, 0.05)
        Gamma[high_field] = (np.sqrt(np.pi) / (2.0 * Kt_hf)) * np.exp(
            1.0 / (4.0 * Kt_hf**2)
        )

    devsim.node_model(device=device, region=region, name="Gamma_n", equation="1.0")
    devsim.set_node_values(
        device=device, region=region, name="Gamma_n", values=Gamma.tolist()
    )
    devsim.node_model(device=device, region=region, name="Gamma_p", equation="1.0")
    devsim.set_node_values(
        device=device, region=region, name="Gamma_p", values=Gamma.tolist()
    )


def setup_surface_recombination(device_info, S_n=None, S_p=None, contact="cathode"):
    """Add surface recombination current at a contact.

    Computes SRV current density at the contact and stores it for extraction.
    The ohmic contact BC in devsim already enforces equilibrium carriers;
    the SRV model provides an additional analytical current component.

    Parameters
    ----------
    device_info : dict
        Device info dict from create_dd_device().
    S_n : float, optional
        Electron surface recombination velocity (cm/s). Default from params.
    S_p : float, optional
        Hole surface recombination velocity (cm/s). Default from params.
    contact : str
        Contact name at which to apply SRV. Default "cathode".
    """
    device = device_info["device_name"]
    params = device_info["params"]
    n_i = device_info["n_i"]

    if S_n is None:
        S_n = params.S_n
    if S_p is None:
        S_p = params.S_p

    ni2 = n_i**2

    devsim.set_parameter(device=device, name=f"S_n_{contact}", value=S_n)
    devsim.set_parameter(device=device, name=f"S_p_{contact}", value=S_p)

    srv_model = (
        f"ElectronCharge * (Electrons * Holes - {ni2}) / "
        f"((Electrons + {n_i}) / {S_p} + (Holes + {n_i}) / {S_n})"
    )

    model_name = f"srv_{contact}"
    devsim.contact_node_model(
        device=device, contact=contact, name=model_name, equation=srv_model
    )

    for var in ("Electrons", "Holes"):
        deriv_name = f"{model_name}:{var}"
        devsim.contact_node_model(
            device=device,
            contact=contact,
            name=deriv_name,
            equation="0",
        )

    device_info["srv_initialized"] = True
    device_info["srv_contact"] = contact
    device_info["srv_S_n"] = S_n
    device_info["srv_S_p"] = S_p
    logger.info(
        f"Surface recombination setup at {contact}: S_n={S_n:.1e}, S_p={S_p:.1e} cm/s"
    )


def extract_dark_current_components(device_info, area=0.05):
    """Extract dark current with component decomposition.

    Returns total current and its SRH, TAT, and SRV contributions.

    Parameters
    ----------
    device_info : dict
        Device info dict with TAT model set up.
    area : float
        Device area (cm^2). Default 0.05 cm^2 (5 mm^2).

    Returns
    -------
    dict
        Dictionary with keys: J_total, J_SRH, J_TAT, J_SRV,
        I_total, I_SRH, I_TAT, I_SRV.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]

    J_total = extract_contact_current(device_info, contact="cathode")

    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))

    q = _Q_SI
    if device_info.get("tat_initialized", False):
        U_SRH_only = np.array(
            devsim.get_node_model_values(
                device=device, region=region, name="U_SRH_only"
            )
        )
        U_TAT = np.array(
            devsim.get_node_model_values(device=device, region=region, name="U_TAT")
        )
        J_SRH = q * np.trapz(U_SRH_only, x)
        J_TAT_bulk = q * np.trapz(U_TAT - U_SRH_only, x)
    else:
        J_SRH = J_total
        J_TAT_bulk = 0.0

    J_SRV = 0.0
    if device_info.get("srv_initialized", False):
        contact = device_info.get("srv_contact", "cathode")
        S_n = device_info.get("srv_S_n", 0.0)
        S_p = device_info.get("srv_S_p", 0.0)
        n_i = device_info["n_i"]

        n_vals = np.array(
            devsim.get_node_model_values(device=device, region=region, name="Electrons")
        )
        p_vals = np.array(
            devsim.get_node_model_values(device=device, region=region, name="Holes")
        )

        if contact == "cathode":
            n_c, p_c = n_vals[-1], p_vals[-1]
        else:
            n_c, p_c = n_vals[0], p_vals[0]

        ni2 = n_i**2
        if S_n > 0 and S_p > 0:
            J_SRV = q * (n_c * p_c - ni2) / ((n_c + n_i) / S_p + (p_c + n_i) / S_n)

    return {
        "J_total": J_total,
        "J_SRH": J_SRH,
        "J_TAT": J_TAT_bulk,
        "J_SRV": J_SRV,
        "I_total": J_total * area,
        "I_SRH": J_SRH * area,
        "I_TAT": J_TAT_bulk * area,
        "I_SRV": J_SRV * area,
    }


def dark_current_sweep(device_info, V_range, area=0.05, V_step=0.5):
    """Sweep reverse voltage and extract dark current components.

    Parameters
    ----------
    device_info : dict
        Device info dict with TAT model set up.
    V_range : array_like
        Array of reverse voltages (V). Negative = reverse bias on anode.
    area : float
        Device area (cm^2). Default 0.05 cm^2.
    V_step : float
        Maximum voltage step for ramping (V).

    Returns
    -------
    dict
        Dictionary with arrays: voltages, I_total, I_SRH, I_TAT, I_SRV,
        J_total, J_SRH, J_TAT, J_SRV.
    """
    import devsim.python_packages.simple_physics as simple_physics

    device = device_info["device_name"]
    region = device_info["region_name"]
    V_range = np.asarray(V_range, dtype=float)

    results = {
        "voltages": [],
        "I_total": [],
        "I_SRH": [],
        "I_TAT": [],
        "I_SRV": [],
        "J_total": [],
        "J_SRH": [],
        "J_TAT": [],
        "J_SRV": [],
    }

    bias_name = simple_physics.GetContactBiasName("anode")

    try:
        current_V = devsim.get_parameter(device=device, name=bias_name)
    except devsim.error:
        current_V = 0.0

    for V_target in V_range:
        delta = V_target - current_V
        if abs(delta) > 1e-12:
            n_steps = max(1, int(np.ceil(abs(delta) / V_step)))
            V_intermediates = np.linspace(current_V, V_target, n_steps + 1)[1:]

            for V_int in V_intermediates:
                devsim.set_parameter(device=device, name=bias_name, value=V_int)
                _compute_node_efield(device, region)
                _compute_gamma_factors(device, region)
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
                        logger.warning(
                            f"dark_current_sweep: failed at V={V_int:.3f}V: {e}"
                        )
                        for key in results:
                            results[key] = np.array(results[key])
                        return results

            current_V = V_target

        _compute_node_efield(device, region)
        _compute_gamma_factors(device, region)

        components = extract_dark_current_components(device_info, area=area)
        results["voltages"].append(V_target)
        for key in (
            "I_total",
            "I_SRH",
            "I_TAT",
            "I_SRV",
            "J_total",
            "J_SRH",
            "J_TAT",
            "J_SRV",
        ):
            results[key].append(components[key])

    for key in results:
        results[key] = np.array(results[key])

    return results


def create_dark_current_device(T=300, N_t=None, S_n=None, S_p=None, **kwargs):
    """Convenience function to create a DD device with TAT and SRV models.

    Parameters
    ----------
    T : float
        Temperature (K).
    N_t : float, optional
        Effective generation rate (cm^-3 s^-1). Default from SiC4H_Parameters.
    S_n : float, optional
        Electron surface recombination velocity (cm/s).
    S_p : float, optional
        Hole surface recombination velocity (cm/s).
    **kwargs
        Additional arguments passed to create_dd_device().

    Returns
    -------
    device_info : dict
        Device info dict with TAT and SRV models initialized.
    """
    device_info = create_dd_device(T=T, **kwargs)

    setup_tat_model(device_info, N_t=N_t)
    setup_surface_recombination(device_info, S_n=S_n, S_p=S_p)

    return device_info


def sensitivity_sweep(
    param_name, param_values, V_eval=-30.0, T=300, area=0.05, base_kwargs=None
):
    """Sweep a single parameter and evaluate dark current at a fixed voltage.

    Creates a fresh device for each parameter value, ramps to V_eval, extracts
    dark current components, and cleans up the device.

    Parameters
    ----------
    param_name : str
        Parameter to sweep. One of:
        - "epi_thickness_cm", "N_D" : passed to create_dd_device via kwargs
        - "N_t", "S_n", "S_p" : passed to create_dark_current_device
    param_values : array_like
        Values of the swept parameter.
    V_eval : float
        Reverse voltage at which to evaluate dark current (V). Default -30V.
    T : float
        Temperature (K). Default 300.
    area : float
        Device area (cm^2). Default 0.05.
    base_kwargs : dict or None
        Baseline device parameters passed to create_dark_current_device.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: param_name column, I_total, I_SRH, I_TAT, I_SRV.
    """
    import pandas as pd

    import devsim.python_packages.simple_physics as simple_physics

    if base_kwargs is None:
        base_kwargs = {}

    param_values = np.asarray(param_values)

    # Classify the parameter target
    dd_device_params = {"epi_thickness_cm", "N_D"}
    dark_current_params = {"N_t", "S_n", "S_p"}

    if param_name not in dd_device_params | dark_current_params:
        raise ValueError(
            f"Unknown param_name '{param_name}'. "
            f"Must be one of {sorted(dd_device_params | dark_current_params)}"
        )

    records = []
    for i, val in enumerate(param_values):
        device_name = f"sweep_{param_name}_{i}"
        kwargs = dict(base_kwargs)
        kwargs["T"] = T

        # Route parameter to the correct function
        dc_kwargs = {}
        if param_name in dd_device_params:
            kwargs[param_name] = val
        else:
            dc_kwargs[param_name] = val

        kwargs["device_name"] = device_name

        try:
            device_info = create_dark_current_device(**dc_kwargs, **kwargs)
            device = device_info["device_name"]
            region = device_info["region_name"]

            # Ramp to V_eval in steps
            bias_name = simple_physics.GetContactBiasName("anode")
            V_step = 0.5
            current_V = 0.0
            n_steps = max(1, int(np.ceil(abs(V_eval) / V_step)))
            V_intermediates = np.linspace(0, V_eval, n_steps + 1)[1:]

            ramp_ok = True
            for V_int in V_intermediates:
                devsim.set_parameter(device=device, name=bias_name, value=V_int)
                _compute_node_efield(device, region)
                _compute_gamma_factors(device, region)
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
                        logger.warning(
                            f"sensitivity_sweep: ramp failed at V={V_int:.3f} "
                            f"for {param_name}={val}: {e}"
                        )
                        ramp_ok = False
                        break

            if ramp_ok:
                _compute_node_efield(device, region)
                _compute_gamma_factors(device, region)
                components = extract_dark_current_components(device_info, area=area)
                records.append(
                    {
                        param_name: val,
                        "I_total": components["I_total"],
                        "I_SRH": components["I_SRH"],
                        "I_TAT": components["I_TAT"],
                        "I_SRV": components["I_SRV"],
                    }
                )
            else:
                records.append(
                    {
                        param_name: val,
                        "I_total": np.nan,
                        "I_SRH": np.nan,
                        "I_TAT": np.nan,
                        "I_SRV": np.nan,
                    }
                )
        finally:
            try:
                devsim.delete_device(device=device_name)
            except Exception:
                pass

    return pd.DataFrame(records)


def plot_dark_current_decomposition(sweep_result, ax=None, title=None):
    """Plot dark current decomposition from a voltage sweep.

    Plots |I_total|, |I_SRH|, |I_TAT|, |I_SRV| vs voltage on semilogy axes
    with publication-quality styling matching plotting.py conventions.

    Parameters
    ----------
    sweep_result : dict
        Output of dark_current_sweep() with keys: voltages, I_total,
        I_SRH, I_TAT, I_SRV.
    ax : matplotlib.axes.Axes or None
        Axes to plot on. If None, creates new figure.
    title : str or None
        Plot title. Default: "Dark Current Decomposition".

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    V = np.asarray(sweep_result["voltages"])

    components = [
        ("I_total", "Total", "k", 2.0, "-"),
        ("I_SRH", "SRH (bulk)", "C0", 1.5, "--"),
        ("I_TAT", "TAT (effective)", "C1", 1.5, "-."),
        ("I_SRV", "SRV (surface)", "C2", 1.5, ":"),
    ]

    for key, label, color, lw, ls in components:
        I = np.abs(np.asarray(sweep_result[key]))
        # Only plot if nonzero somewhere
        if np.any(I > 0):
            ax.semilogy(V, I, color=color, linewidth=lw, linestyle=ls, label=label)

    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("|Dark Current| (A)")
    ax.set_title(title or "Dark Current Decomposition")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    return ax
