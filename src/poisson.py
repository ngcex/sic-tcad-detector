"""devsim Poisson equation solver for 4H-SiC p+/n- diode.

Provides functions to set up and solve the Poisson equation on a 4H-SiC
device created by device.py. Includes voltage ramping for reverse bias
sweeps and extraction of electric field and depletion width.

The solver uses a Poisson-only formulation with clamped exponentials
to handle SiC's extremely low n_i (~5e-9 cm^-3), which causes overflow
in standard devsim formulations under reverse bias.

For depletion width computation under bias, uses the analytical one-sided
depletion formula (validated against numerical equilibrium at 0V) with
punch-through clamping at the epitaxial layer thickness. This approach
avoids the numerical artifacts inherent in applying Boltzmann carrier
statistics under non-equilibrium bias in a Poisson-only framework.

All units CGS (cm, V, V/cm) per devsim convention.

References:
    - devsim simple_physics module pattern
    - devsim diode example: https://devsim.net/examples_diode.html
"""

import logging

import devsim
import devsim.python_packages.simple_physics as simple_physics
from devsim.python_packages.model_create import (
    CreateSolution,
    CreateNodeModel,
    CreateNodeModelDerivative,
    CreateEdgeModel,
    CreateEdgeModelDerivatives,
    InNodeModelList,
)
import numpy as np

from src.analytical import (
    built_in_potential,
    depletion_width as analytical_depletion_width,
)

logger = logging.getLogger(__name__)

# Exponential clamp to prevent overflow with SiC's extreme potential/V_t ratios
_EXP_CLAMP = 700


def _create_sic_potential_only(device, region):
    """Create Poisson equation models for 4H-SiC.

    Custom version of simple_physics.CreateSiliconPotentialOnly with
    clamped exponentials for wide-bandgap stability. Uses direct
    n_i*exp(+/-Potential/V_t) formulations with the argument clamped
    to +/-700 to prevent overflow (exp(700) ~ 1e304, safely within
    double precision range).
    """
    if not InNodeModelList(device, region, "Potential"):
        logger.info("Creating Node Solution Potential")
        CreateSolution(device, region, "Potential")

    elec_arg = f"min(max(Potential/V_t, -{_EXP_CLAMP}), {_EXP_CLAMP})"
    hole_arg = f"min(max(-Potential/V_t, -{_EXP_CLAMP}), {_EXP_CLAMP})"

    elec_i = f"n_i*exp({elec_arg})"
    hole_i = f"n_i*exp({hole_arg})"
    charge_i = "kahan3(IntrinsicHoles, -IntrinsicElectrons, NetDoping)"
    pcharge_i = "-ElectronCharge * IntrinsicCharge"

    for name, expr in (
        ("IntrinsicElectrons", elec_i),
        ("IntrinsicHoles", hole_i),
        ("IntrinsicCharge", charge_i),
        ("PotentialIntrinsicCharge", pcharge_i),
    ):
        CreateNodeModel(device, region, name, expr)
        CreateNodeModelDerivative(device, region, name, expr, "Potential")

    for name, expr in (
        ("ElectricField", "(Potential@n0-Potential@n1)*EdgeInverseLength"),
        ("PotentialEdgeFlux", "Permittivity * ElectricField"),
    ):
        CreateEdgeModel(device, region, name, expr)
        CreateEdgeModelDerivatives(device, region, name, expr, "Potential")

    devsim.equation(
        device=device,
        region=region,
        name="PotentialEquation",
        variable_name="Potential",
        node_model="PotentialIntrinsicCharge",
        edge_model="PotentialEdgeFlux",
        variable_update="log_damp",
    )


def setup_poisson(device_info):
    """Set up the Poisson equation on the device.

    Uses a custom SiC-adapted Poisson model with clamped exponentials
    to prevent overflow at reverse biases up to -60V.

    Parameters
    ----------
    device_info : dict
        Device info dict returned by create_sic_device().
    """
    device = device_info["device_name"]
    region = device_info["region_name"]

    _create_sic_potential_only(device, region)

    for contact in ("anode", "cathode"):
        bias_name = simple_physics.GetContactBiasName(contact)
        devsim.set_parameter(device=device, name=bias_name, value=0.0)
        simple_physics.CreateSiliconPotentialOnlyContact(device, region, contact)

    logger.info("Poisson equation setup complete")


def solve_equilibrium(device_info):
    """Solve the Poisson equation at 0V (thermal equilibrium).

    Parameters
    ----------
    device_info : dict
        Device info dict returned by create_sic_device().

    Raises
    ------
    RuntimeError
        If the solver fails to converge.
    """
    try:
        devsim.solve(
            type="dc",
            absolute_error=1e10,
            relative_error=1e-10,
            maximum_iterations=40,
        )
        logger.info("Equilibrium solve converged")
    except devsim.error:
        try:
            devsim.solve(
                type="dc",
                absolute_error=1e12,
                relative_error=1e-8,
                maximum_iterations=100,
            )
            logger.info("Equilibrium solve converged (relaxed)")
        except devsim.error as e:
            raise RuntimeError(f"Equilibrium solve failed to converge: {e}") from e


def ramp_voltage(
    device_info,
    contact_name="cathode",
    V_start=0.0,
    V_end=-60.0,
    V_step=-0.5,
    abs_err=1e10,
    rel_err=1e-10,
    max_iter=40,
):
    """Ramp bias voltage on a contact in small steps.

    For reverse bias on the cathode, use negative V_end and V_step.

    Parameters
    ----------
    device_info : dict
        Device info dict.
    contact_name : str
        Contact to ramp ("anode" or "cathode").
    V_start : float
        Starting voltage (V).
    V_end : float
        Ending voltage (V).
    V_step : float
        Voltage step (V). Negative for reverse bias ramp.
    abs_err : float
        Absolute error tolerance for solver.
    rel_err : float
        Relative error tolerance for solver.
    max_iter : int
        Maximum Newton iterations.

    Returns
    -------
    results : list of tuple
        List of (voltage, x_array, E_field_array) at each solved bias point.
    """
    device = device_info["device_name"]
    bias_name = simple_physics.GetContactBiasName(contact_name)

    results = []
    V = V_start + V_step

    if V_step < 0:
        condition = lambda v: v >= V_end
    else:
        condition = lambda v: v <= V_end

    while condition(V):
        devsim.set_parameter(device=device, name=bias_name, value=V)
        try:
            devsim.solve(
                type="dc",
                absolute_error=abs_err,
                relative_error=rel_err,
                maximum_iterations=max_iter,
            )
        except devsim.error:
            logger.warning(f"Solve failed at V={V:.2f}V, trying relaxed settings")
            try:
                devsim.solve(
                    type="dc",
                    absolute_error=abs_err * 100,
                    relative_error=rel_err * 10,
                    maximum_iterations=max_iter * 3,
                )
            except devsim.error as e:
                logger.error(
                    f"Solve failed at V={V:.2f}V even with relaxed settings: {e}"
                )
                break

        x, E = extract_electric_field(device_info)
        results.append((V, x, E))
        logger.debug(f"Solved at V={V:.2f}V, E_max={np.max(np.abs(E)):.2e} V/cm")

        V += V_step
        V = round(V, 10)

    return results


def extract_electric_field(device_info):
    """Extract electric field from the solved device.

    Returns positions at edge centers and E-field values.

    Parameters
    ----------
    device_info : dict
        Device info dict.

    Returns
    -------
    x_centers : ndarray
        Edge center positions (cm).
    E_field : ndarray
        Electric field at edge centers (V/cm).
    """
    device = device_info["device_name"]
    region = device_info["region_name"]

    x_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    E_edges = np.array(
        devsim.get_edge_model_values(device=device, region=region, name="ElectricField")
    )
    x_centers = 0.5 * (x_nodes[:-1] + x_nodes[1:])
    return x_centers, E_edges


def extract_depletion_width(device_info, V_applied=0.0):
    """Extract depletion width using analytical formula.

    Uses the one-sided depletion approximation with built-in potential
    computed from the device's ionized doping concentrations and n_i.
    The analytical result is validated against the numerical equilibrium
    solution at 0V and provides accurate punch-through behavior.

    Parameters
    ----------
    device_info : dict
        Device info dict.
    V_applied : float
        Applied bias voltage (V). Negative for reverse bias.

    Returns
    -------
    W : float
        Depletion width (cm), measured from junction into the n-side.

    Note
    ----
    This returns W computed from the **analytical** one-sided depletion
    approximation (W = sqrt(2*eps*(V_bi - V)/(q*N_D))), NOT a W extracted
    from the numerical Poisson electric field solution. The analytical
    formula was validated against the numerical equilibrium solution at
    0V (agreement within ~20%), but under reverse bias it does not account
    for carrier redistribution or non-abrupt junction profiles.

    Phase 2 drift-diffusion code should NOT assume this function returns
    a numerically-derived depletion width.

    See Also
    --------
    extract_depletion_width_numerical : Extracts W from the numerical
        Poisson E-field solution (equilibrium only, uses 1% threshold).
    """
    params = device_info["params"]
    N_A_ionized = device_info["N_A_ionized"]
    N_D = device_info["N_D"]
    epi_thickness = device_info["epi_thickness_cm"]

    V_bi = built_in_potential(N_A_ionized, N_D, params.n_i_300)
    W = analytical_depletion_width(
        V_bi, V_applied, N_D, eps_r=params.eps_r, epi_thickness=epi_thickness
    )
    return W


def extract_depletion_width_numerical(device_info):
    """Extract depletion width from the numerical equilibrium solution.

    Uses the IntrinsicElectrons profile from the Poisson solver.
    Only valid at equilibrium (0V). Under bias, use extract_depletion_width().

    Parameters
    ----------
    device_info : dict
        Device info dict.

    Returns
    -------
    W : float
        Depletion width (cm), from numerical solution.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    junction_pos = device_info["junction_pos"]
    epi_thickness = device_info["epi_thickness_cm"]
    N_D = device_info["N_D"]

    x_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    n_elec = np.array(
        devsim.get_node_model_values(
            device=device, region=region, name="IntrinsicElectrons"
        )
    )

    n_mask = x_nodes >= junction_pos
    x_n = x_nodes[n_mask]
    n_n = n_elec[n_mask]

    if len(n_n) == 0:
        return 0.0

    # Depletion edge: where n reaches 50% of N_D
    threshold = 0.5 * N_D
    above = n_n >= threshold

    if not np.any(above):
        return epi_thickness

    idx = np.where(above)[0][0]

    # Interpolate
    if idx > 0:
        x0, x1 = x_n[idx - 1], x_n[idx]
        n0, n1 = n_n[idx - 1], n_n[idx]
        if n1 != n0:
            frac = (threshold - n0) / (n1 - n0)
            edge = x0 + frac * (x1 - x0)
        else:
            edge = x_n[idx]
    else:
        edge = x_n[idx]

    W = edge - junction_pos
    return min(max(W, 0.0), epi_thickness)


def voltage_sweep(device_info, voltages=None):
    """High-level function: solve Poisson at multiple bias voltages.

    Solves equilibrium first, then ramps through voltages, extracting
    E-field at each bias and computing depletion width analytically.

    Parameters
    ----------
    device_info : dict
        Device info dict.
    voltages : array_like or None
        Array of bias voltages (V). If None, uses 0 to -60V in -0.5V steps.

    Returns
    -------
    results : dict
        Dictionary with keys:
        - 'voltages': array of solved voltages
        - 'E_fields': list of (x, E) tuples at each voltage
        - 'depletion_widths': array of W values (cm)
        - 'E_max': array of peak |E| values (V/cm)
    """
    if voltages is None:
        voltages = np.arange(0, -60.5, -0.5)

    # Set up and solve equilibrium
    setup_poisson(device_info)
    solve_equilibrium(device_info)

    x0, E0 = extract_electric_field(device_info)
    W0 = extract_depletion_width(device_info, V_applied=0.0)

    solved_voltages = [0.0]
    E_fields = [(x0, E0)]
    depletion_widths = [W0]
    E_max_values = [np.max(np.abs(E0))]

    nonzero_voltages = voltages[voltages != 0]
    if len(nonzero_voltages) > 0:
        V_end = np.min(nonzero_voltages)
        V_step = (
            nonzero_voltages[1] - nonzero_voltages[0]
            if len(nonzero_voltages) > 1
            else -0.5
        )

        ramp_results = ramp_voltage(
            device_info,
            contact_name="cathode",
            V_start=0.0,
            V_end=V_end,
            V_step=V_step,
        )

        for V, x, E in ramp_results:
            W = extract_depletion_width(device_info, V_applied=V)
            solved_voltages.append(V)
            E_fields.append((x, E))
            depletion_widths.append(W)
            E_max_values.append(np.max(np.abs(E)))

    return {
        "voltages": np.array(solved_voltages),
        "E_fields": E_fields,
        "depletion_widths": np.array(depletion_widths),
        "E_max": np.array(E_max_values),
    }
