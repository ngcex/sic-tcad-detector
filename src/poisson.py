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
from src.sic_material import intrinsic_concentration

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

    Uses a three-stage fallback: tight → relaxed → aggressive tolerances.
    The aggressive stage is needed for 2D meshes in parametric sweeps where
    prior device state may affect convergence.

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
        except devsim.error:
            try:
                devsim.solve(
                    type="dc",
                    absolute_error=1e14,
                    relative_error=1e-6,
                    maximum_iterations=200,
                )
                logger.info("Equilibrium solve converged (aggressive)")
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

    # Use T-dependent n_i for V_bi calculation
    T = device_info.get("T", 300)
    n_i_T = intrinsic_concentration(T, params)[0]
    V_bi = built_in_potential(N_A_ionized, N_D, n_i_T)
    W = analytical_depletion_width(
        V_bi, V_applied, N_D, eps_r=params.eps_r, epi_thickness=epi_thickness
    )
    return W


def extract_depletion_width_numerical(device_info):
    """Extract depletion width from the numerical solution.

    Compares the carrier concentration (IntrinsicElectrons or Electrons if
    available from DD solver) against the local Donors profile to find
    where the n-side becomes quasi-neutral. Works for both uniform and
    graded doping profiles, and under bias.

    For graded doping under reverse bias, the Poisson-only solver may create
    a non-physical second depletion region near the cathode. This function
    finds only the first depletion-to-neutral transition from the junction
    side, ignoring any cathode-side artifacts.

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

    x_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )

    # Prefer Electrons (DD solution) over IntrinsicElectrons (Poisson-only)
    try:
        n_elec = np.array(
            devsim.get_node_model_values(device=device, region=region, name="Electrons")
        )
    except Exception:
        n_elec = np.array(
            devsim.get_node_model_values(
                device=device, region=region, name="IntrinsicElectrons"
            )
        )

    # Get local donor concentration for position-dependent threshold
    donors = np.array(
        devsim.get_node_model_values(device=device, region=region, name="Donors")
    )

    n_mask = x_nodes >= junction_pos
    x_n = x_nodes[n_mask]
    n_n = n_elec[n_mask]
    donors_n = donors[n_mask]

    if len(n_n) == 0:
        return 0.0

    # Depletion edge: where n reaches 50% of local Donors (position-dependent)
    # Use ratio n/Donors to handle graded doping correctly
    ratio = np.where(donors_n > 0, n_n / donors_n, 0.0)

    # Find first point where ratio >= 0.5 (junction-side depletion-to-neutral)
    above = ratio >= 0.5

    if not np.any(above):
        return epi_thickness

    idx = np.where(above)[0][0]

    # Interpolate
    if idx > 0:
        x0, x1 = x_n[idx - 1], x_n[idx]
        r0, r1 = ratio[idx - 1], ratio[idx]
        if r1 != r0:
            frac = (0.5 - r0) / (r1 - r0)
            edge = x0 + frac * (x1 - x0)
        else:
            edge = x_n[idx]
    else:
        edge = x_n[idx]

    W = edge - junction_pos
    return min(max(W, 0.0), epi_thickness)


def extract_depletion_width_2d_center(device_info, center_x_tol_cm=1e-6):
    """Extract the center-column depletion width from a 2D device.

    The 1D ``extract_depletion_width_numerical`` reads ``x`` as the depth axis
    (1D convention) and therefore returns the WRONG W for a 2D device, where
    ``x`` is lateral and ``y`` is depth. This function mirrors the 1D
    depletion-edge logic but along the ``y`` axis, restricted to the
    symmetry-plane center column (``|x| < center_x_tol_cm``, since ``x = 0`` is
    the symmetry plane in the 2D modules).

    The depletion width is measured from the metallurgical junction
    (``y = junction_pos``) into the n-epi (``y > junction_pos``), as the first
    point where the electron concentration recovers to 50% of the local donor
    concentration (handling the graded profile correctly).

    Parameters
    ----------
    device_info : dict
        2D device info dict. Must have ``dimension == 2`` and the keys
        ``device_name``, ``region_name``, ``junction_pos``, ``epi_thickness_cm``.
    center_x_tol_cm : float
        Half-width of the center-column selection band around ``x = 0`` (cm).

    Returns
    -------
    W : float
        Center-column depletion width (cm), clamped to [0, epi_thickness].

    Raises
    ------
    ValueError
        If ``device_info["dimension"]`` is not 2, or if no nodes fall within
        ``center_x_tol_cm`` of the symmetry plane.
    """
    dimension = device_info.get("dimension")
    if dimension != 2:
        raise ValueError(
            "extract_depletion_width_2d_center requires a 2D device "
            f"(dimension == 2); got dimension={dimension!r}"
        )

    device = device_info["device_name"]
    region = device_info["region_name"]
    junction_pos = device_info["junction_pos"]
    epi_thickness = device_info["epi_thickness_cm"]

    x_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    y_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="y")
    )
    donors = np.array(
        devsim.get_node_model_values(device=device, region=region, name="Donors")
    )

    # Prefer Electrons (DD solution) over IntrinsicElectrons (Poisson-only),
    # mirroring the 1D extractor.
    try:
        n_elec = np.array(
            devsim.get_node_model_values(device=device, region=region, name="Electrons")
        )
    except Exception:
        n_elec = np.array(
            devsim.get_node_model_values(
                device=device, region=region, name="IntrinsicElectrons"
            )
        )

    # Select the center column (symmetry plane at x = 0).
    center_mask = np.abs(x_nodes) < center_x_tol_cm
    if not np.any(center_mask):
        raise ValueError(
            "no nodes within center_x_tol_cm of the symmetry plane "
            f"(min |x| = {np.min(np.abs(x_nodes)):.3e} cm, "
            f"tol = {center_x_tol_cm:.3e} cm)"
        )

    y_c = y_nodes[center_mask]
    n_c = n_elec[center_mask]
    d_c = donors[center_mask]

    # Sort by depth (y ascending).
    order = np.argsort(y_c)
    y_c = y_c[order]
    n_c = n_c[order]
    d_c = d_c[order]

    # Restrict to the n-side (y >= junction_pos).
    n_mask = y_c >= junction_pos
    y_n = y_c[n_mask]
    n_n = n_c[n_mask]
    donors_n = d_c[n_mask]

    if len(n_n) == 0:
        return 0.0

    # Depletion edge: first point where n recovers to 50% of local Donors.
    ratio = np.where(donors_n > 0, n_n / donors_n, 0.0)
    above = ratio >= 0.5

    if not np.any(above):
        return float(epi_thickness)

    idx = np.where(above)[0][0]

    if idx > 0:
        y0, y1 = y_n[idx - 1], y_n[idx]
        r0, r1 = ratio[idx - 1], ratio[idx]
        if r1 != r0:
            frac = (0.5 - r0) / (r1 - r0)
            edge = y0 + frac * (y1 - y0)
        else:
            edge = y_n[idx]
    else:
        edge = y_n[idx]

    W = float(edge - junction_pos)
    return min(max(W, 0.0), float(epi_thickness))


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
