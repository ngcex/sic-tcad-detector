"""Drift-diffusion solver for 4H-SiC p+/n- diode.

Upgrades the Phase 1 Poisson-only solver to coupled Poisson + electron/hole
continuity equations with SRH recombination. Enables current computation
and I-V characterization.

The solver:
1. Assumes Poisson equilibrium already solved (Phase 1 setup)
2. Creates Electrons and Holes solution variables
3. Initializes carriers from equilibrium IntrinsicElectrons/IntrinsicHoles
4. Updates Poisson equation to use actual carriers instead of Boltzmann approx
5. Creates SRH recombination model with proper Jacobian derivatives
6. Creates electron/hole continuity equations with Scharfetter-Gummel currents
7. Sets up drift-diffusion contact equations

All units CGS (cm, cm^-3, A/cm^2, V) per devsim convention.

References:
    - devsim simple_physics.py: CreateSiliconDriftDiffusion pattern
    - devsim simple_dd.py: Bernoulli, CreateElectronCurrent, CreateHoleCurrent
    - Phase 1 clamped exponential approach for SiC n_i stability
"""

from etna.core.solver_tol import DC_REL_ERROR
import logging

import devsim
from devsim.python_packages.model_create import (
    CreateSolution,
    CreateNodeModel,
    CreateNodeModelDerivative,
)
from devsim.python_packages.simple_dd import (
    CreateBernoulli,
    CreateElectronCurrent,
    CreateHoleCurrent,
)
import devsim.python_packages.simple_physics as simple_physics
import numpy as np

from etna.core.device import create_sic_device
from etna.core.poisson import setup_poisson, solve_equilibrium

logger = logging.getLogger(__name__)


def setup_sic_drift_diffusion(device_info):
    """Set up full drift-diffusion equations for 4H-SiC diode.

    Upgrades from Poisson-only to coupled Poisson + electron/hole continuity
    equations with SRH recombination. Must be called after setup_poisson()
    and solve_equilibrium() have been run.

    Parameters
    ----------
    device_info : dict
        Device info dict from create_sic_device(). Must have Poisson
        equilibrium already solved.

    Notes
    -----
    Following the devsim simple_physics.CreateSiliconDriftDiffusion pattern,
    adapted for 4H-SiC with clamped exponentials and SiC material parameters.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]

    # --- Step 1: Create carrier solution variables ---
    CreateSolution(device, region, "Electrons")
    CreateSolution(device, region, "Holes")

    # --- Step 2: Initialize carriers from equilibrium Poisson solution ---
    # IntrinsicElectrons and IntrinsicHoles computed by _create_sic_potential_only
    # with clamped exponentials for SiC stability
    init_n = devsim.get_node_model_values(
        device=device, region=region, name="IntrinsicElectrons"
    )
    init_p = devsim.get_node_model_values(
        device=device, region=region, name="IntrinsicHoles"
    )
    devsim.set_node_values(
        device=device, region=region, name="Electrons", values=init_n
    )
    devsim.set_node_values(device=device, region=region, name="Holes", values=init_p)

    # --- Step 3: Update Poisson to use actual Electrons/Holes ---
    # Replace Boltzmann approximation with actual carrier concentrations
    charge = "kahan3(Holes, -Electrons, NetDoping)"
    pcharge = "-ElectronCharge * IntrinsicCharge2"

    CreateNodeModel(device, region, "IntrinsicCharge2", charge)
    CreateNodeModelDerivative(device, region, "IntrinsicCharge2", charge, "Electrons")
    CreateNodeModelDerivative(device, region, "IntrinsicCharge2", charge, "Holes")

    CreateNodeModel(device, region, "PotentialIntrinsicCharge2", pcharge)
    CreateNodeModelDerivative(
        device, region, "PotentialIntrinsicCharge2", pcharge, "Electrons"
    )
    CreateNodeModelDerivative(
        device, region, "PotentialIntrinsicCharge2", pcharge, "Holes"
    )

    # Update Poisson equation to reference actual carriers
    devsim.equation(
        device=device,
        region=region,
        name="PotentialEquation",
        variable_name="Potential",
        node_model="PotentialIntrinsicCharge2",
        edge_model="PotentialEdgeFlux",
        variable_update="log_damp",
    )

    # --- Step 4: SRH recombination ---
    # USRH = (n*p - ni^2) / (taup*(n + n1) + taun*(p + p1))
    # n1 = n_i, p1 = n_i for midgap traps (already set as parameters)
    USRH = (
        "(Electrons * Holes - n_i^2) / "
        "(taup * (Electrons + n1) + taun * (Holes + p1))"
    )
    CreateNodeModel(device, region, "USRH", USRH)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "USRH", USRH, var)

    # Generation terms for continuity equations
    Gn = "-ElectronCharge * USRH"
    Gp = "+ElectronCharge * USRH"
    CreateNodeModel(device, region, "ElectronGeneration", Gn)
    CreateNodeModel(device, region, "HoleGeneration", Gp)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn, var)
        CreateNodeModelDerivative(device, region, "HoleGeneration", Gp, var)

    # --- Step 5: Bernoulli function and carrier currents ---
    CreateBernoulli(device, region)
    CreateElectronCurrent(device, region, "mu_n")
    CreateHoleCurrent(device, region, "mu_p")

    # --- Step 6: Electron continuity equation ---
    NCharge = "-ElectronCharge * Electrons"
    CreateNodeModel(device, region, "NCharge", NCharge)
    CreateNodeModelDerivative(device, region, "NCharge", NCharge, "Electrons")

    devsim.equation(
        device=device,
        region=region,
        name="ElectronContinuityEquation",
        variable_name="Electrons",
        time_node_model="NCharge",
        edge_model="ElectronCurrent",
        variable_update="positive",
        node_model="ElectronGeneration",
    )

    # --- Step 7: Hole continuity equation ---
    PCharge = "ElectronCharge * Holes"
    CreateNodeModel(device, region, "PCharge", PCharge)
    CreateNodeModelDerivative(device, region, "PCharge", PCharge, "Holes")

    devsim.equation(
        device=device,
        region=region,
        name="HoleContinuityEquation",
        variable_name="Holes",
        time_node_model="PCharge",
        edge_model="HoleCurrent",
        variable_update="positive",
        node_model="HoleGeneration",
    )

    # --- Step 8: Contact equations for drift-diffusion ---
    for contact in ("anode", "cathode"):
        simple_physics.CreateSiliconDriftDiffusionAtContact(
            device, region, contact, is_circuit=False
        )

    # --- Step 9: Initial coupled solve at equilibrium ---
    logger.info("Performing initial coupled DD solve at equilibrium (0V)")
    devsim.solve(
        type="dc",
        absolute_error=1e10,
        relative_error=DC_REL_ERROR,
        maximum_iterations=40,
    )
    logger.info("DD equilibrium solve converged")


def extract_contact_current(device_info, contact="cathode"):
    """Extract total current at a contact (A/cm^2 for 1D).

    Returns the sum of electron and hole contact currents.

    Parameters
    ----------
    device_info : dict
        Device info dict.
    contact : str
        Contact name ("anode" or "cathode").

    Returns
    -------
    I_total : float
        Total contact current (A/cm^2). Positive = conventional current
        flowing into the device at the contact.
    """
    device = device_info["device_name"]
    I_e = devsim.get_contact_current(
        device=device, contact=contact, equation="ElectronContinuityEquation"
    )
    I_h = devsim.get_contact_current(
        device=device, contact=contact, equation="HoleContinuityEquation"
    )
    return I_e + I_h


def create_dd_device(doping_profile="graded", **kwargs):
    """Convenience function to create a device with full DD setup.

    Creates device, sets up Poisson, solves equilibrium, and sets up
    drift-diffusion equations in one call.

    Parameters
    ----------
    doping_profile : str
        "uniform" or "graded". Default "graded".
    **kwargs
        Additional keyword arguments passed to create_sic_device()
        (e.g., device_name, N_D_junction, N_D_bulk, L_transition, etc.).

    Returns
    -------
    device_info : dict
        Device info dict augmented with dd_initialized=True.
    """
    device_info = create_sic_device(doping_profile=doping_profile, **kwargs)

    setup_poisson(device_info)
    solve_equilibrium(device_info)

    setup_sic_drift_diffusion(device_info)

    device_info["dd_initialized"] = True
    return device_info


# _RAMP_SOLVE_COARSE (intermediate ramp points) is currently IDENTICAL to _RAMP_SOLVE
# (target point): same tolerance and the same 40-iteration cap. Kept as a separate name
# so the two can be tuned independently.
_RAMP_SOLVE = dict(type="dc", absolute_error=1e10, relative_error=DC_REL_ERROR, maximum_iterations=40)
_RAMP_SOLVE_COARSE = dict(
    type="dc", absolute_error=1e10, relative_error=DC_REL_ERROR, maximum_iterations=40
)
_RAMP_SOLVE_RELAXED = dict(
    type="dc", absolute_error=1e12, relative_error=1e-8, maximum_iterations=100
)
_RAMP_MODELS = ("Potential", "Electrons", "Holes")


def _snapshot_solution(device, region):
    return {
        m: devsim.get_node_model_values(device=device, region=region, name=m)
        for m in _RAMP_MODELS
    }


def _restore_solution(device, region, snap):
    for m, values in snap.items():
        devsim.set_node_values(device=device, region=region, name=m, values=values)


def _note_relaxed_retry(device_info, V):
    """Log and count a ramp point that only converged with relaxed tolerances."""
    device_info["ramp_relaxed_retries"] = device_info.get("ramp_relaxed_retries", 0) + 1
    logger.warning(
        f"ramp_bias: V={V:.3f}V converged only with relaxed tolerances "
        f"(abs 1e12, rel 1e-8); device {device_info.get('device_name')!r} "
        f"ramp_relaxed_retries={device_info['ramp_relaxed_retries']}"
    )


def _ramp_bias_adaptive(device_info, bias_name, V_from, V_to, h_min, h_max):
    """Adaptive-step bias ramp (see ``ramp_bias``, ``V_step_max``).

    From the second step on, the Newton start for Potential is extrapolated from
    the two previous converged states. The update of Potential is capped per
    Newton iteration (``log_damp``), so a start that is already shifted by the
    bias step saves most of the iterations (2D 50 um / -50 V: 28 -> ~8 per 2 V).
    Carriers keep the previous solution: extrapolating them too gained nothing.

    When a step at the minimum size ``h_min`` still fails, the point is retried once
    with relaxed tolerances (_RAMP_SOLVE_RELAXED: absolute 1e12, relative 1e-8). Such a
    point is less tightly converged; each use is logged (warning) and counted in
    ``device_info["ramp_relaxed_retries"]``.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    sign = 1.0 if V_to >= V_from else -1.0
    V, h = V_from, h_max
    hist = [(V, _snapshot_solution(device, region))]
    while abs(V_to - V) > 1e-12:
        step = min(h, abs(V_to - V))
        V_next = V_to if step >= abs(V_to - V) - 1e-12 else V + sign * step
        final = V_next == V_to
        snap = hist[-1][1]
        if len(hist) >= 2:
            (V1, s1), (V0, s0) = hist[-1], hist[-2]
            ratio = (V_next - V1) / (V1 - V0)
            psi1, psi0 = np.asarray(s1["Potential"]), np.asarray(s0["Potential"])
            devsim.set_node_values(
                device=device,
                region=region,
                name="Potential",
                values=(psi1 + (psi1 - psi0) * ratio).tolist(),
            )
        devsim.set_parameter(device=device, name=bias_name, value=V_next)
        try:
            # _RAMP_SOLVE and _RAMP_SOLVE_COARSE are currently identical (same tolerance
            # and iteration cap at the target and at intermediate points).
            devsim.solve(**(_RAMP_SOLVE if final else _RAMP_SOLVE_COARSE))
        except devsim.error:
            _restore_solution(device, region, snap)
            if step > h_min * (1 + 1e-9):
                h = max(h_min, step / 2)
                continue
            devsim.set_parameter(device=device, name=bias_name, value=V_next)
            try:
                devsim.solve(**_RAMP_SOLVE_RELAXED)
            except devsim.error as e:
                raise RuntimeError(
                    f"ramp_bias: failed to converge at V={V_next:.3f}V: {e}"
                )
            _note_relaxed_retry(device_info, V_next)
        V = V_next
        hist = [hist[-1], (V, _snapshot_solution(device, region))]
        h = min(h * 2, h_max)


def ramp_bias(device_info, V_target, contact="anode", V_step=0.5, V_step_max=None):
    """Ramp contact bias to target voltage incrementally.

    Similar to what iv_sweep does internally but without collecting I-V data.
    Used to establish a bias point before adding generation (for convergence
    stability).

    Parameters
    ----------
    device_info : dict
        Device info dict with DD equations already set up (dd_initialized=True).
    V_target : float
        Target voltage to ramp to (V). Negative = reverse bias for anode contact.
    contact : str
        Contact at which to apply bias ("anode" or "cathode").
    V_step : float
        Maximum voltage step for ramping (V).
    V_step_max : float, optional
        If given, the ramp is adaptive: it starts with steps of ``V_step_max``,
        doubles after each converged step, and on a Newton failure restores the
        last converged solution and halves the step down to ``V_step``, where the
        relaxed-tolerance retry below applies. Intermediate and target points use
        the same tolerance (``DC_REL_ERROR``) and iteration cap. The ramp only chooses the path to
        the target; the converged state at ``V_target`` is the same (to solver
        tolerance). ``None`` keeps fixed steps of ``V_step``.

    Every point that converges only on the relaxed-tolerance retry is logged as a
    warning and counted in ``device_info["ramp_relaxed_retries"]`` (absent = none).

    Returns
    -------
    success : bool
        True if ramping converged at all steps.

    Raises
    ------
    RuntimeError
        If ramping fails to converge even with relaxed tolerances.
    """
    device = device_info["device_name"]
    bias_name = simple_physics.GetContactBiasName(contact)

    # Get current bias value
    try:
        current_V = devsim.get_parameter(device=device, name=bias_name)
    except devsim.error:
        current_V = 0.0

    delta = V_target - current_V
    if abs(delta) < 1e-12:
        return True

    if V_step_max is not None and V_step_max > V_step:
        _ramp_bias_adaptive(device_info, bias_name, current_V, V_target, V_step, V_step_max)
        logger.debug(f"ramp_bias: ramped {contact} to {V_target:.3f}V (adaptive)")
        return True

    n_steps = max(1, int(np.ceil(abs(delta) / V_step)))
    V_intermediates = np.linspace(current_V, V_target, n_steps + 1)[1:]

    for V_int in V_intermediates:
        devsim.set_parameter(device=device, name=bias_name, value=V_int)
        try:
            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=DC_REL_ERROR,
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
                _note_relaxed_retry(device_info, V_int)
            except devsim.error as e:
                raise RuntimeError(
                    f"ramp_bias: failed to converge at V={V_int:.3f}V: {e}"
                )

    logger.debug(f"ramp_bias: ramped {contact} to {V_target:.3f}V")
    return True


def iv_sweep(
    device_info,
    V_range,
    contact="anode",
    V_step_forward=0.1,
    V_step_reverse=0.5,
):
    """Sweep voltage and extract I-V curve from drift-diffusion solution.

    Ramps the contact bias through V_range, solving the coupled DD equations
    at each point and extracting the total contact current.

    For convergence stability, voltages are applied incrementally from the
    current bias state using small internal steps (V_step_forward for V>0,
    V_step_reverse for V<0).

    Parameters
    ----------
    device_info : dict
        Device info dict with DD equations already set up (dd_initialized=True).
    V_range : array_like
        Array of voltages to sweep (V). Should be sorted.
    contact : str
        Contact at which to extract current ("anode" or "cathode").
    V_step_forward : float
        Maximum voltage step for forward bias ramping (V).
    V_step_reverse : float
        Maximum voltage step for reverse bias ramping (V).

    Returns
    -------
    result : dict
        Dictionary with:
        - "voltages": numpy array of applied voltages (V)
        - "currents": numpy array of currents (A/cm^2)
    """
    device = device_info["device_name"]
    V_range = np.asarray(V_range, dtype=float)

    voltages_out = []
    currents_out = []

    # Track current bias state
    current_V = 0.0
    bias_name = simple_physics.GetContactBiasName(contact)

    for V_target in V_range:
        # Determine step size based on direction
        if V_target >= 0:
            step = V_step_forward
        else:
            step = V_step_reverse

        # Ramp from current_V to V_target in small steps
        delta = V_target - current_V
        if abs(delta) < 1e-12:
            # Already at target, just extract current
            pass
        else:
            n_steps = max(1, int(np.ceil(abs(delta) / step)))
            V_intermediates = np.linspace(current_V, V_target, n_steps + 1)[1:]

            for V_int in V_intermediates:
                devsim.set_parameter(device=device, name=bias_name, value=V_int)
                try:
                    devsim.solve(
                        type="dc",
                        absolute_error=1e10,
                        relative_error=DC_REL_ERROR,
                        maximum_iterations=40,
                    )
                except devsim.error:
                    # Retry with relaxed tolerances
                    try:
                        devsim.solve(
                            type="dc",
                            absolute_error=1e12,
                            relative_error=1e-8,
                            maximum_iterations=100,
                        )
                        logger.info(
                            f"iv_sweep: converged at V={V_int:.3f}V with relaxed tolerances"
                        )
                    except devsim.error as e:
                        logger.warning(
                            f"iv_sweep: failed to converge at V={V_int:.3f}V: {e}"
                        )
                        # Return what we have so far
                        return {
                            "voltages": np.array(voltages_out),
                            "currents": np.array(currents_out),
                        }

        # Extract current at the target voltage
        I_total = extract_contact_current(device_info, contact=contact)
        voltages_out.append(V_target)
        currents_out.append(I_total)
        current_V = V_target

        logger.debug(f"iv_sweep: V={V_target:.3f}V, I={I_total:.4e} A/cm^2")

    return {
        "voltages": np.array(voltages_out),
        "currents": np.array(currents_out),
    }
