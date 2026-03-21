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

from src.device import create_sic_device
from src.poisson import setup_poisson, solve_equilibrium

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
        relative_error=1e-10,
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
