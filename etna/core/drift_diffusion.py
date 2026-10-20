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

from etna.core.solver_tol import (
    AGGRESSIVE_ABS_ERROR,
    AGGRESSIVE_REL_ERROR,
    DC_ABS_ERROR,
    DC_REL_ERROR,
    RELAXED_ABS_ERROR,
    RELAXED_REL_ERROR,
)
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
from etna.core.poisson import (
    create_ionized_dopant_models,
    region_contacts,
    setup_poisson,
    solve_equilibrium,
)

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
    # opt-in carrier-dependent ionization (audit-v7 A-02): N_D+(Electrons) and
    # N_A-(Holes) replace the static NetDoping in the Poisson charge
    if create_ionized_dopant_models(
        device, region, "Electrons", "Holes", "", ("Electrons", "Holes")
    ):
        charge = "kahan4(Holes, -Electrons, IonizedDonors, -IonizedAcceptors)"

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
    for contact in region_contacts(device_info):
        simple_physics.CreateSiliconDriftDiffusionAtContact(
            device, region, contact, is_circuit=False
        )

    # --- Step 9: Initial coupled solve at equilibrium ---
    logger.info("Performing initial coupled DD solve at equilibrium (0V)")
    devsim.solve(
        type="dc",
        absolute_error=DC_ABS_ERROR,
        relative_error=DC_REL_ERROR,
        maximum_iterations=40,
    )
    logger.info("DD equilibrium solve converged")


def extract_contact_current(device_info, contact="cathode"):
    """Extract the CONDUCTION current at a contact (A/cm^2 for 1D, A/cm for 2D).

    Returns the sum of electron and hole contact currents, read from the continuity
    equations. It does NOT include the displacement current dQ_contact/dt (the
    Poisson / PotentialEquation term, see :func:`extract_contact_charge` and
    :func:`displacement_current`). In DC that term is zero and this is the terminal
    current. In a transient the terminal current is ``conduction + displacement``
    (Mj-4); the time integral of the two agrees once the device is back in its
    initial state (the displacement term telescopes to Q_end - Q_start), so collected
    charge / CCE from this function is right, the pulse SHAPE may not be: the
    displacement share is large when the contact node is depleted (a contact on the
    depleted epi, e.g. the legacy no-n+ stack at full depletion, where the conduction
    current is zero until the carriers arrive) or when a step is not long compared with
    the dielectric relaxation time of the contact layer (~0.1 ps for the 1e18 p+).
    The transient drivers (single_particle, transient) record both components.

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


def extract_contact_charge(device_info, contact="cathode"):
    """Charge on a contact from the Poisson (PotentialEquation) contact equation.

    ``devsim.get_contact_charge(equation="PotentialEquation")``: the electric
    displacement flux into the contact (the ``contactcharge_edge`` model,
    Permittivity * ElectricField, of CreateSiliconPotentialOnlyContact), C/cm^2 in 1D
    and C/cm in 2D. Its time derivative is the displacement part of the terminal
    current (:func:`displacement_current`).
    """
    return devsim.get_contact_charge(
        device=device_info["device_name"], contact=contact, equation="PotentialEquation"
    )


def displacement_current(Q_contact, Q_contact_prev, dt):
    """Displacement current of one implicit (BDF1) step, (Q_k - Q_{k-1}) / dt_k.

    With this backward difference the TOTAL terminal current of a step is
    ``extract_contact_current + displacement_current``, in the sign convention of
    devsim's contact current (Kirchhoff: in 1D the totals at the two contacts sum to
    zero at every step, the conduction parts alone do not). Summed over steps,
    ``sum(I_disp * dt) = Q_end - Q_start`` exactly, so the total and the conduction
    current carry the same charge once the contact charge is back to its start value.
    Verified against the 1D Shockley-Ramo identity (:func:`ramo_current_1d`) and the
    parallel-plate result I = q N v / d (tests/test_displacement_current.py).
    """
    return (Q_contact - Q_contact_prev) / dt


def ramo_current_1d(device_info, contact="anode"):
    """1D Shockley-Ramo terminal current (A/cm^2): (1/L) * integral (J_n + J_p) dx.

    At fixed contact biases, in 1D the weighting field of either contact is 1/L over
    the whole device (contacts at both ends), so the induced (total) terminal current
    is the conduction current density averaged over the device length. This is an
    independent route to ``extract_contact_current + displacement_current``: it uses
    the edge currents of the whole mesh, not the contact equations. The sign is that
    of devsim's contact current at ``contact`` (anode at the left end, x = x_min).
    1D devices only (edge models ``ElectronCurrent`` / ``HoleCurrent`` / ``EdgeLength``).
    """
    device, region = device_info["device_name"], device_info["region_name"]
    J = np.asarray(
        devsim.get_edge_model_values(
            device=device, region=region, name="ElectronCurrent"
        )
    ) + np.asarray(
        devsim.get_edge_model_values(device=device, region=region, name="HoleCurrent")
    )
    length = np.asarray(
        devsim.get_edge_model_values(device=device, region=region, name="EdgeLength")
    )
    x = np.asarray(devsim.get_node_model_values(device=device, region=region, name="x"))
    I = float(np.sum(J * length) / (x.max() - x.min()))
    if contact == "anode":
        return I
    if contact == "cathode":
        return -I
    raise ValueError(f"contact must be 'anode' or 'cathode', got {contact!r}")


# Sign conventions measured on a 2D device (probe_edges.py):
#   interior node residual: sum_{edges n0=i} F - sum_{edges n1=i} F + NodeVolume_i*Gen_i = 0,
#   F = J_edge * EdgeCouple (flux in the n0->n1 direction)
#   devsim.get_contact_current(c) = -(net F from interior into the contact nodes)
# Hence summing over a node set S (no contact nodes):
#   I_c = sum_S vol*gen + sum_{edges leaving S, not into c} F_out
# and with S = all interior nodes this is the old formula I_maj = vol_term - sum_other I_maj(other).

# Continuity equation, its node (generation) model and its edge (current) model per carrier.
_CARRIER_MODELS = {
    "electrons": (
        "ElectronContinuityEquation",
        "ElectronGeneration",
        "ElectronCurrent",
    ),
    "holes": ("HoleContinuityEquation", "HoleGeneration", "HoleCurrent"),
}


def _contact_node_indices(device, region, contact):
    nodes = set()
    for element in devsim.get_element_node_list(
        device=device, region=region, contact=contact
    ):
        nodes.update(int(k) for k in element)
    return nodes


def _global_model_name(parameter, default):
    """Name of the geometry model devsim assembles with (Cartesian default, or the
    cylindrical one while ``_activate_cylindrical_coords`` is in effect)."""
    try:
        return devsim.get_parameter(name=parameter)
    except devsim.error:
        return default


def node_volume_model_name():
    """Node volume model the equations are assembled with ("NodeVolume", or
    "CylindricalNodeVolume" for an axisymmetric device)."""
    return _global_model_name("node_volume_model", "NodeVolume")


def edge_couple_model_name():
    """Edge coupling model (dual edge length in 2D) the equations are assembled with."""
    return _global_model_name("edge_couple_model", "EdgeCouple")


def contacts_on_region(device, region):
    """Contacts of ``device`` that lie on ``region`` (a multi-region device has others)."""
    out = []
    for c in devsim.get_contact_list(device=device):
        if region in devsim.get_region_list(device=device, contact=c):
            out.append(c)
    return out


def contact_majority_carrier(device, region, contact):
    """ "electrons" if the net doping at the contact nodes is n-type, else "holes".

    The majority carrier is the one whose direct contact current is a small difference
    of huge Scharfetter-Gummel terms (see :func:`extract_steady_state_current`); it is
    read from the doping rather than from the contact's name so that a guard ring
    ("guard_ring_anode", p+) or the 3D electrode (n+ column "cathode" on top) are handled.
    """
    nodes = sorted(_contact_node_indices(device, region, contact))
    if not nodes:
        raise ValueError(f"contact {contact!r} has no nodes on region {region!r}")
    net = np.array(
        devsim.get_node_model_values(device=device, region=region, name="NetDoping")
    )[nodes]
    return "electrons" if float(np.mean(net)) > 0 else "holes"


def _edge_node_indices(device, region):
    """(n0, n1) node index arrays of every edge of the region (devsim edge order)."""
    if "node_index@n0" not in devsim.get_edge_model_list(device=device, region=region):
        devsim.edge_from_node_model(
            device=device, region=region, node_model="node_index"
        )
    n0 = np.array(
        devsim.get_edge_model_values(device=device, region=region, name="node_index@n0")
    ).astype(int)
    n1 = np.array(
        devsim.get_edge_model_values(device=device, region=region, name="node_index@n1")
    ).astype(int)
    return n0, n1


def extract_steady_state_current(device_info, contact="cathode", node_set=None):
    """Steady-state terminal current at ``contact`` in the balanced (well-conditioned) form.

    Same quantity as :func:`extract_contact_current` (A/cm^2 in 1D, A/cm in 2D, same
    sign convention: positive = conventional current from the contact into the device) at
    a converged DC state, evaluated without its round-off floor. Works for any number of
    contacts and for axisymmetric devices (the node volume / edge coupling models are the
    ones devsim assembles with).

    devsim's contact current of the contact's MAJORITY carrier (electrons at an n contact,
    holes at a p+ contact; read from NetDoping, :func:`contact_majority_carrier`) is the
    edge flux into the contact nodes, a small difference of two huge SG drift/diffusion
    terms when the carrier density there is high (n ~ N_D at the cathode; with an n+
    layer under the epi, n ~ 1e19). Its absolute noise floor is ~3e-13 A/cm^2 on the
    default 1D device (measured: independent of N_t, random sign), i.e. 2e-3 of the -10 V
    dark current, and it makes the n+ cathode current unusable. Summing the discretised
    continuity equation of that carrier over a set S of non-contact nodes (internal edge
    fluxes cancel) gives the exact discrete identity

        I_maj(contact) = sum_{S} NodeVolume * Generation_maj
                         + sum_{edges leaving S, not into `contact`} F_maj

    with F = J_edge * EdgeCouple the flux leaving S. S must contain every non-contact
    neighbour of ``contact`` (its edges are the unknown). Both terms are well conditioned
    as long as the edges leaving S lie where that carrier is a minority or the region is
    depleted. The minority carrier at ``contact`` is taken directly.

    ``node_set=None`` takes S = all non-contact nodes: the edges leaving S then end on the
    OTHER contacts and their flux sum is minus that contact's direct current, which is
    well conditioned only if the carrier is a minority there. That is the two-terminal
    diode (p+ anode, n cathode) and the guard-ring cathode (anode and guard are both p+),
    but NOT a p+ contact next to another p+ contact (anode vs guard ring: only the SUM of
    their hole currents is well conditioned this way), so that case raises ValueError.

    ``node_set`` (boolean mask over the region nodes, or node indices) chooses S
    explicitly: for the guard-ring anode, the nodes above a cut through the depleted
    n- epi on the anode side of the gap; for the guard, the other side. The cut edges
    carry a minority/depleted-region flux and the split between the two p+ contacts is
    recovered to round-off (test: the three contact currents sum to zero). Contact nodes
    are always removed from S.

    Valid for steady state only (no dn/dt term): transient extraction keeps using
    :func:`extract_contact_current`.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    contacts = contacts_on_region(device, region)
    if contact not in contacts:
        raise ValueError(
            f"contact must be one of {contacts} on region {region!r} (got {contact!r})"
        )
    carrier = contact_majority_carrier(device, region, contact)
    eq_maj, gen_maj, cur_maj = _CARRIER_MODELS[carrier]
    eq_min = _CARRIER_MODELS["holes" if carrier == "electrons" else "electrons"][0]

    contact_nodes = {c: _contact_node_indices(device, region, c) for c in contacts}
    all_contact_nodes = set().union(*contact_nodes.values())
    vol = np.array(
        devsim.get_node_model_values(
            device=device, region=region, name=node_volume_model_name()
        )
    )
    gen = np.array(
        devsim.get_node_model_values(device=device, region=region, name=gen_maj)
    )
    S = np.ones(vol.size, dtype=bool)
    if node_set is not None:
        node_set = np.asarray(node_set)
        if node_set.dtype == bool:
            if node_set.size != vol.size:
                raise ValueError("node_set mask must have one entry per region node")
            S = node_set.copy()
        else:
            S = np.zeros(vol.size, dtype=bool)
            S[node_set.astype(int)] = True
    S[list(all_contact_nodes)] = False
    volume_term = float(np.sum((vol * gen)[S]))

    n0, n1 = _edge_node_indices(device, region)
    own = np.zeros(vol.size, dtype=bool)
    own[list(contact_nodes[contact])] = True
    # every non-contact neighbour of `contact` must be in S
    neighbours = np.concatenate([n1[own[n0]], n0[own[n1]]])
    neighbours = neighbours[~np.isin(neighbours, list(all_contact_nodes))]
    if not np.all(S[neighbours]):
        raise ValueError(
            f"node_set must contain every non-contact neighbour of contact {contact!r}"
        )

    if node_set is None:
        # Edges leaving S end on the other contacts: their flux sum is minus the direct
        # contact current, well conditioned only where the carrier is a minority.
        same = [
            c
            for c in contacts
            if c != contact and contact_majority_carrier(device, region, c) == carrier
        ]
        if same:
            raise ValueError(
                f"contact {contact!r} and {same} share the majority carrier ({carrier}): "
                "the balanced form needs an explicit node_set (cut through the depleted "
                "region separating them)"
            )
        I_maj = volume_term - sum(
            devsim.get_contact_current(device=device, contact=c, equation=eq_maj)
            for c in contacts
            if c != contact
        )
    else:
        J = np.array(
            devsim.get_edge_model_values(device=device, region=region, name=cur_maj)
        )
        couple = np.array(
            devsim.get_edge_model_values(
                device=device, region=region, name=edge_couple_model_name()
            )
        )
        F = J * couple  # flux in the n0 -> n1 direction (A/cm in 2D)
        leaving = S[n0] & ~S[n1] & ~own[n1]  # out of S, not into `contact`
        entering = S[n1] & ~S[n0] & ~own[n0]
        I_maj = volume_term + float(np.sum(F[leaving]) - np.sum(F[entering]))

    I_min = devsim.get_contact_current(device=device, contact=contact, equation=eq_min)
    return I_maj + I_min


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
_RAMP_SOLVE = dict(
    type="dc",
    absolute_error=DC_ABS_ERROR,
    relative_error=DC_REL_ERROR,
    maximum_iterations=40,
)
_RAMP_SOLVE_COARSE = dict(
    type="dc",
    absolute_error=DC_ABS_ERROR,
    relative_error=DC_REL_ERROR,
    maximum_iterations=40,
)
_RAMP_SOLVE_RELAXED = dict(
    type="dc",
    absolute_error=RELAXED_ABS_ERROR,
    relative_error=RELAXED_REL_ERROR,
    maximum_iterations=100,
)
# Third tier, same as solve_equilibrium's. Needed where the n- epi reaches a zero-flux
# surface beyond a p+ implant edge (guard-ring gap, epi overhang, 3D-electrode top
# gap): there the hole continuity equation limit-cycles at RelError ~5e-8 (AbsError
# ~1e2, the normal p+ round-off floor, p ~ 1e-35 in the deep-depleted strip), so
# devsim.solve at rel 1e-8 "fails" between 5 and 48 V depending on the mesh although
# the physics is converged to ~1e-7. A point accepted here is converged to 1e-6
# relative; check it a posteriori with the balanced contact-current sum (should be
# ~0 relative to the generated current).
_RAMP_SOLVE_AGGRESSIVE = dict(
    type="dc",
    absolute_error=AGGRESSIVE_ABS_ERROR,
    relative_error=AGGRESSIVE_REL_ERROR,
    maximum_iterations=200,
)
# Tier order of the fallback chain; the index is what the ramp bookkeeping stores.
_RAMP_TIERS = (
    ("strict", _RAMP_SOLVE),
    ("relaxed", _RAMP_SOLVE_RELAXED),
    ("aggressive", _RAMP_SOLVE_AGGRESSIVE),
)
_AGGRESSIVE = len(_RAMP_TIERS) - 1
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


def _note_aggressive_retry(device_info, V):
    """Log and count a ramp point that only converged with aggressive tolerances."""
    device_info["ramp_aggressive_retries"] = (
        device_info.get("ramp_aggressive_retries", 0) + 1
    )
    logger.warning(
        f"ramp_bias: V={V:.3f}V converged only with aggressive tolerances "
        f"(abs 1e14, rel 1e-6); device {device_info.get('device_name')!r} "
        f"ramp_aggressive_retries={device_info['ramp_aggressive_retries']}"
    )


def _note_tier(device_info, tier, V):
    if tier == 1:
        _note_relaxed_retry(device_info, V)
    elif tier >= 2:
        _note_aggressive_retry(device_info, V)


def _final_errors(result):
    """(relative, absolute) error of the last Newton iteration of a
    ``devsim.solve(info=True)`` result, maximised over the loaded devices."""
    last = result["iterations"][-1]["devices"]
    return (
        max(d["relative_error"] for d in last),
        max(d["absolute_error"] for d in last),
    )


def solve_tier(tier):
    """Run the DC solve of ``_RAMP_TIERS[tier]`` and return the index of the tier the
    resulting state satisfies, or None when it satisfies none.

    A solve that hits its iteration cap with the final update already below the
    aggressive tolerance (relative < 1e-6, absolute < 1e14) is ACCEPTED as an
    aggressive-tier solution (index ``_AGGRESSIVE``): that is the tolerance limit
    cycle of the lateral p+ edges (RelError ~5e-8 for 40 iterations), and the state
    is exactly what the aggressive tier would stop at, without re-solving. Any other
    failure returns None and the caller restores the last converged state.
    """
    name, kw = _RAMP_TIERS[tier]
    try:
        result = devsim.solve(info=True, **kw)
    except devsim.error:
        return None
    if result.get("converged"):
        return tier
    try:
        rel, ab = _final_errors(result)
    except (KeyError, IndexError, ValueError):
        return None
    agg = _RAMP_SOLVE_AGGRESSIVE
    if rel < agg["relative_error"] and ab < agg["absolute_error"]:
        return _AGGRESSIVE
    return None


def _climb_tiers(device_info, bias_name, V, snap, start):
    """Solve the ramp point ``V`` with the tiers from index ``start`` upward, each from
    the last converged state ``snap``. Returns the accepted tier index or None."""
    device = device_info["device_name"]
    region = device_info["region_name"]
    for tier in range(start, len(_RAMP_TIERS)):
        devsim.set_parameter(device=device, name=bias_name, value=V)
        done = solve_tier(tier)
        if done is not None:
            return done
        _restore_solution(device, region, snap)
    return None


def _ramp_bias_adaptive(device_info, bias_name, V_from, V_to, h_min, h_max):
    """Adaptive-step bias ramp (see ``ramp_bias``, ``V_step_max``).

    From the second step on, the Newton start for Potential is extrapolated from
    the two previous converged states. The update of Potential is capped per
    Newton iteration (``log_damp``), so a start that is already shifted by the
    bias step saves most of the iterations (2D 50 um / -50 V: 28 -> ~8 per 2 V).
    Carriers keep the previous solution: extrapolating them too gained nothing.

    Fallback tiers. A strict solve that fails at a step above ``h_min`` halves the
    step (a genuine Newton failure). At ``h_min`` the point is retried with relaxed
    tolerances (_RAMP_SOLVE_RELAXED: absolute 1e12, relative 1e-8), then aggressive
    ones (_RAMP_SOLVE_AGGRESSIVE: absolute 1e14, relative 1e-6, 200 iterations),
    always from the last converged state. A solve that stalls with its final update
    already below the aggressive tolerance is accepted as an aggressive-tier point
    without halving or re-solving (``solve_tier``: the tolerance limit cycle of a
    lateral p+ edge is not cured by a smaller step). Once a point has needed a
    fallback tier, the following INTERMEDIATE points are solved at that tier first
    (the ramp only chooses the path; the target point always starts strict), and
    each such point is counted like a retry. Each use is logged (warning) and
    counted in ``device_info["ramp_relaxed_retries"]`` /
    ``["ramp_aggressive_retries"]``.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    sign = 1.0 if V_to >= V_from else -1.0
    V, h = V_from, h_max
    hist = [(V, _snapshot_solution(device, region))]
    sticky = 0  # tier intermediate points start at (raised once a fallback was needed)
    while abs(V_to - V) > 1e-12:
        step = min(h, abs(V_to - V))
        V_next = V_to if step >= abs(V_to - V) - 1e-12 else V + sign * step
        final = V_next == V_to
        at_h_min = step <= h_min * (1 + 1e-9)
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
        first = 0 if final else sticky
        done = solve_tier(first)
        if done is None:
            _restore_solution(device, region, snap)
            # A device known to need a fallback tier, or a step already minimal: climb
            # the tiers above the one that failed (never below the sticky one).
            if sticky > 0 or at_h_min:
                done = _climb_tiers(
                    device_info, bias_name, V_next, snap, max(first + 1, sticky)
                )
            if done is None:
                if not at_h_min:
                    h = max(h_min, step / 2)
                    continue
                raise RuntimeError(
                    f"ramp_bias: failed to converge at V={V_next:.3f}V "
                    f"(all tolerance tiers)"
                )
        _note_tier(device_info, done, V_next)
        sticky = max(sticky, done)
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
    warning and counted in ``device_info["ramp_relaxed_retries"]`` (absent = none);
    likewise the aggressive third tier (abs 1e14, rel 1e-6, 200 iterations) in
    ``device_info["ramp_aggressive_retries"]``. Both retries start from the last
    converged state; see ``_ramp_bias_adaptive`` for the acceptance of a stalled
    solve whose final update is already within the aggressive tolerance and for the
    sticky fallback tier of intermediate points (the fixed-step ramp does the same).

    Returns
    -------
    success : bool
        True if ramping converged at all steps.

    Raises
    ------
    RuntimeError
        If ramping fails to converge even with aggressive tolerances.
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
        _ramp_bias_adaptive(
            device_info, bias_name, current_V, V_target, V_step, V_step_max
        )
        logger.debug(f"ramp_bias: ramped {contact} to {V_target:.3f}V (adaptive)")
        return True

    n_steps = max(1, int(np.ceil(abs(delta) / V_step)))
    V_intermediates = np.linspace(current_V, V_target, n_steps + 1)[1:]

    region = device_info["region_name"]
    sticky = 0
    for i, V_int in enumerate(V_intermediates):
        final = i == len(V_intermediates) - 1
        snap = _snapshot_solution(device, region)
        devsim.set_parameter(device=device, name=bias_name, value=V_int)
        first = 0 if final else sticky
        done = solve_tier(first)
        if done is None:
            _restore_solution(device, region, snap)
            done = _climb_tiers(
                device_info, bias_name, V_int, snap, max(first + 1, sticky)
            )
            if done is None:
                raise RuntimeError(
                    f"ramp_bias: failed to converge at V={V_int:.3f}V "
                    f"(all tolerance tiers)"
                )
        _note_tier(device_info, done, V_int)
        sticky = max(sticky, done)

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
    at each point and extracting the total contact current (balanced steady-state
    form, see :func:`extract_steady_state_current`; the direct anode edge current is
    round-off at 300 K reverse bias).

    For convergence stability, voltages are applied incrementally from the
    current bias state using small internal steps (V_step_forward for V>0,
    V_step_reverse for V<0). The ramp starts from the contact's live bias
    parameter, as ramp_bias does, so a device already biased (ramp_bias, an
    earlier sweep) is not jumped back from 0 V in one step (audit-v7 A-13).

    Sign: V is the bias of ``contact``. On the default "anode" (p+ side) a
    positive V is forward bias and a negative V reverse bias; on "cathode" the
    signs are the opposite (cv_sweep and ramp_bias work on the cathode).

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
        - "voltages": numpy array of the voltages reached (V)
        - "currents": numpy array of currents (A/cm^2), same length
        - "truncated": True if a ramp step failed (strict and relaxed solve)
          before the end of V_range; the converged prefix is returned
    """
    device = device_info["device_name"]
    V_range = np.asarray(V_range, dtype=float)

    voltages_out = []
    currents_out = []

    # Track the bias state, starting from the live contact bias (A-13)
    bias_name = simple_physics.GetContactBiasName(contact)
    try:
        current_V = float(devsim.get_parameter(device=device, name=bias_name))
    except devsim.error:
        current_V = 0.0

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
                        absolute_error=DC_ABS_ERROR,
                        relative_error=DC_REL_ERROR,
                        maximum_iterations=40,
                    )
                except devsim.error:
                    # Retry with relaxed tolerances
                    try:
                        devsim.solve(
                            type="dc",
                            absolute_error=RELAXED_ABS_ERROR,
                            relative_error=RELAXED_REL_ERROR,
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
                            "truncated": True,
                        }

        # Extract current at the target voltage
        I_total = extract_steady_state_current(device_info, contact=contact)
        voltages_out.append(V_target)
        currents_out.append(I_total)
        current_V = V_target

        logger.debug(f"iv_sweep: V={V_target:.3f}V, I={I_total:.4e} A/cm^2")

    return {
        "voltages": np.array(voltages_out),
        "currents": np.array(currents_out),
        "truncated": False,
    }
