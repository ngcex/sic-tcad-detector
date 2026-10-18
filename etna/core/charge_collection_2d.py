"""2D charge collection efficiency (CCE) computation utilities.

Provides:
- Area integration over 2D triangular mesh (replaces 1D trapezoid)
- 2D CCE computation from DD simulation
- Lateral CCE scanning (center to edge)
- 2D CCE heatmap generation
- 2D-vs-1D CCE comparison with active volume ratio

All units in CGS (cm, cm^-3, V, s, A/cm) per devsim convention.

Note on symmetry:
    The 2D device represents a half-device (x=0 to x=half_width) by
    exploiting mirror symmetry at x=0.  Both I_collected (contact current)
    and I_generated (volume integral of generation) are computed over the
    half-device, so the factor of 2 cancels in the CCE ratio.  The CCE
    values are therefore identical to what a full-device simulation would
    produce, with half the computational cost.

References:
    - Phase 19: 2D mesh & Poisson validation
    - Phase 20 research: 2D drift-diffusion CCE methodology
    - reference experimental microdosimeter geometry
"""

import logging
import uuid

import devsim
import numpy as np

from etna.core.device2d import create_sic_2d_device
from etna.core.poisson import setup_poisson, solve_equilibrium
from etna.core.drift_diffusion import (
    _AGGRESSIVE,
    _RAMP_TIERS,
    _contact_node_indices,
    setup_sic_drift_diffusion,
    contacts_on_region,
    extract_steady_state_current,
    node_volume_model_name,
    _restore_solution,
    _snapshot_solution,
    ramp_bias,
    solve_tier,
)
from etna.core.charge_collection import add_generation_to_dd
from etna.core.generation_profiles import alpha_depth_shape, confine_generation_to_epi
from etna.core.devsim_reset import delete_device_and_mesh

logger = logging.getLogger(__name__)

# Physical constants (CGS)
Q = 1.602e-19  # C, elementary charge

# Default alpha particle range in SiC
_ALPHA_RANGE_CM = 15e-4  # 15 um


def in_mesh_generation_fraction(
    epi_thickness_cm, alpha_range_cm=_ALPHA_RANGE_CM, n=20001
):
    """Fraction of the alpha's e-h pairs (depth profile) generated inside the simulated epi.

    The lateral scans generate pairs only between the junction and the cathode, i.e.
    inside the epi; the alpha_depth_shape profile extends to ~alpha_range (15 um) and
    beyond, so for a thinner epi part of the deposit is outside the mesh and is not
    counted in I_generated. CCE is therefore relative to in-epi generation. This
    returns int_0^epi shape / int_0^inf shape (1D quadrature on a dense grid, depth
    only: the lateral Gaussian stripe, which is also cut at x = half_width, is ignored).
    ~0.42 / 0.83 / 1.00 for 5 / 10 / 20 um epi at the default 15 um range.
    """
    y_end = max(3.0 * alpha_range_cm, epi_thickness_cm)
    y = np.linspace(0.0, y_end, n)
    total = np.trapezoid(alpha_depth_shape(y, alpha_range_cm=alpha_range_cm), y)
    if total <= 0:
        return 1.0
    y_in = np.linspace(0.0, min(epi_thickness_cm, y_end), n)
    inside = np.trapezoid(alpha_depth_shape(y_in, alpha_range_cm=alpha_range_cm), y_in)
    return float(min(1.0, inside / total))


def _robust_dc_solve(device_info=None):
    """Solve DC with fallback to relaxed, then aggressive tolerances.

    Same three tiers as ``solve_equilibrium`` and ``ramp_bias`` (strict: abs 1e10 /
    rel DC_REL_ERROR / 40 it; relaxed: 1e12 / 1e-8 / 100; aggressive: 1e14 / 1e-6 /
    200), with ``drift_diffusion.solve_tier``'s acceptance of a stalled solve whose
    final update is already within the aggressive tolerance (the lateral p+ edge
    limit cycle). The third tier's use is logged as a warning, and a state accepted
    there should be checked with the balanced contact-current sum
    (:func:`contact_current_sharing`).

    ``device_info`` is optional: a device whose bias ramp already needed the
    aggressive tier (``device_info["ramp_aggressive_retries"]``) is solved at that
    tier directly, which is what the strict attempt would end up accepting anyway
    after 40 cycling iterations. Raises devsim.error when no tier converges.
    """
    start = 0
    if device_info is not None and device_info.get("ramp_aggressive_retries", 0) > 0:
        start = _AGGRESSIVE
    for tier in range(start, len(_RAMP_TIERS)):
        done = solve_tier(tier)
        if done is not None:
            break
    else:
        raise devsim.error("DC solve failed to converge (all tolerance tiers)")
    if done == 1:
        logger.info("DC solve converged with relaxed tolerances")
    elif done >= 2:
        logger.warning(
            "DC solve converged only with aggressive tolerances (abs 1e14, rel 1e-6)"
        )


def integrate_nodes_box_2d(device_info, node_values):
    """Integrate a per-node field with devsim's box rule: sum(NodeVolume * value).

    This is the quadrature the solver itself applies to a node source term, so
    I_generated and the collected current are consistent. The triangle rule of
    ``integrate_over_mesh_2d`` differs from it by O(h_x h_y dG/dx dG/dy), which
    made the edge/centre CCE ratio a mesh artefact (audit PERF-03). The volume model
    is the one the equations are assembled with (``CylindricalNodeVolume`` on the
    axisymmetric 3D-electrode device, where the Cartesian ``NodeVolume`` would weight
    a shell at radius r as if it had no circumference).
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    volume = np.array(
        devsim.get_node_model_values(
            device=device, region=region, name=node_volume_model_name()
        )
    )
    return float(np.sum(volume * np.asarray(node_values, dtype=float)))


def integrate_over_mesh_2d(device_info, node_values):
    """Integrate a per-node scalar field over the 2D triangular mesh.

    Uses vectorized triangle-area integration: for each triangle, compute
    the area via cross product and multiply by the average of the three
    vertex values.

    Parameters
    ----------
    device_info : dict
        Device info dict from create_sic_2d_device.
    node_values : array_like
        Scalar field values at each mesh node.  Length must match the
        number of nodes in the device region.

    Returns
    -------
    integral : float
        Area integral of the field over the mesh (units: [field] * cm^2).
    """
    device = device_info["device_name"]
    region = device_info["region_name"]

    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    y = np.array(devsim.get_node_model_values(device=device, region=region, name="y"))
    node_vals = np.asarray(node_values, dtype=float)

    # Get element (triangle) connectivity -- flat list, reshape to (N_tri, 3)
    elements = devsim.get_element_node_list(device=device, region=region)
    triangles = np.array(elements, dtype=int).reshape(-1, 3)

    # Vertex coordinates for all triangles
    x0, x1, x2 = x[triangles[:, 0]], x[triangles[:, 1]], x[triangles[:, 2]]
    y0, y1, y2 = y[triangles[:, 0]], y[triangles[:, 1]], y[triangles[:, 2]]

    # Triangle areas via cross product: 0.5 * |det|
    areas = 0.5 * np.abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))

    # Average node values per triangle
    v0 = node_vals[triangles[:, 0]]
    v1 = node_vals[triangles[:, 1]]
    v2 = node_vals[triangles[:, 2]]
    v_avg = (v0 + v1 + v2) / 3.0

    return float(np.sum(v_avg * areas))


def create_2d_dd_device(
    half_width_um=50.0, V_bias=50.0, device_name=None, **device_kwargs
):
    """Create a 2D device with full drift-diffusion setup and bias applied.

    Convenience wrapper that calls create_sic_2d_device, sets up Poisson,
    solves equilibrium, sets up DD, and ramps cathode to reverse bias.

    Parameters
    ----------
    half_width_um : float
        Half-width of sensitive volume in micrometers (default 50 um for
        100 um SV).
    V_bias : float
        Reverse bias voltage to apply at cathode (V, positive value).
    device_name : str or None
        Device name.  If None, generates a unique name via uuid4.
    **device_kwargs
        Additional keyword arguments forwarded to create_sic_2d_device
        (e.g., epi_thickness_cm, N_D_bulk, N_D_junction, L_transition, T).

    Returns
    -------
    device_info : dict
        Device info dict with dd_initialized=True and bias applied.
    """
    if device_name is None:
        device_name = f"cce2d_{uuid.uuid4().hex[:8]}"

    device_info = create_sic_2d_device(
        device_name=device_name,
        half_width_um=half_width_um,
        **device_kwargs,
    )
    setup_poisson(device_info)
    solve_equilibrium(device_info)
    setup_sic_drift_diffusion(device_info)
    device_info["dd_initialized"] = True

    # Ramp cathode to reverse bias (positive V on cathode = reverse bias
    # for p+/n- diode, same convention as Phase 19 validation)
    ramp_bias(
        device_info, V_target=V_bias, contact="cathode", V_step=0.5, V_step_max=5.0
    )
    logger.info(
        f"Created 2D DD device '{device_name}': "
        f"half_width={half_width_um:.0f} um, V_bias={V_bias:.1f} V"
    )

    return device_info


def injected_generation(device_info, generation_values):
    """``generation_values`` with the contact nodes zeroed: the generation the solver
    actually injects. On a contact node the contact equation replaces the continuity
    equation (Dirichlet carrier densities), so a source term assigned there never enters
    the device."""
    device = device_info["device_name"]
    region = device_info["region_name"]
    g = np.array(generation_values, dtype=float, copy=True)
    for c in contacts_on_region(device, region):
        g[sorted(_contact_node_indices(device, region, c))] = 0.0
    return g


def compute_cce_2d(device_info, generation_values, contact="cathode", node_set=None):
    """Compute 2D charge collection efficiency from DD simulation.

    CCE = |I_collected| / I_generated, where I_generated is the box-rule integral
    of the generation rate over the NON-CONTACT nodes (``injected_generation``): the
    generation the solver injects. Since 2026-10-18; before, the integral ran over all
    nodes and counted the never-injected contact-node generation, a mesh-dependent CCE
    deficit (1.26e-3 on the legacy planar diode, 50 nm mesh, at the cathode; 0 for an
    alpha stripe on the n+ stack, whose contact nodes carry no generation; 6.1e-4 for a
    uniform generation through the n+).

    The symmetry factor of 2 (half-device) cancels in the ratio since
    both I_collected and I_generated are for the half-device.

    The contact currents devsim reports are those of the LAST solve: the caller must
    inject ``generation_values`` (``add_generation_to_dd``) and solve before calling
    this, or use :func:`cce_with_generation`, which does both and restores the
    device afterwards.

    Parameters
    ----------
    device_info : dict
        Device info dict with DD solved and generation added.
    generation_values : array_like
        Generation rate G(x,y) at each mesh node (cm^-3 s^-1).
    contact : str
        Contact at which to extract collected current.
    node_set : array_like, optional
        Passed to :func:`extract_steady_state_current` (needed for a p+ contact next
        to another p+ contact, e.g. the guard-ring anode).

    Returns
    -------
    cce : float
        Charge collection efficiency, clipped to [0, 1].
    """
    # I_generated = Q * integral(G dA) over the injected (non-contact) nodes [A/cm for 2D]
    I_generated = Q * integrate_nodes_box_2d(
        device_info, injected_generation(device_info, generation_values)
    )

    # I_collected at contact [A/cm for 2D]
    I_collected = abs(
        extract_steady_state_current(device_info, contact, node_set=node_set)
    )

    if I_generated <= 0:
        return 0.0

    cce = I_collected / I_generated
    return float(np.clip(cce, 0.0, 1.0))


def alpha_stripe_generation(
    device_info,
    x_pos_cm,
    gen_rate=1e18,
    stripe_sigma_cm=2e-4,
    depth_from=None,
    zero_on_contacts=True,
):
    """Alpha-like generation stripe entering the epi at lateral position ``x_pos_cm``.

    G(x, y) = gen_rate * depth_shape(y - y_entry) * exp(-0.5 ((x - x_pos) / sigma)^2),
    the profile :func:`cce_lateral_scan` uses at every lateral position: the 1D alpha
    depth shape (``alpha_depth_shape``, flat box + erfc roll-off, peak renormalised to
    ``gen_rate``) below the junction, zero in the p+ layer, and zero beyond the epi end
    when the device has an n+ layer (``confine_generation_to_epi``). ``depth_from``
    overrides the entry depth (default: the junction, ``substrate_thickness_cm``). For
    an axisymmetric device x is the radius and the stripe is a cylindrical shell.

    ``zero_on_contacts=True`` (default since 2026-10-18) zeroes G on the contact nodes.
    The contact equation replaces the continuity equation there (Dirichlet carrier
    densities), so pairs assigned to those nodes are never injected into the device;
    the returned array is then the generation actually injected. :func:`compute_cce_2d`
    normalises by the injected generation either way (since 2026-10-18), so the flag
    no longer moves a CCE; before, ``False`` counted the contact nodes in
    ``I_generated`` and the CCE carried a mesh-dependent deficit equal to their share
    (legacy planar diode at 50 V: 1.26e-3 with the 50 nm device2d mesh at the cathode,
    4.6e-3 with the coarser alternative-structure mesh; 0 on the n+ stack, where the
    stripe is confined to the epi). ``False`` keeps the raw stripe on every node.

    AUDIT Mj-5 (v5): the scan previously used a monotonically DECAYING exponential
    exp(-y/alpha), backwards for an alpha (energy loss per unit path RISES toward the
    end of range) and inconsistent with the 1D alpha profile; both now share
    ``alpha_depth_shape`` (qualitative, not Bragg-peaked by default).
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    x_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    y_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="y")
    )
    y_entry = (
        device_info["substrate_thickness_cm"] if depth_from is None else depth_from
    )
    y_from_entry = y_nodes - y_entry
    # alpha_depth_shape is a pure function of depth (no sorting requirement), safe on
    # unordered 2D node coordinates; clip the p+ side to 0 so the erfc roll-off only
    # sees non-negative depth, then renormalise to a peak of 1 (the contract here is a
    # peak rate, not a total per-alpha dose).
    shape = alpha_depth_shape(
        np.clip(y_from_entry, 0.0, None), alpha_range_cm=_ALPHA_RANGE_CM
    )
    peak = np.max(shape)
    if peak > 0:
        shape = shape / peak
    generation = (
        gen_rate * shape * np.exp(-0.5 * ((x_nodes - x_pos_cm) / stripe_sigma_cm) ** 2)
    )
    generation[y_from_entry < 0] = 0.0
    generation = confine_generation_to_epi(device_info, y_nodes, generation)
    if zero_on_contacts:
        generation = np.array(generation, dtype=float, copy=True)
        for c in contacts_on_region(device, region):
            generation[sorted(_contact_node_indices(device, region, c))] = 0.0
    return generation


def cce_with_generation(device_info, generation, contacts=("cathode",), node_sets=None):
    """Inject ``generation``, re-solve, and return the CCE at each of ``contacts``.

    The device is left as it was: the generation is reset to zero and the pre-injection
    solution (Potential, Electrons, Holes) restored, then one DC solve from that exact
    state (1-2 Newton iterations) so the contact currents devsim reports are consistent
    for the next user. ``node_sets`` maps a contact name to the ``node_set`` of
    :func:`compute_cce_2d` (guard ring: anode and guard share the majority carrier).

    Returns ``{contact: cce}``. A bare ``compute_cce_2d`` call after changing the
    generation array WITHOUT this inject-and-solve step reads the contact currents of
    whatever was solved last (the dark state): that was the notebook-19 bug.
    """
    node_sets = node_sets or {}
    return _solved_with_generation(
        device_info,
        generation,
        lambda G: {
            c: compute_cce_2d(device_info, G, c, node_set=node_sets.get(c))
            for c in contacts
        },
    )


def contact_current_sharing(device_info, generation, node_sets=None):
    """Signed balanced current of EVERY contact of the region under ``generation``,
    as a fraction of the generated current, plus the conservation residual.

    Returns ``(fractions, residual)``: ``fractions[contact] = I_contact / I_gen`` with
    devsim's sign (positive into the device; the cathode of a reverse-biased diode is
    positive, the p+ contacts negative), and ``residual = sum(I_contact) / I_gen``,
    which is zero at a converged steady state up to the solver tolerance: it is the
    a-posteriori check of a bias point accepted on a relaxed/aggressive tier. For a
    guard-ring device ``fractions`` tells how an edge event splits between the
    anode and the guard. ``node_sets`` as in :func:`cce_with_generation`. The device
    is restored afterwards.
    """
    node_sets = node_sets or {}
    device = device_info["device_name"]
    region = device_info["region_name"]
    contacts = contacts_on_region(device, region)
    # injected generation only (contact nodes are Dirichlet), as compute_cce_2d
    I_gen = Q * integrate_nodes_box_2d(
        device_info, injected_generation(device_info, generation)
    )
    currents = _solved_with_generation(
        device_info,
        generation,
        lambda G: {
            c: extract_steady_state_current(device_info, c, node_set=node_sets.get(c))
            for c in contacts
        },
    )
    if I_gen <= 0:
        raise ValueError("generation integrates to zero")
    fractions = {c: cur / I_gen for c, cur in currents.items()}
    return fractions, sum(currents.values()) / I_gen


def _solved_with_generation(device_info, generation, extract):
    """Inject ``generation``, solve, return ``extract(generation)``; then zero the
    generation, restore the pre-injection solution and solve once from it."""
    device = device_info["device_name"]
    region = device_info["region_name"]
    baseline = _snapshot_solution(device, region)
    generation = np.asarray(generation, dtype=float)
    try:
        add_generation_to_dd(device_info, generation)
        _robust_dc_solve(device_info)
        return extract(generation)
    finally:
        add_generation_to_dd(device_info, np.zeros_like(generation))
        _restore_solution(device, region, baseline)
        _robust_dc_solve(device_info)


def cce_lateral_scan(
    device_info,
    n_points=20,
    gen_rate=1e18,
    stripe_sigma_cm=2e-4,
    contact="cathode",
    node_set=None,
    x_max_cm=None,
    zero_on_contacts=True,
):
    """Compute CCE at multiple lateral positions from center to edge.

    At each lateral position x_pos the generation is the alpha stripe of
    :func:`alpha_stripe_generation` (flat-box + erfc depth shape below the junction,
    Gaussian lateral stripe), injected and solved, and the CCE read at ``contact``.

    Parameters
    ----------
    device_info : dict
        Device info dict with DD setup and bias applied.
    n_points : int
        Number of lateral positions to scan.
    gen_rate : float
        Peak generation rate (cm^-3 s^-1) -- the depth profile's plateau
        value, not a total dose.
    stripe_sigma_cm : float
        Gaussian width of the lateral generation stripe (cm).
    contact : str
        Contact for current extraction.
    node_set : array_like, optional
        ``node_set`` of :func:`extract_steady_state_current` for ``contact`` (the
        guard-ring anode needs one; two-terminal devices do not).
    x_max_cm : float, optional
        Last lateral position (default ``device_info["half_width_cm"]``).
    zero_on_contacts : bool
        Passed to :func:`alpha_stripe_generation` (default True since 2026-10-18). The
        CCE does not depend on it: :func:`compute_cce_2d` normalises by the injected
        (non-contact) generation.

    Returns
    -------
    result : dict
        Dictionary with keys:
        - "x_positions_cm": lateral positions scanned (cm)
        - "x_positions_um": lateral positions in micrometers
        - "cce_values": CCE at each position, relative to the pairs generated
          inside the simulated epi (I_collected / I_generated in the mesh)
        - "edge_to_center_ratio": CCE(edge) / CCE(center)
        - "in_mesh_generation_fraction": fraction of the alpha's depth profile
          inside the epi (in_mesh_generation_fraction)
        - "q_collected_over_deposited": cce_values * that fraction, i.e. the
          collected charge over the total pairs the alpha would deposit if its
          range were not cut by the epi/cathode
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    half_width_cm = device_info["half_width_cm"] if x_max_cm is None else x_max_cm
    alpha_range_cm = _ALPHA_RANGE_CM

    # Lateral positions: center (x=0) to edge (x=half_width)
    x_positions = np.linspace(0, half_width_cm, n_points)
    cce_values = []
    baseline = None

    for i, x_pos in enumerate(x_positions):
        generation = alpha_stripe_generation(
            device_info,
            x_pos,
            gen_rate=gen_rate,
            stripe_sigma_cm=stripe_sigma_cm,
            zero_on_contacts=zero_on_contacts,
        )

        # Inject generation, solve, extract CCE
        add_generation_to_dd(device_info, generation)
        _robust_dc_solve(device_info)
        cce = compute_cce_2d(device_info, generation, contact, node_set=node_set)
        cce_values.append(cce)

        logger.info(
            f"Lateral scan [{i + 1}/{n_points}]: x={x_pos * 1e4:.1f} um, CCE={cce:.4f}"
        )

        # Reset generation to zero and restore the baseline. The first reset is
        # solved (that is the baseline); later ones restore its snapshot instead of
        # re-solving (~9-11 Newton iterations each, ~75% of the scan time).
        add_generation_to_dd(device_info, np.zeros_like(generation))
        if baseline is None:
            _robust_dc_solve(device_info)
            baseline = _snapshot_solution(device, region)
        else:
            _restore_solution(device, region, baseline)

    # Contact currents are evaluated during assembly: one solve from the exact
    # baseline (1-2 iterations) leaves the device consistent for later users.
    _robust_dc_solve(device_info)

    cce_values = np.array(cce_values)
    edge_to_center = cce_values[-1] / cce_values[0] if cce_values[0] > 0 else 0.0
    # depth over which pairs are generated: the epi, plus the n+ layer only if the device
    # was built with nplus_generation=True (confine_generation_to_epi)
    gen_depth = device_info["epi_thickness_cm"]
    if device_info.get("nplus_thickness_cm") is not None and device_info.get(
        "nplus_generation"
    ):
        gen_depth += device_info["nplus_thickness_cm"]
    f_in = in_mesh_generation_fraction(gen_depth, alpha_range_cm=alpha_range_cm)

    return {
        "x_positions_cm": x_positions,
        "x_positions_um": x_positions * 1e4,
        "cce_values": cce_values,
        "edge_to_center_ratio": float(edge_to_center),
        "in_mesh_generation_fraction": f_in,
        "q_collected_over_deposited": cce_values * f_in,
    }


def cce_heatmap_2d(device_info, cce_lateral, cce_1d_depth=None):
    """Construct a 2D CCE map using factored form CCE(x,y) = f_edge(x) * f_depth(y).

    f_edge(x) is interpolated from the lateral scan, normalized by the
    center value.  f_depth(y) is either provided or defaults to 1 in the
    depletion region and 0 outside.

    Parameters
    ----------
    device_info : dict
        Device info dict from create_sic_2d_device.
    cce_lateral : dict
        Result from cce_lateral_scan().
    cce_1d_depth : array_like or None
        Per-node depth CCE profile.  If None, uses step function:
        1.0 in the epi (junction_pos <= y <= epi_end), 0.0 in the p+ layer and in
        the n+ layer (default stack), which is not active volume.

    Returns
    -------
    result : dict
        Dictionary with keys:
        - "x_nodes_cm": x-coordinates of mesh nodes (cm)
        - "y_nodes_cm": y-coordinates of mesh nodes (cm)
        - "cce_map": per-node CCE values
        - "active_fraction": fraction of epi nodes (n+ excluded) with CCE > 0.5
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    junction_pos = device_info["substrate_thickness_cm"]

    x_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    y_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="y")
    )

    # f_edge(x): lateral efficiency normalized by center value
    x_lat = cce_lateral["x_positions_cm"]
    cce_lat = cce_lateral["cce_values"]
    cce_center = cce_lat[0] if cce_lat[0] > 0 else 1.0
    f_edge_norm = cce_lat / cce_center

    # Interpolate f_edge to mesh node x-positions
    f_edge = np.interp(x_nodes, x_lat, f_edge_norm)

    # Epi = [junction_pos, epi_end]. With an n+ layer (default stack) the mesh goes on
    # to the cathode below epi_end; y >= junction_pos alone counted the n+ as active
    # epi (same bug class as compare_cce_2d_vs_1d, commit dbc01b0).
    epi_end = device_info.get("epi_end", junction_pos + device_info["epi_thickness_cm"])
    epi_mask = (y_nodes >= junction_pos) & (y_nodes <= epi_end + 1e-12)

    # f_depth(y): depth collection efficiency
    if cce_1d_depth is not None:
        f_depth = np.asarray(cce_1d_depth, dtype=float)
    else:
        # Default: 1.0 in the epi, 0.0 in the p+ and n+ layers
        f_depth = np.where(epi_mask, 1.0, 0.0)

    # Factored 2D CCE map
    cce_map = f_edge * f_depth

    # Active fraction: fraction of epi nodes with CCE > 0.5
    if np.sum(epi_mask) > 0:
        active_fraction = float(np.sum(cce_map[epi_mask] > 0.5) / np.sum(epi_mask))
    else:
        active_fraction = 0.0

    return {
        "x_nodes_cm": x_nodes,
        "y_nodes_cm": y_nodes,
        "cce_map": cce_map,
        "active_fraction": active_fraction,
    }


def _engine_scan_or_none(half_width_um, V_bias, n_points, gen_rate, full_area=False):
    """``lateral_cce_scan_2d`` for ``solver="auto"``; None when the 2D engine declines
    (outside its validated box, not converged), so the caller runs devsim."""
    from etna.core.engine2d import (
        Engine2DFallback,
        Engine2DOutsideValidity,
        engine2d_unsupported_reason,
        lateral_cce_scan_2d,
    )

    reason = engine2d_unsupported_reason(
        half_width_um, V_bias, n_lateral_points=n_points
    )
    if reason is not None:
        logger.info(f"engine2d unsupported ({reason}); using devsim")
        return None
    try:
        return lateral_cce_scan_2d(
            half_width_um,
            V_bias,
            n_lateral_points=n_points,
            gen_rate=gen_rate,
            full_area=full_area,
        )
    except (Engine2DFallback, Engine2DOutsideValidity) as e:
        logger.info(f"engine2d declined ({e}); using devsim")
        return None


def _devsim_scan_keys(scan):
    """The keys ``cce_lateral_scan`` returns, from an engine result."""
    keys = (
        "x_positions_cm",
        "x_positions_um",
        "cce_values",
        "edge_to_center_ratio",
        "in_mesh_generation_fraction",
        "q_collected_over_deposited",
    )
    return {k: scan[k] for k in keys}


def _compare_2d_devsim(half_width_um, V_bias, gen_rate):
    """devsim 2D half of compare_cce_2d_vs_1d: (full-area CCE, lateral scan, half width cm)."""
    device_info_2d = None
    try:
        device_info_2d = create_2d_dd_device(half_width_um=half_width_um, V_bias=V_bias)
        device_2d = device_info_2d["device_name"]
        region_2d = device_info_2d["region_name"]
        junction_pos = device_info_2d["substrate_thickness_cm"]

        # Get 2D mesh coordinates
        y_nodes_2d = np.array(
            devsim.get_node_model_values(device=device_2d, region=region_2d, name="y")
        )

        # Uniform epi generation: G = gen_rate for y > junction, 0 otherwise
        gen_2d = np.where(y_nodes_2d >= junction_pos, gen_rate, 0.0)
        gen_2d = confine_generation_to_epi(device_info_2d, y_nodes_2d, gen_2d)

        # Compute 2D CCE (full area, uniform generation)
        add_generation_to_dd(device_info_2d, gen_2d)
        _robust_dc_solve()
        cce_2d_full = compute_cce_2d(device_info_2d, gen_2d)

        # Reset generation before lateral scan
        x_nodes_2d = np.array(
            devsim.get_node_model_values(device=device_2d, region=region_2d, name="x")
        )
        zero_gen = np.zeros_like(x_nodes_2d)
        add_generation_to_dd(device_info_2d, zero_gen)
        _robust_dc_solve()

        # Lateral scan for edge effects
        lateral_scan = cce_lateral_scan(device_info_2d, n_points=10, gen_rate=gen_rate)

        # Delete the 2D device before compare_cce_2d_vs_1d creates the 1D device:
        # devsim.solve() is global across all loaded devices, so the 2D device
        # (with residual DD state from lateral scan) causes convergence failure
        # when the 1D device tries its initial equilibrium solve.
        half_width_cm = device_info_2d["half_width_cm"]
        delete_device_and_mesh(device_info_2d["device_name"])
        device_info_2d = None  # Prevent double-delete in finally
        return cce_2d_full, lateral_scan, half_width_cm

    finally:
        if device_info_2d is not None:
            try:
                delete_device_and_mesh(device_info_2d["device_name"])
            except Exception:
                pass


def compare_cce_2d_vs_1d(
    half_width_um=50.0, V_bias=50.0, gen_rate=1e18, solver="devsim"
):
    """Compare 2D and 1D CCE for the same device parameters.

    Creates both a 2D device and a 1D device, applies uniform epi
    generation, and computes CCE for each.  Also runs a lateral scan
    on the 2D device to quantify edge effects and compute the
    active-to-geometric volume ratio.

    Parameters
    ----------
    half_width_um : float
        Half-width of 2D SV in micrometers.
    V_bias : float
        Reverse bias voltage (V, positive).
    gen_rate : float
        Uniform generation rate in epi (cm^-3 s^-1).
    solver : {"devsim", "auto"}
        ``"auto"``: the 2D half (full-area CCE and lateral scan) runs on the 2D engine
        (``lateral_cce_scan_2d(..., full_area=True)``, parity with devsim ~1e-14) when
        the device is inside its box, else on devsim. The 1D half is always devsim.

    Returns
    -------
    result : dict
        Dictionary with keys:
        - "solver_2d": "engine2d" or "devsim"
        - "cce_1d": 1D CCE value
        - "cce_2d_center": 2D CCE of the centre column (lateral scan, x = 0)
        - "cce_2d_full_area": 2D CCE with uniform epi generation over the full area
          (before 2026-10-18 "cce_2d_center" returned this value too: audit B-14)
        - "active_to_geometric_ratio": effective active / geometric area
        - "lateral_scan": result from cce_lateral_scan
        - "half_width_um": half-width used
        - "V_bias": bias voltage used
    """
    from etna.core.drift_diffusion import create_dd_device
    from etna.core.charge_collection import compute_cce_from_dd

    if solver not in ("devsim", "auto"):
        raise ValueError(f"solver must be 'devsim' or 'auto', got {solver!r}")
    device_info_1d = None
    engine = (
        _engine_scan_or_none(half_width_um, V_bias, 10, gen_rate, full_area=True)
        if solver == "auto"
        else None
    )

    if engine is not None:
        cce_2d_full = engine["cce_full_area"]
        lateral_scan = _devsim_scan_keys(engine)
        half_width_cm = engine["device_info"]["half_width_cm"]
    else:
        # devsim.solve() is global across all loaded devices: the 2D device is deleted
        # (inside the helper) before the 1D device is created.
        cce_2d_full, lateral_scan, half_width_cm = _compare_2d_devsim(
            half_width_um, V_bias, gen_rate
        )

    try:
        # --- 1D device ---
        dev_id_1d = uuid.uuid4().hex[:8]
        device_info_1d = create_dd_device(
            device_name=f"cce1d_{dev_id_1d}",
            doping_profile="graded",
        )

        # Ramp 1D bias (anode, negative for reverse bias)
        ramp_bias(device_info_1d, V_target=-V_bias, contact="anode", V_step=0.5)

        # 1D generation: uniform in epi
        device_1d = device_info_1d["device_name"]
        region_1d = device_info_1d["region_name"]
        x_nodes_1d = np.array(
            devsim.get_node_model_values(device=device_1d, region=region_1d, name="x")
        )
        junction_1d = device_info_1d["junction_pos"]
        gen_1d = np.where(x_nodes_1d >= junction_1d, gen_rate, 0.0)
        # same confinement as the 2D side (the default stack has an n+ layer)
        gen_1d = confine_generation_to_epi(device_info_1d, x_nodes_1d, gen_1d)

        add_generation_to_dd(device_info_1d, gen_1d)
        _robust_dc_solve()
        cce_1d = compute_cce_from_dd(device_info_1d, gen_1d)

        # Active-to-geometric ratio from lateral scan
        # = integral(CCE(x) dx) / (half_width * CCE_center)
        x_lat = lateral_scan["x_positions_cm"]
        cce_lat = lateral_scan["cce_values"]
        # half_width_cm was saved before 2D device deletion
        cce_center = cce_lat[0]

        if cce_center > 0 and half_width_cm > 0:
            integral_cce = np.trapezoid(cce_lat, x_lat)
            active_to_geometric = integral_cce / (half_width_cm * cce_center)
        else:
            active_to_geometric = 0.0

        return {
            "solver_2d": "devsim" if engine is None else "engine2d",
            "cce_1d": float(cce_1d),
            "cce_2d_center": float(cce_center),
            "cce_2d_full_area": float(cce_2d_full),
            "active_to_geometric_ratio": float(active_to_geometric),
            "lateral_scan": lateral_scan,
            "half_width_um": half_width_um,
            "V_bias": V_bias,
        }

    finally:
        # Clean up devices
        if device_info_1d is not None:
            try:
                delete_device_and_mesh(device_info_1d["device_name"])
            except Exception:
                pass


def cce_vs_bias_lateral(
    half_width_um=50.0, biases=None, n_points=5, gen_rate=1e18, solver="devsim"
):
    """Compute lateral CCE profiles at multiple bias voltages.

    Shows how edge effects evolve from partial to full depletion.
    At low bias (partial depletion), carriers in the undepleted region
    must diffuse to the depletion edge for collection, creating
    position-dependent CCE.  At full depletion, drift dominates
    and CCE is uniform.

    Parameters
    ----------
    half_width_um : float
        Half-width of 2D SV in micrometers.
    biases : list of float or None
        Bias voltages to scan (V, positive). Default: [5, 10, 20, 30, 50].
    n_points : int
        Number of lateral positions per scan.
    gen_rate : float
        Generation rate (cm^-3 s^-1).
    solver : {"devsim", "auto"}
        ``"auto"``: each bias runs on the 2D engine (``lateral_cce_scan_2d``, parity
        with devsim ~1e-14) when inside its box, else on devsim.

    Returns
    -------
    result : dict
        Dictionary with keys:
        - "biases": list of bias voltages
        - "solvers": per-bias "engine2d" / "devsim"
        - "scans": list of cce_lateral_scan results (one per bias)
        - "half_width_um": half-width used
    """
    if biases is None:
        biases = [5.0, 10.0, 20.0, 30.0, 50.0]

    if solver not in ("devsim", "auto"):
        raise ValueError(f"solver must be 'devsim' or 'auto', got {solver!r}")
    scans = []
    solvers = []
    for V in biases:
        logger.info(f"Bias scan: V={V:.0f}V, half_width={half_width_um:.0f}um")
        if solver == "auto":
            engine = _engine_scan_or_none(half_width_um, V, n_points, gen_rate)
            if engine is not None:
                scans.append(_devsim_scan_keys(engine))
                solvers.append("engine2d")
                continue
        solvers.append("devsim")
        device_info = None
        try:
            device_info = create_2d_dd_device(half_width_um=half_width_um, V_bias=V)
            scan = cce_lateral_scan(device_info, n_points=n_points, gen_rate=gen_rate)
            scans.append(scan)
        finally:
            if device_info is not None:
                try:
                    delete_device_and_mesh(device_info["device_name"])
                except Exception:
                    pass

    return {
        "biases": biases,
        "scans": scans,
        "solvers": solvers,
        "half_width_um": half_width_um,
    }
