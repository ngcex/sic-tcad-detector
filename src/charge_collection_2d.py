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
    - Petringa et al. experimental microdosimeter geometry
"""

import logging
import uuid

import devsim
import numpy as np

from src.device2d import create_sic_2d_device
from src.poisson import setup_poisson, solve_equilibrium
from src.drift_diffusion import (
    setup_sic_drift_diffusion,
    extract_contact_current,
    ramp_bias,
)
from src.charge_collection import add_generation_to_dd

logger = logging.getLogger(__name__)

# Physical constants (CGS)
Q = 1.602e-19  # C, elementary charge

# Default alpha particle range in SiC
_ALPHA_RANGE_CM = 15e-4  # 15 um


def _robust_dc_solve():
    """Solve DC with fallback to relaxed tolerances on convergence failure."""
    try:
        devsim.solve(
            type="dc",
            absolute_error=1e10,
            relative_error=1e-10,
            maximum_iterations=40,
        )
    except devsim.error:
        devsim.solve(
            type="dc",
            absolute_error=1e12,
            relative_error=1e-8,
            maximum_iterations=100,
        )
        logger.info("DC solve converged with relaxed tolerances")


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
    ramp_bias(device_info, V_target=V_bias, contact="cathode", V_step=0.5)
    logger.info(
        f"Created 2D DD device '{device_name}': "
        f"half_width={half_width_um:.0f} um, V_bias={V_bias:.1f} V"
    )

    return device_info


def compute_cce_2d(device_info, generation_values, contact="cathode"):
    """Compute 2D charge collection efficiency from DD simulation.

    CCE = |I_collected| / I_generated, where I_generated is computed via
    area integration of the generation rate over the 2D mesh.

    The symmetry factor of 2 (half-device) cancels in the ratio since
    both I_collected and I_generated are for the half-device.

    Parameters
    ----------
    device_info : dict
        Device info dict with DD solved and generation added.
    generation_values : array_like
        Generation rate G(x,y) at each mesh node (cm^-3 s^-1).
    contact : str
        Contact at which to extract collected current.

    Returns
    -------
    cce : float
        Charge collection efficiency, clipped to [0, 1].
    """
    # I_generated = Q * integral(G dA)  [A/cm for 2D]
    I_generated = Q * integrate_over_mesh_2d(device_info, generation_values)

    # I_collected at contact [A/cm for 2D]
    I_collected = abs(extract_contact_current(device_info, contact))

    if I_generated <= 0:
        return 0.0

    cce = I_collected / I_generated
    return float(np.clip(cce, 0.0, 1.0))


def cce_lateral_scan(
    device_info,
    n_points=20,
    gen_rate=1e18,
    stripe_sigma_cm=2e-4,
    contact="cathode",
):
    """Compute CCE at multiple lateral positions from center to edge.

    At each lateral position x_pos, creates a Gaussian-stripe generation
    profile G(x,y) = gen_rate * exp(-y/alpha) * exp(-0.5*((x-x_pos)/sigma)^2)
    where y is measured from the junction (epi entrance).  Generation is
    zero in the substrate (y < junction_pos).

    Parameters
    ----------
    device_info : dict
        Device info dict with DD setup and bias applied.
    n_points : int
        Number of lateral positions to scan.
    gen_rate : float
        Peak generation rate (cm^-3 s^-1).
    stripe_sigma_cm : float
        Gaussian width of the lateral generation stripe (cm).
    contact : str
        Contact for current extraction.

    Returns
    -------
    result : dict
        Dictionary with keys:
        - "x_positions_cm": lateral positions scanned (cm)
        - "x_positions_um": lateral positions in micrometers
        - "cce_values": CCE at each position
        - "edge_to_center_ratio": CCE(edge) / CCE(center)
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    half_width_cm = device_info["half_width_cm"]
    substrate_thickness_cm = device_info["substrate_thickness_cm"]
    alpha = _ALPHA_RANGE_CM

    # Get mesh node coordinates
    x_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    y_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="y")
    )

    # Depth measured from junction (epi entrance)
    junction_pos = substrate_thickness_cm
    y_from_junction = y_nodes - junction_pos

    # Lateral positions: center (x=0) to edge (x=half_width)
    x_positions = np.linspace(0, half_width_cm, n_points)
    cce_values = []

    for i, x_pos in enumerate(x_positions):
        # Build generation profile: exponential depth * Gaussian lateral stripe
        depth_profile = gen_rate * np.exp(-y_from_junction / alpha)
        lateral_profile = np.exp(-0.5 * ((x_nodes - x_pos) / stripe_sigma_cm) ** 2)
        generation = depth_profile * lateral_profile

        # Zero out generation in substrate (y < junction_pos)
        generation[y_from_junction < 0] = 0.0

        # Inject generation, solve, extract CCE
        add_generation_to_dd(device_info, generation)
        _robust_dc_solve()
        cce = compute_cce_2d(device_info, generation, contact)
        cce_values.append(cce)

        logger.info(
            f"Lateral scan [{i+1}/{n_points}]: " f"x={x_pos*1e4:.1f} um, CCE={cce:.4f}"
        )

        # Reset generation to zero and re-solve to restore baseline
        zero_gen = np.zeros_like(x_nodes)
        add_generation_to_dd(device_info, zero_gen)
        _robust_dc_solve()

    cce_values = np.array(cce_values)
    edge_to_center = cce_values[-1] / cce_values[0] if cce_values[0] > 0 else 0.0

    return {
        "x_positions_cm": x_positions,
        "x_positions_um": x_positions * 1e4,
        "cce_values": cce_values,
        "edge_to_center_ratio": float(edge_to_center),
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
        1.0 in the epi (y > junction_pos), 0.0 in substrate.

    Returns
    -------
    result : dict
        Dictionary with keys:
        - "x_nodes_cm": x-coordinates of mesh nodes (cm)
        - "y_nodes_cm": y-coordinates of mesh nodes (cm)
        - "cce_map": per-node CCE values
        - "active_fraction": fraction of epi area with CCE > 0.5
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

    # f_depth(y): depth collection efficiency
    if cce_1d_depth is not None:
        f_depth = np.asarray(cce_1d_depth, dtype=float)
    else:
        # Default: 1.0 in epi (y >= junction_pos), 0.0 in substrate
        f_depth = np.where(y_nodes >= junction_pos, 1.0, 0.0)

    # Factored 2D CCE map
    cce_map = f_edge * f_depth

    # Active fraction: fraction of epi area with CCE > 0.5
    epi_mask = y_nodes >= junction_pos
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


def compare_cce_2d_vs_1d(half_width_um=50.0, V_bias=50.0, gen_rate=1e18):
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

    Returns
    -------
    result : dict
        Dictionary with keys:
        - "cce_1d": 1D CCE value
        - "cce_2d_center": 2D CCE with uniform epi generation
        - "cce_2d_full_area": 2D CCE integrated over full area
        - "active_to_geometric_ratio": effective active / geometric area
        - "lateral_scan": result from cce_lateral_scan
        - "half_width_um": half-width used
        - "V_bias": bias voltage used
    """
    from src.drift_diffusion import create_dd_device
    from src.charge_collection import compute_cce_from_dd

    device_info_2d = None
    device_info_1d = None

    try:
        # --- 2D device ---
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

        # Save 2D results and delete 2D device before creating 1D device.
        # devsim.solve() is global across all loaded devices, so the 2D device
        # (with residual DD state from lateral scan) causes convergence failure
        # when the 1D device tries its initial equilibrium solve.
        half_width_cm = device_info_2d["half_width_cm"]
        devsim.delete_device(device=device_info_2d["device_name"])
        device_info_2d = None  # Prevent double-delete in finally

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
            "cce_1d": float(cce_1d),
            "cce_2d_center": float(cce_2d_full),
            "cce_2d_full_area": float(cce_2d_full),
            "active_to_geometric_ratio": float(active_to_geometric),
            "lateral_scan": lateral_scan,
            "half_width_um": half_width_um,
            "V_bias": V_bias,
        }

    finally:
        # Clean up devices
        if device_info_2d is not None:
            try:
                devsim.delete_device(device=device_info_2d["device_name"])
            except Exception:
                pass
        if device_info_1d is not None:
            try:
                devsim.delete_device(device=device_info_1d["device_name"])
            except Exception:
                pass


def cce_vs_bias_lateral(half_width_um=50.0, biases=None, n_points=5, gen_rate=1e18):
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

    Returns
    -------
    result : dict
        Dictionary with keys:
        - "biases": list of bias voltages
        - "scans": list of cce_lateral_scan results (one per bias)
        - "half_width_um": half-width used
    """
    if biases is None:
        biases = [5.0, 10.0, 20.0, 30.0, 50.0]

    scans = []
    for V in biases:
        logger.info(f"Bias scan: V={V:.0f}V, half_width={half_width_um:.0f}um")
        device_info = None
        try:
            device_info = create_2d_dd_device(half_width_um=half_width_um, V_bias=V)
            scan = cce_lateral_scan(device_info, n_points=n_points, gen_rate=gen_rate)
            scans.append(scan)
        finally:
            if device_info is not None:
                try:
                    devsim.delete_device(device=device_info["device_name"])
                except Exception:
                    pass

    return {
        "biases": biases,
        "scans": scans,
        "half_width_um": half_width_um,
    }
