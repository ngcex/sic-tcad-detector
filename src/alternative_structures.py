"""Alternative SiC microdosimeter structure mesh builders.

Provides geometry-specific mesh builders for non-planar microdosimeter
designs: mesa-etched, 3D electrode (axisymmetric), stacked delta-E/E
telescope, and guard ring structures.

Each builder returns a device_info dict compatible with the existing
pipeline (Poisson, DD, transient, CCE, microdosimetry), with an
additional ``structure_type`` key to identify the geometry.

Coordinate convention (same as device2d.py):
    x = lateral direction (or radial for 3D electrode)
    y = depth direction (0 = top surface, total_depth = bottom)

All units CGS (cm, cm^-3, F/cm, eV, s) per devsim convention.
"""

import logging

import devsim
import devsim.python_packages.simple_physics as simple_physics

from src.sic_material import (
    SiC4H_Parameters,
    intrinsic_concentration,
    mobility_caughey_thomas_T,
    srh_lifetime,
)
from src.incomplete_ionization import ionized_acceptor_concentration
from src.device2d import (
    set_doping_profile_2d,
    set_graded_doping_2d,
    EPS_0,
    Q,
    K_B_EV,
    K_B_J,
    DEFAULT_N_D,
    _N_D_JUNCTION_DEFAULT,
    _N_D_BULK_DEFAULT,
    _L_TRANSITION_DEFAULT,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _set_sic_material_params(device_name, region_name, N_D, N_A, T):
    """Replicate 4H-SiC material parameter setup from device2d.py.

    Sets Permittivity, ElectronCharge, n_i, T, kT, V_t, mu_n, mu_p,
    and SRH parameters on the given device/region.

    Returns a partial dict with material-related keys.
    """
    params = SiC4H_Parameters()
    kT_eV = K_B_EV * T
    kT_J = K_B_J * T
    V_t = kT_J / Q

    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="Permittivity",
        value=params.eps_r * EPS_0,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="ElectronCharge",
        value=Q,
    )

    n_i_T, NC_T, NV_T, E_g_T = intrinsic_concentration(T, params)
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="n_i",
        value=n_i_T,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="T",
        value=T,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="kT",
        value=kT_J,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="V_t",
        value=V_t,
    )

    mu_n = mobility_caughey_thomas_T(N_D, T, "electron", params)
    mu_p = mobility_caughey_thomas_T(N_A, T, "hole", params)
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="mu_n",
        value=mu_n,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="mu_p",
        value=mu_p,
    )

    # SRH parameters
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="n1",
        value=n_i_T,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="p1",
        value=n_i_T,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="taun",
        value=srh_lifetime(T, "electron", params),
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="taup",
        value=srh_lifetime(T, "hole", params),
    )

    return {
        "params": params,
        "n_i": n_i_T,
        "E_g": E_g_T,
        "mu_n": mu_n,
        "mu_p": mu_p,
        "T": T,
    }


def _resolve_doping_defaults(N_D, doping_profile, N_D_junction, N_D_bulk, L_transition):
    """Resolve doping parameter defaults matching device2d.py logic."""
    if N_D is None:
        N_D = DEFAULT_N_D
    _N_D_junction = N_D_junction if N_D_junction is not None else _N_D_JUNCTION_DEFAULT
    _N_D_bulk = N_D_bulk if N_D_bulk is not None else _N_D_BULK_DEFAULT
    _L_transition = L_transition if L_transition is not None else _L_TRANSITION_DEFAULT
    return N_D, _N_D_junction, _N_D_bulk, _L_transition


def _apply_doping(
    device_name,
    region_name,
    junction_pos,
    N_A_ionized,
    N_D,
    doping_profile,
    N_D_junction,
    N_D_bulk,
    L_transition,
):
    """Apply doping profile (graded or uniform) on a region."""
    N_D, _ndj, _ndb, _lt = _resolve_doping_defaults(
        N_D,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )
    if doping_profile == "graded":
        set_graded_doping_2d(
            device_name,
            region_name,
            junction_pos,
            N_A_ionized,
            _ndj,
            _ndb,
            _lt,
        )
        logger.info(
            f"Graded 2D doping on '{region_name}': N_D_junction={_ndj:.2e}, "
            f"N_D_bulk={_ndb:.2e}, L_transition={_lt:.2e} cm"
        )
    else:
        set_doping_profile_2d(
            device_name,
            region_name,
            junction_pos,
            N_A_ionized,
            N_D,
        )
    return N_D


def _build_device_info(
    device_name,
    region_name,
    junction_pos,
    epi_thickness_cm,
    substrate_thickness_cm,
    N_A,
    N_A_ionized,
    N_D,
    doping_profile,
    N_D_junction,
    N_D_bulk,
    L_transition,
    half_width_cm,
    material_info,
    extra_keys=None,
):
    """Assemble a complete device_info dict."""
    num_nodes = len(
        devsim.get_node_model_values(
            device=device_name,
            region=region_name,
            name="x",
        )
    )
    info = {
        "device_name": device_name,
        "region_name": region_name,
        "junction_pos": junction_pos,
        "epi_thickness_cm": epi_thickness_cm,
        "substrate_thickness_cm": substrate_thickness_cm,
        "total_length": substrate_thickness_cm + epi_thickness_cm,
        "N_D": N_D,
        "N_A": N_A,
        "N_A_ionized": N_A_ionized,
        "T": material_info["T"],
        "n_i": material_info["n_i"],
        "E_g": material_info["E_g"],
        "params": material_info["params"],
        "mu_n": material_info["mu_n"],
        "mu_p": material_info["mu_p"],
        "num_nodes": num_nodes,
        "doping_profile": doping_profile,
        "N_D_junction": N_D_junction,
        "N_D_bulk": N_D_bulk,
        "L_transition": L_transition,
        "half_width_cm": half_width_cm,
        "dimension": 2,
    }
    if extra_keys:
        info.update(extra_keys)
    return info


def _add_depth_mesh_lines(mesh_name, junction_pos, total_depth, air_buffer):
    """Add standard depth (y) mesh lines matching device2d.py pattern."""
    # Top air buffer
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=-air_buffer, ps=air_buffer)

    # Anode surface (y=0)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=0, ps=1e-5)

    # Approach junction from p-side
    if junction_pos - 5e-6 > 1e-8:
        devsim.add_2d_mesh_line(
            mesh=mesh_name,
            dir="y",
            pos=junction_pos - 5e-6,
            ps=1e-6,
        )

    # Junction: very fine spacing
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=junction_pos, ps=1e-7)

    # Epi intermediate points
    epi_mid1 = junction_pos + 2e-4
    epi_mid2 = junction_pos + 5e-4
    if epi_mid1 < total_depth - 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=epi_mid1, ps=5e-6)
    if epi_mid2 < total_depth - 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=epi_mid2, ps=5e-6)

    # Cathode surface
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=total_depth, ps=1e-5)

    # Bottom air buffer
    devsim.add_2d_mesh_line(
        mesh=mesh_name,
        dir="y",
        pos=total_depth + air_buffer,
        ps=air_buffer,
    )


# ---------------------------------------------------------------------------
# Public: Mesa-etched structure
# ---------------------------------------------------------------------------


def create_mesa_device(
    device_name="mesa2d",
    region_name="sic",
    pillar_half_width_um=50.0,
    trench_width_um=10.0,
    trench_depth_um=10.0,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=1e-4,
    N_A=1e19,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
):
    """Create a mesa-etched SiC microdosimeter (half-device with symmetry).

    Structure: SiC pillar (active SV) with trenches on both sides.
    Uses mirror symmetry at x=0.

    Lateral layout (x):
        [0, pillar_hw]                     = SiC pillar (active)
        [pillar_hw, pillar_hw + trench_w]  = trench fill (zero doping)

    Trench extends from y=0 down to min(trench_depth, total_depth).
    If trench is shallower than total depth, substrate below trench is
    a separate region with bulk substrate doping.

    Contacts: anode at y=0 on pillar only, cathode at y=total_depth full width.

    Returns device_info with ``structure_type="mesa"``.
    """
    N_D_resolved, _, _, _ = _resolve_doping_defaults(
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )
    params = SiC4H_Parameters()
    pillar_hw_cm = pillar_half_width_um * 1e-4
    trench_w_cm = trench_width_um * 1e-4
    trench_d_cm = trench_depth_um * 1e-4
    junction_pos = substrate_thickness_cm
    total_depth = substrate_thickness_cm + epi_thickness_cm
    air_buffer = 1e-8
    total_width = pillar_hw_cm + trench_w_cm

    # Clamp trench depth to total device depth
    trench_d_cm = min(trench_d_cm, total_depth)

    mesh_name = f"{device_name}_mesh"
    devsim.create_2d_mesh(mesh=mesh_name)

    # -- Lateral (x) mesh lines --
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0, ps=5e-5)
    # Fine spacing at pillar/trench boundary
    if pillar_hw_cm > 2e-4:
        devsim.add_2d_mesh_line(
            mesh=mesh_name,
            dir="x",
            pos=pillar_hw_cm * 0.5,
            ps=2e-4,
        )
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=pillar_hw_cm, ps=1e-5)
    devsim.add_2d_mesh_line(
        mesh=mesh_name,
        dir="x",
        pos=pillar_hw_cm + trench_w_cm * 0.5,
        ps=5e-5,
    )
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=total_width, ps=5e-5)

    # -- Depth (y) mesh lines --
    _add_depth_mesh_lines(mesh_name, junction_pos, total_depth, air_buffer)

    # Add mesh line at trench bottom if it doesn't coincide with existing lines
    if trench_d_cm > 1e-8 and trench_d_cm < total_depth - 1e-8:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=trench_d_cm, ps=1e-6)

    # -- Regions --
    # Top air buffer (full width)
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_top",
        xl=0,
        xh=total_width,
        yl=-air_buffer,
        yh=0,
    )

    # Main SiC pillar region
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region=region_name,
        xl=0,
        xh=pillar_hw_cm,
        yl=0,
        yh=total_depth,
    )

    # Trench fill region (SiC material, will set zero doping)
    if trench_d_cm >= total_depth - 1e-10:
        # Trench goes full depth
        devsim.add_2d_region(
            mesh=mesh_name,
            material="SiC",
            region="trench",
            xl=pillar_hw_cm,
            xh=total_width,
            yl=0,
            yh=total_depth,
        )
    else:
        # Partial trench
        devsim.add_2d_region(
            mesh=mesh_name,
            material="SiC",
            region="trench",
            xl=pillar_hw_cm,
            xh=total_width,
            yl=0,
            yh=trench_d_cm,
        )
        # Substrate below trench
        devsim.add_2d_region(
            mesh=mesh_name,
            material="SiC",
            region="sub_trench",
            xl=pillar_hw_cm,
            xh=total_width,
            yl=trench_d_cm,
            yh=total_depth,
        )

    # Bottom air buffer (full width)
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_bot",
        xl=0,
        xh=total_width,
        yl=total_depth,
        yh=total_depth + air_buffer,
    )

    # -- Contacts --
    # Anode at y=0 on pillar only
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="anode",
        material="metal",
        region=region_name,
        yl=0,
        yh=0,
        bloat=1e-10,
    )
    # Cathode at y=total_depth (on pillar region)
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="cathode",
        material="metal",
        region=region_name,
        yl=total_depth,
        yh=total_depth,
        bloat=1e-10,
    )

    devsim.finalize_mesh(mesh=mesh_name)
    devsim.create_device(mesh=mesh_name, device=device_name)

    # -- Material params on main SiC region --
    N_A_ionized = ionized_acceptor_concentration(N_A, T)
    mat_info = _set_sic_material_params(device_name, region_name, N_D_resolved, N_A, T)

    # -- Doping on main pillar --
    N_D_used = _apply_doping(
        device_name,
        region_name,
        junction_pos,
        N_A_ionized,
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )

    # -- Zero doping on trench fill --
    devsim.node_model(
        device=device_name,
        region="trench",
        name="Donors",
        equation="0",
    )
    devsim.node_model(
        device=device_name,
        region="trench",
        name="Acceptors",
        equation="0",
    )
    devsim.node_model(
        device=device_name,
        region="trench",
        name="NetDoping",
        equation="0",
    )
    # Set minimal material params on trench for solver stability
    _set_sic_material_params(device_name, "trench", N_D_resolved, N_A, T)

    # Sub-trench region gets substrate doping if it exists
    if trench_d_cm < total_depth - 1e-10:
        _set_sic_material_params(device_name, "sub_trench", N_D_resolved, N_A, T)
        # p+ substrate doping in sub_trench
        devsim.node_model(
            device=device_name,
            region="sub_trench",
            name="Acceptors",
            equation=f"{N_A_ionized}",
        )
        devsim.node_model(
            device=device_name,
            region="sub_trench",
            name="Donors",
            equation="0",
        )
        devsim.node_model(
            device=device_name,
            region="sub_trench",
            name="NetDoping",
            equation="Donors - Acceptors",
        )

    device_info = _build_device_info(
        device_name,
        region_name,
        junction_pos,
        epi_thickness_cm,
        substrate_thickness_cm,
        N_A,
        N_A_ionized,
        N_D_used,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
        pillar_hw_cm,
        mat_info,
        extra_keys={"structure_type": "mesa"},
    )

    logger.info(
        f"Created mesa device '{device_name}': {device_info['num_nodes']} nodes, "
        f"pillar_hw={pillar_half_width_um:.0f} um, trench={trench_width_um:.0f} um x "
        f"{trench_depth_um:.0f} um"
    )
    return device_info


# ---------------------------------------------------------------------------
# Public: 3D electrode (axisymmetric)
# ---------------------------------------------------------------------------


def _activate_cylindrical_coords(device_name, region_names):
    """Switch devsim from Cartesian to cylindrical (r,z) integration.

    Must be called AFTER create_device but BEFORE any equation setup.
    x-coordinate becomes radial (r), y-coordinate becomes axial (z).
    """
    devsim.set_parameter(name="raxis_zero", value=0.0)
    devsim.set_parameter(name="raxis_variable", value="x")
    for region in region_names:
        devsim.cylindrical_node_volume(device=device_name, region=region)
        devsim.cylindrical_edge_couple(device=device_name, region=region)
        devsim.cylindrical_surface_area(device=device_name, region=region)
    devsim.set_parameter(name="node_volume_model", value="CylindricalNodeVolume")
    devsim.set_parameter(name="edge_couple_model", value="CylindricalEdgeCouple")
    devsim.set_parameter(
        name="element_edge_couple_model",
        value="ElementCylindricalEdgeCouple",
    )
    devsim.set_parameter(
        name="element_node0_volume_model",
        value="ElementCylindricalNodeVolume@en0",
    )
    devsim.set_parameter(
        name="element_node1_volume_model",
        value="ElementCylindricalNodeVolume@en1",
    )


def restore_cartesian_coords():
    """Restore devsim global parameters to Cartesian defaults.

    Call after deleting a cylindrical device and before creating a new
    Cartesian device, to avoid solver interference.
    """
    devsim.set_parameter(name="node_volume_model", value="NodeVolume")
    devsim.set_parameter(name="edge_couple_model", value="EdgeCouple")
    devsim.set_parameter(
        name="element_edge_couple_model",
        value="ElementEdgeCouple",
    )
    devsim.set_parameter(
        name="element_node0_volume_model",
        value="ElementNodeVolume@en0",
    )
    devsim.set_parameter(
        name="element_node1_volume_model",
        value="ElementNodeVolume@en1",
    )


def create_3d_electrode_device(
    device_name="elec3d",
    region_name="sic",
    outer_radius_um=50.0,
    column_radius_um=5.0,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=1e-4,
    N_A=1e19,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
):
    """Create a 3D electrode SiC microdosimeter (axisymmetric).

    Uses cylindrical coordinates: x = radial (r), y = depth (z).
    Central n+ column at r=0..column_radius as cathode.
    Outer p+ at r=outer_radius as anode.

    Parameters
    ----------
    outer_radius_um : float
        Outer radius of the device (um).
    column_radius_um : float
        Radius of the central n+ electrode column (um).

    Returns device_info with ``structure_type="3d_electrode"``,
    ``coordinate_system="cylindrical"``.
    """
    N_D_resolved, _, _, _ = _resolve_doping_defaults(
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )
    outer_r_cm = outer_radius_um * 1e-4
    col_r_cm = column_radius_um * 1e-4
    junction_pos = substrate_thickness_cm
    total_depth = substrate_thickness_cm + epi_thickness_cm
    air_buffer = 1e-8

    mesh_name = f"{device_name}_mesh"
    devsim.create_2d_mesh(mesh=mesh_name)

    # -- Radial (x) mesh lines: fine near column --
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0, ps=1e-6)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=col_r_cm, ps=5e-6)
    # Intermediate radial points
    r_mid = (col_r_cm + outer_r_cm) * 0.5
    if r_mid > col_r_cm + 1e-5:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=r_mid, ps=2e-4)
    # AUDIT Mj-2 (v5): mesh line at the inner edge of the radial p+ anode wall
    # so the 1 um shell is actually resolved (otherwise it falls between nodes).
    wall_thickness_cm = 1e-4  # 1 um radial p+ shell at the outer anode wall
    devsim.add_2d_mesh_line(
        mesh=mesh_name, dir="x", pos=outer_r_cm - wall_thickness_cm, ps=1e-5
    )
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=outer_r_cm, ps=5e-6)

    # -- Depth (y) mesh lines --
    _add_depth_mesh_lines(mesh_name, junction_pos, total_depth, air_buffer)

    # -- Regions --
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_top",
        yl=-air_buffer,
        yh=0,
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region=region_name,
        yl=0,
        yh=total_depth,
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_bot",
        yl=total_depth,
        yh=total_depth + air_buffer,
    )

    # -- Contacts --
    # Cathode: top of central column (y=0, x=0..column_radius)
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="cathode",
        material="metal",
        region=region_name,
        xl=0,
        xh=col_r_cm,
        yl=0,
        yh=0,
        bloat=1e-10,
    )
    # Anode: outer surface (x=outer_radius, full y range)
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="anode",
        material="metal",
        region=region_name,
        xl=outer_r_cm,
        xh=outer_r_cm,
        yl=0,
        yh=total_depth,
        bloat=1e-10,
    )

    devsim.finalize_mesh(mesh=mesh_name)
    devsim.create_device(mesh=mesh_name, device=device_name)

    # CRITICAL: Activate cylindrical coordinates BEFORE physics
    all_regions = ["air_top", region_name, "air_bot"]
    _activate_cylindrical_coords(device_name, all_regions)

    # -- Material params --
    N_A_ionized = ionized_acceptor_concentration(N_A, T)
    mat_info = _set_sic_material_params(device_name, region_name, N_D_resolved, N_A, T)

    # -- Doping: standard p+/n- profile --
    N_D_used = _apply_doping(
        device_name,
        region_name,
        junction_pos,
        N_A_ionized,
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )

    # -- Central n+ column: high donor doping at x < column_radius --
    # AUDIT Mj-1 (v5): the previous code redefined Donors as
    #   max(Donors_column, N_D_used * step(y - junction_pos))
    # using the SCALAR N_D_used (=DEFAULT_N_D=1.07e15), which DISCARDED the
    # graded epi field installed by _apply_doping (true bulk ~8.5e13, ~12.6x
    # lower) -- corrupting depletion width, field shape, and V_fd for this
    # variant. Fix: rebuild the epi Donors field as Donors_epi from an EXPLICIT
    # expression (NOT `equation="Donors"`, which devsim evaluates lazily and
    # would create a Donors->Donors_epi->Donors cycle), then keep
    # max(column, epi). For uniform doping this reduces to N_D*step(y-junction);
    # for graded it preserves the real exponential profile.
    n_plus = 1e19
    _N_D_res, _N_D_j, _N_D_b, _L_t = _resolve_doping_defaults(
        None, doping_profile, N_D_junction, N_D_bulk, L_transition
    )
    if doping_profile == "graded":
        donors_epi_expr = (
            f"({_N_D_b} + ({_N_D_j} - {_N_D_b}) * "
            f"exp(-max(y - {junction_pos}, 0) / {_L_t})) "
            f"* step(y - {junction_pos})"
        )
    else:
        donors_epi_expr = f"{_N_D_res} * step(y - {junction_pos})"
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Donors_epi",
        equation=donors_epi_expr,
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Donors_column",
        equation=f"{n_plus} * step({col_r_cm} - x)",
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Donors",
        equation="max(Donors_column, Donors_epi)",
    )
    # AUDIT Mj-2 (v5): add a FULL-DEPTH radial p+ shell at the outer anode wall.
    # Previously acceptors were depth-only (step(junction_pos - y)), so the
    # "anode" contact at x=outer_r_cm was metal-on-n- over ~91% of its length
    # and the intended radial p+/n- junction did not exist. The p+ wall over the
    # full y-range creates the real radial junction the 3D-electrode geometry
    # needs. Use max() (not +) to avoid a doubled-doping corner where the wall
    # meets the horizontal p+ substrate. Must precede the NetDoping definition.
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors_wall",
        equation=f"{N_A_ionized} * step(x - {outer_r_cm - wall_thickness_cm})",
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors",
        equation=f"max({N_A_ionized} * step({junction_pos} - y), Acceptors_wall)",
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="NetDoping",
        equation="Donors - Acceptors",
    )

    device_info = _build_device_info(
        device_name,
        region_name,
        junction_pos,
        epi_thickness_cm,
        substrate_thickness_cm,
        N_A,
        N_A_ionized,
        N_D_used,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
        outer_r_cm,
        mat_info,
        extra_keys={
            "structure_type": "3d_electrode",
            "coordinate_system": "cylindrical",
            "outer_radius_cm": outer_r_cm,
            "column_radius_cm": col_r_cm,
        },
    )

    logger.info(
        f"Created 3D electrode device '{device_name}': "
        f"{device_info['num_nodes']} nodes, "
        f"outer_r={outer_radius_um:.0f} um, column_r={column_radius_um:.0f} um, "
        f"cylindrical coordinates active"
    )
    return device_info


# ---------------------------------------------------------------------------
# Public: Stacked delta-E/E telescope
# ---------------------------------------------------------------------------


def _setup_poisson_region(device_name, region_name, contacts):
    """Set up Poisson equation on a single region with given contacts.

    This is a minimal helper that replicates poisson.setup_poisson logic
    for arbitrary contact names, without modifying the poisson module.
    """
    from src.poisson import _create_sic_potential_only

    _create_sic_potential_only(device_name, region_name)
    for contact in contacts:
        bias_name = simple_physics.GetContactBiasName(contact)
        devsim.set_parameter(device=device_name, name=bias_name, value=0.0)
        simple_physics.CreateSiliconPotentialOnlyContact(
            device_name,
            region_name,
            contact,
        )


def create_delta_e_e_device(
    device_name="dee2d",
    region_name_de="delta_e",
    region_name_e="e_stop",
    half_width_um=50.0,
    delta_e_thickness_um=2.0,
    e_stop_thickness_um=50.0,
    substrate_thickness_cm=1e-4,
    N_A=1e19,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
):
    """Create a stacked delta-E/E telescope SiC microdosimeter.

    Two SiC layers with a devsim interface between them:
      - delta-E layer (thin, for particle identification)
      - E-stop layer (thick, for total energy measurement)

    Each layer has its own contact pair for independent readout:
      de_anode (top), estop_cathode (bottom).
    The interface between layers provides current continuity. No contacts
    are placed at the interface boundary (contacts would prevent devsim
    interface creation).

    Returns device_info with ``structure_type="delta_e_e"``, plus
    ``region_name_de``, ``region_name_e``, ``delta_e_thickness_cm``,
    ``e_stop_thickness_cm``.
    """
    N_D_resolved, _, _, _ = _resolve_doping_defaults(
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )
    half_width_cm = half_width_um * 1e-4
    de_thick_cm = delta_e_thickness_um * 1e-4
    e_thick_cm = e_stop_thickness_um * 1e-4
    air_buffer = 1e-8

    # Delta-E layer: y=0 to de_thick_cm
    # Junction in delta-E at substrate_thickness_cm from top
    de_junction = substrate_thickness_cm
    de_total = de_thick_cm

    # E-stop layer: y=de_thick_cm to de_thick_cm + e_thick_cm
    e_start = de_thick_cm
    e_junction = e_start + substrate_thickness_cm
    e_total = e_start + e_thick_cm

    total_depth = e_total

    mesh_name = f"{device_name}_mesh"
    devsim.create_2d_mesh(mesh=mesh_name)

    # -- Lateral (x) mesh lines --
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0, ps=5e-5)
    x_mid = half_width_cm * 0.5
    if x_mid > 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=x_mid, ps=3e-4)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=half_width_cm, ps=5e-4)

    # -- Depth (y) mesh lines --
    # Top air buffer
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=-air_buffer, ps=air_buffer)

    # Delta-E layer surface
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=0, ps=1e-5)

    # Delta-E junction
    if de_junction > 1e-8:
        if de_junction - 5e-6 > 1e-8:
            devsim.add_2d_mesh_line(
                mesh=mesh_name,
                dir="y",
                pos=de_junction - 5e-6,
                ps=1e-6,
            )
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=de_junction, ps=1e-7)

    # Delta-E / E-stop interface
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=de_total, ps=1e-6)

    # E-stop junction
    if e_junction - 5e-6 > e_start + 1e-8:
        devsim.add_2d_mesh_line(
            mesh=mesh_name,
            dir="y",
            pos=e_junction - 5e-6,
            ps=1e-6,
        )
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=e_junction, ps=1e-7)

    # E-stop intermediate points
    e_mid1 = e_junction + 2e-4
    e_mid2 = e_junction + 5e-4
    if e_mid1 < e_total - 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=e_mid1, ps=5e-6)
    if e_mid2 < e_total - 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=e_mid2, ps=5e-6)

    # E-stop bottom
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=e_total, ps=1e-5)

    # Bottom air buffer
    devsim.add_2d_mesh_line(
        mesh=mesh_name,
        dir="y",
        pos=e_total + air_buffer,
        ps=air_buffer,
    )

    # -- Regions --
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_top",
        yl=-air_buffer,
        yh=0,
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region=region_name_de,
        yl=0,
        yh=de_total,
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region=region_name_e,
        yl=de_total,
        yh=e_total,
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_bot",
        yl=e_total,
        yh=e_total + air_buffer,
    )

    # -- Interface between delta-E and E-stop --
    devsim.add_2d_interface(
        mesh=mesh_name,
        name="de_interface",
        region0=region_name_de,
        region1=region_name_e,
        xl=0,
        xh=half_width_cm,
        yl=de_total,
        yh=de_total,
        bloat=1e-10,
    )

    # -- Contacts --
    # Delta-E anode at top surface
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="de_anode",
        material="metal",
        region=region_name_de,
        yl=0,
        yh=0,
        bloat=1e-10,
    )
    # E-stop cathode at bottom surface
    # (no contacts at interface -- would prevent devsim interface creation)
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="estop_cathode",
        material="metal",
        region=region_name_e,
        yl=e_total,
        yh=e_total,
        bloat=1e-10,
    )

    devsim.finalize_mesh(mesh=mesh_name)
    devsim.create_device(mesh=mesh_name, device=device_name)

    # -- Material params on both regions --
    N_A_ionized = ionized_acceptor_concentration(N_A, T)
    mat_info_de = _set_sic_material_params(
        device_name,
        region_name_de,
        N_D_resolved,
        N_A,
        T,
    )
    _set_sic_material_params(device_name, region_name_e, N_D_resolved, N_A, T)

    # -- Doping on delta-E region --
    N_D_used = _apply_doping(
        device_name,
        region_name_de,
        de_junction,
        N_A_ionized,
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )

    # -- Doping on E-stop region (junction relative to e_start) --
    _apply_doping(
        device_name,
        region_name_e,
        e_junction,
        N_A_ionized,
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )

    # Count nodes across primary region
    num_nodes_de = len(
        devsim.get_node_model_values(
            device=device_name,
            region=region_name_de,
            name="x",
        )
    )
    num_nodes_e = len(
        devsim.get_node_model_values(
            device=device_name,
            region=region_name_e,
            name="x",
        )
    )

    device_info = {
        "device_name": device_name,
        "region_name": region_name_de,  # primary region for pipeline compat
        "region_name_de": region_name_de,
        "region_name_e": region_name_e,
        "junction_pos": de_junction,
        "epi_thickness_cm": de_thick_cm - substrate_thickness_cm,
        "substrate_thickness_cm": substrate_thickness_cm,
        "total_length": total_depth,
        "N_D": N_D_used,
        "N_A": N_A,
        "N_A_ionized": N_A_ionized,
        "T": T,
        "n_i": mat_info_de["n_i"],
        "E_g": mat_info_de["E_g"],
        "params": mat_info_de["params"],
        "mu_n": mat_info_de["mu_n"],
        "mu_p": mat_info_de["mu_p"],
        "num_nodes": num_nodes_de + num_nodes_e,
        "doping_profile": doping_profile,
        "N_D_junction": N_D_junction,
        "N_D_bulk": N_D_bulk,
        "L_transition": L_transition,
        "half_width_cm": half_width_cm,
        "dimension": 2,
        "structure_type": "delta_e_e",
        "delta_e_thickness_cm": de_thick_cm,
        "e_stop_thickness_cm": e_thick_cm,
    }

    logger.info(
        f"Created delta-E/E device '{device_name}': "
        f"{num_nodes_de + num_nodes_e} nodes "
        f"(dE={num_nodes_de}, E={num_nodes_e}), "
        f"delta_e={delta_e_thickness_um:.0f} um, e_stop={e_stop_thickness_um:.0f} um"
    )
    return device_info


# ---------------------------------------------------------------------------
# Public: Guard ring
# ---------------------------------------------------------------------------


def create_guard_ring_device(
    device_name="gr2d",
    region_name="sic",
    sv_half_width_um=50.0,
    guard_ring_width_um=5.0,
    guard_ring_gap_um=3.0,
    guard_ring_depth_um=1.0,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=1e-4,
    N_A=1e19,
    N_A_guard=5e18,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
):
    """Create a planar SiC microdosimeter with p+ guard ring.

    Extended planar device with a guard ring beyond the SV edge.

    Lateral layout (half-device, x-direction):
        [0, sv_hw]                   = SV active area
        [sv_hw, sv_hw+gap]           = transition region
        [sv_hw+gap, sv_hw+gap+gr_w]  = guard ring

    Guard ring: additional p+ implant at guard ring position/depth,
    modeled via acceptor doping overlay.

    Contacts: anode at y=0 (on SV), cathode at y=total_depth,
    guard_ring_anode at y=0 on guard ring position.

    Returns device_info with ``structure_type="guard_ring"``.
    """
    N_D_resolved, _, _, _ = _resolve_doping_defaults(
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )
    sv_hw_cm = sv_half_width_um * 1e-4
    gr_w_cm = guard_ring_width_um * 1e-4
    gr_gap_cm = guard_ring_gap_um * 1e-4
    gr_depth_cm = guard_ring_depth_um * 1e-4
    junction_pos = substrate_thickness_cm
    total_depth = substrate_thickness_cm + epi_thickness_cm
    air_buffer = 1e-8

    total_hw_cm = sv_hw_cm + gr_gap_cm + gr_w_cm
    gr_inner = sv_hw_cm + gr_gap_cm
    gr_outer = gr_inner + gr_w_cm

    mesh_name = f"{device_name}_mesh"
    devsim.create_2d_mesh(mesh=mesh_name)

    # -- Lateral (x) mesh lines --
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0, ps=5e-5)
    # SV midpoint
    if sv_hw_cm * 0.5 > 1e-6:
        devsim.add_2d_mesh_line(
            mesh=mesh_name,
            dir="x",
            pos=sv_hw_cm * 0.5,
            ps=2e-4,
        )
    # SV edge (fine)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=sv_hw_cm, ps=1e-5)
    # Gap region
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=gr_inner, ps=1e-5)
    # Guard ring midpoint
    gr_mid = gr_inner + gr_w_cm * 0.5
    if gr_w_cm > 2e-5:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=gr_mid, ps=5e-5)
    # Guard ring outer edge
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=gr_outer, ps=5e-5)

    # -- Depth (y) mesh lines --
    _add_depth_mesh_lines(mesh_name, junction_pos, total_depth, air_buffer)

    # Add mesh line at guard ring depth if distinct
    if gr_depth_cm > 1e-8 and abs(gr_depth_cm - junction_pos) > 1e-8:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=gr_depth_cm, ps=1e-6)

    # -- Regions --
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_top",
        xl=0,
        xh=total_hw_cm,
        yl=-air_buffer,
        yh=0,
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region=region_name,
        xl=0,
        xh=total_hw_cm,
        yl=0,
        yh=total_depth,
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_bot",
        xl=0,
        xh=total_hw_cm,
        yl=total_depth,
        yh=total_depth + air_buffer,
    )

    # -- Contacts --
    # Main anode on SV (x=0..sv_hw, y=0)
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="anode",
        material="metal",
        region=region_name,
        xl=0,
        xh=sv_hw_cm,
        yl=0,
        yh=0,
        bloat=1e-10,
    )
    # Cathode at bottom (full width)
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="cathode",
        material="metal",
        region=region_name,
        yl=total_depth,
        yh=total_depth,
        bloat=1e-10,
    )
    # Guard ring contact
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="guard_ring_anode",
        material="metal",
        region=region_name,
        xl=gr_inner,
        xh=gr_outer,
        yl=0,
        yh=0,
        bloat=1e-10,
    )

    devsim.finalize_mesh(mesh=mesh_name)
    devsim.create_device(mesh=mesh_name, device=device_name)

    # -- Material params --
    N_A_ionized = ionized_acceptor_concentration(N_A, T)
    mat_info = _set_sic_material_params(device_name, region_name, N_D_resolved, N_A, T)

    # -- Standard doping profile --
    N_D_used = _apply_doping(
        device_name,
        region_name,
        junction_pos,
        N_A_ionized,
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )

    # -- Guard ring: additional p+ doping overlay --
    # Guard ring acceptor contribution at (gr_inner < x < gr_outer, y < gr_depth)
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors_GR",
        equation=(
            f"{N_A_guard} * step(x - {gr_inner}) * step({gr_outer} - x) "
            f"* step({gr_depth_cm} - y)"
        ),
    )
    # Redefine Acceptors: substrate doping + guard ring doping (no self-reference)
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors",
        equation=f"{N_A_ionized} * step({junction_pos} - y) + Acceptors_GR",
    )
    # Recalculate NetDoping
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="NetDoping",
        equation="Donors - Acceptors",
    )

    device_info = _build_device_info(
        device_name,
        region_name,
        junction_pos,
        epi_thickness_cm,
        substrate_thickness_cm,
        N_A,
        N_A_ionized,
        N_D_used,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
        sv_hw_cm,
        mat_info,
        extra_keys={
            "structure_type": "guard_ring",
            "guard_ring_contact": "guard_ring_anode",
            "sv_half_width_cm": sv_hw_cm,
            "total_half_width_cm": total_hw_cm,
        },
    )

    logger.info(
        f"Created guard ring device '{device_name}': "
        f"{device_info['num_nodes']} nodes, "
        f"sv_hw={sv_half_width_um:.0f} um, "
        f"gr: gap={guard_ring_gap_um:.0f} um, width={guard_ring_width_um:.0f} um, "
        f"depth={guard_ring_depth_um:.0f} um, N_A_guard={N_A_guard:.1e}"
    )
    return device_info
