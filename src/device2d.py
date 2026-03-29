"""devsim 2D device setup for 4H-SiC p+/n- diode.

Creates a 2D triangular mesh with non-uniform spacing for planar SiC
microdosimeter geometries. Applies 4H-SiC material parameters and graded
(or uniform) epi doping profiles.

Coordinate convention:
    x = lateral direction (0 to half_width_cm), exploiting symmetry
    y = depth direction (0 = anode/top surface, total_depth = cathode/bottom)

The 1D device.py uses x for depth; here depth is mapped to y.
Lateral direction is uniform in doping by construction.

Air buffer regions are created at the top (y < 0) and bottom (y > total_depth)
to enable devsim contact node detection.  Without these buffers, contacts
have zero nodes and the solver fails.

All units CGS (cm, cm^-3, F/cm, eV, s) per devsim convention.

References:
    - devsim documentation: https://devsim.net/
    - devsim mos_2d_create.py example pattern
    - Petringa et al. experimental device structure
"""

import logging

import devsim

from src.sic_material import (
    SiC4H_Parameters,
    intrinsic_concentration,
    mobility_caughey_thomas_T,
    srh_lifetime,
)
from src.incomplete_ionization import ionized_acceptor_concentration

logger = logging.getLogger(__name__)

# Calibrated graded doping defaults (same as device.py)
_N_D_JUNCTION_DEFAULT = 2.9e15  # cm^-3, near junction
_N_D_BULK_DEFAULT = 8.5e13  # cm^-3, deep epi bulk
_L_TRANSITION_DEFAULT = 1e-4  # cm, characteristic decay length

# Physical constants (CGS)
EPS_0 = 8.854e-14  # F/cm, vacuum permittivity
Q = 1.602e-19  # C, elementary charge
K_B_EV = 8.617e-5  # eV/K, Boltzmann constant
K_B_J = 1.3806503e-23  # J/K, Boltzmann constant

# Default uniform donor concentration (same as device.py)
DEFAULT_N_D = 1.07e15  # cm^-3


def set_doping_profile_2d(device_name, region_name, junction_pos, N_A_ionized, N_D):
    """Set step-function (uniform) doping profile on a 2D device.

    Acceptors in p+ substrate (y < junction_pos), donors in n- epi (y >= junction_pos).
    Doping varies only with y (depth); laterally uniform by construction.

    Parameters
    ----------
    device_name : str
        devsim device name.
    region_name : str
        devsim region name.
    junction_pos : float
        Position of metallurgical junction along y-axis (cm).
    N_A_ionized : float
        Ionized acceptor concentration (cm^-3).
    N_D : float
        Donor concentration in epi (cm^-3).
    """
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors",
        equation=f"{N_A_ionized}*step({junction_pos}-y)",
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Donors",
        equation=f"{N_D}*step(y-{junction_pos})",
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="NetDoping",
        equation="Donors-Acceptors",
    )


def set_graded_doping_2d(
    device_name,
    region_name,
    junction_pos,
    N_A_ionized,
    N_D_junction,
    N_D_bulk,
    L_transition,
):
    """Set exponentially graded N_D(y) doping profile on a 2D device.

    N_D(y) = N_D_bulk + (N_D_junction - N_D_bulk) * exp(-(y - y_j) / L)
    for y > junction_pos (n-side).  Acceptors are a step function in the
    p+ substrate.  Doping is laterally uniform (no x-dependence).

    Parameters
    ----------
    device_name : str
        devsim device name.
    region_name : str
        devsim region name.
    junction_pos : float
        Position of metallurgical junction along y-axis (cm).
    N_A_ionized : float
        Ionized acceptor concentration (cm^-3).
    N_D_junction : float
        Donor concentration near junction (cm^-3).
    N_D_bulk : float
        Donor concentration deep in epi bulk (cm^-3).
    L_transition : float
        Characteristic decay length of the doping gradient (cm).
    """
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors",
        equation=f"{N_A_ionized}*step({junction_pos}-y)",
    )

    donor_expr = (
        f"({N_D_bulk} + ({N_D_junction} - {N_D_bulk}) * "
        f"exp(-max(y - {junction_pos}, 0) / {L_transition})) "
        f"* step(y - {junction_pos})"
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Donors",
        equation=donor_expr,
    )

    devsim.node_model(
        device=device_name,
        region=region_name,
        name="NetDoping",
        equation="Donors - Acceptors",
    )


def create_sic_2d_device(
    device_name="sic2d",
    region_name="sic",
    half_width_um=50.0,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=1e-4,
    N_A=1e19,
    N_D=None,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
):
    """Create a 2D p+/n- 4H-SiC diode device in devsim.

    Structure uses symmetry: x = 0 (center/symmetry plane) to half_width_cm.
    Depth: y = 0 (anode, top) to total_depth (cathode, bottom).
    Air buffer regions are added at top and bottom for contact detection.

    Parameters
    ----------
    device_name : str
        Name for the devsim device.
    region_name : str
        Name for the semiconductor region.
    half_width_um : float
        Half-width of the sensitive volume in micrometers. Default 50 um
        (for 100 um SV).
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default 10 um = 10e-4 cm.
    substrate_thickness_cm : float
        p+ substrate thickness (cm). Default 1 um = 1e-4 cm.
    N_A : float
        Total acceptor concentration in p+ substrate (cm^-3).
    N_D : float or None
        Donor concentration in n- epi (cm^-3).  For uniform profile only.
        If None, uses DEFAULT_N_D.
    T : float
        Temperature (K).
    doping_profile : str
        "graded" (default) for exponentially graded N_D(y), or
        "uniform" for step-function N_D.
    N_D_junction : float or None
        For graded profile: donor concentration near junction (cm^-3).
    N_D_bulk : float or None
        For graded profile: donor concentration in bulk epi (cm^-3).
    L_transition : float or None
        For graded profile: characteristic decay length (cm).

    Returns
    -------
    device_info : dict
        Dictionary with device metadata, material parameters, and mesh info.
        Includes all 1D keys plus half_width_cm and dimension=2.
    """
    if N_D is None:
        N_D = DEFAULT_N_D

    params = SiC4H_Parameters()
    half_width_cm = half_width_um * 1e-4  # um -> cm
    junction_pos = substrate_thickness_cm
    total_depth = substrate_thickness_cm + epi_thickness_cm
    air_buffer = 1e-8  # cm, thin buffer for contact detection

    # --- Create 2D mesh ---
    mesh_name = f"{device_name}_mesh"
    devsim.create_2d_mesh(mesh=mesh_name)

    # -- Lateral (x) mesh lines: fine at center, coarser toward edge --
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0, ps=5e-5)
    # Intermediate lateral points for smoother mesh
    x_quarter = half_width_cm * 0.25
    x_half = half_width_cm * 0.5
    x_three_quarter = half_width_cm * 0.75
    if x_quarter > 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=x_quarter, ps=2e-4)
    if x_half > 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=x_half, ps=3e-4)
    if x_three_quarter > 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=x_three_quarter, ps=4e-4)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=half_width_cm, ps=5e-4)

    # -- Depth (y) mesh lines: mirrors 1D spacing from device.py --
    # Top air buffer
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=-air_buffer, ps=air_buffer)

    # Anode surface (y=0)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=0, ps=1e-5)

    # Approach junction from p-side
    if junction_pos - 5e-6 > 1e-8:
        devsim.add_2d_mesh_line(
            mesh=mesh_name, dir="y", pos=junction_pos - 5e-6, ps=1e-6
        )

    # Junction: very fine spacing
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=junction_pos, ps=1e-7)

    # Epi intermediate points (same as 1D)
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
        mesh=mesh_name, dir="y", pos=total_depth + air_buffer, ps=air_buffer
    )

    # -- Regions --
    # Air buffer top (same material as main region per research finding)
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_top",
        yl=-air_buffer,
        yh=0,
    )
    # Main SiC region
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region=region_name,
        yl=0,
        yh=total_depth,
    )
    # Air buffer bottom
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_bot",
        yl=total_depth,
        yh=total_depth + air_buffer,
    )

    # -- Contacts --
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="anode",
        material="metal",
        region=region_name,
        yl=0,
        yh=0,
        bloat=1e-10,
    )
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

    # --- Set 4H-SiC material parameters ---
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

    # Doping- and temperature-dependent mobility
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

    # --- Set doping profile ---
    N_A_ionized = ionized_acceptor_concentration(N_A, T)

    # Resolve graded doping defaults
    _N_D_junction = N_D_junction if N_D_junction is not None else _N_D_JUNCTION_DEFAULT
    _N_D_bulk = N_D_bulk if N_D_bulk is not None else _N_D_BULK_DEFAULT
    _L_transition = L_transition if L_transition is not None else _L_TRANSITION_DEFAULT

    if doping_profile == "graded":
        set_graded_doping_2d(
            device_name,
            region_name,
            junction_pos,
            N_A_ionized,
            _N_D_junction,
            _N_D_bulk,
            _L_transition,
        )
        logger.info(
            f"Graded 2D doping: N_D_junction={_N_D_junction:.2e}, "
            f"N_D_bulk={_N_D_bulk:.2e}, L_transition={_L_transition:.2e} cm"
        )
    else:
        set_doping_profile_2d(
            device_name,
            region_name,
            junction_pos,
            N_A_ionized,
            N_D,
        )

    num_nodes = len(
        devsim.get_node_model_values(device=device_name, region=region_name, name="x")
    )
    logger.info(
        f"Created 2D device '{device_name}': {num_nodes} nodes, "
        f"half_width={half_width_um:.0f} um, "
        f"junction at y={junction_pos*1e4:.1f} um, "
        f"N_A_ionized={N_A_ionized:.2e}, profile={doping_profile}"
    )

    return {
        "device_name": device_name,
        "region_name": region_name,
        "junction_pos": junction_pos,
        "epi_thickness_cm": epi_thickness_cm,
        "substrate_thickness_cm": substrate_thickness_cm,
        "total_length": total_depth,
        "N_D": N_D,
        "N_A": N_A,
        "N_A_ionized": N_A_ionized,
        "T": T,
        "n_i": n_i_T,
        "E_g": E_g_T,
        "params": params,
        "mu_n": mu_n,
        "mu_p": mu_p,
        "num_nodes": num_nodes,
        "doping_profile": doping_profile,
        "N_D_junction": N_D_junction,
        "N_D_bulk": N_D_bulk,
        "L_transition": L_transition,
        "half_width_cm": half_width_cm,
        "dimension": 2,
    }
