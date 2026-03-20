"""devsim 1D device setup for 4H-SiC p+/n- diode.

Creates a 1D mesh with non-uniform spacing, applies 4H-SiC material parameters,
and sets up doping profiles for the p+ substrate / n- epitaxial layer structure.

All units CGS (cm, cm^-3, F/cm, eV, s) per devsim convention.

References:
    - devsim documentation: https://devsim.net/
    - devsim diode example pattern
    - Petringa et al. experimental device structure
"""

import logging

import devsim
import numpy as np

from src.sic_material import SiC4H_Parameters, mobility_caughey_thomas
from src.incomplete_ionization import ionized_acceptor_concentration

logger = logging.getLogger(__name__)

# Calibrated donor concentration from Plan 01 (W(0V) = 1.7 um target)
DEFAULT_N_D = 1.07e15  # cm^-3

# Physical constants (CGS)
EPS_0 = 8.854e-14  # F/cm, vacuum permittivity
Q = 1.602e-19  # C, elementary charge
K_B_EV = 8.617e-5  # eV/K, Boltzmann constant
K_B_J = 1.3806503e-23  # J/K, Boltzmann constant


def create_sic_device(
    device_name="sic_diode",
    region_name="sic",
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=1e-4,
    N_A=1e19,
    N_D=None,
    T=300,
):
    """Create a 1D p+/n- 4H-SiC diode device in devsim.

    Structure: p+ substrate (left, anode) | n- epitaxial layer (right, cathode)

    Parameters
    ----------
    device_name : str
        Name for the devsim device.
    region_name : str
        Name for the semiconductor region.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default 10 um = 10e-4 cm.
    substrate_thickness_cm : float
        p+ substrate thickness (cm). Default 1 um = 1e-4 cm.
    N_A : float
        Total acceptor concentration in p+ substrate (cm^-3).
    N_D : float or None
        Donor concentration in n- epi (cm^-3). If None, uses DEFAULT_N_D.
    T : float
        Temperature (K).

    Returns
    -------
    device_info : dict
        Dictionary with device_name, region_name, junction_pos, epi_thickness_cm,
        N_D, N_A_ionized, params (SiC4H_Parameters instance), and other metadata.
    """
    if N_D is None:
        N_D = DEFAULT_N_D

    params = SiC4H_Parameters()
    junction_pos = substrate_thickness_cm
    total_length = substrate_thickness_cm + epi_thickness_cm

    # --- Create 1D mesh with non-uniform spacing ---
    mesh_name = f"{device_name}_mesh"
    devsim.create_1d_mesh(mesh=mesh_name)

    # p+ substrate: x = 0 (anode contact) to junction_pos
    # Coarse spacing in bulk (~100 nm = 1e-5 cm)
    devsim.add_1d_mesh_line(mesh=mesh_name, pos=0.0, ps=1e-5, tag="top")

    # Approach junction from p-side: finer spacing
    devsim.add_1d_mesh_line(mesh=mesh_name, pos=junction_pos - 5e-6, ps=1e-6)

    # Junction vicinity: very fine spacing (~1 nm = 1e-7 cm)
    devsim.add_1d_mesh_line(mesh=mesh_name, pos=junction_pos, ps=1e-7)

    # Depletion region in n-side: fine spacing near junction, coarser in bulk
    # First ~2 um from junction: fine (~10 nm)
    devsim.add_1d_mesh_line(mesh=mesh_name, pos=junction_pos + 2e-4, ps=5e-6)

    # Bulk n- epi: medium spacing (~50 nm = 5e-6 cm)
    devsim.add_1d_mesh_line(mesh=mesh_name, pos=junction_pos + 5e-4, ps=5e-6)

    # Far end of epi: coarser spacing
    devsim.add_1d_mesh_line(mesh=mesh_name, pos=total_length, ps=1e-5, tag="bot")

    # Contacts
    devsim.add_1d_contact(mesh=mesh_name, name="anode", tag="top", material="metal")
    devsim.add_1d_contact(mesh=mesh_name, name="cathode", tag="bot", material="metal")

    # Region spanning entire device
    devsim.add_1d_region(
        mesh=mesh_name, material="SiC", region=region_name, tag1="top", tag2="bot"
    )

    devsim.finalize_mesh(mesh=mesh_name)
    devsim.create_device(mesh=mesh_name, device=device_name)

    # --- Set 4H-SiC material parameters ---
    kT_eV = K_B_EV * T  # eV
    kT_J = K_B_J * T  # J
    V_t = kT_J / Q  # V (thermal voltage in Volts for devsim)

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
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="n_i",
        value=params.n_i_300,
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

    # Doping-dependent mobility
    mu_n = mobility_caughey_thomas(N_D, carrier="electron")
    mu_p = mobility_caughey_thomas(N_A, carrier="hole")
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

    # SRH parameters (for potential later use)
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="n1",
        value=params.n_i_300,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="p1",
        value=params.n_i_300,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="taun",
        value=params.tau_n,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="taup",
        value=params.tau_p,
    )

    # --- Set doping profile ---
    N_A_ionized = ionized_acceptor_concentration(N_A, T)
    set_doping_profile(device_name, region_name, junction_pos, N_A_ionized, N_D)

    num_nodes = len(
        devsim.get_node_model_values(device=device_name, region=region_name, name="x")
    )
    logger.info(
        f"Created device '{device_name}': {num_nodes} nodes, "
        f"junction at x={junction_pos*1e4:.1f} um, "
        f"N_A_ionized={N_A_ionized:.2e}, N_D={N_D:.2e}"
    )

    return {
        "device_name": device_name,
        "region_name": region_name,
        "junction_pos": junction_pos,
        "epi_thickness_cm": epi_thickness_cm,
        "substrate_thickness_cm": substrate_thickness_cm,
        "total_length": total_length,
        "N_D": N_D,
        "N_A": N_A,
        "N_A_ionized": N_A_ionized,
        "T": T,
        "params": params,
        "mu_n": mu_n,
        "mu_p": mu_p,
        "num_nodes": num_nodes,
    }


def set_doping_profile(device_name, region_name, junction_pos, N_A_ionized, N_D):
    """Set step-function doping profile on the device.

    p-side (x < junction_pos): Acceptors = N_A_ionized
    n-side (x >= junction_pos): Donors = N_D
    NetDoping = Donors - Acceptors

    Parameters
    ----------
    device_name : str
        devsim device name.
    region_name : str
        devsim region name.
    junction_pos : float
        Position of metallurgical junction (cm).
    N_A_ionized : float
        Ionized acceptor concentration (cm^-3).
    N_D : float
        Donor concentration (cm^-3).
    """
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors",
        equation=f"{N_A_ionized}*step({junction_pos}-x)",
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Donors",
        equation=f"{N_D}*step(x-{junction_pos})",
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="NetDoping",
        equation="Donors-Acceptors",
    )
