"""devsim 1D device setup for 4H-SiC p+/n- diode.

Creates a 1D mesh with non-uniform spacing, applies 4H-SiC material parameters,
and sets up doping profiles for the p+ substrate / n- epitaxial layer structure.
Supports both uniform and graded (exponential) epitaxial doping profiles.

All units CGS (cm, cm^-3, F/cm, eV, s) per devsim convention.

References:
    - devsim documentation: https://devsim.net/
    - devsim diode example pattern
    - reference experimental device structure
"""

import logging

import devsim
import numpy as np
import scipy.optimize

from etna.core.sic_material import (
    SiC4H_Parameters,
    mobility_caughey_thomas,
    intrinsic_concentration,
    mobility_caughey_thomas_T,
    srh_lifetime,
)
from etna.core.incomplete_ionization import ionized_acceptor_concentration

logger = logging.getLogger(__name__)

# Calibrated donor concentration from Plan 01 (W(0V) = 1.7 um target)
DEFAULT_N_D = 1.07e15  # cm^-3
# n+ layer under the epi (when enabled): conductive 4H-SiC substrate doping, fully
# ionized. 1e18 and 1e19 give the same W/CCE/dark current to <1e-3 (lab spike).
DEFAULT_N_D_NPLUS = 1e19  # cm^-3

# DEFAULT DEVICE STACK (single source of truth): p+ anode 0.3 um / 10 um graded n- epi /
# n+ substrate (1 um modelled, cathode on its far end), confirmed by the device owner.
# Read by DeviceConfig, create_sic_device / create_sic_2d_device (graded-doping fallbacks),
# cce_vs_bias and the 2D engine gate; SiC4H_Parameters holds N_t / E_t (read back here).
# Every core default reads it: the builders' geometry (p+ 0.3 um, n+ 1 um) and graded
# doping, and the sweep helpers (transient, flash, temperature, C-V, dark current, CCE).
# Explicit arguments win; nplus_thickness_cm=None (explicit) builds the no-n+ legacy
# geometry. Before 2026-10-17 the builders kept p+ 1 um / no n+ and the helpers hardcoded
# 2.90e15 / 8.50e13 / 1.0e-4 (a third doping set), so direct core calls mixed models.
# Calibration (lab/solver2d/bench/nplus_calibrate.py; pinned by
# tests/test_nplus_calibration.py; notes in lab/solver2d/notes/LOG.md):
# C-V: the three W targets (1.7 / 9.5 / 9.73 um at 0 / -10 / -30 V, A = 4 mm2; Photons
# Fig. 4b labels 1.68 / 9.52 / 9.73 um) are W = eps*A/C of a MEASURED C, so the model
# observable is eps*A/C of the small-signal C = dQ/dV that run_cv returns (NOT the
# n/N_D = 0.5 depletion edge; before 2026-10-17 the calibration matched that edge and
# missed its own -10 V point by +2.8 % in eps*A/C: lab/microdosimeter/notes/LOG.md
# sections 4 and 6). The three targets have NO exact solution for (N_D_junction,
# N_D_bulk, L) in this 1D model (3-parameter least squares: residual ~3 %). Hence:
# L_transition held at 0.987 um and (N_D_junction, N_D_bulk) solved EXACTLY on eps*A/C
# at 0 V and -10 V (2 equations, 2 unknowns; residuals < 1e-14 engine, < 3e-8 devsim).
# These two points are therefore CALIBRATION INPUTS, not reproduced validation points.
# The -30 V point is out of fit: +2.56 % (9.979 vs 9.73 um). Three points cannot tell
# apart a ~0.5 um graded n-/n+ buffer, an epi ~2.5 % thinner than nominal (9.75 um, inside
# typical tolerance) and lateral edge/fringe capacitance falling with bias (LOG.md
# section 4: each alone closes the gap); the C-V knee (-5..-15 V) or a SIMS profile could.
# Not picked, not refitted. Previous record (edge calibration): 2.974271e15 / 6.335936e13.
# N_t: rescaled (leakage linear in N_t) so that I(300 K, -60 V, 4 mm2) = 18 pA, the bound
# "< 18 pA up to -60 V" used as equality; N_t is the fitted magnitude of the
# f_depl = min(E/3.5e4, 1) effective-generation PROXY, not a trap density.
# E_t = 0.67 eV: Z1/2 level of radiation_damage (Burin et al. Table I).
DEFAULT_STACK = {
    "substrate_thickness_cm": 0.3e-4,  # p+ anode layer (Petringa 2025 Fig. 1)
    "nplus_thickness_cm": 1e-4,
    "N_D_nplus": DEFAULT_N_D_NPLUS,
    "N_D_junction": 2.958237e15,
    "N_D_bulk": 8.252732e13,
    "L_transition": 0.987e-4,  # held, not fitted
    "N_t": SiC4H_Parameters.N_t,  # 2.805336e12
    "E_t": SiC4H_Parameters.E_t,  # 0.67
}
# Historical name of the same record (it was the opt-in calibration before the flip).
NPLUS_CALIBRATION = DEFAULT_STACK

# LEGACY model (the default before the n+ flip), kept so old numbers stay reproducible:
# cathode contact directly on the n- epi (no n+; non-convergence wall at deep reverse bias,
# a contact artefact), p+ 1 um, the old DeviceConfig doping, N_t 2.8e12, E_t 0.65.
# etna.api.device.legacy_config() builds the DeviceConfig; pinned by tests/test_default_stack.py.
LEGACY_STACK = {
    "substrate_thickness_cm": 1e-4,
    "nplus_thickness_cm": None,
    "N_D_nplus": DEFAULT_N_D_NPLUS,
    "N_D_junction": 2.93e15,
    "N_D_bulk": 8.82e13,
    "L_transition": 0.987e-4,
    "N_t": 2.8e12,
    "E_t": 0.65,
}

# Physical constants (CGS)
EPS_0 = 8.854e-14  # F/cm, vacuum permittivity
Q = 1.602e-19  # C, elementary charge
K_B_EV = 8.617e-5  # eV/K, Boltzmann constant
K_B_J = 1.3806503e-23  # J/K, Boltzmann constant


def build_diode_mesh(
    device_name, region_name, junction_pos, total_length, nplus_thickness_cm=None
):
    """Build and finalize the graded 1D mesh of the p+/n- diode (mesh ``<device_name>_mesh``).

    Anode contact at x = 0, one region spanning the device. Without an n+ layer
    (``nplus_thickness_cm=None``, the historical model) the cathode is at
    ``total_length`` (= the epi end). With an n+ layer the epi end ``total_length``
    is refined (1 nm, abrupt n-/n+ step) and the mesh continues through the n+ to
    the cathode at ``total_length + nplus_thickness_cm``. Shared by
    ``create_sic_device`` and the batched engine's mesh probe so both use the
    identical discretisation.
    """
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
    # Only add intermediate mesh points if they fall strictly within the epi layer
    epi_mid1 = junction_pos + 2e-4  # ~2 um from junction
    epi_mid2 = junction_pos + 5e-4  # ~5 um from junction

    if epi_mid1 < total_length - 1e-6:
        devsim.add_1d_mesh_line(mesh=mesh_name, pos=epi_mid1, ps=5e-6)

    if epi_mid2 < total_length - 1e-6:
        devsim.add_1d_mesh_line(mesh=mesh_name, pos=epi_mid2, ps=5e-6)

    if nplus_thickness_cm is None:
        # Far end of epi: coarser spacing
        devsim.add_1d_mesh_line(mesh=mesh_name, pos=total_length, ps=1e-5, tag="bot")
    else:
        # n-/n+ step at the epi end, then the n+ layer up to the cathode
        devsim.add_1d_mesh_line(mesh=mesh_name, pos=total_length, ps=1e-7)
        devsim.add_1d_mesh_line(
            mesh=mesh_name, pos=total_length + nplus_thickness_cm, ps=1e-5, tag="bot"
        )

    # Contacts
    devsim.add_1d_contact(mesh=mesh_name, name="anode", tag="top", material="metal")
    devsim.add_1d_contact(mesh=mesh_name, name="cathode", tag="bot", material="metal")

    # Region spanning entire device
    devsim.add_1d_region(
        mesh=mesh_name, material="SiC", region=region_name, tag1="top", tag2="bot"
    )

    devsim.finalize_mesh(mesh=mesh_name)


def create_sic_device(
    device_name="sic_diode",
    region_name="sic",
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=DEFAULT_STACK["substrate_thickness_cm"],
    N_A=1e19,
    N_D=None,
    T=300,
    doping_profile="uniform",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
    nplus_thickness_cm=DEFAULT_STACK["nplus_thickness_cm"],
    N_D_nplus=DEFAULT_STACK["N_D_nplus"],
    nplus_generation=False,
):
    """Create a 1D p+/n- 4H-SiC diode device in devsim.

    Structure: p+ anode layer (left, anode) | n- epitaxial layer | [n+ layer] (right, cathode)

    Parameters
    ----------
    device_name : str
        Name for the devsim device.
    region_name : str
        Name for the semiconductor region.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default 10 um = 10e-4 cm.
    substrate_thickness_cm : float
        p+ ANODE layer thickness (cm). Default DEFAULT_STACK (0.3 um); the legacy model
        used 1 um = 1e-4 cm. (Historical name: this is not the n+ substrate, see
        ``nplus_thickness_cm``.)
    N_A : float
        Total acceptor concentration in p+ substrate (cm^-3).
    N_D : float or None
        Donor concentration in n- epi (cm^-3). If None, uses DEFAULT_N_D.
        For graded profile, used as fallback if N_D_junction not specified.
    T : float
        Temperature (K).
    doping_profile : str
        "uniform" (default) for step-function N_D, or "graded" for
        exponentially graded N_D(x) in the epi layer.
    N_D_junction : float or None
        For graded profile: donor concentration near junction (cm^-3).
    N_D_bulk : float or None
        For graded profile: donor concentration in bulk epi (cm^-3).
    L_transition : float or None
        For graded profile: characteristic decay length (cm).
    nplus_thickness_cm : float or None
        Thickness (cm) of an n+ layer under the n- epi, cathode contact on its far end
        (the conductive n+ substrate of the real device, truncated). Default
        DEFAULT_STACK (1 um). None = the historical (legacy) model: cathode contact
        directly on the n- epi, which produces a non-convergence wall at deep reverse
        bias (a contact artefact).
    N_D_nplus : float
        n+ donor concentration (cm^-3), fully ionized. The n+ inherits the device's
        scalar mobility and SRH lifetimes (like the rest of the model).
    nplus_generation : bool
        Source-side option read by the generation-profile builders: False (default)
        = radiation generation confined to the epi (the deposit beyond the epi end is
        dropped, as without an n+ layer); True = the profile continues into the n+.
        Whether the alpha deposit in the n+ counts is an open physics question.

    Returns
    -------
    device_info : dict
        Dictionary with device_name, region_name, junction_pos, epi_thickness_cm,
        N_D, N_A_ionized, params (SiC4H_Parameters instance), and other metadata.
        ``epi_end`` is the x of the n-/n+ boundary (== ``total_length`` without n+).
    """
    if N_D is None:
        N_D = DEFAULT_N_D

    params = SiC4H_Parameters()
    junction_pos = substrate_thickness_cm
    epi_end = substrate_thickness_cm + epi_thickness_cm
    _check_nplus(nplus_thickness_cm, N_D_nplus)
    nplus = None if nplus_thickness_cm is None else (epi_end, N_D_nplus)
    total_length = epi_end + (nplus_thickness_cm or 0.0)

    # --- Create 1D mesh with non-uniform spacing ---
    build_diode_mesh(
        device_name, region_name, junction_pos, epi_end, nplus_thickness_cm=nplus_thickness_cm
    )
    mesh_name = f"{device_name}_mesh"
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
    # T-dependent intrinsic concentration (calibrated: n_i(300) == 5e-9 exactly)
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

    # Doping- and temperature-dependent mobility, both carriers at the EPI doping N_D (the
    # sensitive volume; holes are the minority carriers collected there). Before 2026-10-17
    # mu_p was evaluated at the p+ doping N_A (1e19: 75.7 vs 121 cm2/Vs at 1.07e15;
    # .issues/audit-v7 A-03). Scalar per region: the p+ and n+ inherit the epi values.
    # For graded profiles N_D is the uniform fallback (DEFAULT_N_D) when not given.
    mu_n = mobility_caughey_thomas_T(N_D, T, "electron", params)
    mu_p = mobility_caughey_thomas_T(N_D, T, "hole", params)
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

    # SRH parameters — T-dependent trap concentrations and lifetimes
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

    if doping_profile == "graded":
        # graded-doping fallbacks: the default stack (DEFAULT_STACK)
        _N_D_junction = N_D_junction if N_D_junction is not None else DEFAULT_STACK["N_D_junction"]
        _N_D_bulk = N_D_bulk if N_D_bulk is not None else DEFAULT_STACK["N_D_bulk"]
        _L_transition = L_transition if L_transition is not None else DEFAULT_STACK["L_transition"]
        set_graded_doping_profile(
            device_name,
            region_name,
            junction_pos,
            N_A_ionized,
            _N_D_junction,
            _N_D_bulk,
            _L_transition,
            nplus=nplus,
        )
        logger.info(
            f"Graded doping: N_D_junction={_N_D_junction:.2e}, "
            f"N_D_bulk={_N_D_bulk:.2e}, L_transition={_L_transition:.2e} cm"
        )
    else:
        set_doping_profile(device_name, region_name, junction_pos, N_A_ionized, N_D, nplus=nplus)

    num_nodes = len(
        devsim.get_node_model_values(device=device_name, region=region_name, name="x")
    )
    logger.info(
        f"Created device '{device_name}': {num_nodes} nodes, "
        f"junction at x={junction_pos*1e4:.1f} um, "
        f"N_A_ionized={N_A_ionized:.2e}, N_D={N_D:.2e}, profile={doping_profile}"
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
        "epi_end": epi_end,
        "nplus_thickness_cm": nplus_thickness_cm,
        "N_D_nplus": N_D_nplus if nplus_thickness_cm is not None else None,
        "nplus_generation": bool(nplus_generation),
    }


def _check_nplus(nplus_thickness_cm, N_D_nplus):
    if nplus_thickness_cm is None:
        return
    if not (nplus_thickness_cm > 0):
        raise ValueError(f"nplus_thickness_cm must be > 0 or None (got {nplus_thickness_cm!r})")
    if not (N_D_nplus > 0):
        raise ValueError(f"N_D_nplus must be > 0 (got {N_D_nplus!r})")


def _nplus_term(nplus, coord="x"):
    """Donor expression term of the n+ layer (empty string without n+)."""
    if nplus is None:
        return ""
    x_end, N_np = nplus
    return f" + {N_np} * step({coord} - {x_end})"


def set_doping_profile(device_name, region_name, junction_pos, N_A_ionized, N_D, nplus=None):
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
    nplus : tuple or None
        (x_epi_end, N_D_nplus): adds N_D_nplus * step(x - x_epi_end). None = no n+.
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
        equation=f"{N_D}*step(x-{junction_pos})" + _nplus_term(nplus),
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="NetDoping",
        equation="Donors-Acceptors",
    )


def set_graded_doping_profile(
    device_name,
    region_name,
    junction_pos,
    N_A_ionized,
    N_D_junction,
    N_D_bulk,
    L_transition,
    nplus=None,
):
    """Set exponentially graded N_D(x) doping profile in the epi layer.

    N_D(x) = N_D_bulk + (N_D_junction - N_D_bulk) * exp(-(x - x_j) / L)
    for x > x_junction (n-side). Acceptors remain a step function in the
    p+ substrate.

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
    N_D_junction : float
        Donor concentration near the junction (cm^-3). Higher value (~1e15).
    N_D_bulk : float
        Donor concentration deep in epi bulk (cm^-3). Lower value (~1e13-1e14).
    L_transition : float
        Characteristic decay length of the doping gradient (cm).
    nplus : tuple or None
        (x_epi_end, N_D_nplus): adds N_D_nplus * step(x - x_epi_end). None = no n+.
    """
    # Acceptors: step function in p+ substrate (same as uniform)
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors",
        equation=f"{N_A_ionized}*step({junction_pos}-x)",
    )

    # Donors: exponential grading in n-side epi layer
    # N_D(x) = N_D_bulk + (N_D_junction - N_D_bulk) * exp(-(x - x_j) / L)
    # Only active for x > junction_pos (n-side)
    donor_expr = (
        f"({N_D_bulk} + ({N_D_junction} - {N_D_bulk}) * "
        f"exp(-max(x - {junction_pos}, 0) / {L_transition})) "
        f"* step(x - {junction_pos})"
    ) + _nplus_term(nplus)
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


def apply_damaged_params(device_info, damaged_params):
    """Override SRH lifetimes and doping profile with radiation-damaged values.

    Must be called AFTER create_sic_device() but BEFORE setup_poisson() and
    solve_equilibrium(). The equilibrium solution must be computed with the
    damaged doping profile, not the pristine one -- otherwise the built-in
    potential and carrier distributions will be wrong.

    Staged device creation sequence:
        1. create_sic_device()       -- mesh + pristine doping
        2. apply_damaged_params()    -- override lifetimes + doping
        3. setup_poisson()           -- Poisson with damaged doping
        4. solve_equilibrium()       -- equilibrium at damaged state
        5. setup_sic_drift_diffusion() -- DD equations

    Parameters
    ----------
    device_info : dict
        Device info dict from create_sic_device (before Poisson setup).
    damaged_params : dict
        Output of compute_damaged_params(). Must contain keys:
        - tau_n : float, damaged electron lifetime (s)
        - tau_p : float, damaged hole lifetime (s)
        - N_D_profile : np.ndarray, damaged donor profile in epi region (cm^-3)
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    junction_pos = device_info["junction_pos"]

    # --- Override SRH lifetimes ---
    devsim.set_parameter(
        device=device, region=region, name="taun", value=damaged_params["tau_n"]
    )
    devsim.set_parameter(
        device=device, region=region, name="taup", value=damaged_params["tau_p"]
    )

    # --- Override doping profile ---
    x_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )

    # Get current Acceptors profile
    acceptors = np.array(
        devsim.get_node_model_values(device=device, region=region, name="Acceptors")
    )

    # Build full Donors array: zeros in p+ substrate, damaged profile in epi
    epi_mask = x_nodes >= junction_pos
    damaged_N_D = damaged_params["N_D_profile"]

    donors = np.zeros_like(x_nodes)
    # damaged_N_D corresponds to epi-only nodes
    n_epi_nodes = int(np.sum(epi_mask))
    if len(damaged_N_D) == n_epi_nodes:
        donors[epi_mask] = damaged_N_D
    elif len(damaged_N_D) == len(x_nodes):
        # Full-length profile: use epi portion only
        donors[epi_mask] = damaged_N_D[epi_mask]
    else:
        raise ValueError(
            f"damaged_params['N_D_profile'] length ({len(damaged_N_D)}) does not match "
            f"epi nodes ({n_epi_nodes}) or total nodes ({len(x_nodes)})"
        )

    # Set Donors via set_node_values (overrides the model-computed values)
    devsim.set_node_values(
        device=device, region=region, name="Donors", values=list(donors)
    )

    # Recompute NetDoping = Donors - Acceptors
    net_doping = donors - acceptors
    devsim.set_node_values(
        device=device, region=region, name="NetDoping", values=list(net_doping)
    )

    # epi_mask (x >= junction_pos) also covers the n+ layer, consistently with the
    # reference profiles the callers extract, so the n+ doping is carried through; the
    # logged epi mean must exclude it (1e19 would dominate).
    epi_only = epi_mask & (x_nodes <= device_info.get("epi_end", np.inf) + 1e-12)
    logger.info(
        f"Applied damaged params: tau_n={damaged_params['tau_n']:.3e}s, "
        f"tau_p={damaged_params['tau_p']:.3e}s, "
        f"N_D_epi_mean={np.mean(donors[epi_only]):.3e} cm^-3"
    )


def calibrate_graded_doping(
    target_W_data=None,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=1e-4,
    N_A=1e19,
    T=300,
    x0=None,
    maxiter=80,
):
    """Calibrate graded doping parameters to match experimental W(V) data.

    Fits {N_D_junction, N_D_bulk, L_transition} by minimizing the sum of
    squared relative errors between simulated and experimental depletion
    widths at multiple bias voltages.

    Uses the full drift-diffusion solver for proper carrier transport under
    bias. The target voltages use reverse-bias convention (negative V means
    reverse bias on the diode); internally these are applied as positive
    cathode bias in devsim.

    Parameters
    ----------
    target_W_data : dict or None
        Mapping of voltage (V) to target depletion width (cm).
        Negative voltages = reverse bias on diode.
        Default: {0.0: 1.7e-4, -10.0: 9.5e-4, -30.0: 9.73e-4}.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm).
    substrate_thickness_cm : float
        p+ substrate thickness (cm).
    N_A : float
        Total acceptor concentration (cm^-3).
    T : float
        Temperature (K).
    x0 : array_like or None
        Initial guess [N_D_junction, N_D_bulk, L_transition].
        Default: [2.9e15, 8.5e13, 1e-4].
    maxiter : int
        Maximum optimizer iterations.

    Returns
    -------
    result : dict
        Dictionary with keys: N_D_junction, N_D_bulk, L_transition,
        final_cost, success, W_simulated (dict of V: W pairs).
    """
    from etna.core.drift_diffusion import create_dd_device
    from etna.core.poisson import extract_depletion_width_numerical
    import devsim.python_packages.simple_physics as simple_physics

    if target_W_data is None:
        target_W_data = {0.0: 1.7e-4, -10.0: 9.5e-4, -30.0: 9.73e-4}

    if x0 is None:
        x0 = [2.9e15, 8.5e13, 1e-4]

    voltages = sorted(target_W_data.keys(), reverse=True)  # 0, -10, -30
    W_exp = np.array([target_W_data[v] for v in voltages])
    # Convert reverse bias convention to cathode voltage (positive = reverse)
    cathode_voltages = [-v for v in voltages]  # 0, 10, 30

    import uuid as _uuid

    _run_id = _uuid.uuid4().hex[:8]
    _trial_counter = [0]

    def objective(params_vec):
        """Compute cost function for a trial doping parameter set."""
        N_D_j, N_D_b, L_t = params_vec

        # Enforce bounds via penalty
        if (
            N_D_j < 1e14
            or N_D_j > 1e16
            or N_D_b < 1e12
            or N_D_b > 1e15
            or L_t < 0.5e-4
            or L_t > 5e-4
            or N_D_b >= N_D_j
        ):
            return 1e6

        trial_name = f"cal_{_run_id}_{_trial_counter[0]}"
        mesh_name = f"{trial_name}_mesh"
        _trial_counter[0] += 1

        try:
            device_info = create_dd_device(
                device_name=trial_name,
                epi_thickness_cm=epi_thickness_cm,
                substrate_thickness_cm=substrate_thickness_cm,
                N_A=N_A,
                T=T,
                doping_profile="graded",
                N_D_junction=N_D_j,
                N_D_bulk=N_D_b,
                L_transition=L_t,
            )

            device = device_info["device_name"]
            bias_name = simple_physics.GetContactBiasName("cathode")

            W_sim = []
            current_V = 0.0
            for v_cath in cathode_voltages:
                if abs(v_cath) < 1e-12:
                    W = extract_depletion_width_numerical(device_info)
                else:
                    # Ramp cathode in 0.5V steps
                    V = current_V + 0.5
                    while V <= v_cath + 1e-10:
                        devsim.set_parameter(device=device, name=bias_name, value=V)
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
                        V += 0.5
                    current_V = v_cath
                    W = extract_depletion_width_numerical(device_info)
                W_sim.append(W)

            W_sim = np.array(W_sim)
            cost = np.sum(((W_sim - W_exp) / W_exp) ** 2)

        except Exception as e:
            logger.warning(f"Trial {trial_name} failed: {e}")
            cost = 1e6
        finally:
            try:
                devsim.delete_device(device=trial_name)
            except Exception:
                pass
            try:
                devsim.delete_mesh(mesh=mesh_name)
            except Exception:
                pass

        return cost

    result = scipy.optimize.minimize(
        objective,
        x0,
        method="Nelder-Mead",
        options={"maxiter": maxiter, "xatol": 1e-10, "fatol": 1e-6},
    )

    N_D_j_opt, N_D_b_opt, L_t_opt = result.x

    logger.info(
        f"Calibration complete: N_D_junction={N_D_j_opt:.3e}, "
        f"N_D_bulk={N_D_b_opt:.3e}, L_transition={L_t_opt:.3e} cm, "
        f"cost={result.fun:.6f}, success={result.success}"
    )

    return {
        "N_D_junction": N_D_j_opt,
        "N_D_bulk": N_D_b_opt,
        "L_transition": L_t_opt,
        "final_cost": result.fun,
        "success": result.success,
    }
