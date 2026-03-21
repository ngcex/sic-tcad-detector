"""devsim 1D device setup for 4H-SiC p+/n- diode.

Creates a 1D mesh with non-uniform spacing, applies 4H-SiC material parameters,
and sets up doping profiles for the p+ substrate / n- epitaxial layer structure.
Supports both uniform and graded (exponential) epitaxial doping profiles.

All units CGS (cm, cm^-3, F/cm, eV, s) per devsim convention.

References:
    - devsim documentation: https://devsim.net/
    - devsim diode example pattern
    - Petringa et al. experimental device structure
"""

import logging

import devsim
import numpy as np
import scipy.optimize

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
    doping_profile="uniform",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
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

    if doping_profile == "graded":
        # Defaults for graded profile if not provided
        _N_D_junction = N_D_junction if N_D_junction is not None else 1e15
        _N_D_bulk = N_D_bulk if N_D_bulk is not None else 5e13
        _L_transition = L_transition if L_transition is not None else 2e-4
        set_graded_doping_profile(
            device_name,
            region_name,
            junction_pos,
            N_A_ionized,
            _N_D_junction,
            _N_D_bulk,
            _L_transition,
        )
        logger.info(
            f"Graded doping: N_D_junction={_N_D_junction:.2e}, "
            f"N_D_bulk={_N_D_bulk:.2e}, L_transition={_L_transition:.2e} cm"
        )
    else:
        set_doping_profile(device_name, region_name, junction_pos, N_A_ionized, N_D)

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
        "params": params,
        "mu_n": mu_n,
        "mu_p": mu_p,
        "num_nodes": num_nodes,
        "doping_profile": doping_profile,
        "N_D_junction": N_D_junction,
        "N_D_bulk": N_D_bulk,
        "L_transition": L_transition,
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


def set_graded_doping_profile(
    device_name,
    region_name,
    junction_pos,
    N_A_ionized,
    N_D_junction,
    N_D_bulk,
    L_transition,
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

    Parameters
    ----------
    target_W_data : dict or None
        Mapping of voltage (V) to target depletion width (cm).
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
        Default: [5e14, 5e13, 2e-4].
    maxiter : int
        Maximum optimizer iterations.

    Returns
    -------
    result : dict
        Dictionary with keys: N_D_junction, N_D_bulk, L_transition,
        final_cost, success, W_simulated (dict of V: W pairs).
    """
    from src.poisson import (
        setup_poisson,
        solve_equilibrium,
        ramp_voltage,
        extract_depletion_width_numerical,
    )

    if target_W_data is None:
        target_W_data = {0.0: 1.7e-4, -10.0: 9.5e-4, -30.0: 9.73e-4}

    if x0 is None:
        x0 = [5e14, 5e13, 2e-4]

    voltages = sorted(target_W_data.keys(), reverse=True)  # 0, -10, -30
    W_exp = np.array([target_W_data[v] for v in voltages])

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

        trial_name = f"cal_trial_{_trial_counter[0]}"
        _trial_counter[0] += 1

        try:
            device_info = create_sic_device(
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

            setup_poisson(device_info)
            solve_equilibrium(device_info)

            W_sim = []
            for v in voltages:
                if v == 0.0:
                    W = extract_depletion_width_numerical(device_info)
                else:
                    ramp_voltage(
                        device_info,
                        contact_name="cathode",
                        V_start=0.0,
                        V_end=v,
                        V_step=-0.5 if v < 0 else 0.5,
                    )
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
