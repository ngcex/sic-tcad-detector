"""Capacitance-voltage (C-V) analysis for 4H-SiC p+/n- diode.

Computes junction capacitance from depletion width using the parallel-plate
approximation C = eps * A / W, and provides a C-V sweep function that
uses the numerical Poisson/DD solver to extract W(V) at each bias point.

The 1/C^2 vs V (Mott-Schottky) analysis enables extraction of doping
concentration and built-in potential from the linear fit slope/intercept.

All units CGS (cm, F/cm^2, V) per devsim convention unless otherwise noted.

References:
    - Reference: C-V measurements at 1 kHz, W(0V)=1.7um, W(-10V)=9.5um
    - Sze & Ng, "Physics of Semiconductor Devices", Ch. 2: Junction capacitance
"""

from etna.core.solver_tol import DC_REL_ERROR
import logging
import uuid

import matplotlib.pyplot as plt
import numpy as np

import devsim
import devsim.python_packages.simple_physics as simple_physics

from etna.core.device import apply_damaged_params, create_sic_device
from etna.core.drift_diffusion import setup_sic_drift_diffusion
from etna.core.poisson import (
    extract_depletion_width_numerical,
    setup_poisson,
    solve_equilibrium,
)
from etna.core.radiation_damage import compute_damaged_params, compute_phi_crit
from etna.core.sic_material import srh_lifetime
from etna.core.device import DEFAULT_STACK as _STACK

logger = logging.getLogger(__name__)

# Physical constants (CGS)
EPS_0 = 8.854e-14  # F/cm (vacuum permittivity)
Q = 1.602e-19  # C (elementary charge)


def junction_capacitance(W, eps_r=9.7, area=1.0):
    """Compute junction capacitance from depletion width.

    Uses the parallel-plate approximation: C = eps_r * eps_0 * area / W.

    Parameters
    ----------
    W : float or array_like
        Depletion width (cm).
    eps_r : float
        Relative permittivity of 4H-SiC. Default 9.7.
    area : float
        Junction area (cm^2). Default 1.0 (gives F/cm^2).

    Returns
    -------
    C : float or ndarray
        Capacitance (F if area given, F/cm^2 if area=1.0).
    """
    W = np.asarray(W, dtype=float)
    return eps_r * EPS_0 * area / W


def depletion_width_from_capacitance(C, eps_r=9.7, area=1.0):
    """Compute depletion width from capacitance (inverse of junction_capacitance).

    Parameters
    ----------
    C : float or array_like
        Capacitance (F or F/cm^2).
    eps_r : float
        Relative permittivity of 4H-SiC. Default 9.7.
    area : float
        Junction area (cm^2). Default 1.0.

    Returns
    -------
    W : float or ndarray
        Depletion width (cm).
    """
    C = np.asarray(C, dtype=float)
    return eps_r * EPS_0 * area / C


def compute_cv_from_depletion(voltages, depletion_widths, eps_r=9.7, area=1.0):
    """Compute C-V curve from voltage and depletion width arrays.

    Parameters
    ----------
    voltages : array_like
        Applied bias voltages (V).
    depletion_widths : array_like
        Depletion widths at each voltage (cm).
    eps_r : float
        Relative permittivity. Default 9.7.
    area : float
        Junction area (cm^2). Default 1.0.

    Returns
    -------
    result : dict
        Dictionary with:
        - "voltages": numpy array of voltages (V)
        - "capacitance": numpy array of C values (F or F/cm^2)
        - "one_over_C_squared": numpy array of 1/C^2 values
    """
    voltages = np.asarray(voltages, dtype=float)
    W = np.asarray(depletion_widths, dtype=float)

    C = junction_capacitance(W, eps_r=eps_r, area=area)

    return {
        "voltages": voltages,
        "capacitance": C,
        "one_over_C_squared": 1.0 / C**2,
    }


def cv_sweep(device_info, V_range, eps_r=9.7, area=1.0):
    """Sweep reverse bias and compute C-V from numerical depletion widths.

    At each voltage, ramps bias using the DD/Poisson solver, extracts the
    depletion width numerically from the carrier profile, and computes
    capacitance via the parallel-plate approximation.

    Reverse bias convention: negative V_range values mean reverse bias on
    the diode. Internally, these are applied as positive cathode bias in
    devsim (positive V on cathode = reverse bias for p+/n- diode).

    Parameters
    ----------
    device_info : dict
        Device info dict (DD equations must be set up via create_dd_device).
    V_range : array_like
        Array of reverse bias voltages (V, should be <= 0 for reverse bias).
    eps_r : float
        Relative permittivity. Default 9.7.
    area : float
        Junction area (cm^2). Default 1.0.

    Returns
    -------
    result : dict
        Dictionary with:
        - "voltages": numpy array of voltages actually reached (V), in
          conventional form — shorter than the requested V_range if the DD
          solver hit a convergence wall (typically full depletion /
          punch-through) partway through the sweep; see "truncated".
        - "depletion_widths": numpy array of W values (cm), same length as
          "voltages"
        - "capacitance": numpy array of C values (F or F/cm^2), same length
          as "voltages"
        - "truncated": True if the solver failed to converge before reaching
          the end of V_range (a partial result is returned instead of
          raising); False if the full requested sweep converged
        - "requested_v_stop": the most extreme (deepest reverse bias) voltage
          in V_range, for comparison against voltages.min() when "truncated"
          is True
    """
    V_range = np.asarray(V_range, dtype=float)
    device = device_info["device_name"]
    bias_name = simple_physics.GetContactBiasName("cathode")

    depletion_widths = []
    solved_voltages = []
    truncated = False

    # Extract W at 0V first (current equilibrium state)
    W0 = extract_depletion_width_numerical(device_info)
    current_V_cathode = 0.0

    for V_target in V_range:
        if abs(V_target) < 1e-12:
            # Use equilibrium W
            depletion_widths.append(W0)
            solved_voltages.append(0.0)
            continue

        # Convert conventional reverse bias to cathode voltage
        # Negative V_target (reverse bias) -> positive cathode voltage
        V_cathode_target = -V_target

        # Ramp cathode from the current state in equal steps of at most 0.5 V
        # that end exactly on the target (an off-lattice target used to stop at
        # the last multiple of 0.5 V below it while being labelled V_target).
        n_steps = max(
            1, int(np.ceil(abs(V_cathode_target - current_V_cathode) / 0.5 - 1e-9))
        )
        V_steps = np.linspace(current_V_cathode, V_cathode_target, n_steps + 1)[1:]
        converged = True

        for V in V_steps:
            devsim.set_parameter(device=device, name=bias_name, value=V)
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
                except devsim.error as e:
                    logger.warning(f"cv_sweep: failed at V_cathode={V:.2f}V: {e}")
                    converged = False
                    break

        if converged:
            W = extract_depletion_width_numerical(device_info)
            depletion_widths.append(W)
            solved_voltages.append(V_target)
            current_V_cathode = V_cathode_target
        else:
            # Typically full depletion / punch-through: the depletion edge
            # reaches the back contact and the DD solver's Newton iteration
            # stops converging. This is a physical sweep endpoint, not a
            # bug — stop here and return the converged portion.
            logger.warning(
                f"cv_sweep: stopped sweep at V={V_target:.2f}V "
                f"(requested down to {V_range.min():.2f}V)"
            )
            truncated = True
            break

    voltages = np.array(solved_voltages)
    W_arr = np.array(depletion_widths)
    C_arr = junction_capacitance(W_arr, eps_r=eps_r, area=area)

    return {
        "voltages": voltages,
        "depletion_widths": W_arr,
        "capacitance": C_arr,
        "truncated": truncated,
        "requested_v_stop": float(V_range.min()) if len(V_range) else 0.0,
    }


# Batched-engine scope of cv_at_fluence (see charge_collection._fluence_engine_reason):
# graded/uniform epi, the T box of etna.api.simulation._ENGINE_BOX, reverse biases.
_CV_ENGINE_T = (250.0, 400.0)
# Largest cathode step between two engine solves (V); finer grids are solved as given.
_CV_ENGINE_STEP_V = 10.0


def _cv_engine_reason(dev_kw, V_range):
    """Why the batched engine cannot serve this C-V sweep (None = it can)."""
    profile = dev_kw.get("doping_profile", "graded")
    if profile not in ("graded", "uniform"):
        return f"doping_profile={profile!r}"
    T = dev_kw.get("T", 300.0)
    if not _CV_ENGINE_T[0] <= T <= _CV_ENGINE_T[1]:
        return f"T={T:g} outside [{_CV_ENGINE_T[0]:g}, {_CV_ENGINE_T[1]:g}]"
    if V_range.size == 0 or not np.all(np.isfinite(V_range)):
        return "empty or non-finite bias grid"
    if V_range.max() > 0:
        return "forward bias requested"
    from etna.core.charge_collection import _legacy_depth_reason

    return _legacy_depth_reason(dev_kw, float(-V_range.min()))


def _cv_engine(x, donors, acceptors, info, V_range, area):
    """C-V of one (damaged) device with the batched engine; None if not converged.

    Electrostatics as the devsim path (reverse-bias currents do not move psi):
    quasi-Fermi Poisson on the devsim mesh with the given donors/acceptors, walked
    up in |V| (steps <= _CV_ENGINE_STEP_V), W extracted as
    extract_depletion_width_numerical (n/Donors >= 0.5) and C as cv_sweep.
    """
    from etna.core.batched.extract import depletion_width
    from etna.core.batched.material import thermal_voltage
    from etna.core.batched.poisson import PoissonBatch
    from etna.core.batched.sweeps import _Ramp

    pb = PoissonBatch(
        x, donors[None, :], acceptors[None, :], info["n_i"], thermal_voltage(info["T"])
    )
    W_of = lambda n: float(
        depletion_width(
            x, n, pb.donors, info["junction_pos"], info["epi_thickness_cm"]
        )[0]
    )
    ramp = _Ramp(pb, 1e-8)
    targets = np.unique(-V_range)  # cathode voltages >= 0, increasing
    W_at = {}
    if targets[0] == 0.0:
        W_at[0.0] = W_of(pb.solve(ramp.psi, tol=1e-8).n)
    Vc_prev = 0.0
    for Vc in targets[targets > 0]:
        n_sub = max(1, int(np.ceil((Vc - Vc_prev) / _CV_ENGINE_STEP_V - 1e-9)))
        for V in np.linspace(Vc_prev, Vc, n_sub + 1)[1:]:
            res = ramp.to(V)
        W_at[Vc] = W_of(res.n)
        Vc_prev = Vc
    if not ramp.converged.all():
        return None
    W = np.array([W_at[-v + 0.0] for v in V_range])
    return {
        "voltages": V_range.copy(),
        "depletion_widths": W,
        "capacitance": junction_capacitance(W, area=area),
        "truncated": False,
        "requested_v_stop": float(V_range.min()),
        "solver": "batched",
    }


def cv_at_fluence(
    fluence,
    V_range,
    area=1.0,
    epi_thickness_cm=10e-4,
    N_D_junction=_STACK["N_D_junction"],
    N_D_bulk=_STACK["N_D_bulk"],
    L_transition=_STACK["L_transition"],
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
    phi_crit_threshold=0.90,
    device_kwargs=None,
    solver="auto",
):
    """Compute C-V curve for a device at a given proton fluence.

    The damaged donor profile is applied before the Poisson equilibrium. By
    default the sweep is solved by the batched numba engine (same mesh, doping,
    contacts and W extraction as the devsim path; W agrees with devsim to ~1e-14,
    tests/test_fluence_engine.py). Out of the engine scope, or if it does not
    converge, a fresh damaged devsim DD device runs ``cv_sweep`` (staged device
    creation) and is cleaned up.

    Parameters
    ----------
    fluence : float
        Proton fluence (protons/cm^2).
    V_range : array_like
        Array of reverse bias voltages (V, should be <= 0 for reverse bias).
    area : float
        Junction area (cm^2). Default 1.0.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default 10 um.
    N_D_junction : float
        Junction-side donor concentration (cm^-3). Default: DEFAULT_STACK (2.974271e15).
    N_D_bulk : float
        Bulk-side donor concentration (cm^-3). Default: DEFAULT_STACK (6.335936e13).
    L_transition : float
        Doping transition length (cm). Default: DEFAULT_STACK (0.987 um).
    energy_MeV : float
        Proton energy (MeV). Default 62.0.
    lifetime_model : str
        "linear" or "logarithmic". Default "linear".
    damage_params : RadiationDamageParams or None
        Custom damage parameters. Default: RadiationDamageParams().
    phi_crit_threshold : float
        Fraction of Phi_crit at which to warn. Default 0.90.
    device_kwargs : dict or None
        Extra create_sic_device keywords (geometry: substrate_thickness_cm,
        nplus_thickness_cm, N_D_nplus, ...; or doping), overriding the defaults.
        Default None = the default stack (DEFAULT_STACK); pass
        ``etna.api.device.device_kwargs(legacy_config())`` minus epi for the legacy model.
    solver : {"auto", "engine", "devsim"}
        "auto" (default): batched engine in its scope (graded/uniform epi,
        250-400 K, reverse biases), devsim otherwise or if it does not converge.
        "engine": raise if out of scope. "devsim": the devsim path only.

    Returns
    -------
    dict or None
        C-V result dict with keys: voltages, depletion_widths, capacitance,
        truncated, requested_v_stop, fluence, solver ("batched" or "devsim").
        Returns None if fluence >= Phi_crit (solver would diverge).
    """
    if solver not in ("auto", "engine", "devsim"):
        raise ValueError(
            f"solver must be 'auto', 'engine' or 'devsim' (got {solver!r})"
        )
    V_range = np.asarray(V_range, dtype=float)
    dev_kw = dict(
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=N_D_junction,
        N_D_bulk=N_D_bulk,
        L_transition=L_transition,
    )
    dev_kw.update(device_kwargs or {})

    # Extract pristine lifetimes
    pristine_tau_n = srh_lifetime(300.0, "electron")
    pristine_tau_p = srh_lifetime(300.0, "hole")

    # Create a reference device to extract pristine N_D profile
    ref_id = uuid.uuid4().hex[:8]
    ref_name = f"cv_ref_{ref_id}"
    ref_info = create_sic_device(
        device_name=ref_name,
        **dev_kw,
    )
    try:
        ref_x = np.array(
            devsim.get_node_model_values(
                device=ref_name, region=ref_info["region_name"], name="x"
            )
        )
        ref_donors = np.array(
            devsim.get_node_model_values(
                device=ref_name, region=ref_info["region_name"], name="Donors"
            )
        )
        ref_acceptors = np.array(
            devsim.get_node_model_values(
                device=ref_name, region=ref_info["region_name"], name="Acceptors"
            )
        )
        epi_mask = ref_x >= ref_info["junction_pos"]
        pristine_N_D_profile = ref_donors[epi_mask]
    finally:
        try:
            devsim.delete_device(device=ref_name)
        except Exception:
            pass

    # Check Phi_crit
    phi_crit_info = compute_phi_crit(
        pristine_N_D_profile,
        energy_MeV=energy_MeV,
    )
    phi_crit = phi_crit_info["phi_crit_proton"]

    if fluence >= phi_crit:
        logger.error(
            "Fluence %.3e >= Phi_crit %.3e protons/cm^2: full compensation, "
            "cannot solve. Returning None.",
            fluence,
            phi_crit,
        )
        return None

    if fluence >= phi_crit_threshold * phi_crit:
        logger.warning(
            "Fluence %.3e is >= %.0f%% of Phi_crit (%.3e): approaching "
            "full carrier compensation.",
            fluence,
            phi_crit_threshold * 100,
            phi_crit,
        )

    # Create device for this fluence point
    dev_id = uuid.uuid4().hex[:8]
    dev_name = f"cv_fluence_{dev_id}"

    try:
        # Compute damaged parameters
        damaged = compute_damaged_params(
            pristine_tau_n=pristine_tau_n,
            pristine_tau_p=pristine_tau_p,
            N_D_profile=pristine_N_D_profile,
            fluence=fluence,
            energy_MeV=energy_MeV,
            damage_params=damage_params,
            lifetime_model=lifetime_model,
        )

        if solver != "devsim":
            reason = _cv_engine_reason(dev_kw, V_range)
            if reason is None:
                donors = np.zeros_like(ref_x)
                donors[epi_mask] = damaged["N_D_profile"]
                cv_result = _cv_engine(
                    ref_x, donors, ref_acceptors, ref_info, V_range, area
                )
                if cv_result is not None:
                    cv_result["fluence"] = fluence
                    return cv_result
                logger.warning(
                    "cv_at_fluence: batched engine did not converge, using devsim"
                )
            elif solver == "engine":
                raise ValueError(
                    f"cv_at_fluence: batched engine unavailable ({reason})"
                )
            else:
                logger.info("cv_at_fluence: using devsim (%s)", reason)

        # Staged device creation
        device_info = create_sic_device(
            device_name=dev_name,
            **dev_kw,
        )

        # Apply damage BEFORE Poisson setup
        if fluence > 0:
            apply_damaged_params(device_info, damaged)

        # Continue staged setup
        setup_poisson(device_info)
        solve_equilibrium(device_info)
        setup_sic_drift_diffusion(device_info)

        # Run C-V sweep
        cv_result = cv_sweep(device_info, V_range, area=area)
        cv_result["fluence"] = fluence
        cv_result["solver"] = "devsim"

        logger.info(
            "cv_at_fluence: fluence=%.3e, %d voltage points solved",
            fluence,
            len(cv_result["voltages"]),
        )

        return cv_result

    except Exception as e:
        logger.error("cv_at_fluence failed at fluence=%.3e: %s", fluence, e)
        raise
    finally:
        try:
            devsim.delete_device(device=dev_name)
        except Exception:
            pass


def plot_cv_evolution(cv_results, fluences, ax=None, title=None):
    """Overlay C-V curves at different fluence levels.

    Uses a viridis colormap gradient to distinguish fluence levels, with
    lower fluence in darker colors and higher fluence in lighter colors.

    Parameters
    ----------
    cv_results : list of dict
        List of C-V result dicts (from cv_at_fluence or cv_sweep).
        None entries are skipped (above Phi_crit).
    fluences : array_like
        Fluence values corresponding to each cv_result (protons/cm^2).
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    title : str, optional
        Plot title. Default: "C-V Evolution with Fluence".

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    fluences = np.asarray(fluences, dtype=float)
    cmap = plt.cm.viridis

    # Normalize fluence range for colormap
    valid_fluences = [f for f, r in zip(fluences, cv_results) if r is not None]
    if len(valid_fluences) == 0:
        logger.warning("plot_cv_evolution: no valid C-V results to plot")
        return ax

    f_min = min(valid_fluences)
    f_max = max(valid_fluences)
    if f_max > f_min:
        norm = plt.Normalize(vmin=f_min, vmax=f_max)
    else:
        norm = plt.Normalize(vmin=0, vmax=max(f_max, 1.0))

    for fluence, cv_result in zip(fluences, cv_results):
        if cv_result is None:
            continue
        color = cmap(norm(fluence))
        if fluence == 0:
            label = "Pristine"
        else:
            label = f"{fluence:.1e} p/cm$^2$"
        ax.plot(
            cv_result["voltages"], cv_result["capacitance"], color=color, label=label
        )

    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("Capacitance (F/cm$^2$)")
    ax.set_title(title or "C-V Evolution with Fluence")
    ax.legend(fontsize=8)

    return ax
