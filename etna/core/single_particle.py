"""Single-particle transient charge collection for 4H-SiC microdosimeter.

Provides:
- Ion track generation profile (2D Gaussian along vertical track)
- Single-particle transient simulation (instantaneous injection + BDF1 collection)
- Current pulse analysis (peak, collection time, integrated charge)
- CCE(LET) lookup table construction, save/load with interpolation

The single-particle module bridges TCAD simulation with Monte Carlo particle
transport.  For each ion event, the LET determines the charge deposited along
the track.  A transient drift-diffusion simulation computes the fraction of
charge collected at the contacts (CCE).  By sweeping LET values, a CCE(LET)
lookup table is built once per geometry/bias, then applied to thousands of
MC events without re-running TCAD.

Units: CGS (cm, cm^-3, s, A/cm, C/cm) per devsim convention.
2D current is A/cm (per unit z-depth), charge is C/cm.

References:
    - Phase 20: 2D CCE computation, device creation, mesh integration
    - Phase 21 research: ion track physics, instantaneous injection method
    - TransientSolver pattern (transient.py): BDF1 solve API reference
"""

import json
import logging

import devsim
import numpy as np
import pandas as pd

from etna.core.charge_collection_2d import (
    create_2d_dd_device,
    injected_generation,
    integrate_nodes_box_2d,
)
from etna.core.charge_collection import add_generation_to_dd
from etna.core.drift_diffusion import extract_contact_current
from etna.core.solver_tol import DC_REL_ERROR

logger = logging.getLogger(__name__)

# Physical constants
Q_ELECTRON = 1.602e-19  # C, elementary charge

# Headroom over the physical CCE ceiling (1.0) before the table builder warns. The
# backward-Euler collected charge (collected_charge) has no quadrature error on top
# of the time stepping, and the validated runs give CCE_raw <= 1 (default stack
# 0.99997 at 50 V); an excess beyond this tolerance signals a bug, not numerics.
# (The old 5 % headroom covered the post-injection trapezoid, audit C4.)
_CCE_EXCESS_TOL = 1e-3


def ion_track_generation_2d(
    device_info, LET_keV_um, x_ion_cm=0.0, track_sigma_cm=1e-4, E_pair_eV=7.8
):
    """Create 2D generation profile for a single ion track.

    Converts LET (linear energy transfer) to a spatially-resolved charge
    generation profile on the 2D mesh.  The ion track is vertical (along y)
    through the epi layer, with a Gaussian lateral profile centered at
    x_ion_cm.

    Parameters
    ----------
    device_info : dict
        2D device info from create_2d_dd_device.
    LET_keV_um : float
        Linear energy transfer (keV/um).
    x_ion_cm : float
        Lateral position of ion track (cm). Default: 0.0 (center).
    track_sigma_cm : float
        Gaussian radial width of track (cm). Default: 1e-4 (1 um).
    E_pair_eV : float
        Electron-hole pair creation energy (eV). Default: 7.8 (4H-SiC).

    Returns
    -------
    generation : ndarray
        Generation density at each mesh node (pairs/cm^3).
    Q_generated_C_per_cm : float
        Total generated charge per unit z-depth (C/cm): the box-rule integral of the
        injected generation (contact nodes zeroed) times q_electron, the convention of
        compute_cce_2d.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    junction_pos = device_info["substrate_thickness_cm"]
    epi_thickness_cm = device_info["epi_thickness_cm"]

    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    y = np.array(devsim.get_node_model_values(device=device, region=region, name="y"))

    # LET to linear pair density: pairs/cm
    # LET (keV/um) -> eV/um = LET*1e3 -> eV/cm = LET*1e3*1e4
    # pairs/cm = (LET * 1e3 * 1e4) / E_pair_eV = LET * 1e7 / E_pair_eV
    pairs_per_cm = LET_keV_um * 1e3 / E_pair_eV * 1e4

    # Gaussian lateral profile, normalized so integral over x = 1/cm
    lateral = np.exp(-0.5 * ((x - x_ion_cm) / track_sigma_cm) ** 2)
    lateral /= track_sigma_cm * np.sqrt(2 * np.pi)

    # Track exists only in epi region
    in_epi = (y >= junction_pos) & (y <= junction_pos + epi_thickness_cm)

    # Generation: pairs/cm^3
    generation = np.where(in_epi, pairs_per_cm * lateral, 0.0)
    # Contact nodes carry no injected generation (Dirichlet rows): zero them so that
    # Q_generated counts only the pairs the solver injects, the convention of
    # compute_cce_2d (since 2026-10-18). Matters on the legacy diode, whose cathode row
    # is the epi end (cathode share of the track ~1e-3); with the n+ layer the track
    # touches no contact (anode in the p+, cathode at the end of the n+).
    generation = injected_generation(device_info, generation)

    # Q_generated = integrate(generation * dA) * q  [C/cm for 2D], with the solver's box
    # rule (the quadrature that injects the source term; compute_cce_2d's I_generated).
    # Was the triangle rule (integrate_over_mesh_2d) before 2026-10-18.
    Q_generated_C_per_cm = integrate_nodes_box_2d(device_info, generation) * Q_ELECTRON

    logger.info(
        f"Ion track: LET={LET_keV_um:.1f} keV/um, x_ion={x_ion_cm*1e4:.1f} um, "
        f"pairs/cm={pairs_per_cm:.2e}, Q_gen={Q_generated_C_per_cm:.3e} C/cm"
    )

    return generation, Q_generated_C_per_cm


def simulate_single_particle(
    device_info,
    generation_profile,
    dt_inject=1e-12,
    t_collect=200e-9,
    dt_min=1e-13,
    dt_max=1e-9,
    contact=None,
    rtol_converged=1e-7,
    n_converged=3,
):
    """Simulate single-particle charge collection transient.

    Uses the generation-pulse method (Approach A from research):
    1. Initialize transient state with zero generation
    2. Inject charge as a short pulse (generation/dt_inject for one timestep)
    3. Zero generation and run BDF1 time-stepping until the collected charge
       has converged

    Q_collected is the backward-Euler sum over every step, the injection step
    included (collected_charge); the run stops when the charge added by a step is
    below ``rtol_converged`` of the running total for ``n_converged`` consecutive
    steps, or at ``t_collect`` (then ``converged`` is False and a warning is logged).

    Parameters
    ----------
    device_info : dict
        2D device with DD and bias applied.
    generation_profile : ndarray
        Ion track charge density (pairs/cm^3) from ion_track_generation_2d.
    dt_inject : float
        Injection pulse duration (s). Default: 1e-12 (1 ps).
    t_collect : float
        Maximum collection time to simulate (s). Default: 200e-9 (200 ns).
    dt_min : float
        Minimum time step during collection (s). Default: 1e-13.
    dt_max : float
        Maximum time step during collection (s). Default: 1e-9.
    contact : {"cathode", "anode"} or None
        Contact whose (direct, conduction) current is integrated. None = "anode"
        when the device has an n+ layer (``device_info["nplus_thickness_cm"]`` not
        None), else "cathode" (the historical choice; legacy numbers unchanged).
        With an n+ cathode the direct cathode current is the majority-carrier current
        in 1e19 material, a small difference of huge Scharfetter-Gummel terms
        dominated by round-off (same finding as transient_cce_vs_dose_rate). The
        balanced steady-state form (extract_steady_state_current) is NOT used: it has
        no dQ/dt term and is wrong in a transient.
    rtol_converged : float
        Relative charge increment below which a step counts as converged.
    n_converged : int
        Consecutive converged steps required to stop.

    Returns
    -------
    result : dict
        - "times": end time of each step from the injection start (s), ndarray;
          the first entry is the injection step (t = dt_inject)
        - "currents": contact current at each time (A/cm for 2D), ndarray
        - "Q_collected": collected charge (C/cm for 2D), float
        - "I_dark": baseline dark current (A/cm), float
        - "I_peak": peak signal current (A/cm), float
        - "current_contact": the contact actually integrated, str
        - "converged": False if t_collect / max_steps / a failed step ended the run
    """
    if contact is None:
        contact = (
            "anode" if device_info.get("nplus_thickness_cm") is not None else "cathode"
        )
    if contact not in ("cathode", "anode"):
        raise ValueError(f"contact must be 'cathode' or 'anode', got {contact!r}")
    n_nodes = len(generation_profile)

    # Step 1: Initialize transient state with zero generation
    zeros = np.zeros(n_nodes)
    add_generation_to_dd(device_info, zeros)
    try:
        # the shared DC tolerance (was 1e-10, which on 2D devices at high bias
        # stagnates near 1e-9 and can only reach the relaxed retry; solver_tol)
        devsim.solve(
            type="transient_dc",
            absolute_error=1e10,
            relative_error=DC_REL_ERROR,
            maximum_iterations=30,
        )
    except devsim.error:
        # 2D devices may need relaxed tolerances for transient_dc init
        logger.info("transient_dc failed at strict tolerances, retrying relaxed")
        devsim.solve(
            type="transient_dc",
            absolute_error=1e12,
            relative_error=1e-8,
            maximum_iterations=100,
        )
    logger.debug("Transient state initialized (transient_dc with zero gen)")

    # Step 2: Record dark current
    I_dark = extract_contact_current(device_info, contact)
    logger.debug(f"Dark current: {I_dark:.4e} A/cm")

    # Step 3: Inject charge as short pulse
    G_rate = generation_profile / dt_inject  # pairs/cm^3/s
    add_generation_to_dd(device_info, G_rate)

    try:
        devsim.solve(
            type="transient_bdf1",
            tdelta=dt_inject,
            absolute_error=1e10,
            relative_error=1e-8,
            maximum_iterations=100,
            charge_error=1e10,
        )
    except devsim.error:
        logger.warning("Injection step failed at standard tolerances, retrying relaxed")
        devsim.solve(
            type="transient_bdf1",
            tdelta=dt_inject,
            absolute_error=1e12,
            relative_error=1e-6,
            maximum_iterations=200,
            charge_error=1e10,
        )

    I_after_inject = extract_contact_current(device_info, contact)
    logger.debug(f"Current after injection: {I_after_inject:.4e} A/cm")

    # Step 4: Zero generation and time-step loop.
    # The injection step IS recorded (t = dt_inject) and counted. Each BDF1 step
    # returns the contact current at the END of the step, and the backward-Euler
    # charge it moved through the contact during that step is I_k * dt_k; summing
    # those right-endpoint terms over every step, the injection step included, is
    # the charge the discretised equations actually conserve (collected_charge).
    # The previous estimator (trapezoid over the post-injection samples, stopped at
    # 1 % of the peak) dropped the injection step and the tail and mixed in a
    # quadrature the solver does not use: legacy CCE 0.9776 / 1.0066 at 10 / 50 V
    # (>1), default stack 0.982 / 0.977 (non-monotone). This sum gives legacy
    # 0.9615 / 0.9985 (steady-state CCE 0.9609 / 0.9983) and default 0.99978 / 0.99997.
    # NOTE: extract_contact_current reads ONLY the continuity-equation
    # (conduction) contact current -- it does NOT include displacement current
    # (dD/dt lives in the Poisson/PotentialEquation, which is never queried).
    add_generation_to_dd(device_info, zeros)

    times = [dt_inject]
    currents = [I_after_inject]
    t = dt_inject
    Q_run = (abs(I_after_inject) - abs(I_dark)) * dt_inject
    n_small = 0
    converged = False

    I_peak = I_after_inject
    max_steps = 5000
    step_count = 0

    while t < t_collect and step_count < max_steps:
        # Adaptive dt: 10% rule
        dt = max(dt_min, min(dt_max, t * 0.1))

        try:
            devsim.solve(
                type="transient_bdf1",
                tdelta=dt,
                absolute_error=1e10,
                relative_error=1e-8,
                maximum_iterations=100,
                charge_error=1e10,
            )
        except devsim.error:
            # Retry with relaxed tolerances
            try:
                devsim.solve(
                    type="transient_bdf1",
                    tdelta=dt,
                    absolute_error=1e12,
                    relative_error=1e-6,
                    maximum_iterations=200,
                    charge_error=1e10,
                )
            except devsim.error:
                logger.warning(f"BDF1 step failed at t={t:.3e} s, terminating early")
                break

        t += dt
        I = extract_contact_current(device_info, contact)
        times.append(t)
        currents.append(I)
        step_count += 1

        if abs(I) > abs(I_peak):
            I_peak = I

        # Convergence stop: the charge still being added is negligible against the
        # charge collected so far for n_converged consecutive steps. Validated
        # against 200 ns runs without a stop (legacy and default stack, 10 / 50 V,
        # LET 0.5-500): the truncated sum agrees to <= 1.4e-7 in CCE.
        dQ = (abs(I) - abs(I_dark)) * dt
        Q_run += dQ
        if Q_run != 0.0 and abs(dQ) < rtol_converged * abs(Q_run):
            n_small += 1
            if n_small >= n_converged:
                converged = True
                logger.debug(f"Collection converged at t={t:.3e} s")
                break
        else:
            n_small = 0

    times = np.array(times)
    currents = np.array(currents)
    if not converged:
        logger.warning(
            f"Single-particle collection not converged at t={t:.3e} s "
            f"({step_count} steps, t_collect={t_collect:.1e} s): Q_collected is "
            "truncated"
        )

    Q_collected = collected_charge(times, currents, I_dark)

    logger.info(
        f"Transient complete: {len(times)} steps, t_final={times[-1]:.3e} s, "
        f"Q_collected={Q_collected:.3e} C/cm"
    )

    return {
        "times": times,
        "currents": currents,
        "Q_collected": Q_collected,
        "I_dark": float(I_dark),
        "I_peak": float(I_peak),
        "current_contact": contact,
        "converged": converged,
    }


def collected_charge(times, currents, I_dark, t0=0.0):
    """Backward-Euler collected charge: sum_k (|I_k| - |I_dark|) * (t_k - t_{k-1}).

    ``times`` are the END times of consecutive implicit (BDF1) steps that start at
    ``t0`` (the start of the injection step), ``currents`` the contact current the
    solver returns at each of them. Each term is the charge the discretised
    equations moved through the contact during that step, so the sum is exact for
    the time-discrete system (no quadrature error on top of the time-stepping).
    """
    times = np.asarray(times, dtype=float)
    currents = np.asarray(currents, dtype=float)
    if times.size == 0:
        return 0.0
    dt = np.diff(times, prepend=t0)
    return float(np.sum((np.abs(currents) - abs(I_dark)) * dt))


def analyze_current_pulse(times, currents, I_dark):
    """Extract pulse characteristics from transient simulation.

    Uses the same backward-Euler estimator as simulate_single_particle
    (collected_charge): ``times`` are step end times measured from the injection
    start (t = 0), so Q_collected here equals simulate_single_particle's.

    Parameters
    ----------
    times : array_like
        Step end times (s), from t = 0 at the injection start.
    currents : array_like
        Contact current at each time point (A/cm for 2D).
    I_dark : float
        Baseline dark current (A/cm).

    Returns
    -------
    result : dict
        - "Q_collected": collected charge (C/cm for 2D)
        - "I_peak": peak signal current magnitude (A/cm)
        - "t_peak": time of peak current (s)
        - "t_collection": time to collect 95% of total charge (s)
    """
    times = np.asarray(times, dtype=float)
    currents = np.asarray(currents, dtype=float)

    I_signal = np.abs(currents) - np.abs(I_dark)
    Q_total = collected_charge(times, currents, I_dark)

    I_peak = float(np.max(I_signal))
    t_peak = float(times[np.argmax(I_signal)])

    # Collection time: end of the step at which the cumulative charge reaches 95%
    Q_cumulative = np.cumsum(I_signal * np.diff(times, prepend=0.0))
    idx_95 = int(np.searchsorted(Q_cumulative, 0.95 * Q_total))
    t_collection = float(times[min(idx_95, len(times) - 1)])

    return {
        "Q_collected": Q_total,
        "I_peak": I_peak,
        "t_peak": t_peak,
        "t_collection": t_collection,
    }


def build_cce_let_table(
    half_width_um=50.0,
    V_bias=50.0,
    n_let_points=40,
    let_min=0.1,
    let_max=1000.0,
    x_ion_cm=0.0,
    contact=None,
    device_kwargs=None,
):
    """Build CCE(LET) lookup table from transient simulations.

    Creates a fresh 2D device for each LET value, injects an ion track at
    the center, runs a transient simulation, and extracts CCE.  Device is
    deleted after each simulation to avoid global solver coupling.

    Parameters
    ----------
    half_width_um : float
        Half-width of sensitive volume (um).
    V_bias : float
        Reverse bias voltage (V).
    n_let_points : int
        Number of LET values to simulate.
    let_min : float
        Minimum LET (keV/um).
    let_max : float
        Maximum LET (keV/um).
    x_ion_cm : float
        Lateral position of ion track (cm).
    contact : {"cathode", "anode"} or None
        Forwarded to simulate_single_particle (None = anode with n+, else cathode).
    device_kwargs : dict or None
        Extra keyword arguments for create_2d_dd_device (geometry / doping stack).

    Returns
    -------
    table : pandas.DataFrame
        Columns: LET_keV_um, Q_generated_fC, Q_collected_fC, CCE (clipped to
        [0, 1]), CCE_raw (unclipped Q_collected / Q_generated), t_collection_ns,
        converged.
    """
    LET_values = np.logspace(np.log10(let_min), np.log10(let_max), n_let_points)

    results = []
    for i, LET in enumerate(LET_values):
        device_info = None
        try:
            device_info = create_2d_dd_device(
                half_width_um=half_width_um, V_bias=V_bias, **(device_kwargs or {})
            )

            generation, Q_gen = ion_track_generation_2d(
                device_info, LET, x_ion_cm=x_ion_cm
            )

            sim_result = simulate_single_particle(
                device_info, generation, contact=contact
            )

            Q_col = sim_result["Q_collected"]
            cce_raw = Q_col / Q_gen if Q_gen > 0 else float("nan")
            # Physical ceiling: CCE cannot exceed 1 in a linear DD detector. The
            # backward-Euler estimator gives CCE_raw <= 1 on the validated cases
            # (the old post-injection trapezoid gave 1.0066, audit C4); the clip
            # stays as a guard, with a warning if the excess is beyond numerics.
            if not np.isnan(cce_raw) and cce_raw > 1.0 + _CCE_EXCESS_TOL:
                logger.warning(
                    f"LET={LET:.1f}: CCE_raw={cce_raw:.6f} exceeds the physical "
                    f"ceiling by >{_CCE_EXCESS_TOL:g}; likely a bug "
                    f"(Q_col={Q_col:.3e}, Q_gen={Q_gen:.3e})."
                )
            cce = float(np.clip(cce_raw, 0.0, 1.0))

            pulse_info = analyze_current_pulse(
                sim_result["times"], sim_result["currents"], sim_result["I_dark"]
            )
            t_coll_ns = pulse_info["t_collection"] * 1e9

            results.append(
                {
                    "LET_keV_um": LET,
                    "Q_generated_fC": Q_gen * 1e15,
                    "Q_collected_fC": Q_col * 1e15,
                    "CCE": cce,
                    "CCE_raw": float(cce_raw),
                    "converged": bool(sim_result["converged"]),
                    "t_collection_ns": t_coll_ns,
                }
            )

        except Exception as exc:
            logger.warning(f"Simulation failed for LET={LET:.2f} keV/um: {exc}")
            results.append(
                {
                    "LET_keV_um": LET,
                    "Q_generated_fC": float("nan"),
                    "Q_collected_fC": float("nan"),
                    "CCE": float("nan"),
                    "CCE_raw": float("nan"),
                    "converged": False,
                    "t_collection_ns": float("nan"),
                }
            )

        finally:
            if device_info is not None:
                try:
                    devsim.delete_device(device=device_info["device_name"])
                except Exception:
                    pass

        # Log progress every 5 LET values
        if (i + 1) % 5 == 0 or i == 0:
            logger.info(
                f"CCE(LET) progress: {i+1}/{n_let_points} " f"(LET={LET:.2f} keV/um)"
            )

    return pd.DataFrame(results)


def save_cce_let_table(table_df, filepath, geometry_info=None):
    """Save CCE(LET) table to JSON for Phase 22 consumption.

    Parameters
    ----------
    table_df : pandas.DataFrame
        CCE table from build_cce_let_table.
    filepath : str or Path
        Output file path (.json).
    geometry_info : dict or None
        Optional geometry metadata to include. If None, uses defaults.
    """
    if geometry_info is None:
        geometry_info = {}

    data = {
        "geometry": {
            "half_width_um": geometry_info.get("half_width_um", 50.0),
            "epi_um": geometry_info.get("epi_um", 10.0),
        },
        "bias_V": geometry_info.get("bias_V", 50.0),
        "x_ion_um": geometry_info.get("x_ion_um", 0.0),
        "LET_keV_um": table_df["LET_keV_um"].tolist(),
        "CCE": table_df["CCE"].tolist(),
        "Q_generated_fC": table_df["Q_generated_fC"].tolist(),
        "Q_collected_fC": table_df["Q_collected_fC"].tolist(),
        "t_collection_ns": table_df["t_collection_ns"].tolist(),
    }
    for optional in ("CCE_raw", "converged"):  # build_cce_let_table columns
        if optional in table_df:
            data[optional] = table_df[optional].tolist()

    with open(filepath, "w") as f:
        json.dump(data, f, indent=2)

    logger.info(f"CCE(LET) table saved to {filepath} ({len(table_df)} points)")


def load_cce_let_table(filepath):
    """Load CCE(LET) table and provide interpolation function.

    Parameters
    ----------
    filepath : str or Path
        Path to JSON file from save_cce_let_table.

    Returns
    -------
    cce_interp : callable
        Function cce_interp(LET_keV_um) -> CCE, using log-linear
        interpolation on log10(LET).
    metadata : dict
        Full JSON data including geometry, bias, and raw arrays.
    """
    with open(filepath) as f:
        data = json.load(f)

    LET = np.array(data["LET_keV_um"])
    CCE = np.array(data["CCE"])

    # Filter out NaN values for interpolation
    valid = np.isfinite(LET) & np.isfinite(CCE)
    LET_valid = LET[valid]
    CCE_valid = CCE[valid]

    def cce_interp(let_value):
        """Interpolate CCE at given LET using log-linear interpolation."""
        return float(np.interp(np.log10(let_value), np.log10(LET_valid), CCE_valid))

    logger.info(
        f"Loaded CCE(LET) table from {filepath}: "
        f"{len(LET_valid)} valid points, "
        f"LET range [{LET_valid.min():.2f}, {LET_valid.max():.2f}] keV/um"
    )

    return cce_interp, data
