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

import contextlib
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
from etna.core.solver_tol import (
    DC_ABS_ERROR,
    DC_REL_ERROR,
    RELAXED_ABS_ERROR,
    RELAXED_REL_ERROR,
    TRANSIENT_ABS_ERROR,
    TRANSIENT_RELAXED_ABS_ERROR,
    TRANSIENT_RELAXED_REL_ERROR,
    TRANSIENT_REL_ERROR,
)
from etna.core.devsim_reset import delete_device_and_mesh

logger = logging.getLogger(__name__)

# Physical constants
Q_ELECTRON = 1.602e-19  # C, elementary charge

# Headroom over the physical CCE ceiling (1.0) before the table builder warns. The
# backward-Euler collected charge (collected_charge) has no quadrature error on top
# of the time stepping, and the validated runs give CCE_raw <= 1 (default stack
# 0.99997 at 50 V); an excess beyond this tolerance signals a bug, not numerics.
# (The old 5 % headroom covered the post-injection trapezoid, audit C4.)
_CCE_EXCESS_TOL = 1e-3

# build_cce_let_table(method="linear"): a LET is synthesised from the let_max transient
# once a simulated LET below it agrees with its synthesis to _LINEAR_ANCHOR_TOL
# (CCE_raw); the simulated middle LET must agree to _LINEAR_CHECK_TOL (table tolerance).
_LINEAR_ANCHOR_TOL = 2.5e-7
_LINEAR_CHECK_TOL = 1e-6
# rtol_converged of the let_max reference run: continues past the default stop so the
# synthesised low-LET series (which stop later) are covered; its own row applies the
# default stop to the same series (identical prefix, so that row is unchanged).
_LINEAR_REF_RTOL = 1e-10


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
            absolute_error=DC_ABS_ERROR,
            relative_error=DC_REL_ERROR,
            maximum_iterations=30,
        )
    except devsim.error:
        # 2D devices may need relaxed tolerances for transient_dc init
        logger.info("transient_dc failed at strict tolerances, retrying relaxed")
        devsim.solve(
            type="transient_dc",
            absolute_error=RELAXED_ABS_ERROR,
            relative_error=RELAXED_REL_ERROR,
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
            absolute_error=TRANSIENT_ABS_ERROR,
            relative_error=TRANSIENT_REL_ERROR,
            maximum_iterations=100,
            charge_error=1e10,
        )
    except devsim.error:
        logger.warning("Injection step failed at standard tolerances, retrying relaxed")
        devsim.solve(
            type="transient_bdf1",
            tdelta=dt_inject,
            absolute_error=TRANSIENT_RELAXED_ABS_ERROR,
            relative_error=TRANSIENT_RELAXED_REL_ERROR,
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
                absolute_error=TRANSIENT_ABS_ERROR,
                relative_error=TRANSIENT_REL_ERROR,
                maximum_iterations=100,
                charge_error=1e10,
            )
        except devsim.error:
            # Retry with relaxed tolerances
            try:
                devsim.solve(
                    type="transient_bdf1",
                    tdelta=dt,
                    absolute_error=TRANSIENT_RELAXED_ABS_ERROR,
                    relative_error=TRANSIENT_RELAXED_REL_ERROR,
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
    method="linear",
    n_jobs=1,
):
    """Build CCE(LET) lookup table from transient simulations.

    Builds ONE 2D device at (half_width_um, V_bias), snapshots its dark state, and
    for each LET restores that state, injects an ion track, runs a transient
    simulation and extracts CCE. The device is deleted at the end (devsim solves are
    global: no other device may coexist during the sweep).

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
    method : {"linear", "full"}
        "full": one transient per LET. "linear" (default, n_let_points > 3): the
        let_max transient is synthesised for the other LETs (see Notes); the lowest
        LETs, where the synthesis is least accurate, are simulated bottom-up until one
        agrees with its synthesis to 2.5e-7, and the middle LET is simulated as a
        check (1e-6). Falls back to simulating the remaining LETs (with a warning) if
        the check fails or a transient fails.
    n_jobs : int
        Worker processes (default 1: serial). With n_jobs > 1 the LETs that will (or
        may) be simulated run concurrently in spawned workers, each on its own fresh
        device (same builder and keywords; a fresh device equals the shared one, see
        TestBuildCCELetTableSmallMesh): method="linear" prefetches let_max, the middle
        LET and the lowest n_jobs - 2 anchors, "full" every LET. The table is assembled
        by exactly the serial logic: a prefetched anchor above the bottom-up stop is
        discarded, an anchor not prefetched is simulated in the parent. Parity with
        n_jobs=1: tests/test_cce_let_parallel.py.

    Returns
    -------
    table : pandas.DataFrame
        Columns: LET_keV_um, Q_generated_fC, Q_collected_fC, CCE (clipped to
        [0, 1]), CCE_raw (unclipped Q_collected / Q_generated), t_collection_ns,
        converged, simulated (False for the rows synthesised by method="linear").

    Notes
    -----
    Why method="linear" is enough. The table's LETs are all low injection (peak
    excess density 2.6e9-2.6e12 cm^-3 for 0.5-500 keV/um against N_D_bulk 3.4e13, no
    plasma effects in the model), so the transient is linear in the deposited charge:
    on the fixed time grid of simulate_single_particle (dt depends on t only) the
    contact current at LET L is I_dark + (I_ref(t) - I_dark) * Q_gen(L) / Q_gen(ref).
    The synthesis applies simulate_single_particle's own estimator and convergence
    stop to that series, so it reproduces the per-LET runs, including their residual
    ~1e-5 LET dependence. That dependence is numerical, from collected_charge's
    sum (|I| - |I_dark|) dt: at the anode of the default stack I_dark (+2.0e-13 A/cm
    at hw 50, 50 V) has the opposite sign to the signal current, an offset of order
    I_dark * t_end / Q_gen ~ 1/LET (the sign-aware sum (I_dark - I) dt gives
    0.9999715 / 0.9999703 / 0.9999703 at LET 0.5 / 15.8 / 500 against 0.9999624 /
    0.9999699 / 0.9999703). The Newton tolerances are not the cause
    (TRANSIENT_ABS_ERROR 1e10 -> 1e7 bit-identical, TRANSIENT_REL_ERROR 1e-8 -> 1e-11
    moves CCE by 3e-9; lab/engine_box/bench/let_table_abs_tol.py). Synthesis vs
    simulation (fresh runs, hw 50 / 150): 7.9e-8 / 7.4e-8 at LET 15.8, ~1e-6 at LET
    0.5 (there the solver's tail stalls at I_dark, which a linear model does not
    reproduce: hence the bottom-up simulated anchors). Parity on the full notebook-16
    grids: lab/engine_box/bench/let_table_linear.py.
    """
    if method not in ("linear", "full"):
        raise ValueError(f"method must be 'linear' or 'full', got {method!r}")
    LET_values = np.logspace(np.log10(let_min), np.log10(let_max), n_let_points)

    # One device for the whole table: the build (mesh, Poisson, equilibrium, DD setup,
    # bias ramp; ~25 s at hw 50, a quarter of a point) does not depend on the LET. The
    # dark DC state (Potential / Electrons / Holes at V_bias) is snapshotted once and
    # restored before each LET; simulate_single_particle's transient_dc then resets
    # the time-integration history, as a fresh device would (the transient_cce_vs_
    # dose_rate pattern, WP-F). Results are identical to the former build-per-LET
    # loop (lab/engine_box/bench/let_table_build_once.py).
    pool_stack, prefetch = contextlib.ExitStack(), {}
    with pool_stack:
        if n_jobs is not None and n_jobs > 1:
            # workers start before the parent's own build, so the builds overlap
            prefetch = _prefetch_let_points(
                pool_stack,
                _prefetch_plan(method, n_let_points, n_jobs),
                n_jobs,
                (half_width_um, V_bias, device_kwargs, LET_values, x_ion_cm, contact),
            )
        return _build_cce_let_table_on_device(
            half_width_um,
            V_bias,
            device_kwargs,
            LET_values,
            x_ion_cm,
            contact,
            method,
            prefetch,
        )


def _build_cce_let_table_on_device(
    half_width_um,
    V_bias,
    device_kwargs,
    LET_values,
    x_ion_cm,
    contact,
    method,
    prefetch,
):
    """build_cce_let_table's serial driver; ``prefetch`` maps a grid index to
    (sim_kw, future) of a worker that simulates that LET on its own device."""
    n_let_points = len(LET_values)
    device_info = create_2d_dd_device(
        half_width_um=half_width_um, V_bias=V_bias, **(device_kwargs or {})
    )
    try:
        dark_state = _snapshot_dark_state(device_info)

        def point(i, return_series=False, **sim_kw):
            job = prefetch.pop(i, None)
            if job is not None and job[0] == sim_kw:
                try:
                    row, series = job[1].result()
                    return (row, series) if return_series else row
                except Exception as exc:  # a dead worker: simulate it here instead
                    logger.warning(f"CCE(LET) worker failed for LET index {i}: {exc!r}")
            return _cce_let_point(
                device_info,
                dark_state,
                LET_values[i],
                x_ion_cm,
                contact,
                i,
                n_let_points,
                return_series=return_series,
                **sim_kw,
            )

        rows, done = None, {}
        if method == "linear" and n_let_points > 3:
            rows, done = _linear_table(device_info, LET_values, x_ion_cm, point)
        if rows is None:  # "full", or the linear check failed: simulate the rest
            rows = [done[i] if i in done else point(i) for i in range(n_let_points)]
    finally:
        try:
            delete_device_and_mesh(device_info["device_name"])
        except Exception:
            pass
        for (
            _,
            item,
        ) in prefetch.values():  # speculative anchors the serial logic skipped
            item.cancel()

    return pd.DataFrame(rows)


def _prefetch_plan(method, n_let_points, n_jobs):
    """[(grid index, sim_kw)] that build_cce_let_table(n_jobs > 1) simulates in workers:
    the points the serial logic is certain (let_max, middle) or likely (lowest
    anchors) to ask for; "full" (and n <= 3, which runs full) asks for every LET."""
    if method != "linear" or n_let_points <= 3:
        return [(i, {}) for i in range(n_let_points)]
    i_hi, i_mid = n_let_points - 1, n_let_points // 2
    plan = [(i_hi, {"rtol_converged": _LINEAR_REF_RTOL}), (i_mid, {})]
    anchors = [i for i in range(i_hi) if i != i_mid][: max(n_jobs - 2, 1)]
    return plan + [(i, {}) for i in anchors]


def _prefetch_let_points(stack, plan, n_jobs, build_args):
    """Start ``plan`` on a spawn pool entered on ``stack``; {index: (sim_kw, item)}.
    The jobs are dealt round-robin into at most n_jobs chunks, one device build per
    chunk (the plan's order puts let_max, the longest run, first)."""
    from etna.core.optimization import _spawn_pool

    n_chunks = min(n_jobs, len(plan))
    chunks = [plan[k::n_chunks] for k in range(n_chunks)]
    pool = stack.enter_context(_spawn_pool(n_chunks))
    out = {}
    for chunk in chunks:
        fut = pool.submit(_cce_let_chunk, build_args, chunk)
        for j, (i, kw) in enumerate(chunk):
            out[i] = (kw, _ChunkItem(fut, j))
    return out


class _ChunkItem:
    """Future-like view of job ``j`` of a chunk future (result / cancel)."""

    def __init__(self, fut, j):
        self._fut, self._j = fut, j

    def result(self):
        return self._fut.result()[self._j]

    def cancel(self):
        self._fut.cancel()


def _cce_let_chunk(build_args, jobs):
    """Worker of build_cce_let_table(n_jobs > 1): one fresh device, then the serial
    point routine for each (i, sim_kw); returns [(row, series)]."""
    from etna.core.devsim_lock import quiet_devsim

    half_width_um, V_bias, device_kwargs, LET_values, x_ion_cm, contact = build_args
    with quiet_devsim():
        info = create_2d_dd_device(
            half_width_um=half_width_um, V_bias=V_bias, **(device_kwargs or {})
        )
        try:
            dark = _snapshot_dark_state(info)
            return [
                _cce_let_point(
                    info,
                    dark,
                    LET_values[i],
                    x_ion_cm,
                    contact,
                    i,
                    len(LET_values),
                    return_series=True,
                    **kw,
                )
                for i, kw in jobs
            ]
        finally:
            delete_device_and_mesh(info["device_name"])


def synthesize_linear_response(
    times, currents, I_dark, scale, rtol_converged=1e-7, n_converged=3
):
    """simulate_single_particle's (Q_collected, times, currents) for a generation
    ``scale`` times that of a run (``times``, ``currents``, ``I_dark``), assuming the
    transient is linear in the generation: I = I_dark + scale * (I_run - I_dark) on the
    same time grid, with the same estimator and convergence stop. None if the series
    ends before the stop (the run being scaled is too short to cover it).
    """
    times = np.asarray(times, dtype=float)
    I = I_dark + scale * (np.asarray(currents, dtype=float) - I_dark)
    dt = np.diff(times, prepend=0.0)
    Q_run = (abs(I[0]) - abs(I_dark)) * dt[0]
    n_small = 0
    for k in range(1, len(times)):
        dQ = (abs(I[k]) - abs(I_dark)) * dt[k]
        Q_run += dQ
        if Q_run != 0.0 and abs(dQ) < rtol_converged * abs(Q_run):
            n_small += 1
            if n_small >= n_converged:
                t, Ik = times[: k + 1], I[: k + 1]
                return collected_charge(t, Ik, I_dark), t, Ik
        else:
            n_small = 0
    return None


def _linear_table(device_info, LET_values, x_ion_cm, point):
    """(rows, simulated) for build_cce_let_table(method="linear"); rows is None when the
    let_max transient fails or the middle-LET check misses (the caller then simulates
    the rest). ``simulated`` maps grid index -> simulated row."""
    n = len(LET_values)
    i_hi, i_mid = n - 1, n // 2
    row_hi, ref = point(i_hi, rtol_converged=_LINEAR_REF_RTOL, return_series=True)
    if ref is None or not np.isfinite(row_hi["CCE_raw"]):
        logger.warning("linear CCE(LET) table: let_max transient failed")
        return None, {}
    q_gen_hi = row_hi["Q_generated_fC"] * 1e-15

    def synth(i, Q_gen=None):
        if Q_gen is None:
            _, Q_gen = ion_track_generation_2d(
                device_info, LET_values[i], x_ion_cm=x_ion_cm
            )
        out = synthesize_linear_response(
            ref["times"], ref["currents"], ref["I_dark"], Q_gen / q_gen_hi
        )
        if out is None:
            return None
        Q_col, t, I = out
        t95 = analyze_current_pulse(t, I, ref["I_dark"])["t_collection"]
        row = _table_row(LET_values[i], Q_gen, Q_col, True, t95)
        return {**row, "simulated": False}

    # the let_max row: the default stop applied to its own (longer) series
    hi = synth(i_hi, q_gen_hi)
    if hi is None:
        logger.warning("linear CCE(LET) table: let_max run did not converge")
        return None, {}
    sim = {i_hi: {**hi, "simulated": True}}

    def miss(i):
        s_row = synth(i, sim[i]["Q_generated_fC"] * 1e-15)
        if s_row is None or not np.isfinite(sim[i]["CCE_raw"]):
            return np.inf
        return abs(s_row["CCE_raw"] - sim[i]["CCE_raw"])

    # bottom-up anchors: simulate until a LET agrees with its synthesis
    i_first_synth = i_hi
    for i in range(i_hi):
        sim[i] = point(i)
        if miss(i) <= _LINEAR_ANCHOR_TOL:
            i_first_synth = i + 1
            break
    if i_mid not in sim:
        sim[i_mid] = point(i_mid)
    d_mid = miss(i_mid)
    if not d_mid <= _LINEAR_CHECK_TOL:
        logger.warning(
            f"linear CCE(LET) table: middle LET {LET_values[i_mid]:.3g} keV/um misses "
            f"its synthesis by {d_mid:.2e} in CCE (> {_LINEAR_CHECK_TOL:g}): "
            "simulating the other LETs"
        )
        return None, sim
    rows = []
    for i in range(n):
        row = sim.get(i)
        if row is None:
            row = synth(i) if i >= i_first_synth else None
        if row is None:
            row = sim[i] = point(i)
        rows.append(row)
    logger.info(
        f"linear CCE(LET) table: {len(sim)} of {n} LETs simulated, middle-LET check "
        f"|dCCE| = {d_mid:.2e}"
    )
    return rows, sim


# DC solution variables of the 2D DD device: the complete state create_2d_dd_device
# leaves after the bias ramp (everything else is derived from them at each solve).
_DARK_STATE_SOLUTIONS = ("Potential", "Electrons", "Holes")


def _snapshot_dark_state(device_info):
    device, region = device_info["device_name"], device_info["region_name"]
    return {
        name: list(
            devsim.get_node_model_values(device=device, region=region, name=name)
        )
        for name in _DARK_STATE_SOLUTIONS
    }


def _restore_dark_state(device_info, dark_state):
    device, region = device_info["device_name"], device_info["region_name"]
    for name, values in dark_state.items():
        devsim.set_node_values(device=device, region=region, name=name, values=values)


def _cce_let_point(
    device_info,
    dark_state,
    LET,
    x_ion_cm,
    contact,
    i,
    n_let_points,
    return_series=False,
    **sim_kw,
):
    """One CCE(LET) table row on the shared device, from the restored dark state
    (with ``return_series``: (row, simulate_single_particle result or None))."""
    sim_result = None
    try:
        _restore_dark_state(device_info, dark_state)
        generation, Q_gen = ion_track_generation_2d(device_info, LET, x_ion_cm=x_ion_cm)
        sim_result = simulate_single_particle(
            device_info, generation, contact=contact, **sim_kw
        )
        pulse_info = analyze_current_pulse(
            sim_result["times"], sim_result["currents"], sim_result["I_dark"]
        )
        row = _table_row(
            LET,
            Q_gen,
            sim_result["Q_collected"],
            sim_result["converged"],
            pulse_info["t_collection"],
        )
        row["simulated"] = True
    except Exception as exc:
        logger.warning(f"Simulation failed for LET={LET:.2f} keV/um: {exc}")
        row = {
            "LET_keV_um": LET,
            "Q_generated_fC": float("nan"),
            "Q_collected_fC": float("nan"),
            "CCE": float("nan"),
            "CCE_raw": float("nan"),
            "converged": False,
            "t_collection_ns": float("nan"),
            "simulated": True,
        }

    # Log progress every 5 LET values
    if (i + 1) % 5 == 0 or i == 0:
        logger.info(f"CCE(LET) progress: {i+1}/{n_let_points} (LET={LET:.2f} keV/um)")
    return (row, sim_result) if return_series else row


def _table_row(LET, Q_gen, Q_col, converged, t_collection_s):
    """CCE(LET) table row (charges in C/cm, t_collection in s)."""
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
    return {
        "LET_keV_um": LET,
        "Q_generated_fC": Q_gen * 1e15,
        "Q_collected_fC": Q_col * 1e15,
        "CCE": float(np.clip(cce_raw, 0.0, 1.0)),
        "CCE_raw": float(cce_raw),
        "converged": bool(converged),
        "t_collection_ns": t_collection_s * 1e9,
    }


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

    log_let_valid = np.log10(LET_valid)

    def cce_interp(let_value):
        """Interpolate CCE at given LET using log-linear interpolation (a scalar LET
        gives a float, an array LET an array of the same per-element values)."""
        cce = np.interp(np.log10(let_value), log_let_valid, CCE_valid)
        return float(cce) if np.ndim(cce) == 0 else cce

    # accepts arrays: process_mc_ensemble then looks the events up in one call
    cce_interp.vectorized = True

    # np.interp clamps outside the table: expose the range so callers can flag it
    # (process_mc_ensemble warns on clamped events; audit-v7 B-11)
    cce_interp.let_min = float(LET_valid.min())
    cce_interp.let_max = float(LET_valid.max())

    logger.info(
        f"Loaded CCE(LET) table from {filepath}: "
        f"{len(LET_valid)} valid points, "
        f"LET range [{LET_valid.min():.2f}, {LET_valid.max():.2f}] keV/um"
    )

    return cce_interp, data
