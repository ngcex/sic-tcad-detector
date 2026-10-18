"""Transient drift-diffusion solver for FLASH pulse dynamics in 4H-SiC.

Provides time-resolved simulation of carrier dynamics during pulsed proton
irradiation with adaptive time-stepping that spans the 6-order timescale
gap between microsecond pulse rise times and millisecond pulse durations.

Components:
- pulse_envelope: trapezoidal pulse shape function
- adaptive_dt: time-step selection based on pulse phase
- TransientSolver: main class wrapping devsim transient solve API

The solver uses BDF1 (backward differentiation, first order) for
unconditional stability with stiff semiconductor equations. Generation
rate is modulated at each time step via set_node_values on RadGenRate.

All units CGS (cm, cm^-3, s, A/cm^2) per devsim convention.

References:
    - devsim tran_diode.py example for transient solve pattern
    - Phase 7 decision: transient_dc initialization before BDF1
    - Phase 3 decision: RadGenRate set via set_node_values (data, not expression)
"""

import logging
import uuid

import devsim
import numpy as np
import pandas as pd

from etna.core.charge_collection import add_generation_to_dd
from etna.core.drift_diffusion import (
    create_dd_device,
    extract_contact_current,
    ramp_bias,
)
from etna.core.flash_recombination import add_auger_recombination
from etna.core.generation_profiles import (
    confine_generation_to_epi,
    proton_generation_profile,
)
from etna.core.device import DEFAULT_STACK as _STACK, graded_stack_kwargs
from etna.core.solver_tol import (
    DC_ABS_ERROR,
    DC_REL_ERROR,
    RELAXED_ABS_ERROR,
    RELAXED_REL_ERROR,
)

logger = logging.getLogger(__name__)

# Physical constants (CGS)
Q_ELECTRON = 1.602e-19  # C (elementary charge)

# DC solution variables of the 1D DD device: the complete state that ramp_bias leaves
# (transient_cce_vs_dose_rate restores them per dose rate instead of rebuilding).
_DARK_STATE_SOLUTIONS = ("Potential", "Electrons", "Holes")


def pulse_envelope(t, t_rise, t_duration, t_fall):
    """Trapezoidal pulse envelope returning value in [0, 1].

    Stages:
    - [0, t_rise): linear ramp 0 -> 1
    - [t_rise, t_rise + t_duration): hold at 1.0
    - [t_rise + t_duration, t_rise + t_duration + t_fall): linear ramp 1 -> 0
    - [t_rise + t_duration + t_fall, ...): 0

    Parameters
    ----------
    t : float
        Time (s).
    t_rise : float
        Rise time (s).
    t_duration : float
        Plateau duration (s).
    t_fall : float
        Fall time (s).

    Returns
    -------
    envelope : float
        Envelope value in [0, 1].
    """
    if t < 0:
        return 0.0
    elif t < t_rise:
        return t / t_rise if t_rise > 0 else 1.0
    elif t < t_rise + t_duration:
        return 1.0
    elif t < t_rise + t_duration + t_fall:
        if t_fall > 0:
            return 1.0 - (t - t_rise - t_duration) / t_fall
        else:
            return 0.0
    else:
        return 0.0


def generated_charge_trapezoidal_pulse(
    G_total, t_rise, t_duration, t_fall, q=1.602e-19
):
    """Charge generated over a full trapezoidal pulse envelope.

    Audit C4 fix. The generation rate is modulated by pulse_envelope(t), a
    trapezoid: linear 0->1 over t_rise, flat 1 over t_duration, linear 1->0 over
    t_fall. The time-integral of that envelope is the trapezoid area

        integral(envelope dt) = t_rise/2 + t_duration + t_fall/2

    so the total generated charge is

        Q_gen = q * G_total * (t_rise/2 + t_duration + t_fall/2)

    The previous implementation used only `q * G_total * t_duration` (plateau
    only), omitting the charge generated during the rise and fall ramps. Because
    the collected-charge integral spans the WHOLE waveform, that mismatch made
    CCE = Q_collected / Q_generated exceed 1 (the artifact the old [0, 2] clip
    was masking). At default pulses (t_rise = t_fall = 1 us, t_duration = 1 ms)
    the error is ~0.1%, but for fast pulses with t_rise ~ t_duration it can
    approach +50-100%, so this is a genuine normalization bug, not quadrature
    noise.

    Parameters
    ----------
    G_total : float
        Peak generation rate integrated over space (per-area rate; the same
        quantity the collected-current integral is normalized against).
    t_rise, t_duration, t_fall : float
        Pulse rise time, plateau duration, and fall time (s).
    q : float
        Elementary charge (C). Default 1.602e-19.

    Returns
    -------
    float
        Generated charge over the full envelope (C per unit area).
    """
    envelope_integral = 0.5 * t_rise + t_duration + 0.5 * t_fall
    return q * G_total * envelope_integral


def adaptive_dt(t, t_rise, t_duration, t_fall, dt_min=1e-8, dt_max=1e-4):
    """Select time step based on pulse phase.

    During rise/fall transitions, uses small dt (t_rise/10 or dt_min).
    During plateau and post-pulse, uses dt_max.

    Parameters
    ----------
    t : float
        Current time (s).
    t_rise : float
        Pulse rise time (s).
    t_duration : float
        Pulse plateau duration (s).
    t_fall : float
        Pulse fall time (s).
    dt_min : float
        Minimum time step (s). Default: 1e-8.
    dt_max : float
        Maximum time step (s). Default: 1e-4.

    Returns
    -------
    dt : float
        Selected time step (s), always in [dt_min, dt_max].
    """
    t_plateau_end = t_rise + t_duration
    t_end_fall = t_plateau_end + t_fall

    # During rise phase
    if t < t_rise:
        dt_transition = t_rise / 10 if t_rise > 0 else dt_min
        return max(dt_min, min(dt_transition, dt_max))

    # During plateau
    if t < t_plateau_end:
        return dt_max

    # During fall phase
    if t < t_end_fall:
        dt_transition = t_fall / 10 if t_fall > 0 else dt_min
        return max(dt_min, min(dt_transition, dt_max))

    # Post-pulse
    return dt_max


def pulse_time_step(t, t_rise, t_duration, t_fall, t_end, dt_min=1e-8, dt_max=1e-4):
    """Next step (dt, t_next) of the pulse stepper: adaptive_dt clipped to the next
    phase breakpoint (rise end, plateau end, fall end, t_end).

    ``t_next`` is the breakpoint itself (not ``t + dt`` with its round-off) when the step
    was clipped, so the stepper lands exactly on every phase edge. Before 2026-10-18 only
    t_end clipped: a plateau step of dt_max overran the plateau end whenever t_duration
    was not a multiple of it, the fall ramp was then taken at 100 % generation and never
    resolved (.issues/audit-v7 A-04: +1.0 % generated charge at t_duration 0.55 ms,
    +10.7 % at 50 us; the 1 ms default happened to be commensurate).
    """
    breakpoints = (t_rise, t_rise + t_duration, t_rise + t_duration + t_fall, t_end)
    tol = 1e-6 * dt_min  # round-off of accumulated steps (1e-14 s at the default)
    for b in breakpoints:  # a t within round-off of an edge IS on the edge
        if abs(t - b) <= tol:
            t = b
    dt = adaptive_dt(t, t_rise, t_duration, t_fall, dt_min, dt_max)
    for b in breakpoints:
        if b > t:
            if t + dt >= b - tol:
                return b - t, b
            return dt, t + dt
    return dt, t + dt  # past t_end (callers stop before)


class TransientSolver:
    """Transient drift-diffusion solver for FLASH pulse simulation.

    Wraps the devsim transient solve API with adaptive time-stepping
    and pulse envelope modulation of the generation rate.

    Parameters
    ----------
    device_info : dict
        Device info dict from create_dd_device(), already biased and
        with Auger recombination added. Must have dd_initialized=True.
    contact : {"cathode", "anode"} or None
        Contact at which to extract current. None (default) = "anode" when the device
        has an n+ layer (the direct cathode current on the n+ is round-off, see
        transient_cce_vs_dose_rate), else "cathode" (the historical choice).
    method : str
        Transient solve method. Default: "transient_bdf1".
    """

    def __init__(self, device_info, contact=None, method="transient_bdf1"):
        self.device_info = device_info
        if contact is None:
            has_nplus = device_info.get("nplus_thickness_cm") is not None
            contact = "anode" if has_nplus else "cathode"
        self.contact = contact
        self.method = method

    def initialize(self):
        """Initialize transient state with zero generation.

        CRITICAL: Must be called before simulate_pulse(). Sets zero
        generation and performs transient_dc solve to establish the
        initial charge state needed by BDF methods.
        """
        device = self.device_info["device_name"]
        region = self.device_info["region_name"]

        # Get mesh size for zero array
        x_nodes = np.array(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        zeros = np.zeros_like(x_nodes)

        # Set zero generation
        add_generation_to_dd(self.device_info, zeros)

        # transient_dc solve to initialize time data (Phase 7 pattern); the shared DC
        # tolerance (was 1e-10: Newton is quadratic, the 1D outputs are identical)
        devsim.solve(
            type="transient_dc",
            absolute_error=DC_ABS_ERROR,
            relative_error=DC_REL_ERROR,
            maximum_iterations=30,
        )
        # dark baseline of the state just solved (simulate_pulse(skip_init=True))
        self._I_dark = extract_contact_current(self.device_info, contact=self.contact)
        logger.info("Transient state initialized (transient_dc solve complete)")

    def simulate_pulse(
        self,
        G_spatial,
        t_rise=1e-6,
        t_duration=1e-3,
        t_fall=1e-6,
        t_post=None,
        dt_min=1e-8,
        dt_max=1e-4,
        dose_rate_Gy_s=None,
        skip_init=False,
    ):
        """Simulate a single FLASH pulse with adaptive time-stepping.

        Parameters
        ----------
        G_spatial : array_like
            Spatial generation profile at mesh nodes (cm^-3 s^-1).
            This is the peak generation rate; modulated by envelope(t).
        t_rise : float
            Pulse rise time (s). Default: 1e-6 (1 us).
        t_duration : float
            Pulse plateau duration (s). Default: 1e-3 (1 ms).
        t_fall : float
            Pulse fall time (s). Default: 1e-6 (1 us).
        t_post : float or None
            Post-pulse decay time (s). Default: 5 * t_fall.
        dt_min : float
            Minimum time step (s). Default: 1e-8.
        dt_max : float
            Maximum time step (s). Default: 1e-4.
        dose_rate_Gy_s : float or None
            Dose rate for metadata (Gy/s). Not used in computation.
        skip_init : bool
            If True, do not re-measure the dark current: the baseline is the one
            measured by initialize() (or by the previous pulse). Used by
            simulate_pulse_train, where the device still carries the tail of the
            previous pulse, which a fresh measurement would take as "dark" current.
            Default: False. (Before 2026-10-18 the flag was ignored, A-05.)

        Returns
        -------
        result : dict
            Dictionary with:
            - "times": np.array of step END times (s, from the pulse start t = 0):
              the backward-Euler (BDF1) step solves the state at t + dt, so its
              current is stamped t + dt (before 2026-10-18 it was stamped t, one step
              early, A-04)
            - "currents": np.array of contact currents (A/cm^2) at "times"
            - "I_dark": float, dark current before pulse (A/cm^2)
            - "generation_integral": float, sum over the steps of envelope(t + dt) * dt
              (s): the generation the discrete system actually received, per unit
              G_spatial
            - "t_rise": float
            - "t_duration": float
            - "t_fall": float
            - "dose_rate_Gy_s": float or None

        Notes
        -----
        Each step is clipped to the next phase breakpoint (pulse_time_step), and the
        generation of the implicit step is the envelope at the step END, t + dt.
        """
        G_spatial = np.asarray(G_spatial, dtype=float)

        if t_post is None:
            t_post = 5 * t_fall

        t_end = t_rise + t_duration + t_fall + t_post

        # Dark current before the pulse (skip_init: keep the stored baseline)
        if skip_init and getattr(self, "_I_dark", None) is not None:
            I_dark = self._I_dark
        else:
            I_dark = extract_contact_current(self.device_info, contact=self.contact)
            self._I_dark = I_dark
        logger.info(f"Dark current: {I_dark:.4e} A/cm^2")

        # Time-stepping loop
        t = 0.0
        times = []
        currents = []
        gen_integral = 0.0
        step_count = 0

        while t < t_end:
            dt, t_next = pulse_time_step(
                t, t_rise, t_duration, t_fall, t_end, dt_min, dt_max
            )
            # Implicit (BDF1) step to t_next: generation at the step end
            env = pulse_envelope(t_next, t_rise, t_duration, t_fall)
            add_generation_to_dd(self.device_info, G_spatial * env)

            # Transient solve with retry at relaxed tolerances.
            # charge_error=1e10 effectively disables automatic step rejection
            # based on charge conservation error -- we manage time steps
            # ourselves via pulse_time_step based on the pulse envelope phase.
            # The first attempt uses the shared DC tolerance (was 1e-10); the retry
            # differs only in absolute_error and the iteration cap.
            try:
                devsim.solve(
                    type=self.method,
                    absolute_error=DC_ABS_ERROR,
                    relative_error=DC_REL_ERROR,
                    maximum_iterations=40,
                    tdelta=dt,
                    charge_error=1e10,
                )
            except devsim.error:
                # Retry with relaxed tolerances and more iterations
                devsim.solve(
                    type=self.method,
                    absolute_error=RELAXED_ABS_ERROR,
                    relative_error=RELAXED_REL_ERROR,
                    maximum_iterations=100,
                    tdelta=dt,
                    charge_error=1e10,
                )
                logger.debug(
                    f"Step {step_count}: converged with relaxed tolerances "
                    f"at t={t:.4e}, dt={dt:.4e}"
                )

            # Extract current: the state (and current) at the step end
            I_t = extract_contact_current(self.device_info, contact=self.contact)
            times.append(t_next)
            currents.append(I_t)
            gen_integral += env * dt

            t = t_next
            step_count += 1

            if step_count % 50 == 0:
                logger.info(
                    f"Step {step_count}: t={t:.4e} s, env={env:.3f}, "
                    f"I={I_t:.4e} A/cm^2, dt={dt:.4e} s"
                )

        logger.info(f"Pulse simulation complete: {step_count} steps, t_end={t:.4e} s")

        return {
            "times": np.array(times),
            "currents": np.array(currents),
            "I_dark": float(I_dark),
            "generation_integral": float(gen_integral),
            "t_rise": t_rise,
            "t_duration": t_duration,
            "t_fall": t_fall,
            "dose_rate_Gy_s": dose_rate_Gy_s,
        }

    def compute_transient_cce(self, result, G_spatial, x_nodes, return_raw=False):
        """Compute charge collection efficiency from transient simulation.

        Backward-Euler (BDF1) estimator, the time-discrete charge balance of the
        solved system (the estimator of single_particle.collected_charge):

            Q_collected = sum_k (|I_k| - |I_dark|) * dt_k
            Q_generated = q * integral(G dx) * sum_k envelope(t_k) * dt_k

        with t_k the step END times ("times") and dt_k their spacing from t = 0. Both
        sums run over the same steps, so the ratio is exact for the discrete system
        (no quadrature error on top of the time stepping). The analytic trapezoid
        area ``generated_charge_trapezoidal_pulse`` (t_rise/2 + t_duration + t_fall/2)
        is logged as a cross-check: with ramps of 10 equal clipped steps the discrete
        rise and fall sums are 0.55 t_rise and 0.45 t_fall, equal to it when
        t_rise == t_fall.

        Before 2026-10-18 this was the trapezoid of I over start-of-step stamps,
        normalised by the analytic area, on a stepper that overran the plateau end
        (A-04): 0.0098 instead of ~1 at t_duration 50 us, raw 1.00028 (> 1, clipped)
        at the 1 ms default and 20 Gy/s. Kept differences from single_particle: a
        trapezoidal pulse envelope instead of a one-step impulse, and the 1D
        cathode/anode contact choice of TransientSolver.

        Parameters
        ----------
        result : dict
            Result dictionary from simulate_pulse() (a single pulse; times from 0).
        G_spatial : array_like
            Spatial generation profile at mesh nodes (cm^-3 s^-1).
        x_nodes : array_like
            Mesh node positions (cm).
        return_raw : bool
            Also return the unclipped ratio Q_collected / Q_generated. Default False.

        Returns
        -------
        cce : float
            Transient charge collection efficiency in [0, 1] (the raw ratio clipped to
            the physical ceiling 1; a raw value above 1.05 is logged as a bug).
        cce_raw : float
            Only with ``return_raw=True``: the unclipped ratio, so a clipped 1.0 can be
            told apart from a measured one.
        """
        times = np.asarray(result["times"], dtype=float)
        currents = np.asarray(result["currents"], dtype=float)
        I_dark = result["I_dark"]
        t_duration = result["t_duration"]
        t_rise = result.get("t_rise", 0.0)
        t_fall = result.get("t_fall", 0.0)

        # Collected charge: right-endpoint (backward-Euler) sum over every step of the
        # WHOLE waveform (rise + plateau + fall + post-pulse decay).
        dts = np.diff(times, prepend=0.0)
        I_signal = np.abs(currents) - np.abs(I_dark)
        Q_collected = float(np.sum(I_signal * dts))

        # Total generation rate per unit area: integral of G(x) dx
        G_total = np.trapezoid(
            np.asarray(G_spatial, dtype=float), np.asarray(x_nodes, dtype=float)
        )

        # Generated charge: the envelope the solver applied, summed over the same
        # steps (simulate_pulse "generation_integral"); analytic area as cross-check.
        Q_analytic = generated_charge_trapezoidal_pulse(
            G_total, t_rise, t_duration, t_fall, q=Q_ELECTRON
        )
        if "generation_integral" in result:
            Q_generated = Q_ELECTRON * G_total * result["generation_integral"]
        else:  # a result dict from elsewhere: the analytic envelope area
            Q_generated = Q_analytic

        if Q_generated <= 0:
            logger.warning("Q_generated <= 0, returning CCE = 0")
            return (0.0, float("nan")) if return_raw else 0.0

        cce_raw = Q_collected / Q_generated

        # Physical ceiling: in a linear DD detector CCE cannot exceed 1. The discrete
        # balance is exact up to the solver tolerance, so a larger excess means a real
        # bug (window/normalisation mismatch) and is reported, not silently clipped.
        _QUADRATURE_TOL = 0.05  # 5% headroom over the physical ceiling
        if cce_raw > 1.0 + _QUADRATURE_TOL:
            logger.warning(
                f"Transient CCE_raw={cce_raw:.4f} exceeds physical ceiling by "
                f">{_QUADRATURE_TOL:.0%}; likely a normalization/window bug "
                f"(Q_collected={Q_collected:.4e}, Q_generated={Q_generated:.4e})."
            )
        cce = float(np.clip(cce_raw, 0.0, 1.0))

        logger.info(
            f"Transient CCE: {cce:.4f} (raw={cce_raw:.6f}, "
            f"Q_collected={Q_collected:.4e}, Q_generated={Q_generated:.4e}, "
            f"analytic envelope {Q_analytic:.4e})"
        )
        return (cce, float(cce_raw)) if return_raw else cce


def simulate_pulse_train(
    device_info,
    G_spatial,
    n_pulses=10,
    t_rise=1e-6,
    t_duration=1e-3,
    t_fall=1e-6,
    t_gap=1e-3,
    dt_min=1e-8,
    dt_max=1e-4,
    contact=None,
    dose_rate_Gy_s=None,
):
    """Simulate N consecutive FLASH pulses with inter-pulse carrier memory.

    Creates a TransientSolver, initializes transient state, then loops over
    ``n_pulses`` calls to ``simulate_pulse``. Device state persists between
    pulses (no re-initialization), so inter-pulse carrier memory effects are
    captured naturally by devsim.

    Parameters
    ----------
    device_info : dict
        Device info dict from create_dd_device(), already biased and with
        Auger recombination added.
    G_spatial : array_like
        Spatial generation profile at mesh nodes (cm^-3 s^-1).
    n_pulses : int
        Number of consecutive pulses. Default: 10.
    t_rise : float
        Pulse rise time (s). Default: 1e-6.
    t_duration : float
        Pulse plateau duration (s). Default: 1e-3.
    t_fall : float
        Pulse fall time (s). Default: 1e-6.
    t_gap : float
        Inter-pulse gap duration (s). Default: 1e-3.
    dt_min : float
        Minimum time step (s). Default: 1e-8.
    dt_max : float
        Maximum time step (s). Default: 1e-4.
    contact : {"cathode", "anode"} or None
        Contact for current extraction. None (default) = anode with n+, else cathode
        (see TransientSolver).
    dose_rate_Gy_s : float or None
        Dose rate for metadata. Not used in computation.

    Returns
    -------
    result : dict
        Dictionary with:
        - "times": np.array of concatenated time points (s), continuous
        - "currents": np.array of concatenated contact currents (A/cm^2)
        - "n_pulses": int, number of pulses simulated
        - "pulse_times": list of per-pulse result dicts from simulate_pulse
        - "I_dark": float, dark current before first pulse (A/cm^2)
        - "t_rise", "t_duration", "t_fall", "t_gap": pulse parameters
        - "dose_rate_Gy_s": dose rate metadata
    """
    G_spatial = np.asarray(G_spatial, dtype=float)

    solver = TransientSolver(device_info, contact=contact)
    solver.initialize()

    # Dark current before the first pulse (measured by initialize); every pulse uses
    # it as its baseline (skip_init), not the tail of the previous pulse (A-05)
    I_dark = solver._I_dark

    all_times = []
    all_currents = []
    pulse_results = []
    t_offset = 0.0

    for i in range(n_pulses):
        logger.info(f"Pulse {i + 1}/{n_pulses}: t_offset={t_offset:.4e} s")

        # Use t_gap as post-pulse decay time (inter-pulse gap)
        result_i = solver.simulate_pulse(
            G_spatial,
            t_rise=t_rise,
            t_duration=t_duration,
            t_fall=t_fall,
            t_post=t_gap,
            dt_min=dt_min,
            dt_max=dt_max,
            dose_rate_Gy_s=dose_rate_Gy_s,
            skip_init=True,
        )

        pulse_results.append(result_i)

        # Offset times to make continuous timeline
        all_times.append(result_i["times"] + t_offset)
        all_currents.append(result_i["currents"])

        # Advance offset by this pulse's total duration
        t_offset += t_rise + t_duration + t_fall + t_gap

    logger.info(
        f"Pulse train complete: {n_pulses} pulses, " f"total time = {t_offset:.4e} s"
    )

    return {
        "times": np.concatenate(all_times),
        "currents": np.concatenate(all_currents),
        "n_pulses": n_pulses,
        "pulse_times": pulse_results,
        "I_dark": float(I_dark),
        "t_rise": t_rise,
        "t_duration": t_duration,
        "t_fall": t_fall,
        "t_gap": t_gap,
        "dose_rate_Gy_s": dose_rate_Gy_s,
    }


def transient_cce_vs_dose_rate(
    V_bias=-30.0,
    dose_rates=None,
    t_rise=1e-6,
    t_duration=1e-3,
    t_fall=1e-6,
    dt_min=1e-8,
    dt_max=1e-4,
    epi_thickness_cm=10e-4,
    device_kwargs=None,
    contact=None,
):
    """Sweep dose rates and compute transient CCE at each.

    Builds ONE device (Auger + bias ramp to V_bias), snapshots its dark DC state, and
    for each dose rate restores that state, simulates a single FLASH pulse and extracts
    the transient CCE; the device is deleted once at the end (also on error). Results
    are returned as a pandas DataFrame for easy comparison with steady-state CCE.

    Parameters
    ----------
    V_bias : float
        Reverse bias voltage (V, negative). Default: -30.
    dose_rates : array_like or None
        Dose rates to sweep (Gy/s). Default: [20, 50, 100, 150, 200, 230].
    t_rise : float
        Pulse rise time (s). Default: 1e-6.
    t_duration : float
        Pulse plateau duration (s). Default: 1e-3.
    t_fall : float
        Pulse fall time (s). Default: 1e-6.
    dt_min : float
        Minimum time step (s). Default: 1e-8.
    dt_max : float
        Maximum time step (s). Default: 1e-4.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default: 10e-4 (10 um).
    device_kwargs : dict or None
        Additional keyword arguments for create_dd_device (e.g. N_A, T,
        substrate_thickness_cm, doping_profile, N_D); override the defaults.
    contact : {"cathode", "anode"} or None
        Contact whose (direct) current is integrated. None = "anode" when the device has
        an n+ layer, else "cathode". The direct cathode current on the n+ is the
        majority-carrier current in 1e19 material, dominated by round-off: on the
        pre-2026-10-18 stack (-30 V, 20-230 Gy/s) it scattered 0.983-1.0 (clipped) vs the
        anode's 0.99969-0.99987. Measured 2026-10-18 on the current default stack (local
        mobility, step-end BDF1 estimator; anode): 0.999725 / 0.999868 / 0.999869 /
        0.999870 / 0.999801 / 0.999822 at 20 / 50 / 100 / 150 / 200 / 230 Gy/s (spread
        1.5e-4: the direct anode current carries round-off of ~1e-4 of q*integral(G);
        the balanced steady-state CCE at 20 Gy/s is 0.999809, the direct steady-state
        anode value 0.999729). With the scalar epi mobility (before the A-03 fix) the
        20 Gy/s point read 1.00028 (clipped to 1). Without n+ the cathode is kept (spread 3e-8;
        legacy numbers unchanged). The balanced steady-state form is not used: it has
        no dQ/dt term.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with columns ["dose_rate_Gy_s", "transient_cce",
        "transient_cce_raw"] (the unclipped Q_collected / Q_generated; transient_cce
        is it clipped to [0, 1]); ``df.attrs["contact"]`` is the integrated contact.
    """
    if dose_rates is None:
        dose_rates = np.array([20, 50, 100, 150, 200, 230], dtype=float)
    else:
        dose_rates = np.asarray(dose_rates, dtype=float)

    cce_values = []
    raw_values = []
    dev_kw = dict(
        epi_thickness_cm=epi_thickness_cm,
        **graded_stack_kwargs(),
    )
    dev_kw.update(device_kwargs or {})
    if contact is None:
        # (a missing key = create_sic_device's default, the stack's n+ layer)
        has_nplus = (
            dev_kw.get("nplus_thickness_cm", _STACK["nplus_thickness_cm"]) is not None
        )
        contact = "anode" if has_nplus else "cathode"
    if contact not in ("cathode", "anode"):
        raise ValueError(f"contact must be 'cathode' or 'anode', got {contact!r}")

    # One device for the whole sweep: build, Auger and the bias ramp do not depend on
    # the dose rate. The dark DC state (Potential / Electrons / Holes) at V_bias is
    # snapshotted once and restored before each dose rate's initialize(), whose
    # transient_dc solve also resets the time-integration history; results are
    # bit-identical to the former build-per-dose-rate loop (6 rates: 7.37 -> 2.13 s,
    # lab/engine_box/bench/transient_build_once.py).
    dev_id = uuid.uuid4().hex[:8]
    device_info = create_dd_device(device_name=f"transient_sweep_{dev_id}", **dev_kw)
    device = device_info["device_name"]
    region = device_info["region_name"]
    try:
        add_auger_recombination(device_info)
        # Ramp bias via anode (reverse bias)
        ramp_bias(device_info, V_bias, contact="anode")
        dark_state = {
            name: list(
                devsim.get_node_model_values(device=device, region=region, name=name)
            )
            for name in _DARK_STATE_SOLUTIONS
        }
        x_nodes = np.array(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        for dose_rate in dose_rates:
            for name, values in dark_state.items():
                devsim.set_node_values(
                    device=device, region=region, name=name, values=values
                )
            G_spatial = proton_generation_profile(
                x_nodes, E_MeV=62, dose_rate_Gy_s=dose_rate
            )

            # Zero generation in p+ region
            junction_pos = device_info["junction_pos"]
            G_spatial[x_nodes < junction_pos] = 0.0
            G_spatial = confine_generation_to_epi(device_info, x_nodes, G_spatial)

            # Transient simulation
            solver = TransientSolver(device_info, contact=contact)
            solver.initialize()

            result = solver.simulate_pulse(
                G_spatial,
                t_rise=t_rise,
                t_duration=t_duration,
                t_fall=t_fall,
                dt_min=dt_min,
                dt_max=dt_max,
                dose_rate_Gy_s=dose_rate,
            )

            cce, cce_raw = solver.compute_transient_cce(
                result, G_spatial, x_nodes, return_raw=True
            )
            cce_values.append(cce)
            raw_values.append(cce_raw)

            logger.info(
                f"transient_cce_vs_dose_rate: {dose_rate:.0f} Gy/s "
                f"-> CCE={cce:.4f} (raw {cce_raw:.6f})"
            )
    finally:
        try:
            devsim.delete_device(device=device)
        except Exception:
            pass

    df = pd.DataFrame(
        {
            "dose_rate_Gy_s": dose_rates,
            "transient_cce": cce_values,
            "transient_cce_raw": raw_values,
        }
    )
    df.attrs["contact"] = contact
    return df
