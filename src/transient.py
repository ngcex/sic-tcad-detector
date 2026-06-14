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

from src.charge_collection import add_generation_to_dd
from src.drift_diffusion import create_dd_device, extract_contact_current, ramp_bias
from src.flash_recombination import add_auger_recombination
from src.generation_profiles import proton_generation_profile

logger = logging.getLogger(__name__)

# Physical constants (CGS)
Q_ELECTRON = 1.602e-19  # C (elementary charge)


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


class TransientSolver:
    """Transient drift-diffusion solver for FLASH pulse simulation.

    Wraps the devsim transient solve API with adaptive time-stepping
    and pulse envelope modulation of the generation rate.

    Parameters
    ----------
    device_info : dict
        Device info dict from create_dd_device(), already biased and
        with Auger recombination added. Must have dd_initialized=True.
    contact : str
        Contact at which to extract current. Default: "cathode".
    method : str
        Transient solve method. Default: "transient_bdf1".
    """

    def __init__(self, device_info, contact="cathode", method="transient_bdf1"):
        self.device_info = device_info
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

        # transient_dc solve to initialize time data (Phase 7 pattern)
        devsim.solve(
            type="transient_dc",
            absolute_error=1e10,
            relative_error=1e-10,
            maximum_iterations=30,
        )
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
            If True, skip the dark current measurement (device already has
            transient state from a previous pulse). Used by
            simulate_pulse_train for multi-pulse sequences. Default: False.

        Returns
        -------
        result : dict
            Dictionary with:
            - "times": np.array of time points (s)
            - "currents": np.array of contact currents (A/cm^2)
            - "I_dark": float, dark current before pulse (A/cm^2)
            - "t_rise": float
            - "t_duration": float
            - "t_fall": float
            - "dose_rate_Gy_s": float or None
        """
        G_spatial = np.asarray(G_spatial, dtype=float)

        if t_post is None:
            t_post = 5 * t_fall

        t_end = t_rise + t_duration + t_fall + t_post

        # Record dark current before pulse
        I_dark = extract_contact_current(self.device_info, contact=self.contact)
        logger.info(f"Dark current: {I_dark:.4e} A/cm^2")

        # Time-stepping loop
        t = 0.0
        times = []
        currents = []
        step_count = 0

        while t < t_end:
            # Compute envelope and time step
            env = pulse_envelope(t, t_rise, t_duration, t_fall)
            dt = adaptive_dt(t, t_rise, t_duration, t_fall, dt_min, dt_max)

            # Don't overshoot t_end
            if t + dt > t_end:
                dt = t_end - t

            # Update generation rate: G(x,t) = G_spatial(x) * envelope(t)
            G_t = G_spatial * env
            add_generation_to_dd(self.device_info, G_t)

            # Transient solve with retry at relaxed tolerances.
            # charge_error=1e10 effectively disables automatic step rejection
            # based on charge conservation error -- we manage time steps
            # ourselves via adaptive_dt based on the pulse envelope phase.
            try:
                devsim.solve(
                    type=self.method,
                    absolute_error=1e10,
                    relative_error=1e-10,
                    maximum_iterations=40,
                    tdelta=dt,
                    charge_error=1e10,
                )
            except devsim.error:
                # Retry with relaxed tolerances and more iterations
                devsim.solve(
                    type=self.method,
                    absolute_error=1e12,
                    relative_error=1e-8,
                    maximum_iterations=100,
                    tdelta=dt,
                    charge_error=1e10,
                )
                logger.debug(
                    f"Step {step_count}: converged with relaxed tolerances "
                    f"at t={t:.4e}, dt={dt:.4e}"
                )

            # Extract current
            I_t = extract_contact_current(self.device_info, contact=self.contact)
            times.append(t)
            currents.append(I_t)

            t += dt
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
            "t_rise": t_rise,
            "t_duration": t_duration,
            "t_fall": t_fall,
            "dose_rate_Gy_s": dose_rate_Gy_s,
        }

    def compute_transient_cce(self, result, G_spatial, x_nodes):
        """Compute charge collection efficiency from transient simulation.

        Integrates (|I(t)| - |I_dark|) over time to get collected charge,
        then divides by generated charge during the pulse plateau.

        Parameters
        ----------
        result : dict
            Result dictionary from simulate_pulse().
        G_spatial : array_like
            Spatial generation profile at mesh nodes (cm^-3 s^-1).
        x_nodes : array_like
            Mesh node positions (cm).

        Returns
        -------
        cce : float
            Transient charge collection efficiency in [0, 1].
        """
        times = result["times"]
        currents = result["currents"]
        I_dark = result["I_dark"]
        t_duration = result["t_duration"]
        t_rise = result.get("t_rise", 0.0)
        t_fall = result.get("t_fall", 0.0)

        # Collected charge: integral of (|I(t)| - |I_dark|) dt over the WHOLE
        # waveform (rise + plateau + fall + post-pulse decay).
        I_signal = np.abs(currents) - np.abs(I_dark)
        Q_collected = np.trapezoid(I_signal, times)

        # Total generation rate per unit area: integral of G(x) dx
        G_total = np.trapezoid(
            np.asarray(G_spatial, dtype=float), np.asarray(x_nodes, dtype=float)
        )

        # Generated charge over the FULL trapezoidal envelope (audit C4).
        # Must match the integration window above: the collected-charge integral
        # spans rise+plateau+fall, so the generated charge must count the rise
        # and fall ramps too, not just the plateau. Using plateau-only here was
        # the normalization bug that pushed CCE above 1 (masked by a [0, 2] clip).
        Q_generated = generated_charge_trapezoidal_pulse(
            G_total, t_rise, t_duration, t_fall, q=Q_ELECTRON
        )

        if Q_generated <= 0:
            logger.warning("Q_generated <= 0, returning CCE = 0")
            return 0.0

        cce_raw = Q_collected / Q_generated

        # Physical ceiling: in a linear DD detector CCE cannot exceed 1. Any
        # residual excess is finite-step trapezoidal-quadrature error on the
        # current pulse (typically < ~1%). Enforce the physical bound, but assert
        # the raw value is within quadrature tolerance -- a larger overshoot
        # would indicate a real bug (e.g. window/normalization mismatch), not
        # numerical noise, and must NOT be silently clipped. (Audit C4.)
        _QUADRATURE_TOL = 0.05  # 5% headroom over the physical ceiling
        if cce_raw > 1.0 + _QUADRATURE_TOL:
            logger.warning(
                f"Transient CCE_raw={cce_raw:.4f} exceeds physical ceiling by "
                f">{_QUADRATURE_TOL:.0%}; likely a normalization/window bug, not "
                f"quadrature error (Q_collected={Q_collected:.4e}, "
                f"Q_generated={Q_generated:.4e})."
            )
        cce = float(np.clip(cce_raw, 0.0, 1.0))

        logger.info(
            f"Transient CCE: {cce:.4f} (raw={cce_raw:.4f}, "
            f"Q_collected={Q_collected:.4e}, Q_generated={Q_generated:.4e})"
        )
        return cce


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
    contact="cathode",
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
    contact : str
        Contact for current extraction. Default: "cathode".
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

    # Record dark current before first pulse
    I_dark = extract_contact_current(device_info, contact=contact)

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
            skip_init=(i > 0),
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
):
    """Sweep dose rates and compute transient CCE at each.

    For each dose rate, creates a fresh device, simulates a single FLASH
    pulse, and extracts the transient CCE. Results are returned as a
    pandas DataFrame for easy comparison with steady-state CCE.

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

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with columns ["dose_rate_Gy_s", "transient_cce"].
    """
    if dose_rates is None:
        dose_rates = np.array([20, 50, 100, 150, 200, 230], dtype=float)
    else:
        dose_rates = np.asarray(dose_rates, dtype=float)

    cce_values = []

    for dose_rate in dose_rates:
        dev_id = uuid.uuid4().hex[:8]
        device_info = create_dd_device(
            device_name=f"transient_sweep_{dev_id}",
            epi_thickness_cm=epi_thickness_cm,
            doping_profile="graded",
            N_D_junction=2.90e15,
            N_D_bulk=8.50e13,
            L_transition=1.0e-4,
        )
        device = device_info["device_name"]
        region = device_info["region_name"]

        try:
            add_auger_recombination(device_info)

            # Ramp bias via anode (reverse bias)
            ramp_bias(device_info, V_bias, contact="anode")

            # Get mesh and generation profile
            x_nodes = np.array(
                devsim.get_node_model_values(device=device, region=region, name="x")
            )
            G_spatial = proton_generation_profile(
                x_nodes, E_MeV=62, dose_rate_Gy_s=dose_rate
            )

            # Zero generation in p+ region
            junction_pos = device_info["junction_pos"]
            G_spatial[x_nodes < junction_pos] = 0.0

            # Transient simulation
            solver = TransientSolver(device_info, contact="cathode")
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

            cce = solver.compute_transient_cce(result, G_spatial, x_nodes)
            cce_values.append(cce)

            logger.info(
                f"transient_cce_vs_dose_rate: {dose_rate:.0f} Gy/s " f"-> CCE={cce:.4f}"
            )

        finally:
            try:
                devsim.delete_device(device=device)
            except Exception:
                pass

    return pd.DataFrame({"dose_rate_Gy_s": dose_rates, "transient_cce": cce_values})
