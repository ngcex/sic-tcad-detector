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

import devsim
import numpy as np

from src.charge_collection import add_generation_to_dd
from src.drift_diffusion import extract_contact_current

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

            # Transient solve
            devsim.solve(
                type=self.method,
                absolute_error=1e10,
                relative_error=1e-10,
                maximum_iterations=30,
                tdelta=dt,
                charge_error=1e-2,
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
            Transient charge collection efficiency, clipped to [0, 2].
        """
        times = result["times"]
        currents = result["currents"]
        I_dark = result["I_dark"]
        t_duration = result["t_duration"]

        # Collected charge: integral of (|I(t)| - |I_dark|) dt
        I_signal = np.abs(currents) - np.abs(I_dark)
        Q_collected = np.trapezoid(I_signal, times)

        # Total generation rate per unit area: integral of G(x) dx
        G_total = np.trapezoid(
            np.asarray(G_spatial, dtype=float), np.asarray(x_nodes, dtype=float)
        )

        # Generated charge during plateau: q * G_total * t_duration
        Q_generated = Q_ELECTRON * G_total * t_duration

        if Q_generated <= 0:
            logger.warning("Q_generated <= 0, returning CCE = 0")
            return 0.0

        cce = Q_collected / Q_generated

        # Clip to [0, 2] -- allow slight overshoot for transit-time effects
        cce = float(np.clip(cce, 0.0, 2.0))

        logger.info(
            f"Transient CCE: {cce:.4f} "
            f"(Q_collected={Q_collected:.4e}, Q_generated={Q_generated:.4e})"
        )
        return cce
