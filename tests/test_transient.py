"""Tests for transient FLASH pulse simulation.

Validates:
- Trapezoidal pulse envelope shape at key time points
- Adaptive time-step selection during pulse phases
- Single-pulse simulation produces physically reasonable I(t)
- Transient CCE converges toward steady-state value
"""

import uuid

import devsim
import numpy as np
import pytest

from src.drift_diffusion import create_dd_device, ramp_bias, extract_contact_current
from src.charge_collection import add_generation_to_dd
from src.flash_recombination import add_auger_recombination, cce_vs_dose_rate
from src.generation_profiles import proton_generation_profile
from src.transient import (
    TransientSolver,
    pulse_envelope,
    adaptive_dt,
    generated_charge_trapezoidal_pulse,
)


# ---------------------------------------------------------------------------
# Pulse envelope unit tests
# ---------------------------------------------------------------------------


def test_pulse_envelope_stages():
    """Verify pulse_envelope returns correct values at each stage."""
    t_rise = 1e-6
    t_duration = 1e-3
    t_fall = 1e-6

    # Before pulse
    assert pulse_envelope(0.0, t_rise, t_duration, t_fall) == 0.0

    # Midpoint of rise
    val_mid_rise = pulse_envelope(t_rise / 2, t_rise, t_duration, t_fall)
    assert abs(val_mid_rise - 0.5) < 1e-10

    # Top of rise
    val_top = pulse_envelope(t_rise, t_rise, t_duration, t_fall)
    assert abs(val_top - 1.0) < 1e-10

    # During plateau (midpoint)
    t_mid_plateau = t_rise + t_duration / 2
    val_plateau = pulse_envelope(t_mid_plateau, t_rise, t_duration, t_fall)
    assert abs(val_plateau - 1.0) < 1e-10

    # Midpoint of fall
    t_mid_fall = t_rise + t_duration + t_fall / 2
    val_mid_fall = pulse_envelope(t_mid_fall, t_rise, t_duration, t_fall)
    assert abs(val_mid_fall - 0.5) < 1e-10

    # After fall
    t_after = t_rise + t_duration + t_fall + 1e-6
    val_after = pulse_envelope(t_after, t_rise, t_duration, t_fall)
    assert val_after == 0.0


def test_pulse_envelope_edges():
    """Verify pulse_envelope at exact transition boundaries."""
    t_rise = 1e-6
    t_duration = 1e-3
    t_fall = 1e-6

    # Start
    assert pulse_envelope(0.0, t_rise, t_duration, t_fall) == 0.0

    # End of rise = start of plateau
    assert abs(pulse_envelope(t_rise, t_rise, t_duration, t_fall) - 1.0) < 1e-10

    # End of plateau = start of fall
    t_plateau_end = t_rise + t_duration
    assert abs(pulse_envelope(t_plateau_end, t_rise, t_duration, t_fall) - 1.0) < 1e-10

    # End of fall
    t_fall_end = t_rise + t_duration + t_fall
    assert abs(pulse_envelope(t_fall_end, t_rise, t_duration, t_fall) - 0.0) < 1e-10


# ---------------------------------------------------------------------------
# Adaptive dt unit tests
# ---------------------------------------------------------------------------


def test_adaptive_dt_transitions():
    """Verify adaptive_dt returns small dt during transitions, large during plateau."""
    t_rise = 1e-6
    t_duration = 1e-3
    t_fall = 1e-6
    dt_min = 1e-8
    dt_max = 1e-4

    # During rise: should use t_rise/10 = 1e-7
    dt_rise = adaptive_dt(t_rise / 2, t_rise, t_duration, t_fall, dt_min, dt_max)
    assert dt_rise == pytest.approx(t_rise / 10, rel=1e-10)

    # During plateau: should use dt_max
    t_mid_plateau = t_rise + t_duration / 2
    dt_plateau = adaptive_dt(t_mid_plateau, t_rise, t_duration, t_fall, dt_min, dt_max)
    assert dt_plateau == dt_max

    # During fall: should use t_fall/10 = 1e-7
    t_mid_fall = t_rise + t_duration + t_fall / 2
    dt_fall = adaptive_dt(t_mid_fall, t_rise, t_duration, t_fall, dt_min, dt_max)
    assert dt_fall == pytest.approx(t_fall / 10, rel=1e-10)

    # Post-pulse: should use dt_max
    t_post = t_rise + t_duration + t_fall + 1e-3
    dt_post = adaptive_dt(t_post, t_rise, t_duration, t_fall, dt_min, dt_max)
    assert dt_post == dt_max


def test_adaptive_dt_bounds():
    """Verify dt is always within [dt_min, dt_max] bounds."""
    t_rise = 1e-6
    t_duration = 1e-3
    t_fall = 1e-6
    dt_min = 1e-8
    dt_max = 1e-4

    # Test at many time points spanning the full pulse
    t_end = t_rise + t_duration + t_fall + 1e-3
    test_times = np.linspace(0, t_end, 200)

    for t in test_times:
        dt = adaptive_dt(t, t_rise, t_duration, t_fall, dt_min, dt_max)
        assert dt >= dt_min, f"dt={dt} < dt_min={dt_min} at t={t}"
        assert dt <= dt_max, f"dt={dt} > dt_max={dt_max} at t={t}"


# ---------------------------------------------------------------------------
# Generated-charge bookkeeping (audit C4: CCE>1 normalization)
# ---------------------------------------------------------------------------


def test_generated_charge_counts_full_envelope():
    """Generated charge must integrate the WHOLE trapezoidal envelope.

    Audit C4: the old code normalized CCE by q*G_total*t_duration (plateau
    only), omitting the charge generated during rise and fall. The correct
    generated charge is q*G_total*(t_rise/2 + t_duration + t_fall/2), the area
    under the trapezoidal envelope. Omitting rise/fall makes the denominator too
    small and inflates CCE above 1 (the bug the [0,2] clip was masking).
    """
    q = 1.602e-19
    G_total = 1.0e15  # cm^-3 s^-1 integrated over x -> cm^-2 s^-1 (per-area rate)
    t_rise, t_duration, t_fall = 2e-6, 1e-3, 3e-6

    Q = generated_charge_trapezoidal_pulse(G_total, t_rise, t_duration, t_fall, q=q)
    expected = q * G_total * (t_rise / 2 + t_duration + t_fall / 2)
    assert Q == pytest.approx(expected, rel=1e-12)


def test_generated_charge_exceeds_plateau_only():
    """Full-envelope charge must be strictly larger than the plateau-only value."""
    q = 1.602e-19
    G_total = 1.0e15
    t_rise, t_duration, t_fall = 5e-4, 1e-3, 5e-4  # rise/fall comparable to plateau
    Q_full = generated_charge_trapezoidal_pulse(
        G_total, t_rise, t_duration, t_fall, q=q
    )
    Q_plateau_only = q * G_total * t_duration
    # rise+fall add (t_rise+t_fall)/2 = 5e-4 s of generation on top of 1e-3 plateau
    assert Q_full > Q_plateau_only
    assert Q_full == pytest.approx(Q_plateau_only * 1.5, rel=1e-12)


def test_generated_charge_reduces_to_plateau_when_no_ramps():
    """With zero rise/fall the envelope is a pure rectangle (plateau only)."""
    q = 1.602e-19
    G_total = 2.0e15
    Q = generated_charge_trapezoidal_pulse(G_total, 0.0, 1e-3, 0.0, q=q)
    assert Q == pytest.approx(q * G_total * 1e-3, rel=1e-12)


# ---------------------------------------------------------------------------
# Integration tests (slow -- require devsim device simulation)
# ---------------------------------------------------------------------------


def _make_transient_device(name_prefix="test_transient"):
    """Create a DD device with Auger, biased for transient simulation."""
    dev_id = uuid.uuid4().hex[:8]
    device_info = create_dd_device(
        device_name=f"{name_prefix}_{dev_id}",
        doping_profile="graded",
        N_D_junction=2.90e15,
        N_D_bulk=8.50e13,
        L_transition=1.0e-4,
    )
    add_auger_recombination(device_info)
    ramp_bias(device_info, -30.0, contact="anode")
    return device_info


@pytest.mark.slow
def test_single_pulse_simulation():
    """Integration test: single FLASH pulse produces valid I(t) waveform."""
    device_info = _make_transient_device()
    device = device_info["device_name"]
    region = device_info["region_name"]

    try:
        # Get mesh and generation profile
        x_nodes = np.array(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        G_spatial = proton_generation_profile(x_nodes, E_MeV=62, dose_rate_Gy_s=100.0)

        # Zero generation in p+ region
        junction_pos = device_info["junction_pos"]
        G_spatial[x_nodes < junction_pos] = 0.0

        # Create solver, initialize, simulate
        solver = TransientSolver(device_info, contact="cathode")
        solver.initialize()

        result = solver.simulate_pulse(
            G_spatial,
            t_rise=1e-6,
            t_duration=1e-3,
            t_fall=1e-6,
            dt_min=1e-8,
            dt_max=5e-5,
        )

        # Validate result structure
        assert "times" in result
        assert "currents" in result
        assert "I_dark" in result

        times = result["times"]
        currents = result["currents"]
        I_dark = result["I_dark"]

        # Enough time steps taken
        assert len(times) > 20, f"Only {len(times)} steps taken"

        # Current increases during pulse (signal above dark current)
        assert np.max(np.abs(currents)) > np.abs(I_dark), (
            f"Max current {np.max(np.abs(currents)):.4e} not above "
            f"dark current {np.abs(I_dark):.4e}"
        )

        # Compute transient CCE -- should be physically reasonable
        cce = solver.compute_transient_cce(result, G_spatial, x_nodes)
        assert 0.5 < cce < 1.5, f"Transient CCE={cce:.4f} outside [0.5, 1.5]"

    finally:
        try:
            devsim.delete_device(device=device)
        except Exception:
            pass


@pytest.mark.slow
def test_transient_cce_matches_steady_state():
    """Validation: transient CCE converges toward steady-state value."""
    device_info = _make_transient_device("test_cce_match")
    device = device_info["device_name"]
    region = device_info["region_name"]

    try:
        # Get mesh and generation profile
        x_nodes = np.array(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        G_spatial = proton_generation_profile(x_nodes, E_MeV=62, dose_rate_Gy_s=100.0)

        junction_pos = device_info["junction_pos"]
        G_spatial[x_nodes < junction_pos] = 0.0

        # Transient CCE with longer pulse for better convergence
        solver = TransientSolver(device_info, contact="cathode")
        solver.initialize()

        result = solver.simulate_pulse(
            G_spatial,
            t_rise=1e-6,
            t_duration=5e-3,  # longer pulse
            t_fall=1e-6,
            dt_min=1e-8,
            dt_max=1e-4,
        )

        transient_cce = solver.compute_transient_cce(result, G_spatial, x_nodes)

    finally:
        try:
            devsim.delete_device(device=device)
        except Exception:
            pass

    # Steady-state CCE at same conditions
    ss_result = cce_vs_dose_rate(
        dose_rates_Gy_s=[100.0],
        V_bias=-30.0,
        E_MeV=62,
    )
    steady_state_cce = ss_result["cce_values"][0]

    # Physical ceiling (audit C4): CCE must not exceed 1 once the generated
    # charge is normalized over the full pulse envelope. Guards against a
    # re-regression of the plateau-only normalization bug.
    assert 0.0 <= transient_cce <= 1.0, (
        f"Transient CCE ({transient_cce:.4f}) outside physical [0, 1] range "
        f"-- normalization/window regression (C4)?"
    )

    # Transient should be within 20% of steady-state
    deviation = abs(transient_cce - steady_state_cce)
    assert deviation < 0.2, (
        f"Transient CCE ({transient_cce:.4f}) differs from steady-state "
        f"({steady_state_cce:.4f}) by {deviation:.4f} (> 0.2 tolerance)"
    )
