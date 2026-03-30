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

from src.charge_collection_2d import create_2d_dd_device, integrate_over_mesh_2d
from src.charge_collection import add_generation_to_dd
from src.drift_diffusion import extract_contact_current

logger = logging.getLogger(__name__)

# Physical constants
Q_ELECTRON = 1.602e-19  # C, elementary charge


def ion_track_generation_2d(
    device_info, LET_keV_um, x_ion_cm=0.0, track_sigma_cm=1e-4, E_pair_eV=8.4
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
        Electron-hole pair creation energy (eV). Default: 8.4 (4H-SiC).

    Returns
    -------
    generation : ndarray
        Generation density at each mesh node (pairs/cm^3).
    Q_generated_C_per_cm : float
        Total generated charge per unit z-depth (C/cm), computed via
        mesh area integration: integrate(generation) * q_electron.
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

    # Q_generated = integrate(generation * dA) * q  [C/cm for 2D]
    Q_generated_C_per_cm = integrate_over_mesh_2d(device_info, generation) * Q_ELECTRON

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
    contact="cathode",
):
    """Simulate single-particle charge collection transient.

    Uses the generation-pulse method (Approach A from research):
    1. Initialize transient state with zero generation
    2. Inject charge as a short pulse (generation/dt_inject for one timestep)
    3. Zero generation and run BDF1 time-stepping until current returns to
       baseline

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
    contact : str
        Contact for current extraction.

    Returns
    -------
    result : dict
        - "times": time points (s), ndarray
        - "currents": contact current at each time (A/cm for 2D), ndarray
        - "Q_collected": integrated collected charge (C/cm for 2D), float
        - "I_dark": baseline dark current (A/cm), float
        - "I_peak": peak signal current (A/cm), float
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    n_nodes = len(generation_profile)

    # Step 1: Initialize transient state with zero generation
    zeros = np.zeros(n_nodes)
    add_generation_to_dd(device_info, zeros)
    devsim.solve(
        type="transient_dc",
        absolute_error=1e10,
        relative_error=1e-10,
        maximum_iterations=30,
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
        )
    except devsim.error:
        logger.warning("Injection step failed at standard tolerances, retrying relaxed")
        devsim.solve(
            type="transient_bdf1",
            tdelta=dt_inject,
            absolute_error=1e12,
            relative_error=1e-6,
            maximum_iterations=200,
        )

    I_after_inject = extract_contact_current(device_info, contact)
    logger.debug(f"Current after injection: {I_after_inject:.4e} A/cm")

    # Step 4: Zero generation and time-step loop
    add_generation_to_dd(device_info, zeros)

    times = [dt_inject]
    currents = [I_after_inject]
    t = dt_inject

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
                )
            except devsim.error:
                logger.warning(f"BDF1 step failed at t={t:.3e} s, terminating early")
                break

        t += dt
        I = extract_contact_current(device_info, contact)
        times.append(t)
        currents.append(I)

        if abs(I) > abs(I_peak):
            I_peak = I

        # Early termination: signal decayed to 1% of peak
        if abs(I - I_dark) < 0.01 * abs(I_peak - I_dark) and t > 10 * dt_inject:
            logger.debug(f"Early termination at t={t:.3e} s (signal < 1% peak)")
            break

        step_count += 1

    times = np.array(times)
    currents = np.array(currents)

    # Q_collected = integral(|I| - |I_dark|) dt, using absolute values
    Q_collected = float(np.trapezoid(np.abs(currents) - np.abs(I_dark), times))

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
    }


def analyze_current_pulse(times, currents, I_dark):
    """Extract pulse characteristics from transient simulation.

    Parameters
    ----------
    times : array_like
        Time points (s).
    currents : array_like
        Contact current at each time point (A/cm for 2D).
    I_dark : float
        Baseline dark current (A/cm).

    Returns
    -------
    result : dict
        - "Q_collected": integrated charge (C/cm for 2D)
        - "I_peak": peak signal current magnitude (A/cm)
        - "t_peak": time of peak current (s)
        - "t_collection": time to collect 95% of total charge (s)
    """
    times = np.asarray(times)
    currents = np.asarray(currents)

    I_signal = np.abs(currents) - np.abs(I_dark)
    Q_total = float(np.trapezoid(I_signal, times))

    I_peak = float(np.max(I_signal))
    t_peak = float(times[np.argmax(I_signal)])

    # Collection time: when cumulative charge reaches 95%
    if len(times) > 1:
        Q_cumulative = np.cumsum(0.5 * (I_signal[:-1] + I_signal[1:]) * np.diff(times))
        idx_95 = np.searchsorted(Q_cumulative, 0.95 * Q_total)
        t_collection = float(times[min(idx_95 + 1, len(times) - 1)])
    else:
        t_collection = float(times[0])

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

    Returns
    -------
    table : pandas.DataFrame
        Columns: LET_keV_um, Q_generated_fC, Q_collected_fC, CCE,
        t_collection_ns.
    """
    LET_values = np.logspace(np.log10(let_min), np.log10(let_max), n_let_points)

    results = []
    for i, LET in enumerate(LET_values):
        device_info = None
        try:
            device_info = create_2d_dd_device(
                half_width_um=half_width_um, V_bias=V_bias
            )

            generation, Q_gen = ion_track_generation_2d(
                device_info, LET, x_ion_cm=x_ion_cm
            )

            sim_result = simulate_single_particle(device_info, generation)

            Q_col = sim_result["Q_collected"]
            cce = Q_col / Q_gen if Q_gen > 0 else float("nan")

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
