"""Design optimization utilities for 4H-SiC microdosimeter feasibility study.

Provides:
- Parametric sweep of SV geometry (half-width, epi thickness, doping, bias)
  with CCE uniformity ranking
- Noise floor estimation from dark current (shot noise limit)
- Multi-criteria structure comparison scoring (planar, mesa, 3D, etc.)
- Dark current extraction helper for 2D devices

All units follow project conventions:
    - Lengths: cm (internal), um (user-facing parameters)
    - Currents: A/cm (2D) or A (total)
    - Energies: eV, keV
    - Lineal energy: keV/um

References:
    - Phase 20: charge_collection_2d.py CCE methodology
    - Phase 14: dark_current.py TAT model
    - Phase 23: microdosimetry.py mean_chord_length
    - ICRU Report 36: microdosimetric definitions
"""

import itertools
import logging
import uuid
import warnings

import devsim
import numpy as np
import pandas as pd

from src.analytical import built_in_potential, full_depletion_voltage_graded
from src.charge_collection_2d import create_2d_dd_device, cce_lateral_scan
from src.dark_current import setup_tat_model
from src.devsim_reset import reset_devsim_fully
from src.drift_diffusion import extract_contact_current
from src.microdosimetry import mean_chord_length
from src.sic_material import SiC4H_Parameters

# AUDIT Mj-3 (v5): minimum CCE for a config to be a trustworthy microdosimeter.
# Below this, edge/center uniformity is meaningless (both-low CCE looks uniform).
CCE_FLOOR = 0.95

logger = logging.getLogger(__name__)

# Physical constants
Q = 1.602e-19  # C, elementary charge
# E-h pair creation energy: single source of truth in SiC4H_Parameters (7.8 eV,
# measured 4H-SiC W-value). Audit M6: was hardcoded 8.4 eV here, ~7.7% high.
E_PAIR_EV = (
    SiC4H_Parameters().E_pair_eV
)  # eV, electron-hole pair creation energy in 4H-SiC


def microdosimetric_sweep(
    half_widths_um=(25, 50, 100, 150),
    epi_thicknesses_cm=(5e-4, 10e-4, 20e-4),
    N_D_bulks=(5e13, 8.5e13, 5e14),
    V_biases=(20, 40, 60),
    n_lateral_points=10,
    max_configs=None,
):
    """Sweep SV geometry and doping parameters, ranking by CCE uniformity.

    Creates a full combinatorial grid of (half_width, epi_thickness, N_D_bulk,
    V_bias) and for each configuration runs a lateral CCE scan to quantify
    edge-to-center uniformity.

    Parameters
    ----------
    half_widths_um : tuple of float
        SV half-widths to sweep (micrometers).
    epi_thicknesses_cm : tuple of float
        Epitaxial layer thicknesses to sweep (cm).
    N_D_bulks : tuple of float
        Bulk donor concentrations to sweep (cm^-3).
    V_biases : tuple of float
        Reverse bias voltages to sweep (V).
    n_lateral_points : int
        Number of lateral positions per CCE scan.
    max_configs : int or None
        If set, truncate the parameter grid to this many configurations
        (with a warning).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: half_width_um, epi_thickness_cm, N_D_bulk,
        V_bias, center_cce, edge_cce, edge_center_ratio, cce_std, plus the
        Mj-3 validity gate: V_fd (full-depletion voltage from the graded
        profile), is_fully_depleted (|V_bias| >= V_fd), passes_cce_floor
        (both center & edge CCE >= CCE_FLOOR=0.95), and is_valid (both).
        Sorted so all valid configs rank above all invalid ones, then by
        edge_center_ratio descending within each group. Invalid (partially
        depleted or low-CCE) configs are RETAINED and flagged, not dropped --
        do not trust the uniformity ratio of a row with is_valid=False.
    """
    grid = list(
        itertools.product(half_widths_um, epi_thicknesses_cm, N_D_bulks, V_biases)
    )
    total = len(grid)

    if max_configs is not None and total > max_configs:
        warnings.warn(
            f"Parameter grid has {total} configurations, truncating to {max_configs}. "
            f"Reduce sweep ranges for exhaustive coverage.",
            UserWarning,
            stacklevel=2,
        )
        grid = grid[:max_configs]
        total = len(grid)

    records = []
    for i, (hw, epi, nd, vb) in enumerate(grid):
        logger.info(f"Sweep {i + 1}/{total}: hw={hw}, epi={epi}, nd={nd:.1e}, vb={vb}")
        device_info = None
        try:
            device_info = create_2d_dd_device(
                half_width_um=hw,
                V_bias=vb,
                epi_thickness_cm=epi,
                N_D_bulk=nd,
            )
            scan = cce_lateral_scan(device_info, n_points=n_lateral_points)
            cce_vals = np.array(scan["cce_values"])

            # AUDIT Mj-3: full-depletion gate. A partially-depleted device with
            # both-low CCE shows a deceptively uniform edge/center ratio and
            # would rank artificially high. Compute V_fd from the device's
            # actual graded profile (NOT the uniform N_D_bulk estimate, which is
            # anti-conservative) and require both full depletion AND an absolute
            # CCE floor before the uniformity ratio is trustworthy.
            V_bi = built_in_potential(
                device_info["N_A_ionized"],
                device_info["N_D_junction"],
                device_info["n_i"],
            )
            V_fd = full_depletion_voltage_graded(
                epi_thickness=device_info["epi_thickness_cm"],
                N_D_bulk=device_info["N_D_bulk"],
                N_D_junction=device_info["N_D_junction"],
                L_transition=device_info["L_transition"],
                V_bi=V_bi,
            )
            center_cce = float(cce_vals[0])
            edge_cce = float(cce_vals[-1])
            is_fully_depleted = abs(vb) >= V_fd
            passes_cce_floor = center_cce >= CCE_FLOOR and edge_cce >= CCE_FLOOR
            records.append(
                {
                    "half_width_um": hw,
                    "epi_thickness_cm": epi,
                    "N_D_bulk": nd,
                    "V_bias": vb,
                    "center_cce": center_cce,
                    "edge_cce": edge_cce,
                    "edge_center_ratio": float(scan["edge_to_center_ratio"]),
                    "cce_std": float(np.std(cce_vals)),
                    "V_fd": float(V_fd),
                    "is_fully_depleted": bool(is_fully_depleted),
                    "passes_cce_floor": bool(passes_cce_floor),
                    "is_valid": bool(is_fully_depleted and passes_cce_floor),
                }
            )
        except Exception as e:
            logger.warning(
                f"Sweep {i + 1}/{total} failed (hw={hw}, epi={epi}, "
                f"nd={nd:.1e}, vb={vb}): {e}"
            )
            records.append(
                {
                    "half_width_um": hw,
                    "epi_thickness_cm": epi,
                    "N_D_bulk": nd,
                    "V_bias": vb,
                    "center_cce": np.nan,
                    "edge_cce": np.nan,
                    "edge_center_ratio": np.nan,
                    "cce_std": np.nan,
                    "V_fd": np.nan,
                    "is_fully_depleted": False,
                    "passes_cce_floor": False,
                    "is_valid": False,
                }
            )
            # Reset devsim global state after failure to prevent
            # contamination of subsequent configurations (a failed
            # ramp_bias can leave solver state that causes all later
            # equilibrium solves to diverge). reset_devsim_fully also clears
            # any cylindrical-axis globals (P03/P30) and preserves the solver.
            reset_devsim_fully(preserve_solver=True)
        finally:
            if device_info is not None:
                try:
                    devsim.delete_device(device=device_info["device_name"])
                except Exception:
                    pass

    return _rank_sweep_results(pd.DataFrame(records))


def _rank_sweep_results(df):
    """Rank microdosimetric-sweep results with the Mj-3 validity gate.

    Valid configs (fully-depleted AND above the CCE floor) rank above all
    invalid ones, then by uniformity (edge_center_ratio) within each group.
    This guarantees a deceptive partially-depleted config can never outrank a
    genuinely uniform fully-depleted one. Invalid rows are retained (flagged
    via is_valid=False), not dropped. Pure function for testability.
    """
    df = df.sort_values(
        ["is_valid", "edge_center_ratio"],
        ascending=[False, False],
        na_position="last",
    )
    return df.reset_index(drop=True)


def estimate_noise_floor(
    I_dark_A,
    t_shaping_s=1e-6,
    sv_thickness_um=10.0,
    sv_width_um=100.0,
    sv_depth_um=None,
    k_sigma=3.0,
):
    """Estimate detector-intrinsic noise floor from dark current.

    Computes the minimum detectable charge, energy, and lineal energy
    from shot noise on the dark current. This represents the fundamental
    detector-intrinsic limit; real systems will have higher noise from
    readout electronics, electromagnetic interference, and other sources.

    Parameters
    ----------
    I_dark_A : float
        Dark current in Amperes (total, not density).
    t_shaping_s : float
        Shaping time of the readout electronics (seconds).
    sv_thickness_um : float
        Sensitive volume thickness (micrometers).
    sv_width_um : float
        Sensitive volume width (micrometers).
    sv_depth_um : float or None
        Sensitive volume depth (micrometers). If None, uses slab
        approximation for mean chord length.
    k_sigma : float
        Number of standard deviations for detection threshold
        (default 3 for 99.7% confidence).

    Returns
    -------
    dict
        Dictionary with keys:
        - I_dark_A: input dark current
        - sigma_shot_C: shot noise RMS charge (Coulombs)
        - Q_min_fC: minimum detectable charge (femtocoulombs)
        - E_min_keV: minimum detectable energy (keV)
        - y_min_keV_um: minimum detectable lineal energy (keV/um)
        - l_bar_um: mean chord length (um)
        - t_shaping_s: shaping time used
        - k_sigma: detection threshold used
    """
    sigma_shot = np.sqrt(2.0 * Q * abs(I_dark_A) * t_shaping_s)
    Q_min = k_sigma * sigma_shot

    # Energy threshold from pair creation energy
    E_min_eV = Q_min * E_PAIR_EV / Q
    E_min_keV = E_min_eV / 1e3

    # Mean chord length for lineal energy conversion
    l_bar = mean_chord_length(sv_thickness_um, sv_width_um, sv_depth_um)

    # Minimum detectable lineal energy
    y_min_keV_um = E_min_keV / l_bar

    return {
        "I_dark_A": I_dark_A,
        "sigma_shot_C": float(sigma_shot),
        "Q_min_fC": float(Q_min * 1e15),
        "E_min_keV": float(E_min_keV),
        "y_min_keV_um": float(y_min_keV_um),
        "l_bar_um": float(l_bar),
        "t_shaping_s": t_shaping_s,
        "k_sigma": k_sigma,
    }


def score_structures(
    metrics_dict,
    weights=None,
):
    """Score and rank alternative microdosimeter structures.

    Applies multi-criteria weighted scoring with min-max normalization
    across structures. Metrics are classified as "higher is better"
    (cce_uniformity) or "lower is better" (noise_floor, spectral_resolution,
    fabrication_complexity).

    Parameters
    ----------
    metrics_dict : dict
        {structure_name: {metric_name: value}} where metrics are:
        - cce_uniformity: edge/center CCE ratio (higher = better)
        - noise_floor: y_min in keV/um (lower = better)
        - spectral_resolution: y_D/y_F ratio or d(y) peak std
          (lower = better)
        - fabrication_complexity: 1-4 scale (lower = better)
    weights : dict or None
        {metric_name: weight}. Default: cce_uniformity=0.30,
        noise_floor=0.20, spectral_resolution=0.20,
        fabrication_complexity=0.30. Weights must sum to 1.0.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: structure, each raw metric, each
        normalized metric (suffixed _norm), weighted_score.
        Sorted by weighted_score descending (best first).
    """
    if weights is None:
        weights = {
            "cce_uniformity": 0.30,
            "noise_floor": 0.20,
            "spectral_resolution": 0.20,
            "fabrication_complexity": 0.30,
        }

    # Metrics where higher is better
    higher_is_better = {"cce_uniformity"}

    # Build raw DataFrame
    structures = list(metrics_dict.keys())
    metric_names = list(weights.keys())

    rows = []
    for name in structures:
        row = {"structure": name}
        for m in metric_names:
            row[m] = metrics_dict[name].get(m, 0.0)
        rows.append(row)

    df = pd.DataFrame(rows)

    # Normalize each metric to [0, 1]
    for m in metric_names:
        vals = df[m].values.astype(float)
        v_min = vals.min()
        v_max = vals.max()

        if v_max == v_min:
            # All equal -- assign 0.5
            df[f"{m}_norm"] = 0.5
        else:
            normalized = (vals - v_min) / (v_max - v_min)
            if m not in higher_is_better:
                # Lower is better: invert
                normalized = 1.0 - normalized
            df[f"{m}_norm"] = normalized

    # Compute weighted score
    df["weighted_score"] = sum(weights[m] * df[f"{m}_norm"] for m in metric_names)

    df = df.sort_values("weighted_score", ascending=False).reset_index(drop=True)
    return df


def get_dark_current_2d(half_width_um=50.0, V_bias=50.0, **device_kwargs):
    """Extract total dark current from a 2D DD device with TAT model.

    Creates a 2D drift-diffusion device, applies the trap-assisted tunneling
    model, re-solves, and extracts the dark current from the cathode contact.
    The 2D current (A/cm) is multiplied by the SV depth to get total current
    in Amperes.

    Parameters
    ----------
    half_width_um : float
        SV half-width (micrometers).
    V_bias : float
        Reverse bias voltage (V, positive).
    **device_kwargs
        Additional arguments for create_2d_dd_device (e.g.,
        epi_thickness_cm, N_D_bulk).

    Returns
    -------
    float
        Dark current in Amperes.
    """
    device_info = None
    try:
        device_info = create_2d_dd_device(
            half_width_um=half_width_um,
            V_bias=V_bias,
            **device_kwargs,
        )

        # Apply TAT model and re-solve
        setup_tat_model(device_info)
        devsim.solve(
            type="dc",
            absolute_error=1e10,
            relative_error=1e-10,
            maximum_iterations=40,
        )

        # Extract 2D current (A/cm)
        I_2d = extract_contact_current(device_info, contact="cathode")

        # Convert to total current: multiply by SV depth
        # Default SV depth = 2 * half_width (square cross-section)
        sv_depth_cm = 2.0 * half_width_um * 1e-4
        I_total = I_2d * sv_depth_cm

        logger.info(
            f"Dark current 2D: I_2d={I_2d:.3e} A/cm, "
            f"sv_depth={sv_depth_cm * 1e4:.0f} um, I_total={I_total:.3e} A"
        )

        return float(I_total)

    finally:
        if device_info is not None:
            try:
                devsim.delete_device(device=device_info["device_name"])
            except Exception:
                pass
