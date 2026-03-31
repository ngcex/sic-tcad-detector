"""Microdosimetric spectra computation for 4H-SiC microdosimeter.

Transforms pulse height distributions (energy deposited per event) into
standard microdosimetric observables following ICRU Report 36 conventions:

- Lineal energy: y = epsilon / l_bar (keV/um)
- Frequency distribution f(y): normalized so integral f(y)*dy = 1
- Dose distribution d(y): d(y) = y * f(y) / y_F
- Frequency-mean lineal energy: y_F = integral(y * f(y) * dy)
- Dose-mean lineal energy: y_D = integral(y * d(y) * dy)
- Tissue-equivalence correction via kappa = S_tissue / S_SiC

Also provides publication-quality plotting functions for y*d(y) vs log(y)
and y*f(y) vs log(y) spectra.

References:
    - ICRU Report 36 (1983): Microdosimetry definitions and conventions
    - Cauchy theorem: l_bar = 4V/S for convex bodies under isotropic irradiation
    - Phase 22: mc_coupling.py (process_mc_ensemble output)
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Mean chord length
# ---------------------------------------------------------------------------


def mean_chord_length(sv_thickness_um, sv_width_um=None, sv_depth_um=None):
    """Compute mean chord length for the sensitive volume geometry.

    For a rectangular parallelepiped: l_bar = 4V/S (Cauchy's theorem,
    valid for convex bodies under isotropic irradiation).
    If only thickness given: slab approximation l_bar = 2*thickness
    (valid when lateral dimensions >> thickness).

    Parameters
    ----------
    sv_thickness_um : float
        SV thickness in micrometers.
    sv_width_um : float, optional
        SV width in micrometers.
    sv_depth_um : float, optional
        SV depth in micrometers.

    Returns
    -------
    float
        Mean chord length in micrometers.
    """
    t = sv_thickness_um
    if sv_width_um is not None and sv_depth_um is not None:
        w, d = sv_width_um, sv_depth_um
        V = w * t * d
        S = 2.0 * (w * t + w * d + t * d)
        l_bar = 4.0 * V / S
        logger.debug(
            "3D mean chord length: %.3f um (w=%.1f, t=%.1f, d=%.1f)",
            l_bar,
            w,
            t,
            d,
        )
        return l_bar
    else:
        l_bar = 2.0 * t
        logger.debug("Slab approximation mean chord length: %.3f um (t=%.1f)", l_bar, t)
        return l_bar


# ---------------------------------------------------------------------------
# Log-spaced y-bins (ICRU 36 convention)
# ---------------------------------------------------------------------------


def make_y_bins(y_min=0.01, y_max=1e4, bins_per_decade=50):
    """Create log-spaced lineal energy bins per ICRU Report 36 convention.

    Parameters
    ----------
    y_min : float
        Minimum lineal energy in keV/um.
    y_max : float
        Maximum lineal energy in keV/um.
    bins_per_decade : int
        Number of bins per decade of lineal energy.

    Returns
    -------
    dict
        Keys: bin_edges, bin_centers, bin_widths (all ndarray).
    """
    n_decades = np.log10(y_max) - np.log10(y_min)
    n_bins = int(round(n_decades * bins_per_decade))
    bin_edges = np.logspace(np.log10(y_min), np.log10(y_max), n_bins + 1)
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])  # geometric mean
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    return {
        "bin_edges": bin_edges,
        "bin_centers": bin_centers,
        "bin_widths": bin_widths,
    }


# ---------------------------------------------------------------------------
# Lineal energy spectrum: f(y), d(y), y_F, y_D
# ---------------------------------------------------------------------------


def lineal_energy_spectrum(
    collected_energies_keV,
    l_bar_um,
    y_min=0.01,
    y_max=1e4,
    bins_per_decade=50,
):
    """Compute lineal energy spectrum from collected energies.

    Converts energy depositions to lineal energy y = epsilon / l_bar,
    then computes the frequency distribution f(y) and dose distribution
    d(y) on log-spaced bins following ICRU Report 36 conventions.

    Parameters
    ----------
    collected_energies_keV : array_like
        Collected energy per event in keV.
    l_bar_um : float
        Mean chord length in micrometers.
    y_min : float
        Minimum lineal energy bin edge in keV/um.
    y_max : float
        Maximum lineal energy bin edge in keV/um.
    bins_per_decade : int
        Number of bins per decade.

    Returns
    -------
    dict
        bin_edges, bin_centers, bin_widths : ndarray
        f_y, d_y : ndarray -- frequency and dose distributions
        y_F, y_D : float -- frequency-mean and dose-mean lineal energy
        n_events : int -- number of events used
        y_values : ndarray -- raw lineal energy per event
    """
    energies = np.asarray(collected_energies_keV, dtype=np.float64)
    y_values = energies / l_bar_um  # keV / um

    bins = make_y_bins(y_min=y_min, y_max=y_max, bins_per_decade=bins_per_decade)
    bin_edges = bins["bin_edges"]
    bin_centers = bins["bin_centers"]
    bin_widths = bins["bin_widths"]

    # Histogram into log-spaced bins
    counts, _ = np.histogram(y_values, bins=bin_edges)
    n_total = len(y_values)

    # f(y): normalized frequency distribution, integral f(y)*dy = 1
    f_y = np.zeros_like(bin_centers)
    nonzero = bin_widths > 0
    f_y[nonzero] = counts[nonzero] / (n_total * bin_widths[nonzero])

    # Validate normalization
    integral_f = np.sum(f_y * bin_widths)
    if abs(integral_f - 1.0) > 0.01:
        logger.warning(
            "f(y) normalization: integral = %.4f (expected 1.0, "
            "deviation %.2f%%). Some events may fall outside [%.2e, %.2e] keV/um.",
            integral_f,
            abs(integral_f - 1.0) * 100,
            y_min,
            y_max,
        )

    # y_F: frequency-mean lineal energy
    y_F = np.sum(bin_centers * f_y * bin_widths)

    # d(y): dose distribution
    d_y = np.zeros_like(bin_centers)
    if y_F > 0:
        d_y = bin_centers * f_y / y_F

    # y_D: dose-mean lineal energy
    y_D = np.sum(bin_centers * d_y * bin_widths)

    # Validate y_D >= y_F (Jensen's inequality)
    if y_D < y_F and y_F > 0:
        logger.warning(
            "y_D (%.4f) < y_F (%.4f): possible normalization issue. "
            "Jensen's inequality requires y_D >= y_F.",
            y_D,
            y_F,
        )

    logger.info(
        "Lineal energy spectrum: %d events, y_F=%.3f keV/um, y_D=%.3f keV/um",
        n_total,
        y_F,
        y_D,
    )

    return {
        "bin_edges": bin_edges,
        "bin_centers": bin_centers,
        "bin_widths": bin_widths,
        "f_y": f_y,
        "d_y": d_y,
        "y_F": y_F,
        "y_D": y_D,
        "n_events": n_total,
        "y_values": y_values,
    }


# ---------------------------------------------------------------------------
# Microdosimetric means (standalone convenience)
# ---------------------------------------------------------------------------


def compute_microdosimetric_means(bin_centers, f_y, bin_widths):
    """Compute y_F and y_D from a pre-computed f(y) distribution.

    Parameters
    ----------
    bin_centers : ndarray
        Bin center values (keV/um).
    f_y : ndarray
        Frequency distribution values.
    bin_widths : ndarray
        Bin widths (dy).

    Returns
    -------
    dict
        y_F : float -- frequency-mean lineal energy
        y_D : float -- dose-mean lineal energy
    """
    y_F = np.sum(bin_centers * f_y * bin_widths)
    if y_F > 0:
        y_D = np.sum(bin_centers**2 * f_y * bin_widths) / y_F
    else:
        y_D = 0.0
    return {"y_F": y_F, "y_D": y_D}


# ---------------------------------------------------------------------------
# Stopping power loading and kappa computation
# ---------------------------------------------------------------------------


def _load_stopping_powers(filepath):
    """Load stopping power CSV file.

    Parameters
    ----------
    filepath : str or Path
        Path to CSV with columns: energy_MeV, stopping_power_MeV_cm2_per_g.

    Returns
    -------
    tuple
        (energy_MeV, stopping_power) as ndarrays.
    """
    df = pd.read_csv(filepath)
    return df["energy_MeV"].values, df["stopping_power_MeV_cm2_per_g"].values


def compute_kappa_table(water_csv_path=None, sic_csv_path=None):
    """Compute tissue-equivalence correction factor kappa(E) = S_water / S_SiC.

    Parameters
    ----------
    water_csv_path : str or Path, optional
        Path to water stopping power CSV. Defaults to data/stopping_power_water.csv.
    sic_csv_path : str or Path, optional
        Path to SiC stopping power CSV. Defaults to data/stopping_power_sic.csv.

    Returns
    -------
    dict
        energy_MeV : ndarray -- common energy grid
        kappa : ndarray -- tissue-equivalence correction factor
    """
    data_dir = Path(__file__).parent.parent / "data"
    if water_csv_path is None:
        water_csv_path = data_dir / "stopping_power_water.csv"
    if sic_csv_path is None:
        sic_csv_path = data_dir / "stopping_power_sic.csv"

    e_water, s_water = _load_stopping_powers(water_csv_path)
    e_sic, s_sic = _load_stopping_powers(sic_csv_path)

    # Interpolate to common energy grid (use water grid as reference)
    s_sic_interp = np.interp(e_water, e_sic, s_sic)

    kappa = s_water / s_sic_interp

    logger.info(
        "Kappa table computed: %d points, range [%.3f, %.3f]",
        len(kappa),
        np.min(kappa),
        np.max(kappa),
    )

    return {"energy_MeV": e_water, "kappa": kappa}


# ---------------------------------------------------------------------------
# Tissue-equivalence correction
# ---------------------------------------------------------------------------


def tissue_equivalence_correction(
    y_values, event_energies_keV, kappa_table=None, kappa_constant=0.58
):
    """Apply tissue-equivalence correction to lineal energy values.

    Converts SiC detector lineal energy to tissue-equivalent values using
    kappa = S_tissue / S_SiC.

    Parameters
    ----------
    y_values : array_like
        Lineal energy values in keV/um (SiC response).
    event_energies_keV : array_like
        Energy per event in keV (used for energy-dependent kappa lookup).
    kappa_table : dict, optional
        Dictionary with energy_MeV and kappa arrays from compute_kappa_table.
        If None, uses constant kappa approximation.
    kappa_constant : float
        Constant kappa value used when kappa_table is None. Default 0.58.

    Returns
    -------
    ndarray
        Tissue-equivalent lineal energy values in keV/um.
    """
    y_values = np.asarray(y_values, dtype=np.float64)
    event_energies_keV = np.asarray(event_energies_keV, dtype=np.float64)

    if kappa_table is not None:
        # Energy-dependent kappa: interpolate per event
        event_energies_MeV = event_energies_keV / 1e3
        kappa_per_event = np.interp(
            event_energies_MeV,
            kappa_table["energy_MeV"],
            kappa_table["kappa"],
        )
        y_tissue = kappa_per_event * y_values
        logger.info(
            "Energy-dependent tissue correction: kappa range [%.3f, %.3f]",
            np.min(kappa_per_event),
            np.max(kappa_per_event),
        )
    else:
        logger.warning(
            "Using constant kappa = %.3f for tissue-equivalence correction. "
            "This is an approximation; kappa varies ~20-30%% across the "
            "energy range. Use compute_kappa_table() for energy-dependent "
            "correction.",
            kappa_constant,
        )
        y_tissue = kappa_constant * y_values

    return y_tissue


# ---------------------------------------------------------------------------
# Publication-quality plotting
# ---------------------------------------------------------------------------


def plot_yd_spectrum(
    bin_centers, d_y, y_F=None, y_D=None, ax=None, label=None, color=None
):
    """Plot y*d(y) vs y (semilog-x) following microdosimetric conventions.

    The standard microdosimetric dose spectrum: area under curve is
    proportional to dose in each lineal energy interval.

    Parameters
    ----------
    bin_centers : ndarray
        Lineal energy bin centers in keV/um.
    d_y : ndarray
        Dose distribution values.
    y_F : float, optional
        Frequency-mean lineal energy for vertical line annotation.
    y_D : float, optional
        Dose-mean lineal energy for vertical line annotation.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    label : str, optional
        Legend label for the spectrum curve.
    color : str, optional
        Color for the spectrum curve.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    yd_y = bin_centers * d_y
    plot_kwargs = {}
    if label is not None:
        plot_kwargs["label"] = label
    if color is not None:
        plot_kwargs["color"] = color

    ax.semilogx(bin_centers, yd_y, **plot_kwargs)

    if y_F is not None:
        ax.axvline(
            y_F,
            color="gray",
            linestyle="--",
            linewidth=0.8,
            label=f"$y_F$ = {y_F:.2f} keV/$\\mu$m",
        )
    if y_D is not None:
        ax.axvline(
            y_D,
            color="gray",
            linestyle="-",
            linewidth=0.8,
            label=f"$y_D$ = {y_D:.2f} keV/$\\mu$m",
        )

    ax.set_xlabel("y (keV/$\\mu$m)")
    ax.set_ylabel("y $\\cdot$ d(y)")
    ax.set_title("Dose-weighted lineal energy spectrum")
    ax.legend(fontsize=8)

    return ax


def plot_yf_spectrum(bin_centers, f_y, y_F=None, ax=None, label=None, color=None):
    """Plot y*f(y) vs y (semilog-x) following microdosimetric conventions.

    Parameters
    ----------
    bin_centers : ndarray
        Lineal energy bin centers in keV/um.
    f_y : ndarray
        Frequency distribution values.
    y_F : float, optional
        Frequency-mean lineal energy for vertical line annotation.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    label : str, optional
        Legend label for the spectrum curve.
    color : str, optional
        Color for the spectrum curve.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    yf_y = bin_centers * f_y
    plot_kwargs = {}
    if label is not None:
        plot_kwargs["label"] = label
    if color is not None:
        plot_kwargs["color"] = color

    ax.semilogx(bin_centers, yf_y, **plot_kwargs)

    if y_F is not None:
        ax.axvline(
            y_F,
            color="gray",
            linestyle="--",
            linewidth=0.8,
            label=f"$y_F$ = {y_F:.2f} keV/$\\mu$m",
        )

    ax.set_xlabel("y (keV/$\\mu$m)")
    ax.set_ylabel("y $\\cdot$ f(y)")
    ax.set_title("Frequency-weighted lineal energy spectrum")
    ax.legend(fontsize=8)

    return ax
