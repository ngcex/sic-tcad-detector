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


# 4H-SiC mass fractions (Si 28.0855, C 12.011 g/mol; 1:1 stoichiometry).
SIC_MASS_FRACTION_SI = 0.7004
SIC_MASS_FRACTION_C = 0.2996


def sic_stopping_power_bragg(
    e_si,
    s_si,
    e_c,
    s_c,
    energy_grid=None,
    w_si=SIC_MASS_FRACTION_SI,
    w_c=SIC_MASS_FRACTION_C,
):
    """Compose the SiC mass stopping power via Bragg additivity.

    Bragg's rule: the mass stopping power of a compound is the mass-fraction-
    weighted sum of its elemental mass stopping powers,

        (S/rho)_SiC(E) = w_Si * (S/rho)_Si(E) + w_C * (S/rho)_C(E)

    with w_Si = 0.7004, w_C = 0.2996 for 4H-SiC. This lets SiC stopping power be
    built from NIST PSTAR/ASTAR elemental proton data (Si and C) rather than from
    a single hand-entered SiC table -- the path the Phase 27 spec calls for.

    Both elemental inputs are interpolated (log-log) onto a common grid.

    Parameters
    ----------
    e_si, s_si : array_like
        Silicon proton energy (MeV) and mass stopping power (MeV cm^2/g).
    e_c, s_c : array_like
        Carbon proton energy (MeV) and mass stopping power (MeV cm^2/g).
    energy_grid : array_like, optional
        Common output energy grid (MeV). Defaults to the silicon energy grid.
    w_si, w_c : float
        SiC mass fractions (default 4H-SiC values).

    Returns
    -------
    tuple
        (energy_MeV, sic_stopping_power MeV cm^2/g) ndarrays.
    """
    e_si = np.asarray(e_si, float)
    e_c = np.asarray(e_c, float)
    grid = np.asarray(energy_grid, float) if energy_grid is not None else e_si

    def _loglog_interp(x, xp, fp):
        return np.exp(np.interp(np.log(x), np.log(xp), np.log(fp)))

    s_si_g = _loglog_interp(grid, e_si, np.asarray(s_si, float))
    s_c_g = _loglog_interp(grid, e_c, np.asarray(s_c, float))
    return grid, w_si * s_si_g + w_c * s_c_g


def compute_kappa_table(water_csv_path=None, sic_csv_path=None, source="legacy"):
    """Compute tissue-equivalence correction factor kappa(E) = S_water / S_SiC.

    .. danger::
       AUDIT C-1 (v5) -- the ``source="legacy"`` path reads
       ``data/stopping_power_{water,sic}.csv`` which are FABRICATED placeholders
       (kappa ~0.58, flat, sign-inverted). The physics requires kappa > 1 (water
       Z/A 0.555 vs SiC 0.499; I_water ~78 eV vs I_SiC ~136 eV); real NIST PSTAR
       gives ~1.24 at 1 MeV decreasing to ~1.13 at 100 MeV.

       Phase 27 machinery (``source="bragg"``) is now wired: it builds SiC by
       Bragg additivity from elemental Si + C proton stopping powers and divides
       real water PSTAR by it. **It is still DATA-BLOCKED**: the input files
       ``data/srim/{water,si,c}_proton.csv`` are placeholders that must be filled
       with real NIST PSTAR (water, Si, C) tabulated proton mass stopping powers
       (see ``data/srim/README.md``). Until then ``source="bragg"`` raises a clear
       error rather than returning a fabricated number. Network fetch is not used
       (offline/CI environment). See ``.planning/PHYSICS_AUDIT_v5.md`` (C-1/C-2).

    Parameters
    ----------
    water_csv_path, sic_csv_path : str or Path, optional
        Override input paths (legacy source only).
    source : {"legacy", "bragg"}
        "legacy" -- divide the (placeholder) water CSV by the SiC CSV (the v3.0
        flat-kappa artefact; kept for regression continuity).
        "bragg" -- compute SiC via Bragg additivity from elemental Si + C PSTAR
        data in ``data/srim/`` and divide real water PSTAR by it (Phase 27).

    Returns
    -------
    dict
        energy_MeV : ndarray, kappa : ndarray, source : str
    """
    data_dir = Path(__file__).parent.parent / "data"

    if source == "bragg":
        srim = data_dir / "srim"
        paths = {
            "water": srim / "water_proton.csv",
            "si": srim / "si_proton.csv",
            "c": srim / "c_proton.csv",
        }
        for name, p in paths.items():
            if not p.exists() or _is_placeholder_stopping_csv(p):
                raise FileNotFoundError(
                    f"compute_kappa_table(source='bragg') needs real PSTAR data at "
                    f"{p} ({name}). It is a placeholder/missing -- DO NOT fabricate. "
                    f"See data/srim/README.md for the exact NIST PSTAR files to drop in."
                )
        e_w, s_w = _load_stopping_powers(paths["water"])
        e_si, s_si = _load_stopping_powers(paths["si"])
        e_c, s_c = _load_stopping_powers(paths["c"])
        grid, s_sic = sic_stopping_power_bragg(e_si, s_si, e_c, s_c, energy_grid=e_w)
        kappa = s_w / s_sic
        logger.info(
            "Kappa(bragg): %d points, range [%.3f, %.3f]",
            len(kappa),
            float(np.min(kappa)),
            float(np.max(kappa)),
        )
        return {"energy_MeV": e_w, "kappa": kappa, "source": "bragg"}

    # legacy (v3.0 flat-kappa artefact)
    if water_csv_path is None:
        water_csv_path = data_dir / "stopping_power_water.csv"
    if sic_csv_path is None:
        sic_csv_path = data_dir / "stopping_power_sic.csv"
    e_water, s_water = _load_stopping_powers(water_csv_path)
    e_sic, s_sic = _load_stopping_powers(sic_csv_path)
    s_sic_interp = np.interp(e_water, e_sic, s_sic)
    kappa = s_water / s_sic_interp
    logger.info(
        "Kappa(legacy): %d points, range [%.3f, %.3f] (PLACEHOLDER -- see C-1)",
        len(kappa),
        float(np.min(kappa)),
        float(np.max(kappa)),
    )
    return {"energy_MeV": e_water, "kappa": kappa, "source": "legacy"}


def _is_placeholder_stopping_csv(path):
    """True if a stopping-power CSV is the documented placeholder (header marker)."""
    try:
        first = Path(path).read_text().splitlines()[0]
    except (OSError, IndexError):
        return True
    return "PLACEHOLDER" in first.upper()


# ---------------------------------------------------------------------------
# Tissue-equivalence correction
# ---------------------------------------------------------------------------


def tissue_equivalence_correction(
    y_values,
    event_energies_keV=None,
    kappa_table=None,
    kappa_constant=0.58,
    particle_energy_MeV=None,
):
    """Apply tissue-equivalence correction to lineal energy values.

    Converts SiC detector lineal energy to tissue-equivalent values using
    kappa = S_tissue / S_SiC.

    AUDIT C-2 (v5): kappa(E) is physically a function of the particle KINETIC
    energy, NOT the energy deposited in the thin SV (for a penetrating particle
    the deposited energy is ~1e4x below the kinetic energy). Pass the kinetic
    energy via ``particle_energy_MeV``:

    - scalar -> mono-energetic beam (e.g. 62.0 for a 62 MeV proton beam): exact,
      and the recommended path for clinical beams;
    - array (len == y_values) -> per-event primary kinetic energy from MC truth.

    If ``particle_energy_MeV`` is omitted and a ``kappa_table`` is given, the
    function falls back to the **energy-averaged** kappa (documented approximation)
    rather than silently using the wrong (deposited) energy variable.

    Parameters
    ----------
    y_values : array_like
        Lineal energy values in keV/um (SiC response).
    event_energies_keV : array_like, optional
        DEPRECATED for kappa lookup (this is deposited energy, the C-2 bug). Kept
        only for backward-compatible call signatures; ignored for the lookup.
    kappa_table : dict, optional
        From ``compute_kappa_table``. If None, uses ``kappa_constant``.
    kappa_constant : float
        Constant kappa when no table is given. Default 0.58 (legacy placeholder;
        set to ~1.2 once real PSTAR data is in place, see C-1).
    particle_energy_MeV : float or array_like, optional
        Particle KINETIC energy (MeV) for the kappa(E) lookup -- scalar (beam) or
        per-event array. Correct variable per C-2.

    Returns
    -------
    ndarray
        Tissue-equivalent lineal energy values in keV/um.
    """
    y_values = np.asarray(y_values, dtype=np.float64)

    if kappa_table is None:
        logger.warning(
            "Constant kappa = %.3f (approximation; real kappa(E) varies ~10-30%%. "
            "Use compute_kappa_table(source='bragg') once PSTAR data is in place).",
            kappa_constant,
        )
        return kappa_constant * y_values

    E_grid = np.asarray(kappa_table["energy_MeV"], float)
    K_grid = np.asarray(kappa_table["kappa"], float)

    if particle_energy_MeV is not None:
        # C-2 correct path: lookup by KINETIC energy (scalar beam or per-event).
        Ek = np.asarray(particle_energy_MeV, dtype=np.float64)
        kappa = np.interp(Ek, E_grid, K_grid)
        y_tissue = kappa * y_values
        kmin, kmax = float(np.min(kappa)), float(np.max(kappa))
        logger.info(
            "Energy-dependent tissue correction (kinetic E): kappa [%.3f, %.3f]",
            kmin,
            kmax,
        )
        return y_tissue

    # No kinetic energy supplied: fall back to ENERGY-AVERAGED kappa, NOT the
    # deposited-energy lookup (which was the C-2 defect). Documented approximation.
    kappa_avg = float(np.mean(K_grid))
    logger.warning(
        "tissue_equivalence_correction: no particle_energy_MeV given; using "
        "energy-averaged kappa = %.3f as a documented approximation. Pass the "
        "particle KINETIC energy (scalar beam or per-event MC truth) for the "
        "correct energy-dependent correction (AUDIT C-2).",
        kappa_avg,
    )
    return kappa_avg * y_values


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
