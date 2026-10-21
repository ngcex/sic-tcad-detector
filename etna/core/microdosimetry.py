"""Microdosimetric spectra computation for 4H-SiC microdosimeter.

Transforms pulse height distributions (energy deposited per event) into
standard microdosimetric observables following ICRU Report 36 conventions:

- Lineal energy: y = epsilon / l_bar (keV/um)
- Frequency distribution f(y): normalized so integral f(y)*dy = 1
- Dose distribution d(y): d(y) = y * f(y) / y_F
- Frequency-mean lineal energy: y_F = integral(y * f(y) * dy)
- Dose-mean lineal energy: y_D = integral(y * d(y) * dy)
- Tissue-equivalence correction via kappa(E) = L_tissue / L_SiC (LET ratio at the
  particle KINETIC energy, real NIST PSTAR/ASTAR data, see compute_kappa_table)

Also provides publication-quality plotting functions for y*d(y) vs log(y)
and y*f(y) vs log(y) spectra.

References:
    - ICRU Report 36 (1983): Microdosimetry definitions and conventions
    - Cauchy theorem: l_bar = 4V/S for convex bodies under isotropic irradiation
    - Phase 22: mc_coupling.py (process_mc_ensemble output)
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Mean chord length
# ---------------------------------------------------------------------------


def mean_chord_length(sv_thickness_um, sv_width_um=None, sv_depth_um=None):
    """Compute mean chord length for the sensitive volume geometry.

    For a rectangular parallelepiped: l_bar = 4V/S (Cauchy's theorem,
    valid for convex bodies under isotropic irradiation).
    If width or depth is missing: slab approximation l_bar = 2*thickness
    (valid when lateral dimensions >> thickness). run_microdosimetry passes both
    (depth defaults to the width there, audit-v7 D-15 follow-up); a width alone
    still gives the slab value here (estimate_noise_floor relies on that).

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
        n_events : int -- number of input events
        n_in_range, n_zero, n_out_of_range : int -- events in f(y) (y in range),
            with epsilon <= 0, with y > 0 outside the bins (the last two excluded)
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
    # f(y) is the distribution of events that deposit energy (ICRU 36): epsilon = 0
    # events (and positive y outside [y_min, y_max], which cannot be binned) are not
    # part of the normalisation. Dividing by all events scaled y_F by the non-zero
    # fraction while leaving y_D unchanged (audit v8 finding 5).
    n_in = int(counts.sum())
    n_zero = int(np.sum(y_values <= 0))
    n_out = n_total - n_in - n_zero
    if n_out > 0:
        logger.warning(
            "%d of %d events with y > 0 fall outside [%.2e, %.2e] keV/um and are "
            "excluded from f(y)",
            n_out,
            n_total,
            y_min,
            y_max,
        )

    # f(y): normalized frequency distribution, integral f(y)*dy = 1
    f_y = np.zeros_like(bin_centers)
    nonzero = bin_widths > 0
    if n_in > 0:
        f_y[nonzero] = counts[nonzero] / (n_in * bin_widths[nonzero])

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
        "n_in_range": n_in,
        "n_zero": n_zero,
        "n_out_of_range": n_out,
        "y_values": y_values,
    }


# ---------------------------------------------------------------------------
# Microdosimetric means (standalone convenience)
# ---------------------------------------------------------------------------


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
    import pandas as pd  # lazy: keeps `import etna` light

    df = pd.read_csv(filepath)
    return df["energy_MeV"].values, df["stopping_power_MeV_cm2_per_g"].values


# 4H-SiC mass fractions (Si 28.0855, C 12.011 g/mol; 1:1 stoichiometry). Single
# source: etna.core.stopping_power (re-exported here for backward compatibility).
from etna.core.stopping_power import (  # noqa: E402
    SIC_MASS_FRACTION_C,
    SIC_MASS_FRACTION_SI,
)


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

    with w_Si = 0.7004, w_C = 0.2996 for 4H-SiC. Generic array helper; the bundled
    NIST data path is ``etna.core.stopping_power.sic_mass_stopping`` (used by
    ``compute_kappa_table(source="bragg")``).

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


def compute_kappa_table(
    water_csv_path=None,
    sic_csv_path=None,
    source=None,
    particle="proton",
    tissue="water",
    ratio="linear",
    carbon="c_amorphous",
):
    """SiC -> tissue lineal-energy conversion factor kappa(E) vs KINETIC energy.

    ``source="bragg"`` (the physics path) uses REAL NIST data (PSTAR for protons,
    ASTAR for alphas; ``data/nist_star/``, downloaded 2026-10-20 with provenance).
    NIST has no SiC table: SiC is composed by Bragg additivity of NIST Si + C
    (``etna.core.stopping_power``; caveat: binding/phase effects ignored, roughly a few
    % for protons > 1 MeV and up to ~10% near the stopping maximum -- an unverified
    literature estimate).

    Convention (``ratio``):

    - ``"linear"`` (default): kappa = L_tissue / L_SiC, the ratio of LINEAR
      electronic stopping powers (LET ratio). For a track crossing the SV, lineal
      energy y ~ LET, so y_tissue = kappa * y_SiC. This is the convention of the Si
      microdosimetry literature (Si -> muscle ~0.58 for protons; from these same
      NIST files: 0.584 at 10 MeV). For SiC -> water: 0.426 @1 MeV, 0.391 @10 MeV,
      0.378 @62 MeV, 0.375 @150 MeV (protons).
    - ``"mass"``: (S/rho)_tissue / (S/rho)_SiC (1.368 @1 MeV -> 1.204 @150 MeV for
      protons): the Bragg-Gray DOSE conversion factor, NOT a lineal-energy factor.

    History (C-1, v5 audit): the old ``data/stopping_power_{water,sic}.csv`` tables
    are fabricated (water 274 MeV cm^2/g at 0.1 MeV vs PSTAR 817) and give a flat
    ~0.58. The v5/v6 text called that value "sign-inverted" against an expected mass
    ratio ~1.13-1.24; with the linear (lineal-energy) convention the real SiC value
    is ~0.38-0.43, so 0.58 is wrong in magnitude and flatness, and the quoted
    1.13-1.24 range matched neither convention. The legacy path is kept only as the
    explicit ``source="legacy_unsafe"`` opt-in.

    .. important::
       There is deliberately NO silent default for ``source``: omitting it raises
       ``ValueError`` (v6 physics-honesty fix).

    Parameters
    ----------
    water_csv_path, sic_csv_path : str or Path, optional
        Override input paths (``source="legacy_unsafe"`` only).
    source : {"bragg", "legacy_unsafe"}
        Required.
    particle : {"proton", "alpha"}
        kappa depends on the particle: never apply a proton table to alpha events.
    tissue : {"water", "muscle_icru"}
        Tissue reference (NIST liquid water or ICRU striated muscle).
    ratio : {"linear", "mass"}
        See above. Default "linear" (lineal-energy / LET conversion).
    carbon : {"c_amorphous", "c_graphite"}
        NIST carbon used in the Bragg composition (<= 1% effect above 1 MeV).

    Returns
    -------
    dict
        energy_MeV, kappa (ndarrays on the NIST kinetic-energy grid), source,
        particle, tissue, ratio.

    Raises
    ------
    ValueError
        If ``source`` is omitted or not one of ``"bragg"``/``"legacy_unsafe"``.
    """
    from etna._paths import data_root

    data_dir = data_root()

    if source not in ("bragg", "legacy_unsafe"):
        raise ValueError(
            "compute_kappa_table() requires an explicit source= kwarg -- there "
            "is no silent default. Pass source='bragg' for the real NIST "
            "PSTAR/ASTAR physics path, or source='legacy_unsafe' to "
            "deliberately opt into the FABRICATED flat ~0.58 placeholder "
            "(real SiC->water LET ratio is ~0.38-0.43 for protons) for pipeline "
            f"illustration only. Got source={source!r}."
        )

    if source == "bragg":
        from etna.core.stopping_power import kappa_table

        res = kappa_table(particle=particle, tissue=tissue, ratio=ratio, carbon=carbon)
        logger.info(
            "Kappa(bragg, %s, %s, %s): %d points, range [%.3f, %.3f]",
            particle,
            tissue,
            ratio,
            len(res["kappa"]),
            float(np.min(res["kappa"])),
            float(np.max(res["kappa"])),
        )
        return res

    # legacy_unsafe (v3.0 flat-kappa artefact) -- explicit opt-in only, see above.
    if water_csv_path is None:
        water_csv_path = data_dir / "stopping_power_water.csv"
    if sic_csv_path is None:
        sic_csv_path = data_dir / "stopping_power_sic.csv"
    e_water, s_water = _load_stopping_powers(water_csv_path)
    e_sic, s_sic = _load_stopping_powers(sic_csv_path)
    s_sic_interp = np.interp(e_water, e_sic, s_sic)
    kappa = s_water / s_sic_interp
    logger.warning(
        "Kappa(legacy_unsafe): %d points, range [%.3f, %.3f] -- FABRICATED "
        "PLACEHOLDER, NOT REAL PHYSICS (real SiC->water LET ratio ~0.38-0.43 "
        "for protons, see C-1). Caller explicitly opted in via "
        "source='legacy_unsafe'.",
        len(kappa),
        float(np.min(kappa)),
        float(np.max(kappa)),
    )
    return {"energy_MeV": e_water, "kappa": kappa, "source": "legacy_unsafe"}


# ---------------------------------------------------------------------------
# Tissue-equivalence correction
# ---------------------------------------------------------------------------


_UNSET = object()


def tissue_equivalence_correction(
    y_values,
    event_energies_keV=None,
    kappa_table=None,
    kappa_constant=_UNSET,
    particle_energy_MeV=None,
    particle=None,
):
    """Apply tissue-equivalence correction to lineal energy values.

    Converts SiC detector lineal energy to tissue-equivalent values,
    y_tissue = kappa(E) * y_SiC, with kappa the SiC -> tissue LET ratio from
    ``compute_kappa_table(source="bragg")`` (real NIST PSTAR/ASTAR data).

    AUDIT C-2 (v5, fixed): kappa(E) is a function of the particle KINETIC energy,
    NOT the energy deposited in the thin SV (for a penetrating particle the
    deposited energy is ~1e4x below the kinetic energy). Pass the kinetic energy via
    ``particle_energy_MeV``:

    - scalar -> mono-energetic beam (e.g. 62.0 for a 62 MeV proton beam);
    - array (len == y_values) -> per-event primary kinetic energy from MC truth, e.g.
      ``process_mc_ensemble(...)["event_kinetic_MeV"]`` (needs an ``ekin_MeV``
      column in the MC events, see ``etna.core.mc_coupling``).

    The lookup is log-log interpolation on the NIST grid and raises outside it. If
    ``particle_energy_MeV`` is omitted and a ``kappa_table`` is given, the function
    falls back to the mean of the tabulated kappa (crude documented approximation,
    logged as a warning) rather than using the deposited energy.

    .. important::
       There is deliberately NO silent default for ``kappa_constant`` (v6
       physics-honesty fix). If ``kappa_table`` is None, pass ``kappa_constant``
       explicitly: a real value (SiC -> water LET ratio ~0.38-0.43 for protons, see
       ``compute_kappa_table``) or the literal string ``"legacy_unsafe"`` to opt
       into the fabricated flat 0.58 placeholder. Omitting it raises ``ValueError``.

    Parameters
    ----------
    y_values : array_like
        Lineal energy values in keV/um (SiC response).
    event_energies_keV : array_like, optional
        DEPRECATED for kappa lookup (this is deposited energy, the C-2 bug). Kept
        only for backward-compatible call signatures; ignored for the lookup.
    kappa_table : dict, optional
        From ``compute_kappa_table``. If None, uses ``kappa_constant``.
    kappa_constant : float or "legacy_unsafe"
        Constant kappa when no table is given (required in that case).
    particle_energy_MeV : float or array_like, optional
        Particle KINETIC energy (MeV) for the kappa(E) lookup.
    particle : {"proton", "alpha"}, optional
        If given, must match ``kappa_table["particle"]`` (kappa is particle-specific).

    Returns
    -------
    ndarray
        Tissue-equivalent lineal energy values in keV/um.

    Raises
    ------
    ValueError
        If ``kappa_table`` is None and ``kappa_constant`` is omitted, if the particle
        does not match the table, or if a kinetic energy is outside the table.
    """
    y_values = np.asarray(y_values, dtype=np.float64)

    if kappa_table is None:
        if kappa_constant is _UNSET:
            raise ValueError(
                "tissue_equivalence_correction() requires an explicit "
                "kappa_constant= when kappa_table is None -- there is no "
                "silent default. Pass a real value (SiC->water LET ratio "
                "~0.38-0.43 for protons, see compute_kappa_table) or "
                "kappa_constant='legacy_unsafe' to deliberately opt into the "
                "FABRICATED flat 0.58 placeholder for pipeline illustration only."
            )
        if kappa_constant == "legacy_unsafe":
            kappa_constant = 0.58
            logger.warning(
                "Constant kappa = %.3f -- FABRICATED PLACEHOLDER, NOT REAL PHYSICS "
                "(real SiC->water LET ratio ~0.38-0.43 for protons, see C-1). "
                "Caller explicitly opted in via kappa_constant='legacy_unsafe'.",
                kappa_constant,
            )
        else:
            logger.info(
                "Constant kappa = %.3f (caller-supplied; the real SiC->water LET "
                "ratio varies with kinetic energy -- use "
                "compute_kappa_table(source='bragg') for kappa(E)).",
                kappa_constant,
            )
        return kappa_constant * y_values

    table_particle = kappa_table.get("particle")
    if (
        particle is not None
        and table_particle is not None
        and particle != table_particle
    ):
        raise ValueError(
            f"kappa_table is for {table_particle!r} but particle={particle!r}: kappa "
            "is particle-specific, build the table with compute_kappa_table("
            f"source='bragg', particle={particle!r})."
        )

    E_grid = np.asarray(kappa_table["energy_MeV"], float)
    K_grid = np.asarray(kappa_table["kappa"], float)

    if particle_energy_MeV is not None:
        # C-2 correct path: lookup by KINETIC energy (scalar beam or per-event).
        Ek = np.asarray(particle_energy_MeV, dtype=np.float64)
        if (
            np.any(~np.isfinite(Ek))
            or np.any(Ek < E_grid[0])
            or np.any(Ek > E_grid[-1])
        ):
            raise ValueError(
                "particle_energy_MeV must be finite and inside the kappa table "
                f"[{E_grid[0]:g}, {E_grid[-1]:g}] MeV"
            )
        kappa = np.exp(np.interp(np.log(Ek), np.log(E_grid), np.log(K_grid)))
        y_tissue = kappa * y_values
        logger.info(
            "Energy-dependent tissue correction (kinetic E): kappa [%.3f, %.3f]",
            float(np.min(kappa)),
            float(np.max(kappa)),
        )
        return y_tissue

    # No kinetic energy supplied: fall back to the table-MEAN kappa, NOT the
    # deposited-energy lookup (which was the C-2 defect). Crude approximation.
    kappa_avg = float(np.mean(K_grid))
    logger.warning(
        "tissue_equivalence_correction: no particle_energy_MeV given; using the "
        "mean of the tabulated kappa = %.3f as a crude approximation. Pass the "
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
    import matplotlib.pyplot as plt  # lazy: keeps `import etna` free of pyplot

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
    if ax.get_legend_handles_labels()[0]:  # nothing labelled -> no empty legend
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
    import matplotlib.pyplot as plt  # lazy: keeps `import etna` free of pyplot

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
    if ax.get_legend_handles_labels()[0]:  # nothing labelled -> no empty legend
        ax.legend(fontsize=8)

    return ax
