"""Radiation generation rate profiles for 4H-SiC detector simulation.

Provides analytical models for carrier generation from:
- Alpha particles (Am-241, 5.486 MeV): by default the real Bragg curve in SiC from
  NIST ASTAR (Si + C, Bragg additivity), legacy flat box as opt-out
- Proton beams: dose-rate in SiC (default; energy-independent, accepted for 30-150 MeV
  therapeutic energies only) or, opt-in, water-dose or fluence-rate input with the NIST
  PSTAR energy dependence (any energy of the PSTAR grid)

Energy ranges (2026-10-21): the default mode's 30-150 MeV check is a scope choice of the
flat entrance-region profile, not a data limit. PSTAR covers 1 keV - 10 GeV and the NIEL
hardness table (etna.core.niel_tables / radiation_damage) 1-1000 MeV, but below ~10 MeV
a flat profile is wrong for a 10 um epi: L_e,SiC(10 um) / L_e,SiC(0) = 1.009 at 10 MeV,
1.031 at 5 MeV, 1.17 at 2 MeV, and a 1 MeV proton stops at 10.8 um (NIST PSTAR Si + C,
csda_range_cm / residual_energy_MeV). The range was therefore not widened to the NIEL
table's; use ``fluence_rate_cm2_s`` (depth-resolved) for other energies.
- Dose rate to generation rate conversion

All units in CGS (cm, cm^-3, eV, s) per devsim convention.

References:
    - NIST PSTAR / ASTAR (data/nist_star/, downloaded 2026-10-20)
    - Ioffe NSM Archive for 4H-SiC material properties
"""

import numpy as np

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

from etna.core.sic_material import SiC4H_Parameters
from etna.core.stopping_power import (
    AM241_ALPHA_E0_MEV,
    alpha_bragg_curve_sic,
    csda_range_cm,
    linear_stopping,
    mass_stopping,
    residual_energy_MeV,
)

_params = SiC4H_Parameters()
RHO_SIC = _params.rho  # g/cm^3, 4H-SiC density
E_PAIR_SIC_EV = _params.E_pair_eV  # eV, electron-hole pair creation energy in 4H-SiC

# NIST PSTAR CSDA ranges in water (mm), read from data/nist_star/pstar_water.csv
# (downloaded 2026-10-20). The 62 MeV entry was a hand-typed 31.0 before (PSTAR, log-log interpolated: 32.81).
PROTON_RANGE_WATER_MM = {
    E: round(float(csda_range_cm("proton", "water", float(E))) * 10.0, 3)
    for E in (30, 62, 70, 150)
}

# CSDA range of a 5.486 MeV (Am-241) alpha in 4H-SiC, NIST ASTAR Si + C by Bragg
# additivity (18.43 um). Was a hand-set 15 um before 2026-10-20.
ALPHA_RANGE_SIC_CM = float(csda_range_cm("alpha", "sic", AM241_ALPHA_E0_MEV))


# ---------------------------------------------------------------------------
# Dose rate conversion
# ---------------------------------------------------------------------------


def dose_rate_to_generation(dose_rate_Gy_s, rho_g_cm3=RHO_SIC, E_pair_eV=E_PAIR_SIC_EV):
    """Convert dose rate (Gy/s) to carrier generation rate (cm^-3 s^-1).

    Uses CGS conversion: 1 Gy = 1e4 erg/g.

    Parameters
    ----------
    dose_rate_Gy_s : float or array
        Dose rate in Gray per second.
    rho_g_cm3 : float
        Material density (g/cm^3). Default: 3.21 (4H-SiC).
    E_pair_eV : float
        Electron-hole pair creation energy (eV). Default: 7.8 (4H-SiC).

    Returns
    -------
    G : float or array
        Generation rate (cm^-3 s^-1).

    Notes
    -----
    G = dose_rate * rho * 1e4 / (E_pair_eV * 1.602e-12)

    where 1e4 converts Gy -> erg/g and 1.602e-12 converts eV -> erg.
    """
    return dose_rate_Gy_s * rho_g_cm3 * 1e4 / (E_pair_eV * 1.602e-12)


# ---------------------------------------------------------------------------
# Alpha particle generation profile
# ---------------------------------------------------------------------------


def _nist_alpha_curve_unit(E0_MeV=AM241_ALPHA_E0_MEV):
    """NIST (ASTAR, SiC by Bragg additivity) alpha ionization curve on a unit depth axis.

    Returns ``(u, f)`` with ``u = depth / CSDA range`` (uniform, 0..1.3) and ``f`` the
    electronic dE/dx normalised to 1.0 at the entrance, Bohr range straggling included.
    """
    s, f, R, _ = alpha_bragg_curve_sic(E0_MeV)
    return s / R, f / f[0]


def alpha_depth_shape(x_cm, alpha_range_cm=ALPHA_RANGE_SIC_CM, bragg_peak_ratio="nist"):
    """Unnormalized spatial shape of the alpha e-h pair generation profile.

    Pure function of `x_cm` — no normalization, no assumption that `x_cm` is sorted or
    evenly spaced, safe to evaluate on unordered mesh node arrays (e.g. a 2D device's
    node coordinates).

    Parameters
    ----------
    x_cm : array_like
        Depth positions in SiC (cm), measured from detector entrance.
        Values are NOT required to be sorted.
    alpha_range_cm : float
        Alpha particle range in SiC (cm). Default: ``ALPHA_RANGE_SIC_CM``, the CSDA
        range of a 5.486 MeV (Am-241) alpha in 4H-SiC from NIST ASTAR Si + C by Bragg
        additivity (18.43 um; it was a hand-set 15 um before 2026-10-20).
    bragg_peak_ratio : "nist" (default), None or float
        - ``"nist"``: the REAL Bragg curve. Electronic stopping power S_e(E(x)) of a
          5.486 MeV alpha along its CSDA path in SiC (NIST ASTAR Si + C, Bragg
          additivity), with the computed Bohr range straggling (0.165 um), depth axis
          scaled so the CSDA range equals ``alpha_range_cm``. Entrance = 1.0, peak =
          2.49 at 0.92 R (data-derived, see ``stopping_power.alpha_bragg_summary``).
        - ``None``: opt-out, the pre-2026-10 flat box (1.0 to 0.8 R, erfc roll-off of
          width 0.1 R; no Bragg peak).
        - float: opt-out, the old qualitative linear ramp from 1.0 to this ratio over
          0..0.8 R, then the erfc roll-off.

    Returns
    -------
    shape : ndarray
        Unnormalized profile shape, same length as x_cm, 1.0 at the entrance; callers
        normalize as needed (to total deposited energy, to a peak generation rate,
        etc.).
    """
    x = np.asarray(x_cm, dtype=float)
    R = alpha_range_cm

    if isinstance(bragg_peak_ratio, str):
        if bragg_peak_ratio != "nist":
            raise ValueError(
                f"bragg_peak_ratio must be 'nist', None or a float (got {bragg_peak_ratio!r})"
            )
        u, f = _nist_alpha_curve_unit()
        return np.interp(x / R, u, f, left=1.0, right=0.0)

    # Transition region parameters (legacy box / ramp shapes)
    x_flat = 0.8 * R  # flat/rising region ends at 80% of range
    sigma = 0.1 * R  # roll-off width (10% of range)

    if bragg_peak_ratio is None:
        shape_factor = np.ones_like(x)
    else:
        # Linear ramp from 1.0 (entrance) to bragg_peak_ratio (end of the flat/rising
        # region, just before the erfc roll-off) -- qualitative legacy option.
        shape_factor = 1.0 + (bragg_peak_ratio - 1.0) * np.clip(x / x_flat, 0.0, 1.0)

    # Smooth roll-off to zero beyond the flat/rising region.
    # erfc((x - x_flat) / (sqrt(2) * sigma)) transitions from 2 to 0.
    # Divide by 2 to get transition from 1 to 0.
    from scipy.special import erfc  # lazy: keeps `import etna` light

    rolloff = 0.5 * erfc((x - x_flat) / (np.sqrt(2) * sigma))
    return shape_factor * rolloff


def alpha_generation_profile(
    x_cm,
    alpha_range_cm=ALPHA_RANGE_SIC_CM,
    total_energy_eV=5.486e6,
    E_pair_eV=E_PAIR_SIC_EV,
    bragg_peak_ratio="nist",
):
    """Spatial e-h pair generation profile for Am-241 alpha particles in SiC.

    Evaluates `alpha_depth_shape` and normalizes it so the total number
    of generated e-h pairs matches the alpha's deposited energy. Requires
    `x_cm` to be SORTED and evenly-representative of the integration
    domain (uses `np.diff`/trapezoidal integration) — for an unordered
    array (e.g. raw mesh node coordinates), use `alpha_depth_shape`
    directly and normalize separately (see `charge_collection_2d.py`'s
    `cce_lateral_scan` for an example: peak-normalized, not
    energy-normalized).

    Parameters
    ----------
    x_cm : array_like
        Depth positions in SiC (cm), measured from detector entrance.
        Must be sorted (ascending or descending) for the trapezoidal
        normalization integral to be meaningful.
    alpha_range_cm : float
        Alpha particle range in SiC (cm). Default: ``ALPHA_RANGE_SIC_CM`` (18.43 um,
        NIST-derived, see `alpha_depth_shape`).
    total_energy_eV : float
        Total kinetic energy of alpha particle (eV). Default: 5.486e6 (Am-241).
    E_pair_eV : float
        Electron-hole pair creation energy (eV). Default: 7.8 (4H-SiC).
    bragg_peak_ratio : "nist" (default), None or float
        See `alpha_depth_shape` (None = legacy flat box, opt-out).

    Returns
    -------
    G : ndarray
        Generation profile (cm^-1): pairs per cm depth per incident alpha.
        Multiply by fluence rate (alphas/cm^2/s) to get cm^-3 s^-1.

    Notes
    -----
    The profile is normalized so that integral(G(x) dx) = total_energy_eV / E_pair_eV
    over the given ``x_cm`` grid, REGARDLESS of the shape — the shape changes where the
    pairs are, not how many. If ``x_cm`` stops before the end of the track, the full
    pair count is redistributed over the grid (the historical convention; for the
    in-epi fraction see ``charge_collection_2d.in_mesh_generation_fraction``).
    """
    x = np.asarray(x_cm, dtype=float)
    profile = alpha_depth_shape(x, alpha_range_cm, bragg_peak_ratio)

    # Normalize so integral = total pairs per alpha (independent of shape)
    total_pairs = total_energy_eV / E_pair_eV
    dx = np.diff(x)
    # Trapezoidal integration for normalization
    integral = np.sum(0.5 * (profile[:-1] + profile[1:]) * dx)

    if integral > 0:
        G = profile * (total_pairs / integral)
    else:
        G = np.zeros_like(x)

    return G


# ---------------------------------------------------------------------------
# Tilted alpha track (microdosimeter alpha test, Petringa 2025 Phys. Med. 140)
# ---------------------------------------------------------------------------

# Parameters of the LEGACY qualitative two-power-law curve (shape="bragg_powerlaw").
# They were literature orders of magnitude, not data; the NIST ASTAR curve used by the
# default shape="bragg" gives S_e peak/entrance 2.50 (5.5 MeV, before straggling) and
# the stopping maximum at ~0.50 MeV (stopping_power.alpha_bragg_summary).
ALPHA_PEAK_TO_ENTRANCE = 2.3  # legacy, not data
ALPHA_E_PEAK_EV = 0.75e6  # eV, legacy, not data
ALPHA_STRAGGLE_FRAC = 0.02  # legacy range straggling sigma / range (NIST+Bohr: 0.009)


def alpha_path_stopping(
    s_cm,
    alpha_range_cm=18e-4,
    shape="bragg",
    E0_eV=5.5e6,
    peak_to_entrance=ALPHA_PEAK_TO_ENTRANCE,
    E_peak_eV=ALPHA_E_PEAK_EV,
    straggle_frac=ALPHA_STRAGGLE_FRAC,
):
    """Energy loss per unit PATH length dE/ds (eV/cm) along an alpha track.

    ``s_cm`` must be a uniform, ascending grid starting at 0 (the convolution below
    assumes it). The result integrates to ``E0_eV`` when the grid covers the track.

    shape="bragg" (default): the REAL curve, electronic S_e(E(s)) of an ``E0_eV``
    alpha along its CSDA path in 4H-SiC from NIST ASTAR Si + C (Bragg additivity),
    with the computed Bohr range straggling, path axis scaled so the CSDA range equals
    ``alpha_range_cm`` (NIST: 18.50 um at 5.5 MeV; the API default 18 um compresses it
    by 3%). Peak/entrance ~2.49. ``peak_to_entrance``/``E_peak_eV``/``straggle_frac``
    are ignored.
    shape="bragg_powerlaw": the legacy QUALITATIVE curve (before 2026-10-20): S(E) ~
    (E/E_pk)^-a above E_pk and ~ (E/E_pk)^0.5 below, a fixed by ``peak_to_entrance``,
    Gaussian straggling sigma = straggle_frac * range; peak/entrance ~2.1 after
    smoothing. Not data.
    shape="box": the legacy flat box (``alpha_depth_shape(..., bragg_peak_ratio=None)``,
    flat to 0.8 R, erfc roll-off) along the path.
    """
    s = np.asarray(s_cm, dtype=float)
    R = float(alpha_range_cm)
    if shape == "box":
        f = alpha_depth_shape(s, alpha_range_cm=R, bragg_peak_ratio=None)
    elif shape == "bragg":
        u, g = _nist_alpha_curve_unit(E0_eV / 1e6)
        f = np.interp(s / R, u, g, right=0.0)
    elif shape == "bragg_powerlaw":
        a = np.log(peak_to_entrance) / np.log(E0_eV / E_peak_eV)
        E = np.geomspace(1.0, E0_eV, 20001)
        S = np.where(E > E_peak_eV, (E / E_peak_eV) ** (-a), (E / E_peak_eV) ** 0.5)
        inv = 1.0 / S
        from_end = np.concatenate(
            [[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(E))]
        )
        scale = R / from_end[-1]
        s_of_E = (
            R - from_end * scale
        )  # path travelled when the energy is E (descending)
        f = np.interp(s, s_of_E[::-1], (S / scale)[::-1], right=0.0)
        f[s > R] = 0.0
        if len(s) > 1:
            sig = straggle_frac * R
            ds = s[1] - s[0]
            k = np.arange(-int(5 * sig / ds), int(5 * sig / ds) + 1) * ds
            ker = np.exp(-0.5 * (k / sig) ** 2)
            ker /= ker.sum()
            m = len(k) // 2
            # straggling smears the track END only: pad the entrance side with its value
            padded = np.concatenate([np.full(m, f[0]), f, np.zeros(m)])
            f = np.convolve(padded, ker, mode="valid")
    else:
        raise ValueError(
            f"shape must be 'bragg', 'bragg_powerlaw' or 'box' (got {shape!r})"
        )
    total = np.trapezoid(f, s)
    return f * (E0_eV / total) if total > 0 else f


def tilted_alpha_generation(
    x_cm,
    theta_deg=65.0,
    alpha_range_cm=18e-4,
    dead_layer_cm=0.0,
    shape="bragg",
    E0_eV=5.5e6,
    E_pair_eV=E_PAIR_SIC_EV,
    n_path=20001,
):
    """e-h pairs per cm of DEPTH per alpha, for a track tilted ``theta_deg`` from the normal.

    Depth ``x_cm`` is measured from the SiC surface (x = 0, the anode side of the p+);
    any order is accepted. A front dead layer (metal contact) of SiC-equivalent normal
    thickness ``dead_layer_cm`` is crossed first: depth x corresponds to path
    s = (x + dead_layer)/cos(theta), and G(x) = (dE/ds)(s) / cos(theta) / E_pair.
    The track ends at depth range*cos(theta) - dead_layer (18 um at 65 deg: 7.61 um).

    Returns ``(G, dead_layer_fraction)``: G (cm^-1, same shape as x; integrates to
    the pairs generated in SiC) and the fraction of E0 lost in the dead layer.
    """
    if not (0.0 <= theta_deg < 90.0):
        raise ValueError(f"theta_deg must be in [0, 90) (got {theta_deg!r})")
    if dead_layer_cm < 0:
        raise ValueError(f"dead_layer_cm must be >= 0 (got {dead_layer_cm!r})")
    c = np.cos(np.radians(theta_deg))
    s = np.linspace(0.0, 1.3 * alpha_range_cm, n_path)
    f = alpha_path_stopping(s, alpha_range_cm, shape, E0_eV)
    path = (np.asarray(x_cm, dtype=float) + dead_layer_cm) / c
    G = np.interp(path, s, f, right=0.0) / c / E_pair_eV
    s_dead = dead_layer_cm / c
    if s_dead > 0:
        sd = np.linspace(0.0, s_dead, 2001)
        lost = float(np.trapezoid(np.interp(sd, s, f), sd) / E0_eV)
    else:
        lost = 0.0
    return G, lost


# ---------------------------------------------------------------------------
# Proton Bragg peak generation profile
# ---------------------------------------------------------------------------


def _interpolate_proton_range_water(E_MeV):
    """Interpolate proton CSDA range in water for arbitrary energy.

    Uses log-log interpolation of tabulated NIST PSTAR data.

    Parameters
    ----------
    E_MeV : float
        Proton kinetic energy (MeV).

    Returns
    -------
    R_water_mm : float
        CSDA range in water (mm).

    Raises
    ------
    ValueError
        If energy is outside the tabulated range [30, 150] MeV.
    """
    energies = sorted(PROTON_RANGE_WATER_MM.keys())
    ranges = [PROTON_RANGE_WATER_MM[e] for e in energies]

    E_min, E_max = energies[0], energies[-1]
    if E_MeV < E_min or E_MeV > E_max:
        raise ValueError(
            f"Energy {E_MeV} MeV outside tabulated range " f"[{E_min}, {E_max}] MeV"
        )

    # Log-log interpolation (range ~ E^p power law)
    log_E = np.log(energies)
    log_R = np.log(ranges)
    R_water_mm = np.exp(np.interp(np.log(E_MeV), log_E, log_R))

    return float(R_water_mm)


def proton_generation_profile(
    x_cm, E_MeV, dose_rate_Gy_s=1.0, *, dose_medium="sic", fluence_rate_cm2_s=None
):
    """Carrier generation rate profile for proton beams in SiC detector.

    Three input conventions (audit-v7 B-10):

    - default, ``dose_medium="sic"``: ``dose_rate_Gy_s`` is the dose rate IN SiC. The
      generation G = D' rho / E_pair is then flat and ENERGY-INDEPENDENT by
      construction (at equal SiC dose rate, 30/70/150 MeV give the same G); ``E_MeV``
      is only validated against the 30-150 MeV range. Unchanged since v1.
    - ``dose_medium="water"`` (opt-in): ``dose_rate_Gy_s`` is the dose rate to WATER
      (the clinical quantity). It is converted to SiC dose with the NIST PSTAR mass
      electronic stopping-power ratio (S/rho)_SiC/(S/rho)_water at ``E_MeV``
      (Bragg-Gray, no secondary-electron transport; 0.823 at 62 MeV, 0.815 at 30 MeV,
      0.831 at 150 MeV). Flat in depth.
    - ``fluence_rate_cm2_s`` given (opt-in): proton fluence rate (cm^-2 s^-1); then
      G(x) = Phi' L_e,SiC(E(x)) / E_pair with E(x) the residual CSDA energy at depth x
      (NIST PSTAR Si + C, Bragg additivity) and G = 0 beyond the CSDA range. This
      carries the full energy dependence (L_SiC is ~3.4x larger at 30 MeV than at
      150 MeV) and the depth dependence for low-energy protons. ``dose_rate_Gy_s`` and
      ``dose_medium`` are ignored.

    For the two opt-in modes ``E_MeV`` may be anywhere in the PSTAR grid
    (1 keV - 10 GeV).

    Parameters
    ----------
    x_cm : array_like
        Depth positions in SiC (cm), relative to detector entrance.
    E_MeV : float
        Proton kinetic energy at the detector entrance (MeV).
    dose_rate_Gy_s : float
        Dose rate at detector entrance (Gy/s), in the medium ``dose_medium``.
    dose_medium : {"sic", "water"}
        Medium of ``dose_rate_Gy_s``.
    fluence_rate_cm2_s : float, optional
        Proton fluence rate; selects the fluence mode.

    Returns
    -------
    G : ndarray
        Generation rate (cm^-3 s^-1) at each position.
    """
    x = np.asarray(x_cm, dtype=float)

    if fluence_rate_cm2_s is not None:
        E = residual_energy_MeV("proton", "sic", float(E_MeV), x)
        L = np.zeros_like(E)
        alive = E > 1e-3  # PSTAR grid starts at 1 keV
        L[alive] = linear_stopping("proton", "sic", E[alive], "electronic")  # MeV/cm
        return float(fluence_rate_cm2_s) * L * 1e6 / E_PAIR_SIC_EV

    if dose_medium == "water":
        ratio = float(
            mass_stopping("proton", "sic", float(E_MeV))
            / mass_stopping("proton", "water", float(E_MeV))
        )
        return np.full_like(
            x, dose_rate_to_generation(dose_rate_Gy_s * ratio), dtype=float
        )
    if dose_medium != "sic":
        raise ValueError(f"dose_medium must be 'sic' or 'water' (got {dose_medium!r})")

    # Default: SiC dose rate. The range lookup only validates E_MeV (raises outside the
    # 30-150 MeV table: the scope of the flat profile, see the module docstring; PSTAR
    # itself covers 1 keV - 10 GeV); the profile is the flat entrance-region value.
    if E_MeV not in PROTON_RANGE_WATER_MM:
        _interpolate_proton_range_water(E_MeV)

    G_entrance = dose_rate_to_generation(dose_rate_Gy_s)
    return np.full_like(x, G_entrance, dtype=float)


def confine_generation_to_epi(device_info, depth, generation):
    """Zero ``generation`` beyond the epi end when the device has an n+ layer.

    Without an n+ layer (``device_info["nplus_thickness_cm"]`` None or absent) the mesh
    ends at the epi end and ``generation`` is returned unchanged (same object). With an
    n+ layer the deposit in the n+ is kept only if ``device_info["nplus_generation"]``
    is True. Default False keeps the historical normalisation (the part of the track
    beyond the epi end is not counted), which is a modelling choice, not a physics claim:
    whether the alpha deposit in the n+ substrate should count is an open question.

    ``depth`` is the node coordinate along the stack (x in 1D, y in 2D), same order as
    ``generation``. The epi-end node itself is kept.
    """
    if device_info.get("nplus_thickness_cm") is None or device_info.get(
        "nplus_generation", False
    ):
        return generation
    out = np.array(generation, dtype=float, copy=True)
    out[np.asarray(depth, dtype=float) > device_info["epi_end"] + 1e-12] = 0.0
    return out
