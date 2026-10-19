"""Radiation generation rate profiles for 4H-SiC detector simulation.

Provides analytical models for carrier generation from:
- Alpha particles (Am-241, 5.486 MeV) with smooth spatial profile
- Proton Bragg peak beams (30-150 MeV therapeutic energies)
- Dose rate to generation rate conversion

All units in CGS (cm, cm^-3, eV, s) per devsim convention.

References:
    - NIST PSTAR database for proton CSDA ranges in water
    - Bortfeld 1997, Med Phys 24(12) for Bragg curve analytical model
    - Ioffe NSM Archive for 4H-SiC material properties
"""

import numpy as np
from scipy.special import erfc

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

from etna.core.sic_material import SiC4H_Parameters

_params = SiC4H_Parameters()
RHO_SIC = _params.rho  # g/cm^3, 4H-SiC density
E_PAIR_SIC_EV = _params.E_pair_eV  # eV, electron-hole pair creation energy in 4H-SiC

# NIST PSTAR CSDA ranges in water (mm)
PROTON_RANGE_WATER_MM = {
    30: 8.85,
    62: 31.0,
    70: 40.8,
    150: 157.7,
}


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


def alpha_depth_shape(x_cm, alpha_range_cm=15e-4, bragg_peak_ratio=None):
    """Unnormalized spatial shape of the alpha e-h pair generation profile.

    Approximately uniform (or, if `bragg_peak_ratio` is set, linearly
    rising) from 0 to ~0.8*range, then smooth erfc roll-off (width
    ~0.1*range) to avoid numerical ringing in the DD solver. Pure
    function of `x_cm` — no normalization, no assumption that `x_cm` is
    sorted or evenly spaced, safe to evaluate on unordered mesh node
    arrays (e.g. a 2D device's node coordinates).

    Parameters
    ----------
    x_cm : array_like
        Depth positions in SiC (cm), measured from detector entrance.
        Values are NOT required to be sorted.
    alpha_range_cm : float
        Alpha particle range in SiC (cm). Default: 15e-4 (15 um).
    bragg_peak_ratio : float or None
        If given, linearly ramps the flat region from 1.0 (entrance) to
        this ratio (just before the erfc roll-off), approximating the
        real Bragg-curve rise in stopping power near end-of-range instead
        of a flat box. Default None (flat box, original behavior).

        QUALITATIVE shape correction, not a calibrated physical value: a
        real peak/entrance stopping-power ratio for a 5.486 MeV alpha in
        4H-SiC requires ASTAR/SRIM data this project does not have
        offline access to (same constraint as the kappa(E) NIST PSTAR
        data gap — see `data/srim/README.md`). Order-of-magnitude
        literature reasoning (Bethe-Bloch shape) suggests ~2.5-3x, well
        below a therapeutic proton's Bragg ratio (~30-50x, per Bortfeld
        1997) since an Am-241 alpha starts already well up its own
        stopping-power curve. Do NOT treat any specific value as
        validated — if precision matters, source a real ASTAR-Si or
        SRIM-4H-SiC stopping power table first (audit finding Mj-5,
        PHYSICS_AUDIT_v5.md).

    Returns
    -------
    shape : ndarray
        Unnormalized profile shape, same length as x_cm. Peaks at 1.0
        (bragg_peak_ratio=None) or bragg_peak_ratio (otherwise); callers
        normalize as needed (to total deposited energy, to a peak
        generation rate, etc.).
    """
    x = np.asarray(x_cm, dtype=float)
    R = alpha_range_cm

    # Transition region parameters
    x_flat = 0.8 * R  # flat/rising region ends at 80% of range
    sigma = 0.1 * R  # roll-off width (10% of range)

    if bragg_peak_ratio is None:
        # Flat-box shape (default: see docstring on why a calibrated peak
        # shape isn't available offline).
        shape_factor = np.ones_like(x)
    else:
        # Linear ramp from 1.0 (entrance) to bragg_peak_ratio (end of the
        # flat/rising region, just before the erfc roll-off) — a
        # QUALITATIVE stand-in for the real (nonlinear) Bragg-Kleeman rise
        # in stopping power near end-of-range. See docstring: not calibrated.
        shape_factor = 1.0 + (bragg_peak_ratio - 1.0) * np.clip(x / x_flat, 0.0, 1.0)

    # Smooth roll-off to zero beyond the flat/rising region.
    # erfc((x - x_flat) / (sqrt(2) * sigma)) transitions from 2 to 0.
    # Divide by 2 to get transition from 1 to 0.
    rolloff = 0.5 * erfc((x - x_flat) / (np.sqrt(2) * sigma))
    return shape_factor * rolloff


def alpha_generation_profile(
    x_cm,
    alpha_range_cm=15e-4,
    total_energy_eV=5.486e6,
    E_pair_eV=E_PAIR_SIC_EV,
    bragg_peak_ratio=None,
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
        Alpha particle range in SiC (cm). Default: 15e-4 (15 um).
    total_energy_eV : float
        Total kinetic energy of alpha particle (eV). Default: 5.486e6 (Am-241).
    E_pair_eV : float
        Electron-hole pair creation energy (eV). Default: 7.8 (4H-SiC).
    bragg_peak_ratio : float or None
        See `alpha_depth_shape`. Default None (flat box, original
        behavior — kept as default; see that function's docstring for
        why a calibrated peak shape isn't available offline).

    Returns
    -------
    G : ndarray
        Generation profile (cm^-1): pairs per cm depth per incident alpha.
        Multiply by fluence rate (alphas/cm^2/s) to get cm^-3 s^-1.

    Notes
    -----
    The profile is normalized so that integral(G(x) dx) = total_energy_eV / E_pair_eV,
    giving the total number of e-h pairs generated per alpha particle,
    REGARDLESS of `bragg_peak_ratio` — a local peak/entrance ratio changes
    the SHAPE, not the total charge. This matters for depth-dependent
    quantities (track containment near end-of-range, depth-resolved
    collection efficiency) but is a second-order correction for
    integrated CCE in a fully-depleted, uniform-collection detector.
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

# Qualitative Bragg-curve parameters for a 5.5 MeV alpha (see alpha_path_stopping).
ALPHA_PEAK_TO_ENTRANCE = (
    2.3  # S(E_peak)/S(E0): ASTAR-Si-like ORDER OF MAGNITUDE, not data
)
ALPHA_E_PEAK_EV = (
    0.75e6  # eV, energy of the stopping-power maximum (order of magnitude)
)
ALPHA_STRAGGLE_FRAC = 0.02  # range straggling sigma / range


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

    shape="bragg": two-power-law stopping power S(E) ~ (E/E_pk)^-a above E_pk and
    ~ (E/E_pk)^0.5 below, with a fixed by ``peak_to_entrance`` = S(E_pk)/S(E0);
    s(E) = int_E^E0 dE'/S scaled to ``alpha_range_cm``; Gaussian range straggling
    (sigma = straggle_frac * range) smears the track end. QUALITATIVE: the 2.3 ratio
    and 0.75 MeV peak are literature orders of magnitude for alphas in Si, not a
    tabulated ASTAR/SRIM 4H-SiC curve (none is available offline; data/srim/README.md).
    After smoothing the peak/entrance ratio is ~2.1.
    shape="box": ``alpha_depth_shape`` (flat to 0.8 R, erfc roll-off) along the path,
    the profile ``run_cce`` uses at normal incidence.
    """
    s = np.asarray(s_cm, dtype=float)
    R = float(alpha_range_cm)
    if shape == "box":
        f = alpha_depth_shape(s, alpha_range_cm=R)
    elif shape == "bragg":
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
        raise ValueError(f"shape must be 'bragg' or 'box' (got {shape!r})")
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


def proton_generation_profile(x_cm, E_MeV, dose_rate_Gy_s=1.0):
    """Carrier generation rate profile for proton beams in SiC detector.

    For therapeutic proton energies (30-150 MeV), the CSDA range far
    exceeds the thin SiC detector (~10 um), so the generation profile
    within the detector is approximately flat (entrance dose region,
    well before the Bragg peak).

    Energy-independent by construction: at equal dose rate the profile is the same
    for every E_MeV, which is only validated against the supported range (audit-v7
    B-10). Curves "at 30/70/150 MeV" at equal dose rate are therefore identical.

    Parameters
    ----------
    x_cm : array_like
        Depth positions in SiC (cm), relative to detector entrance.
    E_MeV : float
        Proton kinetic energy (MeV). Supported range: 30-150 MeV.
        Values between tabulated energies are interpolated.
    dose_rate_Gy_s : float
        Dose rate at detector entrance (Gy/s). Default: 1.0.

    Returns
    -------
    G : ndarray
        Generation rate (cm^-3 s^-1) at each position.

    Notes
    -----
    Proton ranges in SiC are scaled from water via density ratio:
    R_SiC = R_water * (rho_water / rho_SiC).

    For all therapeutic energies, range >> detector thickness,
    so the profile is essentially uniform within the detector.
    """
    x = np.asarray(x_cm, dtype=float)

    # The range lookup only validates E_MeV (raises outside the tabulated energies). The
    # SiC range (R_water * rho_water / rho_SiC, first order) is not used: the profile
    # below is the flat entrance-region value for every energy (audit B-10, open:
    # energy independence of this profile).
    if E_MeV not in PROTON_RANGE_WATER_MM:
        _interpolate_proton_range_water(E_MeV)

    # Convert dose rate to generation rate
    G_entrance = dose_rate_to_generation(dose_rate_Gy_s)

    # For entrance region (x << range), profile is flat
    G = np.full_like(x, G_entrance, dtype=float)

    return G


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
