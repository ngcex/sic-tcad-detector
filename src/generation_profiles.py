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

from src.sic_material import SiC4H_Parameters

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
        Electron-hole pair creation energy (eV). Default: 8.4 (4H-SiC).

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


def alpha_generation_profile(
    x_cm, alpha_range_cm=15e-4, total_energy_eV=5.486e6, E_pair_eV=E_PAIR_SIC_EV
):
    """Spatial e-h pair generation profile for Am-241 alpha particles in SiC.

    Uses a smooth profile: approximately uniform from 0 to ~0.8*range,
    then smooth erfc roll-off (width ~0.1*range) to avoid numerical ringing
    in the DD solver.

    Parameters
    ----------
    x_cm : array_like
        Depth positions in SiC (cm), measured from detector entrance.
    alpha_range_cm : float
        Alpha particle range in SiC (cm). Default: 15e-4 (15 um).
    total_energy_eV : float
        Total kinetic energy of alpha particle (eV). Default: 5.486e6 (Am-241).
    E_pair_eV : float
        Electron-hole pair creation energy (eV). Default: 8.4 (4H-SiC).

    Returns
    -------
    G : ndarray
        Generation profile (cm^-1): pairs per cm depth per incident alpha.
        Multiply by fluence rate (alphas/cm^2/s) to get cm^-3 s^-1.

    Notes
    -----
    The profile is normalized so that integral(G(x) dx) = total_energy_eV / E_pair_eV,
    giving the total number of e-h pairs generated per alpha particle.
    """
    x = np.asarray(x_cm, dtype=float)
    R = alpha_range_cm

    # Transition region parameters
    x_flat = 0.8 * R  # flat region ends at 80% of range
    sigma = 0.1 * R  # roll-off width (10% of range)

    # Smooth profile: erfc roll-off from flat region
    # erfc((x - x_flat) / (sqrt(2) * sigma)) transitions from 2 to 0
    # Divide by 2 to get transition from 1 to 0
    profile = 0.5 * erfc((x - x_flat) / (np.sqrt(2) * sigma))

    # Normalize so integral = total pairs per alpha
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

    # Get range in water (mm), interpolating if needed
    if E_MeV in PROTON_RANGE_WATER_MM:
        R_water_mm = PROTON_RANGE_WATER_MM[E_MeV]
    else:
        R_water_mm = _interpolate_proton_range_water(E_MeV)

    # Scale to SiC via density ratio (first-order approximation)
    # R_SiC = R_water * (rho_water / rho_SiC)
    _R_sic_mm = R_water_mm * (1.0 / RHO_SIC)  # noqa: F841

    # Convert dose rate to generation rate
    G_entrance = dose_rate_to_generation(dose_rate_Gy_s)

    # For entrance region (x << range), profile is flat
    G = np.full_like(x, G_entrance, dtype=float)

    return G
