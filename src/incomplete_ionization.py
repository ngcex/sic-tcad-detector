"""Incomplete ionization model for Al acceptors in 4H-SiC.

Al acceptors in 4H-SiC have a deep ionization energy (E_A ~ 220 meV),
far above kT ~ 26 meV at 300K. The simple Gibbs formula predicts < 0.01%
ionization, but experimental data shows 10-30% ionization at N_A ~ 1e19
due to impurity band formation and screening at high doping.

This module implements a hybrid model:
  - For N_A < 1e18 cm^-3: Gibbs formula with doping-dependent E_A correction
  - For N_A >= 1e18 cm^-3: Empirical model calibrated to literature values

The transition between regimes uses smooth interpolation to avoid
discontinuities in downstream calculations.

All units CGS (cm^-3, eV, K) per devsim convention.

References:
    - TU Wien Ayalew thesis, Section 3.6
    - Pensl & Choyke, Physica B 185 (1993) 264
    - Kimoto & Cooper, "Fundamentals of Silicon Carbide Technology"
"""

import numpy as np


def _gibbs_ionization(N_A, T, E_A, g_A, k_B):
    """Simplified Gibbs distribution for ionized acceptor fraction.

    f_A = 1 / (1 + g_A * exp(E_A / (k_B * T)))

    Valid in the low-doping regime where impurity band effects are negligible.
    """
    kT = k_B * T
    if kT <= 0:
        return 0.0
    return 1.0 / (1.0 + g_A * np.exp(E_A / kT))


def _doping_dependent_EA(N_A, E_A0=0.220):
    """Compute doping-dependent ionization energy.

    E_A(N_A) = E_A0 - alpha * N_A^(1/3)

    At high doping, screening and band-tailing reduce the effective
    ionization energy. The coefficient alpha is calibrated so that
    the Gibbs formula produces physically reasonable ionization
    fractions in the moderate-doping regime.

    Parameters
    ----------
    N_A : float
        Total acceptor concentration (cm^-3).
    E_A0 : float
        Isolated acceptor ionization energy (eV).

    Returns
    -------
    E_A_eff : float
        Effective ionization energy (eV), clamped >= 0.
    """
    alpha = 3.4e-8  # eV*cm, calibrated coefficient
    E_A_eff = E_A0 - alpha * N_A ** (1.0 / 3.0)
    return max(E_A_eff, 0.0)


def _empirical_high_doping_ionization(N_A, T, T_ref=300.0):
    """Empirical ionization fraction for high-doping regime (N_A >= 1e18).

    At very high doping (N_A > 1e18 cm^-3), impurity band formation
    causes ionization fractions much higher than the Gibbs formula predicts.
    Literature reports ~10-30% ionization at N_A ~ 1e19, 300K.

    This empirical model interpolates between:
      - f ~ 0.05 at N_A = 1e18
      - f ~ 0.15 at N_A = 1e19
      - f ~ 0.30 at N_A = 1e20

    Temperature scaling: ionization increases with temperature due to
    thermal activation from the impurity band.

    Parameters
    ----------
    N_A : float
        Total acceptor concentration (cm^-3).
    T : float
        Temperature (K).
    T_ref : float
        Reference temperature for calibration (K).

    Returns
    -------
    f_A : float
        Ionized fraction (0 to 1).
    """
    # Log-linear interpolation in the high-doping regime
    # Calibrated: f(1e18) ~ 0.05, f(1e19) ~ 0.15, f(1e20) ~ 0.30
    log_NA = np.log10(max(N_A, 1e18))
    # Linear in log10(N_A): f = a + b * (log10(N_A) - 18)
    # f(18) = 0.05, f(19) = 0.15 => b = 0.10, a = 0.05
    f_ref = 0.05 + 0.10 * (log_NA - 18.0)
    f_ref = min(max(f_ref, 0.01), 0.95)  # clamp to physical range

    # Temperature scaling: mild increase with T
    # At 600K, ionization roughly doubles from impurity band broadening
    T_factor = 1.0 + 0.5 * (T - T_ref) / T_ref
    T_factor = max(T_factor, 0.5)  # don't let it go negative/zero

    f_A = f_ref * T_factor
    return min(f_A, 0.95)  # cap at 95%


def ionized_acceptor_fraction(N_A, T=300, E_A=None, g_A=4, k_B=8.617e-5):
    """Compute fraction of ionized Al acceptors in 4H-SiC.

    Hybrid model combining:
      1. Gibbs formula with doping-dependent E_A for N_A < 1e18 cm^-3
      2. Empirical impurity-band model for N_A >= 1e18 cm^-3
      3. Smooth interpolation in the transition region

    The 10-30% ionization at N_A ~ 1e19, 300K reported in literature
    cannot be reproduced by the Gibbs formula alone (which gives < 0.01%).
    At such high doping, impurity band formation, screening, and band-tailing
    dominate over the single-acceptor Gibbs statistics.

    Parameters
    ----------
    N_A : float
        Total Al acceptor concentration (cm^-3).
    T : float
        Temperature (K).
    E_A : float or None
        Acceptor ionization energy (eV). If None, uses doping-dependent
        correction starting from 0.220 eV.
    g_A : int
        Degeneracy factor (4 for Al in 4H-SiC).
    k_B : float
        Boltzmann constant (eV/K).

    Returns
    -------
    f_A : float
        Ionized fraction, between 0 and 1.

    References
    ----------
    TU Wien Ayalew thesis, Section 3.6.
    Kimoto & Cooper, "Fundamentals of Silicon Carbide Technology".
    """
    if N_A <= 0:
        return 0.0

    # Compute Gibbs-based fraction
    if E_A is None:
        E_A_eff = _doping_dependent_EA(N_A)
    else:
        E_A_eff = E_A
    f_gibbs = _gibbs_ionization(N_A, T, E_A_eff, g_A, k_B)

    # Compute empirical high-doping fraction
    f_empirical = _empirical_high_doping_ionization(N_A, T)

    # Smooth transition between regimes using a logistic blend
    # Below 1e17: pure Gibbs; above 1e19: pure empirical
    # Transition centered at ~1e18
    log_NA = np.log10(N_A)
    # Sigmoid centered at log10(N_A) = 18, width ~0.5 decades
    weight = 1.0 / (1.0 + np.exp(-(log_NA - 18.0) / 0.5))

    f_A = (1.0 - weight) * f_gibbs + weight * f_empirical

    # Ensure physical bounds
    return float(min(max(f_A, 0.0), 1.0))


def ionized_acceptor_concentration(N_A, T=300, **kwargs):
    """Compute ionized acceptor concentration N_A^-.

    N_A^- = N_A * ionized_acceptor_fraction(N_A, T, ...)

    This is the effective p-side doping used in all downstream
    calculations (built-in potential, depletion width, etc.).

    Parameters
    ----------
    N_A : float
        Total Al acceptor concentration (cm^-3).
    T : float
        Temperature (K).
    **kwargs
        Additional keyword arguments passed to ionized_acceptor_fraction.

    Returns
    -------
    N_A_ionized : float
        Ionized acceptor concentration (cm^-3).
    """
    f = ionized_acceptor_fraction(N_A, T, **kwargs)
    return N_A * f
