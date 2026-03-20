"""Analytical electrostatics formulas for p-n junction analysis.

Provides reference analytical solutions for built-in potential, depletion width,
and electric field profile. These serve as validation benchmarks for the
numerical devsim solver and as standalone tools for quick estimates.

All units CGS (cm, cm^-3, V, V/cm, F/cm) per devsim convention.

References:
    - Sze & Ng, "Physics of Semiconductor Devices", 3rd ed.
    - Neamen, "Semiconductor Physics and Devices"
    - PVEducation: https://www.pveducation.org/pvcdrom/pn-junctions/
"""

import numpy as np

# Physical constants (CGS)
EPS_0 = 8.854e-14  # F/cm, vacuum permittivity
Q = 1.602e-19  # C, elementary charge
K_B = 8.617e-5  # eV/K, Boltzmann constant


def built_in_potential(N_A_ionized, N_D, n_i, T=300, k_B=K_B):
    """Compute built-in potential for a p-n junction.

    V_bi = kT * ln(N_A^- * N_D / n_i^2)

    Since kT is in eV and q = 1 in eV units, V_bi comes out in Volts directly.

    IMPORTANT: Use IONIZED acceptor concentration (N_A^-), not total N_A.
    For 4H-SiC with N_A = 1e19 and ~15% ionization: N_A^- ~ 1.5e18.

    Parameters
    ----------
    N_A_ionized : float
        Ionized acceptor concentration (cm^-3). Must be > 0.
    N_D : float
        Donor concentration (cm^-3). Must be > 0.
    n_i : float
        Intrinsic carrier concentration (cm^-3). Must be > 0.
    T : float
        Temperature (K).
    k_B : float
        Boltzmann constant (eV/K).

    Returns
    -------
    V_bi : float
        Built-in potential (V).

    References
    ----------
    Sze & Ng, Ch. 2; Neamen, Ch. 7.
    """
    if N_A_ionized <= 0 or N_D <= 0 or n_i <= 0:
        raise ValueError("All concentrations must be positive")

    kT = k_B * T  # eV = V (since q=1 in eV units)
    V_bi = kT * np.log(N_A_ionized * N_D / n_i**2)
    return V_bi


def depletion_width(V_bi, V_applied, N_D, eps_r=9.7, epi_thickness=None):
    """Compute one-sided depletion width for asymmetric p+/n- junction.

    W = sqrt(2 * eps * (V_bi - V_applied) / (q * N_D))

    Valid when N_A >> N_D, so the depletion region extends almost entirely
    into the lightly-doped n- side.

    Parameters
    ----------
    V_bi : float
        Built-in potential (V).
    V_applied : float
        Applied bias voltage (V). Negative for reverse bias.
        V_bi - V_applied increases under reverse bias.
    N_D : float
        Donor concentration in n- region (cm^-3).
    eps_r : float
        Relative permittivity (9.7 for 4H-SiC).
    epi_thickness : float or None
        Epitaxial layer thickness (cm). If provided, W is clamped
        to this value (punch-through condition).

    Returns
    -------
    W : float
        Depletion width (cm).

    Notes
    -----
    For the Petringa device:
      - W(0V) ~ 1.7 um = 1.7e-4 cm
      - W(-10V) ~ 9.5 um (approaches punch-through)
      - W(-30V) ~ 9.73 um (essentially fully depleted)
    """
    eps = eps_r * EPS_0  # F/cm

    V_total = V_bi - V_applied  # increases for reverse bias (V_applied < 0)
    if V_total < 0:
        # Forward bias exceeding V_bi -- no depletion
        return 0.0

    W = np.sqrt(2.0 * eps * V_total / (Q * N_D))

    if epi_thickness is not None and W > epi_thickness:
        W = epi_thickness

    return W


def electric_field_profile(x_cm, V_bi, V_applied, N_D, eps_r=9.7):
    """Compute electric field profile across a one-sided p+/n- junction.

    In the depletion region (0 to W):
        E(x) = -q * N_D * (W - x) / eps

    This gives a triangular profile with maximum magnitude at the
    junction (x=0) and zero at the depletion edge (x=W).
    Outside the depletion region, E = 0.

    Parameters
    ----------
    x_cm : array_like
        Position array (cm). x=0 is the metallurgical junction.
    V_bi : float
        Built-in potential (V).
    V_applied : float
        Applied bias voltage (V). Negative for reverse bias.
    N_D : float
        Donor concentration (cm^-3).
    eps_r : float
        Relative permittivity.

    Returns
    -------
    E_field : ndarray
        Electric field array (V/cm). Negative values indicate field
        pointing from n-side toward p-side (opposing diffusion).
    """
    x = np.asarray(x_cm, dtype=float)
    eps = eps_r * EPS_0

    W = depletion_width(V_bi, V_applied, N_D, eps_r)

    E_field = np.zeros_like(x)
    in_depletion = (x >= 0) & (x <= W)
    E_field[in_depletion] = -Q * N_D * (W - x[in_depletion]) / eps

    return E_field


def depletion_width_vs_bias(V_biases, V_bi, N_D, eps_r=9.7, epi_thickness=None):
    """Compute depletion width for an array of bias voltages.

    Convenience function for generating W vs V curves.

    Parameters
    ----------
    V_biases : array_like
        Array of applied bias voltages (V).
    V_bi : float
        Built-in potential (V).
    N_D : float
        Donor concentration (cm^-3).
    eps_r : float
        Relative permittivity.
    epi_thickness : float or None
        Epitaxial layer thickness (cm) for punch-through clamping.

    Returns
    -------
    W_array : ndarray
        Array of depletion widths (cm), same length as V_biases.
    """
    V_biases = np.asarray(V_biases, dtype=float)
    W_array = np.array(
        [depletion_width(V_bi, V, N_D, eps_r, epi_thickness) for V in V_biases]
    )
    return W_array
