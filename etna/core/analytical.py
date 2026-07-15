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


def full_depletion_voltage_graded(
    epi_thickness,
    N_D_bulk,
    N_D_junction,
    L_transition,
    V_bi,
    eps_r=9.7,
    n_d_uniform=None,
):
    """Reverse bias required to fully deplete a graded-doped epi layer.

    For a one-sided p+/n- junction the potential to deplete from the junction
    out to depth ``t`` with a depth-dependent donor density N_D(y) is the
    Poisson double integral:

        V_fd = (q / eps) * integral_0^t [ y * N_D(y) ] dy  -  V_bi

    The project's graded profile (set_graded_doping_2d) is

        N_D(y) = N_D_bulk + (N_D_junction - N_D_bulk) * exp(-y / L)

    for which the moment integral has the closed form

        integral_0^t y*N_D(y) dy =
            N_D_bulk * t^2 / 2
          + (N_D_junction - N_D_bulk) * ( L^2 - exp(-t/L) * (L*t + L^2) )

    AUDIT Mj-3 (v5): the uniform-N_D_bulk estimate (q*N_D_bulk*t^2/(2 eps) - V_bi)
    is NOT a valid gate for the graded sweep -- it is anti-conservative, because
    the near-junction doping (up to N_D_junction) is much higher than the bulk
    asymptote, so the real V_fd is LARGER. Using N_D_bulk would mark
    not-fully-depleted configs as depleted (gate fails open). This function
    integrates the actual profile.

    Parameters
    ----------
    epi_thickness : float
        Epitaxial layer thickness t (cm).
    N_D_bulk : float
        Deep-bulk donor asymptote (cm^-3).
    N_D_junction : float
        Near-junction donor concentration (cm^-3).
    L_transition : float
        Exponential grading length L (cm).
    V_bi : float
        Built-in potential (V).
    eps_r : float
        Relative permittivity (9.7 for 4H-SiC).
    n_d_uniform : float or None
        If given, ignore the graded profile and use a uniform donor density
        n_d_uniform, collapsing to V_fd = q*N_D*t^2/(2 eps) - V_bi. Provided
        as an explicit, documented fallback for uniform-doping devices.

    Returns
    -------
    V_fd : float
        Reverse bias magnitude (V, positive) needed to fully deplete the epi.
        Compare against the applied reverse-bias magnitude: fully depleted iff
        |V_bias| >= V_fd.
    """
    eps = eps_r * EPS_0  # F/cm
    t = epi_thickness
    if n_d_uniform is not None:
        moment = n_d_uniform * t**2 / 2.0
    else:
        L = L_transition
        moment = N_D_bulk * t**2 / 2.0 + (N_D_junction - N_D_bulk) * (
            L**2 - np.exp(-t / L) * (L * t + L**2)
        )
    return (Q / eps) * moment - V_bi


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
    For the reference device:
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
