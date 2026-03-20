"""4H-SiC material parameters for TCAD simulation.

All units in CGS (cm, cm^-3, F/cm, eV, s) per devsim convention.
Parameters sourced from literature with citations in field comments.

References:
    - Ioffe NSM Archive: https://www.ioffe.ru/SVA/NSM/Semicond/SiC/
    - TU Wien Ayalew thesis: https://www.iue.tuwien.ac.at/phd/ayalew/
    - Sze & Ng, "Physics of Semiconductor Devices", 3rd ed.
"""

from dataclasses import dataclass
import numpy as np


@dataclass
class SiC4H_Parameters:
    """4H-SiC material parameters at 300K.

    All values sourced from literature with citations.
    Units are CGS (cm, cm^-3, eV, s, F/cm) for devsim compatibility.
    """

    # --- Bandgap ---
    E_g: float = 3.26  # eV at 300K, Ioffe NSM Archive
    E_g_0: float = 3.265  # eV at 0K, Varshni extrapolation
    E_g_alpha: float = 6.5e-4  # eV/K, Varshni parameter (Ioffe NSM)
    E_g_beta: float = 1300.0  # K, Varshni parameter (Ioffe NSM)

    # --- Dielectric constant ---
    eps_r: float = 9.7  # relative permittivity, Ioffe NSM Archive

    # --- Effective masses (density of states) ---
    m_e_dos: float = 0.77  # m0, electron DOS effective mass (Ioffe NSM)
    m_h_dos: float = 1.0  # m0, hole DOS effective mass (Ioffe NSM, approx)
    M_c: int = 3  # equivalent conduction band minima (4H-SiC)

    # --- Effective density of states at 300K ---
    NC_300: float = 1.69e19  # cm^-3, computed: 2*Mc*(2pi*m_e*kT/h^2)^1.5
    NV_300: float = 2.49e19  # cm^-3, computed: 2*(2pi*m_h*kT/h^2)^1.5

    # --- Intrinsic carrier concentration at 300K ---
    # n_i = sqrt(NC*NV) * exp(-Eg/(2*kT))
    # ~5e-9 cm^-3 (18 orders below Si's 1e10)
    n_i_300: float = 5.0e-9  # cm^-3, Ioffe/TU Wien

    # --- Electron mobility: Caughey-Thomas model ---
    # Source: TU Wien Ayalew thesis, Table 3.5
    mu_n_max: float = 950.0  # cm^2/Vs at 300K, low-doping limit
    mu_n_min: float = 40.0  # cm^2/Vs, high-doping limit
    N_ref_n: float = 1.94e17  # cm^-3, reference doping concentration
    alpha_n: float = 0.61  # doping exponent

    # --- Hole mobility: Caughey-Thomas model ---
    # Source: TU Wien Ayalew thesis, Table 3.5
    mu_p_max: float = 125.0  # cm^2/Vs at 300K, low-doping limit
    mu_p_min: float = 15.9  # cm^2/Vs, high-doping limit
    N_ref_p: float = 1.76e19  # cm^-3, reference doping concentration
    alpha_p: float = 0.34  # doping exponent

    # --- SRH recombination lifetimes ---
    tau_n: float = 1.0e-9  # s, electron lifetime in p-type (Ioffe NSM)
    tau_p: float = 6.0e-7  # s, hole lifetime in n-type (Ioffe NSM)

    # --- Auger recombination coefficients ---
    C_n: float = 5.0e-31  # cm^6/s, electron Auger (Ioffe NSM)
    C_p: float = 2.0e-31  # cm^6/s, hole Auger (Ioffe NSM)

    # --- Radiative recombination ---
    B: float = 1.5e-12  # cm^3/s, radiative coefficient (estimation)

    # --- Al acceptor incomplete ionization ---
    E_A: float = 0.220  # eV, Al acceptor ionization energy (TU Wien)
    g_A: int = 4  # degeneracy factor for Al in SiC

    # --- N donor ionization energies ---
    E_D_hex: float = 0.050  # eV, nitrogen on hexagonal site (TU Wien)
    E_D_cub: float = 0.092  # eV, nitrogen on cubic site (TU Wien)


def compute_ni(T=300):
    """Compute intrinsic carrier concentration for 4H-SiC from first principles.

    Uses Varshni bandgap model and Boltzmann statistics for NC, NV.

    Parameters
    ----------
    T : float
        Temperature in Kelvin.

    Returns
    -------
    n_i : float
        Intrinsic carrier concentration (cm^-3).
    NC : float
        Effective density of states in conduction band (cm^-3).
    NV : float
        Effective density of states in valence band (cm^-3).
    E_g : float
        Bandgap at temperature T (eV).

    References
    ----------
    Ioffe NSM Archive, TU Wien Ayalew thesis.
    """
    k_B = 8.617e-5  # eV/K
    h = 6.626e-34  # J*s
    m0 = 9.109e-31  # kg

    # 4H-SiC parameters
    m_e = 0.77 * m0  # DOS effective mass, electrons
    m_h = 1.0 * m0  # DOS effective mass, holes
    M_c = 3  # conduction band minima

    # Varshni bandgap: E_g(T) = E_g(0) - alpha*T^2/(T+beta)
    E_g = 3.265 - 6.5e-4 * T**2 / (T + 1300.0)

    # Effective density of states
    # NC = 2 * Mc * (2*pi*m_e*kB_J*T / h^2)^(3/2)  [m^-3, then convert to cm^-3]
    kT_J = k_B * T * 1.602e-19  # eV -> J
    NC = 2 * M_c * (2 * np.pi * m_e * kT_J / h**2) ** 1.5 * 1e-6  # m^-3 -> cm^-3
    NV = 2 * (2 * np.pi * m_h * kT_J / h**2) ** 1.5 * 1e-6

    n_i = np.sqrt(NC * NV) * np.exp(-E_g / (2 * k_B * T))
    return n_i, NC, NV, E_g


def mobility_caughey_thomas(N_total, carrier="electron"):
    """Doping-dependent mobility using Caughey-Thomas model.

    mu(N) = mu_min + (mu_max - mu_min) / (1 + (N/N_ref)^alpha)

    Parameters
    ----------
    N_total : float
        Total doping concentration (cm^-3).
    carrier : str
        'electron' or 'hole'.

    Returns
    -------
    mu : float
        Mobility (cm^2/Vs).

    References
    ----------
    TU Wien Ayalew thesis, Table 3.5.
    """
    if carrier == "electron":
        mu_max, mu_min = 950.0, 40.0  # cm^2/Vs
        N_ref, alpha = 1.94e17, 0.61
    elif carrier == "hole":
        mu_max, mu_min = 125.0, 15.9  # cm^2/Vs
        N_ref, alpha = 1.76e19, 0.34
    else:
        raise ValueError(f"carrier must be 'electron' or 'hole', got '{carrier}'")

    mu = mu_min + (mu_max - mu_min) / (1.0 + (N_total / N_ref) ** alpha)
    return mu
