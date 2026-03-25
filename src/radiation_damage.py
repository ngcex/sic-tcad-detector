"""Radiation damage physics for 4H-SiC proton irradiation.

Pure-Python module (NO devsim dependency) implementing:
- Defect introduction rates (Z1/2, EH4, EH6/7) linear in fluence
- Carrier lifetime degradation (linear and logarithmic models)
- Carrier removal reducing effective doping
- NIEL hardness factor scaling across proton energies
- High-level compute_damaged_params interface
- Thermal annealing kinetics with per-defect Arrhenius recovery

All units in CGS (cm, cm^-3, eV, s) per project convention.
Defect parameters from Burin et al., arXiv:2407.16710 (2024).

NIEL hardness factors are placeholders pending SR-NIEL calculator lookup.
The module works with any NIEL table -- values are data, not code.

References:
    - Burin et al., arXiv:2407.16710 (2024): Defect introduction rates,
      trap energy levels, capture cross-sections for Z1/2, EH4, EH6/7
    - SR-NIEL Calculator (sr-niel.org): NIEL values for protons in SiC
"""

from dataclasses import dataclass
import logging

import numpy as np

from src.sic_material import SiC4H_Parameters

logger = logging.getLogger(__name__)

# Physical constants
_k_B_J = 1.3806e-23  # J/K, Boltzmann constant
_K_B_EV = 8.617e-5  # eV/K, Boltzmann constant in eV
_m0_kg = 9.109e-31  # kg, electron rest mass

# Material constants from SiC4H_Parameters (avoid duplicating values)
_sic = SiC4H_Parameters()
_M_E_DOS = _sic.m_e_dos  # 0.77 m0
_M_H_DOS = _sic.m_h_dos  # 1.0 m0


# ---------------------------------------------------------------------------
# NIEL hardness factors for protons in SiC
# kappa(E) = NIEL_proton(E) / NIEL_neutron(1 MeV)
# PLACEHOLDER values -- must be replaced with SR-NIEL calculator data
# before production use. See STATE.md blockers.
# ---------------------------------------------------------------------------
NIEL_HARDNESS_PROTON_SIC: dict[float, float] = {
    30: 0.50,  # placeholder -- obtain from SR-NIEL
    62: 0.35,  # placeholder -- obtain from SR-NIEL
    70: 0.33,  # placeholder -- obtain from SR-NIEL
    150: 0.22,  # placeholder -- obtain from SR-NIEL
}


@dataclass
class RadiationDamageParams:
    """Radiation damage constants for 4H-SiC.

    All introduction rates in cm^-1 (concentration per unit fluence).
    All capture cross-sections in cm^2.
    All energy levels in eV (below Ec unless noted).

    Reference energy: 1 MeV neutron equivalent (neq).
    Scale to proton energies via NIEL hardness factors.

    Source: Burin et al., arXiv:2407.16710 (2024), Table I.
    """

    # --- Z1/2 center (carbon vacancy) ---
    eta_Z12: float = 5.0  # cm^-1, introduction rate
    E_Z12: float = 0.67  # eV below Ec
    sigma_n_Z12: float = 2e-14  # cm^2, electron capture cross-section
    sigma_p_Z12: float = 3.5e-14  # cm^2, hole capture cross-section
    type_Z12: str = "acceptor"

    # --- EH6/7 center ---
    eta_EH67: float = 1.6  # cm^-1, introduction rate
    E_EH67: float = 1.60  # eV below Ec
    sigma_n_EH67: float = 9e-12  # cm^2, electron capture cross-section
    sigma_p_EH67: float = 3.8e-14  # cm^2, hole capture cross-section
    type_EH67: str = "donor"

    # --- EH4 center ---
    eta_EH4: float = 2.4  # cm^-1, introduction rate
    E_EH4: float = 1.03  # eV below Ec
    sigma_n_EH4: float = 5e-13  # cm^2, electron capture cross-section
    sigma_p_EH4: float = 5e-14  # cm^2, hole capture cross-section
    type_EH4: str = "acceptor"

    # --- Carrier removal ---
    eta_removal: float = 5.0  # cm^-1, effective carrier removal rate

    # --- Provenance ---
    source: str = "Burin et al., arXiv:2407.16710 (2024)"
    reference_particle: str = "1 MeV neutron equivalent"

    def __post_init__(self):
        """Validate that all physical parameters are positive."""
        for name in ("eta_Z12", "eta_EH67", "eta_EH4", "eta_removal"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive, got {getattr(self, name)}")
        for name in (
            "sigma_n_Z12",
            "sigma_p_Z12",
            "sigma_n_EH67",
            "sigma_p_EH67",
            "sigma_n_EH4",
            "sigma_p_EH4",
        ):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive, got {getattr(self, name)}")


@dataclass
class AnnealingParams:
    """Thermal annealing kinetics parameters for each defect type in 4H-SiC.

    Each defect has an activation energy (E_a) for thermally-driven recovery
    and an attempt frequency (nu_0). The recovery fraction follows first-order
    Arrhenius kinetics: f = 1 - exp(-nu_0 * t * exp(-E_a / (k_B * T))).

    Activation energies are calibrated so that:
    - Z1/2 is thermally stable below ~1000C on hour timescales (E_a=4.5 eV)
    - EH6/7 is stable below ~1200C (E_a=3.2 eV)
    - EH4 anneals significantly at ~400-600C (E_a=1.8 eV)

    See 17-RESEARCH.md open question 4 for Z1/2 calibration rationale:
    lower E_a values (e.g., 3.5 eV) give unwanted annealing at 1000C.

    Confidence: Z1/2 E_a HIGH (calibrated to stability data); EH4, EH6/7
    E_a MEDIUM-LOW (estimated from thermal stability observations, not
    directly measured via isothermal Arrhenius analysis).

    Source: Estimated from thermal stability observations in literature.
        Z1/2 migration barrier 3.7-4.2 eV (PRB 86, 075205), effective
        E_a raised to 4.5 eV for recombination kinetics. EH4 ~1.8 eV
        from 400-600C annealing onset. EH6/7 ~3.2 eV from 1200C stability.
        nu_0 = 1e13 Hz (standard solid-state Debye frequency).
    """

    # Z1/2: Carbon vacancy -- extremely thermally stable
    E_a_Z12: float = 4.5  # eV, activation energy for annealing
    nu_0_Z12: float = 1e13  # Hz, attempt frequency

    # EH6/7: Also very stable (related to carbon vacancy configurations)
    E_a_EH67: float = 3.2  # eV
    nu_0_EH67: float = 1e13  # Hz

    # EH4: Moderate stability, anneals at ~400-600C
    E_a_EH4: float = 1.8  # eV
    nu_0_EH4: float = 1e13  # Hz

    # Provenance
    source: str = (
        "Estimated from thermal stability observations. "
        "Z1/2: PRB 86, 075205 (migration barrier ~3.7-4.2 eV, "
        "effective E_a=4.5 eV for recombination kinetics). "
        "EH4/EH6/7: stability onset temperatures. "
        "nu_0=1e13 Hz (standard Debye frequency). "
        "Confidence: Z1/2 HIGH, EH4/EH6/7 MEDIUM-LOW."
    )

    def __post_init__(self):
        """Validate that all physical parameters are positive."""
        for name in ("E_a_Z12", "E_a_EH67", "E_a_EH4"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive, got {getattr(self, name)}")
        for name in ("nu_0_Z12", "nu_0_EH67", "nu_0_EH4"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive, got {getattr(self, name)}")


# ---------------------------------------------------------------------------
# Pure functions -- all stateless, taking explicit parameters
# ---------------------------------------------------------------------------


def defect_concentration(eta: float, fluence: float) -> float:
    """Compute defect concentration from introduction rate and fluence.

    N_defect = eta * Phi [cm^-3].

    Parameters
    ----------
    eta : float
        Introduction rate (cm^-1).
    fluence : float
        Fluence (neq/cm^2).

    Returns
    -------
    float
        Defect concentration (cm^-3).
    """
    return eta * fluence


def defect_concentrations(params: RadiationDamageParams, fluence_neq: float) -> dict:
    """Compute defect concentrations for all three defect types.

    Parameters
    ----------
    params : RadiationDamageParams
        Damage constants (introduction rates in cm^-1).
    fluence_neq : float
        Fluence in 1-MeV neutron equivalent (neq/cm^2).

    Returns
    -------
    dict
        Keys: N_Z12, N_EH67, N_EH4 (cm^-3).
    """
    return {
        "N_Z12": defect_concentration(params.eta_Z12, fluence_neq),
        "N_EH67": defect_concentration(params.eta_EH67, fluence_neq),
        "N_EH4": defect_concentration(params.eta_EH4, fluence_neq),
    }


def degraded_lifetime(
    tau_0: float,
    K_tau: float,
    fluence: float,
    model: str = "linear",
    alpha: float = 0.8,
) -> float:
    """Compute radiation-degraded carrier lifetime.

    Parameters
    ----------
    tau_0 : float
        Pristine (pre-irradiation) lifetime (s).
    K_tau : float
        Lifetime damage constant (cm^2/s).
    fluence : float
        Fluence (neq/cm^2).
    model : str
        "linear" or "logarithmic".
    alpha : float
        Exponent for logarithmic model (default 0.8, empirical).

    Returns
    -------
    float
        Degraded lifetime (s). Always >= 1e-15.

    Raises
    ------
    ValueError
        If model is not "linear" or "logarithmic".

    Notes
    -----
    Linear model: 1/tau = 1/tau_0 + K_tau * Phi
    Logarithmic model: tau = tau_0 / (1 + K_tau * tau_0 * Phi)^alpha
        More gradual saturation at high fluence.
    """
    if model == "linear":
        inv_tau = 1.0 / tau_0 + K_tau * fluence
        return max(1.0 / inv_tau, 1e-15)
    elif model == "logarithmic":
        factor = (1.0 + K_tau * tau_0 * fluence) ** alpha
        return max(tau_0 / factor, 1e-15)
    else:
        raise ValueError(
            f"Unknown lifetime model '{model}'. Use 'linear' or 'logarithmic'."
        )


def effective_doping(N_D: float, eta: float, fluence: float) -> float:
    """Compute effective doping after carrier removal.

    N_eff = max(N_D - eta * Phi, 0) [cm^-3].
    Floor at zero prevents unphysical negative doping.

    Parameters
    ----------
    N_D : float
        Original donor concentration (cm^-3).
    eta : float
        Carrier removal rate (cm^-1).
    fluence : float
        Fluence (neq/cm^2).

    Returns
    -------
    float
        Effective doping (cm^-3), >= 0.
    """
    return max(N_D - eta * fluence, 0.0)


def apply_carrier_removal(
    N_D_profile: np.ndarray, eta_removal: float, fluence_neq: float, floor: float = 0.0
) -> np.ndarray:
    """Apply carrier removal position-dependently to doping profile.

    N_D_damaged(x) = max(N_D(x) - eta * Phi, floor)

    Parameters
    ----------
    N_D_profile : array
        Original donor concentration at each node (cm^-3).
    eta_removal : float
        Carrier removal rate (cm^-1).
    fluence_neq : float
        Fluence in 1 MeV neutron equivalent (neq/cm^2).
    floor : float
        Minimum N_D value (cm^-3). Default 0.

    Returns
    -------
    np.ndarray
        Damaged doping profile (cm^-3).
    """
    removal = eta_removal * fluence_neq
    result = np.maximum(N_D_profile - removal, floor)

    # Warn if any position approaches compensation
    if np.any(result < 1e10):
        logger.warning(
            "Effective doping below 1e10 cm^-3 at %d position(s) -- "
            "approaching full compensation. Fluence_neq=%.2e",
            int(np.sum(result < 1e10)),
            fluence_neq,
        )

    return result


def compute_K_tau(
    params: RadiationDamageParams, carrier: str = "electron", T: float = 300.0
) -> float:
    """Compute lifetime damage constant from defect capture cross-sections.

    K_tau = sum_i (eta_i * sigma_i * v_th)

    Units: [cm^-1 * cm^2 * cm/s] = [cm^2/s]
    So 1/tau = 1/tau_0 + K_tau * Phi has units [1/s] = [1/s] + [cm^2/s * 1/cm^2].

    Parameters
    ----------
    params : RadiationDamageParams
        Defect parameters.
    carrier : str
        "electron" or "hole".
    T : float
        Temperature (K). Default 300.

    Returns
    -------
    float
        K_tau (cm^2/s).

    Raises
    ------
    ValueError
        If carrier is not "electron" or "hole".
    """
    if carrier == "electron":
        m_eff = _M_E_DOS * _m0_kg
        v_th = np.sqrt(3 * _k_B_J * T / m_eff) * 100  # m/s -> cm/s
        K_tau = (
            params.eta_Z12 * params.sigma_n_Z12 * v_th
            + params.eta_EH67 * params.sigma_n_EH67 * v_th
            + params.eta_EH4 * params.sigma_n_EH4 * v_th
        )
    elif carrier == "hole":
        m_eff = _M_H_DOS * _m0_kg
        v_th = np.sqrt(3 * _k_B_J * T / m_eff) * 100  # m/s -> cm/s
        K_tau = (
            params.eta_Z12 * params.sigma_p_Z12 * v_th
            + params.eta_EH67 * params.sigma_p_EH67 * v_th
            + params.eta_EH4 * params.sigma_p_EH4 * v_th
        )
    else:
        raise ValueError(f"carrier must be 'electron' or 'hole', got '{carrier}'")
    return K_tau


def get_hardness_factor(energy_MeV: float, niel_table: dict | None = None) -> float:
    """Get NIEL hardness factor for a given proton energy.

    Uses linear interpolation between table entries for intermediate energies.

    Parameters
    ----------
    energy_MeV : float
        Proton energy (MeV).
    niel_table : dict, optional
        Mapping energy (MeV) -> hardness factor. Defaults to
        NIEL_HARDNESS_PROTON_SIC.

    Returns
    -------
    float
        Hardness factor (dimensionless).
    """
    if niel_table is None:
        niel_table = NIEL_HARDNESS_PROTON_SIC
    energies = sorted(niel_table.keys())
    factors = [niel_table[e] for e in energies]
    return float(np.interp(energy_MeV, energies, factors))


def scale_to_proton_energy(
    damage_constant: float, energy_MeV: float, niel_table: dict | None = None
) -> float:
    """Scale a 1-MeV-neq damage constant to a specific proton energy.

    scaled = damage_constant * kappa(energy)

    Parameters
    ----------
    damage_constant : float
        Damage constant calibrated to 1 MeV neutron equivalent.
    energy_MeV : float
        Proton energy (MeV).
    niel_table : dict, optional
        Custom NIEL table. Defaults to NIEL_HARDNESS_PROTON_SIC.

    Returns
    -------
    float
        Scaled damage constant.
    """
    kappa = get_hardness_factor(energy_MeV, niel_table)
    return damage_constant * kappa


def compute_phi_crit(
    N_D_profile: np.ndarray,
    eta_removal: float = 5.0,
    energy_MeV: float = 62.0,
    niel_table: dict | None = None,
) -> dict:
    """Compute critical fluence for full carrier compensation.

    Phi_crit is the proton fluence at which the minimum doping position
    reaches full compensation (N_eff = 0). Beyond Phi_crit, the device
    cannot sustain a depletion region and the solver will diverge.

    Phi_crit_neq = min(N_D where N_D > 0) / eta_removal
    Phi_crit_proton = Phi_crit_neq / kappa(energy)

    Parameters
    ----------
    N_D_profile : np.ndarray
        Donor concentration profile (cm^-3). Only positive values are
        considered (p-side zeros are ignored).
    eta_removal : float
        Carrier removal rate (cm^-1). Default 5.0.
    energy_MeV : float
        Proton energy (MeV). Default 62.0.
    niel_table : dict, optional
        Custom NIEL hardness table. Defaults to NIEL_HARDNESS_PROTON_SIC.

    Returns
    -------
    dict
        Keys:
        - phi_crit_proton : float, critical proton fluence (protons/cm^2)
        - phi_crit_neq : float, critical 1-MeV-neq fluence (neq/cm^2)
        - N_D_min : float, minimum positive doping in profile (cm^-3)
        - kappa : float, NIEL hardness factor used

    Raises
    ------
    ValueError
        If N_D_profile has no positive values.
    """
    positive_mask = N_D_profile > 0
    if not np.any(positive_mask):
        raise ValueError("N_D_profile has no positive values; cannot compute Phi_crit")

    N_D_min = float(np.min(N_D_profile[positive_mask]))
    kappa = get_hardness_factor(energy_MeV, niel_table)

    phi_crit_neq = N_D_min / eta_removal
    phi_crit_proton = phi_crit_neq / kappa

    logger.info(
        "Phi_crit: %.3e protons/cm^2 (%.3e neq/cm^2) at %g MeV "
        "(N_D_min=%.3e, eta=%.1f, kappa=%.3f)",
        phi_crit_proton,
        phi_crit_neq,
        energy_MeV,
        N_D_min,
        eta_removal,
        kappa,
    )

    return {
        "phi_crit_proton": phi_crit_proton,
        "phi_crit_neq": phi_crit_neq,
        "N_D_min": N_D_min,
        "kappa": kappa,
    }


def annealing_fraction(
    T: float,
    t: float,
    E_a: float,
    nu_0: float = 1e13,
) -> float:
    """Compute fraction of defects removed by thermal annealing.

    First-order Arrhenius kinetics:
        N(t) = N_0 * exp(-k * t)
        k = nu_0 * exp(-E_a / (k_B * T))
        f = 1 - N(t)/N_0 = 1 - exp(-k * t)

    Parameters
    ----------
    T : float
        Annealing temperature (K). Must be > 0.
    t : float
        Annealing time (s). Returns 0.0 if t <= 0.
    E_a : float
        Activation energy (eV).
    nu_0 : float
        Attempt frequency (Hz). Default 1e13 (Debye frequency).

    Returns
    -------
    float
        Recovery fraction in [0, 1]. 0 = no recovery, 1 = full recovery.
    """
    if t <= 0:
        return 0.0
    rate = nu_0 * np.exp(-E_a / (_K_B_EV * T))
    exponent = rate * t
    # Clip to avoid overflow in exp for very large exponents
    if exponent > 700:
        return 1.0
    return float(1.0 - np.exp(-exponent))


def defect_recovery_fractions(
    T_anneal: float,
    t_anneal: float,
    anneal_params: AnnealingParams | None = None,
) -> dict:
    """Compute recovery fraction for each defect type at given thermal treatment.

    Parameters
    ----------
    T_anneal : float
        Annealing temperature (K).
    t_anneal : float
        Annealing time (s).
    anneal_params : AnnealingParams, optional
        Per-defect activation energies and attempt frequencies.
        Defaults to AnnealingParams().

    Returns
    -------
    dict
        Keys: f_Z12, f_EH67, f_EH4 -- each in [0, 1].
    """
    if anneal_params is None:
        anneal_params = AnnealingParams()
    return {
        "f_Z12": annealing_fraction(
            T_anneal, t_anneal, anneal_params.E_a_Z12, anneal_params.nu_0_Z12
        ),
        "f_EH67": annealing_fraction(
            T_anneal, t_anneal, anneal_params.E_a_EH67, anneal_params.nu_0_EH67
        ),
        "f_EH4": annealing_fraction(
            T_anneal, t_anneal, anneal_params.E_a_EH4, anneal_params.nu_0_EH4
        ),
    }


def compute_damaged_params(
    pristine_tau_n: float,
    pristine_tau_p: float,
    N_D_profile: np.ndarray,
    fluence: float,
    energy_MeV: float = 62.0,
    damage_params: RadiationDamageParams | None = None,
    lifetime_model: str = "linear",
    T: float = 300.0,
) -> dict:
    """Compute all radiation-degraded device parameters for a given fluence.

    CRITICAL: Short-circuits at fluence <= 0 -- returns pristine values
    with NO arithmetic operations (regression safety per Pitfall 1).

    Parameters
    ----------
    pristine_tau_n : float
        Pre-irradiation electron lifetime (s).
    pristine_tau_p : float
        Pre-irradiation hole lifetime (s).
    N_D_profile : np.ndarray
        Original donor concentration profile (cm^-3).
    fluence : float
        Proton fluence (protons/cm^2).
    energy_MeV : float
        Proton energy (MeV). Default 62.0.
    damage_params : RadiationDamageParams, optional
        Damage constants. Defaults to RadiationDamageParams().
    lifetime_model : str
        "linear" or "logarithmic". Default "linear".
    T : float
        Temperature (K). Default 300.

    Returns
    -------
    dict
        Keys: tau_n, tau_p, N_D_profile, N_Z12, N_EH67, N_EH4,
        fluence, fluence_neq, energy_MeV, lifetime_model.
    """
    # Short-circuit: zero or negative fluence returns pristine values unchanged
    if fluence <= 0:
        return {
            "tau_n": pristine_tau_n,
            "tau_p": pristine_tau_p,
            "N_D_profile": N_D_profile,
            "N_Z12": 0.0,
            "N_EH67": 0.0,
            "N_EH4": 0.0,
            "fluence": 0.0,
            "fluence_neq": 0.0,
            "energy_MeV": energy_MeV,
            "lifetime_model": lifetime_model,
        }

    if damage_params is None:
        damage_params = RadiationDamageParams()

    # Convert proton fluence to 1-MeV-neq equivalent
    kappa = get_hardness_factor(energy_MeV)
    fluence_neq = fluence * kappa

    # Lifetime degradation
    K_tau_n = compute_K_tau(damage_params, carrier="electron", T=T)
    K_tau_p = compute_K_tau(damage_params, carrier="hole", T=T)
    tau_n = degraded_lifetime(
        pristine_tau_n, K_tau_n, fluence_neq, model=lifetime_model
    )
    tau_p = degraded_lifetime(
        pristine_tau_p, K_tau_p, fluence_neq, model=lifetime_model
    )

    # Carrier removal (position-dependent)
    N_D_damaged = apply_carrier_removal(
        N_D_profile, damage_params.eta_removal, fluence_neq
    )

    # Defect concentrations
    defects = defect_concentrations(damage_params, fluence_neq)

    return {
        "tau_n": tau_n,
        "tau_p": tau_p,
        "N_D_profile": N_D_damaged,
        "N_Z12": defects["N_Z12"],
        "N_EH67": defects["N_EH67"],
        "N_EH4": defects["N_EH4"],
        "fluence": fluence,
        "fluence_neq": fluence_neq,
        "energy_MeV": energy_MeV,
        "lifetime_model": lifetime_model,
    }


def compute_annealed_params(
    pristine_tau_n: float,
    pristine_tau_p: float,
    N_D_profile: np.ndarray,
    fluence: float,
    energy_MeV: float = 62.0,
    T_anneal: float | None = None,
    t_anneal: float | None = None,
    damage_params: RadiationDamageParams | None = None,
    anneal_params: AnnealingParams | None = None,
    lifetime_model: str = "linear",
    T_device: float = 300.0,
) -> dict:
    """Compute device parameters after irradiation + thermal annealing.

    Composes radiation damage with per-defect thermal recovery:
    1. Compute irradiated state via compute_damaged_params
    2. Apply per-defect Arrhenius recovery fractions
    3. Recompute lifetimes from reduced defect concentrations (not interpolation)
    4. Scale carrier removal recovery by Z1/2 fraction (Z1/2-dominated)

    If T_anneal or t_anneal is None, returns compute_damaged_params output
    unchanged (no annealing = pass-through).

    Parameters
    ----------
    pristine_tau_n : float
        Pre-irradiation electron lifetime (s).
    pristine_tau_p : float
        Pre-irradiation hole lifetime (s).
    N_D_profile : np.ndarray
        Original donor concentration profile (cm^-3).
    fluence : float
        Proton fluence (protons/cm^2).
    energy_MeV : float
        Proton energy (MeV). Default 62.0.
    T_anneal : float, optional
        Annealing temperature (K). None = no annealing.
    t_anneal : float, optional
        Annealing time (s). None = no annealing.
    damage_params : RadiationDamageParams, optional
        Radiation damage constants. Defaults to RadiationDamageParams().
    anneal_params : AnnealingParams, optional
        Annealing kinetics parameters. Defaults to AnnealingParams().
    lifetime_model : str
        "linear" or "logarithmic". Default "linear".
    T_device : float
        Device operating temperature (K) for v_th in K_tau. Default 300.

    Returns
    -------
    dict
        Same keys as compute_damaged_params plus: f_Z12, f_EH67, f_EH4,
        T_anneal, t_anneal. Defect concentrations and lifetimes reflect
        post-anneal recovery.
    """
    # Pass-through: no annealing if T_anneal or t_anneal not specified
    if T_anneal is None or t_anneal is None:
        return compute_damaged_params(
            pristine_tau_n,
            pristine_tau_p,
            N_D_profile,
            fluence,
            energy_MeV,
            damage_params,
            lifetime_model,
            T_device,
        )

    if damage_params is None:
        damage_params = RadiationDamageParams()

    # Step 1: Get irradiated state
    damaged = compute_damaged_params(
        pristine_tau_n,
        pristine_tau_p,
        N_D_profile,
        fluence,
        energy_MeV,
        damage_params,
        lifetime_model,
        T_device,
    )

    # Short-circuit: zero fluence means nothing to anneal
    if fluence <= 0:
        damaged["f_Z12"] = 0.0
        damaged["f_EH67"] = 0.0
        damaged["f_EH4"] = 0.0
        damaged["T_anneal"] = T_anneal
        damaged["t_anneal"] = t_anneal
        return damaged

    # Step 2: Per-defect recovery fractions
    fractions = defect_recovery_fractions(T_anneal, t_anneal, anneal_params)
    f_Z12 = fractions["f_Z12"]
    f_EH67 = fractions["f_EH67"]
    f_EH4 = fractions["f_EH4"]

    # Step 3: Reduce defect concentrations
    N_Z12_annealed = damaged["N_Z12"] * (1.0 - f_Z12)
    N_EH67_annealed = damaged["N_EH67"] * (1.0 - f_EH67)
    N_EH4_annealed = damaged["N_EH4"] * (1.0 - f_EH4)

    # Step 4: Recompute lifetimes from annealed defect concentrations
    # CRITICAL: Do NOT interpolate between pristine and damaged lifetimes
    # (Pitfall 5 from RESEARCH.md). Instead, recompute K_tau from reduced etas.
    # We compute K_tau directly rather than using dataclasses.replace on
    # RadiationDamageParams, because full recovery (f=1.0) would produce
    # eta=0 which fails RadiationDamageParams validation.
    eta_Z12_eff = damage_params.eta_Z12 * (1.0 - f_Z12)
    eta_EH67_eff = damage_params.eta_EH67 * (1.0 - f_EH67)
    eta_EH4_eff = damage_params.eta_EH4 * (1.0 - f_EH4)

    # Electron K_tau
    m_eff_e = _M_E_DOS * _m0_kg
    v_th_e = np.sqrt(3 * _k_B_J * T_device / m_eff_e) * 100  # cm/s
    K_tau_n = (
        eta_Z12_eff * damage_params.sigma_n_Z12 * v_th_e
        + eta_EH67_eff * damage_params.sigma_n_EH67 * v_th_e
        + eta_EH4_eff * damage_params.sigma_n_EH4 * v_th_e
    )

    # Hole K_tau
    m_eff_h = _M_H_DOS * _m0_kg
    v_th_h = np.sqrt(3 * _k_B_J * T_device / m_eff_h) * 100  # cm/s
    K_tau_p = (
        eta_Z12_eff * damage_params.sigma_p_Z12 * v_th_h
        + eta_EH67_eff * damage_params.sigma_p_EH67 * v_th_h
        + eta_EH4_eff * damage_params.sigma_p_EH4 * v_th_h
    )

    fluence_neq = damaged["fluence_neq"]
    tau_n_annealed = degraded_lifetime(
        pristine_tau_n, K_tau_n, fluence_neq, model=lifetime_model
    )
    tau_p_annealed = degraded_lifetime(
        pristine_tau_p, K_tau_p, fluence_neq, model=lifetime_model
    )

    # Step 5: Carrier removal recovery -- Z1/2-dominated
    # Scale eta_removal by (1 - f_Z12) since carrier removal tracks Z1/2
    eta_removal_annealed = damage_params.eta_removal * (1.0 - f_Z12)
    N_D_annealed = apply_carrier_removal(N_D_profile, eta_removal_annealed, fluence_neq)

    return {
        "tau_n": tau_n_annealed,
        "tau_p": tau_p_annealed,
        "N_D_profile": N_D_annealed,
        "N_Z12": N_Z12_annealed,
        "N_EH67": N_EH67_annealed,
        "N_EH4": N_EH4_annealed,
        "fluence": fluence,
        "fluence_neq": fluence_neq,
        "energy_MeV": energy_MeV,
        "lifetime_model": lifetime_model,
        "f_Z12": f_Z12,
        "f_EH67": f_EH67,
        "f_EH4": f_EH4,
        "T_anneal": T_anneal,
        "t_anneal": t_anneal,
    }
