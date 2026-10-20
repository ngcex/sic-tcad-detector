"""Incomplete ionization of dopants in 4H-SiC: the LEGACY scalar Al fraction and the
OPT-IN carrier-dependent model (``IonizationModel``, audit-v7 A-02).

Two models live here.

1. LEGACY (the default, ``DeviceConfig(incomplete_ionization=False)``):
   ``ionized_acceptor_fraction`` / ``ionized_acceptor_concentration``, ONE scalar
   ionized fraction for the p+ Al (0.132 at 1e19, 300 K), donors fully ionized. It is
   kept bit-identical so no baseline moves. Known defect (audit-v7 A-02): its low-doping
   branch is the occupancy with E_F pinned at E_V, so the fraction RISES with N_A below
   ~1e18 (1.8e-4 at 1e15 where charge neutrality gives ~0.67). DeviceConfig rejects
   N_A < 1e17, which limits the exposure; the high-doping branch is an uncited empirical
   blend (comments below).

2. OPT-IN carrier-dependent model (``incomplete_ionization=True`` or an
   ``IonizationModel``), steady-state Boltzmann form of the Gibbs distribution (Burin et
   al., "TCAD Parameters for 4H-SiC: A Review", arXiv:2410.06798 v8, Eq. 82):
       N_D+ = sum_l w_l N_D / (1 + g_D n / n1_l),   n1_l = N_C exp(-dE_l / kT)
       N_A- = sum_l w_l N_A / (1 + g_A p / p1_l),   p1_l = N_V exp(-dE_l / kT)
   with the Pearson-Bardeen doping dependence dE(N) = dE0 - alpha N^(1/3) (Burin Eq. 85,
   N = the local density of that dopant, the "dop" reading of Tables XXIV/XXV), clamped
   at 0 (Eq. 85 can push a level into the band, Burin p.107). The ionized densities
   follow the LOCAL n and p inside the Poisson equation of every solver, so dopants in
   a depleted region are fully ionized and only the neutral regions freeze out; the
   site split uses Burin Eq. 88 (equal hexagonal / cubic shares). Parameters (all
   exposed in ``IonizationModel``):
     - N donor, Kagamihara et al., J. Appl. Phys. 96, 5601 (2004) (Burin Table XXV,
       FCCS, "dop"): hexagonal 70.9 meV, alpha 3.38e-5 meV cm; cubic 123.7 meV,
       alpha 4.65e-5 meV cm; shares 1/2 and 1/2. Burin: this fit "fits the
       measurements best" for both sites.
     - Al acceptor, Achatz et al., Appl. Phys. Lett. 92, 072103 (2008) (Burin Table
       XXIV, "dop"): 220 meV, alpha 2.32e-5 meV cm (Burin: dE0 within 220 +- 20 meV,
       alpha (3.2 +- 1.7)e-5 meV cm across the literature). One effective level.
     - Degeneracy g_D = 2, g_A = 4, the values "commonly used" (Burin Sec. IX.B and its
       references; Ayalew thesis TU Wien 2004, Sec. 3.6).
     - N_C(T), N_V(T): etna.core.sic_material.effective_dos.
   ASSUMED, not from the sources: the same dE(N) law at every temperature (no T
   dependence of g, Burin Eq. 81 not used), no compensation term in N (Burin p.106
   discusses N = N_D - N_K), no impurity-band / deep Al level (Burin p.110 "d" level),
   quasi-static dopant occupancy (no capture/emission time constants: the small-signal
   C is the low-frequency one). At 300 K this gives a neutral ionized fraction of
   0.029 for Al at 1e19 (the legacy scalar says 0.132), 0.96 for N at 8.2e15, 0.87 at
   4.5e16 and 0.51 at 1e19.

The rest of this docstring describes the LEGACY scalar model.

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

from dataclasses import dataclass

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


# ---------------------------------------------------------------------------------------
# OPT-IN carrier-dependent model (audit-v7 A-02). See the module docstring for the
# equations, sources and assumptions. Everything below is numpy; the solvers (devsim,
# batched 1D, engine2d) all take the per-node fields from ``ionization_fields`` so they
# discretise the identical charge (parity), and the static equilibrium net doping from
# ``neutral_net_doping`` (ohmic contacts, Newton start).
# ---------------------------------------------------------------------------------------

MEV_CM = 1e-3  # meV cm -> eV cm (Burin quotes alpha in meV cm)


@dataclass(frozen=True)
class DopantLevel:
    """One ionization level: share ``fraction`` of the dopant sits on it, with
    dE(N) = max(dE0_eV - alpha_eVcm * N^(1/3), 0) (Burin Eq. 85)."""

    fraction: float
    dE0_eV: float
    alpha_eVcm: float

    def energy(self, N):
        N = np.asarray(N, dtype=float)
        return np.maximum(self.dE0_eV - self.alpha_eVcm * np.cbrt(N), 0.0)


# N donor, hexagonal and cubic sites, Kagamihara et al., J. Appl. Phys. 96, 5601 (2004),
# as tabulated in Burin et al. arXiv:2410.06798 Table XXV; equal site shares (Eq. 88).
NITROGEN_KAGAMIHARA_2004 = (
    DopantLevel(0.5, 70.9e-3, 3.38e-5 * MEV_CM),  # hexagonal
    DopantLevel(0.5, 123.7e-3, 4.65e-5 * MEV_CM),  # cubic
)
# Al acceptor, one effective level, Achatz et al., Appl. Phys. Lett. 92, 072103 (2008),
# as tabulated in Burin et al. Table XXIV.
ALUMINUM_ACHATZ_2008 = (DopantLevel(1.0, 220e-3, 2.32e-5 * MEV_CM),)


@dataclass(frozen=True)
class IonizationModel:
    """Parameters of the carrier-dependent incomplete-ionization model (hashable).

    Defaults: N (Kagamihara 2004, two sites) and Al (Achatz 2008), g_D = 2, g_A = 4.
    The level fractions of each species should sum to 1.
    """

    donor_levels: tuple = NITROGEN_KAGAMIHARA_2004
    g_D: float = 2.0
    acceptor_levels: tuple = ALUMINUM_ACHATZ_2008
    g_A: float = 4.0

    def __post_init__(self):
        for name in ("donor_levels", "acceptor_levels"):
            levels = tuple(getattr(self, name))
            if not levels:
                raise ValueError(f"IonizationModel.{name} must not be empty")
            if not all(isinstance(lv, DopantLevel) for lv in levels):
                raise TypeError(f"IonizationModel.{name} must be DopantLevel entries")
            if abs(sum(lv.fraction for lv in levels) - 1.0) > 1e-12:
                raise ValueError(f"IonizationModel.{name} fractions must sum to 1")
            if any(lv.fraction < 0 or lv.dE0_eV < 0 or lv.alpha_eVcm < 0 for lv in levels):
                raise ValueError(f"IonizationModel.{name}: negative parameter")
            object.__setattr__(self, name, levels)
        if not (self.g_D > 0 and self.g_A > 0):
            raise ValueError("IonizationModel degeneracies must be > 0")


DEFAULT_IONIZATION_MODEL = IonizationModel()


def resolve_ionization(spec):
    """``incomplete_ionization`` keyword -> IonizationModel or None (legacy model).

    None / False -> None; True -> DEFAULT_IONIZATION_MODEL; an IonizationModel as is.
    """
    if spec is None or spec is False:
        return None
    if spec is True:
        return DEFAULT_IONIZATION_MODEL
    if isinstance(spec, IonizationModel):
        return spec
    raise TypeError(
        "incomplete_ionization must be a bool, None or an IonizationModel "
        f"(got {type(spec).__name__})"
    )


@dataclass(frozen=True)
class IonizationFields:
    """Per-node inputs of the ionized charge, arrays of shape (L, *node_shape), cm^-3.

    ``d_dens[l]`` = share of the donors on level l, ``d_ref[l]`` = n1_l / g_D, so
    N_D+ = sum_l d_dens[l] / (1 + n / d_ref[l]); the same for acceptors with p.
    """

    d_dens: np.ndarray
    d_ref: np.ndarray
    a_dens: np.ndarray
    a_ref: np.ndarray


def _thermal_voltage(T):
    from etna.core.device import K_B_J, Q  # the solvers' constants (V_t = kT/q)

    return K_B_J * float(T) / Q


def ionization_fields(donors, acceptors_total, T, model=DEFAULT_IONIZATION_MODEL):
    """IonizationFields of a doping profile (TOTAL donors and acceptors, cm^-3) at T (K)."""
    from etna.core.sic_material import effective_dos

    donors = np.asarray(donors, dtype=float)
    acceptors_total = np.asarray(acceptors_total, dtype=float)
    NC, NV = effective_dos(float(T))
    V_t = _thermal_voltage(T)

    def build(levels, N, Neff, g):
        dens = np.stack([lv.fraction * N for lv in levels])
        ref = np.stack([Neff * np.exp(-lv.energy(N) / V_t) / g for lv in levels])
        return dens, ref

    d_dens, d_ref = build(model.donor_levels, donors, NC, model.g_D)
    a_dens, a_ref = build(model.acceptor_levels, acceptors_total, NV, model.g_A)
    return IonizationFields(d_dens, d_ref, a_dens, a_ref)


def ionized_densities(fields, n, p):
    """(N_D+, N_A-) in cm^-3 at carrier densities n, p (arrays of the node shape)."""
    n = np.asarray(n, dtype=float)
    p = np.asarray(p, dtype=float)
    Dp = (fields.d_dens / (1.0 + n / fields.d_ref)).sum(axis=0)
    Am = (fields.a_dens / (1.0 + p / fields.a_ref)).sum(axis=0)
    return Dp, Am


def neutral_state(fields, n_i, iterations=200):
    """Local charge neutrality at every node: (n, p, N_D+, N_A-), each of the node shape.

    Solves p - n + N_D+(n) - N_A-(p) = 0 with n p = n_i^2 by bisection in
    u = ln(n / n_i) (the residual is strictly decreasing in u), 200 halvings of a
    bracket of width < 600: u to machine precision, deterministic, the same numbers for
    every solver. This is the equilibrium state of a neutral region (ohmic contacts).
    """
    n_i = float(n_i)
    N_max = float(
        max(np.max(fields.d_dens.sum(axis=0)), np.max(fields.a_dens.sum(axis=0)), 1.0)
    )
    half = np.log(N_max / n_i) + 5.0  # |u| of any neutral state is below this
    shape = fields.d_dens.shape[1:]
    lo = np.full(shape, -half)
    hi = np.full(shape, half)

    def residual(u):
        n = n_i * np.exp(u)
        p = n_i * np.exp(-u)
        Dp, Am = ionized_densities(fields, n, p)
        return p - n + Dp - Am

    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        pos = residual(mid) > 0.0  # root lies above mid
        lo = np.where(pos, mid, lo)
        hi = np.where(pos, hi, mid)
    u = 0.5 * (lo + hi)
    n = n_i * np.exp(u)
    p = n_i * np.exp(-u)
    Dp, Am = ionized_densities(fields, n, p)
    return n, p, Dp, Am


def neutral_net_doping(fields, n_i):
    """Equilibrium ionized net doping N_D+ - N_A- of the neutral state at every node
    (cm^-3): the ``NetDoping`` of a device with the carrier-dependent model."""
    _, _, Dp, Am = neutral_state(fields, n_i)
    return Dp - Am


def neutral_ionized_acceptor_concentration(
    N_A, T, model=DEFAULT_IONIZATION_MODEL, N_D=0.0
):
    """Ionized acceptors (cm^-3) of a neutral region with N_A acceptors (and N_D donors)
    at T under the carrier-dependent model (the p+ anode's equilibrium value)."""
    from etna.core.sic_material import intrinsic_concentration, SiC4H_Parameters

    n_i = intrinsic_concentration(float(T), SiC4H_Parameters())[0]
    f = ionization_fields(np.array([N_D]), np.array([N_A]), T, model)
    return float(neutral_state(f, n_i)[3][0])
