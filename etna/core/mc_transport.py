"""Fast condensed-history energy-deposition Monte Carlo for thin SiC sensitive volumes.

What it is: an event-by-event sampler of the energy a light ion (proton, alpha, C-12)
deposits in a rectangular sensitive volume (SV), with

- the MEAN energy loss tied to the bundled NIST PSTAR/ASTAR electronic stopping powers
  (SiC by Bragg additivity of NIST Si + C, ``etna.core.stopping_power``);
- energy-loss FLUCTUATIONS (straggling) from a line-by-line numba port of the Geant4
  11.3 models (``G4UniversalFluctuation`` / ``G4UrbanFluctuation`` / the Gaussian
  branch of ``G4IonFluctuations``; L. Urban, Geant4 Physics Reference Manual), which
  cover the thin-layer (Landau / Vavilov-like) and thick-layer (Bohr Gaussian) regimes;
- explicit DELTA RAYS above a cut T_cut (restricted energy loss below it, the Geant4
  scheme): 1/T^2 (1 - beta^2 T/T_max) spectrum, two-body emission angle, straight-line
  transport with the NIST ESTAR CSDA range (``etna.core.electron_range``), so delta-ray
  energy can escape the SV (or enter it from a SiC surround);
- adaptive stepping (energy-loss fraction per step) so particles that stop inside the SV
  (alpha Bragg peak) end at the CSDA range;
- per-event deterministic counter-based random streams (SplitMix64 seeded from
  (seed, event index)): results do not depend on the numba thread count.

Calibrated vs predicted: nothing here is fitted. Mean stopping and electron ranges are
NIST data (SiC composed); the fluctuation model is Geant4's (its parameters are Geant4's
defaults, quoted from the 11.3 source). C-12 stopping is MODELLED as z_eff^2 x PSTAR at
the same velocity (Ziegler effective charge as in G4ionEffectiveCharge, with the Fermi
velocity set to 1 v_Bohr; < 0.5 % effect on z_eff above 2 MeV/u, unverified below).

Deliberately NOT modelled (quantified in docs/USAGE.md and lab/mc_transport):
nuclear interactions / fragmentation, multiple Coulomb scattering of the primary
(straight tracks), electron scattering (delta rays go straight: an upper bound on
escape distance, ``delta_range_factor`` scales it), Yang's charge-exchange factor of
G4IonFluctuations, the Lindhard partition of nuclear stopping (nuclear loss slows the
ion but is not deposited as ionisation).

Units in the API: MeV (kinetic energies), keV (deposits), um (geometry).
"""

from __future__ import annotations

import functools
import math
from dataclasses import dataclass, field

import numpy as np
from numba import njit, prange

from etna.core import stopping_power as sp
from etna.core.electron_range import electron_csda_table

# --- physical constants (CODATA 2018 / PDG) -------------------------------------------
ME_C2 = 0.51099895  # MeV
R_E_CM = 2.8179403262e-13  # classical electron radius, cm
TWOPI_MC2_RCL2 = 2.0 * math.pi * ME_C2 * R_E_CM**2  # MeV cm^2 (Geant4 twopi_mc2_rcl2)
AMU_C2 = 931.49410242  # MeV
PROTON_MC2 = 938.27208816  # MeV
N_A = 6.02214076e23

PARTICLES = {
    # name: (rest mass MeV, charge z, mass number A, PDG code)
    "proton": (PROTON_MC2, 1, 1, 2212),
    "alpha": (3727.3794066, 2, 4, 1000020040),
    "c12": (12.0 * AMU_C2 - 6.0 * ME_C2, 6, 12, 1000060120),
}

# --- materials ------------------------------------------------------------------------
# Elemental NIST data used for the composed values: Z, A (g/mol), I (eV) as printed on
# the saved NIST composition pages (data/nist_star/README.md).
_ELEM = {"Si": (14, 28.0855, 173.0), "C": (6, 12.011, 81.0), "H": (1, 1.008, 19.2)}


@dataclass(frozen=True)
class Material:
    """Target medium of the SV (and of its optional surround)."""

    name: str
    density: float  # g/cm^3
    z_over_a: float  # mol/g
    I_eV: float  # mean excitation energy
    zeff_mass: float  # sum w_i Z_i (G4IonisParamMat::ComputeFluctModel 'Zeff')
    zeff_atom: float  # atom-density-weighted Z (G4IonisParamMat fZeff)
    star_name: str  # key for etna.core.stopping_power ("sic", "si", "water")
    estar_name: str  # key for etna.core.electron_range

    @property
    def electron_density_cm3(self):
        return self.density * N_A * self.z_over_a


def _sic():
    w = {"Si": sp.SIC_MASS_FRACTION_SI, "C": sp.SIC_MASS_FRACTION_C}
    zoa = {k: _ELEM[k][0] / _ELEM[k][1] for k in w}
    z_over_a = sum(w[k] * zoa[k] for k in w)
    # Bragg rule for I: ln I = sum w_i (Z/A)_i ln I_i / (Z/A)
    ln_i = sum(w[k] * zoa[k] * math.log(_ELEM[k][2]) for k in w) / z_over_a
    return Material(
        name="sic",
        density=sp.density("sic"),
        z_over_a=z_over_a,
        I_eV=math.exp(ln_i),
        zeff_mass=sum(w[k] * _ELEM[k][0] for k in w),
        zeff_atom=0.5 * (14 + 6),
        star_name="sic",
        estar_name="sic",
    )


@functools.lru_cache(maxsize=None)
def get_material(name):
    """``"sic"`` (4H-SiC, rho from SiC4H_Parameters), ``"si"`` or ``"water"``."""
    if name == "sic":
        return _sic()
    if name == "si":
        return Material("si", 2.33, 14 / 28.0855, 173.0, 14.0, 14.0, "si", "si")
    if name == "water":
        wh = 2 * 1.008 / (2 * 1.008 + 15.999)
        return Material(
            "water",
            1.0,
            0.55509,  # NIST ESTAR water Z/A
            75.0,
            wh * 1 + (1 - wh) * 8,
            (2 * 1 + 8) / 3.0,
            "water",
            "water",
        )
    raise ValueError(f"material must be 'sic', 'si' or 'water' (got {name!r})")


# --- projectile tables ----------------------------------------------------------------


def ion_effective_charge(particle, T_MeV, material):
    """Effective charge (units of e) as in Geant4 11.3 ``G4ionEffectiveCharge``.

    Protons: 1. Alpha: Ziegler's He formula. C-12: Ziegler heavy-ion formula with the
    material Fermi velocity set to 1 v_Bohr (Geant4 uses tabulated element values;
    the full charge is returned above 20 MeV x z per proton mass in either case).
    """
    mass, z, _, _ = PARTICLES[particle]
    T = np.asarray(T_MeV, dtype=float)
    if z <= 1:
        return np.ones_like(T)
    red = T * PROTON_MC2 / mass  # MeV, "reduced" energy
    out = np.full_like(T, float(z))
    low = red <= z * 20.0
    if not np.any(low):
        return out
    r = np.maximum(red[low], 1e-3)
    zm = material.zeff_atom
    if z == 2:
        c = (0.2865, 0.1266, -0.001429, 0.02402, -0.01135, 0.001475)
        Q = np.maximum(0.0, np.log(r * 1e3 * AMU_C2 / PROTON_MC2))
        x = c[0] + c[1] * Q + c[2] * Q**2 + c[3] * Q**3 + c[4] * Q**4 + c[5] * Q**5
        ex = np.where(x < 0.2, x * (1 - 0.5 * x), 1.0 - np.exp(-x))
        tq2 = (7.6 - Q) ** 2
        tt = (0.007 + 0.00005 * zm) * np.where(
            tq2 < 0.2, 1.0 - tq2 + 0.5 * tq2 * tq2, np.exp(-tq2)
        )
        out[low] = z * (1.0 + tt) * np.sqrt(ex)
        return out
    zi13 = z ** (1.0 / 3.0)
    zi23 = zi13 * zi13
    eF = 0.025  # MeV: 25 keV * vF^2 with vF = 1 (assumption, see docstring)
    v1sq = r / eF
    vF = math.sqrt(eF / 0.025)
    y = np.where(
        v1sq > 1.0,
        vF * np.sqrt(v1sq) * (1.0 + 0.2 / v1sq) / zi23,
        0.692308 * vF * (1.0 + 0.666666 * v1sq + v1sq * v1sq / 15.0) / zi23,
    )
    y3 = np.exp(0.3 * np.log(y))
    q = np.maximum(
        1.0 - np.exp(0.803 * y3 - 1.3167 * y3 * y3 - 0.38157 * y - 0.008983 * y * y),
        1.0 / z,
    )
    out[low] = z * q
    return out


@dataclass(frozen=True)
class ProjectileTables:
    """Log-uniform kinetic-energy grid of linear stopping powers for one ion+material."""

    particle: str
    material: Material
    lnE0: float
    dlnE: float
    L_el: np.ndarray  # MeV/cm electronic
    L_nuc: np.ndarray  # MeV/cm nuclear
    zeff2: np.ndarray  # effective charge squared (fluctuations, delta production)
    range_cm: np.ndarray  # CSDA range from L_el + L_nuc
    E_min: float
    E_max: float
    stopping_source: str


@functools.lru_cache(maxsize=None)
def projectile_tables(particle, material_name, n=1601):
    """Stopping tables for ``particle`` in ``material_name`` (NIST, SiC Bragg-composed;
    C-12 = z_eff^2 x PSTAR at the same velocity)."""
    if particle not in PARTICLES:
        raise ValueError(f"particle must be one of {sorted(PARTICLES)}")
    mat = get_material(material_name)
    mass, z, A, _ = PARTICLES[particle]
    if particle in ("proton", "alpha"):
        e_tab, _ = sp.mass_stopping(particle, "si", component="electronic")
        lo, hi = e_tab[0], e_tab[-1]
        if particle == "alpha":
            hi = min(hi, 1000.0)
        E = np.geomspace(lo, hi, n)
        L_el = sp.mass_stopping(particle, mat.star_name, E, "electronic") * mat.density
        L_nuc = sp.mass_stopping(particle, mat.star_name, E, "nuclear") * mat.density
        src = f"NIST {'PSTAR' if particle == 'proton' else 'ASTAR'} ({mat.star_name})"
        zeff = ion_effective_charge(particle, E, mat)
    else:
        # same velocity: proton energy Ep = E * m_p / M
        e_tab, _ = sp.mass_stopping("proton", "si", component="electronic")
        lo = e_tab[0] * mass / PROTON_MC2 * 25.0  # >= ~25 keV/u: Ziegler q range
        hi = min(e_tab[-1], 5000.0) * mass / PROTON_MC2
        E = np.geomspace(lo, hi, n)
        Ep = E * PROTON_MC2 / mass
        zeff = ion_effective_charge(particle, E, mat)
        L_el = (
            zeff**2
            * sp.mass_stopping("proton", mat.star_name, Ep, "electronic")
            * mat.density
        )
        L_nuc = np.zeros_like(E)
        src = "MODEL: z_eff^2 x NIST PSTAR at equal velocity (Ziegler z_eff)"
    L_tot = L_el + L_nuc
    p = -np.log(L_tot[1] / L_tot[0]) / np.log(E[1] / E[0])
    r0 = E[0] / (L_tot[0] * (1.0 + max(p, -0.9)))
    inv = 1.0 / L_tot
    rng = r0 + np.concatenate(
        [[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(E))]
    )
    for a in (L_el, L_nuc, zeff, rng):
        a.setflags(write=False)
    return ProjectileTables(
        particle=particle,
        material=mat,
        lnE0=float(np.log(E[0])),
        dlnE=float(np.log(E[1] / E[0])),
        L_el=L_el,
        L_nuc=L_nuc,
        zeff2=zeff**2,
        range_cm=rng,
        E_min=float(E[0]),
        E_max=float(E[-1]),
        stopping_source=src,
    )


@functools.lru_cache(maxsize=None)
def _electron_tables(estar_name, density):
    e, r = electron_csda_table(estar_name)
    lnT = np.log(e)
    lnR = np.log(r / density)
    return float(lnT[0]), float(lnT[1] - lnT[0]), lnR


def electron_range_um(T_keV, material="sic"):
    """Electron CSDA range (um) in ``material`` (NIST ESTAR, SiC composed)."""
    mat = get_material(material)
    lnT0, dlnT, lnR = _electron_tables(mat.estar_name, mat.density)
    T = np.asarray(T_keV, dtype=float) * 1e-3
    x = (np.log(T) - lnT0) / dlnT
    return np.exp(np.interp(x, np.arange(lnR.size), lnR)) * 1e4


# =====================================================================================
# numba kernel
# =====================================================================================

_GOLD = np.uint64(0x9E3779B97F4A7C15)
_M1 = np.uint64(0xBF58476D1CE4E5B9)
_M2 = np.uint64(0x94D049BB133111EB)


@njit(cache=True, inline="always")
def _next_u64(st):
    st[0] = st[0] + _GOLD
    z = st[0]
    z = (z ^ (z >> np.uint64(30))) * _M1
    z = (z ^ (z >> np.uint64(27))) * _M2
    return z ^ (z >> np.uint64(31))


@njit(cache=True, inline="always")
def _uniform(st):
    # (0, 1): 53 random bits, never exactly 0
    return ((_next_u64(st) >> np.uint64(11)) + 0.5) * (1.0 / 9007199254740992.0)


@njit(cache=True)
def _seed_state(seed, event):
    st = np.empty(1, dtype=np.uint64)
    st[0] = np.uint64(seed) * _M1 + np.uint64(event) * _GOLD + _M2
    for _ in range(4):
        _next_u64(st)
    return st


@njit(cache=True, inline="always")
def _gauss(st):
    u1 = _uniform(st)
    u2 = _uniform(st)
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


@njit(cache=True)
def _poisson(st, mean):
    # G4Poisson (Geant4 global/HEPRandom/include/G4Poisson.hh): direct sum for
    # mean <= 16, Gaussian approximation above.
    if mean <= 0.0:
        return 0
    if mean <= 16.0:
        pos = _uniform(st)
        pv = math.exp(-mean)
        ps = pv
        n = 0
        while ps <= pos and n < 1000:
            n += 1
            pv *= mean / n
            ps += pv
        return n
    v = mean + _gauss(st) * math.sqrt(mean) + 0.5
    if v <= 0.0:
        return 0
    if v > 2.0e9:
        return 2000000000
    return int(v)


@njit(cache=True)
def _gamma(st, a):
    # Marsaglia-Tsang; shape a < 1 by boosting
    if a < 1.0:
        return _gamma(st, a + 1.0) * _uniform(st) ** (1.0 / a)
    d = a - 1.0 / 3.0
    c = 1.0 / math.sqrt(9.0 * d)
    while True:
        x = _gauss(st)
        v = 1.0 + c * x
        if v <= 0.0:
            continue
        v = v * v * v
        u = _uniform(st)
        if math.log(u) < 0.5 * x * x + d - d * v + d * math.log(v):
            return d * v


@njit(cache=True)
def _sample_gauss_into(st, eav, esig2):
    # G4UniversalFluctuation::SampleGauss
    sig = math.sqrt(esig2)
    if eav < 0.25 * sig:
        return eav + (2.0 * _uniform(st) - 1.0) * eav
    for _ in range(1000):
        x = eav + sig * _gauss(st)
        if 0.0 <= x <= 2.0 * eav:
            return x
    return eav


# fluctuation model codes
FLUCT_NONE, FLUCT_UNIVERSAL, FLUCT_URBAN, FLUCT_ION, FLUCT_BOHR = 0, 1, 2, 3, 4
FLUCT_CODES = {
    "none": FLUCT_NONE,
    "universal": FLUCT_UNIVERSAL,
    "urban": FLUCT_URBAN,
    "ion": FLUCT_ION,
    "bohr": FLUCT_BOHR,
}


@njit(cache=True)
def _fluct_gauss_bohr(st, mean, siga2):
    siga = math.sqrt(siga2)
    sn = mean / siga
    if sn >= 2.0:
        for _ in range(1000):
            x = mean + siga * _gauss(st)
            if 0.0 <= x <= 2.0 * mean:
                return x
        return mean
    if sn > 0.1:
        neff = sn * sn
        return mean * _gamma(st, neff) / neff
    return 2.0 * mean * _uniform(st)


@njit(cache=True)
def _sample_fluct(
    st,
    model,
    mean,
    tcut,
    tmax,
    length,
    T,
    mass,
    q2,
    n_el,
    I,
    lnI,
    e0,
    f1,
    f2,
    e1f,
    e2f,
    lnE1,
    lnE2,
    ion_thr,
):
    """One step's sampled (restricted) energy loss, MeV. Port of Geant4 11.3."""
    if model == FLUCT_NONE or mean < 1.0e-5:  # minLoss = 10 eV
        return mean
    gam = T / mass + 1.0
    gam2 = gam * gam
    beta2 = 1.0 - 1.0 / gam2
    siga2 = (tmax / beta2 - 0.5 * tcut) * TWOPI_MC2_RCL2 * length * q2 * n_el
    if model == FLUCT_BOHR:
        return _fluct_gauss_bohr(st, mean, siga2)
    if model == FLUCT_ION and T <= ion_thr:
        # G4IonFluctuations low-velocity branch WITHOUT Yang's factor (documented),
        # with its slowing-down correction of the dispersion
        if mean > 0.2 * T:
            g2 = (T - mean) / mass + 1.0
            b2 = 1.0 - 1.0 / (g2 * g2)
            if b2 < 0.2 * beta2:
                b2 = 0.2 * beta2
            x = b2 / beta2
            x3 = 1.0 / (x * x * x)
            siga2 *= 0.25 * (1.0 + x) * (x3 + (1.0 / b2 - 0.5) / (1.0 / beta2 - 0.5))
        return _fluct_gauss_bohr(st, mean, siga2)
    # G4UniversalFluctuation::SampleFluctuations
    if mean >= 10.0 * tcut and tmax <= 2.0 * tcut:
        siga = math.sqrt(siga2)
        sn = mean / siga
        if sn >= 2.0:
            for _ in range(1000):
                x = mean + siga * _gauss(st)
                if 0.0 <= x <= 2.0 * mean:
                    return x
            return mean
        neff = sn * sn
        return mean * _gamma(st, neff) / neff
    if tcut <= e0:
        return mean
    scaling = min(1.0 + 0.5e-3 / tcut, 1.5)
    meanl = mean / scaling
    w2 = math.log(2.0 * ME_C2 * beta2 * gam2) - beta2 if tcut > I else 0.0
    rate, fw, a0, nmax = 0.56, 4.0, 42.0, 8.0
    a1 = 0.0
    a2 = 0.0
    loss = 0.0
    if model == FLUCT_URBAN:
        # G4UrbanFluctuation::SampleGlandz (two excitation levels)
        e1 = e1f
        e2 = e2f
        if tcut > I and w2 > lnI:
            if w2 > lnE2:
                C = meanl * (1.0 - rate) / (w2 - lnI)
                a1 = C * f1 * (w2 - lnE1) / e1f
                a2 = C * f2 * (w2 - lnE2) / e2f
            else:
                a1 = meanl * (1.0 - rate) / e1
            if a1 < a0:
                fwnow = 0.5 + (fw - 0.5) * math.sqrt(a1 / a0)
                a1 /= fwnow
                e1 *= fwnow
            else:
                a1 /= fw
                e1 *= fw
    else:
        # G4UniversalFluctuation::SampleGlandz (one excitation level, E1 = I)
        e1 = I
        e2 = 0.0
        if tcut > e1:
            a1 = meanl * (1.0 - rate) / e1
            if a1 < a0:
                fwnow = 0.1 + (fw - 0.1) * math.sqrt(a1 / a0)
                a1 /= fwnow
                e1 *= fwnow
            else:
                a1 /= fw
                e1 *= fw
    w1 = tcut / e0
    a3 = rate * meanl * (tcut - e0) / (e0 * tcut * math.log(w1))
    if a1 + a2 <= 0.0:
        a3 /= rate
    emean = 0.0
    sig2e = 0.0
    for lev in range(2):
        ax = a1 if lev == 0 else a2
        ex = e1 if lev == 0 else e2
        if ax > 0.0:
            if ax > nmax:
                emean += ax * ex
                sig2e += ax * ex * ex
            else:
                p = _poisson(st, ax)
                if p > 0:
                    loss += ((p + 1) - 2.0 * _uniform(st)) * ex
    if sig2e > 0.0:
        loss += _sample_gauss_into(st, emean, sig2e)
    if a3 > 0.0:
        emean = 0.0
        sig2e = 0.0
        p3 = a3
        alfa = 1.0
        if a3 > nmax:
            alfa = w1 * (nmax + a3) / (w1 * nmax + a3)
            alfa1 = alfa * math.log(alfa) / (alfa - 1.0)
            namean = a3 * w1 * (alfa - 1.0) / ((w1 - 1.0) * alfa)
            emean += namean * e0 * alfa1
            sig2e += e0 * e0 * namean * (alfa - alfa1 * alfa1)
            p3 = a3 - namean
        w3 = alfa * e0
        if tcut > w3:
            w = (tcut - w3) / tcut
            nnb = _poisson(st, p3)
            for _ in range(nnb):
                loss += w3 / (1.0 - w * _uniform(st))
        if sig2e > 0.0:
            loss += _sample_gauss_into(st, emean, sig2e)
    return loss * scaling


@njit(cache=True, inline="always")
def _interp_lin(tab, lnE0, dlnE, E):
    x = (math.log(E) - lnE0) / dlnE
    if x <= 0.0:
        return tab[0]
    n = tab.size
    if x >= n - 1:
        return tab[n - 1]
    i = int(x)
    f = x - i
    return tab[i] * (1.0 - f) + tab[i + 1] * f


@njit(cache=True, inline="always")
def _interp_loglog(tab, lnE0, dlnE, E):
    x = (math.log(E) - lnE0) / dlnE
    n = tab.size
    if x <= 0.0:
        return tab[0]
    if x >= n - 1:
        return tab[n - 1]
    i = int(x)
    f = x - i
    return math.exp(math.log(tab[i]) * (1.0 - f) + math.log(tab[i + 1]) * f)


@njit(cache=True)
def _energy_from_range(rng_tab, lnE0, dlnE, R):
    # inverse of the (increasing) projectile range table
    n = rng_tab.size
    if R <= rng_tab[0]:
        return math.exp(lnE0) * R / rng_tab[0]
    if R >= rng_tab[n - 1]:
        return math.exp(lnE0 + dlnE * (n - 1))
    i = np.searchsorted(rng_tab, R) - 1
    f = (R - rng_tab[i]) / (rng_tab[i + 1] - rng_tab[i])
    return math.exp(lnE0 + dlnE * (i + f))


@njit(cache=True)
def _e_range(lnR, lnT0, dlnT, T):
    # electron range (cm) at T (MeV), log-log
    if T <= 0.0:
        return 0.0
    x = (math.log(T) - lnT0) / dlnT
    n = lnR.size
    if x <= 0.0:
        # below the table (1 keV): R ~ T^(1+p) power law continuation
        p1 = (lnR[1] - lnR[0]) / dlnT
        return math.exp(lnR[0] + p1 * (math.log(T) - lnT0))
    if x >= n - 1:
        return math.exp(lnR[n - 1])
    i = int(x)
    f = x - i
    return math.exp(lnR[i] * (1.0 - f) + lnR[i + 1] * f)


@njit(cache=True)
def _e_energy(lnR, lnT0, dlnT, R):
    # inverse: electron energy (MeV) with residual range R (cm)
    if R <= 0.0:
        return 0.0
    lr = math.log(R)
    n = lnR.size
    if lr <= lnR[0]:
        p1 = (lnR[1] - lnR[0]) / dlnT
        return math.exp(lnT0 + (lr - lnR[0]) / p1)
    if lr >= lnR[n - 1]:
        return math.exp(lnT0 + dlnT * (n - 1))
    i = np.searchsorted(lnR, lr) - 1
    f = (lr - lnR[i]) / (lnR[i + 1] - lnR[i])
    return math.exp(lnT0 + dlnT * (i + f))


@njit(cache=True, inline="always")
def _ray_box(px, py, pz, dx, dy, dz, lo, hi):
    """Slab method: (t_in, t_out) of the ray p + t d with the box; t_in > t_out if miss."""
    t0 = -1e300
    t1 = 1e300
    for k in range(3):
        p = px if k == 0 else (py if k == 1 else pz)
        d = dx if k == 0 else (dy if k == 1 else dz)
        if abs(d) < 1e-300:
            if p < lo[k] or p > hi[k]:
                return 1.0, -1.0
        else:
            a = (lo[k] - p) / d
            b = (hi[k] - p) / d
            if a > b:
                a, b = b, a
            if a > t0:
                t0 = a
            if b < t1:
                t1 = b
    return t0, t1


@njit(cache=True)
def _sample_entry(st, mode, bdir, plo, phi):
    """Entry point on the padded box surface and inward unit direction.

    mode 0: parallel beam along ``bdir``, uniform over the projected area;
    mode 1: isotropic (mu-randomness): face ~ area, cosine-law inward direction.
    """
    a = np.empty(6)
    L = phi - plo
    areas = (
        L[1] * L[2],
        L[1] * L[2],
        L[0] * L[2],
        L[0] * L[2],
        L[0] * L[1],
        L[0] * L[1],
    )
    # faces: 0 x-lo (n=-x), 1 x-hi, 2 y-lo, 3 y-hi, 4 z-lo, 5 z-hi
    for f in range(6):
        ax = f // 2
        out_sign = -1.0 if f % 2 == 0 else 1.0
        if mode == 0:
            c = -out_sign * bdir[ax]  # inward cos with beam
            a[f] = areas[f] * c if c > 0 else 0.0
        else:
            a[f] = areas[f]
    tot = a.sum()
    u = _uniform(st) * tot
    face = 5
    acc = 0.0
    for f in range(6):
        acc += a[f]
        if u <= acc:
            face = f
            break
    ax = face // 2
    out_sign = -1.0 if face % 2 == 0 else 1.0
    p = np.empty(3)
    for k in range(3):
        p[k] = plo[k] + _uniform(st) * L[k]
    p[ax] = plo[ax] if face % 2 == 0 else phi[ax]
    d = np.empty(3)
    if mode == 0:
        d[0] = bdir[0]
        d[1] = bdir[1]
        d[2] = bdir[2]
    else:
        ct = math.sqrt(_uniform(st))
        stt = math.sqrt(max(0.0, 1.0 - ct * ct))
        ph = 2.0 * math.pi * _uniform(st)
        # local frame: n_in = -out_sign e_ax
        e1 = (ax + 1) % 3
        e2 = (ax + 2) % 3
        d[ax] = -out_sign * ct
        d[e1] = stt * math.cos(ph)
        d[e2] = stt * math.sin(ph)
    return p, d


@njit(cache=True, inline="always")
def _rotate_dir(dx, dy, dz, ct, ph):
    """Direction at polar angle acos(ct), azimuth ph around unit vector d."""
    st_ = math.sqrt(max(0.0, 1.0 - ct * ct))
    # build orthonormal basis (u, v, d)
    if abs(dz) < 0.9:
        ux, uy, uz = -dy, dx, 0.0
    else:
        ux, uy, uz = 0.0, -dz, dy
    nu = math.sqrt(ux * ux + uy * uy + uz * uz)
    ux /= nu
    uy /= nu
    uz /= nu
    vx = dy * uz - dz * uy
    vy = dz * ux - dx * uz
    vz = dx * uy - dy * ux
    c = math.cos(ph) * st_
    s = math.sin(ph) * st_
    return (
        c * ux + s * vx + ct * dx,
        c * uy + s * vy + ct * dy,
        c * uz + s * vz + ct * dz,
    )


@njit(cache=True)
def _add_track_hits(hits, n_hits, i, nh, max_hits, p, d, s0, length, energy, hit_step):
    """Record a primary deposit spread evenly over ceil(length / hit_step) points of the
    path [s0, s0 + length] (recording only: the sampled step is unchanged)."""
    nseg = 1
    if length > 0.0 and hit_step > 0.0:
        nseg = max(1, int(math.ceil(length / hit_step - 1e-9)))
    e = energy / nseg
    for k in range(nseg):
        sk = s0 + (k + 0.5) * length / nseg
        if nh < max_hits:
            hits[i, nh, 0] = p[0] + sk * d[0]
            hits[i, nh, 1] = p[1] + sk * d[1]
            hits[i, nh, 2] = p[2] + sk * d[2]
            hits[i, nh, 3] = e
            hits[i, nh, 4] = 1.0
            nh += 1
        else:
            hits[i, max_hits - 1, 3] += e
            n_hits[i, 1] += 1
    return nh


@njit(cache=True, parallel=True)
def _transport_kernel(
    n,
    event_offset,
    seed,
    E0,
    mode,
    bdir,
    svlo,
    svhi,
    margin,
    lnE0,
    dlnE,
    L_el,
    L_nuc,
    zeff2,
    rng_tab,
    E_min,
    mass,
    n_el,
    I,
    lnI,
    e0,
    f1,
    f2,
    e1f,
    e2f,
    lnE1,
    lnE2,
    ion_thr,
    tcut,
    fluct,
    delta_on,
    lnR,
    lnT0,
    dlnT,
    rfac,
    step_frac,
    max_hits,
    hit_step,
    hits,
    n_hits,
    edep,
    edep_p,
    edep_d,
    chord,
    e_in,
    entry,
    direc,
    stop_s,
    ndelta,
):
    plo = svlo - margin
    phi = svhi + margin
    for i in prange(n):
        ev = event_offset + i
        st = _seed_state(seed, ev)
        E = E0[ev % E0.size] if E0.size > 1 else E0[0]
        p, d = _sample_entry(st, mode, bdir, plo, phi)
        t0, t1 = _ray_box(p[0], p[1], p[2], d[0], d[1], d[2], plo, phi)
        s_end = t1
        a, b = _ray_box(p[0], p[1], p[2], d[0], d[1], d[2], svlo, svhi)
        hit_sv = b > a and b > 0.0
        s_in = max(a, 0.0) if hit_sv else 1e300
        s_out = b if hit_sv else -1e300
        chord[i] = (s_out - s_in) if hit_sv else 0.0
        for k in range(3):
            direc[i, k] = d[k]
            entry[i, k] = p[k] + s_in * d[k] if hit_sv else np.nan
        e_in[i] = np.nan
        dep_p = 0.0
        dep_d = 0.0
        nd = 0
        nh = 0
        s = 0.0
        stop_s[i] = np.nan
        alive = True
        guard = 0
        while alive and s < s_end - 1e-15 and guard < 100000:
            guard += 1
            if hit_sv and s >= s_in - 1e-15 and np.isnan(e_in[i]):
                e_in[i] = E
            # next geometric boundary
            nb = s_end
            if hit_sv:
                if s < s_in - 1e-15:
                    nb = s_in
                elif s < s_out - 1e-15:
                    nb = s_out
            inside = hit_sv and (s >= s_in - 1e-15) and (s < s_out - 1e-15)
            if E <= E_min:
                if inside:
                    dep_p += E
                    if max_hits > 0:
                        nh = _add_track_hits(
                            hits, n_hits, i, nh, max_hits, p, d, s, 0.0, E, hit_step
                        )
                stop_s[i] = s
                E = 0.0
                alive = False
                break
            Lt = _interp_loglog(L_el, lnE0, dlnE, E) + _interp_lin(L_nuc, lnE0, dlnE, E)
            ds = min(nb - s, step_frac * E / Lt)
            R = _interp_lin(rng_tab, lnE0, dlnE, E)
            if R <= ds * (1.0 + 1e-12):
                # stops within this step: deposit the rest (electronic), end of range
                if inside:
                    dep_p += E
                    if max_hits > 0:
                        nh = _add_track_hits(
                            hits, n_hits, i, nh, max_hits, p, d, s, R, E, hit_step
                        )
                stop_s[i] = s + R
                E = 0.0
                alive = False
                break
            # mid-step energy for the mean loss (2nd order)
            Em = _energy_from_range(rng_tab, lnE0, dlnE, R - 0.5 * ds)
            Lel = _interp_loglog(L_el, lnE0, dlnE, Em)
            Lnu = _interp_lin(L_nuc, lnE0, dlnE, Em)
            q2 = _interp_lin(zeff2, lnE0, dlnE, Em)
            Eafter = _energy_from_range(rng_tab, lnE0, dlnE, R - ds)
            mean_tot = E - Eafter  # CSDA loss over the step (el + nuc)
            frac_el = Lel / (Lel + Lnu)
            mean_el = mean_tot * frac_el
            gam = Em / mass + 1.0
            beta2 = 1.0 - 1.0 / (gam * gam)
            mr = ME_C2 / mass
            tmax = 2.0 * ME_C2 * beta2 * gam * gam / (1.0 + 2.0 * gam * mr + mr * mr)
            tc = tcut if (delta_on and tcut < tmax) else tmax
            mean_res = mean_el
            nd_step = 0
            if delta_on and tmax > tcut:
                pref = TWOPI_MC2_RCL2 * n_el * q2 / beta2
                sub = pref * (math.log(tmax / tcut) - beta2 * (1.0 - tcut / tmax)) * ds
                if sub < mean_el:
                    mean_res = mean_el - sub
                else:
                    mean_res = 0.0
                sig = pref * (
                    (1.0 / tcut - 1.0 / tmax) - beta2 / tmax * math.log(tmax / tcut)
                )
                nd_step = _poisson(st, sig * ds)
            loss = _sample_fluct(
                st,
                fluct,
                mean_res,
                tc,
                tmax,
                ds,
                Em,
                mass,
                q2,
                n_el,
                I,
                lnI,
                e0,
                f1,
                f2,
                e1f,
                e2f,
                lnE1,
                lnE2,
                ion_thr,
            )
            dE = loss + (mean_tot - mean_el)
            if inside:
                dep_p += loss
                if max_hits > 0:
                    nh = _add_track_hits(
                        hits, n_hits, i, nh, max_hits, p, d, s, ds, loss, hit_step
                    )
            for _j in range(nd_step):
                # T from 1/T^2 between tcut and tmax, accept (1 - beta2 T / tmax)
                T = tcut
                for _r in range(1000):
                    u = _uniform(st)
                    T = tcut * tmax / (tcut + u * (tmax - tcut))
                    if _uniform(st) <= 1.0 - beta2 * T / tmax:
                        break
                if T >= E - dE:
                    continue
                dE += T
                nd += 1
                sd = s + _uniform(st) * ds
                qx = p[0] + sd * d[0]
                qy = p[1] + sd * d[1]
                qz = p[2] + sd * d[2]
                # two-body kinematics (M >> m): cos^2 = T (tmax + 2m) / (tmax (T + 2m))
                ct = math.sqrt(
                    min(1.0, T * (tmax + 2.0 * ME_C2) / (tmax * (T + 2.0 * ME_C2)))
                )
                ux, uy, uz = _rotate_dir(
                    d[0], d[1], d[2], ct, 2.0 * math.pi * _uniform(st)
                )
                Rd = rfac * _e_range(lnR, lnT0, dlnT, T)
                ta, tb = _ray_box(qx, qy, qz, ux, uy, uz, svlo, svhi)
                if tb > ta and tb > 0.0:
                    ta = max(ta, 0.0)
                    if ta < Rd:
                        tb = min(tb, Rd)
                        e_a = (
                            _e_energy(lnR, lnT0, dlnT, (Rd - ta) / rfac)
                            if ta > 0.0
                            else T
                        )
                        e_b = (
                            _e_energy(lnR, lnT0, dlnT, (Rd - tb) / rfac)
                            if tb < Rd
                            else 0.0
                        )
                        dd = e_a - e_b
                        if dd > 0.0:
                            dep_d += dd
                            if max_hits > 0:
                                tm = 0.5 * (ta + tb)
                                if nh < max_hits:
                                    hits[i, nh, 0] = qx + tm * ux
                                    hits[i, nh, 1] = qy + tm * uy
                                    hits[i, nh, 2] = qz + tm * uz
                                    hits[i, nh, 3] = dd
                                    hits[i, nh, 4] = 1.0 + nd
                                    nh += 1
                                else:
                                    hits[i, max_hits - 1, 3] += dd
                                    n_hits[i, 1] += 1
            E -= dE
            s += ds
            if E <= 0.0:
                stop_s[i] = s
                alive = False
        edep_p[i] = dep_p * 1e3
        edep_d[i] = dep_d * 1e3
        edep[i] = (dep_p + dep_d) * 1e3
        ndelta[i] = nd
        n_hits[i, 0] = nh
        if max_hits > 0:
            for k in range(nh):
                hits[i, k, 3] *= 1e3


# =====================================================================================
# Python API
# =====================================================================================


@dataclass
class TrackEvents:
    """Per-event results of :func:`sample_track_events` (arrays of length n_events).

    edep_keV = edep_primary_keV (restricted, fluctuated, local) + edep_delta_keV (delta
    rays above T_cut that end inside the SV). chord_um is the primary's geometric chord
    in the SV (0: the track missed the SV; such events can still score delta energy
    with a surround). Positions in um, SV frame: x, y in [-w/2, w/2], z in [0, t]
    (z = 0 is the beam-entrance face for the normal beam).
    """

    edep_keV: np.ndarray
    edep_primary_keV: np.ndarray
    edep_delta_keV: np.ndarray
    chord_um: np.ndarray
    E_sv_entry_MeV: np.ndarray
    entry_um: np.ndarray
    direction: np.ndarray
    stop_path_um: np.ndarray
    n_delta: np.ndarray
    hits: np.ndarray | None
    n_hits: np.ndarray | None
    meta: dict = field(default_factory=dict)

    @property
    def n_events(self):
        return int(self.edep_keV.size)

    def lineal_energy(self, l_bar_um):
        """y = epsilon / l_bar (keV/um) of the events with epsilon > 0."""
        e = self.edep_keV
        return e[e > 0] / l_bar_um


def auto_tcut_keV(sv_um, material="sic", fraction=0.05):
    """Delta-ray production cut: electron CSDA range = ``fraction`` x the smallest SV
    dimension, clamped to [1, 1000] keV (the ESTAR table starts at 1 keV)."""
    target = fraction * float(min(sv_um))
    T = np.geomspace(1.0, 1000.0, 2001)
    r = electron_range_um(T, material)
    return float(np.clip(np.interp(target, r, T), 1.0, 1000.0))


def _beam_direction(irradiation):
    if irradiation == "normal":
        return 0, np.array([0.0, 0.0, 1.0])
    if irradiation == "isotropic":
        return 1, np.array([0.0, 0.0, 1.0])
    th = math.radians(float(irradiation))
    if not (0.0 <= th < math.pi / 2):
        raise ValueError("irradiation angle must be in [0, 90) degrees")
    return 0, np.array([math.sin(th), 0.0, math.cos(th)])


def default_fluctuation_model(particle):
    """Geant4 11.3 G4EmStandardPhysics_option4 choice: protons G4UrbanFluctuation,
    alpha / ions G4IonFluctuations (-> G4UniversalFluctuation at high velocity)."""
    return "urban" if particle == "proton" else "ion"


def sample_track_events(
    particle,
    energy_MeV,
    sv_um=(100.0, 100.0, 10.0),
    material="sic",
    irradiation="normal",
    n_events=10000,
    seed=0,
    tcut_keV=None,
    delta_rays=True,
    fluctuations=None,
    surround_um=0.0,
    delta_range_factor=1.0,
    step_fraction=0.05,
    max_hits=0,
    event_offset=0,
    hit_step_um=1.0,
):
    """Sample ``n_events`` primary tracks through a box SV and score the deposits.

    Parameters
    ----------
    particle : {"proton", "alpha", "c12"}
    energy_MeV : float or array
        Kinetic energy (total, MeV; C-12: 12 x MeV/u) at the entrance of the scoring
        box (the SV plus ``surround_um``). An array is a beam spectrum: event k uses
        ``energy_MeV[k % len]`` (shuffle it first if it is ordered).
    sv_um : (wx, wy, t)
        SV size (um); the SV is [-wx/2, wx/2] x [-wy/2, wy/2] x [0, t].
    irradiation : "normal" | "isotropic" | angle in degrees
        Parallel beam along +z (uniform over the face), mu-randomness isotropic chords,
        or a parallel beam tilted by the angle in the x-z plane.
    tcut_keV : float or None
        Delta-ray production cut (None: :func:`auto_tcut_keV`, range = 5 % of the
        smallest SV dimension). Below it the loss is local (restricted + fluctuations).
    delta_rays : bool
        False: no explicit delta rays (all loss local, T_cut -> T_max, i.e. the full
        Landau/Vavilov-like loss deposited at the track).
    fluctuations : {"urban", "universal", "ion", "bohr", "none"} or None
        None = Geant4 option4 choice (:func:`default_fluctuation_model`).
        "none" + delta_rays=False is the deterministic CSDA baseline (the deposit is the
        mean loss along the chord, as ETNA's LET x chord estimate).
    surround_um : float
        Thickness of the same material around the SV (delta rays produced there can
        enter the SV; tracks start at its outer face, so it also degrades the energy).
        0 = bare SV in vacuum.
    delta_range_factor : float
        Straight-line delta path = factor x CSDA range (1 = maximal escape distance).
    step_fraction : float
        Maximum fractional energy loss per step (adaptive steps; boundaries split steps).
    max_hits : int
        > 0 records up to this many energy deposits per event (positions + energy) for
        the hit-list exchange (:mod:`etna.core.mc_exchange`).
    hit_step_um : float
        Recorded primary deposits are spread over points <= this far apart along the
        path (default 1 um, Allpix DepositionGeant4's ``max_step_length``); recording
        only, the sampled steps (and every distribution) are unchanged.
    """
    if fluctuations is None:
        fluctuations = default_fluctuation_model(particle)
    if fluctuations not in FLUCT_CODES:
        raise ValueError(f"fluctuations must be one of {sorted(FLUCT_CODES)}")
    sv = np.asarray(sv_um, dtype=float)
    if sv.shape != (3,) or np.any(~np.isfinite(sv)) or np.any(sv <= 0):
        raise ValueError("sv_um must be three positive finite sizes (um)")
    if not (0.0 < step_fraction <= 0.5):
        raise ValueError("step_fraction must be in (0, 0.5]")
    if not (hit_step_um > 0):
        raise ValueError("hit_step_um must be > 0")
    if surround_um < 0 or delta_range_factor <= 0:
        raise ValueError("surround_um must be >= 0 and delta_range_factor > 0")
    n_events = int(n_events)
    if n_events < 1:
        raise ValueError("n_events must be >= 1")
    tab = projectile_tables(particle, material)
    mat = tab.material
    E0 = np.atleast_1d(np.asarray(energy_MeV, dtype=float))
    if np.any(~np.isfinite(E0)) or np.any(E0 <= tab.E_min) or np.any(E0 > tab.E_max):
        raise ValueError(
            f"{particle} energy must be in ({tab.E_min:g}, {tab.E_max:g}] MeV for the "
            "stopping table"
        )
    if tcut_keV is None:
        tcut_keV = auto_tcut_keV(sv, material)
    if not (1.0 <= tcut_keV <= 1000.0):
        raise ValueError("tcut_keV must be in [1, 1000] keV (ESTAR table range)")
    mode, bdir = _beam_direction(irradiation)
    svlo = np.array([-sv[0] / 2, -sv[1] / 2, 0.0]) * 1e-4
    svhi = np.array([sv[0] / 2, sv[1] / 2, sv[2]]) * 1e-4
    lnT0, dlnT, lnR = _electron_tables(mat.estar_name, mat.density)
    # Urban two-level parameters (G4IonisParamMat::ComputeFluctModel)
    zeff = mat.zeff_mass
    I = mat.I_eV * 1e-6
    lnI = math.log(I)
    if zeff > 2.1:
        f2 = 2.0 / zeff
        f1 = 1.0 - f2
        e2 = 10.0 * zeff * zeff * 1e-6
        lnE2 = math.log(e2)
        lnE1 = (lnI - f2 * lnE2) / f1
    else:
        f2, f1 = 0.0, 1.0
        e2 = 40e-6 if zeff > 1.1 else 10e-6
        lnE2 = math.log(e2)
        lnE1 = lnI
    e1 = math.exp(lnE1)
    mass = PARTICLES[particle][0]
    nh = max(int(max_hits), 0)
    hits = np.zeros((n_events if nh else 1, max(nh, 1), 5))
    n_hits = np.zeros((n_events, 2), dtype=np.int64)
    out = {
        k: np.empty(n_events)
        for k in ("edep", "edep_p", "edep_d", "chord", "e_in", "stop_s")
    }
    entry = np.empty((n_events, 3))
    direc = np.empty((n_events, 3))
    ndelta = np.empty(n_events, dtype=np.int64)
    _transport_kernel(
        n_events,
        int(event_offset),
        int(seed) & 0xFFFFFFFF,
        E0,
        mode,
        bdir,
        svlo,
        svhi,
        surround_um * 1e-4,
        tab.lnE0,
        tab.dlnE,
        np.asarray(tab.L_el),
        np.asarray(tab.L_nuc),
        np.asarray(tab.zeff2),
        np.asarray(tab.range_cm),
        tab.E_min,
        mass,
        mat.electron_density_cm3,
        I,
        lnI,
        10e-6,
        f1,
        f2,
        e1,
        e2,
        lnE1,
        lnE2,
        10.0 / PROTON_MC2 * PARTICLES[particle][1] * mass,  # G4IonFluctuations limit
        tcut_keV * 1e-3,
        FLUCT_CODES[fluctuations],
        bool(delta_rays),
        lnR,
        lnT0,
        dlnT,
        float(delta_range_factor),
        float(step_fraction),
        nh,
        float(hit_step_um) * 1e-4,
        hits,
        n_hits,
        out["edep"],
        out["edep_p"],
        out["edep_d"],
        out["chord"],
        out["e_in"],
        entry,
        direc,
        out["stop_s"],
        ndelta,
    )
    meta = {
        "particle": particle,
        "material": material,
        "sv_um": tuple(float(v) for v in sv),
        "irradiation": irradiation,
        "tcut_keV": float(tcut_keV),
        "delta_rays": bool(delta_rays),
        "fluctuations": fluctuations,
        "surround_um": float(surround_um),
        "delta_range_factor": float(delta_range_factor),
        "step_fraction": float(step_fraction),
        "seed": int(seed),
        "event_offset": int(event_offset),
        "stopping_source": tab.stopping_source,
        "label": "PREDICTED (NIST stopping + Geant4-port fluctuations; not fitted)",
    }
    if nh:
        hits[:, :, 0:3] *= 1e4  # cm -> um
    return TrackEvents(
        edep_keV=out["edep"],
        edep_primary_keV=out["edep_p"],
        edep_delta_keV=out["edep_d"],
        chord_um=out["chord"] * 1e4,
        E_sv_entry_MeV=out["e_in"],
        entry_um=entry * 1e4,
        direction=direc,
        stop_path_um=out["stop_s"] * 1e4,
        n_delta=ndelta,
        hits=hits if nh else None,
        n_hits=n_hits if nh else None,
        meta=meta,
    )


# --- analytic references (PDG "Passage of particles through matter") ---------------

K_PDG = 0.307075  # MeV mol^-1 cm^2


def kinematics(particle, T_MeV):
    """(beta^2, gamma, T_max MeV) of the ion at kinetic energy T."""
    mass = PARTICLES[particle][0]
    gam = T_MeV / mass + 1.0
    b2 = 1.0 - 1.0 / gam**2
    mr = ME_C2 / mass
    tmax = 2 * ME_C2 * b2 * gam**2 / (1 + 2 * gam * mr + mr**2)
    return b2, gam, tmax


def landau_xi_keV(particle, T_MeV, thickness_um, material="sic", z=None):
    """xi = (K/2) (Z/A) z^2 (rho x / beta^2) in keV."""
    mat = get_material(material)
    b2, _, _ = kinematics(particle, T_MeV)
    zz = PARTICLES[particle][1] if z is None else z
    return (
        0.5
        * K_PDG
        * mat.z_over_a
        * zz**2
        * mat.density
        * thickness_um
        * 1e-4
        / b2
        * 1e3
    )


def vavilov_kappa(particle, T_MeV, thickness_um, material="sic"):
    """kappa = xi / T_max (Landau << 0.01 < Vavilov < ~10 < Gaussian)."""
    _, _, tmax = kinematics(particle, T_MeV)
    return landau_xi_keV(particle, T_MeV, thickness_um, material) / (tmax * 1e3)


def landau_mpv_keV(particle, T_MeV, thickness_um, material="sic"):
    """PDG most probable loss Delta_p = xi [ln(2mc^2 b^2 g^2 / I) + ln(xi/I) + j - b^2 - delta],
    j = 0.200, delta (density effect) neglected (< 1e-2 relative below ~1 GeV/u)."""
    mat = get_material(material)
    b2, gam, _ = kinematics(particle, T_MeV)
    xi = landau_xi_keV(particle, T_MeV, thickness_um, material) * 1e-3
    I = mat.I_eV * 1e-6
    v = xi * (math.log(2 * ME_C2 * b2 * gam**2 / I) + math.log(xi / I) + 0.200 - b2)
    return v * 1e3


def landau_fwhm_keV(particle, T_MeV, thickness_um, material="sic"):
    """Landau FWHM = 4.018 xi (PDG)."""
    return 4.018 * landau_xi_keV(particle, T_MeV, thickness_um, material)


def bohr_variance_keV2(particle, T_MeV, thickness_um, material="sic"):
    """Vavilov/Bohr variance of the total loss: xi T_max (1 - beta^2/2) (keV^2)."""
    b2, _, tmax = kinematics(particle, T_MeV)
    return (
        landau_xi_keV(particle, T_MeV, thickness_um, material)
        * tmax
        * 1e3
        * (1 - 0.5 * b2)
    )


def mean_csda_deposit_keV(particle, E_MeV, path_um, material="sic"):
    """Mean electronic CSDA energy loss (keV) over a straight path of ``path_um`` from
    entry energy ``E_MeV``: the integral of L_el / L_tot dE between the entry and exit
    energies of the range table (the deterministic baseline; a particle that stops
    deposits its whole residual energy below the table minimum, as the sampler)."""
    tab = projectile_tables(particle, material)
    E = np.asarray(E_MeV, dtype=float)
    lnE = np.log(tab.E_min) + tab.dlnE * np.arange(tab.L_el.size)
    Egrid = np.exp(lnE)
    rng = np.asarray(tab.range_cm)
    R1 = np.interp(np.log(E), lnE, rng) - np.asarray(path_um, dtype=float) * 1e-4
    stops = R1 <= rng[0]
    E1 = np.exp(np.interp(np.maximum(R1, rng[0]), rng, lnE))
    frac = tab.L_el / (tab.L_el + tab.L_nuc)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (frac[1:] + frac[:-1]) * np.diff(Egrid))])
    hi = np.interp(E, Egrid, cum)
    lo = np.where(stops, -tab.E_min, np.interp(E1, Egrid, cum))
    return (hi - lo) * 1e3
