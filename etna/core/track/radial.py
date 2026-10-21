"""Radial distribution of the pairs around an ion path (track structure).

Every model is a cumulative fraction F(r) of the deposited energy (= of the pairs,
at a constant pair energy) within radius r of the path. The engine lumps the pairs
of an annulus [r_a, r_b] into that finite-volume cell exactly, F(r_b) - F(r_a), so
the total is conserved on any mesh and a core thinner than the first cell is a
sub-grid quantity (the cell carries its mean density). The share beyond the domain
radius R goes into the outermost cell.

Models
------
``gaussian(w)``
    rho(r) ~ exp(-r^2 / w^2), F = 1 - exp(-r^2 / w^2): the generic commercial-TCAD
    heavy-ion form (one radius parameter; the exact form and default of Sentaurus'
    HeavyIon ``Wt_hi`` are NOT verified here, see the issue).
``core_penumbra(r_min, r_max)``
    D = const for r < r_min, ~ 1/r^2 up to r_max, 0 beyond (the Local-Effect-Model
    form, Korcyl arXiv:1410.5250 eq. 1.20, lambda normalising to LET).
``whk(beta)``
    Waligorski, Hamm & Katz, Nucl. Tracks Radiat. Meas. 11, 309 (1986): the Zhang
    et al. (1985) delta-ray dose with the range-energy law R = k omega^alpha and the
    WHK correction (1 + K(r)), constants as given for WATER in Korcyl arXiv:1410.5250
    eqs. 2.3-2.12 (checked on the rendered page). Transferred to SiC by mass-thickness
    scaling (every length x rho_water / rho_SiC), the usual Katz practice; this
    transfer is OUR approximation, not a published SiC fit. Only the shape is used: it
    is normalised to the stopping-power LET (the WHK integral reproduces it in water).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

M_E_C2_MEV = 0.51099895
RHO_WATER = 1.0  # g/cm^3
RHO_SIC = 3.21


def beta_of(E_MeV_per_u):
    """Ion speed beta from the kinetic energy per nucleon (MeV/u, u = 931.494 MeV)."""
    g = 1.0 + float(E_MeV_per_u) / 931.49410242
    return math.sqrt(1.0 - 1.0 / (g * g))


@dataclass(frozen=True)
class RadialModel:
    name: str
    params: tuple

    def cdf(self, r_cm):
        r = np.asarray(r_cm, dtype=float)
        if self.name == "gaussian":
            (w,) = self.params
            return -np.expm1(-((r / w) ** 2))
        if self.name == "core_penumbra":
            a, b = self.params
            tot = 1.0 + 2.0 * math.log(b / a)
            inner = (r / a) ** 2
            outer = 1.0 + 2.0 * np.log(np.clip(r, a, b) / a)
            return np.where(r < a, inner, outer) / tot
        if self.name == "whk":
            rg, F = _whk_table(*self.params)
            return np.interp(r, rg, F, left=0.0, right=1.0)
        if self.name == "thermalized":
            rg, F = _conv_table(*self.params)
            return np.interp(r, rg, F, left=0.0, right=1.0)
        raise ValueError(self.name)

    def annulus_fractions(self, faces_cm):
        """Fraction in each [faces[i], faces[i+1]] (the last absorbs r > faces[-1])."""
        F = self.cdf(np.asarray(faces_cm, dtype=float))
        out = np.diff(F)
        out[-1] += 1.0 - F[-1]
        return out

    def radius_containing(self, frac):
        """Radius (cm) holding ``frac`` of the pairs (bisection on log r)."""
        lo, hi = 1e-12, 1.0
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if float(self.cdf(mid)) < frac:
                lo = mid
            else:
                hi = mid
        return hi


def gaussian(w_cm=50e-7):
    """Gaussian radial profile exp(-r^2/w^2); default w = 50 nm."""
    return RadialModel("gaussian", (float(w_cm),))


def core_penumbra(r_min_cm, r_max_cm):
    if not 0 < r_min_cm < r_max_cm:
        raise ValueError("need 0 < r_min < r_max")
    return RadialModel("core_penumbra", (float(r_min_cm), float(r_max_cm)))


def whk(beta, rho=RHO_SIC):
    """Waligorski-Hamm-Katz radial dose shape at ion speed ``beta`` in a medium of
    density ``rho`` (g/cm^3; mass-thickness transfer from water)."""
    if not 0 < beta < 1:
        raise ValueError("beta must be in (0, 1)")
    return RadialModel("whk", (float(beta), float(rho)))


def thermalized(model, s_cm):
    """``model`` convolved in the plane with the Gaussian exp(-r^2/s^2) / (pi s^2): the
    pairs spread by ``s`` before drift-diffusion applies (hot-carrier thermalization /
    the finite range of the slowest delta electrons; the size of s is NOT known for
    4H-SiC and is a parameter of every result that depends on the core density).
    ``thermalized(m, s).cdf`` -> ``m.cdf`` as s -> 0."""
    if not s_cm > 0:
        raise ValueError("s_cm must be positive")
    return RadialModel("thermalized", (model, float(s_cm)))


@lru_cache(maxsize=64)
def _conv_table(model, s, n=700):
    """CDF of the convolution, from the shell masses of ``model`` on a log grid:
    rho'(r) = sum_j dF_j exp(-(r - r_j)^2 / s^2) i0e(2 r r_j / s^2) / (pi s^2)."""
    from scipy.special import i0e

    r_hi = max(model.radius_containing(1.0 - 1e-9), 10 * s)
    rb = np.geomspace(min(1e-10, s * 1e-3), r_hi, n)
    Fb = model.cdf(rb)
    rc = np.sqrt(rb[1:] * rb[:-1])
    dF = np.diff(Fb)
    dF[0] += Fb[0]
    ro = np.geomspace(min(1e-10, s * 1e-3), r_hi + 6 * s, n)
    arg = 2.0 * ro[:, None] * rc[None, :] / s**2
    ker = np.exp(-((ro[:, None] - rc[None, :]) ** 2) / s**2) * i0e(arg)
    rho = ker @ dF / (np.pi * s**2)
    integrand = 2 * np.pi * ro * rho
    F = np.concatenate(
        [[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(ro))]
    )
    F = F / F[-1]
    rg = np.concatenate([[0.0], ro])
    F = np.concatenate([[0.0], F])
    rg.setflags(write=False)
    F.setflags(write=False)
    return rg, F


def whk_ranges(beta, rho=RHO_SIC):
    """(omega_max MeV, R_max cm, theta cm) of the WHK model."""
    w = 2.0 * M_E_C2_MEV * beta * beta / (1.0 - beta * beta)
    if beta <= 0.031265:
        alpha, k = 1.079, 1.0355e-1  # kg m^-2 MeV^-alpha
    else:
        alpha, k = 1.667, 6.0138
    # kg m^-2 -> g cm^-2 (x 0.1), then / rho (g cm^-3) -> cm
    R = 0.1 * k * w**alpha / rho
    theta = 0.1 * k * (10e-6) ** alpha / rho
    return w, R, theta, alpha


@lru_cache(maxsize=64)
def _whk_table(beta, rho, n=4000):
    _, R, theta, alpha = whk_ranges(beta, rho)
    scale = RHO_WATER / rho  # water lengths -> this medium
    M = 1e-8 * scale  # 1e-10 m in cm
    Nn = (1.5e-7 + beta * 5e-7) * scale  # (1.5e-9 + beta 5e-9) m in cm
    L = (8.0 if beta <= 0.031265 else 19.0) * beta ** (1.0 / 3.0)
    r = np.geomspace(1e-10, R, n)
    base = (
        (1.0 / r**2)
        / alpha
        * np.clip(1.0 - (r + theta) / (R + theta), 0.0, None) ** (1.0 / alpha)
        / (1.0 + theta / r)
    )
    K = np.where(r < M, 0.0, L * ((r - M) / Nn) * np.exp(-(r - M) / Nn))
    dose = base * (1.0 + K)
    integrand = 2 * math.pi * r * dose
    F = np.concatenate(
        [[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r))]
    )
    F /= F[-1]
    rg = np.concatenate([[0.0], r])
    F = np.concatenate([[0.0], F])
    rg.setflags(write=False)
    F.setflags(write=False)
    return rg, F
