"""Analytic energy-loss straggling references: Landau, Vavilov, Bohr (Gaussian).

Independent of the sampler in :mod:`etna.core.mc_transport` (no shared code), used to
benchmark it. Free-electron theories: they ignore atomic binding, so for very thin
absorbers (kappa << 0.01) the real (Bichsel / Urban) distributions are wider than
Landau's and peak slightly higher; for kappa >~ 0.1 Vavilov is the reference and for
kappa >~ 10 the Gaussian of Bohr variance.

Vavilov density (P.V. Vavilov, Sov. Phys. JETP 5 (1957) 749), in the integral form of
the CERNLIB G116 / Schorr implementation:

    phi(lambda) = (1 / pi) exp(kappa (1 + beta^2 C)) *
                  int_0^inf exp(kappa f1(u)) cos(lambda u + kappa f2(u)) du
    f1(u) = beta^2 (ln u - Ci(u)) - cos u - u Si(u)
    f2(u) = u (ln u - Ci(u)) + sin u + beta^2 Si(u)
    lambda = (Delta - <Delta>) / T_max - kappa (1 + beta^2 - C),  C = Euler's constant

(lambda in units of T_max = xi / kappa.) The transcription is validated numerically in
the tests by the first three cumulants of Vavilov's distribution,
k_n = xi T_max^(n-1) (1/(n-1) - beta^2/n): normalisation, mean, variance
xi T_max (1 - beta^2/2) and third cumulant xi T_max^2 (1/2 - beta^2/3) all reproduced
to < 1e-6.
"""

from __future__ import annotations

import numpy as np
from scipy import integrate, special

EULER = 0.5772156649015329


def _f12(u, beta2):
    si, ci = special.sici(u)
    lc = np.log(u) - ci
    f1 = beta2 * lc - np.cos(u) - u * si
    f2 = u * lc + np.sin(u) + beta2 * si
    return f1, f2


def vavilov_pdf_lambda(lam, kappa, beta2, u_max=None, n=20001):
    """Vavilov density in the variable lambda (see module docstring)."""
    lam = np.atleast_1d(np.asarray(lam, dtype=float))
    if u_max is None:
        # exp(kappa f1) with f1 ~ -pi u / 2 for large u: integrate to e^-40
        u_max = 80.0 / (np.pi * kappa) + 50.0
    u = np.linspace(1e-12, u_max, n)
    f1, f2 = _f12(u, beta2)
    w = np.exp(kappa * f1)
    out = np.empty(lam.size)
    for i, lv in enumerate(lam):
        out[i] = integrate.trapezoid(w * np.cos(lv * u + kappa * f2), u)
    return np.exp(kappa * (1.0 + beta2 * EULER)) / np.pi * out


def vavilov_pdf_delta(delta_keV, mean_keV, xi_keV, kappa, beta2, **kw):
    """Vavilov density in the energy loss Delta (1/keV); T_max = xi / kappa."""
    d = np.asarray(delta_keV, dtype=float)
    tmax = xi_keV / kappa
    lam = (d - mean_keV) / tmax - kappa * (1.0 + beta2 - EULER)
    return vavilov_pdf_lambda(lam, kappa, beta2, **kw) / tmax


def landau_pdf_delta(delta_keV, mpv_keV, xi_keV):
    """Landau density with the PDG most probable value Delta_p and width xi.

    Classical Landau variable lambda = (Delta - Delta_p) / xi - 0.22278 (mode at
    -0.22278, FWHM 4.018). scipy.stats.landau uses x = (2/pi) (lambda - ln(pi/2))
    (its mode is at -0.4293 and FWHM 2.558 = 4.018 x 2/pi; checked in the tests).
    """
    from scipy.stats import landau

    lam = (np.asarray(delta_keV, dtype=float) - mpv_keV) / xi_keV - 0.22278
    x = (2.0 / np.pi) * (lam - np.log(np.pi / 2.0))
    return landau.pdf(x) * (2.0 / np.pi) / xi_keV


def gaussian_pdf_delta(delta_keV, mean_keV, var_keV2):
    d = np.asarray(delta_keV, dtype=float)
    return np.exp(-0.5 * (d - mean_keV) ** 2 / var_keV2) / np.sqrt(2 * np.pi * var_keV2)
