"""Analytic straggling references (etna.core.straggling), issue #7."""

import numpy as np
import pytest

from etna.core import straggling as S


@pytest.mark.parametrize(
    "kappa,beta2", [(0.2, 0.3), (0.5, 0.1), (1.0, 0.5), (3.0, 0.2)]
)
def test_vavilov_cumulants(kappa, beta2):
    # k_n = xi T_max^(n-1) (1/(n-1) - beta^2/n); lambda in T_max units
    lam = np.linspace(-kappa * (1 + beta2 - S.EULER) - 8, 25, 6001)
    p = S.vavilov_pdf_lambda(lam, kappa, beta2)
    m0 = np.trapezoid(p, lam)
    m1 = np.trapezoid(p * lam, lam)
    var = np.trapezoid(p * (lam - m1) ** 2, lam)
    m3 = np.trapezoid(p * (lam - m1) ** 3, lam)
    assert m0 == pytest.approx(1.0, abs=1e-6)
    assert m1 == pytest.approx(-kappa * (1 + beta2 - S.EULER), abs=1e-6)
    assert var == pytest.approx(kappa * (1 - beta2 / 2), rel=1e-5)
    assert m3 == pytest.approx(kappa * (0.5 - beta2 / 3), rel=1e-4)


def test_vavilov_delta_normalised():
    d = np.linspace(0, 400, 4001)
    p = S.vavilov_pdf_delta(d, 117.0, 4.0, 0.5, 0.02)
    assert np.trapezoid(p, d) == pytest.approx(1.0, abs=1e-4)
    assert np.trapezoid(p * d, d) == pytest.approx(117.0, rel=1e-4)


def test_landau_mode_and_fwhm():
    d = np.linspace(-20, 120, 140001)
    p = S.landau_pdf_delta(d, 10.0, 2.0)
    i = int(np.argmax(p))
    assert d[i] == pytest.approx(10.0, abs=2e-3)
    half = p[i] / 2
    w = d[i:][p[i:] < half][0] - d[:i][p[:i] < half][-1]
    assert w == pytest.approx(4.018 * 2.0, rel=1e-3)


def test_gaussian():
    d = np.linspace(-50, 50, 10001)
    assert np.trapezoid(S.gaussian_pdf_delta(d, 0, 25.0), d) == pytest.approx(1.0)
