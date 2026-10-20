"""etna.measurements I(T): Arrhenius fit and model E_a, on SYNTHETIC data generated
in-test (an analytic Arrhenius law, or the ETNA dark-current model itself; the small
inline CSV literals are made-up SYNTHETIC numbers, not measurements)."""

from __future__ import annotations

import numpy as np
import pytest

import etna
import etna.measurements as M
from etna.measurements._common import K_B_EV
from tests._measurement_synth import it_csv

T = np.arange(280.0, 345.0, 5.0)


def _arrhenius(T, Ea, p=2.0, lnA=-10.0):
    return T**p * np.exp(lnA - Ea / (K_B_EV * T))


def test_noise_free_fit_is_exact():
    f = M.fit_arrhenius(T, _arrhenius(T, 1.234), T_power=2.0)
    assert f.E_a_eV == pytest.approx(1.234, abs=1e-9)
    assert f.ln_A == pytest.approx(-10.0, abs=1e-7)
    assert f.sigma_E_a_eV < 1e-9 and not f.sigma_from_data


def test_noisy_fit_recovers_E_a_within_its_uncertainty():
    rng = np.random.default_rng(11)
    I0 = _arrhenius(T, 1.0)
    hits, chi2 = 0, []
    for _ in range(200):
        I = I0 * (1.0 + 0.03 * rng.standard_normal(T.size))
        f = M.fit_arrhenius(T, I, sigma_I_A=0.03 * I0, T_power=2.0)
        hits += abs(f.E_a_eV - 1.0) < f.sigma_E_a_eV
        chi2.append(f.chi2 / f.dof)
    assert 0.6 < hits / 200 < 0.76  # 1-sigma coverage ~68 %
    assert 0.8 < np.mean(chi2) < 1.2


def test_prefactor_exponent_shifts_E_a_by_kT():
    I = _arrhenius(T, 1.0, p=0.0)
    f0 = M.fit_arrhenius(T, I, T_power=0.0)
    f2 = M.fit_arrhenius(T, I, T_power=2.0)
    shift = f0.E_a_eV - f2.E_a_eV
    assert shift == pytest.approx(2 * K_B_EV * T.mean(), rel=0.02)


def test_model_round_trip_and_comparison():
    cfg = etna.DeviceConfig()
    Im, solvers = M.model_current_vs_T(cfg, T, -30.0)
    assert set(solvers) == {"batched"}
    text = it_csv(T, Im, meta=["# bias_V: -30"])
    an = M.analyze_it(M.load_it_csv(text), cfg)
    (c,) = an.comparisons
    assert c.delta_eV == pytest.approx(0.0, abs=1e-6)
    # the model's E_a: E_g/2 by assumption, bent by n_i ~ T^1.5, E_g(T), Hurkx Gamma(T)
    assert 1.55 < c.model_fit.E_a_eV < 1.8
    assert c.E_g_half_eV == pytest.approx(1.63, abs=0.01)
    js = an.to_json()
    assert '"E_a_eV"' in js and '"source_sha256"' in js
    assert "activation energy" in M.report_markdown(it_analysis=an)


def test_several_biases_in_one_file_and_sign_convention():
    cfg = etna.DeviceConfig()
    Ts = np.array([290.0, 300.0, 310.0, 320.0])
    rows_T, rows_V, rows_I = [], [], []
    for V in (-10.0, -40.0):
        Im, _ = M.model_current_vs_T(cfg, Ts, V)
        rows_T += list(Ts)
        rows_V += [-V] * len(Ts)  # reverse bias written positive
        rows_I += list(Im)
    d = M.load_it_csv(
        it_csv(rows_T, rows_I, V=rows_V, meta=["# bias_convention: reverse_positive"])
    )
    assert d.biases() == [-10.0, -40.0]
    an = M.analyze_it(d, cfg)
    assert [f.bias_V for f in an.fits] == [-10.0, -40.0]
    assert all(abs(c.delta_eV) < 1e-6 for c in an.comparisons)


def test_celsius_column_and_no_bias():
    lines = ["T_C,I_nA", "10,1.0", "20,2.0", "30,4.5"]
    d = M.load_it_csv("\n".join(lines))
    np.testing.assert_allclose(d.T_K, [283.15, 293.15, 303.15])
    an = M.analyze_it(d, etna.DeviceConfig())
    assert an.comparisons == [] and any("no bias" in n for n in an.notes)


def test_arrhenius_figure():
    from etna.measurements import plots
    from etna.viz import export_figure

    d = M.load_it_csv(it_csv(T, _arrhenius(T, 1.0), meta=["# bias_V: -20"]))
    fig = plots.arrhenius_figure(M.analyze_it(d, None, T_power=2.0))
    assert export_figure(fig, "png", dpi=60)[:4] == b"\x89PNG"
