"""etna.measurements C-V: round trips on SYNTHETIC data generated in-test from the model.

simulate (etna.run_cv) -> write CSV (tmp_path) -> load_cv_csv -> compare / invert ->
recover. Inverting a curve produced by the same model is an "inverse crime" (no model
error, no noise): the tolerances below test the machinery, not the method's accuracy on
a real device.
"""

from __future__ import annotations

import numpy as np
import pytest

import etna
import etna.measurements as M
from etna.core.batched.material import EPS
from tests._measurement_synth import cv_csv, write

TRUTH = etna.measured_cv_tabulated_config  # model preset used as the synthetic truth


@pytest.fixture(scope="module")
def truth_cv():
    cfg = TRUTH()
    r = etna.run_cv(cfg, 0.0, -60.0, 121)  # 0.5 V grid
    return cfg, r


def test_roundtrip_compare_is_exact_on_the_data_grid(truth_cv, tmp_path):
    cfg, r = truth_cv
    p = write(tmp_path, "cv.csv", cv_csv(r.x, r.y, meta=["# temperature_k: 300"]))
    d = M.load_cv_csv(p)
    assert d.area_cm2 == 0.04 and d.area_source == "metadata" and d.T_K == 300.0
    np.testing.assert_allclose(d.W_um, r.metadata["W_from_C"] * 1e4, rtol=1e-14)
    cmp = M.compare_cv(d, cfg)
    assert cmp.exact_grid and cmp.n_points == 121
    assert cmp.rms_W_um < 1e-9 and cmp.rms_C_rel < 1e-12
    assert cmp.model_metadata["solver"] == "batched"


def test_units_and_sign_conventions_are_equivalent(truth_cv):
    _, r = truth_cv
    ref = M.load_cv_csv(cv_csv(r.x, r.y))
    for unit in ("F", "nF"):
        d = M.load_cv_csv(cv_csv(r.x, r.y, unit=unit))
        np.testing.assert_allclose(d.C_per_area, ref.C_per_area, rtol=1e-13)
    d = M.load_cv_csv(cv_csv(r.x, r.y, convention="reverse_positive"))
    np.testing.assert_array_equal(d.V, ref.V)
    # per-area capacitance in nF/cm2, area as a column (mm2)
    lines = ["area_mm2,V_V,C_per_area_nF_cm2"]
    lines += [f"4,{v:.17g},{c / 0.04 * 1e9:.17g}" for v, c in zip(r.x, r.y)]
    d = M.load_cv_csv("\n".join(lines))
    assert d.area_cm2 == pytest.approx(0.04) and d.area_source.startswith("column")
    np.testing.assert_allclose(d.W_um, ref.W_um, rtol=1e-13)
    # area argument used only when the file has none
    d = M.load_cv_csv(cv_csv(r.x, r.y, area_cm2=None), area_cm2=0.04)
    assert d.area_source == "argument"
    np.testing.assert_allclose(d.W_um, ref.W_um, rtol=1e-14)


def test_nonuniform_biases_use_the_interpolated_model(tmp_path):
    cfg = TRUTH()
    # model solved exactly at biases that are NOT on a uniform grid from 0 V
    off = etna.run_cv(cfg, -0.123, -40.0, 300)
    d = M.load_cv_csv(cv_csv(off.x, off.y))
    cmp = M.compare_cv(d, cfg)  # 0.05 V model grid + linear interpolation
    assert not cmp.exact_grid
    # interpolation error of the 0.05 V model grid, largest at the knee
    assert cmp.max_abs_W_um < 0.01 and cmp.rms_W_um < 2e-3
    assert any("interpolated" in n for n in cmp.notes)
    fine_cmp = M.compare_cv(d, cfg, step_V=0.005)
    assert fine_cmp.max_abs_W_um < 1e-4


def test_compare_separates_two_models(truth_cv):
    _, r = truth_cv
    d = M.load_cv_csv(cv_csv(r.x, r.y))
    cmp = M.compare_cv(d, etna.DeviceConfig())  # parametric default vs the table
    assert 0.02 < cmp.rms_W_um < 0.2
    assert cmp.max_abs_W_um > 0.2  # the knee: parametric form misses it


def test_forward_points_are_dropped_and_noted(truth_cv):
    cfg, r = truth_cv
    V = np.r_[0.5, r.x]
    C = np.r_[r.y[0] * 1.2, r.y]
    cmp = M.compare_cv(M.load_cv_csv(cv_csv(V, C)), cfg)
    assert cmp.n_forward_dropped == 1 and cmp.n_points == 121
    assert any("forward" in n for n in cmp.notes)


def test_sigma_of_W_from_sigma_of_C(truth_cv):
    _, r = truth_cv
    d = M.load_cv_csv(cv_csv(r.x, r.y, sigma_rel=0.01))
    np.testing.assert_allclose(d.sigma_W_um(), 0.01 * d.W_um, rtol=1e-12)
    cmp = M.compare_cv(d, TRUTH())
    assert cmp.chi2 == pytest.approx(0.0, abs=1e-12)


def test_depletion_profile_recovers_uniform_doping():
    cfg = etna.DeviceConfig(doping_profile="uniform", N_D=1e15)
    r = etna.run_cv(cfg, 0.0, -40.0, 81)  # 1e15: the 10 um epi is not depleted at 40 V
    dep = M.depletion_profile(M.load_cv_csv(cv_csv(r.x, r.y)))
    mid = (dep.depth_um > 2.0) & (dep.depth_um < 6.0)
    assert mid.sum() > 10
    np.testing.assert_allclose(dep.N_cm3[mid], 1e15, rtol=0.03)
    assert np.all(dep.debye_um[mid] < 0.2)


@pytest.fixture(scope="module")
def inversion(truth_cv):
    _, r = truth_cv
    d = M.load_cv_csv(cv_csv(r.x, r.y))
    # base geometry = the default stack (parametric doping is replaced by the table)
    return d, M.invert_cv(d, etna.DeviceConfig())


@pytest.mark.slow
def test_inversion_recovers_the_synthetic_table(inversion):
    d, inv = inversion
    truth = np.log10([n for _, n in TRUTH().donor_table])
    np.testing.assert_allclose(inv.knots_um, [k for k, _ in TRUTH().donor_table])
    assert inv.converged
    assert inv.rms_W_um < 0.005  # stated tolerance: 5 nm rms on eps*A/C
    # stated tolerance: 0.1 dex at the knots inside the probed depth range
    assert inv.constrained.sum() == 12  # 2 .. 9.75 um (W(0) = 1.68, W(-60) = 9.89 um)
    assert np.max(np.abs(inv.log10_N - truth)[inv.constrained]) < 0.1
    lo, hi = inv.probed_depth_um
    assert lo == pytest.approx(d.W_um[0]) and hi == pytest.approx(d.W_um[-1])
    assert inv.model_provenance["solver"] == "batched"
    assert any("outside the probed depth" in n for n in inv.notes)


@pytest.mark.slow
def test_inversion_json_round_trip(inversion, tmp_path):
    d, inv = inversion
    p = tmp_path / "donor_table.json"
    p.write_text(inv.to_json())
    rec = M.read_donor_table_json(p)
    assert rec["method"].startswith("cv_inversion")
    assert rec["data"]["source_sha256"] == d.sha256
    table = M.load_donor_table_json(p, epi_thickness_um=10.0)
    assert table == inv.donor_table
    cfg = etna.DeviceConfig(donor_table=table)
    cmp = M.compare_cv(d, cfg)
    assert cmp.rms_W_um == pytest.approx(inv.rms_W_um, rel=1e-6, abs=1e-9)
    assert "C-V inversion" in M.report_markdown(cv_comparison=cmp, cv_inversion=inv)


def test_inversion_refuses_a_device_outside_the_engine_box(truth_cv):
    _, r = truth_cv
    d = M.load_cv_csv(cv_csv(r.x, r.y))
    with pytest.raises(ValueError, match="validated domain"):
        M.invert_cv(d, etna.DeviceConfig(T=500.0))


def test_inversion_uses_the_file_temperature(truth_cv):
    _, r = truth_cv
    d = M.load_cv_csv(cv_csv(r.x, r.y, meta=["# temperature_k: 500"]))
    with pytest.raises(ValueError, match="T=500"):
        M.invert_cv(d, etna.DeviceConfig())


def test_W_uses_the_model_permittivity(truth_cv):
    _, r = truth_cv
    d = M.load_cv_csv(cv_csv(r.x, r.y))
    np.testing.assert_allclose(d.W_um, EPS / (r.y / 0.04) * 1e4, rtol=1e-14)


def test_figures_build_and_export(truth_cv):
    from etna.measurements import plots
    from etna.viz import export_figure

    cfg, r = truth_cv
    d = M.load_cv_csv(cv_csv(r.x, r.y))
    fig = plots.cv_comparison_figure(M.compare_cv(d, etna.DeviceConfig()))
    assert export_figure(fig, "png", dpi=60)[:4] == b"\x89PNG"


def test_stacked_batch_matches_run_cv_per_table():
    """The inversion's one-solve Jacobian batch reproduces etna.run_cv for each table."""
    from etna.measurements.cv import _stacked_W, default_knots_um, table_config

    k = default_knots_um(10.0)
    base = etna.DeviceConfig()
    cfgs = [
        table_config(base, k, np.full(len(k), 14.0)),
        table_config(base, k, np.linspace(15.5, 13.3, len(k))),
        table_config(base, k, np.r_[np.full(len(k) - 3, 13.5), 15.0, 16.0, 16.5]),
    ]
    V = np.linspace(0.0, -60.0, 61)
    W = _stacked_W(cfgs, 61, -60.0, V, True)
    for j, c in enumerate(cfgs):
        ref = etna.run_cv(c, 0.0, -60.0, 61).metadata["W_from_C"] * 1e4
        np.testing.assert_allclose(W[j], ref, rtol=0, atol=1e-9)


@pytest.mark.slow
def test_inversion_figure_and_fast_path(inversion):
    from etna.measurements import plots
    from etna.viz import export_figure

    d, inv = inversion
    assert not any("slow path" in n for n in inv.notes)
    fig = plots.inversion_figure(inv, M.depletion_profile(d))
    assert export_figure(fig, "png", dpi=60)[:4] == b"\x89PNG"


def test_inversion_flags_biases_beyond_the_cross_check():
    r = etna.run_cv(TRUTH(), 0.0, -120.0, 121)
    d = M.load_cv_csv(cv_csv(r.x, r.y))
    inv = M.invert_cv(d, etna.DeviceConfig(), max_nfev=2)  # notes only, not the fit
    assert any("beyond cross_checked_to_V = 100 V" in n for n in inv.notes)
