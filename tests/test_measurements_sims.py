"""etna.measurements SIMS/SRP: round trips on SYNTHETIC profiles generated in-test from
the model's own donor profile (etna.core.device numpy twins) plus an erfc Al p+ layer."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import erfc

import etna
import etna.measurements as M
from tests._measurement_synth import sims_csv, write

XJ = 0.3  # um, synthetic junction depth (= the default p+ thickness)


def _synthetic(cfg, step_um=0.01, depth_max=12.0, al=True):
    depth = np.arange(0.0, depth_max + 1e-9, step_um)
    N = M.model_donor_profile(cfg, depth - XJ)
    N = np.where(depth < XJ, N[np.argmin(np.abs(depth - XJ))], N)  # N continues into p+
    Al = 1e19 * 0.5 * erfc((depth - 0.27) / 0.01) + 1e12 if al else None
    return depth, N, Al


def test_interp_round_trip_is_exact():
    cfg = etna.measured_cv_tabulated_config()
    depth, N, _ = _synthetic(cfg, al=False)
    prof = M.load_sims_csv(sims_csv(depth, N, meta=[f"# junction_depth_um: {XJ}"]))
    knots = [k for k, _ in cfg.donor_table]
    t = M.resample_profile(prof, cfg, knots_um=knots, method="interp")
    assert t.junction_source == "file metadata"
    np.testing.assert_allclose(t.N_D, [n for _, n in cfg.donor_table], rtol=1e-9)
    new = t.to_device_config(etna.DeviceConfig())
    np.testing.assert_allclose(
        np.array(new.donor_table), np.array(cfg.donor_table), rtol=1e-9
    )
    assert new.substrate_thickness_um == XJ and new.buffer_thickness_um is None
    assert M.profile_vs_model(t, cfg)["max_abs_log10"] < 1e-9


def test_nm_and_um_depth_columns_agree():
    cfg = etna.DeviceConfig()
    depth, N, Al = _synthetic(cfg, step_um=0.05)
    a = M.load_sims_csv(sims_csv(depth, N, Al, unit="um"))
    b = M.load_sims_csv(sims_csv(depth, N, Al, unit="nm"))
    np.testing.assert_allclose(a.depth_um, b.depth_um, rtol=1e-14)


def test_junction_from_the_al_n_crossing_and_acceptor_plateau():
    cfg = etna.DeviceConfig()
    depth, N, Al = _synthetic(cfg)
    t = M.resample_profile(M.load_sims_csv(sims_csv(depth, N, Al)), cfg)
    assert t.junction_source == "Al/N crossing"
    assert abs(t.junction_depth_um - 0.27) < 0.04  # erfc centre 0.27, width 10 nm
    assert t.N_A_plateau == pytest.approx(1e19, rel=1e-3)
    new = t.to_device_config(cfg)
    assert new.N_A == pytest.approx(1e19, rel=1e-3)
    assert new.substrate_thickness_um == pytest.approx(t.junction_depth_um)


def test_bin_mean_with_noise_recovers_the_profile_within_its_uncertainty():
    # smooth SYNTHETIC profile (abrupt epi end, no table): bin mean is unbiased there
    cfg = etna.nplus_calibrated_config()
    depth, N, _ = _synthetic(cfg, step_um=0.005, al=False)
    rng = np.random.default_rng(3)
    noisy = N * np.exp(rng.normal(0.0, 0.05, N.size))  # 5 % log-normal scatter
    prof = M.load_sims_csv(
        sims_csv(
            depth, noisy, meta=[f"# junction_depth_um: {XJ}", "# rel_sigma_N: 0.1"]
        )
    )
    t = M.resample_profile(prof, cfg, spacing_um=0.5)
    ref = M.model_donor_profile(cfg, t.knots_um)
    inner = (t.knots_um > 0.0) & (t.knots_um < cfg.epi_thickness_um)
    err = np.abs(t.N_D / ref - 1.0)[inner]
    assert np.all(err < 2 * t.rel_uncertainty[inner])  # within 2 sigma
    assert np.all(t.rel_uncertainty[inner] >= 0.1)  # calibration sigma not / sqrt(n)
    assert np.all(t.n_rows[inner] >= 90)


def test_detection_limit_censors_knots():
    cfg = etna.measured_cv_tabulated_config()  # bulk 1.3e13
    depth, N, _ = _synthetic(cfg, al=False)
    prof = M.load_sims_csv(
        sims_csv(
            depth,
            N,
            meta=[f"# junction_depth_um: {XJ}", "# detection_limit_N_cm3: 3e13"],
        )
    )
    t = M.resample_profile(prof, cfg)
    assert t.censored.any()
    assert np.all(t.N_D[t.censored] <= 3e13 * (1 + 1e-9))
    with pytest.raises(ValueError, match="detection limit"):
        t.to_device_config(cfg)
    assert t.to_device_config(cfg, allow_censored=True).donor_table == t.donor_table


def test_short_profile_is_flagged_extrapolated():
    cfg = etna.DeviceConfig()
    depth, N, _ = _synthetic(cfg, depth_max=7.0, al=False)
    t = M.resample_profile(
        M.load_sims_csv(sims_csv(depth, N)), cfg, junction_depth_um=XJ
    )
    assert t.junction_source == "argument"
    assert t.extrapolated.any() and t.knots_um[t.extrapolated].min() > 6.6
    with pytest.raises(ValueError, match="outside the measured depth"):
        t.to_device_config(cfg, allow_extrapolated=False)


def test_json_round_trip_and_figure(tmp_path):
    from etna.measurements import plots
    from etna.viz import export_figure

    cfg = etna.DeviceConfig()
    depth, N, Al = _synthetic(cfg, step_um=0.02)
    prof = M.load_sims_csv(write(tmp_path, "sims.csv", sims_csv(depth, N, Al)))
    t = M.resample_profile(prof, cfg)
    p = tmp_path / "t.json"
    p.write_text(t.to_json(cfg))
    assert M.load_donor_table_json(p, epi_thickness_um=10.0) == t.donor_table
    assert M.read_donor_table_json(p)["data"]["source_file"] == "sims.csv"
    fig = plots.profile_figure(prof, t, config=cfg)
    assert export_figure(fig, "png", dpi=60)[:4] == b"\x89PNG"
    assert "SIMS/SRP profile" in M.report_markdown(profile_table=t)
