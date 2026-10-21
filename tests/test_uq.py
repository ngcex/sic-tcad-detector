"""etna.uq: registry, samplers, Sobol' estimators (Ishigami), engine evaluators, run_uq."""

from __future__ import annotations

import numpy as np
import pytest

from etna import DeviceConfig
from etna.uq import (
    DEFAULT_PARAMS,
    OBSERVABLES,
    ParameterSpace,
    SaltelliDesign,
    UncertainParameter,
    correlated_normal,
    default_parameters,
    indices_table,
    ishigami,
    ishigami_indices,
    run_uq,
    sample,
    sobol_indices,
    tied_removal,
)
from etna.uq import observables as O
from etna.uq.params import LABELS, R_N_NOMINAL, G_LIT_PER_NEQ_SIC_DAMAGE, resolve_space

CFG = DeviceConfig()


def _ishigami_space():
    return ParameterSpace(
        [UncertainParameter(f"x{i + 1}", 0.0, "uniform", -np.pi, np.pi, "model:x")
         for i in range(3)]
    )  # fmt: skip


def _ishigami_run(N, seed=1, n_boot=0):
    d = SaltelliDesign(_ishigami_space(), seed=seed)
    P = d.new_points(N)
    fA, fB, fAB = SaltelliDesign.split(ishigami(P), N, 3)
    return sobol_indices(fA, fB, fAB, n_boot=n_boot, seed=seed)


# ------------------------------------------------------------------ registry


def test_registry_provenance_and_scope():
    reg = default_parameters(CFG)
    for p in reg.values():
        assert p.source and p.label in LABELS and p.applies_to
    for obs, names in DEFAULT_PARAMS.items():
        assert obs in OBSERVABLES
        for n in names:
            assert obs in reg[n].applies_to, (obs, n)
    # DESIGN-1 L tolerance 0.46 +- 0.3 um reaches below the engine box (0.3 um)
    L = reg["L_transition_um"]
    assert L.lower == 0.3
    assert L.truncated_mass == pytest.approx((0.3 - (CFG.L_transition_um - 0.3)) / 0.6)
    assert L.support[0] == pytest.approx(0.3)
    # the d spread is a data spread, said so
    assert "not a fit error" in reg["d_removal"].source
    # tau_scale is for the track facades only
    with pytest.raises(ValueError, match="does not enter"):
        resolve_space("cce", CFG, ["tau_scale"])


def test_distributions_and_truncation():
    p = UncertainParameter(
        "a", 1.0, "normal", 1.0, 0.5, "model:a", lower=0.8, upper=1.1
    )
    u = np.linspace(0, 1, 101)
    x = p.ppf(u)
    assert x.min() >= 0.8 and x.max() <= 1.1 and np.all(np.diff(x) >= 0)
    q = UncertainParameter("b", 1.0, "loguniform", 0.1, 10.0, "model:b")
    assert q.ppf(0.5) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        UncertainParameter("c", 1.0, "uniform", 2.0, 1.0, "model:c")
    with pytest.raises(ValueError):
        UncertainParameter("c", 1.0, "beta", 0.0, 1.0, "model:c")


def test_correlated_draws_reproduce_covariance():
    cov = np.array([[0.09, 0.05], [0.05, 0.04]])
    X = correlated_normal([4.7, 1.0], cov, 4096, method="sobol", seed=3)
    assert np.allclose(np.cov(X.T), cov, rtol=0.05, atol=2e-3)
    base = [UncertainParameter("d", 4.7, "normal", 4.7, 0.3, "model:d"),
            UncertainParameter("s", 1.0, "normal", 1.0, 0.2, "model:s"),
            UncertainParameter("z", 0.0, "uniform", 0.0, 1.0, "model:z")]  # fmt: skip
    sp = ParameterSpace.with_covariance(base, ["d", "s"], [4.7, 1.0], cov)
    assert sp.factors == [[0, 1], [2]]
    assert sp.factor_names == ["d+s", "z"]
    Y = sample(sp, 4096, "sobol", seed=2)
    assert np.allclose(np.cov(Y[:, :2].T), cov, rtol=0.05, atol=2e-3)
    assert abs(np.corrcoef(Y[:, 0], Y[:, 2])[0, 1]) < 0.05
    # the Saltelli AB matrix swaps the whole correlated block
    d = SaltelliDesign(sp, seed=0)
    d.new_points(8)
    assert np.array_equal(d.AB[0][:, :2], d.B[:, :2])
    assert np.array_equal(d.AB[0][:, 2], d.A[:, 2])


def test_samplers_seeded():
    sp = _ishigami_space()
    for m in ("sobol", "lhs", "mc"):
        a, b = sample(sp, 64, m, seed=5), sample(sp, 64, m, seed=5)
        assert np.array_equal(a, b)
        assert not np.array_equal(a, sample(sp, 64, m, seed=6))
        assert a.min() >= -np.pi and a.max() <= np.pi


# ------------------------------------------------------------------ estimators


def test_ishigami_converges_to_closed_form():
    ref = ishigami_indices()
    assert ref["V"] == pytest.approx(7**2 / 8 + 0.1 * np.pi**4 / 5
                                     + 0.1**2 * np.pi**8 / 18 + 0.5)  # fmt: skip
    errs = []
    for N in (256, 1024, 4096):
        s = _ishigami_run(N)
        errs.append(max(np.abs(s.S1 - ref["S1"]).max(), np.abs(s.ST - ref["ST"]).max()))
    assert errs[-1] < 0.01
    assert errs[0] > errs[-1]


def test_ishigami_bootstrap_ci_and_scipy_crosscheck():
    from scipy.stats import sobol_indices as scipy_sobol
    from scipy.stats import uniform

    ref = ishigami_indices()
    s = _ishigami_run(1024, n_boot=300)
    cover = ((s.ST_ci[0] <= ref["ST"]) & (ref["ST"] <= s.ST_ci[1])).sum() + (
        (s.S1_ci[0] <= ref["S1"]) & (ref["S1"] <= s.S1_ci[1])).sum()  # fmt: skip
    assert cover >= 5  # of 6 (95 % intervals)
    assert 0 < s.max_ci_width < 0.3
    r = scipy_sobol(func=lambda x: ishigami(x.T), n=1024,
                    dists=[uniform(-np.pi, 2 * np.pi)] * 3, random_state=1)  # fmt: skip
    assert np.allclose(s.ST, r.total_order, atol=0.03)
    assert np.allclose(s.S1, r.first_order, atol=0.03)


def test_design_extension_is_the_larger_design_and_reproducible():
    sp = _ishigami_space()
    a = SaltelliDesign(sp, seed=4)
    a.new_points(64)
    a.new_points(64)
    b = SaltelliDesign(sp, seed=4)
    b.new_points(128)
    assert np.array_equal(a.A, b.A) and np.array_equal(a.AB, b.AB)
    assert np.array_equal(_ishigami_run(256, seed=9, n_boot=20).ST_ci,
                          _ishigami_run(256, seed=9, n_boot=20).ST_ci)  # fmt: skip


def test_unresolved_and_failed_rows():
    rng = np.random.default_rng(0)
    fA, fB = rng.random(64), rng.random(64)
    fAB = rng.random((2, 64))
    s = sobol_indices(fA * 1e-9, fB * 1e-9, fAB * 1e-9, n_boot=0, resolution=1e-6)
    assert s.below_resolution and np.all(np.isnan(s.ST))
    fA[3] = np.nan
    s = sobol_indices(fA, fB, fAB, n_boot=0)
    assert s.n_used == 63


# ------------------------------------------------------------------ evaluators


def _nominal(obs):
    sp = resolve_space(obs, CFG)
    return sp, sp.nominal()[None, :]


def test_evaluators_match_facades_at_nominal():
    from etna import run_cce, run_cv, run_dark_current, run_transit_time

    sp, v = _nominal("cv")
    r = run_cv(CFG, 0.0, -30.0, 31)
    assert np.allclose(O.eval_cv(CFG, v, sp).Y[0], r.y, rtol=1e-12, atol=0)
    sp, v = _nominal("cce")
    r = run_cce(CFG, -1.0, -30.0, 10)
    assert np.allclose(O.eval_cce(CFG, v, sp).Y[0], r.y, rtol=0, atol=1e-9)
    sp, v = _nominal("dark_current")
    r = run_dark_current(CFG, 0.0, -60.0, 7)
    assert np.allclose(O.eval_dark_current(CFG, v, sp).Y[0], r.y, rtol=1e-10, atol=0)
    sp, v = _nominal("transit_time")
    r = run_transit_time(CFG)
    assert np.allclose(O.eval_transit_time(CFG, v, sp).Y[0], r.metadata["t_e_low_field"],
                       rtol=1e-12, atol=0)  # fmt: skip


def test_v_fd_matches_run_cv_rule():
    from etna import run_cv
    from etna.api.simulation import _reach_through_V
    from etna.core.batched.sweeps import reverse_sweep
    from etna.uq.observables import _batch, fd_threshold_cm

    sp, v = _nominal("v_fd")
    vfd = O.eval_v_fd(CFG, v, sp).Y[0, 0]
    assert abs(vfd - run_cv(CFG, 0.0, -10.0, 11).metadata["reach_through_V"]) <= 0.05
    # same threshold as the facade's rule, on one grid
    V = -np.arange(0.0, 6.0, 0.25)
    pb, _, kws = _batch([CFG], [0])
    W = reverse_sweep(pb, V, kws[0]["substrate_thickness_cm"], kws[0]["epi_thickness_cm"],
                      capacitance=False).W[:, 0]  # fmt: skip
    k = np.flatnonzero(W >= fd_threshold_cm(CFG))[0]
    assert _reach_through_V(W, V, CFG) == pytest.approx(abs(V[k]))


def test_breakdown_compensation_microdosimetry_at_nominal():
    from etna import run_breakdown, run_microdosimetry
    from etna._paths import data_path
    from etna.api.simulation import compensation_fluence

    sp, v = _nominal("breakdown")
    vbd = O.eval_breakdown(CFG, v, sp).Y[0, 0]
    assert vbd == pytest.approx(run_breakdown(CFG, mesh_check=False).metadata["V_BD"],
                                abs=1e-6)  # fmt: skip
    sp, v = _nominal("compensation_fluence")
    phi = O.eval_compensation_fluence(CFG, v, sp).Y[0, 0]
    assert phi == pytest.approx(compensation_fluence(CFG, 62.0), rel=1e-12)
    sp, v = _nominal("microdosimetry")
    yF, yD = O.eval_microdosimetry(CFG, v, sp).Y[0]
    r = run_microdosimetry(CFG, data_path("synthetic_mc_events.csv"))
    assert yF == pytest.approx(r.metadata["y_F"], rel=0.02)  # facade: binned moments
    assert yD == pytest.approx(r.metadata["y_D"], rel=0.02)


def test_mobility_scales_enter_transport_only():
    sp = resolve_space("cce", CFG, ["mu_n_scale", "mu_p_scale"])
    vals = np.array([[1.0, 1.0], [1.08, 0.65]])
    cce = O.eval_cce(CFG, vals, sp, voltages=[-2.0, -20.0]).Y
    assert not np.allclose(cce[0], cce[1], rtol=0, atol=1e-7)
    tt = O.eval_transit_time(CFG, vals, sp, voltages=[-50.0]).Y
    assert tt[1, 0] == pytest.approx(tt[0, 0] / 1.08, rel=1e-9)
    # C-V is electrostatic: the factors are not even read
    spc = resolve_space("cv", CFG, ["N_D_bulk"])
    c = O.eval_cv(CFG, np.array([[CFG.N_D_bulk]] * 2), spc, voltages=[-5.0]).Y
    assert c[0, 0] == c[1, 0]


def test_tied_removal_keeps_the_per_proton_calibration():
    d0 = 4.69
    assert tied_removal(d0, R_N_NOMINAL) == pytest.approx(d0)
    for r_n in (1.02, 1.20):
        s = R_N_NOMINAL / r_n
        g = s * (G_LIT_PER_NEQ_SIC_DAMAGE + tied_removal(d0, r_n))
        assert g == pytest.approx(G_LIT_PER_NEQ_SIC_DAMAGE + d0)


def test_cce_fluence_matches_run_radiation_damage_at_nominal():
    from etna import run_radiation_damage

    sp, v = _nominal("cce_fluence")
    fl = np.array([1e12, 1e14])
    y = O.eval_cce_fluence(CFG, v, sp, fluences=fl).Y[0]
    r = run_radiation_damage(CFG, fluences=fl, V_bias=-40.0, proton_energy_MeV=62.0)
    assert np.allclose(y, r.y, rtol=0, atol=1e-12)


def test_out_of_box_samples_are_nan_not_devsim():
    p = UncertainParameter("N_D_bulk", CFG.N_D_bulk, "uniform", 1e13, 3e13,
                           "config:N_D_bulk", applies_to=("cv",))  # fmt: skip
    sp = ParameterSpace([p])
    out = O.eval_cv(CFG, np.array([[1.5e13], [2.5e13]]), sp, voltages=[0.0, -5.0])
    assert out.info["n_outside_box"] == 1
    assert np.all(np.isnan(out.Y[0])) and np.all(np.isfinite(out.Y[1]))


# ------------------------------------------------------------------ run_uq


def test_run_uq_bands_indices_reproducible():
    r = run_uq(CFG, "cv", params=["N_D_bulk", "L_transition_um"], n=8,
               voltages=[0.0, -2.0, -10.0], n_boot=50, seed=3)  # fmt: skip
    md = r.metadata
    assert r.sim_type == "uq" and md["solver"] == "batched"
    assert md["n_evaluations"] == (2 + 2) * 8 and md["n_failed"] == 0
    assert np.all(md["band_low"] <= r.y) and np.all(r.y <= md["band_high"])
    assert md["numerical_error"].shape == (3,)
    assert md["sobol"]["S1"].shape == (2, 3)
    rows = indices_table(r)
    assert [row["factor"] for row in rows] == ["N_D_bulk", "L_transition_um"]
    r2 = run_uq(CFG, "cv", params=["N_D_bulk", "L_transition_um"], n=8,
                voltages=[0.0, -2.0, -10.0], n_boot=50, seed=3)  # fmt: skip
    assert np.array_equal(md["outputs"], r2.metadata["outputs"])
    assert np.array_equal(
        md["sobol"]["ST_ci"], r2.metadata["sobol"]["ST_ci"], equal_nan=True
    )


def test_run_uq_adaptive_and_plain_designs():
    r = run_uq(
        CFG, "compensation_fluence", n=16, target_ci=0.08, n_max=4096, n_boot=100
    )
    ad = r.metadata["adaptive"]
    ns = [h["n"] for h in ad["history"]]
    assert ns[0] == 16 and all(b == 2 * a for a, b in zip(ns, ns[1:]))
    assert ad["converged"] and ad["achieved_ci"] <= 0.08
    names = r.metadata["factor_names"]
    st = dict(zip(names, r.metadata["sobol"]["ST"]))
    assert st["N_D_bulk"] > st["d_removal"] > 0  # +-20 % spec tolerance dominates
    for m in ("lhs", "mc", "sobol"):
        p = run_uq(CFG, "compensation_fluence", n=32, method=m)
        assert p.metadata["sobol"] is None and p.metadata["n_evaluations"] == 32
    with pytest.raises(ValueError, match="unknown observable"):
        run_uq(CFG, "nope")


def test_run_uq_model_form_breakdown_and_figures():
    from etna.uq.plots import uq_figure

    r = run_uq(CFG, "microdosimetry", n=8, n_boot=20)
    mf = r.metadata["model_form"]
    assert mf["mc_statistics_ci95"]["y_F"][0] < r.metadata["nominal"][0]
    assert uq_figure(r, "band").etna_kind == "uq_band"
    r = run_uq(CFG, "cce", params=["mu_p_scale", "epi_thickness_um"], n=4, n_boot=20,
               voltages=[-2.0, -10.0], numerical_error=False)  # fmt: skip
    assert uq_figure(r, "indices").etna_kind == "uq_indices"
    assert uq_figure(r, "band").etna_kind == "uq_band"
