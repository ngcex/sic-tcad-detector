"""Physical-model library (etna.core.physics_models), impact-ionization numerics
(etna.core.avalanche), the coupled avalanche / high-field DD engine
(etna.core.batched.models_dd) and the facades (etna.api.models).

Gates (measured 2026-10-21 in brackets): the coupled solver with every model off
reproduces ContinuityBatch's current to round-off [<= 1e-15]; DD with avalanche vs the
McIntyre generation-weighted multiplication in (1 - 1/M) <= 1e-4 up to M ~ 400
[4.3e-6]; Hatakeyama's measured V_BD = 1940 (1e16/N_D)^0.8 V within 10 % [-0.6 / +1.8 /
+7.0 % at 3e16 / 1e17 / 2e17: consistency, the coefficients were fitted on these
diodes]; V_BD mesh refinement <= 0.01 % [0.0047 %].
"""

import dataclasses
import json

import numpy as np
import pytest

from etna import DeviceConfig
from etna.api.models import refine_mesh, run_avalanche, run_breakdown, run_transit_time
from etna.core.avalanche import (
    alpha_edges_kernel,
    breakdown_voltage,
    ionization_integrals,
)
from etna.core.batched.mesh import node_positions
from etna.core.batched.models_dd import ModelContinuityBatch, edge_mobility_kernel
from etna.core.batched.sweeps import _Ramp, alpha_generation, make_batch, srh_continuity
from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM
from etna.core.physics_models import (
    HIGH_FIELD_MOBILITY,
    IMPACT_IONIZATION,
    LOW_FIELD_MOBILITY,
    AlphaParams,
    PhysicsModels,
    hatakeyama_anisotropic_alpha,
)
from etna.core.sic_material import SiC4H_Parameters, mobility_caughey_thomas_T

SUB, NPL = 0.3e-4, 1e-4


# ------------------------------------------------------------------- registry


def test_default_selection_is_etna_physics():
    m = PhysicsModels()
    assert m.is_default
    assert m.impact_coefficients(300.0) is None
    assert m.high_field_params(300.0) is None
    assert hash(m) == hash(PhysicsModels())
    json.dumps(m.describe())


@pytest.mark.parametrize(
    "kw, msg",
    [
        ({"impact_ionization": "Hatakeyama2004"}, "not a registered model"),
        ({"high_field_mobility": 3}, "must be None"),
        ({"low_field_mobility": None}, "cannot be None"),
        ({"allow_T_extrapolation": 1}, "must be a bool"),
    ],
)
def test_strict_validation(kw, msg):
    with pytest.raises(ValueError, match=msg):
        PhysicsModels(**kw)


def test_basal_only_set_is_refused_for_the_vertical_diode():
    basal = dataclasses.replace(IMPACT_IONIZATION["hatakeyama2004"], n_par=None)
    with pytest.raises(ValueError, match="parallel-to-c"):
        PhysicsModels(impact_ionization=basal)


def test_temperature_range_is_enforced():
    with pytest.raises(ValueError, match="allow_T_extrapolation"):
        PhysicsModels(impact_ionization="hatakeyama2004").impact_coefficients(350.0)
    m = PhysicsModels(impact_ionization="hatakeyama2004", allow_T_extrapolation=True)
    assert m.impact_coefficients(350.0)[0] == pytest.approx(1.76e8)
    assert m.T_extrapolated(350.0)
    # niwa2015: hole a and b grow with T (c, d > 0), electrons T-independent
    c300 = PhysicsModels(impact_ionization="niwa2015").impact_coefficients(300.0)
    c400 = PhysicsModels(impact_ionization="niwa2015").impact_coefficients(400.0)
    assert c400[0] == c300[0] and c400[1] == c300[1]
    assert c400[4] == pytest.approx(3.14e6 * (1 + 6.3e-3 * 100))
    assert c400[5] == pytest.approx(11.8e6 * (1 + 1.23e-3 * 100))


def test_every_set_has_provenance():
    for reg in (LOW_FIELD_MOBILITY, HIGH_FIELD_MOBILITY, IMPACT_IONIZATION):
        for name, m in reg.items():
            assert m.name == name
            assert m.provenance.status in (
                "verified-primary", "verified-secondary", "assumption", "legacy"
            )  # fmt: skip
            assert m.provenance.source
    assert IMPACT_IONIZATION["ayalew2004_legacy"].provenance.status == "legacy"
    assert IMPACT_IONIZATION["hatakeyama2004"].provenance.doi == "10.1063/1.1784520"


@pytest.mark.parametrize("carrier", ["electron", "hole"])
def test_etna_default_low_field_is_the_solver_mobility(carrier):
    N = np.geomspace(1e13, 1e20, 15)
    for T in (150.0, 300.0, 450.0):
        ref = mobility_caughey_thomas_T(N, T, carrier)
        got = LOW_FIELD_MOBILITY["etna_default"].mobility(N, T, carrier)
        np.testing.assert_array_equal(got, ref)
    # ayalew2004 adds (T/300)^-0.5 on mu_min to the legacy (perpendicular) set: equal
    # at 300 K, different elsewhere
    a = LOW_FIELD_MOBILITY["ayalew2004"]
    leg = SiC4H_Parameters(mobility="legacy_perp")
    np.testing.assert_allclose(
        a.mobility(N, 300.0, carrier),
        mobility_caughey_thomas_T(N, 300.0, carrier, leg),
        rtol=1e-15,
    )
    assert a.mobility(1e20, 450.0, carrier) < mobility_caughey_thomas_T(
        1e20, 450.0, carrier, leg
    )


def test_high_field_kernel_matches_reference():
    hf = HIGH_FIELD_MOBILITY["lades2000"]
    mu0 = np.full((50, 1), 900.0)
    E = np.geomspace(1e2, 1e7, 50)[:, None]
    v, b = hf.vsat_beta(300.0, "electron")
    got = edge_mobility_kernel(mu0, E, np.array([v]), np.array([b]))
    np.testing.assert_allclose(got, hf.mobility(mu0, E, 300.0, "electron"), rtol=1e-14)
    # velocity saturates at v_sat
    assert got[-1, 0] * E[-1, 0] == pytest.approx(v, rel=1e-3)
    # T laws: v_sat ~ T^-0.44, beta ~ T^1
    v4, b4 = hf.vsat_beta(400.0, "electron")
    assert v4 == pytest.approx(2.2e7 * (4 / 3) ** -0.44)
    assert b4 == pytest.approx(1.2 * 4 / 3)


def test_alpha_kernel_matches_reference():
    E = np.geomspace(1e5, 1e7, 40)
    for name, m in IMPACT_IONIZATION.items():
        T = 300.0
        c = m.coefficients(T)
        got = alpha_edges_kernel(E[None, :], *(np.array([v]) for v in c[:4]))[0]
        np.testing.assert_allclose(got, m.alpha(E, T, "electron"), rtol=1e-13)
    # Hatakeyama Table I at 3 MV/cm along c
    a = IMPACT_IONIZATION["hatakeyama2004"].alpha(3e6, 300.0, "hole")
    assert a == pytest.approx(3.41e8 * np.exp(-2.5e7 / 3e6))
    assert AlphaParams(1.0, 1.0).alpha(0.5) == 0.0


def test_hatakeyama_interpolation_reduces_to_axes():
    m = IMPACT_IONIZATION["hatakeyama2004"]
    F = np.array([2e6, 3e6, 4e6])
    for carrier in ("electron", "hole"):
        np.testing.assert_allclose(
            hatakeyama_anisotropic_alpha(m, F, 0 * F, carrier),
            m.alpha(F, carrier=carrier, direction="c_axis"), rtol=1e-12,
        )  # fmt: skip
        np.testing.assert_allclose(
            hatakeyama_anisotropic_alpha(m, 0 * F, F, carrier),
            m.alpha(F, carrier=carrier, direction="basal"), rtol=1e-12,
        )  # fmt: skip
        mid = hatakeyama_anisotropic_alpha(m, F / np.sqrt(2), F / np.sqrt(2), carrier)
        assert np.all(np.isfinite(mid)) and np.all(mid > 0)


# --------------------------------------------------------- ionization integrals


def _const_field(alpha_n, alpha_p, W=10e-4, N=201):
    """Potential of a constant 3 MV/cm field and coefficients giving alpha_n, alpha_p."""
    x = np.linspace(0.0, W, N)
    psi = (-3e6 * x)[None, :]
    # alpha = a exp(-(b/E)^m) with b tiny: alpha = a
    coef = (alpha_n, 1e-9, 0.0, 1.0, alpha_p, 1e-9, 0.0, 1.0)
    return x, psi, coef


def test_integral_identity_and_constant_field_limits():
    x, psi, coef = _const_field(800.0, 2500.0)
    I = ionization_integrals(x, psi, coef)
    phiW = np.sum((I.alpha_n[0] - I.alpha_p[0]) * np.diff(x))
    assert 1 - I.I_n[0] == pytest.approx(np.exp(-phiW) * (1 - I.I_p[0]), rel=1e-12)
    # alpha = beta: I = alpha W and M(x) = 1 / (1 - alpha W) for every x
    a = 0.5 / 10e-4
    x, psi, coef = _const_field(a, a)
    I = ionization_integrals(x, psi, coef)
    assert I.I_n[0] == pytest.approx(0.5, rel=1e-12)
    assert I.I_p[0] == pytest.approx(0.5, rel=1e-12)
    np.testing.assert_allclose(I.M_x[0], 2.0, rtol=1e-12)
    # only electrons ionize: M_n = exp(alpha W), I_p = 0
    x, psi, coef = _const_field(a, 0.0)
    I = ionization_integrals(x, psi, coef)
    assert I.M_n[0] == pytest.approx(np.exp(0.5), rel=1e-12)
    assert I.I_p[0] == 0.0


def test_refine_mesh_bisects():
    x = np.array([0.0, 1.0, 3.0])
    np.testing.assert_array_equal(refine_mesh(x, 1), [0, 0.5, 1, 2, 3])
    assert refine_mesh(x, 2).size == 9


# ----------------------------------------------------------------- breakdown


@pytest.fixture(scope="module")
def hata_batch():
    x = node_positions(10e-4, SUB, NPL)
    ND = np.array([3e16, 1e17, 2e17])
    pb = make_batch(
        x, SUB, 1e19, ND, ND, 1e-4, 300.0, nplus=(SUB + 10e-4, 1e19), buffer=None
    )
    return pb, ND


def test_breakdown_vs_hatakeyama_measured(hata_batch):
    """Consistency with the measured V_BD(N_D) of the diodes the set was fitted on."""
    pb, ND = hata_batch
    coef = PhysicsModels(impact_ionization="hatakeyama2004").impact_coefficients(300.0)
    bd = breakdown_voltage(pb, coef)
    assert bd.converged.all()
    ref = 1940.0 * (1e16 / ND) ** 0.8
    np.testing.assert_allclose(bd.V_BD / ref, 1.0, atol=0.10)
    # the root: I_p = 1 at V_BD; peak field = SISPAD-2005 critical field within 5 %
    I = ionization_integrals(pb.x, bd.psi, coef)
    np.testing.assert_allclose(I.I_p, 1.0, atol=1e-5)
    Fc = 2.70e6 * (1e16 / ND) ** -0.1
    np.testing.assert_allclose(bd.E_max_V_cm / Fc, 1.0, atol=0.05)
    assert np.all(bd.quadrature_error_V < 0.1)


def test_breakdown_nearly_uniform_field_analytic():
    """p+ / 1e12 / n+: nearly constant field E = V/W; with alpha = beta = a exp(-b/E)
    the breakdown is a W exp(-b W / V) = 1, V = b W / ln(a W)."""
    W = 10e-4
    x = node_positions(W, SUB, NPL)
    pb = make_batch(
        x, SUB, 1e19, 1e12, 1e12, 1e-4, 300.0, nplus=(SUB + W, 1e19), buffer=None
    )
    a, b = 3.0e6, 1.5e7
    bd = breakdown_voltage(pb, (a, b, 0.0, 1.0, a, b, 0.0, 1.0))
    # the field also reaches ~0.1 um into the p+ (N_A- 1.3e18) and the n+: W_eff > W
    assert bd.V_BD[0] == pytest.approx(b * W / np.log(a * W), rel=0.03)


# --------------------------------------------------------- coupled DD engine


@pytest.fixture(scope="module")
def default_engine():
    from etna.api.models import _batch

    cfg = DeviceConfig()
    pb, x, kw = _batch(cfg)
    G = alpha_generation(x, SUB, ALPHA_RANGE_SIC_CM, 1e18, epi_end=SUB + 10e-4)
    cb = srh_continuity(pb, 300.0, kw["N_A"])
    return pb, x, G, cb


def test_models_off_reproduce_continuity_batch(default_engine):
    pb, x, G, cb = default_engine
    ramp = _Ramp(pb, 1e-10, 10.0)
    mc = ModelContinuityBatch(pb, cb.mu_n, cb.mu_p, cb.tau_n, cb.tau_p)
    for V in (2.0, 10.0, 60.0):
        ramp.to(V)
        tr = cb.solve(
            ramp.psi, V_cathode=V, G=G, tol=1e-10, poisson_tol=1e-10, max_iter=200
        )
        r = mc.solve_models(ramp.psi, V_cathode=V, G=G)
        assert r.converged.all()
        assert abs(r.J[0]) == pytest.approx(abs(tr.J_cathode[0]), rel=1e-12)
        assert r.J_spread[0] < 1e-12


def test_avalanche_dd_matches_mcintyre(default_engine):
    pb, x, G, cb = default_engine
    coef = PhysicsModels(impact_ionization="hatakeyama2004").impact_coefficients(300.0)
    V_BD = breakdown_voltage(pb, coef).V_BD[0]
    assert V_BD == pytest.approx(2225.78, abs=0.05)
    ramp = _Ramp(pb, 1e-10, 10.0)
    base = ModelContinuityBatch(pb, cb.mu_n, cb.mu_p, cb.tau_n, cb.tau_p)
    mc = ModelContinuityBatch(pb, cb.mu_n, cb.mu_p, cb.tau_n, cb.tau_p, ii=coef)
    for frac in (0.9, 0.99, 0.999):
        V = frac * V_BD
        ramp.to(V)
        r0 = base.solve_models(ramp.psi, V_cathode=V, G=G)
        r = mc.solve_models(ramp.psi, V_cathode=V, G=G)
        assert r.converged.all() and r.newton_iterations[0] <= 5
        gain = r.J[0] / r0.J[0]
        M = ionization_integrals(x, ramp.psi, coef).generation_weighted(x, G)[0]
        assert (1 - 1 / gain) == pytest.approx(1 - 1 / M, rel=1e-4)


def test_past_breakdown_is_not_converged(default_engine):
    pb, x, G, cb = default_engine
    coef = PhysicsModels(impact_ionization="hatakeyama2004").impact_coefficients(300.0)
    ramp = _Ramp(pb, 1e-10, 10.0)
    ramp.to(2300.0)
    mc = ModelContinuityBatch(pb, cb.mu_n, cb.mu_p, cb.tau_n, cb.tau_p, ii=coef)
    r = mc.solve_models(ramp.psi, V_cathode=2300.0, G=G)
    assert not r.converged.any()


# ------------------------------------------------------------------- facades


def test_run_breakdown_default_device():
    r = run_breakdown(DeviceConfig(), n_points=12)
    md = r.metadata
    assert md["V_BD"] == pytest.approx(2225.78, abs=0.05)
    assert md["mesh_error_V"] / md["V_BD"] < 1e-4
    assert md["converged"] and not md["field_extrapolated"]
    assert np.all(np.diff(r.y) > 0) and r.y[0] == 1.0 and r.y[-1] > 1e3
    assert md["physics_models"]["impact_ionization"]["name"] == "hatakeyama2004"
    assert "upper" in md["note"] or "at or below" in md["note"]
    json.dumps(md["physics_models"])


def test_run_avalanche_alpha_and_dark():
    cfg = DeviceConfig()
    m = PhysicsModels(impact_ionization="hatakeyama2004")
    a = run_avalanche(cfg, m, v_start=-100.0, n_points=5, approach_breakdown=True)
    md = a.metadata
    assert all(md["converged"])
    np.testing.assert_allclose(
        1 - 1 / md["gain"][1:], 1 - 1 / md["M_mcintyre"][1:], rtol=1e-4
    )
    assert a.y[-1] > 100.0  # CCE with gain, not clipped
    d = run_avalanche(
        cfg, m, v_start=-60.0, n_points=4, approach_breakdown=True, source="dark"
    )
    assert np.all(np.diff(d.y) > 0) and d.metadata["gain"][-1] > 100
    # points past V_BD are NaN with a reason
    b = run_avalanche(cfg, m, -2000.0, -2400.0, 3)
    assert np.isnan(b.y[-1]) and b.metadata["reasons"][-1] == "at or beyond V_BD"


def test_high_field_mobility_effects_are_small_on_cce_large_on_transit():
    cfg = DeviceConfig()
    hf = PhysicsModels(high_field_mobility="lades2000")
    c = run_avalanche(cfg, hf, -10.0, -100.0, 3)
    dcce = c.y - c.metadata["cce_baseline"]
    assert np.all(dcce < 0) and np.all(np.abs(dcce) < 1e-3)
    t = run_transit_time(cfg, hf, -20.0, -100.0, 3)
    ratio = t.metadata["t_e"] / t.metadata["t_e_low_field"]
    assert np.all(ratio > 1.5) and ratio[-1] > ratio[0]


def test_facades_refuse_what_they_do_not_carry():
    with pytest.raises(NotImplementedError, match="1D"):
        run_breakdown(DeviceConfig(half_width_um=50.0))
    with pytest.raises(NotImplementedError, match="bandgap_narrowing"):
        run_avalanche(DeviceConfig(), PhysicsModels(bandgap_narrowing="lades2000"))
    with pytest.raises(ValueError, match="allow_T_extrapolation"):
        run_breakdown(DeviceConfig(T=350.0))
    with pytest.raises(ValueError, match="impact_ionization"):
        run_breakdown(DeviceConfig(), PhysicsModels())
