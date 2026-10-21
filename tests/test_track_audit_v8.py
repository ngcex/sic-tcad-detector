"""Regression tests of the v8 audit fixes (.planning/AUDIT_v8_2D_TRACK.md M1, m2-m5,
m7; issue #21) and of issue #16 (tilted-alpha plasma slices and DeviceConfig.mobility).
Each test failed before its fix (except the analytic J.E check, which adds the
physics check of the heat source that m2 found missing)."""

import math

import numpy as np
import pytest

from etna import DeviceConfig, microdosimeter_config
from etna.core.batched.material import Q
from etna.core.sic_material import srh_lifetime
from etna.core.track import phd
from etna.core.track.cross_section import CrossSectionDevice, slice_device
from etna.core.track.phd import MESH_TINY, build_track_device, simulate_track
from etna.core.track.radial import gaussian
from etna.core.track.solver import W0, W2, TrackNotConverged, TrackSystem, dark_state
from etna.core.track.thermal import C_V, LatticeHeat
from etna.core.track.tilted import (
    PlasmaTable,
    column_device,
    path_profile,
    tilted_alpha_cce,
)

COLUMN_KW = dict(h_min_cm=1e-6, h_max_cm=1e-4, z_growth=1.6)


# ------------------------------------------------------------ M1 (issue #21)
def test_set_temperature_scales_the_device_lifetimes():
    dev = column_device(DeviceConfig(), tau_scale=0.05)
    s = TrackSystem(dev)
    s.set_temperature(np.full(s.N, dev.T))
    assert np.allclose(s.taun, dev.tau_n, rtol=1e-12)  # was the library 1 ns
    assert np.allclose(s.taup, dev.tau_p, rtol=1e-12)  # was the library 600 ns
    T1 = dev.T + 50.0
    s.set_temperature(np.full(s.N, T1))
    for tau, d, c in ((s.taun, dev.tau_n, "electron"), (s.taup, dev.tau_p, "hole")):
        law = srh_lifetime(T1, c) / srh_lifetime(dev.T, c)
        assert tau[0] / d == pytest.approx(law, rel=1e-12)


def test_coupled_run_at_vanishing_heating_equals_isothermal_with_tau_scale():
    """The audit's reproduction: column, tau_scale 0.05, 50 V, LET 0.1, Gaussian 50 nm.
    The lattice heats by ~1e-10 K, so the coupled run must equal the isothermal one
    (recombined fraction 4.05e-4; it was 2.07e-5, the unscaled-lifetime value)."""
    dev = column_device(DeviceConfig(), tau_scale=0.05)
    dark = dark_state(dev, 50.0)
    kw = dict(device=dev, dark=dark, LET_keV_um=0.1, radial=gaussian(50e-7),
              V_bias=50.0, t_max=5e-6)  # fmt: skip
    iso = simulate_track(**kw, excess=False)
    hot = simulate_track(**kw, thermal={"couple": True})
    assert hot.extra["thermal"]["dT_max_K"] < 1e-8
    assert iso.recombined_fraction == pytest.approx(4.05e-4, rel=0.02)
    assert hot.recombined_fraction == pytest.approx(iso.recombined_fraction, rel=1e-3)


# ------------------------------------------------------------ m5
def test_simulate_track_restores_isothermal_coefficients_on_failure(monkeypatch):
    dev = build_track_device(DeviceConfig(), mesh=MESH_TINY)
    dark = dark_state(dev, 10.0)
    system = dark[2]
    mun0, taun0 = system.mun.copy(), system.taun.copy()

    def failing(sysx, *a, thermal=None, **k):
        sysx.set_temperature(np.full(sysx.N, dev.T + 100.0))  # a hot, coupled step
        raise TrackNotConverged("forced")

    monkeypatch.setattr(phd, "run_transient", failing)
    with pytest.raises(TrackNotConverged):
        simulate_track(device=dev, dark=dark, LET_keV_um=10.0, V_bias=10.0,
                       thermal={"couple": True})  # fmt: skip
    assert np.array_equal(system.mun, mun0)
    assert np.array_equal(system.taun, taun0)


# ------------------------------------------------------------ m2: J.E source physics
def test_joule_source_matches_an_analytic_resistor():
    """Independent check of the J.E heat source (sign, units, edge power): a neutral
    n-type slab (N_D 1e16) under a uniform field E is an exact discrete steady state
    (uniform n, linear psi); the Joule power per unit length is sigma E^2 x area with
    sigma = q (mu_n N_D + mu_p n_i^2 / N_D), and the adiabatic lattice warms uniformly
    by sigma E^2 h / C_v in a step h."""
    d = slice_device(0.0, R_cm=3e-4, r_min_cm=2e-6, growth=1.5, background="neutral",
                     N_D=1e16)  # fmt: skip
    x, _, s = d.dark()
    E = 1e3
    s.set_reference(-E * np.tile(d.u, d.Mv))
    mun, mup = d.node_mobility()
    sigma = Q * (mun[0] * d.N_D + mup[0] * d.n_i**2 / d.N_D)
    P = sigma * E**2 * (d.u[-1] - d.u[0]) * d.v[-1]  # W per cm of slab
    assert float(s.joule(x).sum()) == pytest.approx(P, rel=1e-10)
    heat = LatticeHeat(s, x)
    _, joule, recomb = heat.sources(x)
    assert joule == pytest.approx(P, rel=1e-10) and recomb == 0.0
    h = 1e-9
    heat.advance(s, x, x, x, h, out_weights=(W0, W0, W2))
    dT = heat.T - d.T
    assert np.allclose(dT, sigma * E**2 * h / C_V, rtol=1e-6)


# ------------------------------------------------------------ m3 / m4 / #16 (tilted)
def _toy_table(**settings):
    nodes = np.array([3e-4, 10e-4, 16e-4])
    Eg = np.array([0.0, 1e3, 1e4, 1e5])
    ftab = np.array([[1e-3, 3e-3, 1e-3, 1e-4]] * 3) * np.array([[1.0], [1.5], [2.0]])
    return PlasmaTable(nodes, np.array([220.0, 300.0, 450.0]), np.array([1.2, 0.8, 0.3]),
                       Eg, ftab, np.array([0.0, 1e-4, 1e-5, 0.0]),
                       settings=dict(settings))  # fmt: skip


def test_tilted_validity_flag_per_bias():
    with pytest.warns(RuntimeWarning, match="outside its validity"):
        r = tilted_alpha_cce(microdosimeter_config(), (5.0, 50.0), table=_toy_table(),
                             column_kw=COLUMN_KW)  # fmt: skip
    # audit table (microdosimeter, 65 deg, 18 um): 0.66 of the pairs below 1e4 V/cm at
    # 5 V, none at >= 20 V
    assert r.frac_below_1e4[0] == pytest.approx(0.66, abs=0.05)
    assert r.frac_below_1e4[1] == 0.0
    assert list(r.valid) == [False, True]
    # the flag agrees with the per-field loss split (bands < 2e3, 2e3-1e4, > 1e4)
    lb = r.loss_by_field()
    assert lb[1, :2].sum() == 0.0 and lb[0, :2].sum() > 0.0


def test_precomputed_table_is_checked_against_the_run():
    cfg = microdosimeter_config()
    dev = column_device(cfg, tau_scale=1.0, **COLUMN_KW)
    good = dict(theta_deg=65.0, s_nm=10.0, tau_n=dev.tau_n, tau_p=dev.tau_p,
                dead_layer_um=0.0, mobility="c_axis")  # fmt: skip
    r = tilted_alpha_cce(cfg, (50.0,), table=_toy_table(**good), column_kw=COLUMN_KW)
    assert r.valid.all()
    for bad, call in (
        (dict(tau_n=dev.tau_n * 25), dict(tau_scale=1.0)),  # table at another tau
        ({}, dict(tau_scale=0.04)),  # run at another tau
        (dict(s_nm=1.0), {}),
        (dict(theta_deg=60.0), {}),
        ({}, dict(dead_layer_um=0.5)),
        (dict(mobility="legacy_perp"), {}),
    ):
        with pytest.raises(ValueError, match="does not match"):
            tilted_alpha_cce(cfg, (50.0,), table=_toy_table(**{**good, **bad}),
                             column_kw=COLUMN_KW, **call)  # fmt: skip


def test_slice_device_uses_the_mobility_set():
    mu = {}
    for m in ("c_axis", "legacy_perp"):
        d = CrossSectionDevice(u=np.array([-1.0, 1.0]), v=np.array([0.0, 1.0]),
                               E_perp=0.0, mobility=m)  # fmt: skip
        mu[m] = d.node_mobility()[0][0]
    # issue #16: mu_n at 1e14 is 1162.4 (c_axis) vs 941.1 (legacy_perp)
    assert mu["c_axis"] == pytest.approx(1162.4, abs=0.1)
    assert mu["legacy_perp"] == pytest.approx(941.1, abs=0.1)


def test_tilted_forwards_the_config_mobility_to_the_slices(monkeypatch):
    from etna.core.track import tilted

    seen = []

    class _Fake:
        loss, ledger_residual, converged, wall_s = 1e-4, 0.0, True, 0.0

    def fake_slice_loss(E, LET, rad, **kw):
        seen.append(kw.get("mobility"))
        return _Fake()

    monkeypatch.setattr(tilted, "slice_loss", fake_slice_loss)
    cfg = microdosimeter_config(mobility="legacy_perp")
    r = tilted_alpha_cce(cfg, (50.0,), column_kw=COLUMN_KW, table_kw=dict(n_nodes=2))
    assert seen and set(seen) == {"legacy_perp"}
    assert r.table.settings["mobility"] == "legacy_perp"
    # and the stored mobility guards the table against a c_axis run
    with pytest.raises(ValueError, match="mobility"):
        tilted_alpha_cce(microdosimeter_config(), (50.0,), table=r.table,
                         column_kw=COLUMN_KW)  # fmt: skip


def test_path_profile_unchanged():  # guard: the toy tables above assume this path
    s, f, E = path_profile()
    assert s[-1] == pytest.approx(1.3 * 18e-4) and math.isclose(E[0], 5.5e6)


# ------------------------------------------------------------ m7: CCE clip flagged
def test_fv2d_cce_above_one_is_flagged_not_silently_clipped(monkeypatch):
    from etna.core.fv2d import cce_lateral_scan_fv2d, prepared_solver, workflows
    from etna.core.fv2d.structures import build_mesh, guard_ring_structure

    st = guard_ring_structure(sv_half_width_um=10.0)
    s = prepared_solver(build_mesh(st), {"cathode": 5.0}, dark_check=False)
    cut = st.info["sv_half_width_cm"] + 0.5 * st.info["gap_cm"]
    kw = dict(n_points=2, contact="anode", x_max_cm=st.info["sv_half_width_cm"],
              node_sets=workflows.cut_node_sets(s.mesh, cut))  # fmt: skip
    ref = cce_lateral_scan_fv2d(s, **kw)
    assert np.array_equal(ref["cce_raw"], ref["cce_values"])  # within [0, 1]: unchanged
    assert not ref["cce_clipped"].any()
    orig = workflows._currents_all

    def inflated(*a):  # a current-extraction error: 50 % too much current
        return {c: 1.5 * v for c, v in orig(*a).items()}

    monkeypatch.setattr(workflows, "_currents_all", inflated)
    with pytest.warns(RuntimeWarning, match="CCE above 1"):
        bad = cce_lateral_scan_fv2d(s, **kw)
    assert np.all(bad["cce_values"] <= 1.0)
    assert bad["cce_clipped"][0] and bad["cce_raw"][0] > 1.0


def test_engine2d_cce_above_one_is_flagged_not_silently_clipped(monkeypatch):
    from etna.core.engine2d.device import Device2D
    from etna.core.engine2d.scan import lateral_cce_scan_2d

    orig = Device2D.illuminated_currents
    monkeypatch.setattr(Device2D, "illuminated_currents",
                        lambda self, G: 1.5 * orig(self, G))  # fmt: skip
    with pytest.warns(RuntimeWarning, match="CCE above 1"):
        bad = lateral_cce_scan_2d(25.0, 20.0, n_lateral_points=2)
    assert bad["cce_clipped"].all() and np.all(bad["cce_raw"] > 1.0)
    assert np.array_equal(bad["cce_values"], np.minimum(bad["cce_raw"], 1.0))
    # the true CCE (raw / 1.5) is below 1: values inside [0, 1] are never altered
    assert np.all(bad["cce_raw"] / 1.5 < 1.0)
