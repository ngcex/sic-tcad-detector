"""Axisymmetric single-ion track engine (etna.core.track): radial models, grid, the 1D
limits against the batched engine, ledgers, time-step convergence, opt-in physics and
the lattice heat equation. All on MESH_TINY (a few seconds each)."""

import math

import numpy as np
import pytest

from etna import DeviceConfig
from etna.api.device import device_kwargs
from etna.core.batched.material import Q
from etna.core.batched.sweeps import _Ramp, make_batch, srh_continuity
from etna.core.track.grid import axial_nodes, radial_nodes
from etna.core.track.phd import (
    MESH_TINY,
    build_track_device,
    ion_let_profile,
    simulate_track,
    track_pairs,
    whk_along,
)
from etna.core.track.radial import (
    beta_of,
    core_penumbra,
    gaussian,
    thermalized,
    whk,
    whk_ranges,
)
from etna.core.track.solver import W0, W2, dark_state, run_transient, steady_cce
from etna.core.track.thermal import C_V, LatticeHeat, line_source_peak

V_BIAS = 30.0


@pytest.fixture(scope="module")
def tiny():
    dev = build_track_device(DeviceConfig(), mesh=MESH_TINY)
    return dev, dark_state(dev, V_BIAS)


# ------------------------------------------------------------------ radial models
def test_gaussian_and_core_penumbra_cdfs():
    g = gaussian(50e-7)
    assert float(g.cdf(50e-7)) == pytest.approx(1 - math.exp(-1))
    cp = core_penumbra(1e-7, 1e-4)
    assert float(cp.cdf(1e-4)) == pytest.approx(1.0)
    assert float(cp.cdf(0.5e-7)) == pytest.approx(0.25 / (1 + 2 * math.log(1e3)))
    faces = np.concatenate([[0.0], np.geomspace(1e-8, 1e-5, 30)])
    for m in (g, cp, whk(0.05), thermalized(whk(0.05), 5e-7)):
        fr = m.annulus_fractions(faces)
        assert fr.sum() == pytest.approx(1.0, abs=1e-12)
        assert np.all(fr >= -1e-15)


def test_whk_shape():
    b = beta_of(1.0)  # 1 MeV/u
    assert b == pytest.approx(0.0463, rel=1e-3)
    m = whk(b)
    _, R, _, _ = whk_ranges(b)
    assert float(m.cdf(R)) == pytest.approx(1.0)
    r50, r90 = m.radius_containing(0.5), m.radius_containing(0.9)
    assert 1e-8 < r50 < 1e-7 < r90 < R  # half the energy within ~1 nm, tail to R_max
    # faster ions: wider penumbra
    assert whk(beta_of(10.0)).radius_containing(0.9) > r90


def test_thermalized_gaussian_is_gaussian():
    t = thermalized(gaussian(1e-6), 1e-6)
    r = np.array([5e-7, 1e-6, 2e-6])
    np.testing.assert_allclose(
        t.cdf(r), gaussian(math.sqrt(2) * 1e-6).cdf(r), atol=1e-3
    )


# ------------------------------------------------------------------ grid
def test_grid_volumes_and_nodes(tiny):
    dev, _ = tiny
    L = dev.z[-1]
    assert dev.volumes().sum() == pytest.approx(math.pi * dev.R**2 * L, rel=1e-12)
    assert dev.column().volumes().sum() == pytest.approx(dev.volumes().sum(), rel=1e-12)
    r = radial_nodes(1e-7, 1e-4, 1.5)
    assert r[0] == 0 and r[1] == 1e-7 and r[-1] == 1e-4
    z = axial_nodes([0.3e-4, 10.3e-4], 11.3e-4)
    assert 0.3e-4 in z and 10.3e-4 in z and z[-1] == 11.3e-4


def test_track_pairs_total(tiny):
    dev, _ = tiny
    pairs = track_pairs(dev, 100.0, gaussian(50e-7))
    epi = (dev.z > dev.junction_cm) & (dev.z < dev.epi_end_cm)
    expected = 100.0 * 1e7 / 7.8 * np.dot(dev.axial_dual(), epi)
    assert pairs.sum() == pytest.approx(expected, rel=1e-12)
    pk = pairs.reshape(dev.K, dev.M)
    assert pk[0].sum() == 0 and pk[-1].sum() == 0


# ------------------------------------------------------------------ 1D limits
def _batched(dev, V):
    kw = device_kwargs(DeviceConfig())
    epi_end = kw["substrate_thickness_cm"] + kw["epi_thickness_cm"]
    pb = make_batch(
        dev.z, kw["substrate_thickness_cm"], kw["N_A"], kw["N_D_junction"],
        kw["N_D_bulk"], kw["L_transition"], kw["T"], nplus=(epi_end, kw["N_D_nplus"]),
        buffer=(epi_end, kw["buffer_thickness_cm"], kw["buffer_end_doping"]),
    )  # fmt: skip
    ramp = _Ramp(pb, 1e-10, 10.0)
    ramp.to(V)
    return pb, ramp, srh_continuity(pb, kw["T"], kw["N_A"])


def test_dark_state_matches_batched_engine():
    dev = build_track_device(DeviceConfig(), mesh=MESH_TINY, auger=False).column()
    x, xbc, s = dark_state(dev, V_BIAS)
    pb, ramp, cb = _batched(dev, V_BIAS)
    tr = cb.solve(
        ramp.psi, V_cathode=V_BIAS, tol=1e-12, max_iter=200, poisson_tol=1e-10
    )
    assert np.abs(s.psi_ref + x[0::3] - tr.psi[0]).max() < 1e-10
    for c, ref in ((1, tr.n[0]), (2, tr.p[0])):
        m = ref > 1e5
        assert np.max(np.abs(x[c::3][m] - ref[m]) / ref[m]) < 1e-9


def test_steady_cce_matches_batched_engine():
    dev = build_track_device(DeviceConfig(), mesh=MESH_TINY, auger=False).column()
    x, xbc, s = dark_state(dev, V_BIAS)
    pb, ramp, cb = _batched(dev, V_BIAS)
    G = np.where((dev.z > dev.junction_cm) & (dev.z < dev.epi_end_cm), 1e18, 0.0)
    G[0] = G[-1] = 0.0
    tg = cb.solve(ramp.psi, V_cathode=V_BIAS, G=G, tol=1e-12, max_iter=200,
                  poisson_tol=1e-10)  # fmt: skip
    cce_ref = abs(tg.J_cathode[0]) / (Q * np.dot(dev.axial_dual(), G))
    assert steady_cce(s, x, xbc, G) == pytest.approx(cce_ref, abs=1e-7)


def test_low_injection_axisym_equals_column(tiny):
    """Low injection: the laterally integrated track is the 1D column problem.

    NOT the steady CCE of the same profile: in the 4H-SiC depletion region n_i ~ 5e-9,
    so SRH of the excess, n p / (tau_p n + tau_n p), needs both carriers at the same
    place and time; it is homogeneous but not linear, and a transient (electrons and
    holes separate) recombines less than steady illumination (here 5.9e-5 vs 2.0e-4 of
    the pairs; the steady value itself matches the 1D engine, test above)."""
    dev, dark = tiny
    r2 = simulate_track(device=dev, dark=dark, LET_keV_um=0.1, radial=gaussian(3e-4))
    col = dev.column()
    xc, bc, sc = dark_state(col, V_BIAS)
    pairs = track_pairs(dev, 0.1, gaussian(3e-4)).reshape(dev.K, dev.M).sum(axis=1)
    r1 = run_transient(sc, xc, bc, pairs)
    assert r1.converged and r2.converged
    assert r2.cce == pytest.approx(
        r1.cce, abs=1e-5
    )  # majority float64 floor ~2e-6 here
    G = pairs / sc.vol * (1e18 / np.max(pairs / sc.vol))
    assert steady_cce(sc, xc, bc, G) < r1.cce


# ------------------------------------------------------------------ ledgers / numerics
def test_ledgers_close_at_high_let(tiny):
    dev, dark = tiny
    r = simulate_track(device=dev, dark=dark, LET_keV_um=300.0, radial=gaussian(20e-7),
                       energy=True)  # fmt: skip
    assert r.converged
    assert r.ledger_residual < 1e-6
    assert r.extra["quad_err"] < 1e-6
    assert abs(r.cce - r.cce_anode) < 1e-5
    e = r.extra["energy"]
    assert e["residual"] < 1e-3  # discrete Poynting theorem, O(h^3) time error
    # Joule heat per collected charge = bias + built-in-like drop (V)
    v_eff = e["joule_J"] / (Q * r.cce * r.N_gen)
    assert V_BIAS < v_eff < V_BIAS + 4.0


def test_time_step_convergence(tiny):
    dev, dark = tiny
    kw = dict(device=dev, dark=dark, LET_keV_um=1000.0, radial=gaussian(20e-7))
    lo = simulate_track(rtol=1e-3, **kw)
    hi = simulate_track(rtol=1e-5, **kw)
    assert hi.n_steps > lo.n_steps
    assert abs(lo.cce - hi.cce) < 2e-4
    assert lo.recombined_fraction > 1e-5  # the plasma regime is exercised


def test_generation_pulse_matches_instantaneous(tiny):
    dev, dark = tiny
    kw = dict(device=dev, dark=dark, LET_keV_um=10.0, radial=gaussian(50e-7))
    a = simulate_track(**kw)
    b = simulate_track(pulse_sigma_s=2e-12, **kw)
    assert b.N_gen == pytest.approx(a.N_gen, rel=1e-3)  # scheme-integrated pulse
    assert abs(a.cce - b.cce) < 1e-4


# ------------------------------------------------------------------ opt-in physics
def test_injection_auger_and_carrier_scattering():
    cfg = DeviceConfig()
    base = dict(LET_keV_um=1000.0, V_bias=V_BIAS, radial=whk(beta_of(1.0)))
    out = {}
    for name, kw in (
        ("const", {}),
        ("inj", {"auger_model": "injection"}),
        ("ccs", {"carrier_scattering": True}),
    ):
        from etna.core.track.grid import device_from_config

        m = dict(MESH_TINY)
        r = radial_nodes(m.pop("r_min_cm"), m.pop("R_cm"), m.pop("r_growth"))
        dev = device_from_config(cfg, r=r, **m, **kw)
        out[name] = simulate_track(device=dev, **base)
    assert out["inj"].recombined_fraction < out["const"].recombined_fraction
    assert out["ccs"].recombined_fraction != out["const"].recombined_fraction
    for r in out.values():
        assert r.ledger_residual < 1e-6


def test_ion_profile_alpha():
    let, epu = ion_let_profile("alpha", 5.486)
    z = np.array([0.0, 5e-4, 17e-4, 30e-4])
    L = let(z)
    assert 100 < L[0] < 250 and L[2] > L[0] and L[3] == 0.0  # Bragg rise, then range
    m = whk_along(epu)(5e-4)
    assert m.name == "whk"


# ------------------------------------------------------------------ heat equation
def test_heat_line_source_analytic():
    """Isotropic heat diffusion of an instantaneous Gaussian line source (no electrical
    sources) vs the analytic axis temperature; energy conserved to round-off."""
    from etna.core.track.grid import device_from_config

    dev = device_from_config(
        DeviceConfig(), r=radial_nodes(2e-7, 20e-4, 1.25), h_min_cm=1e-5,
        h_max_cm=2e-4, z_growth=2.0,
    )  # fmt: skip
    x, xbc, s = dark_state(dev, 0.0)
    k = 4.0
    heat = LatticeHeat(s, x, kappa_r=k, kappa_z=k)
    w = 20e-7
    E_line = 1e-9  # J/cm
    rf = dev.radial_faces()
    frac = gaussian(w).annulus_fractions(rf)
    E_node = (E_line * dev.axial_dual())[:, None] * frac[None, :]
    heat.T = heat.T + E_node.ravel() / heat.CV
    E0 = heat.stored()
    t, h = 0.0, 1e-13
    checks = {}
    while t < 2e-10:
        heat.advance(s, x, x, x, h, (W0, W0, W2))
        t += h
        h *= 1.2
        for tc in (2e-11, 2e-10):
            if tc not in checks and t >= tc:
                checks[tc] = (t, heat.T.reshape(dev.K, dev.M)[dev.K // 2, 0] - 300.0)
    for tc, (tt, dT) in checks.items():
        assert dT == pytest.approx(line_source_peak(E_line, w, tt, kappa=k), rel=0.03)
    assert heat.stored() == pytest.approx(E0, rel=1e-9)


def test_electrothermal_coupled_run(tiny):
    dev, dark = tiny
    kw = dict(device=dev, dark=dark, LET_keV_um=300.0, radial=gaussian(20e-7))
    iso = simulate_track(**kw)
    hot = simulate_track(thermal={"couple": True}, **kw)
    th = hot.extra["thermal"]
    assert th["ledger_residual"] < 1e-6
    assert th["dT_max_K"] > 0
    assert th["heat_prompt_J"] == pytest.approx(
        hot.N_gen * (7.8 - 3.26) * Q, rel=0.01
    )  # (E_pair - E_g) per pair
    assert abs(hot.cce - iso.cce) < 1e-4  # 20 nm core: lattice heating is small
    assert C_V == pytest.approx(0.69 * 3.21)
    # cathode heat sink: the sink heat is the integrated conductive flux (independent
    # of the stored heat), so the ledger checks the scheme's conservation
    sink = simulate_track(thermal={"sink": "cathode"}, **kw).extra["thermal"]
    assert sink["heat_sink_J"] > 0
    assert sink["ledger_residual"] < 1e-6


def test_invalid_arguments(tiny):
    dev, _ = tiny
    with pytest.raises(ValueError):
        track_pairs(dev, 1.0, region="bulk")
    with pytest.raises(ValueError):
        gaussian(1e-6) and thermalized(gaussian(1e-6), 0.0)
    with pytest.raises(ValueError):
        build_track_device(DeviceConfig(), mesh=dict(MESH_TINY)) and __import__(
            "etna.core.track.grid", fromlist=["x"]
        ).device_from_config(DeviceConfig(), auger_model="bogus")
