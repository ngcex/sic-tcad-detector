"""Tilted-track follow-up of issue #9: the excess-variable kernel (float64 floor), the
transverse cross-section solver, the 2.5D composition, the tau calibration and the
facades. Small meshes; each test a few seconds."""

import math
from decimal import Decimal, getcontext

import numpy as np
import pytest

from etna import DeviceConfig, microdosimeter_config
from etna.core.track.cross_section import (
    CrossSectionDevice,
    slice_device,
    slice_loss,
    slice_pairs,
)
from etna.core.track.grid import radial_nodes
from etna.core.track.kernels import bern_diff
from etna.core.track.phd import MESH_TINY, build_track_device, simulate_track
from etna.core.track.radial import gaussian, thermalized
from etna.core.track.solver import TrackSystem, dark_state, run_transient
from etna.core.track.tilted import (
    PlasmaTable,
    calibrate_tau_scale,
    path_profile,
    tilted_alpha_cce,
)

XS = dict(R_cm=3e-4, r_min_cm=2e-7, growth=1.6)  # small transverse slice (tests)


# ------------------------------------------------------------ excess unknowns
def test_bern_diff_matches_high_precision_reference():
    getcontext().prec = 50

    def B(x):
        return Decimal(1) if x == 0 else x / (x.exp() - 1)

    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(1500):
        x = float(rng.choice([1, -1]) * 10 ** rng.uniform(-8, 1.9))
        if rng.random() < 0.1:
            x = 0.0
        d = float(rng.choice([1, -1]) * 10 ** rng.uniform(-14, 0.3))
        ref = B(Decimal(x) + Decimal(d)) - B(Decimal(x))
        worst = max(worst, float(abs((Decimal(bern_diff(x, d)) - ref) / ref)))
    assert worst < 1e-13


@pytest.fixture(scope="module")
def tiny0():
    dev = build_track_device(DeviceConfig(), mesh=MESH_TINY)
    return dev, dark_state(dev, 0.0)


@pytest.mark.parametrize("let", [1.0, 1000.0])
def test_excess_unknowns_match_absolute_and_remove_the_floor(tiny0, let):
    dev, dark = tiny0
    r0 = simulate_track(device=dev, dark=dark, LET_keV_um=let, V_bias=0.0, excess=False)
    r1 = simulate_track(device=dev, dark=dark, LET_keV_um=let, V_bias=0.0, excess=True)
    assert r1.extra["excess"] and not r0.extra["excess"]
    assert r1.cce == pytest.approx(r0.cce, abs=max(5 * r0.ledger_residual, 1e-8))
    # ledger_floor is accumulated from the carrier content held in the unknowns (not a
    # constant): excess form ~1e-13 (MESH_TINY), below its a priori n_steps x 4.4e-16
    assert 0.0 < r1.extra["ledger_floor"] <= 4.4e-16 * r1.n_steps
    assert r1.ledger_residual < 1e-7
    if let == 1.0:  # the absolute form sits at the float64 floor of the 1e19 layers
        assert r0.extra["ledger_floor"] > 1e-8
        # the absolute ledger cannot beat its floor; the excess form breaks through it
        assert r0.ledger_residual > 0.5 * r0.extra["ledger_floor"]
        assert r1.ledger_residual < 0.2 * r0.extra["ledger_floor"]
        assert r1.ledger_residual < 0.1 * r0.ledger_residual


def test_excess_default_and_opt_out(tiny0):
    dev, dark = tiny0
    r = simulate_track(device=dev, dark=dark, LET_keV_um=10.0, V_bias=0.0)
    assert r.extra["excess"] is True  # default (no thermal / energy)
    r = simulate_track(device=dev, dark=dark, LET_keV_um=10.0, V_bias=0.0, energy=True)
    assert r.extra["excess"] is False  # energy theorem needs absolute states
    with pytest.raises(ValueError):
        simulate_track(device=dev, dark=dark, LET_keV_um=10.0, V_bias=0.0,
                       energy=True, excess=True)  # fmt: skip


# ------------------------------------------------------------ cross-section
def test_slice_dark_state_is_exact():
    for kw in (dict(E_perp=3e4), dict(E_perp=1e4, N_space=8e13),
               dict(E_perp=0.0, background="neutral")):  # fmt: skip
        d = slice_device(R_cm=3e-4, **kw)
        x, xbc, s = d.dark()
        _, _, phi, rpsi = s.evaluate(x, xbc, jac=False)
        inner = s.is_contact == 0
        # Poisson exact (linear / quadratic psi, FV Laplacian); carrier rates ~ 0
        assert np.max(np.abs(rpsi[inner])) < 1e-6 * np.max(s.vol) * 1e10
        s.set_excess(x)
        _, _, phi, rpsi = s.evaluate(np.zeros_like(x), np.zeros_like(x), jac=False)
        assert np.all(phi == 0.0) and np.all(rpsi == 0.0)
    with pytest.raises(ValueError):
        slice_device(1e4, background="neutral")


def test_slice_lumping_conserves_and_matches_exact_gaussian():
    d = slice_device(0.0, R_cm=2e-4, r_min_cm=2e-7, growth=1.4)
    w1, s = 30e-7, 20e-7
    exact = slice_pairs(d, 1e8, gaussian(math.hypot(w1, s)))
    sampled = slice_pairs(d, 1e8, thermalized(gaussian(w1), s), n_theta=128, n_r=800)
    assert exact.sum() == pytest.approx(0.5e8, rel=1e-12)
    assert sampled.sum() == pytest.approx(0.5e8, rel=1e-9)
    # same spatial distribution to sub-cell accuracy: first moment of r^2
    uu = np.tile(d.u, d.Mv)
    vv = np.repeat(d.v, d.Mu)
    m_e = np.dot(exact, uu**2 + vv**2) / exact.sum()
    m_s = np.dot(sampled, uu**2 + vv**2) / sampled.sum()
    assert m_s == pytest.approx(m_e, rel=0.05)


def test_slice_at_zero_field_matches_the_axisymmetric_radial_problem():
    """Same physics, different mesh: the Cartesian half-plane slice at E = 0 vs an
    axisymmetric (r) ring mesh of the same radius, intrinsic background, absorbing
    outer boundary, Gaussian track. SRH off (tau 1 s): the loss is the early Auger
    loss in the core, independent of the box shape (square vs disc)."""
    R, L = 3e-4, 300.0
    rad = gaussian(40e-7)
    xs = slice_loss(0.0, L, rad, R_cm=R, r_min_cm=2e-7, growth=1.3, tau_n=1.0, tau_p=1.0)
    # axisymmetric 1-row mesh (unit length along the track)
    ref = CrossSectionDevice(
        u=np.array([-1.0, 1.0]), v=np.array([0.0, 1.0]), E_perp=0.0, tau_n=1.0,
        tau_p=1.0,
    )
    r = radial_nodes(2e-7, R, 1.3)
    M = r.size
    rf = np.concatenate([[0.0], 0.5 * (r[:-1] + r[1:]), [R]])
    area = math.pi * (rf[1:] ** 2 - rf[:-1] ** 2)
    nb = -np.ones((M, 2), dtype=np.int64)
    geo = np.zeros((M, 2))
    for i in range(M):
        if i > 0:
            nb[i, 0], geo[i, 0] = i - 1, 2 * math.pi * rf[i] / (r[i] - r[i - 1])
        if i < M - 1:
            nb[i, 1], geo[i, 1] = i + 1, 2 * math.pi * rf[i + 1] / (r[i + 1] - r[i])
    is_c = np.zeros(M, dtype=np.int64)
    is_c[-1] = 1
    sysa = TrackSystem(ref, geometry=dict(vol=area, netdop=np.zeros(M), nb=nb, geo=geo,
                                          is_contact=is_c, anode=np.zeros(0),
                                          cathode=np.array([M - 1])))  # fmt: skip
    ref.node_mobility = lambda T=None: tuple(
        np.full(M, m[0]) for m in CrossSectionDevice.node_mobility(ref)
    )
    sysa.set_temperature(None)
    x = np.zeros(3 * M)
    x[1::3] = ref.n_i
    x[2::3] = ref.n_i
    line = L * 1e7 / 7.8
    pairs = line * rad.annulus_fractions(rf)
    pairs[-2] += pairs[-1]
    pairs[-1] = 0.0
    ra = run_transient(
        sysa, x, x.copy(), pairs, rtol=1e-4, t_max=2e-6, keep_series=False
    )
    assert xs.converged and ra.converged
    assert xs.recombined_fraction == pytest.approx(ra.recombined_fraction, rel=0.03)


def test_slice_linear_response_and_field_dependence():
    rad = gaussian(30e-7)
    lo = slice_loss(1e4, 0.001, rad, **XS).recombined_fraction
    lo2 = slice_loss(1e4, 0.01, rad, **XS).recombined_fraction
    assert lo2 == pytest.approx(lo, rel=0.02)  # linear: independent of LET
    hi = slice_loss(1e4, 300.0, rad, **XS)
    hi2 = slice_loss(4e4, 300.0, rad, **XS)
    assert hi.ledger_residual < 1e-6
    assert hi.recombined_fraction > 10 * lo2  # plasma excess
    assert hi2.recombined_fraction < hi.recombined_fraction  # more field, less loss


# ------------------------------------------------------------ composition / calibration
def _toy_table():
    s, f, E = path_profile()
    nodes = np.array([3e-4, 10e-4, 16e-4])
    Eg = np.array([0.0, 1e3, 1e4, 1e5])
    ftab = np.array([[1e-3, 3e-3, 1e-3, 1e-4]] * 3) * np.array([[1.0], [1.5], [2.0]])
    return PlasmaTable(nodes, np.array([220.0, 300.0, 450.0]), np.array([1.2, 0.8, 0.3]),
                       Eg, ftab, np.array([0.0, 1e-4, 1e-5, 0.0]))  # fmt: skip


def test_plasma_table_interpolation_and_roundtrip():
    t = _toy_table()
    assert t(np.array([3e-4]), np.array([1e4]))[0] == pytest.approx(1e-3 - 1e-5)
    mid = t(np.array([6.5e-4]), np.array([1e5]))[0]
    assert mid == pytest.approx(0.5 * (1e-4 + 1.5e-4))
    t2 = PlasmaTable.from_dict(t.to_dict())
    assert np.allclose(t2.f_ex, t.f_ex) and np.allclose(t2.s_cm, t.s_cm)


def test_tilted_composition_with_a_given_table():
    cfg = microdosimeter_config()
    with pytest.warns(RuntimeWarning, match="outside its validity"):  # 5 V
        r = tilted_alpha_cce(cfg, (5.0, 50.0), table=_toy_table(),
                             column_kw=dict(h_min_cm=1e-6, h_max_cm=1e-4,
                                            z_growth=1.6))  # fmt: skip
    assert np.all(r.cce_column > 0.99)
    assert np.all((r.plasma_loss > 0) & (r.plasma_loss < 3e-3))
    # at 50 V every path element sits in 1e4..1e5 V/cm, where the toy f_ex lies in
    # [1e-4, 2e-3 - 1e-5]: the pair-weighted loss must too (weights sum to 1)
    assert r.valid[1] and 1e-4 <= r.plasma_loss[1] <= 2e-3 - 1e-5
    assert r.plasma_loss[0] > r.plasma_loss[1]  # weaker field at 5 V
    assert np.all(r.column_ledger < 1e-6)


def test_calibrate_tau_scale_recovers_a_synthetic_truth():
    taus = np.array([1.0, 1 / 3, 1 / 10, 1 / 30, 1 / 100])
    V = np.array([0.0, 5, 10, 20, 30, 40, 50])

    def shape(t):
        loss = 0.3 * (1 / t) ** 0.5 / 30 / (1 + V / 4.0)
        y = 1 - loss
        return y / y[-1]

    shapes = np.array([shape(t) for t in taus])
    truth = 1 / 17.0
    cal = calibrate_tau_scale(taus, shapes, shape(truth) * 0.99, sigma=1e-3, T=300.0)
    assert cal.tau_scale == pytest.approx(truth, rel=0.05)
    assert cal.ci68[0] < truth < cal.ci68[1]
    assert cal.free_scale == pytest.approx(0.99, rel=1e-3)
    assert cal.label == "calibrated" and not cal.at_edge
    assert cal.tau_p_ns == pytest.approx(600 * cal.tau_scale, rel=1e-6)


# ------------------------------------------------------------ facades
def test_facades_smoke():
    import etna

    r = etna.run_track_event(DeviceConfig(), LET_keV_um=100.0, V_bias=10.0, mesh="tiny")
    assert r.sim_type == "track_event" and 0.9 < r.metadata["cce"] <= 1.0
    assert r.metadata["excess"] is True and r.metadata["ledger_residual"] < 1e-6
    with pytest.warns(RuntimeWarning, match="outside its validity"):  # 10 V
        t = etna.run_tilted_alpha_track(
            microdosimeter_config(), (10.0, 50.0), table=_toy_table()
        )
    assert t.sim_type == "cce_tilted_track" and np.all(t.x == [-10.0, -50.0])
    assert list(t.metadata["valid"]) == [False, True]
    assert 0.1 < t.metadata["frac_below_1e4"][0] < 0.4  # audit: 0.24 at 10 V
