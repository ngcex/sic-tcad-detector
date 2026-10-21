"""Event-by-event response (etna.core.mc_response), issue #7."""

import numpy as np
import pytest

from etna.core import mc_response as mr
from etna.core.single_particle import load_cce_let_table


def test_lineal_moments_exact():
    e = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    m = mr.lineal_moments(e, 2.0)
    y = np.array([1, 2, 3, 4]) / 2.0
    assert m["n"] == 4
    assert m["y_F"] == pytest.approx(y.mean())
    assert m["y_D"] == pytest.approx((y**2).mean() / y.mean())


def test_lineal_se_matches_spread_of_replicas():
    rng = np.random.default_rng(0)
    ests = []
    ses = []
    for _ in range(200):
        e = rng.gamma(2.0, 3.0, 2000)
        m = mr.lineal_moments(e, 1.0)
        ests.append(m["y_D"])
        ses.append(m["se_y_D"])
    assert np.std(ests) == pytest.approx(np.mean(ses), rel=0.15)


def test_default_l_bar():
    assert mr.default_l_bar_um((100, 100, 10), "normal") == 10.0
    assert mr.default_l_bar_um((100, 100, 10), "isotropic") == pytest.approx(
        16.6667, rel=1e-4
    )
    assert mr.default_l_bar_um((100, 100, 10), 60.0) == pytest.approx(20.0)


def test_adaptive_converges_and_is_batch_independent():
    a = mr.adaptive_lineal_energy("proton", 62.0, rel_tol=0.01, batch=5000, seed=3)
    assert a["converged"]
    assert a["se_y_D"] / a["y_D"] < 0.01
    b = mr.adaptive_lineal_energy(
        "proton",
        62.0,
        rel_tol=1e-9,
        batch=a["n_events"] // 2,
        max_events=a["n_events"],
        seed=3,
    )
    assert not b["converged"]
    assert b["y_D"] == pytest.approx(a["y_D"], rel=1e-12)
    with pytest.raises(ValueError):
        mr.adaptive_lineal_energy("proton", 62.0, rel_tol=0.0)


def test_pulse_heights_ideal_and_noise():
    e = np.full(20000, 50.0)
    ph = mr.pulse_heights(e, np.full(e.size, 10.0))
    assert np.allclose(ph["measured_keV"], 50.0)
    assert np.allclose(ph["LET_keV_um"], 5.0)
    ph = mr.pulse_heights(e, np.full(e.size, 10.0), fano=0.1, W_eV=7.8)
    n = 50e3 / 7.8
    assert ph["n_pairs"].std() == pytest.approx(np.sqrt(0.1 * n), rel=0.05)
    ph = mr.pulse_heights(
        e, np.full(e.size, 10.0), enc_electrons=500, threshold_keV=50.0
    )
    assert ph["charge_e"].std() == pytest.approx(500, rel=0.05)
    assert 0.4 < ph["n_triggered"] / e.size < 0.6


def test_pulse_heights_cce_table_and_let_modes():
    cce, _ = load_cce_let_table("data/cce_let_table_100um.json")
    e = np.array([5.0, 50.0, 500.0])
    c = np.array([10.0, 20.0, 0.0])
    ph = mr.pulse_heights(e, c, cce_interp=cce, sv_thickness_um=10.0)
    assert ph["LET_keV_um"] == pytest.approx([0.5, 2.5, 50.0])
    assert np.allclose(ph["CCE"], cce(ph["LET_keV_um"]))
    th = mr.pulse_heights(e, let_mode="thickness", sv_thickness_um=10.0)
    assert th["LET_keV_um"] == pytest.approx([0.5, 5.0, 50.0])
    for bad in (dict(let_mode="x"), dict(fano=2.0), dict(W_eV=0.0)):
        with pytest.raises(ValueError):
            mr.pulse_heights(e, c, **bad)
    with pytest.raises(ValueError):
        mr.pulse_heights(-e, c)
