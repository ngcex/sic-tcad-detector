"""Condensed-history sampler (etna.core.mc_transport) against NIST / analytic limits.

Issue #7. Tolerances are set from the sample size (heavy-tailed losses: >= 5 sigma).
"""

import numba
import numpy as np
import pytest

from etna.core import mc_transport as mt
from etna.core import straggling as S
from etna.core.electron_range import (
    electron_csda_table,
    estar_default_csda,
    estar_total_stopping,
)
from etna.core.stopping_power import bohr_range_straggling_cm, csda_range_cm


@pytest.mark.parametrize("material", ["si", "c_amorphous", "water"])
def test_electron_csda_matches_nist_column(material):
    e, r = electron_csda_table(material)
    ed, rd = estar_default_csda(material)
    mine = np.exp(np.interp(np.log(ed), np.log(e), np.log(r)))
    sel = ed >= 0.02
    assert np.max(np.abs(mine[sel] / rd[sel] - 1)) < 0.005
    assert np.max(np.abs(mine / rd - 1)) < 0.01


def test_estar_sic_is_bragg_composition():
    e, s = estar_total_stopping("sic")
    _, s_si = estar_total_stopping("si")
    assert e[0] == pytest.approx(1e-3)
    assert np.all(s < np.maximum(s_si, estar_total_stopping("c_amorphous")[1] * 2))
    assert mt.electron_range_um(10.0, "sic") == pytest.approx(1.02, rel=0.05)


def test_material_sic_parameters():
    m = mt.get_material("sic")
    assert m.z_over_a == pytest.approx(0.4988, abs=2e-4)
    assert 130 < m.I_eV < 145  # Bragg log-mix of NIST I(Si)=173, I(C)=81
    assert m.zeff_mass == pytest.approx(11.60, abs=0.01)


def test_effective_charge_limits():
    sic = mt.get_material("sic")
    assert mt.ion_effective_charge("c12", [12 * 300.0], sic)[0] == 6.0
    q = mt.ion_effective_charge("c12", [12 * 10.0], sic)[0]
    assert 5.95 < q < 6.0
    qa = mt.ion_effective_charge("alpha", [5.486], sic)[0]
    assert 1.9 < qa < 2.05  # Ziegler's (1 + tt) factor can exceed 2 (as in Geant4)
    assert mt.ion_effective_charge("proton", [1.0], sic)[0] == 1.0


@pytest.mark.parametrize(
    "particle,E,tol",
    [
        ("proton", 62.0, 0.01),
        ("proton", 10.0, 0.005),
        ("alpha", 5.486, 0.002),
        ("c12", 12 * 290.0, 0.01),
    ],
)
def test_mean_deposit_matches_csda_integral(particle, E, tol):
    ev = mt.sample_track_events(particle, E, n_events=200_000, delta_rays=False, seed=1)
    ref = float(mt.mean_csda_deposit_keV(particle, E, 10.0))
    assert ev.edep_keV.mean() == pytest.approx(ref, rel=tol)


def test_deterministic_baseline_is_csda_mean():
    ev = mt.sample_track_events(
        "proton", 62.0, n_events=1000, delta_rays=False, fluctuations="none"
    )
    ref = float(mt.mean_csda_deposit_keV("proton", 62.0, 10.0))
    assert np.allclose(ev.edep_keV, ref, rtol=1e-3)
    assert np.allclose(ev.chord_um, 10.0)


@pytest.mark.parametrize(
    "particle,E", [("proton", 62.0), ("proton", 250.0), ("c12", 12 * 290.0)]
)
def test_variance_close_to_vavilov(particle, E):
    ev = mt.sample_track_events(particle, E, n_events=400_000, delta_rays=False, seed=2)
    vb = mt.bohr_variance_keV2(particle, E, 10.0)
    assert 0.8 < ev.edep_keV.var() / vb < 1.2


def _mode(x, bins=400):
    lo, hi = np.percentile(x, [0.5, 98.0])
    h, e = np.histogram(x, bins=bins, range=(lo, hi))
    hs = np.convolve(h.astype(float), np.ones(7) / 7, mode="same")
    return 0.5 * (e[1:] + e[:-1])[np.argmax(hs)]


def test_vavilov_regime_most_probable_value():
    # 10 MeV p in 10 um SiC: kappa = 0.54
    part, E, t = "proton", 10.0, 10.0
    ev = mt.sample_track_events(part, E, n_events=400_000, delta_rays=False, seed=3)
    kap = mt.vavilov_kappa(part, E, t)
    b2, _, _ = mt.kinematics(part, E)
    xi = mt.landau_xi_keV(part, E, t)
    d = np.linspace(60, 200, 2001)
    p = S.vavilov_pdf_delta(d, ev.edep_keV.mean(), xi, kap, b2)
    assert _mode(ev.edep_keV) == pytest.approx(d[np.argmax(p)], rel=0.03)


def test_landau_regime_wider_than_landau():
    # PDG RPP 2024 sec. 34: thin-Si distributions are wider than 4 xi, MPV close to Eq. 34.12
    ev = mt.sample_track_events(
        "proton", 150.0, n_events=400_000, delta_rays=False, seed=4
    )
    mpv = _mode(ev.edep_keV)
    ref = mt.landau_mpv_keV("proton", 150.0, 10.0)
    assert 1.0 < mpv / ref < 1.2
    assert np.percentile(ev.edep_keV, 75) - np.percentile(
        ev.edep_keV, 25
    ) > 4.018 * 0.5 * (mt.landau_xi_keV("proton", 150.0, 10.0))


def test_gaussian_regime_bohr_width():
    ev = mt.sample_track_events(
        "alpha", 5.486, sv_um=(100, 100, 2), n_events=200_000, delta_rays=False, seed=5
    )
    sd = np.sqrt(mt.bohr_variance_keV2("alpha", 5.486, 2.0))
    assert ev.edep_keV.std() == pytest.approx(sd, rel=0.1)


def test_isotropic_mean_chord_is_cauchy():
    ev = mt.sample_track_events(
        "proton", 62.0, sv_um=(100, 60, 10), irradiation="isotropic", n_events=200_000
    )
    V = 100 * 60 * 10
    Sa = 2 * (100 * 60 + 100 * 10 + 60 * 10)
    assert np.all(ev.chord_um > 0)
    assert ev.chord_um.mean() == pytest.approx(4 * V / Sa, rel=0.01)


def test_tilted_beam_chord():
    ev = mt.sample_track_events(
        "proton",
        62.0,
        sv_um=(1e4, 1e4, 10),
        irradiation=60.0,
        n_events=2000,
        delta_rays=False,
        fluctuations="none",
    )
    inner = np.abs(ev.entry_um[:, 0]) < 4000
    assert np.allclose(ev.chord_um[inner], 20.0, rtol=1e-9)


def test_delta_escape_and_equilibrium():
    ref = float(mt.mean_csda_deposit_keV("proton", 62.0, 10.0))
    bare = mt.sample_track_events("proton", 62.0, n_events=200_000, seed=6)
    assert bare.edep_keV.mean() / ref == pytest.approx(0.914, abs=0.015)
    assert np.allclose(bare.edep_keV, bare.edep_primary_keV + bare.edep_delta_keV)
    eq = mt.sample_track_events(
        "proton", 62.0, n_events=200_000, seed=6, surround_um=300.0
    )
    hit = eq.chord_um > 0
    # equilibrium recovers the unrestricted mean (+ ~1 % from slowing in the surround)
    assert eq.edep_keV.sum() / hit.sum() / ref == pytest.approx(1.01, abs=0.015)


def test_alpha_range_and_straggling():
    ev = mt.sample_track_events("alpha", 5.486, sv_um=(100, 100, 40), n_events=100_000)
    R = float(csda_range_cm("alpha", "sic", 5.486)) * 1e4
    assert np.nanmean(ev.stop_path_um) == pytest.approx(R, rel=0.005)
    sb = bohr_range_straggling_cm("alpha", 5.486) * 1e4
    assert np.nanstd(ev.stop_path_um) == pytest.approx(sb, rel=0.2)


def test_step_size_convergence():
    kw = dict(
        sv_um=(100, 100, 20), n_events=2000, delta_rays=False, fluctuations="none"
    )
    a = mt.sample_track_events("alpha", 5.486, step_fraction=0.4, **kw)
    b = mt.sample_track_events("alpha", 5.486, step_fraction=0.01, **kw)
    assert np.nanmean(a.stop_path_um) == pytest.approx(
        np.nanmean(b.stop_path_um), rel=2e-3
    )


def test_reproducible_and_thread_independent():
    a = mt.sample_track_events("proton", 62.0, n_events=5000, seed=9)
    b = mt.sample_track_events("proton", 62.0, n_events=5000, seed=9)
    assert np.array_equal(a.edep_keV, b.edep_keV)
    n0 = numba.get_num_threads()
    try:
        numba.set_num_threads(1)
        c = mt.sample_track_events("proton", 62.0, n_events=5000, seed=9)
    finally:
        numba.set_num_threads(n0)
    assert np.array_equal(a.edep_keV, c.edep_keV)
    h1 = mt.sample_track_events("proton", 62.0, n_events=2500, seed=9)
    h2 = mt.sample_track_events(
        "proton", 62.0, n_events=2500, seed=9, event_offset=2500
    )
    assert np.array_equal(np.concatenate([h1.edep_keV, h2.edep_keV]), a.edep_keV)
    d = mt.sample_track_events("proton", 62.0, n_events=5000, seed=10)
    assert not np.array_equal(a.edep_keV, d.edep_keV)


def test_hits_conserve_energy():
    ev = mt.sample_track_events("c12", 12 * 62.0, n_events=2000, max_hits=64, seed=3)
    assert np.allclose(ev.hits[:, :, 3].sum(axis=1), ev.edep_keV, rtol=1e-9, atol=1e-9)
    nh = ev.n_hits[:, 0]
    i = int(np.argmax(nh))
    h = ev.hits[i, : nh[i]]
    assert np.all((h[:, 2] >= 0) & (h[:, 2] <= 10.0))


def test_energy_spectrum_input():
    E = np.array([10.0, 100.0])
    ev = mt.sample_track_events(
        "proton", E, n_events=4, delta_rays=False, fluctuations="none"
    )
    assert np.allclose(ev.edep_keV[0::2], ev.edep_keV[0])
    assert ev.edep_keV[0] > 3 * ev.edep_keV[1]


@pytest.mark.parametrize(
    "kw",
    [
        dict(energy_MeV=0.0),
        dict(energy_MeV=np.nan),
        dict(sv_um=(1, 1)),
        dict(sv_um=(1, -1, 1)),
        dict(tcut_keV=0.5),
        dict(fluctuations="landau"),
        dict(irradiation=95.0),
        dict(step_fraction=0.0),
        dict(n_events=0),
        dict(surround_um=-1.0),
    ],
)
def test_validation(kw):
    args = dict(particle="proton", energy_MeV=62.0)
    args.update(kw)
    with pytest.raises(ValueError):
        mt.sample_track_events(**args)


def test_auto_tcut():
    tc = mt.auto_tcut_keV((100, 100, 10))
    assert mt.electron_range_um(tc) == pytest.approx(0.5, rel=0.02)
    assert mt.auto_tcut_keV((1, 1, 0.01)) == 1.0
