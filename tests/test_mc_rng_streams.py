"""Per-event random streams of etna.core.mc_transport are independent (issue #17).

Before the fix the start state was a linear function of (seed, event): event k+1 drew
event k's numbers shifted by one, so consecutive deposits were correlated (0.72 for
100 MeV p in 10 um) and the iid standard errors were 1.6-2.3x too small.
"""

from __future__ import annotations

import numba
import numpy as np
import pytest

from etna.core import mc_response as mr
from etna.core import mc_transport as mt


def _draws(seed, event, n):
    st = mt._seed_state(seed, event)
    return np.array([mt._uniform(st) for _ in range(n)])


def test_streams_are_not_shifted_copies():
    for seed in (0, 1, 7):
        for k in (0, 5, 1000):
            a, b = _draws(seed, k, 8), _draws(seed, k + 1, 8)
            # no shared draws at any lag (the bug: a[1:] == b[:-1] bit for bit)
            assert np.intersect1d(a, b).size == 0


def test_seed_event_combination_is_not_symmetric():
    # a symmetric hash (seed ^ event) would give (3, 5) and (5, 3) the same stream
    assert not np.array_equal(_draws(3, 5, 4), _draws(5, 3, 4))
    assert not np.array_equal(_draws(0, 0, 4), _draws(1, 1, 4))


def test_normal_beam_entry_points_uncorrelated():
    ev = mt.sample_track_events("proton", 100.0, (100, 100, 10), n_events=20000)
    assert np.max(np.abs(ev.entry_um[1:, 0] - ev.entry_um[:-1, 1])) > 1.0


def test_lag1_autocorrelation_of_deposits():
    n = 20000
    ev = mt.sample_track_events(
        "proton", 100.0, (100, 100, 10), n_events=n, delta_rays=False
    )
    e = ev.edep_keV
    r = np.corrcoef(e[1:], e[:-1])[0, 1]
    assert abs(r) < 4.0 / np.sqrt(n)  # was 0.72


def test_seed_scatter_matches_reported_standard_error():
    """sd over seeds of y_F / y_D == the reported delta-method SE (ratio ~1).

    100 seeds x 20k events (62 MeV p, the USAGE setup): the sd of a sample sd from 100
    values is ~7 %, so [0.75, 1.3] is a > 3 sigma band. Before the fix: 1.99 (y_F),
    1.53 (y_D).
    """
    yF, yD, sF, sD = [], [], [], []
    for s in range(100):
        ev = mt.sample_track_events(
            "proton", 62.0, (100, 100, 10), n_events=20000, seed=s
        )
        m = mr.lineal_moments(ev.edep_keV, 10.0)
        yF.append(m["y_F"])
        yD.append(m["y_D"])
        sF.append(m["se_y_F"])
        sD.append(m["se_y_D"])
    rF = np.std(yF, ddof=1) / np.mean(sF)
    rD = np.std(yD, ddof=1) / np.mean(sD)
    assert 0.75 < rF < 1.3, rF
    assert 0.75 < rD < 1.3, rD


def test_per_event_reproducible_across_threads_and_batches():
    a = mt.sample_track_events("proton", 62.0, n_events=4000, seed=21, max_hits=8)
    n0 = numba.get_num_threads()
    try:
        numba.set_num_threads(1)
        b = mt.sample_track_events("proton", 62.0, n_events=4000, seed=21, max_hits=8)
    finally:
        numba.set_num_threads(n0)
    np.testing.assert_array_equal(a.edep_keV, b.edep_keV)
    np.testing.assert_array_equal(a.hits, b.hits)
    h1 = mt.sample_track_events("proton", 62.0, n_events=2000, seed=21)
    h2 = mt.sample_track_events(
        "proton", 62.0, n_events=2000, seed=21, event_offset=2000
    )
    np.testing.assert_array_equal(
        np.concatenate([h1.edep_keV, h2.edep_keV]), a.edep_keV
    )


@pytest.mark.parametrize("seed", [0, 2**32 - 1])
def test_seed_range(seed):
    ev = mt.sample_track_events("proton", 62.0, n_events=100, seed=seed)
    assert np.all(np.isfinite(ev.edep_keV))
