"""grid_from_devsim memo: a hit equals a fresh build, returns independent copies, and a
monkeypatched builder is still called (the builder is part of the key)."""

from __future__ import annotations

import numpy as np

from etna.api.device import DeviceConfig, device_kwargs
from etna.core import device2d
from etna.core.engine2d import grid as g


def _kw(**over):
    return dict(half_width_um=25.0, **device_kwargs(DeviceConfig(**over)))


def _same(a, b):
    ga, ia = a
    gb, ib = b
    for f in ("xs", "ys", "perm", "donors", "acceptors"):
        np.testing.assert_array_equal(getattr(ga, f), getattr(gb, f))
    assert ia.keys() == ib.keys()


def test_hit_equals_fresh_build_and_is_a_copy():
    g._GRID_CACHE.clear()
    kw = _kw()
    first = g.grid_from_devsim(**kw)
    assert len(g._GRID_CACHE) == 1
    second = g.grid_from_devsim(**kw)
    assert len(g._GRID_CACHE) == 1
    _same(first, second)
    fresh = g._grid_from_devsim(device2d.create_sic_2d_device, kw)
    _same(second, fresh)
    assert second[0] is not first[0] and second[1] is not first[1]
    second[1]["T"] = -1.0  # a caller mutating info does not reach the cache
    assert g.grid_from_devsim(**kw)[1]["T"] == first[1]["T"]
    g.grid_from_devsim(**_kw(T=350.0))  # other keywords: other entry
    assert len(g._GRID_CACHE) == 2


def test_monkeypatched_builder_is_called(monkeypatch):
    g.grid_from_devsim(**_kw())  # warm the cache with the real builder
    calls = []
    real = device2d.create_sic_2d_device

    def spy(*a, **k):
        calls.append(k)
        return real(*a, **k)

    monkeypatch.setattr(device2d, "create_sic_2d_device", spy)
    g.grid_from_devsim(**_kw())
    assert len(calls) == 1


def test_unhashable_keyword_bypasses_the_cache(monkeypatch):
    n = len(g._GRID_CACHE)
    calls = []
    monkeypatch.setattr(
        g, "_grid_from_devsim", lambda b, kw: calls.append(kw) or ("grid", {})
    )
    assert g.grid_from_devsim(half_width_um=25.0, weird={"a": 1}) == ("grid", {})
    assert len(calls) == 1 and len(g._GRID_CACHE) == n
