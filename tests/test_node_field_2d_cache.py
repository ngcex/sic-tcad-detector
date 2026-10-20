"""run_field 2D's node field with cached Delaunay triangulations is bit-identical to the
uncached plotting2d.py pattern, cold and on a cache hit."""

from __future__ import annotations

import numpy as np

from etna.api import simulation as sim
from etna.api.device import DeviceConfig, device_kwargs
from etna.core.engine2d import dark_state_2d


def _reference(x_coords, y_coords, potential):
    """The pre-2026-10-21 _node_field_2d (no cache)."""
    from scipy.interpolate import LinearNDInterpolator

    x_reg = np.linspace(x_coords.min(), x_coords.max(), 100)
    y_reg = np.linspace(y_coords.min(), y_coords.max(), 200)
    X_reg, Y_reg = np.meshgrid(x_reg, y_reg)
    V_grid = LinearNDInterpolator(list(zip(x_coords, y_coords)), potential)(
        X_reg, Y_reg
    )
    Ey, Ex = np.gradient(-V_grid, y_reg[1] - y_reg[0], x_reg[1] - x_reg[0])
    E = np.sqrt(Ex**2 + Ey**2)
    f = LinearNDInterpolator(list(zip(X_reg.ravel(), Y_reg.ravel())), E.ravel())(
        x_coords, y_coords
    )
    f[np.isnan(f)] = 0.0
    return f


def test_cached_node_field_is_bit_identical():
    cfg = DeviceConfig(half_width_um=25.0)
    sim._DELAUNAY_CACHE.clear()
    states = [
        dark_state_2d(cfg.half_width_um, v, **device_kwargs(cfg)) for v in (10.0, 30.0)
    ]
    for st in states:  # cold for the first, cache hit (same mesh) for the second
        x, y, psi = st["x"], st["y"], st["Potential"]
        np.testing.assert_array_equal(
            sim._node_field_2d(x, y, psi), _reference(x, y, psi)
        )
    assert len(sim._DELAUNAY_CACHE) == 2  # one mesh, one regular grid
    # a different mesh (wider device) gets its own entries
    st = dark_state_2d(50.0, 10.0, **device_kwargs(cfg))
    x, y, psi = st["x"], st["y"], st["Potential"]
    np.testing.assert_array_equal(sim._node_field_2d(x, y, psi), _reference(x, y, psi))
    assert len(sim._DELAUNAY_CACHE) == 4


def test_cache_is_bounded():
    sim._DELAUNAY_CACHE.clear()
    rng = np.random.default_rng(0)
    for _ in range(sim._DELAUNAY_CACHE_MAX + 3):
        pts = rng.random((30, 2))
        sim._delaunay(("t", pts.tobytes()), lambda: pts)
    assert len(sim._DELAUNAY_CACHE) == sim._DELAUNAY_CACHE_MAX
    sim._DELAUNAY_CACHE.clear()
