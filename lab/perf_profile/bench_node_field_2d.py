"""Where does run_field 2D's _node_field_2d time go (Delaunay vs evaluation)?

usage: uv run python lab/perf_profile/bench_node_field_2d.py
"""

from __future__ import annotations

import time

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

from etna import set_devsim_quiet
from etna.core.engine2d import dark_state_2d
from etna.api.device import DeviceConfig, device_kwargs
from etna.api.simulation import _node_field_2d

set_devsim_quiet(True)
cfg = DeviceConfig(half_width_um=50.0)
st = dark_state_2d(cfg.half_width_um, 30.0, **device_kwargs(cfg))
x, y, pot = st["x"], st["y"], st["Potential"]
print("nodes", x.size, "unique x", np.unique(x).size, "unique y", np.unique(y).size)

t = time.perf_counter()
ref = _node_field_2d(x, y, pot)
print(f"_node_field_2d total {time.perf_counter() - t:.3f} s")

pts = np.column_stack([x, y])
t = time.perf_counter()
tri1 = Delaunay(pts)
t1 = time.perf_counter() - t
xr = np.linspace(x.min(), x.max(), 100)
yr = np.linspace(y.min(), y.max(), 200)
X, Y = np.meshgrid(xr, yr)
t = time.perf_counter()
V = LinearNDInterpolator(tri1, pot)(X, Y)
t2 = time.perf_counter() - t
t = time.perf_counter()
tri2 = Delaunay(np.column_stack([X.ravel(), Y.ravel()]))
t3 = time.perf_counter() - t
Ey, Ex = np.gradient(-V, yr[1] - yr[0], xr[1] - xr[0])
E = np.sqrt(Ex**2 + Ey**2)
t = time.perf_counter()
f = LinearNDInterpolator(tri2, E.ravel())(x, y)
t4 = time.perf_counter() - t
f[np.isnan(f)] = 0.0
print(
    f"Delaunay(mesh) {t1:.3f}  eval->grid {t2:.3f}  Delaunay(grid) {t3:.3f}  eval->nodes {t4:.3f}"
)
print("bit-identical with precomputed triangulations:", np.array_equal(f, ref))
t = time.perf_counter()
_ = list(zip(x, y))
print(f"list(zip) of nodes {time.perf_counter() - t:.3f}")
