"""Tensor-product grid of ETNA's 2D devsim device, as numba-friendly arrays.

devsim's 2D mesh (etna.core.device2d.create_sic_2d_device) is a tensor grid
xs (nx) x ys (ny) with one diagonal per rectangle whose EdgeCouple is ~1e-33
(notes/LOG.md, STEP 0). The finite-volume scheme is therefore exactly the
5-point stencil on rectangular dual cells:

    horizontal edge (i,j)-(i+1,j): couple dy_dual[j], length hx[i]
    vertical   edge (i,j)-(i,j+1): couple dx_dual[i], length hy[j]
    node volume dx_dual[i] * dy_dual[j]

Unknowns are numbered k = j*nx + i (x fastest -> half-bandwidth nx). devsim's
node order is mapped by coordinates (``perm[devsim_index] = k``), never assumed.
"""

from __future__ import annotations

import uuid
from dataclasses import dataclass
from functools import lru_cache

import numpy as np


@dataclass(frozen=True)
class Grid2D:
    xs: np.ndarray  # (nx,) cm
    ys: np.ndarray  # (ny,) cm
    perm: np.ndarray  # (N,) devsim node index -> k
    tri: np.ndarray  # (T, 3) triangles in k numbering (devsim's elements)

    @property
    def nx(self):
        return len(self.xs)

    @property
    def ny(self):
        return len(self.ys)

    @property
    def N(self):
        return self.nx * self.ny

    @property
    def hx(self):
        return np.diff(self.xs)

    @property
    def hy(self):
        return np.diff(self.ys)

    @staticmethod
    def _dual(h, n):
        d = np.zeros(n)
        d[:-1] += h / 2
        d[1:] += h / 2
        return d

    @property
    def dxd(self):
        return self._dual(self.hx, self.nx)

    @property
    def dyd(self):
        return self._dual(self.hy, self.ny)

    def to_devsim(self, field_k):
        """(N,) array in k order -> devsim node order."""
        return np.asarray(field_k)[self.perm]

    def from_devsim(self, field_dev):
        out = np.empty(self.N)
        out[self.perm] = field_dev
        return out

    def integrate_tri(self, values_k):
        """etna integrate_over_mesh_2d: sum over triangles of area * vertex mean."""
        X = np.tile(self.xs, self.ny)
        Y = np.repeat(self.ys, self.nx)
        t = self.tri
        x0, x1, x2 = X[t[:, 0]], X[t[:, 1]], X[t[:, 2]]
        y0, y1, y2 = Y[t[:, 0]], Y[t[:, 1]], Y[t[:, 2]]
        area = 0.5 * np.abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        v = np.asarray(values_k)
        return float(np.sum((v[t[:, 0]] + v[t[:, 1]] + v[t[:, 2]]) / 3.0 * area))


def grid_from_nodes(x, y, elements):
    xs, ys = np.unique(x), np.unique(y)
    i = np.searchsorted(xs, x)
    j = np.searchsorted(ys, y)
    if not (np.array_equal(xs[i], x) and np.array_equal(ys[j], y)):
        raise ValueError("mesh is not a tensor-product grid")
    perm = j * len(xs) + i
    if len(np.unique(perm)) != len(x) or len(x) != len(xs) * len(ys):
        raise ValueError("mesh is not a full tensor-product grid")
    tri = perm[np.asarray(elements, dtype=int).reshape(-1, 3)]
    return Grid2D(xs=xs, ys=ys, perm=perm, tri=tri)


@lru_cache(maxsize=32)
def grid_from_devsim(
    half_width_um=50.0, epi_thickness_cm=10e-4, substrate_thickness_cm=1e-4
):
    """Build the devsim mesh only (no equations, no solve), read it, delete it."""
    import devsim

    from etna.core.device2d import create_sic_2d_device

    name = f"g2d_{uuid.uuid4().hex[:8]}"
    info = create_sic_2d_device(
        device_name=name,
        half_width_um=half_width_um,
        epi_thickness_cm=epi_thickness_cm,
        substrate_thickness_cm=substrate_thickness_cm,
    )
    reg = info["region_name"]
    try:
        x = np.array(devsim.get_node_model_values(device=name, region=reg, name="x"))
        y = np.array(devsim.get_node_model_values(device=name, region=reg, name="y"))
        el = devsim.get_element_node_list(device=name, region=reg)
    finally:
        devsim.delete_device(device=name)
        devsim.delete_mesh(mesh=f"{name}_mesh")
    return grid_from_nodes(x, y, el)
