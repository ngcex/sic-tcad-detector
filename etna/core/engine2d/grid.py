"""devsim's 2D diode mesh as tensor-grid arrays, read from a mesh-only build.

``create_sic_2d_device`` is called with no equations and no solve (a few tens of
ms), its node coordinates, element list, NodeVolume, doping and contact node sets
are read, and the device is deleted. Everything the engine relies on is then
*asserted* from that device, not assumed (lab/solver2d/notes/LOG.md, STEP 0):

* nodes form a full tensor grid xs (nx) x ys (ny);
* NodeVolume equals the rectangular dual cell dx_dual*dy_dual (so the FV scheme
  is the 5-point stencil; devsim's diagonals carry ~1e-33 coupling);
* Donors and Acceptors depend on y only (identical across every column);
* the anode is exactly the full row y = 0 and the cathode the full row y = max.

Engine numbering is k = j*nx + i (x fastest); ``perm[devsim_index] = k``.
"""

from __future__ import annotations

import uuid
from dataclasses import dataclass

import numpy as np

from etna.core.devsim_lock import serialized


class NotLaterallyUniform(ValueError):
    """The devsim device violates an assumption of the 2D engine."""


@dataclass(frozen=True)
class Grid2D:
    xs: np.ndarray  # (nx,) cm
    ys: np.ndarray  # (ny,) cm
    perm: np.ndarray  # (N,) devsim node index -> k
    donors: np.ndarray  # (ny,) cm^-3, as devsim evaluated them
    acceptors: np.ndarray  # (ny,)

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
    def dual(h, n):
        d = np.zeros(n)
        d[:-1] += h / 2
        d[1:] += h / 2
        return d

    @property
    def dxd(self):
        return self.dual(self.hx, self.nx)

    @property
    def dyd(self):
        return self.dual(self.hy, self.ny)

    @property
    def volume(self):
        """(N,) box volumes in k order == devsim NodeVolume."""
        return np.outer(self.dyd, self.dxd).ravel()

    def to_devsim(self, field_k):
        return np.asarray(field_k)[self.perm]


def grid_from_arrays(x, y, volume, donors, acceptors, anode, cathode):
    """Validate devsim node arrays and build the Grid2D (raises NotLaterallyUniform)."""
    xs, ys = np.unique(x), np.unique(y)
    nx, ny = len(xs), len(ys)
    i = np.searchsorted(xs, x)
    j = np.searchsorted(ys, y)
    perm = j * nx + i
    if len(x) != nx * ny or len(np.unique(perm)) != len(x):
        raise NotLaterallyUniform("mesh is not a full tensor-product grid")
    if nx < 2 or ny < 3:
        raise NotLaterallyUniform(f"degenerate grid {nx} x {ny}")
    vol_rect = Grid2D.dual(np.diff(xs), nx)[i] * Grid2D.dual(np.diff(ys), ny)[j]
    if np.max(np.abs(volume - vol_rect) / vol_rect) > 1e-10:
        raise NotLaterallyUniform("NodeVolume is not the rectangular dual cell")
    rows = []
    for name, v in (("Donors", donors), ("Acceptors", acceptors)):
        A = np.empty((ny, nx))
        A[j, i] = v
        if np.any(A.max(axis=1) != A.min(axis=1)):
            raise NotLaterallyUniform(f"{name} vary laterally")
        rows.append(A[:, 0].copy())
    for name, nodes, row in (("anode", anode, 0), ("cathode", cathode, ny - 1)):
        if not np.array_equal(np.sort(perm[nodes]), np.arange(nx) + row * nx):
            raise NotLaterallyUniform(f"{name} is not exactly the full row j={row}")
    return Grid2D(xs=xs, ys=ys, perm=perm, donors=rows[0], acceptors=rows[1])


def _freeze(v):
    """Hashable form of a keyword value (tuples for sequences/arrays); TypeError if none."""
    if isinstance(v, np.ndarray):
        return ("ndarray", v.dtype.str, v.shape, v.tobytes())
    if isinstance(v, (list, tuple)):
        return tuple(_freeze(e) for e in v)
    hash(v)
    return v


_GRID_CACHE: "dict" = {}  # key -> (Grid2D, info), insertion order = age
_GRID_CACHE_MAX = 16


@serialized
def grid_from_devsim(**device_kwargs):
    """Mesh-only ``create_sic_2d_device(**device_kwargs)`` -> validated Grid2D.

    Memoized on the keywords (and the builder function, so a monkeypatched builder is
    called): the mesh-only devsim build costs 0.1-0.3 s per call and is repeated for
    every bias of one device. A hit returns deep copies (the caller may mutate ``info``);
    exceptions are not cached; unhashable keywords bypass the cache.
    """
    import copy

    from etna.core import device2d

    builder = device2d.create_sic_2d_device
    try:
        key = (
            builder,
            tuple(sorted((k, _freeze(v)) for k, v in device_kwargs.items())),
        )
    except TypeError:
        return _grid_from_devsim(builder, device_kwargs)
    hit = _GRID_CACHE.pop(key, None)
    if hit is None:
        hit = _grid_from_devsim(builder, device_kwargs)
    _GRID_CACHE[key] = hit
    while len(_GRID_CACHE) > _GRID_CACHE_MAX:
        _GRID_CACHE.pop(next(iter(_GRID_CACHE)))
    return copy.deepcopy(hit)


def _grid_from_devsim(create_sic_2d_device, device_kwargs):
    import devsim

    name = f"e2d_mesh_{uuid.uuid4().hex[:8]}"
    info = create_sic_2d_device(device_name=name, **device_kwargs)
    reg = info["region_name"]
    try:
        g = lambda n: np.array(
            devsim.get_node_model_values(device=name, region=reg, name=n)
        )
        contacts = {}
        for c in ("anode", "cathode"):
            # a contact node model is evaluated on that contact's nodes only (0 elsewhere)
            devsim.contact_node_model(
                device=name, contact=c, name=f"is_{c}", equation="1"
            )
            contacts[c] = np.flatnonzero(g(f"is_{c}"))
        arrays = dict(
            x=g("x"),
            y=g("y"),
            volume=g("NodeVolume"),
            donors=g("Donors"),
            acceptors=g("Acceptors"),
        )
    finally:
        devsim.delete_device(device=name)
        devsim.delete_mesh(mesh=f"{name}_mesh")
    return (
        grid_from_arrays(
            anode=contacts["anode"], cathode=contacts["cathode"], **arrays
        ),
        info,
    )
