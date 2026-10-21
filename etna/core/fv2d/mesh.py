"""Tensor-grid finite-volume mesh of a laterally NON-uniform 2D device (Cartesian or (r, z)).

A ``FVMesh`` carries everything the discretisation needs as arrays on the node grid
(ny, nx) (row j = depth y_j, column i = lateral x_i, or radius r_i):

* ``vol``      node (box) volume: cm^2 per cm of device length (Cartesian half-device) or
               cm^3 (cylindrical, includes the 2 pi);
* ``gh, gv``   edge geometric factor EdgeCouple / EdgeLength of the lateral edges
               (ny, nx-1) and of the depth edges (ny-1, nx) (cm per cm, or cm);
* ``donors, acceptors`` (ionized acceptors, as devsim's ``Acceptors``), ``net``;
* ``mun_h/mun_v/mup_h/mup_v`` edge mobilities (cm^2/Vs) on the same edge arrays;
* ``contacts`` name -> boolean (ny, nx) mask of the contact's (Dirichlet) nodes;
* ``active`` boolean (ny, nx) mask of the nodes that carry equations (inactive nodes,
  e.g. outside an L-shaped region, are frozen and decoupled: their edges are zeroed).

Two producers: :mod:`etna.core.fv2d.devsim_import` (devsim's own NodeVolume /
CylindricalNodeVolume and EdgeCouple, verbatim: parity mode) and
:mod:`etna.core.fv2d.structures` (exact box integrals on a native tensor mesh: the
Cartesian dual cell, or 2 pi int r dr over the dual cell in (r, z), see
:func:`box_geometry`).

Solver ordering: the shorter grid dimension is numbered fastest so the half-bandwidth of
the 5-point matrices is min(nx, ny) (``SolverOrder``).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class FVMesh:
    xs: np.ndarray
    ys: np.ndarray
    vol: np.ndarray
    gh: np.ndarray
    gv: np.ndarray
    donors: np.ndarray
    acceptors: np.ndarray
    net: np.ndarray
    mun_h: np.ndarray
    mun_v: np.ndarray
    mup_h: np.ndarray
    mup_v: np.ndarray
    contacts: dict
    coordinate_system: str = "cartesian"
    active: np.ndarray | None = None
    info: dict = field(default_factory=dict)
    # incomplete ionization is not supported here (static doping); see USAGE

    def __post_init__(self):
        ny, nx = len(self.ys), len(self.xs)
        for name, shape in (
            ("vol", (ny, nx)),
            ("gh", (ny, nx - 1)),
            ("gv", (ny - 1, nx)),
            ("donors", (ny, nx)),
            ("acceptors", (ny, nx)),
            ("net", (ny, nx)),
            ("mun_h", (ny, nx - 1)),
            ("mup_h", (ny, nx - 1)),
            ("mun_v", (ny - 1, nx)),
            ("mup_v", (ny - 1, nx)),
        ):
            a = np.ascontiguousarray(getattr(self, name), dtype=float)
            if a.shape != shape:
                raise ValueError(f"FVMesh.{name} has shape {a.shape}, expected {shape}")
            setattr(self, name, a)
        if self.coordinate_system not in ("cartesian", "cylindrical"):
            raise ValueError(f"coordinate_system {self.coordinate_system!r}")
        if not self.contacts:
            raise ValueError("FVMesh needs at least one contact")
        for c, m in self.contacts.items():
            m = np.asarray(m, dtype=bool)
            if m.shape != (ny, nx) or not m.any():
                raise ValueError(f"contact {c!r}: mask shape {m.shape} / empty")
            self.contacts[c] = m
        if self.active is None:
            self.active = np.ones((ny, nx), dtype=bool)
        self.active = np.asarray(self.active, dtype=bool)

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
    def X(self):
        return np.broadcast_to(self.xs[None, :], (self.ny, self.nx))

    @property
    def Y(self):
        return np.broadcast_to(self.ys[:, None], (self.ny, self.nx))

    def contact_any(self):
        m = np.zeros((self.ny, self.nx), dtype=bool)
        for c in self.contacts.values():
            m |= c
        return m


def dual_lengths(h, n):
    d = np.zeros(n)
    d[:-1] += h / 2
    d[1:] += h / 2
    return d


def box_geometry(xs, ys, coordinate_system="cartesian"):
    """Exact box (dual-cell) geometry of a tensor mesh: (vol, gh, gv).

    Cartesian: vol = dx_dual dy_dual, lateral edge factor dy_dual / hx, depth edge
    factor dx_dual / hy (devsim's 2D values on a rectangle mesh, diagonals 0).
    Cylindrical (x = r): vol = 2 pi int r dr over the dual interval * dy_dual, lateral
    (radial) edge factor 2 pi r_mid dy_dual / hx, depth (axial) edge factor
    2 pi int r dr over the dual interval / hy. devsim's CylindricalEdgeCouple equals these
    edge factors to round-off; its CylindricalNodeVolume differs by up to 0.2 % (5/3 on
    the axis) because it is assembled per triangle (checked 2026-10-21 on the 3D-electrode
    mesh): parity runs import devsim's volumes instead (devsim_import).
    """
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    hx, hy = np.diff(xs), np.diff(ys)
    nx, ny = len(xs), len(ys)
    dyd = dual_lengths(hy, ny)
    if coordinate_system == "cartesian":
        dxd = dual_lengths(hx, nx)
        vol = np.outer(dyd, dxd)
        gh = np.outer(dyd, 1.0 / hx)
        gv = np.outer(1.0 / hy, dxd)
    elif coordinate_system == "cylindrical":
        if xs[0] < 0:
            raise ValueError("cylindrical mesh needs r >= 0")
        mid = 0.5 * (xs[:-1] + xs[1:])
        lo = np.r_[xs[0], mid]
        hi = np.r_[mid, xs[-1]]
        ring = np.pi * (hi**2 - lo**2)  # 2 pi int r dr
        vol = np.outer(dyd, ring)
        gh = np.outer(dyd, 2 * np.pi * mid / hx)
        gv = np.outer(1.0 / hy, ring)
    else:
        raise ValueError(f"coordinate_system {coordinate_system!r}")
    return vol, gh, gv


@dataclass(frozen=True)
class SolverOrder:
    """Flat solver numbering with the SHORTER dimension fastest.

    ``transpose`` False: k = j * nx + i (x fastest, m = nx). True: k = i * ny + j
    (y fastest, m = ny). ``flat(A)`` maps a (ny, nx) node array to the (N,) solver
    vector, ``grid(v)`` back. The fast-direction edges (``ef``) are the x edges when
    not transposed, else the y edges.
    """

    nx: int
    ny: int
    transpose: bool

    @staticmethod
    def for_mesh(mesh):
        return SolverOrder(mesh.nx, mesh.ny, mesh.nx > mesh.ny)

    @property
    def m(self):
        return self.ny if self.transpose else self.nx

    @property
    def ns(self):
        return self.nx if self.transpose else self.ny

    def flat(self, A):
        A = np.asarray(A)
        return np.ascontiguousarray((A.T if self.transpose else A).ravel())

    def grid(self, v):
        v = np.asarray(v)
        if self.transpose:
            return v.reshape(self.nx, self.ny).T.copy()
        return v.reshape(self.ny, self.nx).copy()

    def edges(self, eh, ev):
        """(fast, slow) edge arrays in solver layout: fast (ns, m-1), slow (ns-1, m)."""
        if self.transpose:
            return np.ascontiguousarray(ev.T), np.ascontiguousarray(eh.T)
        return np.ascontiguousarray(eh), np.ascontiguousarray(ev)
