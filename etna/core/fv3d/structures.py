"""3D tensor-grid structures: column lattices, square-trench microdosimeter cells, test cells.

A ``Structure3D`` is a box [0, Lx] x [0, Ly] x [0, Lz] (cm; z = depth from the top
surface) on graded node lines, a net-doping function of the bulk (the free nodes) and
ohmic electrodes = bodies of nodes held at Dirichlet values (the n+ / p+ columns,
trench walls, implanted planes). Every box side is a natural (zero-flux) boundary: for a
lattice cell with two mirror planes (square, rectangular) the QUARTER cell with
reflecting sides IS the periodic lattice. Diagonal (eighth-cell) symmetry is not
representable on a tensor grid.

Circular columns are staircased: a node is in a column iff its (x, y) lies inside the
circle, so the column surface converges at first order in the local line spacing (the
builders refine the lines over the whole column radius). Lines are built with
``etna.core.fv2d.structures.graded_lines`` (local spacing grows geometrically by
``ratio`` away from each knot): Poisson only needs ~ the Debye length (0.45 um at
8e13 cm^-3) at the column surfaces, NOT the 0.1 nm contact cells of
``etna.core.contact_mesh`` (those resolve a drift-diffusion minority-carrier sink; here
they would only create cell aspect ratios that slow the multigrid smoother).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from etna.core.fv2d.structures import graded_lines

N_D_DEFAULT = (
    8.252732e13  # cm^-3, the n- epi of lab/geometry2d/bench/axisym_validity.py
)
N_ELECTRODE = 1e19  # cm^-3, n+ / p+ column and wall doping


@dataclass
class Electrode:
    """Ohmic electrode body joined to the bias contact ``contact``, majority ``kind`` "n"
    or "p" with doping ``doping`` (sets the built-in potential of its Dirichlet value).

    ``mask(X, Y, Z)`` returns either a LEVEL SET (float array, <= 0 inside the body) or a
    boolean node mask. With a level set the solver places the Dirichlet value on the
    true surface (cut-cell couplings, :class:`etna.core.fv3d.Poisson3D` ``cut_cell``):
    second-order instead of the first-order staircase of a node mask."""

    contact: str
    kind: str
    mask: Callable
    doping: float = N_ELECTRODE


@dataclass
class Structure3D:
    name: str
    xs: np.ndarray
    ys: np.ndarray
    zs: np.ndarray
    doping: Callable  # (X, Y, Z) -> net doping N_D - N_A of the bulk (cm^-3)
    electrodes: list
    T: float = 300.0
    info: dict = field(default_factory=dict)

    @property
    def shape(self):
        return (len(self.zs), len(self.ys), len(self.xs))

    @property
    def N(self):
        return int(np.prod(self.shape))

    def coords(self):
        return np.meshgrid(self.zs, self.ys, self.xs, indexing="ij")


def lines(knots, lo, hi, ratio=1.15, h_max=None):
    return graded_lines(knots, lo, hi, ratio=ratio, h_max=h_max)


def _cap(lateral, depth, partial):
    """Level set of a body = lateral shape intersected with a depth half-space (max of
    the two; only for a partial-depth body: a through body keeps its exact lateral level
    set on every node, also on the bottom plane). Near the tip of a partial body the
    max-composition places the lateral cut-cell surface only approximately (within one
    line spacing of the tip)."""
    return np.maximum(lateral, depth) if partial else lateral


def _disc_knots(c, r, h, lo, hi):
    """Knots covering a disc edge projection [c - r, c + r] with spacing h (the staircase
    surface is resolved at every angle, not only at x = c +- r)."""
    m = max(int(np.ceil(2 * r / (4 * h))), 1)
    pts = np.linspace(c - r, c + r, 2 * m + 1)
    return [(p, h) for p in pts if lo - 1e-15 <= p <= hi + 1e-15]


def column_cell(
    pitch=50e-4,
    r_n=5e-4,
    r_p=5e-4,
    thickness=50e-4,
    depth_n=None,
    depth_p=None,
    p_from="top",
    back=None,
    N_D=N_D_DEFAULT,
    h_fine=0.25e-4,
    h_bulk=2e-4,
    h_z=None,
    ratio=1.15,
    T=300.0,
):
    """Quarter cell of a square lattice of 3D columns (pitch ``pitch``): n+ cathode
    column (radius ``r_n``) at the cell centre (0, 0), p+ anode columns (``r_p``) at the
    corners -> (a, a), a = pitch / 2, n- bulk ``N_D``, cell thickness ``thickness``.

    ``depth_n`` / ``depth_p``: column depths (None = through the whole thickness). The n+
    column grows from the top surface z = 0; the p+ columns from the top (``p_from =
    "top"``, single-sided) or from the bottom (``"bottom"``, double-sided). ``back``:
    None, "p+" (planar anode at z = thickness) or "n+" (planar cathode there).
    Through columns with no back plane make the problem z-invariant (exactly the 2D
    Cartesian cross-section of lab/geometry2d/bench/axisym_validity.py).
    """
    a = pitch / 2
    dn = thickness if depth_n is None else depth_n
    dp = thickness if depth_p is None else depth_p
    if r_n + r_p >= np.sqrt(2) * a:
        raise ValueError("columns overlap: r_n + r_p >= pitch / sqrt(2)")
    kx = [(0.0, h_fine), (a, h_bulk)]
    kx += _disc_knots(0.0, r_n, h_fine, 0, a) + _disc_knots(a, r_p, h_fine, 0, a)
    xs = lines(kx, 0.0, a, ratio, h_bulk)
    hz = h_z or h_bulk
    kz = [(0.0, h_fine), (thickness, h_fine if back else hz)]
    for d in (dn, thickness - dp if p_from == "bottom" else dp):
        if 0 < d < thickness:
            kz.append((d, h_fine))
    zs = lines(kz, 0.0, thickness, ratio, hz)

    # level sets (<= 0 inside)
    def n_col(X, Y, Z):
        return _cap(np.hypot(X, Y) - r_n, Z - dn, dn < thickness)

    def p_col(X, Y, Z):
        dz = (Z - dp) if p_from == "top" else (thickness - dp - Z)
        return _cap(np.hypot(X - a, Y - a) - r_p, dz, dp < thickness)

    el = [Electrode("cathode", "n", n_col), Electrode("anode", "p", p_col)]
    if back in ("p+", "n+"):
        el.append(Electrode("anode" if back == "p+" else "cathode", back[0],
                            lambda X, Y, Z: thickness - Z))  # fmt: skip
    elif back is not None:
        raise ValueError("back must be None, 'p+' or 'n+'")
    info = dict(kind="column_cell", pitch=pitch, r_n=r_n, r_p=r_p, thickness=thickness,
                depth_n=dn, depth_p=dp, p_from=p_from, back=back, N_D=N_D,
                symmetry="quarter cell (mirror planes x = 0, a and y = 0, a)",
                cell_area_factor=4.0)  # fmt: skip
    return Structure3D("column_cell", xs, ys=xs.copy(), zs=zs,
                       doping=lambda X, Y, Z: np.full(X.shape, N_D), electrodes=el,
                       T=T, info=info)  # fmt: skip


def trench_cell(
    pitch=50e-4,
    r_n=5e-4,
    wall=5e-4,
    thickness=10e-4,
    depth_n=None,
    depth_wall=None,
    back=None,
    N_D=N_D_DEFAULT,
    h_fine=0.25e-4,
    h_bulk=1e-4,
    h_z=None,
    ratio=1.15,
    T=300.0,
):
    """Quarter of a square-cell microdosimeter (3D sensitive volume): n+ column (``r_n``)
    at the cell centre, p+ trench WALL of width ``wall`` along the cell boundary (shared
    with the neighbour cells: half of it, wall / 2, inside this cell's quarter at x = a
    and y = a), depths ``depth_n`` / ``depth_wall`` from the top (None = through),
    optional planar ``back`` ("p+" / "n+") at z = thickness. The axisymmetric proxy of
    this cell (an annular trench) is fv2d's (r, z) model."""
    a = pitch / 2
    dn = thickness if depth_n is None else depth_n
    dw = thickness if depth_wall is None else depth_wall
    xw = a - wall / 2
    kx = [(0.0, h_fine), (xw, h_fine), (a, h_fine)] + _disc_knots(0, r_n, h_fine, 0, a)
    xs = lines(kx, 0.0, a, ratio, h_bulk)
    hz = h_z or h_bulk
    kz = [(0.0, h_fine), (thickness, h_fine if back else hz)]
    for d in (dn, dw):
        if 0 < d < thickness:
            kz.append((d, h_fine))
    zs = lines(kz, 0.0, thickness, ratio, hz)
    el = [
        Electrode("cathode", "n",
                  lambda X, Y, Z: _cap(np.hypot(X, Y) - r_n, Z - dn, dn < thickness)),
        Electrode("anode", "p",
                  lambda X, Y, Z: _cap(np.minimum(xw - X, xw - Y), Z - dw, dw < thickness)),
    ]  # fmt: skip
    if back in ("p+", "n+"):
        el.append(Electrode("anode" if back == "p+" else "cathode", back[0],
                            lambda X, Y, Z: thickness - Z))  # fmt: skip
    elif back is not None:
        raise ValueError("back must be None, 'p+' or 'n+'")
    info = dict(kind="trench_cell", pitch=pitch, r_n=r_n, wall=wall, thickness=thickness,
                depth_n=dn, depth_wall=dw, back=back, N_D=N_D,
                symmetry="quarter cell (mirror planes x = 0, a and y = 0, a)",
                cell_area_factor=4.0)  # fmt: skip
    return Structure3D("trench_cell", xs, xs.copy(), zs,
                       lambda X, Y, Z: np.full(X.shape, N_D), el, T, info)  # fmt: skip


def coax_cell(r_in=5e-4, r_out=25e-4, length=4e-4, N_D=N_D_DEFAULT, h=0.5e-4,
              nz=3, T=300.0):  # fmt: skip
    """z-invariant coaxial test cell (quarter): n+ cathode r <= r_in, p+ anode r >= r_out
    in the square [0, r_out + 2h]^2, uniform lines of spacing ``h`` (level-set circles),
    ``nz`` z lines over ``length``. Reference: the (r, z) fv2d solution / analytic coax.
    """
    L = r_out + 2 * h
    xs = np.linspace(0.0, L, int(round(L / h)) + 1)
    zs = np.linspace(0.0, length, nz)
    el = [Electrode("cathode", "n", lambda X, Y, Z: np.hypot(X, Y) - r_in),
          Electrode("anode", "p", lambda X, Y, Z: r_out - np.hypot(X, Y))]  # fmt: skip
    info = dict(kind="coax_cell", r_in=r_in, r_out=r_out, N_D=N_D, h=h)
    return Structure3D("coax_cell", xs, xs.copy(), zs,
                       lambda X, Y, Z: np.full(X.shape, N_D), el, T, info)  # fmt: skip


def slab(thickness=10e-4, N_D=N_D_DEFAULT, zs=None, nxy=3, T=300.0):
    """x,y-invariant planar p+ (z = 0, anode) / n- / n+ (z = thickness, cathode) diode:
    the 1D reference (analytic depletion / overdepleted parabola)."""
    zs = np.linspace(0.0, thickness, 101) if zs is None else np.asarray(zs, float)
    xs = np.linspace(0.0, 1e-4, nxy)
    el = [Electrode("anode", "p", lambda X, Y, Z: Z + 0.0),
          Electrode("cathode", "n", lambda X, Y, Z: thickness - Z)]  # fmt: skip
    info = dict(kind="slab", thickness=thickness, N_D=N_D)
    return Structure3D("slab", xs, xs.copy(), zs,
                       lambda X, Y, Z: np.full(X.shape, N_D), el, T, info)  # fmt: skip


def extrude_2d(xs, ys, zs, doping2d, electrodes2d, T=300.0, name="extruded", info=None):
    """z-extrusion of a 2D Cartesian cross-section: ``doping2d(X, Y)``,
    ``electrodes2d`` = [(contact, kind, mask2d(X, Y), doping)]; parity with fv2d."""
    el = [Electrode(c, k, (lambda m: lambda X, Y, Z: m(X, Y))(m), d)
          for c, k, m, d in electrodes2d]  # fmt: skip
    return Structure3D(name, np.asarray(xs, float), np.asarray(ys, float),
                       np.asarray(zs, float), lambda X, Y, Z: doping2d(X, Y), el, T,
                       dict(info or {}, kind="extruded"))  # fmt: skip
