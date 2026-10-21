"""Native (devsim-free) tensor meshes of the notebook-19 structures for the fv2d engine.

Each ``*_structure`` function takes the keywords of its devsim twin in
etna.core.alternative_structures (same defaults: the ALT_STACK pin, graded epi) and
returns a :class:`Structure`: the doping as numpy functions of (x, y) that reproduce the
devsim node models (tests/test_fv2d.py checks them against devsim's evaluated Donors /
Acceptors / NetDoping on devsim's own mesh), the contact segments, and the "knot" lines
(position, spacing) of the tensor mesh. :func:`build_mesh` turns a Structure and node
lines into an FVMesh with EXACT box geometry (mesh.box_geometry); the doping is always
re-evaluated on the nodes from the functions, never interpolated (refine.py inserts lines
and rebuilds).

devsim's ``step(x)`` is 1 for x >= 0 (checked on the builders' meshes: implant-edge
nodes carry the implant), and the abrupt implant edges are reproduced with it; a
``lateral_straggle_cm`` > 0 gives the builders' logistic roll-off.
Only ``doping_profile="graded"`` (the builders' default) is implemented here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from etna.core.alternative_structures import (
    ALT_EPI_PS,
    ALT_JUNCTION_PS,
    ALT_LATERAL_EDGE_PS,
    ALT_STACK,
)
from etna.core.device import DEFAULT_N_A, acceptor_total_scale, edge_mobility
from etna.core.incomplete_ionization import ionized_acceptor_concentration
from etna.core.sic_material import (
    SiC4H_Parameters,
    intrinsic_concentration,
    mobility_caughey_thomas_T,
)

from .mesh import FVMesh, box_geometry


@dataclass
class Structure:
    """Geometry + doping of a 2D device, mesh-independent.

    ``doping(X, Y) -> (donors, acceptors_ionized, acceptors_total)`` on node arrays;
    ``contacts``: name -> (side, lo, hi) with side in {"top", "bottom", "left", "right"}
    and [lo, hi] the covered interval of the other coordinate (inclusive), or a callable
    (X, Y) -> boolean node mask (e.g. a whole electrode column as a Dirichlet body);
    ``x_knots`` / ``y_knots``: [(position, spacing)] lines of the initial mesh.
    """

    name: str
    width: float
    depth: float
    doping: Callable
    contacts: dict
    x_knots: list
    y_knots: list
    coordinate_system: str = "cartesian"
    T: float = 300.0
    info: dict = field(default_factory=dict)


def step(x):
    """devsim ``step``: 1 for x >= 0, else 0."""
    return (np.asarray(x) >= 0).astype(float)


def _edge(s, straggle):
    """alternative_structures._edge_expr: 1 where s < 0, 0 where s > 0."""
    if straggle and straggle > 0:
        return 1.0 / (1.0 + np.exp(np.minimum(s / straggle, 700.0)))
    return step(-s)


def _graded(y, yj, Nj, Nb, L):
    """graded_base_expr (exponential epi, clamped at its junction value above yj)."""
    return Nb + (Nj - Nb) * np.exp(-np.maximum(y - yj, 0.0) / L)


def _stack(N_D_junction, N_D_bulk, L_transition):
    return (
        ALT_STACK["N_D_junction"] if N_D_junction is None else N_D_junction,
        ALT_STACK["N_D_bulk"] if N_D_bulk is None else N_D_bulk,
        ALT_STACK["L_transition"] if L_transition is None else L_transition,
    )


def _check_graded(doping_profile):
    if doping_profile != "graded":
        raise NotImplementedError(
            "fv2d native structures implement doping_profile='graded' only"
        )


def _base_info(T, epi, sub, half_width, name, coordinate_system, **extra):
    prm = SiC4H_Parameters()
    n_i = intrinsic_concentration(T, prm)[0]
    info = dict(
        T=T,
        n_i=n_i,
        params=prm,
        structure_type=name,
        coordinate_system=coordinate_system,
        substrate_thickness_cm=sub,
        junction_pos=sub,
        epi_thickness_cm=epi,
        total_length=sub + epi,
        half_width_cm=half_width,
        nplus_thickness_cm=None,
        buffer_thickness_cm=None,
        stack="ALT_STACK",
        dimension=2,
    )
    info.update(extra)
    return info


def _depth_knots(junction, total, extra=()):
    """Depth lines of alternative_structures._add_depth_mesh_lines (region part)."""
    k = [(0.0, 1e-5)]
    if junction - 5e-6 > 1e-8:
        k.append((junction - 5e-6, max(1e-6, ALT_JUNCTION_PS)))
    k.append((junction, ALT_JUNCTION_PS))
    for y in (junction + 2e-4, junction + 5e-4):
        if y < total - 1e-6:
            k.append((y, ALT_EPI_PS))
    k.append((total, 1e-5))
    k.extend(extra)
    return k


# --------------------------------------------------------------------- builders
def guard_ring_structure(
    sv_half_width_um=50.0,
    guard_ring_width_um=5.0,
    guard_ring_gap_um=3.0,
    guard_ring_depth_um=1.0,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=ALT_STACK["substrate_thickness_cm"],
    N_A=DEFAULT_N_A,
    N_A_guard=5e18,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
    with_guard=True,
    lateral_straggle_cm=0.0,
):
    """Twin of create_guard_ring_device (same keywords and doping)."""
    _check_graded(doping_profile)
    Nj, Nb, L = _stack(N_D_junction, N_D_bulk, L_transition)
    sv, w, gap, gd = (
        v * 1e-4
        for v in (
            sv_half_width_um,
            guard_ring_width_um,
            guard_ring_gap_um,
            guard_ring_depth_um,
        )
    )
    yj = substrate_thickness_cm
    total = yj + epi_thickness_cm
    gi, go = sv + gap, sv + gap + w
    NAi = ionized_acceptor_concentration(N_A, T)
    NGi = ionized_acceptor_concentration(N_A_guard, T) if with_guard else 0.0
    sA = acceptor_total_scale(N_A, NAi)
    sG = acceptor_total_scale(N_A_guard, NGi) if with_guard else 0.0

    def doping(X, Y):
        don = _graded(Y, yj, Nj, Nb, L) + 0.0 * X
        a_an = NAi * step(yj - Y) * _edge(X - sv, lateral_straggle_cm)
        a_gr = NGi * step(gd - Y) * _edge(gi - X, lateral_straggle_cm)
        return don, a_an + a_gr, sA * a_an + sG * a_gr

    edge_ps = max(0.5 * lateral_straggle_cm, ALT_LATERAL_EDGE_PS)
    xk = [(0.0, 5e-5)]
    if sv * 0.5 > 1e-6:
        xk.append((sv * 0.5, 2e-4))
    xk += [(sv, edge_ps), (sv + 0.5 * gap, 2e-5), (gi, edge_ps)]
    if w > 2e-5:
        xk.append((gi + 0.5 * w, 5e-5))
    xk.append((go, 1e-5))
    extra = [(gd, 1e-6)] if gd > 1e-8 and abs(gd - yj) > 1e-8 else []
    contacts = {"anode": ("top", 0.0, sv), "cathode": ("bottom", 0.0, go)}
    if with_guard:
        contacts["guard_ring_anode"] = ("top", gi, go)
    info = _base_info(
        T, epi_thickness_cm, substrate_thickness_cm, sv, "guard_ring", "cartesian",
        sv_half_width_cm=sv, gap_cm=gap, guard_inner_cm=gi, guard_outer_cm=go,
        guard_depth_cm=gd, total_half_width_cm=go, N_A=N_A, N_A_ionized=NAi,
        N_A_guard=N_A_guard if with_guard else 0.0,
        guard_ring_contact="guard_ring_anode" if with_guard else None,
        contacts=tuple(contacts), N_D_junction=Nj, N_D_bulk=Nb, L_transition=L,
    )  # fmt: skip
    return Structure(
        "guard_ring",
        go,
        total,
        doping,
        contacts,
        xk,
        _depth_knots(yj, total, extra),
        "cartesian",
        T,
        info,
    )


def mesa_structure(
    pillar_half_width_um=50.0,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=ALT_STACK["substrate_thickness_cm"],
    N_A=DEFAULT_N_A,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
    **_trench,
):
    """Twin of create_mesa_device: the solved region is the pillar (the trench carries no
    equations: its sidewall is a zero-flux boundary), i.e. the planar half-device."""
    _check_graded(doping_profile)
    Nj, Nb, L = _stack(N_D_junction, N_D_bulk, L_transition)
    hw = pillar_half_width_um * 1e-4
    yj = substrate_thickness_cm
    total = yj + epi_thickness_cm
    NAi = ionized_acceptor_concentration(N_A, T)
    sA = acceptor_total_scale(N_A, NAi)

    def doping(X, Y):
        don = _graded(Y, yj, Nj, Nb, L) * step(Y - yj) + 0.0 * X
        acc = NAi * step(yj - Y) + 0.0 * X
        return don, acc, sA * acc

    xk = [(0.0, 5e-5)]
    if hw > 4e-4:
        xk += [(0.25 * hw, 2e-4), (0.5 * hw, 3e-4), (0.75 * hw, 4e-4)]
    xk.append((hw, 2e-4))
    contacts = {"anode": ("top", 0.0, hw), "cathode": ("bottom", 0.0, hw)}
    info = _base_info(
        T, epi_thickness_cm, substrate_thickness_cm, hw, "mesa", "cartesian",
        N_A=N_A, N_A_ionized=NAi, contacts=tuple(contacts),
        N_D_junction=Nj, N_D_bulk=Nb, L_transition=L,
    )  # fmt: skip
    return Structure(
        "mesa",
        hw,
        total,
        doping,
        contacts,
        xk,
        _depth_knots(yj, total),
        "cartesian",
        T,
        info,
    )


def electrode3d_structure(
    outer_radius_um=50.0,
    column_radius_um=5.0,
    top_gap_um=2.0,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=ALT_STACK["substrate_thickness_cm"],
    N_A=DEFAULT_N_A,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
    lateral_straggle_cm=0.0,
):
    """Twin of create_3d_electrode_device (axisymmetric (r, z))."""
    _check_graded(doping_profile)
    Nj, Nb, L = _stack(N_D_junction, N_D_bulk, L_transition)
    ro, rc, tg = outer_radius_um * 1e-4, column_radius_um * 1e-4, top_gap_um * 1e-4
    si = rc + tg
    wall = 1e-4
    wi = ro - wall
    if not (rc < si < wi < ro):
        raise ValueError("3D electrode: column + top gap must fit inside the p+ wall")
    yj = substrate_thickness_cm
    total = yj + epi_thickness_cm
    NAi = ionized_acceptor_concentration(N_A, T)
    sA = acceptor_total_scale(N_A, NAi)

    def doping(X, Y):
        don = np.maximum(1e19 * step(rc - X), _graded(Y, yj, Nj, Nb, L))
        a_sheet = NAi * step(yj - Y) * _edge(si - X, lateral_straggle_cm)
        a_wall = NAi * _edge(wi - X, lateral_straggle_cm) + 0.0 * Y
        acc = np.maximum(a_sheet, a_wall)
        return don, acc, sA * acc

    edge_ps = max(0.5 * lateral_straggle_cm, ALT_LATERAL_EDGE_PS)
    xk = [
        (0.0, 5e-5),
        (rc, 2e-5),
        (si, edge_ps),
        (0.5 * (si + wi), 2e-4),
        (wi, edge_ps),
        (ro, 2e-5),
    ]
    contacts = {"cathode": ("top", 0.0, rc), "anode": ("right", 0.0, total)}
    info = _base_info(
        T, epi_thickness_cm, substrate_thickness_cm, ro, "3d_electrode", "cylindrical",
        outer_radius_cm=ro, column_radius_cm=rc, top_gap_cm=tg, sheet_inner_cm=si,
        wall_thickness_cm=wall, N_A=N_A, N_A_ionized=NAi, contacts=tuple(contacts),
        N_D_junction=Nj, N_D_bulk=Nb, L_transition=L,
    )  # fmt: skip
    return Structure(
        "3d_electrode",
        ro,
        total,
        doping,
        contacts,
        xk,
        _depth_knots(yj, total),
        "cylindrical",
        T,
        info,
    )


# ------------------------------------------------------------------- meshing
def graded_lines(knots, extent_lo, extent_hi, ratio=1.15, h_max=None):
    """Node coordinates through every knot (pos, spacing): local spacing
    h(x) = min_k (ps_k + (ratio - 1) |x - pos_k|) (capped at ``h_max``), each knot
    interval meshed by marching h and rescaled to end exactly on the next knot."""
    kn = sorted(
        {
            (float(p), float(h))
            for p, h in knots
            if extent_lo - 1e-15 <= p <= extent_hi + 1e-15
        }
    )
    pos = np.array([p for p, _ in kn])
    ps = np.array([h for _, h in kn])
    # merge duplicates (same position): finest spacing
    up, idx = np.unique(pos, return_inverse=True)
    hs = np.full(len(up), np.inf)
    np.minimum.at(hs, idx, ps)
    pos, ps = up, hs
    if pos[0] > extent_lo:
        pos, ps = np.r_[extent_lo, pos], np.r_[ps[0], ps]
    if pos[-1] < extent_hi:
        pos, ps = np.r_[pos, extent_hi], np.r_[ps, ps[-1]]

    def h(x):
        v = np.min(ps + (ratio - 1.0) * np.abs(x - pos))
        return min(v, h_max) if h_max else v

    out = [pos[0]]
    for a, b in zip(pos[:-1], pos[1:]):
        xs = [a]
        x = a
        while True:
            # symmetric step: average of the spacing here and one step ahead
            dx = h(x)
            dx = 0.5 * (dx + h(min(x + dx, b)))
            if x + dx >= b - 0.3 * dx:
                break
            x += dx
            xs.append(x)
        xs.append(b)
        xs = np.array(xs)
        xs = a + (xs - a) * (b - a) / (xs[-1] - a)
        out.extend(xs[1:])
    return np.array(out)


def initial_lines(structure, ratio=1.15):
    s = structure
    xs = graded_lines(s.x_knots, 0.0, s.width, ratio)
    ys = graded_lines(s.y_knots, 0.0, s.depth, ratio)
    return xs, ys


def _contact_mask(side, lo, hi, xs, ys, tol=1e-12):
    ny, nx = len(ys), len(xs)
    m = np.zeros((ny, nx), dtype=bool)
    if side in ("top", "bottom"):
        j = 0 if side == "top" else ny - 1
        m[j, (xs >= lo - tol) & (xs <= hi + tol)] = True
    else:
        i = 0 if side == "left" else nx - 1
        m[(ys >= lo - tol) & (ys <= hi + tol), i] = True
    return m


def build_mesh(structure, xs=None, ys=None, params=None):
    """FVMesh of ``structure`` on the node lines xs, ys (default: initial_lines)."""
    s = structure
    if xs is None or ys is None:
        xs0, ys0 = initial_lines(s)
        xs = xs0 if xs is None else xs
        ys = ys0 if ys is None else ys
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    prm = params or s.info.get("params") or SiC4H_Parameters()
    X, Y = np.meshgrid(xs, ys)
    don, acc, acc_tot = (np.asarray(a, float) * np.ones_like(X) for a in s.doping(X, Y))
    vol, gh, gv = box_geometry(xs, ys, s.coordinate_system)
    Ntot = don + acc_tot
    mob = {}
    for c, tag in (("electron", "n"), ("hole", "p")):
        mu = mobility_caughey_thomas_T(Ntot, s.T, c, prm)
        mob[f"mu{tag}_h"] = edge_mobility(mu)  # along x (last axis)
        mob[f"mu{tag}_v"] = edge_mobility(mu.T).T  # along y
    contacts = {
        c: (
            np.asarray(seg(X, Y), bool)
            if callable(seg)
            else _contact_mask(*seg, xs, ys)
        )
        for c, seg in s.contacts.items()
    }
    return FVMesh(
        xs=xs, ys=ys, vol=vol, gh=gh, gv=gv, donors=don, acceptors=acc, net=don - acc,
        contacts=contacts, coordinate_system=s.coordinate_system,
        info=dict(s.info), **mob,
    )  # fmt: skip
