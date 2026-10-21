"""fv3d: 3D Poisson on tensor grids (issue #4). Small grids, each test well under 10 s."""

import math

import numpy as np
import pytest
import scipy.sparse.linalg as spla

from etna.core.device2d import EPS_0, Q
from etna.core.fv3d import (
    Electrode,
    Poisson3D,
    Structure3D,
    coax_cell,
    column_cell,
    extrude_2d,
    slab,
)
from etna.core.fv3d.mg import Multigrid, pcg
from etna.core.sic_material import SiC4H_Parameters

ND = 8.252732e13


def small_cell(**kw):
    g = dict(pitch=20e-4, r_n=2e-4, r_p=2e-4, thickness=4e-4, h_fine=0.5e-4,
             h_bulk=1e-4)  # fmt: skip
    g.update(kw)
    return column_cell(**g)


def test_operator_and_vcycle_symmetric_pcg_matches_direct():
    n = 13
    xs = np.linspace(0, 1, n)
    zs = np.r_[0, np.cumsum(1.2 ** np.arange(n - 1))]
    zs /= zs[-1]
    free = np.ones((n, n, n), np.uint8)
    free[-1] = 0
    free[: n // 2, :3, :3] = 0
    mg = Multigrid(xs, xs, zs, free, coarse_n=200)
    assert len(mg.levels) >= 2
    rng = np.random.default_rng(1)
    c = rng.random(free.shape)
    mg.set_c(c)
    y, z = (rng.standard_normal(free.shape) * (free != 0) for _ in range(2))
    L = mg.fine
    a, b = np.vdot(y, L.apply(z)), np.vdot(z, L.apply(y))
    assert abs(a - b) <= 1e-13 * abs(a)
    a, b = np.vdot(y, mg.vcycle(z)), np.vdot(z, mg.vcycle(y))
    assert abs(a - b) <= 1e-12 * abs(a)
    bvec = rng.standard_normal(free.shape) * (free != 0)
    x, it, rr = pcg(mg, bvec, rtol=1e-12)
    xr = spla.spsolve(L.sparse().tocsc(), bvec.ravel()).reshape(bvec.shape)
    assert rr < 1e-12 and it < 40
    assert np.abs(x - xr).max() <= 1e-10 * np.abs(xr).max()


def test_parity_with_fv2d_on_extruded_cross_section():
    """z-extrusion of a 2D Cartesian column cell = fv2d's Solver2D on the same lines."""
    from etna.core.fv2d import Solver2D
    from etna.core.fv2d.structures import Structure, _base_info, build_mesh, graded_lines

    a, r, NP = 10e-4, 2e-4, 1e19

    def doping(X, Y):
        acc = np.where(np.hypot(X - a, Y - a) <= r, NP, 0.0)
        return np.where(np.hypot(X, Y) <= r, NP, ND), acc, acc

    info = _base_info(300.0, 2 * a, 0.0, a, "sq", "cartesian")
    cont = {"cathode": lambda X, Y: np.hypot(X, Y) <= r,
            "anode": lambda X, Y: np.hypot(X - a, Y - a) <= r}  # fmt: skip
    st = Structure("sq", a, a, doping, cont, [], [], "cartesian", 300.0, info)
    xs = graded_lines([(0.0, 2e-5), (r, 2e-5), (a - r, 2e-5), (a, 1e-4)], 0.0, a, 1.3)
    mesh = build_mesh(st, xs, xs)
    s2 = Solver2D(mesh).equilibrium()
    el = [("cathode", "n", cont["cathode"], NP), ("anode", "p", cont["anode"], NP - ND)]
    p = Poisson3D(
        extrude_2d(xs, xs, [0.0, 1e-4, 2e-4], lambda X, Y: np.full(X.shape, ND), el)
    )
    p.equilibrium()
    for V in (5.0, 40.0):
        s2.set_bias({"cathode": V})
        r3 = p.solve({"cathode": V})
        assert np.abs(r3.psi - s2.potential()[None]).max() < 1e-10


def test_overdepleted_slab_is_the_exact_parabola_on_graded_lines():
    T = 10e-4
    zs = np.unique(np.r_[np.linspace(0, T, 11), np.geomspace(1e-7, T, 15)])
    p = Poisson3D(slab(T, ND, zs))
    V = 150.0
    psi = p.solve({"cathode": V}, carriers=False).psi[:, 1, 1]
    eps = SiC4H_Parameters().eps_r * EPS_0
    a0 = -p.Vt * math.log(1e19 / p.ni)
    b0 = V + p.Vt * math.log(1e19 / p.ni)
    exact = a0 + (b0 - a0) * zs / T + Q * ND / (2 * eps) * zs * (T - zs)
    assert np.abs(psi - exact).max() < 1e-9


def test_cut_cell_coax_converges_at_second_order():
    rin, rout = 5e-4, 15e-4
    errs = []
    for h in (1e-4, 0.5e-4):
        p = Poisson3D(coax_cell(rin, rout, N_D=ND, h=h, nz=2))
        r = p.solve({"cathode": 50.0}, carriers=False)
        eps = SiC4H_Parameters().eps_r * EPS_0
        pin = 50.0 + p.Vt * math.log(1e19 / p.ni)
        pout = -p.Vt * math.log(1e19 / p.ni)
        k = Q * ND / (4 * eps)
        A = (pout - pin + k * (rout**2 - rin**2)) / math.log(rout / rin)
        B = pin + k * rin**2 - A * math.log(rin)
        Z, Y, X = np.meshgrid(p.zs, p.ys, p.xs, indexing="ij")
        R = np.hypot(X, Y)
        m = (p.free != 0) & (R > rin + 1.5e-4) & (R < rout - 1.5e-4)
        errs.append(np.abs(r.psi[m] - (-k * R[m] ** 2 + A * np.log(R[m]) + B)).max())
    assert errs[1] < 0.06  # V, of 50 V
    assert errs[0] / errs[1] > 3.0  # second order (staircase: ~1.4)


def test_quarter_cell_equals_the_mirrored_lattice():
    """Mirror planes = periodic lattice: the quarter cell equals the half cell built
    explicitly with the neighbour's n+ column."""
    q = small_cell()
    a = q.info["pitch"] / 2
    xs2 = np.r_[q.xs, 2 * a - q.xs[::-1][1:]]
    r = q.info["r_n"]
    el = [Electrode("cathode", "n", lambda X, Y, Z: np.minimum(np.hypot(X, Y), np.hypot(X - 2 * a, Y)) - r),
          Electrode("anode", "p", lambda X, Y, Z: np.hypot(X - a, Y - a) - r)]  # fmt: skip
    h = Structure3D("half", xs2, q.ys, q.zs, q.doping, el, 300.0, {})
    pq, ph = Poisson3D(q), Poisson3D(h)
    pq.equilibrium()
    ph.equilibrium()
    sq, sh = pq.solve({"cathode": 30.0}), ph.solve({"cathode": 30.0})
    nx = len(q.xs)
    assert np.abs(sh.psi[:, :, :nx] - sq.psi).max() < 1e-7
    assert np.abs(sh.psi[:, :, nx - 1 :] - sq.psi[:, :, ::-1]).max() < 1e-7


def test_null_point_weighting_sum_and_direct_jump():
    p = Poisson3D(small_cell())
    p.equilibrium()
    s = p.solve({"cathode": 60.0})  # direct jump from 0 V: line search + bisection
    assert p.depleted_fraction(s) > 0.99
    # the |E| minimum sits on the symmetric point between two anode columns
    sad = p.saddles(s, 1, exclude_contact_layers=2)
    a = small_cell().info["pitch"] / 2
    assert sad and {round(sad[0]["x"] / a, 6), round(sad[0]["y"] / a, 6)} == {0.0, 1.0}
    wc, _ = p.weighting_potential("cathode")
    wa, _ = p.weighting_potential("anode")
    assert np.abs(wc + wa - 1.0).max() < 1e-8
    assert wc.min() > -1e-10 and wc.max() < 1 + 1e-10


def test_run_field_3d_facade():
    import etna

    r = etna.run_field_3d(None, "column", V_bias=(0.5, 2.0, 20.0), pitch=20, r_n=2,
                          r_p=2, thickness=4, h_fine=0.5, h_bulk=1.0,
                          weighting=("cathode",))  # fmt: skip
    assert isinstance(r, etna.SimResult) and r.sim_type == "field3d"
    assert np.all(np.diff(r.y) > 0) and r.y[-1] > 0.99
    md = r.metadata
    assert md["V_depletion"] is not None and 2.0 < md["V_depletion"] < 20.0
    assert md["psi"].shape == md["shape"] and "cathode" in md["weighting"]
    assert md["structure_info"]["pitch"] == pytest.approx(20e-4)
    with pytest.raises(ValueError):
        etna.run_field_3d(None, "hexagon")


def test_ramo_signal_pair_collects_one_charge():
    p = Poisson3D(small_cell())
    p.equilibrium()
    s = p.solve({"cathode": 60.0})
    w, _ = p.weighting_potential("cathode")
    from etna.core.fv3d import ramo_signal

    start = (5e-4, 3e-4, 2e-4)
    e = ramo_signal(p, s, w, start, "electron")
    h = ramo_signal(p, s, w, start, "hole")
    assert e["stopped"] == "electrode" and h["stopped"] == "electrode"
    assert e["Q"][-1] + h["Q"][-1] == pytest.approx(-1.0, abs=0.05)


def test_null_point_field_vanishes_under_refinement():
    """Without the face-normal zeroing (raw one-sided np.gradient) the field at the
    point (a, 0) between two anode columns falls with h: a null point, not a saddle of
    finite field (issue #4's 90-360 V/cm)."""
    raw = []
    for h in (1.0, 0.5):
        p = Poisson3D(small_cell(h_fine=0.5e-4 * h, h_bulk=1e-4 * h))
        p.equilibrium()
        s = p.solve({"cathode": 40.0})
        Ez, Ey, Ex = np.gradient(s.psi, p.zs, p.ys, p.xs)
        raw.append(np.hypot(Ex, Ey)[1, 0, -1])
    assert raw[1] < 0.7 * raw[0]


def test_facade_trench_partial_columns_linear_and_secant():
    import etna

    g = dict(pitch=20, thickness=4, h_fine=0.5, h_bulk=1.0, N_D=ND)
    t = etna.run_field_3d(None, "trench", V_bias=(1.0, 10.0), r_n=2, wall=2,
                          find_depletion_voltage=True, weighting=(), **g)  # fmt: skip
    assert t.metadata["structure_info"]["kind"] == "trench_cell"
    assert 1.0 < t.metadata["V_depletion"] < 10.0 and t.y[-1] > 0.99
    c = etna.run_field_3d(None, "column", V_bias=(50.0,), r_n=2, r_p=2, depth_n=2,
                          depth_p=2, back="p+", carriers=False, **g)  # fmt: skip
    assert np.isnan(c.y).all() and c.metadata["max_E"][0] > 0
    b = etna.run_field_3d(None, "column", V_bias=(5.0, 20.0), r_n=2, r_p=2, depth_n=3,
                          depth_p=3, p_from="bottom", weighting=(), **g)  # fmt: skip
    assert b.y[1] > b.y[0]
