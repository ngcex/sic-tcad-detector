"""fv3d validation: parity with fv2d (z-extrusion), 1D analytic slab, coax vs analytic and
vs the axisymmetric fv2d (r, z) solution, manufactured Laplace (order of accuracy).

Usage: uv run python lab/fv3d/bench/validate.py
"""

import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "geometry2d" / "bench"))

from etna.core.device2d import EPS_0, Q  # noqa: E402
from etna.core.fv3d import Poisson3D, coax_cell, extrude_2d, slab  # noqa: E402
from etna.core.fv3d.mg import Multigrid, apply_op, pcg  # noqa: E402
from etna.core.sic_material import SiC4H_Parameters  # noqa: E402

EPS = SiC4H_Parameters().eps_r * EPS_0


def parity_fv2d():
    import axisym_validity as av

    from etna.core.fv2d import Solver2D
    from etna.core.fv2d.structures import build_mesh, graded_lines

    st = av.square_cell()
    mesh = build_mesh(
        st, graded_lines(st.x_knots, 0, st.width), graded_lines(st.y_knots, 0, st.depth)
    )
    s2 = Solver2D(mesh).equilibrium()
    a = av.P / 2
    el = [("cathode", "n", lambda X, Y: np.hypot(X, Y) <= av.RC, av.N_PLUS),
          ("anode", "p", lambda X, Y: np.hypot(X - a, Y - a) <= av.RP, av.N_PLUS - av.ND)]  # fmt: skip
    p3 = Poisson3D(extrude_2d(mesh.xs, mesh.ys, np.linspace(0, 4e-4, 3),
                              lambda X, Y: np.full(X.shape, av.ND), el))  # fmt: skip
    p3.equilibrium()
    print(
        f"[parity] square column cell (axisym_validity), {mesh.nx}x{mesh.ny}x3 = {p3.N} nodes"
    )
    for V in (10, 50, 75, 100):
        s2.set_bias({"cathode": V})
        r = p3.solve({"cathode": V})
        d = np.abs(r.psi - s2.potential()[None]).max()
        print(f"  {V:4d} V: max|psi_3D - psi_fv2d| {d:.1e} V, depleted {p3.depleted_fraction(r):.5f}, "
              f"Newton {r.newton_iterations}, CG {r.cg_iterations}")  # fmt: skip


def slab_1d():
    T, N = 10e-4, 8.252732e13
    print("[slab] p+ (z=0) / n- 8.25e13 / n+ (z=T), x,y-invariant")
    # (1) fully-depleted linear problem on a graded grid: nodal exact parabola
    zs = np.unique(np.r_[np.linspace(0, T, 31), np.geomspace(1e-7, T, 40)])
    p = Poisson3D(slab(T, N, zs))
    V = 200.0
    r = p.solve({"cathode": V}, carriers=False)
    psi = r.psi[:, 1, 1]
    a0 = -p.Vt * math.log(1e19 / p.ni)
    b0 = V + p.Vt * math.log(1e19 / p.ni)
    exact = a0 + (b0 - a0) * zs / T + Q * N / (2 * EPS) * zs * (T - zs)
    print(
        f"  carriers=none, {len(zs)} graded z lines: max|psi - parabola| {np.abs(psi - exact).max():.1e} V"
    )
    # (2) Boltzmann carriers: depletion width vs the abrupt one-sided junction
    T = 30e-4  # thick enough that the depletion stays inside the n- up to 20 V
    print("  Boltzmann carriers, 30 um n-, 1201 uniform z lines")
    p = Poisson3D(slab(T, N, np.linspace(0, T, 1201)))
    p.equilibrium()
    Vbi = p.Vt * math.log(1e19 * N / p.ni**2)
    for V in (0.0, 2.0, 5.0, 10.0, 20.0):
        r = p.solve({"cathode": V})
        n = r.n[:, 1, 1]
        i = np.nonzero(n >= 0.5 * N)[0][0]
        z0, z1 = p.zs[i - 1], p.zs[i]
        W = z0 + (0.5 * N - n[i - 1]) * (z1 - z0) / (n[i] - n[i - 1])
        Wa = math.sqrt(2 * EPS * (Vbi + V - 2 * p.Vt) / (Q * N))
        print(f"  {V:5.1f} V: W(n = N/2) {W * 1e4:.4f} um, analytic sqrt(2 eps (Vbi + V - 2kT/q)/qN) "
              f"{Wa * 1e4:.4f} um, diff {(W - Wa) * 1e4:+.4f} um (L_D {math.sqrt(EPS * p.Vt / (Q * N)) * 1e4:.3f} um)")  # fmt: skip


def coax_analytic(r, rin, rout, psi_in, psi_out, N):
    k = Q * N / (4 * EPS)
    A = (psi_out - psi_in + k * (rout**2 - rin**2)) / math.log(rout / rin)
    B = psi_in + k * rin**2 - A * math.log(rin)
    return -k * r**2 + A * np.log(r) + B


def coax():
    rin, rout, N, V = 5e-4, 20e-4, 8.252732e13, 100.0
    print(
        f"[coax] r_in 5 um, r_out 20 um (staircase circles), N_D 8.25e13, {V} V, carriers=none"
    )
    rows = []
    for h in (1e-4, 0.5e-4, 0.25e-4, 0.125e-4):
        p = Poisson3D(coax_cell(rin, rout, N_D=N, h=h, nz=2))
        r = p.solve({"cathode": V}, carriers=False)
        psi_in = V + p.Vt * math.log(1e19 / p.ni)
        psi_out = -p.Vt * math.log(1e19 / p.ni)
        Z, Y, X = np.meshgrid(p.zs, p.ys, p.xs, indexing="ij")
        R = np.hypot(X, Y)
        m = (p.free != 0) & (R > rin + 2e-4) & (R < rout - 2e-4)
        ex = coax_analytic(R[m], rin, rout, psi_in, psi_out, N)
        err = np.abs(r.psi[m] - ex).max()
        # scalar functional: potential at r = 12.5 um on the x axis
        i = np.argmin(np.abs(p.xs - 12.5e-4))
        rows.append((h, p.N, err, r.psi[0, 0, i]))
        print(f"  h {h * 1e4:5.3f} um: N {p.N:7d}, max|psi - analytic| (r_in+2 < r < r_out-2 um) {err:.3e} V, "
              f"psi(12.5 um) {r.psi[0, 0, i]:.6f} V (analytic {coax_analytic(12.5e-4, rin, rout, psi_in, psi_out, N):.6f})")  # fmt: skip
    e = np.array([r[2] for r in rows])
    print("  observed order of the max error:", np.round(np.log2(e[:-1] / e[1:]), 2))


def coax_vs_axisym():
    """Boltzmann carriers: potential of the 3D (cut-cell) coax vs the fv2d (r, z)
    solution of the same coax (the axisymmetric reference of a truly cylindrical cell).
    """
    from etna.core.fv2d import Solver2D
    from etna.core.fv2d.structures import (
        Structure,
        _base_info,
        build_mesh,
        graded_lines,
    )

    rin, rout, N = 5e-4, 20e-4, 8.252732e13
    dz = 1e-4

    def doping(X, Y):
        return (
            np.where(X <= rin, 1e19, N),
            np.where(X >= rout, 1e19 + N, 0.0),
            np.where(X >= rout, 1e19 + N, 0.0),
        )

    info = _base_info(300.0, dz, 0.0, rout, "coax", "cylindrical")
    st = Structure("coax", rout, dz, doping, {"cathode": lambda X, Y: X <= rin, "anode": lambda X, Y: X >= rout},
                   [(0, 0.05e-4), (rin, 0.02e-4), (rout, 0.02e-4)], [(0, dz / 2), (dz, dz / 2)],
                   "cylindrical", 300.0, info)  # fmt: skip
    mesh = build_mesh(st, graded_lines(st.x_knots, 0, rout, h_max=0.05e-4),
                      graded_lines(st.y_knots, 0, dz))  # fmt: skip
    s2 = Solver2D(mesh).equilibrium()
    print("[coax vs axisym fv2d] Boltzmann carriers: max|psi_3D(x, y) - psi_rz(r)| over the "
          f"annulus (fv2d reference: {mesh.nx} r lines, <= 0.05 um)")  # fmt: skip
    Vs = (5, 20, 50)
    ref = []
    for V in Vs:
        s2.set_bias({"cathode": V})
        ref.append(s2.potential()[0].copy())
    rows = []
    for h in (1e-4, 0.5e-4, 0.25e-4, 0.125e-4):
        p = Poisson3D(coax_cell(rin, rout, N_D=N, h=h, nz=2))
        p.equilibrium()
        Z, Y, X = np.meshgrid(p.zs, p.ys, p.xs, indexing="ij")
        R = np.hypot(X, Y)
        m = (p.free != 0) & (R > rin) & (R < rout)
        errs = []
        for V, pr in zip(Vs, ref):
            r = p.solve({"cathode": V})
            errs.append(np.abs(r.psi[m] - np.interp(R[m], mesh.xs, pr)).max())
        rows.append(errs)
        print(f"  3D h {h * 1e4:.3f} um ({p.N} nodes):",
              ", ".join(f"{V} V {e:.2e} V" for V, e in zip(Vs, errs)))  # fmt: skip
    e = np.array(rows)
    print("  observed order per bias:", np.round(np.log2(e[:-1] / e[1:]), 2).T.tolist())
    print("  NOTE: not convergent by design. The difference is the nm-scale carrier spill-over\n"
          "  layer at the ohmic electrode surfaces (L_D(1e19) ~ 1 nm, a ~0.05-0.1 V drop), which\n"
          "  neither grid resolves: the fv2d reference itself moves by 0.045 V between 0.05 and\n"
          "  0.01 um surface spacing (scratch check). Without carriers (carriers=False) the same\n"
          "  coax converges at second order to the analytic solution (section [coax] above).")  # fmt: skip


def manufactured():
    """Laplace u = cos(pi x / Lx) cos(pi y / Ly) cosh(g z) / cosh(g Lz) (zero normal
    derivative on the x, y sides = the mirror planes), Dirichlet planes z = 0, Lz."""
    Lx, Ly, Lz = 1.0, 1.0, 1.0
    g = math.pi * math.sqrt(1 / Lx**2 + 1 / Ly**2)
    print("[manufactured Laplace] graded lines (ratio 1.1), max nodal error vs h")
    errs = []
    for n in (9, 17, 33, 65):
        t = np.linspace(0, 1, n)
        xs = Lx * (
            t + 0.15 * np.sin(np.pi * t) ** 2 * np.sign(0.5 - t) * 0.2
        )  # smooth grading
        xs = np.sort(xs)
        ys = Ly * t**1.3
        zs = Lz * (1 - (1 - t) ** 1.2)
        Z, Y, X = np.meshgrid(zs, ys, xs, indexing="ij")
        u = (
            np.cos(math.pi * X / Lx)
            * np.cos(math.pi * Y / Ly)
            * np.cosh(g * Z)
            / math.cosh(g * Lz)
        )
        free = np.ones(u.shape, np.uint8)
        free[0] = free[-1] = 0
        mg = Multigrid(xs, ys, zs, free)
        mg.set_c(np.zeros(u.shape))
        ubc = np.where(free == 0, u, 0.0)
        F = np.empty(u.shape)
        apply_op(ubc, F, *mg.fine.geo, np.zeros(u.shape), free)
        d, it, _ = pcg(mg, -F, rtol=1e-12)
        e = np.abs(ubc + d - u).max()
        errs.append(e)
        print(f"  n {n:3d}: max error {e:.3e}  (PCG {it} it)")
    e = np.array(errs)
    print("  observed order:", np.round(np.log2(e[:-1] / e[1:]), 2))


if __name__ == "__main__":
    t0 = time.perf_counter()
    for f in (parity_fv2d, slab_1d, coax, coax_vs_axisym, manufactured):
        f()
    print(f"total {time.perf_counter() - t0:.1f} s")
