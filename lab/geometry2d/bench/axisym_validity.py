"""Validity of the axisymmetric (Wigner-Seitz) proxy of a 3D-column array (fv2d).

Columns through the whole epi (the classic 3D-detector layout): the plane normal to the
columns is a 2D Cartesian problem. Square lattice of pitch P: n+ cathode column (radius
rc) at the cell centre, p+ anode column (radius rp) at the corners, n- bulk N_D. By
symmetry the quarter cell [0, P/2]^2 with the n+ quarter-disc at (0, 0) and the p+
quarter-disc at (P/2, P/2) (zero-flux sides). The axisymmetric proxy (what
create_3d_electrode_device / electrode3d_structure model) is a cylinder of the same
area, R = P / sqrt(pi), with the p+ as an outer shell of the same area
(R_in = sqrt(R^2 - rp^2)). Both columns are Dirichlet bodies (ohmic electrodes).

Compared vs reverse bias: depleted area fraction (n < N_D / 2), the bias of 99 %
depletion, the min / max |E| over the depleted bulk. The square cell has saddle points of
low field at the edge midpoints (between two p+ columns) that a cylinder cannot have.

Usage: uv run python lab/geometry2d/bench/axisym_validity.py
"""

import numpy as np

from etna.core.fv2d import Solver2D
from etna.core.fv2d.structures import Structure, _base_info, build_mesh, graded_lines

P = 50e-4
RC = RP = 5e-4
ND = 8.252732e13
N_PLUS = 1e19


def square_cell():
    a = P / 2

    def doping(X, Y):
        r0 = np.hypot(X, Y)
        r1 = np.hypot(X - a, Y - a)
        don = np.where(r0 <= RC, N_PLUS, ND)
        acc = np.where(r1 <= RP, N_PLUS, 0.0)
        return don, acc, acc

    k = [(0.0, 2e-5), (RC, 2e-6), (a - RP, 2e-6), (a, 2e-5)]
    info = _base_info(300.0, P, 0.0, a, "square_cell", "cartesian")
    contacts = {
        "cathode": lambda X, Y: np.hypot(X, Y) <= RC,
        "anode": lambda X, Y: np.hypot(X - a, Y - a) <= RP,
    }
    return Structure(
        "square_cell", a, a, doping, contacts, k, k, "cartesian", 300.0, info
    )


def cylinder():
    R = P / np.sqrt(np.pi)
    Rin = np.sqrt(R * R - RP * RP)
    dz = 1e-4

    def doping(X, Y):
        don = np.where(X <= RC, N_PLUS, ND)
        acc = np.where(X >= Rin, N_PLUS, 0.0)
        return don, acc, acc

    info = _base_info(300.0, P, 0.0, R, "cylinder", "cylindrical")
    contacts = {"cathode": lambda X, Y: X <= RC, "anode": lambda X, Y: X >= Rin}
    xk = [(0.0, 2e-5), (RC, 2e-6), (Rin, 2e-6), (R, 2e-5)]
    return Structure(
        "cylinder",
        R,
        dz,
        doping,
        contacts,
        xk,
        [(0, dz / 2), (dz, dz / 2)],
        "cylindrical",
        300.0,
        info,
    )


def field(s):
    psi = s.potential()
    m = s.mesh
    Ey, Ex = np.gradient(psi, m.ys, m.xs)
    return np.hypot(Ex, Ey)


def sweep(st, Vs):
    mesh = build_mesh(
        st, graded_lines(st.x_knots, 0, st.width), graded_lines(st.y_knots, 0, st.depth)
    )
    s = Solver2D(mesh).equilibrium()
    bulk = (mesh.contact_any() == 0) & (mesh.donors < 1e18) & (mesh.acceptors == 0)
    w = mesh.vol * bulk
    out = []
    for V in Vs:
        s.set_bias({"cathode": V})
        n, _ = s.density()
        dep = (n < 0.5 * ND) & bulk
        E = field(s)
        Ed = E[dep]
        out.append((V, float(np.sum(w[dep]) / np.sum(w)), float(Ed.min()) if Ed.size else np.nan,
                    float(E[bulk].max())))  # fmt: skip
    return mesh, out


if __name__ == "__main__":
    Vs = [1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 200]
    res = {}
    for st in (square_cell(), cylinder()):
        mesh, out = sweep(st, Vs)
        res[st.name] = out
        print(f"{st.name}: {mesh.nx}x{mesh.ny} nodes")
        for V, f, emin, emax in out:
            print(
                f"  V {V:4.0f}  depleted {f:.4f}  min|E| depleted {emin:9.3e}  max|E| {emax:9.3e} V/cm"
            )
    for name, out in res.items():
        Vfd = next((V for V, f, *_ in out if f >= 0.99), None)
        print(name, "first sweep bias with >= 99 % of the bulk depleted:", Vfd)
