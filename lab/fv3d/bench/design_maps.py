"""Design maps with fv3d (N_D 8.25e13, 300 K):

1. 99 % depletion bias of the square column lattice (through columns, r_n = r_p) vs pitch
   and column radius, against the equal-area axisymmetric proxy (fv2d (r, z): n+ core,
   p+ outer shell of the same area as the corner columns), the proxy issue #4 found
   invalid for column lattices.
2. Square-cell microdosimeter (central n+ column r 5 um, p+ trench wall 5 um, through,
   10 um thick) vs its annular-trench proxy (the case where the (r, z) model IS meant to
   apply), pitch 25 / 50 / 80 um.
3. Low-field map of the true-3D column cell (30 um deep columns + p+ back plane, 50 um,
   pitch 50 um, r 5 um) at 100 V: null points (|E| local minima) and the bulk fraction
   below 0.5 / 1 / 2 kV/cm; figure lab/fv3d/figures/low_field_map.png (min over z of
   |E| in the depleted bulk, and the x = a mirror plane).

Usage: uv run python lab/fv3d/bench/design_maps.py [1|2|3]
"""

import sys
import time
from pathlib import Path

import numpy as np

from etna.core.fv3d import Poisson3D, column_cell, trench_cell
from etna.core.fv3d.analysis import depletion_voltage

ND = 8.252732e13
FIG = Path(__file__).resolve().parents[1] / "figures"


def proxy_v99(R, r_in, R_shell_in, target=0.99, V_max=2000.0):
    """V99 of the axisymmetric proxy: n+ core r <= r_in (cathode), p+ shell
    R_shell_in <= r <= R (anode), fv2d (r, z), 1 V continuation steps + interpolation.
    """
    from etna.core.fv2d import Solver2D
    from etna.core.fv2d.structures import (
        Structure,
        _base_info,
        build_mesh,
        graded_lines,
    )

    dz = 1e-4

    def doping(X, Y):
        acc = np.where(X >= R_shell_in, 1e19 + ND, 0.0)
        return np.where(X <= r_in, 1e19, ND), acc, acc

    info = _base_info(300.0, dz, 0.0, R, "proxy", "cylindrical")
    st = Structure("proxy", R, dz, doping,
                   {"cathode": lambda X, Y: X <= r_in, "anode": lambda X, Y: X >= R_shell_in},
                   [(0, 0.05e-4), (r_in, 0.02e-4), (R_shell_in, 0.02e-4), (R, 0.05e-4)],
                   [(0, dz / 2), (dz, dz / 2)], "cylindrical", 300.0, info)  # fmt: skip
    mesh = build_mesh(
        st,
        graded_lines(st.x_knots, 0, R, h_max=0.1e-4),
        graded_lines(st.y_knots, 0, dz),
    )
    s = Solver2D(mesh).equilibrium()
    bulk = mesh.contact_any() == 0
    w = mesh.vol * bulk
    prev = (0.0, 0.0)
    V = 0.0
    step = 1.0
    while V < V_max:
        V += step
        s.set_bias({"cathode": V})
        n, _ = s.density()
        f = float(np.sum(w[(n < 0.5 * ND) & bulk]) / np.sum(w))
        if f >= target:
            return prev[0] + (target - prev[1]) * (V - prev[0]) / (f - prev[1])
        prev = (V, f)
        step = 1.0 if V < 50 else 5.0
    return None


def map_columns():
    print(
        "1. through-column square lattice: V99 (V), 3D (fv3d) vs equal-area axisymmetric proxy"
    )
    print("   pitch um  r um   nodes   V99 3D   V99 proxy  ratio   wall s")
    for P in (30e-4, 50e-4, 80e-4):
        for r in (2.5e-4, 5e-4):
            t0 = time.perf_counter()
            st = column_cell(pitch=P, r_n=r, r_p=r, thickness=2e-4, h_z=1e-4)
            p = Poisson3D(st)
            p.equilibrium()
            V3, _ = depletion_voltage(
                p, V_step=max(5.0, 0.5 * (P / 1e-4) ** 2 * 0.03), tol=0.05
            )
            R = P / np.sqrt(np.pi)
            Vp = proxy_v99(R, r, np.sqrt(R * R - r * r))
            print(f"   {P * 1e4:6.0f} {r * 1e4:6.1f} {p.N:8d} {V3:8.2f} {Vp:10.2f} {V3 / Vp:6.2f} "
                  f"{time.perf_counter() - t0:7.1f}", flush=True)  # fmt: skip


def map_trench():
    print(
        "2. square-trench microdosimeter cell (r_n 5 um, wall 5 um, through, 10 um thick):"
    )
    print(
        "   pitch um   nodes   V99 3D   V99 annular proxy  ratio  min|E| at 1.5 V99 (V/cm)"
    )
    for P in (25e-4, 50e-4, 80e-4):
        st = trench_cell(pitch=P, r_n=5e-4, wall=5e-4, thickness=2e-4, h_z=1e-4)
        p = Poisson3D(st)
        p.equilibrium()
        V3, _ = depletion_voltage(p, V_step=5.0, tol=0.05)
        s = p.solve({"cathode": 1.5 * V3})
        sad = p.saddles(s, 1, exclude_contact_layers=2)
        # annular proxy: equal area of the SV inside the wall (side P - wall)
        side = P - 5e-4
        R = side / np.sqrt(np.pi)
        Vp = proxy_v99(R + 1e-4, 5e-4, R)  # 1 um p+ shell beyond the SV radius
        print(f"   {P * 1e4:6.0f} {p.N:8d} {V3:8.2f} {Vp:12.2f} {V3 / Vp:12.2f} "
              f"{sad[0]['E'] if sad else float('nan'):10.1f} at "
              f"({sad[0]['x'] * 1e4:.1f}, {sad[0]['y'] * 1e4:.1f}) um", flush=True)  # fmt: skip


def low_field():
    print("3. true-3D column cell (30 um columns, p+ back plane, 50 um), 100 V")
    st = column_cell(depth_n=30e-4, depth_p=30e-4, back="p+")
    p = Poisson3D(st)
    p.equilibrium()
    for V in (50.0, 100.0):
        s = p.solve({"cathode": V})
    print(f"   {st.shape} = {p.N} nodes, depleted {p.depleted_fraction(s):.5f}")
    for thr in (500.0, 1e3, 2e3):
        print(
            f"   bulk fraction with |E| < {thr:.0f} V/cm (or undepleted): {p.low_field_fraction(s, thr):.5f}"
        )
    print("   |E| local minima (x, y, z um; |E| V/cm):")
    for d in p.saddles(s, 8, exclude_contact_layers=2):
        print(
            f"     ({d['x'] * 1e4:5.1f}, {d['y'] * 1e4:5.1f}, {d['z'] * 1e4:5.1f})  {d['E']:9.2f}"
        )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    E = p.field(s.psi)[3]
    dep = p.depleted(s)
    Em = np.where(dep, E, np.nan)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    um = 1e4
    im = ax[0].pcolormesh(
        p.xs * um,
        p.ys * um,
        np.log10(np.nanmin(Em, axis=0)),
        shading="nearest",
        cmap="viridis",
    )
    ax[0].set(
        title="log10 min_z |E| (V/cm), depleted bulk, 100 V",
        xlabel="x (um)",
        ylabel="y (um)",
        aspect="equal",
    )
    fig.colorbar(im, ax=ax[0])
    im = ax[1].pcolormesh(
        p.ys * um, p.zs * um, np.log10(Em[:, :, -1]), shading="nearest", cmap="viridis"
    )
    ax[1].set(
        title="log10 |E| on the mirror plane x = a", xlabel="y (um)", ylabel="z (um)"
    )
    ax[1].invert_yaxis()
    fig.colorbar(im, ax=ax[1])
    fig.tight_layout()
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "low_field_map.png", dpi=110)
    print(f"   figure: {FIG / 'low_field_map.png'}")


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3"]
    for w, f in (("1", map_columns), ("2", map_trench), ("3", low_field)):
        if w in which:
            f()
