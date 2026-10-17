"""Coarser 2D y-mesh around the n-/n+ step (default stack): node counts + engine gate.

The default stack (n+ under the epi) adds a 1 nm mesh line at the epi end. devsim grades
the spacing geometrically between neighbouring y-lines, so the WHOLE interval from the
last epi line (epi_mid2 = junction + 5 um, ps 50 nm; or epi_mid1 for thin epis) to the
epi end is graded 50 nm -> 1 nm at ratio ~1.01: ~330 extra rows on a 10 um epi. The n+
itself is only ~47 rows (1 nm -> 100 nm at ratio 1.1 over 1 um).

Candidates add lines only in the n+ branch (legacy mesh untouched):
  epi_end - d_e (ps p_e)    : grade to 1 nm over d_e instead of microns
  epi_end + d_n (ps p_n)    : same on the n+ side (matters for thick n+)
The lines are injected by patching devsim.add_2d_mesh_line during the build (lab only).

usage:
  uv run python lab/solver2d/bench/nplus_mesh.py count            # node counts per candidate
  uv run python lab/solver2d/bench/nplus_mesh.py gate CAND [out]  # engine2d old vs CAND
"""

from __future__ import annotations

import contextlib
import csv
import sys
import time
from pathlib import Path

import devsim
import numpy as np

from etna.core.batched.extract import depletion_width
from etna.core.charge_collection_2d import _ALPHA_RANGE_CM
from etna.core.device import DEFAULT_STACK
from etna.core.engine2d.device import Device2D
from etna.core.engine2d.grid import grid_from_devsim
from etna.core.generation_profiles import alpha_depth_shape

# name -> list of (offset from epi_end in cm, ps in cm); negative offset = epi side.
CANDIDATES = {
    "fine": [],
    "e05": [(-0.5e-4, 5e-6)],
    "e05n02": [(-0.5e-4, 5e-6), (0.2e-4, 2e-6)],
    "e1n02": [(-1.0e-4, 5e-6), (0.2e-4, 2e-6)],
    "e02n02": [(-0.2e-4, 2e-6), (0.2e-4, 2e-6)],
    "e05n01": [(-0.5e-4, 5e-6), (0.1e-4, 1e-6)],
    "e1n05": [(-1.0e-4, 1e-5), (0.5e-4, 5e-6)],
    # keeps few-nm cells over the last ~0.3 um (where the V >= V_rt depletion edge sits)
    "k3": [(-0.1e-4, 2e-7), (-0.3e-4, 5e-7), (-0.8e-4, 1.5e-6), (-2e-4, 4e-6), (0.2e-4, 2e-6)],
    # refinements, for the self-convergence check (not candidates)
    "sfine": [(-0.5e-4, 2e-7), (-2e-4, 1e-6)],
    "ultra": [(-0.05e-4, 5e-8), (-0.2e-4, 1e-7), (-0.5e-4, 2e-7), (-1e-4, 3e-7),
              (-2e-4, 5e-7), (-3e-4, 1e-6), (0.05e-4, 5e-8), (0.2e-4, 2e-7)],
}


@contextlib.contextmanager
def mesh_variant(extra, epi_end):
    """Inject y-lines at epi_end + offset (ps) into the next create_sic_2d_device mesh."""
    orig = devsim.add_2d_mesh_line
    done = set()

    def patched(**kw):
        orig(**kw)
        mesh = kw["mesh"]
        if kw["dir"] == "y" and mesh not in done and abs(kw["pos"] - epi_end) < 1e-12:
            done.add(mesh)
            for off, ps in extra:
                orig(mesh=mesh, dir="y", pos=epi_end + off, ps=ps)

    devsim.add_2d_mesh_line = patched
    try:
        yield
    finally:
        devsim.add_2d_mesh_line = orig


def build(cand, **kw):
    sub = kw.get("substrate_thickness_cm", DEFAULT_STACK["substrate_thickness_cm"])
    epi_end = sub + kw.get("epi_thickness_cm", 10e-4)
    t0 = time.perf_counter()
    with mesh_variant(CANDIDATES[cand], epi_end):
        grid, info = grid_from_devsim(**kw)
    return grid, info, time.perf_counter() - t0


def generation(grid, info, x_pos, gen_rate=1e18, sigma=2e-4):
    ys, xs = grid.ys, grid.xs
    yj = ys - info["substrate_thickness_cm"]
    sh = alpha_depth_shape(np.clip(yj, 0.0, None), alpha_range_cm=_ALPHA_RANGE_CM)
    row = gen_rate * sh / sh.max()
    row[yj < 0] = 0.0
    row[ys > info["epi_end"] + 1e-12] = 0.0
    lat = np.exp(-0.5 * ((xs - x_pos) / sigma) ** 2)
    return (row[:, None] * lat[None, :]).ravel()


def solve(grid, info, V, npts=5):
    t0 = time.perf_counter()
    dev = Device2D(grid, info)
    dev.dark(float(V))
    nx = grid.nx
    u = dev.u.reshape(grid.ny, nx)[:, 0] * dev.Vt
    n = dev.n.reshape(grid.ny, nx)[:, 0]
    W = float(
        depletion_width(
            grid.ys,
            n[None, :],
            grid.donors[None, :],
            info["substrate_thickness_cm"],
            info["epi_thickness_cm"],
        )[0]
    )
    dev.factor_continuity()
    cce = []
    for xp in np.linspace(0, info["half_width_cm"], npts):
        G = generation(grid, info, xp)
        cce.append(abs(dev.illuminated_current(G)) / float(np.sum(grid.volume * G)))
    I0 = dev.illuminated_current(np.zeros(grid.N))  # dark (SRH) cathode current / q
    return dict(W=W, cce=np.array(cce), psi=u, I0=I0, t=time.perf_counter() - t0)


def epi_compare(ga, ia, sa, gb, sb):
    """psi / field of b interpolated onto a's epi nodes (a = fine reference)."""
    jp, ee = ia["substrate_thickness_cm"], ia["epi_end"]
    m = (ga.ys >= jp) & (ga.ys <= ee)
    pb = np.interp(ga.ys[m], gb.ys, sb["psi"])
    dpsi = float(np.max(np.abs(pb - sa["psi"][m])))
    Ea = -np.diff(sa["psi"]) / np.diff(ga.ys)
    Eb = -np.diff(sb["psi"]) / np.diff(gb.ys)
    ya = 0.5 * (ga.ys[1:] + ga.ys[:-1])
    yb = 0.5 * (gb.ys[1:] + gb.ys[:-1])
    me = (ya > jp + 0.05e-4) & (ya < ee - 0.05e-4)  # epi interior, away from both steps
    Eb_i = np.interp(ya[me], yb, Eb)
    dE = float(np.max(np.abs(Eb_i - Ea[me])) / np.max(np.abs(Ea[me])))
    return dpsi, dE


# (hw um, epi um, n+ um, V)
GATE = [
    (25, 10, 1, 5),
    (25, 10, 1, 10),
    (25, 10, 1, 20),
    (25, 10, 1, 50),
    (25, 10, 1, 100),
    (25, 5, 1, 5),
    (25, 5, 1, 50),
    (25, 5, 3, 100),
    (25, 20, 1, 5),
    (25, 20, 1, 30),
    (25, 20, 1, 100),
    (25, 20, 5, 60),
    (25, 10, 5, 5),
    (25, 10, 5, 50),
    (25, 10, 2, 30),
    (50, 10, 1, 50),
    (50, 5, 5, 100),
    (50, 20, 1, 40),
    (150, 10, 1, 50),
    (150, 20, 5, 100),
    (150, 10, 1, 5),
]
FIELDS = [
    "cand",
    "hw",
    "epi",
    "np",
    "V",
    "ny_f",
    "ny_c",
    "N_f",
    "N_c",
    "W_f_um",
    "relW",
    "relCCE",
    "relI0",
    "I0_f",
    "dpsi_V",
    "relE",
    "t_f",
    "t_c",
]


def cmd_count():
    for cand in CANDIDATES:
        row = [cand]
        for hw, epi, npu in (
            (150, 10, 1),
            (25, 5, 1),
            (25, 20, 1),
            (25, 10, 5),
            (25, 20, 5),
        ):
            g, i, t = build(
                cand,
                half_width_um=hw,
                epi_thickness_cm=epi * 1e-4,
                nplus_thickness_cm=npu * 1e-4,
            )
            row.append(f"hw{hw}/e{epi}/n{npu}: ny {g.ny} N {g.N}")
        print(" | ".join(row), flush=True)
    g, i, t = build("fine", half_width_um=150, nplus_thickness_cm=None)
    print(f"legacy-geometry stack, no n+, hw150: ny {g.ny} N {g.N}")


def cmd_gate(cand, out):
    out = Path(out)
    new = not out.exists()
    with out.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(FIELDS)
        for hw, epi, npu, V in GATE:
            kw = dict(
                half_width_um=hw,
                epi_thickness_cm=epi * 1e-4,
                nplus_thickness_cm=npu * 1e-4,
            )
            ga, ia, _ = build("fine", **kw)
            gb, ib, _ = build(cand, **kw)
            sa, sb = solve(ga, ia, V), solve(gb, ib, V)
            dpsi, dE = epi_compare(ga, ia, sa, gb, sb)
            relW = abs(sb["W"] - sa["W"]) / sa["W"]
            relC = float(np.max(np.abs(sb["cce"] - sa["cce"]) / sa["cce"]))
            relI = (
                abs(sb["I0"] - sa["I0"]) / abs(sa["I0"]) if sa["I0"] else float("nan")
            )
            r = [
                cand,
                hw,
                epi,
                npu,
                V,
                ga.ny,
                gb.ny,
                ga.N,
                gb.N,
                f"{sa['W']*1e4:.6f}",
                f"{relW:.2e}",
                f"{relC:.2e}",
                f"{relI:.2e}",
                f"{sa['I0']:.3e}",
                f"{dpsi:.2e}",
                f"{dE:.2e}",
                f"{sa['t']:.2f}",
                f"{sb['t']:.2f}",
            ]
            w.writerow(r)
            f.flush()
            print(" ".join(map(str, r)), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "count":
        cmd_count()
    else:
        cmd_gate(
            sys.argv[2],
            (
                sys.argv[3]
                if len(sys.argv) > 3
                else "lab/solver2d/notes/nplus_mesh_gate.csv"
            ),
        )
