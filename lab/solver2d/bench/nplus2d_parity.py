"""2D engine WITH an n+ layer vs devsim: dark state (psi, n, p), W, lateral CCE.

Engine side uses the engine internals directly (Device2D on grid_from_devsim with the
n+ keywords; the public gate is bypassed) so this runs before/independently of the
production gate. Generation exactly as cce_lateral_scan: alpha_depth_shape stripe,
zero in the p+, confine_generation_to_epi (unless nplus_generation). devsim side:
create_2d_dd_device(...) + dark dump + cce_lateral_scan (balanced extraction).

usage: uv run python lab/solver2d/bench/nplus2d_parity.py START STOP [out.csv]
appends rows to lab/solver2d/notes/nplus2d_parity.csv (default)
"""

from __future__ import annotations

import csv
import sys
import time
import uuid
from pathlib import Path

import devsim
import numpy as np

from etna.core.batched.extract import depletion_width
from etna.core.charge_collection_2d import (
    _ALPHA_RANGE_CM,
    cce_lateral_scan,
    create_2d_dd_device,
)
from etna.core.engine2d.device import Device2D
from etna.core.engine2d.grid import grid_from_devsim
from etna.core.generation_profiles import alpha_depth_shape

D = 8.823059e13
# (hw um, epi um, N_D_bulk, V, n+ um, N_D_nplus, nplus_generation)
POINTS = [
    (25, 10, D, 20, 1, 1e19, False),
    (25, 10, D, 60, 1, 1e19, False),
    (25, 10, D, 100, 1, 1e19, False),
    (25, 10, 5e13, 80, 5, 1e18, False),
    (25, 5, 5e13, 10, 1, 1e18, False),
    (25, 5, 5e14, 60, 3, 1e19, False),
    (25, 5, D, 100, 5, 1e19, False),
    (25, 7.5, 5e14, 40, 2, 3e18, False),
    (25, 15, 5e13, 50, 1, 1e19, False),
    (25, 20, D, 100, 1, 1e19, False),
    (25, 20, 5e14, 10, 5, 1e18, False),
    (25, 20, 5e14, 80, 1, 1e18, False),
    (25, 10, 5e14, 5, 1, 1e19, False),
    (25, 10, D, 30, 1, 1e19, True),
    (25, 10, D, 80, 5, 1e18, True),
    (50, 10, D, 50, 1, 1e19, False),
    (50, 5, 5e13, 100, 5, 1e19, False),
    (50, 20, 5e13, 40, 2, 1e18, True),
    (100, 10, 5e13, 80, 1, 1e19, False),
    (150, 10, D, 60, 1, 1e19, False),
    (150, 20, 5e14, 100, 5, 1e18, False),
]
NPTS = 10


def generation(ys, xs, junction, epi_end, x_pos, keep_nplus, gen_rate=1e18, sigma=2e-4):
    yj = ys - junction
    sh = alpha_depth_shape(np.clip(yj, 0.0, None), alpha_range_cm=_ALPHA_RANGE_CM)
    sh = sh / sh.max()
    row = gen_rate * sh
    row[yj < 0] = 0.0
    if not keep_nplus:
        row[ys > epi_end + 1e-12] = 0.0
    lat = np.exp(-0.5 * ((xs - x_pos) / sigma) ** 2)
    return (row[:, None] * lat[None, :]).ravel()


def engine(kw, hw, V):
    t0 = time.perf_counter()
    grid, info = grid_from_devsim(half_width_um=hw, **kw)
    dev = Device2D(grid, info)
    dev.dark(float(V))
    dark = dict(u=dev.u * dev.Vt, n=dev.n.copy(), p=dev.p.copy())
    dev.factor_continuity()
    xpos = np.linspace(0, info["half_width_cm"], NPTS)
    raw = []
    for xp in xpos:
        G = generation(
            grid.ys,
            grid.xs,
            info["substrate_thickness_cm"],
            info["epi_end"],
            xp,
            kw.get("nplus_generation", False),
        )
        I = dev.illuminated_current(G)
        raw.append(abs(I) / float(np.sum(grid.volume * G)))
    t = time.perf_counter() - t0
    return grid, info, dark, np.array(raw), t, max(dev.gummel_iterations)


def ref(kw, hw, V):
    t0 = time.perf_counter()
    name = f"np2d_{uuid.uuid4().hex[:8]}"
    info = create_2d_dd_device(half_width_um=hw, V_bias=V, device_name=name, **kw)
    try:
        g = lambda m: np.array(
            devsim.get_node_model_values(
                device=name, region=info["region_name"], name=m
            )
        )
        dark = {
            k: g(m)
            for k, m in (
                ("x", "x"),
                ("y", "y"),
                ("u", "Potential"),
                ("n", "Electrons"),
                ("p", "Holes"),
                ("donors", "Donors"),
            )
        }
        scan = cce_lateral_scan(info, n_points=NPTS)
    finally:
        devsim.delete_device(device=name)
    return dark, np.asarray(scan["cce_values"]), time.perf_counter() - t0


def W_col(ys, n, donors, info):
    return float(
        depletion_width(
            ys,
            n[None, :],
            donors[None, :],
            info["substrate_thickness_cm"],
            info["epi_thickness_cm"],
        )[0]
    )


def logdiff(a, b):
    m = (a > 1.0) & (b > 1.0)
    return float(np.max(np.abs(np.log10(a[m]) - np.log10(b[m]))))


FIELDS = [
    "hw",
    "epi_um",
    "nd",
    "V",
    "np_um",
    "nd_np",
    "np_gen",
    "W_d_um",
    "dW_um",
    "dpsi_V",
    "dlogn",
    "dlogp",
    "center_d",
    "dC",
    "dE",
    "max_dCCE",
    "t_d",
    "t_e",
    "gummel_max",
]


def main():
    a, b = int(sys.argv[1]), int(sys.argv[2])
    out = Path(
        sys.argv[3] if len(sys.argv) > 3 else "lab/solver2d/notes/nplus2d_parity.csv"
    )
    new = not out.exists()
    with out.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(FIELDS)
        for hw, epi, nd, V, npu, ndn, gen in POINTS[a:b]:
            kw = dict(
                epi_thickness_cm=epi * 1e-4,
                N_D_bulk=nd,
                nplus_thickness_cm=npu * 1e-4,
                N_D_nplus=ndn,
                nplus_generation=gen,
            )
            grid, info, de, cce_e, te, git = engine(kw, hw, float(V))
            dd, cce_d, td = ref(kw, hw, float(V))
            # devsim arrays -> engine order k
            inv = np.empty_like(grid.perm)
            inv[grid.perm] = np.arange(grid.N)
            dk = {k: v[inv] for k, v in dd.items()}
            assert np.allclose(dk["y"], np.repeat(grid.ys, grid.nx))
            col = slice(0, None, grid.nx)
            Wd = W_col(grid.ys, dk["n"][col], dk["donors"][col], info)
            We = W_col(grid.ys, de["n"][col], np.repeat(grid.donors, 1), info)
            dCCE = cce_e - np.clip(cce_d, 0, 1)
            row = [
                hw,
                epi,
                nd,
                V,
                npu,
                ndn,
                gen,
                Wd * 1e4,
                (We - Wd) * 1e4,
                float(np.max(np.abs(de["u"] - dk["u"]))),
                logdiff(de["n"], dk["n"]),
                logdiff(de["p"], dk["p"]),
                cce_d[0],
                dCCE[0],
                dCCE[-1],
                float(np.max(np.abs(dCCE))),
                round(td, 2),
                round(te, 3),
                git,
            ]
            w.writerow(row)
            f.flush()
            print("ROW", row, flush=True)


if __name__ == "__main__":
    main()
