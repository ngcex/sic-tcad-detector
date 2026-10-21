"""Grid convergence of the 3D column cell (pitch 50 um, r = 5 um columns, N_D 8.25e13):
bias of 99 % depletion (node-threshold and charge criteria) and the low-field volume
(|E| < 1 kV/cm) at 150 V (the saddle minimum itself is a symmetric NULL point, |E| -> 0), on line sets scaled by s (h_fine = 0.25 s um, h_bulk = 2 s um),
observed order from three levels and Richardson-extrapolated values with error estimates.

A: through columns (z-invariant, the issue #4 case; 2 um thick, compare 2D fv2d 75 V)
B: single-sided 30 um deep columns + planar p+ back contact, 50 um thick (true 3D)

Usage: uv run python lab/fv3d/bench/convergence.py [A|B]
"""

import sys
import time

import numpy as np

from etna.core.fv3d import Poisson3D, column_cell
from etna.core.fv3d.analysis import depletion_voltage

CASES = {
    "A": (dict(thickness=2e-4, h_z=1e-4), (2.0, 1.41421356, 1.0, 0.70710678)),
    "B": (dict(depth_n=30e-4, depth_p=30e-4, back="p+"), (2.0, 1.41421356, 1.0)),
}


def run(case):
    kw, scales = CASES[case]
    rows = []
    for s in scales:
        st = column_cell(h_fine=0.25e-4 * s, h_bulk=2e-4 * s,
                         **{k: (v * s if k == "h_z" else v) for k, v in kw.items()})  # fmt: skip
        t0 = time.perf_counter()
        p = Poisson3D(st)
        p.equilibrium()
        V99, info = depletion_voltage(p, target=0.99, V_step=20.0, tol=0.02)
        V99c, _ = depletion_voltage(p, target=0.99, V_step=20.0, tol=0.02, mode="charge")
        sol = p.solve({"cathode": 150.0})
        sad = p.saddles(sol, 1, exclude_contact_layers=2)
        Emin = sad[0]["E"] if sad else np.nan
        low = p.low_field_fraction(sol, 1e3)
        wall = time.perf_counter() - t0
        rows.append((s, p.N, V99, V99c, low, Emin, wall))
        print(f"{case} s {s:.3f}: grid {st.shape} = {p.N} nodes, cut edges {p.n_cut}, "
              f"V99 {V99:.3f} V (charge criterion {V99c:.3f} V), 150 V: bulk fraction with "
              f"|E| < 1 kV/cm {low:.5f}, min |E| {Emin:.2f} V/cm at "
              f"({sad[0]['x'] * 1e4:.1f}, {sad[0]['y'] * 1e4:.1f}, {sad[0]['z'] * 1e4:.1f}) um, "
              f"{info['newton']} Newton / {info['cg']} CG (first scan), {wall:.1f} s", flush=True)  # fmt: skip
    a = np.array([r[:5] for r in rows])
    for name, col in (("V99", 2), ("V99 (charge)", 3), ("|E| < 1 kV/cm fraction", 4)):
        q = a[:, col]
        for i in range(len(q) - 2):
            q1, q2, q3 = q[i : i + 3]
            r = a[i, 0] / a[i + 1, 0]
            num, den = q1 - q2, q2 - q3
            pobs = np.log(num / den) / np.log(a[i, 0] / a[i + 1, 0]) if den and num / den > 0 else -1
            if not 0.3 < pobs < 4.0:
                print(f"  {name}: levels {i}-{i + 2} not asymptotic ({q1:.4f}, {q2:.4f}, {q3:.4f}): "
                      f"error estimate |q2 - q3| = {abs(q2 - q3):.3g}")  # fmt: skip
                continue
            pobs = np.log(num / den) / np.log(r)
            ext = q3 - den / (r**pobs - 1)
            print(f"  {name}: observed order {pobs:.2f}, Richardson {ext:.4f}, "
                  f"error of the finest {abs(q3 - ext):.3g}")  # fmt: skip


if __name__ == "__main__":
    for c in sys.argv[1:] or ["A", "B"]:
        run(c)
