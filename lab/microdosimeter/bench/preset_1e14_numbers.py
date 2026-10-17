"""Production numbers of the microdosimeter preset at N_D = 1e14 (candidate change after
bench/hypotheses.py): run_cv reach-through and plateau, tilted CCE pins, cold start.

usage: uv run python lab/microdosimeter/bench/preset_1e14_numbers.py [N_D]
"""

import sys

import numpy as np

import etna
from etna.api import simulation as sim

N_D = float(sys.argv[1]) if len(sys.argv) > 1 else 1e14
for p in (0.3, 1.0):
    c = etna.microdosimeter_config(p, N_D=N_D)
    r = etna.run_cv(c, 0.0, -60.0, 121)
    W = np.array(r.metadata["W_from_C"])[[0, 20, 60, 120]] * 1e4
    print(
        f"p+{p}: solver {r.metadata['solver']} reach_through_V {r.metadata['reach_through_V']} "
        f"C(-60) {r.y[-1]:.5e} F  x4/25 {r.y[-1] * 4 / 25:.4e}  W_from_C(0,-10,-30,-60) um {np.round(W, 3)}"
    )
c = etna.microdosimeter_config(N_D=N_D)
r = etna.run_cce_tilted_alpha(c, 0.0, -60.0, 13)
V = np.round(r.x, 6)
n = r.y / r.y[V == -50.0][0]
print("tilted raw:", {float(v): round(float(y), 5) for v, y in zip(V, r.y)})
print("tilted norm -50:", {float(v): round(float(y), 5) for v, y in zip(V, n)})
print(
    "cce_vs_E0(-50)",
    round(float(r.metadata["cce_vs_E0"][V == -50.0][0]), 5),
    "reach",
    r.metadata["reach_through_V"],
    "monotone",
    bool(np.all(np.diff(r.y) >= -1e-6)),
)
p = etna.run_cce_tilted_alpha(c, -50.0, -50.0, 1, generation_in_p_plus=True)
print("p+epi -50:", round(float(p.y[0]), 5))
r2 = etna.run_cce(c, -50.0, -60.0, 3)
d = sim._run_cce_devsim(c, np.array([-50.0, -55.0, -60.0]))
print(
    "run_cce solver",
    r2.metadata["solver"],
    "fallback_reason" in r2.metadata,
    "max|eng-devsim|",
    float(np.max(np.abs(r2.y - d.y))),
    "reach",
    r2.metadata["reach_through_V"],
)
a = etna.run_cce(c, -10.0, -40.0, 7)
b = etna.run_cce_tilted_alpha(
    c, -10.0, -40.0, 7, theta_deg=0.0, alpha_range_um=15.0, shape="box"
)
print("normal-incidence box vs run_cce max diff", float(np.max(np.abs(a.y - b.y))))
