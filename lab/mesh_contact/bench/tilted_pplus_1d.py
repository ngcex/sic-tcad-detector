"""run_cce_tilted_alpha(generation_in_p_plus=True) is the one default-geometry 1D path
that puts generation on a contact node (the anode at x = 0, on the p+). Does its anode
contact cell (100 nm) matter? Microdosimeter preset, 65 deg, Bragg shape, 0 dead layer;
the anode cell bisected 0-6 times (contact-only), batched engine, trapz normalisation as
the facade.

Usage: uv run python lab/mesh_contact/bench/tilted_pplus_1d.py
"""

import numpy as np

from etna import microdosimeter_config
from etna.api.device import device_kwargs
from etna.api.simulation import _engine_batch
from etna.core.batched.mesh import node_positions
from etna.core.batched.sweeps import cce_sweep
from etna.core.generation_profiles import tilted_alpha_generation

V = np.array([-1.0, -5.0, -10.0, -30.0])
cfg = microdosimeter_config()
kw = device_kwargs(cfg)
sub, epi = kw["substrate_thickness_cm"], kw["epi_thickness_cm"]
x0 = node_positions(epi, sub, kw["nplus_thickness_cm"])


def anode_refine(x, times):
    out = set(x)
    o = x[1]
    for _ in range(times):
        o = 0.5 * (x[0] + o)
        out.add(o)
    return np.array(sorted(out))


prev = None
for lev in range(7):
    x = anode_refine(list(x0), lev)
    pb, x = _engine_batch([kw], x)
    G, _ = tilted_alpha_generation(x, theta_deg=65.0, alpha_range_cm=18e-4,
                                   dead_layer_cm=0.0, shape="bragg", E0_eV=5.5e6)  # fmt: skip
    G[x > sub + epi + 1e-12] = 0.0
    G = G * (1e18 / G.max())
    r = cce_sweep(pb, V, sub, G, T=cfg.T, N_A=cfg.N_A)
    c = r.cce[:, 0]
    d = "" if prev is None else np.array2string(c - prev, precision=2)
    print(
        f"anode cell {(x[1] - x[0]) * 1e7:7.2f} nm  CCE {np.array2string(c, precision=6)}  d {d}"
    )
    prev = c
