"""STEP 0: dump devsim's 2D mesh (geometry + couplings + equilibrium) to an npz.

usage: uv run python lab/solver2d/bench/dump_mesh.py HALF_WIDTH_UM [V_BIAS]
writes lab/solver2d/notes/mesh_hw{HW}.npz
"""

import sys
import time
from pathlib import Path

import devsim
import numpy as np

from etna.core.charge_collection_2d import create_2d_dd_device

hw = float(sys.argv[1]) if len(sys.argv) > 1 else 25.0
t0 = time.perf_counter()
info = create_2d_dd_device(half_width_um=hw, V_bias=0.0)
t1 = time.perf_counter()
dev, reg = info["device_name"], info["region_name"]
g = lambda name: np.array(
    devsim.get_node_model_values(device=dev, region=reg, name=name)
)
ge = lambda name: np.array(
    devsim.get_edge_model_values(device=dev, region=reg, name=name)
)
x, y = g("x"), g("y")
devsim.edge_from_node_model(device=dev, region=reg, node_model="node_index")
out = dict(
    x=x,
    y=y,
    NodeVolume=g("NodeVolume"),
    Donors=g("Donors"),
    Acceptors=g("Acceptors"),
    NetDoping=g("NetDoping"),
    Potential=g("Potential"),
    Electrons=g("Electrons"),
    Holes=g("Holes"),
    IntrinsicElectrons=g("IntrinsicElectrons"),
    EdgeCouple=ge("EdgeCouple"),
    EdgeLength=ge("EdgeLength"),
    edges=np.stack([ge("node_index@n0"), ge("node_index@n1")], axis=1).astype(int),
    elements=np.array(devsim.get_element_node_list(device=dev, region=reg)).reshape(
        -1, 3
    ),
)
for c in ("anode", "cathode"):
    # contact node models are evaluated only on the contact's nodes (0 elsewhere)
    devsim.contact_node_model(device=dev, contact=c, name=f"is_{c}", equation="1")
    out[f"contact_{c}"] = np.flatnonzero(g(f"is_{c}"))
out["regions"] = np.array(devsim.get_region_list(device=dev))
for r in devsim.get_region_list(device=dev):
    if r != reg:
        out[f"n_nodes_{r}"] = len(
            devsim.get_node_model_values(device=dev, region=r, name="x")
        )
out["params"] = np.array(
    [
        info["N_A_ionized"],
        info["n_i"],
        info["mu_n"],
        info["mu_p"],
        devsim.get_parameter(device=dev, region=reg, name="taun"),
        devsim.get_parameter(device=dev, region=reg, name="taup"),
        devsim.get_parameter(device=dev, region=reg, name="V_t"),
        devsim.get_parameter(device=dev, region=reg, name="Permittivity"),
    ]
)
Path("lab/solver2d/notes").mkdir(parents=True, exist_ok=True)
np.savez(f"lab/solver2d/notes/mesh_hw{hw:g}.npz", **out)
print(
    f"build {t1 - t0:.2f}s nodes {len(x)} edges {len(out['edges'])} tris {len(out['elements'])}"
)
