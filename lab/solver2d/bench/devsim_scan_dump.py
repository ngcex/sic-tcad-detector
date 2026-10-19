"""Dump devsim's dark + illuminated states and contact-current ingredients for the
lateral CCE scan (convention pinning for the 2D engine).

usage: uv run python lab/solver2d/bench/devsim_scan_dump.py HW_UM V_BIAS N_POINTS
writes lab/solver2d/notes/scan_hw{HW}_V{V}.npz
"""

import sys
import time

import devsim
import numpy as np

from etna.core.charge_collection import add_generation_to_dd
from etna.core.charge_collection_2d import (
    _ALPHA_RANGE_CM,
    _robust_dc_solve,
    compute_cce_2d,
    create_2d_dd_device,
    integrate_over_mesh_2d,
)
from etna.core.drift_diffusion import extract_contact_current
from etna.core.generation_profiles import alpha_depth_shape

hw = float(sys.argv[1])
V = float(sys.argv[2])
npts = int(sys.argv[3])
t0 = time.perf_counter()
info = create_2d_dd_device(half_width_um=hw, V_bias=V)
t1 = time.perf_counter()
dev, reg = info["device_name"], info["region_name"]
g = lambda n: np.array(devsim.get_node_model_values(device=dev, region=reg, name=n))
ge = lambda n: np.array(devsim.get_edge_model_values(device=dev, region=reg, name=n))
out = {f"dark_{k}": g(k) for k in ("Potential", "Electrons", "Holes")}
x, y = g("x"), g("y")
yj = y - info["substrate_thickness_cm"]
shape = alpha_depth_shape(np.clip(yj, 0.0, None), alpha_range_cm=_ALPHA_RANGE_CM)
shape = shape / shape.max()
xpos = np.linspace(0, info["half_width_cm"], npts)
for c in ("anode", "cathode"):
    devsim.contact_node_model(device=dev, contact=c, name=f"is_{c}", equation="1")
    out[f"contact_{c}"] = np.flatnonzero(g(f"is_{c}"))
cce, Ie, Ih, Igen = [], [], [], []
for k, xp in enumerate(xpos):
    G = 1e18 * shape * np.exp(-0.5 * ((x - xp) / 2e-4) ** 2)
    G[yj < 0] = 0.0
    add_generation_to_dd(info, G)
    _robust_dc_solve()
    cce.append(compute_cce_2d(info, G, "cathode"))
    Ie.append(
        devsim.get_contact_current(
            device=dev, contact="cathode", equation="ElectronContinuityEquation"
        )
    )
    Ih.append(
        devsim.get_contact_current(
            device=dev, contact="cathode", equation="HoleContinuityEquation"
        )
    )
    Igen.append(integrate_over_mesh_2d(info, G))
    if k in (0, npts - 1):
        for n in (
            "Potential",
            "Electrons",
            "Holes",
            "ElectronGeneration",
            "HoleGeneration",
            "RadGenRate",
        ):
            out[f"pt{k}_{n}"] = g(n)
        for n in ("ElectronCurrent", "HoleCurrent", "EdgeCouple"):
            out[f"pt{k}_{n}"] = ge(n)
        out[f"pt{k}_Ian"] = np.array(
            [
                devsim.get_contact_current(device=dev, contact="anode", equation=e)
                for e in ("ElectronContinuityEquation", "HoleContinuityEquation")
            ]
        )
    add_generation_to_dd(info, np.zeros_like(x))
    _robust_dc_solve()
t2 = time.perf_counter()
devsim.edge_from_node_model(device=dev, region=reg, node_model="node_index")
out["edges"] = np.stack([ge("node_index@n0"), ge("node_index@n1")], 1).astype(int)
out.update(
    x=x,
    y=y,
    xpos=xpos,
    cce=np.array(cce),
    I_e=np.array(Ie),
    I_h=np.array(Ih),
    I_gen_int=np.array(Igen),
    NodeVolume=g("NodeVolume"),
    t_build=t1 - t0,
    t_scan=t2 - t1,
)
np.savez(f"lab/solver2d/notes/scan_hw{hw:g}_V{V:g}.npz", **out)
print(f"build+ramp {t1 - t0:.2f}s scan {t2 - t1:.2f}s cce {np.array(cce)}")
