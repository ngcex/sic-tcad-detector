"""Prototype: build the transient device once, restore the dark state per dose rate."""

import time
import uuid

import devsim
import numpy as np

from etna.api.device import DeviceConfig, device_kwargs
from etna.core import transient as tr
from etna.core.devsim_lock import quiet_devsim
from etna.core.drift_diffusion import create_dd_device, ramp_bias

cfg = DeviceConfig()
kw = device_kwargs(cfg)
rates = np.array([20.0, 50.0, 100.0, 150.0, 200.0, 230.0])

with quiet_devsim():
    t = time.perf_counter()
    df = tr.transient_cce_vs_dose_rate(
        -30.0, rates, epi_thickness_cm=kw.pop("epi_thickness_cm"), device_kwargs=kw
    )
    t_ref = time.perf_counter() - t

    t = time.perf_counter()
    kw2 = device_kwargs(cfg)
    info = create_dd_device(device_name=f"tr_{uuid.uuid4().hex[:6]}", **kw2)
    dev, reg = info["device_name"], info["region_name"]
    tr.add_auger_recombination(info)
    ramp_bias(info, -30.0, contact="anode")
    t_build = time.perf_counter() - t
    names = ("Potential", "Electrons", "Holes")
    snap = {
        n: list(devsim.get_node_model_values(device=dev, region=reg, name=n))
        for n in names
    }
    x = np.array(devsim.get_node_model_values(device=dev, region=reg, name="x"))
    out = []
    for r in rates:
        for n in names:
            devsim.set_node_values(device=dev, region=reg, name=n, values=snap[n])
        G = tr.proton_generation_profile(x, E_MeV=62, dose_rate_Gy_s=r)
        G[x < info["junction_pos"]] = 0.0
        G = tr.confine_generation_to_epi(info, x, G)
        s = tr.TransientSolver(info, contact=df.attrs["contact"])
        s.initialize()
        res = s.simulate_pulse(G, dose_rate_Gy_s=r)
        out.append(s.compute_transient_cce(res, G, x))
    devsim.delete_device(device=dev)
    t_once = time.perf_counter() - t

print(
    "RES ref",
    round(t_ref, 2),
    "s; build-once",
    round(t_once, 2),
    "s (build+ramp",
    round(t_build, 2),
    "s)",
)
print("RES max|dCCE|", np.max(np.abs(np.array(out) - df["transient_cce"].to_numpy())))
