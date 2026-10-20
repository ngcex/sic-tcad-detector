"""Prototype: seed a devsim DD device from the batched engine instead of ramping.

ramp_bias (adaptive, 0.5..2 V steps, tiered tolerances) costs 0.6 s of run_field 1D,
2.3 s of run_transient and ~25 s of each 2D CCE(LET) device build. The engine solves
the same reverse-bias state (single quasi-Fermi Poisson, 1e-14 vs devsim) in ms.
Seed: set the cathode bias, copy Potential / Electrons / Holes from the engine, and
run ONE strict devsim DC solve (the same tier-0 solve the ramp ends with). Compare the
devsim state with the ramped one.

usage: uv run python lab/perf_profile/proto_seeded_ramp.py [bias_V ...]
"""

from __future__ import annotations

import sys
import time

import numpy as np

import devsim
from devsim.python_packages import simple_physics
from etna import DeviceConfig, set_devsim_quiet
from etna.api.device import build_device, device_kwargs
from etna.api.simulation import _engine_batch
from etna.core.batched.sweeps import _Ramp
from etna.core.devsim_reset import delete_device_and_mesh
from etna.core.drift_diffusion import ramp_bias, solve_tier

set_devsim_quiet(True)
cfg = DeviceConfig()


def state(info):
    d, r = info["device_name"], info["region_name"]
    return {
        m: np.array(devsim.get_node_model_values(device=d, region=r, name=m))
        for m in ("Potential", "Electrons", "Holes")
    }


for V in [float(a) for a in sys.argv[1:]] or [-10.0, -30.0, -100.0]:
    info = build_device(cfg)
    t = time.perf_counter()
    ramp_bias(info, V_target=-V, contact="cathode", V_step=0.5, V_step_max=2.0)
    t_ramp = time.perf_counter() - t
    ref = state(info)
    delete_device_and_mesh(info["device_name"])

    info = build_device(cfg)
    t = time.perf_counter()
    pb, x = _engine_batch([device_kwargs(cfg)])
    res = _Ramp(pb, 1e-8, 10.0).to(-V)
    d, r = info["device_name"], info["region_name"]
    devsim.set_parameter(
        device=d, name=simple_physics.GetContactBiasName("cathode"), value=-V
    )
    for m, arr in (
        ("Potential", res.psi[0]),
        ("Electrons", res.n[0]),
        ("Holes", res.p[0]),
    ):
        # devsim rejects zero carriers (the engine's Boltzmann n/p underflow in the
        # depleted epi): floor them; the potential is copied as is
        vals = arr if m == "Potential" else np.maximum(arr, 1e-300)
        devsim.set_node_values(device=d, region=r, name=m, values=vals.tolist())
    tier = solve_tier(0)
    t_seed = time.perf_counter() - t
    got = state(info)
    delete_device_and_mesh(d)
    dpsi = np.max(np.abs(got["Potential"] - ref["Potential"]))
    rel = {
        m: np.max(
            np.abs(got[m] - ref[m])
            / np.maximum(np.abs(ref[m]), 1e-30 + 1e-12 * np.abs(ref[m]).max())
        )
        for m in ("Electrons", "Holes")
    }
    print(
        f"V={V:6.1f}: ramp {t_ramp:.3f} s, seeded {t_seed:.3f} s (x{t_ramp / t_seed:.1f}), tier {tier}; "
        f"max|dpsi| {dpsi:.2e} V, max rel dn {rel['Electrons']:.2e}, dp {rel['Holes']:.2e}"
    )
