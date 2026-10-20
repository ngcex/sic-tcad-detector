"""Prototype: 1D run_field on the batched engine instead of a devsim DD ramp.

The engine's Poisson-with-quasi-Fermi potential equals devsim's reverse-bias DD
potential to ~1e-14 (lab/gpu_solver notes); the 1D field is the edge difference
quotient of it, the same as devsim's ElectricField edge model, interpolated onto
the nodes exactly like run_field does. Compares node arrays and timing.

usage: uv run python lab/perf_profile/proto_field1d_engine.py [bias_V ...]
"""

from __future__ import annotations

import sys
import time

import numpy as np

from etna import DeviceConfig, run_field, set_devsim_quiet
from etna.api.device import device_kwargs
from etna.api.simulation import _engine_batch
from etna.core.batched import reverse_sweep

set_devsim_quiet(True)


def field_engine(config, bias_V):
    kw = device_kwargs(config)
    pb, x = _engine_batch([kw])
    sw = reverse_sweep(
        pb,
        np.array([bias_V]),
        kw["substrate_thickness_cm"],
        kw["epi_thickness_cm"],
        keep_psi=True,
        capacitance=False,
    )
    if not sw.converged[-1, 0]:
        raise RuntimeError("engine did not converge")
    psi = sw.psi[-1, 0]
    x_c = 0.5 * (x[1:] + x[:-1])
    E_edge = (psi[:-1] - psi[1:]) / np.diff(
        x
    )  # devsim ElectricField = (V0 - V1) / EdgeLength
    return x, psi, np.abs(np.interp(x, x_c, E_edge))


cfg = DeviceConfig()
for V in [float(a) for a in sys.argv[1:]] or [-10.0, -30.0, -50.0, -100.0]:
    field_engine(cfg, V)  # warm
    t = time.perf_counter()
    ref = run_field(cfg, V)
    t_dev = time.perf_counter() - t
    ts = []
    for _ in range(3):
        t = time.perf_counter()
        x, psi, E = field_engine(cfg, V)
        ts.append(time.perf_counter() - t)
    xs = ref.mesh.x_coords
    same_mesh = xs.shape == x.shape and np.allclose(xs, x, rtol=0, atol=1e-16)
    dpsi = np.max(np.abs(psi - ref.metadata["potential"])) if same_mesh else np.nan
    Er = ref.y
    rel = np.max(np.abs(E - Er)) / np.max(Er) if same_mesh else np.nan
    print(
        f"V={V:6.1f}: devsim {t_dev:.3f} s, engine {min(ts):.4f} s (x{t_dev / min(ts):.0f}); same mesh {same_mesh}; "
        f"max|dpsi| {dpsi:.2e} V; max|dE|/max E {rel:.2e}"
    )
