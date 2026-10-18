"""dark_current_sweep ramp step: V_step 2-5 V vs the 0.5 V reference.

Default config (DeviceConfig()), run_dark_current's grid (0..-100 V, 20 points),
250-400 K. Prints max |dI/I| over the grid and the wall time per step size.

usage: uv run python lab/engine_box/bench/dark_vstep.py
"""

from __future__ import annotations

import time

import devsim
import numpy as np

from etna.api.device import DeviceConfig, device_kwargs
from etna.core.dark_current import create_dark_current_device, dark_current_sweep
from etna.core.devsim_lock import quiet_devsim
from etna.core.devsim_reset import reset_devsim_fully

BIAS = np.linspace(0.0, -100.0, 20)
STEPS = (0.5, 1.0, 2.0, 3.0, 5.0, 6.0)


def sweep(cfg, V_step):
    reset_devsim_fully()
    info = create_dark_current_device(N_t=cfg.N_t, E_t=cfg.E_t, **device_kwargs(cfg))
    try:
        t = time.perf_counter()
        r = dark_current_sweep(info, BIAS, area=cfg.area_cm2, V_step=V_step)
        return r, time.perf_counter() - t
    finally:
        devsim.delete_device(device=info["device_name"])


def main():
    for T in (250.0, 300.0, 350.0, 400.0):
        cfg = DeviceConfig(T=T)
        with quiet_devsim():
            ref, t_ref = sweep(cfg, 0.5)
        line = [f"T={T:.0f} ref 0.5 V: {t_ref:.2f}s pts {len(ref['voltages'])}"]
        for s in STEPS[1:]:
            with quiet_devsim():
                r, t = sweep(cfg, s)
            n = min(len(r["voltages"]), len(ref["voltages"]))
            rel = np.abs(r["I_total"][:n] / ref["I_total"][:n] - 1)
            # 0 V point: I ~ round-off, skip it in the relative error
            rel = rel[np.abs(r["voltages"][:n]) > 1e-9]
            line.append(
                f"{s:.0f}V: {t:.2f}s max|dI/I| {rel.max():.2e} pts {len(r['voltages'])}"
            )
        print(" | ".join(line), flush=True)


if __name__ == "__main__":
    main()
