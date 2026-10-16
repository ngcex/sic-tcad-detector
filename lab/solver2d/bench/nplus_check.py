"""n+ layer in the production devsim builders (nplus_thickness_cm option): wall, W, CCE, dark current.

Compares with the monkeypatched spike (lab/solver2d/bench/nplus.py, LOG.md "n+ substrate spike").
Device = create_dd_device graded defaults (N_D_junction 2.9e15, N_D_bulk 8.5e13, L 1 um).

usage: uv run python lab/solver2d/bench/nplus_check.py TASK [NPLUS_UM ...]
  TASK = wall | cv | cce | cce_alpha | dark | wall2d
"""

from __future__ import annotations

import sys
import time
import uuid

import devsim
import numpy as np


def kw(t_um, gen=False):
    if t_um is None:
        return {}
    return dict(nplus_thickness_cm=t_um * 1e-4, N_D_nplus=1e19, nplus_generation=gen)


def wall(t_um):
    from etna.core.devsim_reset import reset_devsim_fully
    from etna.core.drift_diffusion import create_dd_device, ramp_bias

    info = create_dd_device(device_name=f"w{uuid.uuid4().hex[:6]}", **kw(t_um))
    V = 0.0
    try:
        while V < 150.0 - 1e-9:
            ramp_bias(info, V_target=V + 0.5, contact="cathode", V_step=0.5)
            V += 0.5
    except RuntimeError:
        pass
    reset_devsim_fully(preserve_solver=True)
    return f"last converged {V:.1f} V (0.5 V fixed steps)"


def cv(t_um):
    from etna.core.cv_analysis import cv_sweep
    from etna.core.drift_diffusion import create_dd_device

    info = create_dd_device(device_name=f"c{uuid.uuid4().hex[:6]}", **kw(t_um))
    r = cv_sweep(info, -np.array([0.0, 5.0, 10.0, 30.0, 50.0]), area=0.05)
    devsim.delete_device(device=info["device_name"])
    return "W_um " + np.array2string(
        np.asarray(r["depletion_widths"]) * 1e4, precision=6
    )


def cce(t_um, gen=False):
    from etna.core.charge_collection import cce_vs_bias

    r = cce_vs_bias(-np.arange(10.0, 55.0, 10.0), device_kwargs=kw(t_um, gen))
    raw = np.asarray(r["I_collected"]) / r["I_generated"]
    return (
        "CCE "
        + np.array2string(r["cce_values"], precision=6)
        + " raw "
        + np.array2string(raw, precision=6)
        + f" I_gen {r['I_generated']:.6e} truncated {r['truncated']}"
    )


def dark(t_um):
    from etna.core.dark_current import create_dark_current_device, dark_current_sweep

    info = create_dark_current_device(
        N_t=2.8e12, device_name=f"d{uuid.uuid4().hex[:6]}", **kw(t_um)
    )
    r = dark_current_sweep(info, np.array([-10.0, -30.0, -60.0]), area=0.04)
    devsim.delete_device(device=info["device_name"])
    return "I_A " + np.array2string(np.abs(r["I_total"]), precision=7)


def wall2d(t_um, hw=25.0, epi_um=10.0, nd=5e13, v_end=150.0):
    from etna.core.charge_collection_2d import create_2d_dd_device
    from etna.core.devsim_reset import reset_devsim_fully
    from etna.core.drift_diffusion import ramp_bias

    info = create_2d_dd_device(
        half_width_um=hw,
        V_bias=0.0,
        epi_thickness_cm=epi_um * 1e-4,
        N_D_bulk=nd,
        **kw(t_um),
    )
    V = 0.0
    try:
        while V < v_end - 1e-9:
            ramp_bias(info, V_target=V + 0.5, contact="cathode", V_step=0.5)
            V += 0.5
    except RuntimeError:
        pass
    reset_devsim_fully(preserve_solver=True)
    return f"2D hw={hw} epi={epi_um} N_D_bulk={nd:.1e}: last converged {V:.1f} V"


if __name__ == "__main__":
    task = sys.argv[1]
    for arg in sys.argv[2:] or ["none"]:
        t_um = None if arg == "none" else float(arg)
        t0 = time.perf_counter()
        if task == "cce_alpha":
            out = cce(t_um, gen=True)
        else:
            out = globals()[task](t_um)
        print(
            f"RES {task} nplus={arg}um: {out} ({time.perf_counter() - t0:.1f}s)",
            flush=True,
        )
