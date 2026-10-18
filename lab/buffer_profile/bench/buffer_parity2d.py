"""Parity of the 2D engine (etna.core.engine2d) vs devsim WITH the graded n-/n+
transition and at the low N_D_bulk of the measured-C-V preset (box 5e13 -> 2e13).

Dark state: dark_state_2d vs create_2d_dd_device + ramp (Potential on every node);
lateral CCE: lateral_cce_scan_2d vs cce_lateral_scan. hw 25 um (the devsim 2D build is
tens of seconds per config).

usage: uv run python lab/buffer_profile/bench/buffer_parity2d.py
"""

from __future__ import annotations

import time
import uuid

import devsim
import numpy as np

from etna.core.charge_collection_2d import cce_lateral_scan, create_2d_dd_device
from etna.core.device import MEASURED_CV_STACK as M
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.engine2d import (
    dark_state_2d,
    engine2d_unsupported_reason,
    lateral_cce_scan_2d,
)

PRESET = dict(
    N_D_junction=M["N_D_junction"],
    L_transition=M["L_transition"],
    substrate_thickness_cm=M["substrate_thickness_cm"],
    nplus_thickness_cm=M["nplus_thickness_cm"],
    N_D_nplus=M["N_D_nplus"],
    buffer_thickness_cm=M["buffer_thickness_cm"],
    buffer_end_doping=M["buffer_end_doping"],
)
NOBUF = {k: v for k, v in PRESET.items() if not k.startswith("buffer")}
DEF = dict(substrate_thickness_cm=0.3e-4, nplus_thickness_cm=1e-4, N_D_nplus=1e19)

# (label, hw, V_bias, epi, N_D_bulk, extra kwargs)
CASES = [
    ("preset hw25 20V", 25.0, 20.0, 10e-4, M["N_D_bulk"], PRESET),
    ("preset hw25 60V", 25.0, 60.0, 10e-4, M["N_D_bulk"], PRESET),
    ("preset hw25 100V", 25.0, 100.0, 10e-4, M["N_D_bulk"], PRESET),
    ("preset doping no buffer 40V", 25.0, 40.0, 10e-4, M["N_D_bulk"], NOBUF),
    ("default bulk 2e13 40V", 25.0, 40.0, 10e-4, 2e13, DEF),
    (
        "default + buffer 0.5um 1e16 60V",
        25.0,
        60.0,
        10e-4,
        8.252732e13,
        dict(DEF, buffer_thickness_cm=0.5e-4, buffer_end_doping=1e16),
    ),
    (
        "default + buffer 2.5um 1e19 60V",
        25.0,
        60.0,
        10e-4,
        8.252732e13,
        dict(DEF, buffer_thickness_cm=2.5e-4, buffer_end_doping=1e19),
    ),
    (
        "bulk 2e13 epi 5um buffer 0.5um 1e16 20V",
        25.0,
        20.0,
        5e-4,
        2e13,
        dict(DEF, buffer_thickness_cm=0.5e-4, buffer_end_doping=1e16),
    ),
    (
        "bulk 2e13 epi 20um buffer 2.5um 1e19 100V",
        25.0,
        100.0,
        20e-4,
        2e13,
        dict(DEF, buffer_thickness_cm=2.5e-4, buffer_end_doping=1e19),
    ),
]


def main():
    worst = dict(psi=0.0, cce=0.0)
    for label, hw, vb, epi, nd, extra in CASES:
        # no-buffer cases explicit (builder default = the stack's transition since 2026-10-18)
        extra = {"buffer_thickness_cm": None, "buffer_end_doping": None, **extra}
        why = engine2d_unsupported_reason(hw, vb, epi, nd, 3, **extra)
        assert why is None, (label, why)
        t = time.time()
        e = lateral_cce_scan_2d(
            hw, vb, epi_thickness_cm=epi, N_D_bulk=nd, n_lateral_points=3, **extra
        )
        st = dark_state_2d(hw, vb, epi_thickness_cm=epi, N_D_bulk=nd, **extra)
        te = time.time() - t
        t = time.time()
        name = f"b2d_{uuid.uuid4().hex[:6]}"
        info = create_2d_dd_device(
            half_width_um=hw,
            V_bias=vb,
            device_name=name,
            epi_thickness_cm=epi,
            N_D_bulk=nd,
            **extra,
        )
        psi_d = np.array(
            devsim.get_node_model_values(
                device=name, region=info["region_name"], name="Potential"
            )
        )
        d = cce_lateral_scan(info, n_points=3)
        devsim.delete_device(device=name)
        reset_devsim_fully(preserve_solver=True)
        td = time.time() - t
        dpsi = float(np.max(np.abs(st["Potential"] - psi_d)))
        dc = float(
            np.max(np.abs(np.asarray(e["cce_values"]) - np.asarray(d["cce_values"])))
        )
        worst["psi"], worst["cce"] = max(worst["psi"], dpsi), max(worst["cce"], dc)
        print(
            f"RES {label:44s} dpsi {dpsi:.2e} V dCCE {dc:.2e} CCE {e['cce_values'][0]:.6f} "
            f"engine {te:.1f}s devsim {td:.1f}s",
            flush=True,
        )
    print("RES WORST", {k: f"{v:.2e}" for k, v in worst.items()})


if __name__ == "__main__":
    main()
