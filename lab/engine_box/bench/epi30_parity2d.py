"""Parity of the 2D engine vs devsim for epi thicknesses beyond the old 20 um bound.

Dark state: dark_state_2d vs create_2d_dd_device + ramp (Potential on every node, the
run_field path); lateral CCE: lateral_cce_scan_2d vs cce_lateral_scan. Device keywords
come from etna.api.device.device_kwargs(DeviceConfig(...)) exactly as the facades pass
them. The epi bound of ENGINE2D_BOX / ENGINE2D_DARK_BOX is lifted to 30 um for the run.

usage: uv run python lab/engine_box/bench/epi30_parity2d.py   (several minutes)
"""

from __future__ import annotations

import time
import uuid

import devsim
import numpy as np

from etna.api.device import DeviceConfig, device_kwargs, legacy_config
from etna.core.charge_collection_2d import cce_lateral_scan, create_2d_dd_device
from etna.core.devsim_lock import quiet_devsim
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.engine2d import scan as scan_mod
from etna.core.engine2d import dark_state_2d, lateral_cce_scan_2d

scan_mod.ENGINE2D_BOX["epi_thickness_cm"] = (5e-4, 30e-4)
scan_mod.ENGINE2D_DARK_BOX["epi_thickness_cm"] = (5e-4, 30e-4)

# (label, config, V_bias)
CASES = [
    (
        "default epi 22 hw50 50V",
        DeviceConfig(half_width_um=50.0, epi_thickness_um=22.0),
        50.0,
    ),
    (
        "default epi 25 hw25 20V",
        DeviceConfig(half_width_um=25.0, epi_thickness_um=25.0),
        20.0,
    ),
    (
        "default epi 30 hw25 60V",
        DeviceConfig(half_width_um=25.0, epi_thickness_um=30.0),
        60.0,
    ),
    (
        "default epi 30 hw25 100V",
        DeviceConfig(half_width_um=25.0, epi_thickness_um=30.0),
        100.0,
    ),
    (
        "no buffer n+ epi 30 hw25 60V",
        DeviceConfig(
            half_width_um=25.0,
            epi_thickness_um=30.0,
            buffer_thickness_um=None,
            buffer_end_doping=None,
        ),
        60.0,
    ),
    (
        "bulk 2e13 epi 30 hw150 80V",
        DeviceConfig(half_width_um=150.0, epi_thickness_um=30.0, N_D_bulk=2e13),
        80.0,
    ),
    (
        "legacy (no n+) epi 30 hw25 60V",
        legacy_config(half_width_um=25.0, epi_thickness_um=30.0),
        60.0,
    ),
]


def main():
    worst = dict(psi=0.0, cce=0.0)
    for label, cfg, vb in CASES:
        kw = device_kwargs(cfg)
        hw = cfg.half_width_um
        t = time.time()
        st = dark_state_2d(hw, vb, **kw)
        kw_scan = dict(kw)
        epi, nd = kw_scan.pop("epi_thickness_cm"), kw_scan.pop("N_D_bulk")
        e = lateral_cce_scan_2d(
            hw, vb, epi_thickness_cm=epi, N_D_bulk=nd, n_lateral_points=3, **kw_scan
        )
        te = time.time() - t
        t = time.time()
        name = f"e30_{uuid.uuid4().hex[:6]}"
        with quiet_devsim():
            info = create_2d_dd_device(
                half_width_um=hw, V_bias=vb, device_name=name, **kw
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
            f"RES {label:34s} dpsi {dpsi:.2e} V dCCE {dc:.2e} CCE {e['cce_values'][0]:.6f} "
            f"engine {te:.1f}s devsim {td:.1f}s",
            flush=True,
        )
    print("RES WORST", {k: f"{v:.2e}" for k, v in worst.items()})


if __name__ == "__main__":
    main()
