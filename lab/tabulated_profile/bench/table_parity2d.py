"""Parity of the 2D engine (etna.core.engine2d) vs devsim with a TABULATED epi profile.

Dark state: dark_state_2d vs create_2d_dd_device + ramp (Potential on every node);
lateral CCE: lateral_cce_scan_2d vs cce_lateral_scan (3 points). hw 25 um (the devsim 2D
build is tens of seconds per config). Cases: the 19-knot measured-C-V table at 20 / 60 /
100 V, the table scaled by -0.3 / +0.3 dex (ENGINE2D_TABLE_BOX corners) at 40 V, a 3-knot
table at 60 V, and the table stretched x2 on a 20 um epi (n+ 5 um) at 100 V.

usage: uv run python lab/tabulated_profile/bench/table_parity2d.py
"""

from __future__ import annotations

import time
import uuid

import devsim
import numpy as np

from etna.core.charge_collection_2d import cce_lateral_scan, create_2d_dd_device
from etna.core.device import MEASURED_CV_TABLE, MEASURED_CV_TABULATED_STACK as S
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.engine2d import (
    dark_state_2d,
    engine2d_unsupported_reason,
    lateral_cce_scan_2d,
)


def scaled(table, dex):
    return tuple((d, n * 10.0**dex) for d, n in table)


BASE = dict(
    substrate_thickness_cm=S["substrate_thickness_cm"],
    nplus_thickness_cm=S["nplus_thickness_cm"],
    N_D_nplus=S["N_D_nplus"],
    buffer_thickness_cm=None,
    buffer_end_doping=None,
)
T20 = tuple((2 * d, n) for d, n in MEASURED_CV_TABLE)
# (label, hw, V_bias, epi, extra kwargs)
CASES = [
    ("table hw25 20V", 25.0, 20.0, 10e-4, dict(BASE, donor_table=MEASURED_CV_TABLE)),
    ("table hw25 60V", 25.0, 60.0, 10e-4, dict(BASE, donor_table=MEASURED_CV_TABLE)),
    ("table hw25 100V", 25.0, 100.0, 10e-4, dict(BASE, donor_table=MEASURED_CV_TABLE)),
    (
        "table -0.3 dex 40V",
        25.0,
        40.0,
        10e-4,
        dict(BASE, donor_table=scaled(MEASURED_CV_TABLE, -0.3)),
    ),
    (
        "table +0.3 dex 40V",
        25.0,
        40.0,
        10e-4,
        dict(BASE, donor_table=scaled(MEASURED_CV_TABLE, 0.3)),
    ),
    (
        "3-knot table 60V",
        25.0,
        60.0,
        10e-4,
        dict(BASE, donor_table=((0.0, 1e15), (5.0, 5e13), (10.0, 1e16))),
    ),
    (
        "table x2 depth, epi 20um, n+ 5um 100V",
        25.0,
        100.0,
        20e-4,
        dict(BASE, donor_table=T20, nplus_thickness_cm=5e-4),
    ),
]


def main():
    worst = dict(psi=0.0, cce=0.0)
    for label, hw, vb, epi, extra in CASES:
        why = engine2d_unsupported_reason(hw, vb, epi, None, 3, **extra)
        assert why is None, (label, why)
        t = time.time()
        e = lateral_cce_scan_2d(
            hw, vb, epi_thickness_cm=epi, n_lateral_points=3, **extra
        )
        st = dark_state_2d(hw, vb, epi_thickness_cm=epi, **extra)
        te = time.time() - t
        t = time.time()
        name = f"t2d_{uuid.uuid4().hex[:6]}"
        info = create_2d_dd_device(
            half_width_um=hw,
            V_bias=vb,
            device_name=name,
            epi_thickness_cm=epi,
            **extra,
        )
        psi_d = np.array(
            devsim.get_node_model_values(
                device=name, region=info["region_name"], name="Potential"
            )
        )
        don_d = np.array(
            devsim.get_node_model_values(
                device=name, region=info["region_name"], name="Donors"
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
            f"RES {label:40s} dpsi {dpsi:.2e} V dCCE {dc:.2e} CCE {e['cce_values'][0]:.6f} "
            f"donors[max] {don_d.max():.3g} engine {te:.1f}s devsim {td:.1f}s",
            flush=True,
        )
    print("RES WORST", {k: f"{v:.2e}" for k, v in worst.items()})


if __name__ == "__main__":
    main()
