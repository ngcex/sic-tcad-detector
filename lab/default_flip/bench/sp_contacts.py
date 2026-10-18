"""Single-particle CCE vs steady-state 2D CCE, legacy diode and default stack (task 2 check).

usage: uv run python lab/default_flip/bench/sp_contacts.py
"""

from __future__ import annotations

import devsim

import etna.core.single_particle as sp
from etna.core.charge_collection import add_generation_to_dd
from etna.core.charge_collection_2d import (
    _robust_dc_solve,
    compute_cce_2d,
    create_2d_dd_device,
)
from etna.core.devsim_reset import reset_devsim_fully

LEGACY_2D = dict(
    substrate_thickness_cm=1e-4,
    nplus_thickness_cm=None,
    buffer_thickness_cm=None,
    buffer_end_doping=None,
    N_D_junction=2.931417e15,
    N_D_bulk=8.823059e13,
    L_transition=9.872267e-5,
)


def run(label, V, kw):
    reset_devsim_fully()
    info = create_2d_dd_device(half_width_um=25.0, V_bias=V, device_name="spc", **kw)
    try:
        g, q = sp.ion_track_generation_2d(info, LET_keV_um=1.0)
        G = g * 1e3
        add_generation_to_dd(info, G)
        _robust_dc_solve()
        cce_ss = compute_cce_2d(info, G, contact="cathode")
        # contact None: cathode on the legacy diode, anode with n+ (the direct n+ cathode
        # current is round-off dominated, simulate_single_particle docstring)
        r = sp.simulate_single_particle(info, g, contact=None)
        cce_sp = r["Q_collected"] / q
        print(
            f"RES {label}: single-particle {cce_sp:.6f}  steady-state {cce_ss:.6f}  diff {cce_sp - cce_ss:+.2e}",
            flush=True,
        )
    finally:
        devsim.delete_device(device="spc")


if __name__ == "__main__":
    run("legacy hw25 10V", 10.0, LEGACY_2D)
    run("default hw25 50V", 50.0, {})
