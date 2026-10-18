"""2D CCE normalisation: the never-injected contact-node generation (task 2 of the flip).

Prints every 2D CCE number that the contact-node convention can move, devsim and engine,
plus the contact-node share of the box-rule generation integral. Run before and after
the change and diff.

usage: uv run python lab/default_flip/bench/cce2d_contacts.py OUT.json
"""

from __future__ import annotations

import json
import sys
import uuid

import devsim
import numpy as np

from etna.core.charge_collection import add_generation_to_dd
from etna.core.charge_collection_2d import (
    _robust_dc_solve,
    alpha_stripe_generation,
    cce_lateral_scan,
    compute_cce_2d,
    create_2d_dd_device,
    integrate_nodes_box_2d,
)
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.drift_diffusion import _contact_node_indices, contacts_on_region
from etna.core.engine2d import lateral_cce_scan_2d

LEGACY = dict(
    substrate_thickness_cm=1e-4,
    nplus_thickness_cm=None,
    buffer_thickness_cm=None,
    buffer_end_doping=None,
    N_D_junction=2.958237e15,
    L_transition=0.987e-4,
    N_D_bulk=8.252732e13,
)
CASES = {
    "default hw25 20V": (25.0, 20.0, {}),
    "default hw25 50V": (25.0, 50.0, {}),
    "default hw25 50V nplus_generation": (25.0, 50.0, {"nplus_generation": True}),
    "legacy planar hw25 50V": (25.0, 50.0, LEGACY),
    "legacy planar hw50 20V": (50.0, 20.0, LEGACY),
}


def contact_share(info, G):
    dev, reg = info["device_name"], info["region_name"]
    Gc = np.array(G, dtype=float, copy=True)
    mask = np.zeros(Gc.size, bool)
    for c in contacts_on_region(dev, reg):
        mask[sorted(_contact_node_indices(dev, reg, c))] = True
    tot = integrate_nodes_box_2d(info, Gc)
    Gc[~mask] = 0.0
    return integrate_nodes_box_2d(info, Gc) / tot if tot > 0 else 0.0


def main():
    out = {}
    for label, (hw, vb, kw) in CASES.items():
        reset_devsim_fully()
        name = f"c2d_{uuid.uuid4().hex[:6]}"
        info = create_2d_dd_device(half_width_um=hw, V_bias=vb, device_name=name, **kw)
        G = alpha_stripe_generation(info, 0.0, zero_on_contacts=False)
        share = contact_share(info, G)
        d = cce_lateral_scan(info, n_points=3)
        # v3 recipe (uniform generation across epi + n+, contact nodes included)
        y = np.array(devsim.get_node_model_values(device=name, region="sic", name="y"))
        gu = np.where(y >= info["substrate_thickness_cm"], 1e18, 0.0)
        add_generation_to_dd(info, gu)
        _robust_dc_solve()
        v3 = float(compute_cce_2d(info, gu))
        ushare = contact_share(info, gu)
        devsim.delete_device(device=name)
        reset_devsim_fully()
        try:
            e = lateral_cce_scan_2d(hw, vb, n_lateral_points=3, **kw)["cce_values"]
            e = [float(v) for v in e]
        except Exception as exc:  # noqa: BLE001
            e = f"{type(exc).__name__}: {exc}"
        out[label] = dict(
            devsim=[float(v) for v in d["cce_values"]],
            engine=e,
            stripe_contact_share=share,
            uniform_cce=v3,
            uniform_contact_share=ushare,
        )
        print("RES", label, json.dumps(out[label]), flush=True)
    json.dump(out, open(sys.argv[1], "w"), indent=1)


if __name__ == "__main__":
    main()
