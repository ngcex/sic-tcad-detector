"""Adaptive tensor-line refinement (fv2d, native meshes): CCE mesh convergence on the
guard ring (10 V, anode scan, as lab/edge_cce) and the 3D electrode (50 V, radial
cathode CCE), against the devsim-mesh value (devsim_baselines.json).

Usage: uv run python lab/geometry2d/bench/amr.py [gr elec]
"""

import json
import sys
from pathlib import Path

import numpy as np

from etna.core.fv2d import cut_node_sets
from etna.core.fv2d.refine import adaptive_cce_scan
from etna.core.fv2d.structures import electrode3d_structure, guard_ring_structure

NOTES = Path(__file__).resolve().parents[1] / "notes"
REF = json.loads((NOTES / "devsim_baselines.json").read_text())


def gr(tol=1e-5):
    st = guard_ring_structure()
    info = st.info
    xcut = info["sv_half_width_cm"] + 0.5 * info["gap_cm"]
    ref = REF["gr"]
    kw = dict(contact="anode", x_positions_cm=ref["x_cm"],
              node_sets_fn=lambda m: cut_node_sets(m, xcut))  # fmt: skip
    return adaptive_cce_scan(st, {"cathode": ref["V"]}, kw, tol=tol), np.array(
        ref["cce_anode"]
    )


def elec(tol=1e-5):
    st = electrode3d_structure()
    ref = REF["elec"]
    kw = dict(contact="cathode", x_positions_cm=ref["r_cm"])
    return adaptive_cce_scan(st, {"cathode": ref["V"]}, kw, tol=tol), np.array(
        ref["cce_cathode"]
    )


if __name__ == "__main__":
    for tag in sys.argv[1:] or ["gr", "elec"]:
        res, ref = {"gr": gr, "elec": elec}[tag]()
        for L in res["levels"]:
            print(
                tag,
                f"level {L['level']}: {L['nx']}x{L['ny']} = {L['nodes']} nodes, "
                f"{L['seconds']:.2f} s, dCCE vs previous {L.get('dcce', float('nan')):.2e}, "
                f"max |CCE - devsim mesh| {np.max(np.abs(L['cce'] - ref)):.2e}, "
                f"new lines x {L.get('n_new_x', 0)} y {L.get('n_new_y', 0)}",
            )
        print(
            tag, "error estimate", res["error_estimate"], "converged", res["converged"]
        )
        print(tag, "final cce", np.round(res["scan"]["cce_values"], 7).tolist())
