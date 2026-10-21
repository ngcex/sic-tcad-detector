"""fv2d vs devsim on the notebook-19 structures (same devsim mesh, imported verbatim).

Reads lab/geometry2d/notes/devsim_baselines.json / devsim_psi.npz (devsim_baselines.py)
and prints, per structure: wall times (mesh import, ramp, scan), max |dCCE|, max |dpsi|,
the conservation residual and the dark-check result.

Usage: uv run python lab/geometry2d/bench/parity.py [mesa gr elec]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)

from etna.core.alternative_structures import (  # noqa: E402
    create_3d_electrode_device,
    create_guard_ring_device,
    create_mesa_device,
)
from etna.core.fv2d import (  # noqa: E402
    alpha_stripe_generation_fv2d,
    cce_lateral_scan_fv2d,
    contact_sharing_fv2d,
    cut_node_sets,
    mesh_from_builder,
    prepared_solver,
)

NOTES = Path(__file__).resolve().parents[1] / "notes"
REF = json.loads((NOTES / "devsim_baselines.json").read_text())
PSI = np.load(NOTES / "devsim_psi.npz")


def psi_diff(solver, tag):
    m = solver.mesh
    x, y, p = PSI[f"{tag}_x"], PSI[f"{tag}_y"], PSI[f"{tag}_psi"]
    i = np.searchsorted(m.xs, x)
    j = np.searchsorted(m.ys, y)
    return float(np.max(np.abs(solver.potential()[j, i] - p)))


def run(tag):
    t0 = time.perf_counter()
    if tag == "mesa":
        mesh, _ = mesh_from_builder(
            create_mesa_device,
            pillar_half_width_um=50.0,
            trench_width_um=10.0,
            trench_depth_um=10.0,
        )
    elif tag == "gr":
        mesh, _ = mesh_from_builder(
            create_guard_ring_device,
            sv_half_width_um=50.0,
            guard_ring_width_um=5.0,
            guard_ring_gap_um=3.0,
            guard_ring_depth_um=1.0,
        )
    else:
        mesh, _ = mesh_from_builder(
            create_3d_electrode_device,
            outer_radius_um=50.0,
            column_radius_um=5.0,
            top_gap_um=2.0,
        )
    t_mesh = time.perf_counter() - t0
    ref = REF[tag]
    s = prepared_solver(mesh, {"cathode": ref["V"]})
    t_ramp = time.perf_counter() - t0 - t_mesh
    out = dict(
        nodes=mesh.N,
        grid=(mesh.nx, mesh.ny),
        t_mesh=t_mesh,
        t_ramp=t_ramp,
        ramp_steps=len(s.ramp_steps),
        newton=s.newton_iterations,
        dark_check=s.dark_check_result,
        dpsi=psi_diff(s, tag),
    )
    t1 = time.perf_counter()
    if tag == "mesa":
        sc = cce_lateral_scan_fv2d(s, n_points=len(ref["cce"]), contact="cathode")
        out["dcce"] = float(np.max(np.abs(sc["cce_values"] - ref["cce"])))
    elif tag == "gr":
        info = mesh.info
        ns = cut_node_sets(mesh, info["sv_half_width_cm"] + 0.5 * info["gap_cm"])
        sc = cce_lateral_scan_fv2d(
            s,
            n_points=len(ref["cce_anode"]),
            contact="anode",
            node_sets=ns,
            x_max_cm=info["sv_half_width_cm"] + 1.5e-4,
        )
        out["dcce"] = float(np.max(np.abs(sc["cce_values"] - ref["cce_anode"])))
        ds = []
        for row in ref["sharing"]:
            fr, res = contact_sharing_fv2d(
                s, alpha_stripe_generation_fv2d(mesh, row["x_um"] * 1e-4), ns
            )
            ds.append(max(abs(fr[c] - row[c]) for c in fr))
        out["dsharing"] = max(ds)
    else:
        sc = cce_lateral_scan_fv2d(s, contact="cathode", x_positions_cm=ref["r_cm"])
        out["dcce"] = float(np.max(np.abs(sc["cce_values"] - ref["cce_cathode"])))
        sa = cce_lateral_scan_fv2d(s, contact="anode", x_positions_cm=ref["r_cm"])
        out["dcce_anode"] = float(np.max(np.abs(sa["cce_values"] - ref["cce_anode"])))
    out["t_scan"] = time.perf_counter() - t1
    out["max_residual"] = float(np.max(np.abs(sc["residuals"])))
    out["max_direct_mismatch"] = float(np.max(sc["direct_mismatch"]))
    out["cce"] = sc["cce_values"].tolist()
    out["t_total"] = time.perf_counter() - t0
    out["devsim_t_total"] = ref["t_total_s"]
    return out


if __name__ == "__main__":
    for tag in sys.argv[1:] or ["mesa", "gr", "elec"]:
        print(tag, json.dumps(run(tag), default=str))
