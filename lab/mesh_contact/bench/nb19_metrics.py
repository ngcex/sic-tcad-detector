"""Notebook-19 metrics on the legacy and the graded depth mesh (fv2d on devsim's own
meshes, imported verbatim: parity with devsim 1e-15, so these ARE the devsim numbers of
each mesh, without devsim's 30-110 s guard-ring ramp).

Settings as notebook 19 (scripts/create_notebook_19.py): SV half-width 50 um, 50 V
(V_BIAS), 10 V (V_BIAS_REFERENCE), N_SCAN 9, stripe sigma 2 um; planar = the legacy
stack 2D diode (create_sic_2d_device, nplus None); guard ring 5 um / gap 3 um / depth 1
um, anode read-out inside the SV and edge-event sharing at x = 48 / 50 / 51.5 um;
overhang reference (no ring) at 10 V. The 3D electrode has no contact on the epi (its
depth mesh is unchanged) and the telescope has no fv2d twin (lab/mesh_contact
bench/telescope_devsim.py).

Usage: uv run python lab/mesh_contact/bench/nb19_metrics.py
"""

import json
import time
from pathlib import Path

import numpy as np

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.alternative_structures import (  # noqa: E402
    create_guard_ring_device,
    create_mesa_device,
)
from etna.core.device2d import create_sic_2d_device  # noqa: E402
from etna.core.fv2d import (  # noqa: E402
    alpha_stripe_generation_fv2d,
    cce_lateral_scan_fv2d,
    contact_sharing_fv2d,
    cut_node_sets,
    mesh_from_builder,
    prepared_solver,
)

SV = 50.0
V_BIAS, V_REF, N_SCAN, SIG = 50.0, 10.0, 9, 2e-4
LEGACY_STACK = dict(
    substrate_thickness_cm=1e-4,
    nplus_thickness_cm=None,
    buffer_thickness_cm=None,
    buffer_end_doping=None,
    N_D_junction=2.958237e15,
    N_D_bulk=8.252732e13,
    L_transition=0.987e-4,
)
GR_KW = dict(
    sv_half_width_um=SV,
    guard_ring_width_um=5.0,
    guard_ring_gap_um=3.0,
    guard_ring_depth_um=1.0,
)
EDGE_UM = (SV - 2.0, SV, SV + 1.5)


def sharing(s, mesh, ns):
    rows = {}
    for x in EDGE_UM:
        frac, res = contact_sharing_fv2d(
            s, alpha_stripe_generation_fv2d(mesh, x * 1e-4), ns
        )
        rows[x] = dict(
            anode=-frac["anode"],
            guard=-frac.get("guard_ring_anode", 0.0),
            lost=1.0 - frac["cathode"],
            residual=res,
        )
    return rows


def run(depth_mesh):
    out = {}
    t = time.perf_counter()
    m, _ = mesh_from_builder(
        create_sic_2d_device, half_width_um=SV, depth_mesh=depth_mesh, **LEGACY_STACK
    )
    for V in (V_REF, V_BIAS):
        sc = cce_lateral_scan_fv2d(
            prepared_solver(m, {"cathode": V}, dark_check=False), n_points=N_SCAN
        )
        out[f"planar_{V:g}V"] = sc["cce_values"].tolist()
    m, _ = mesh_from_builder(
        create_mesa_device,
        pillar_half_width_um=SV,
        trench_width_um=10.0,
        trench_depth_um=10.0,
        depth_mesh=depth_mesh,
    )
    for V in (V_REF, V_BIAS):
        sc = cce_lateral_scan_fv2d(
            prepared_solver(m, {"cathode": V}, dark_check=False), n_points=N_SCAN
        )
        out[f"mesa_{V:g}V"] = sc["cce_values"].tolist()
    m, _ = mesh_from_builder(create_guard_ring_device, depth_mesh=depth_mesh, **GR_KW)
    info = m.info
    ns = cut_node_sets(m, info["sv_half_width_cm"] + 0.5 * info["gap_cm"])
    for V in (V_REF, V_BIAS):
        s = prepared_solver(m, {"cathode": V}, dark_check=False)
        sc = cce_lateral_scan_fv2d(
            s,
            n_points=N_SCAN,
            contact="anode",
            node_sets=ns,
            x_max_cm=info["sv_half_width_cm"] - 3 * SIG,
        )
        out[f"gr_anode_{V:g}V"] = sc["cce_values"].tolist()
        # lab/edge_cce map (11 points to mid-gap, 10 V) and its f ratios
        if V == V_REF:
            e = cce_lateral_scan_fv2d(
                s,
                n_points=11,
                contact="anode",
                node_sets=ns,
                x_max_cm=info["sv_half_width_cm"] + 1.5e-4,
            )
            out["edge_map_10V"] = e["cce_values"].tolist()
            out["edge_map_x_um"] = e["x_positions_um"].tolist()
        out[f"gr_sharing_{V:g}V"] = {str(k): v for k, v in sharing(s, m, ns).items()}
    m, _ = mesh_from_builder(
        create_guard_ring_device, with_guard=False, depth_mesh=depth_mesh, **GR_KW
    )
    s = prepared_solver(m, {"cathode": V_REF}, dark_check=False)
    out["overhang_sharing_10V"] = {str(k): v for k, v in sharing(s, m, None).items()}
    out["seconds"] = time.perf_counter() - t
    return out


if __name__ == "__main__":
    res = {dm: run(dm) for dm in ("legacy", "graded")}
    Path(__file__).with_name("nb19_metrics.json").write_text(json.dumps(res, indent=1))
    L, G = res["legacy"], res["graded"]
    for k in L:
        if k.startswith(("planar", "mesa", "gr_anode")):
            a, b = np.array(L[k]), np.array(G[k])
            print(f"{k:16s} centre {a[0]:.6f} -> {b[0]:.6f} ({b[0] - a[0]:+.2e})  "
                  f"edge {a[-1]:.6f} -> {b[-1]:.6f}  edge/ctr {a[-1] / a[0]:.6f} -> "
                  f"{b[-1] / b[0]:.6f}")  # fmt: skip
    a, b = np.array(L["edge_map_10V"]), np.array(G["edge_map_10V"])
    print("edge map 10 V  x_um:", np.round(L["edge_map_x_um"], 2).tolist())
    print("   f legacy:", np.round(a / a[0], 6).tolist())
    print("   f graded:", np.round(b / b[0], 6).tolist())
    for k in ("gr_sharing_10V", "gr_sharing_50V", "overhang_sharing_10V"):
        for x in L[k]:
            la, ga = L[k][x], G[k][x]
            print(f"{k:21s} x {float(x):5.1f}: anode {la['anode']:.4f}->{ga['anode']:.4f} "
                  f"guard {la['guard']:.4f}->{ga['guard']:.4f} lost {la['lost']:.4f}->"
                  f"{ga['lost']:.4f} resid {ga['residual']:+.0e}")  # fmt: skip
    print(f"time legacy {L['seconds']:.1f} s, graded {G['seconds']:.1f} s")
