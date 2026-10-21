"""devsim reference numbers (and wall time) for the alternative structures of notebook 19.

Writes lab/geometry2d/notes/devsim_baselines.json (a few floats per structure) and
lab/geometry2d/notes/devsim_psi.npz (potential at the target bias, for field parity).
The fv2d parity bench (bench/parity.py) and tests/test_fv2d_parity.py read the JSON.

Usage: uv run python lab/geometry2d/bench/devsim_baselines.py [mesa gr elec]
Stack: ALT_STACK (etna.core.alternative_structures), recorded in the JSON.
"""

import json
import sys
import time
from pathlib import Path

import devsim
import numpy as np

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)

from etna.core.alternative_structures import (  # noqa: E402
    ALT_STACK,
    create_3d_electrode_device,
    create_guard_ring_device,
    create_mesa_device,
    guard_ring_node_sets,
    restore_cartesian_coords,
)
from etna.core.charge_collection_2d import (  # noqa: E402
    alpha_stripe_generation,
    cce_lateral_scan,
    cce_with_generation,
    contact_current_sharing,
)
from etna.core.drift_diffusion import ramp_bias, setup_sic_drift_diffusion  # noqa: E402
from etna.core.poisson import setup_poisson, solve_equilibrium  # noqa: E402

NOTES = Path(__file__).resolve().parents[1] / "notes"
OUT = NOTES / "devsim_baselines.json"
PSI = NOTES / "devsim_psi.npz"


def _psi(info):
    d, r = info["device_name"], info["region_name"]
    g = lambda n: np.array(
        devsim.get_node_model_values(device=d, region=r, name=n)
    )  # noqa: E731
    return g("x"), g("y"), g("Potential")


def _prepare(info, V, V_step_max=2.5):
    setup_poisson(info)
    solve_equilibrium(info)
    setup_sic_drift_diffusion(info)
    ramp_bias(info, V_target=V, contact="cathode", V_step=0.5, V_step_max=V_step_max)


def run_mesa(V=50.0, n=9):
    t0 = time.time()
    info = create_mesa_device(
        pillar_half_width_um=50.0,
        trench_width_um=10.0,
        trench_depth_um=10.0,
        depth_mesh="legacy",  # references stay on the legacy depth mesh (README)
    )
    _prepare(info, V)
    t_ramp = time.time() - t0
    scan = cce_lateral_scan(info, n_points=n, contact="cathode")
    both = cce_with_generation(
        info, alpha_stripe_generation(info, 0.0), contacts=("cathode", "anode")
    )
    t = time.time() - t0
    out = dict(
        V=V,
        nodes=info["num_nodes"],
        x_cm=scan["x_positions_cm"].tolist(),
        cce=scan["cce_values"].tolist(),
        cce_anode_center=both["anode"],
        t_ramp_s=t_ramp,
        t_total_s=t,
    )
    x, y, p = _psi(info)
    devsim.delete_device(device=info["device_name"])
    return out, (x, y, p)


def run_gr(V=10.0, n=11):
    t0 = time.time()
    info = create_guard_ring_device(
        sv_half_width_um=50.0,
        guard_ring_width_um=5.0,
        guard_ring_gap_um=3.0,
        guard_ring_depth_um=1.0,
        depth_mesh="legacy",
    )
    _prepare(info, V)
    t_ramp = time.time() - t0
    ns = guard_ring_node_sets(info)
    x_max = info["sv_half_width_cm"] + 1.5e-4
    scan = cce_lateral_scan(
        info, n_points=n, contact="anode", node_set=ns["anode"], x_max_cm=x_max
    )
    t_scan = time.time() - t0
    sharing = []
    for x_um in (48.0, 50.0, 51.5):
        frac, resid = contact_current_sharing(
            info, alpha_stripe_generation(info, x_um * 1e-4), node_sets=ns
        )
        sharing.append(dict(x_um=x_um, residual=resid, **frac))
    t = time.time() - t0
    out = dict(
        V=V,
        nodes=info["num_nodes"],
        x_cm=scan["x_positions_cm"].tolist(),
        cce_anode=scan["cce_values"].tolist(),
        sharing=sharing,
        t_ramp_s=t_ramp,
        t_scan_s=t_scan,
        t_total_s=t,
    )
    x, y, p = _psi(info)
    devsim.delete_device(device=info["device_name"])
    return out, (x, y, p)


def run_elec(V=50.0):
    t0 = time.time()
    info = create_3d_electrode_device(
        outer_radius_um=50.0, column_radius_um=5.0, top_gap_um=2.0
    )
    try:
        _prepare(info, V)
        t_ramp = time.time() - t0
        wall_inner = info["outer_radius_cm"] - info["wall_thickness_cm"]
        rs = np.linspace(info["sheet_inner_cm"] + 1e-4, wall_inner - 3 * 2e-4, 6)
        cath, an = [], []
        for r in rs:
            both = cce_with_generation(
                info, alpha_stripe_generation(info, r), contacts=("cathode", "anode")
            )
            cath.append(both["cathode"])
            an.append(both["anode"])
        t = time.time() - t0
        out = dict(
            V=V,
            nodes=info["num_nodes"],
            r_cm=rs.tolist(),
            cce_cathode=cath,
            cce_anode=an,
            t_ramp_s=t_ramp,
            t_total_s=t,
        )
        x, y, p = _psi(info)
        devsim.delete_device(device=info["device_name"])
    finally:
        restore_cartesian_coords()
    return out, (x, y, p)


if __name__ == "__main__":
    which = sys.argv[1:] or ["mesa", "gr", "elec"]
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    psi = dict(np.load(PSI)) if PSI.exists() else {}
    res["stack"] = {k: v for k, v in ALT_STACK.items()}
    for w in which:
        out, (x, y, p) = {"mesa": run_mesa, "gr": run_gr, "elec": run_elec}[w]()
        res[w] = out
        psi[f"{w}_x"], psi[f"{w}_y"], psi[f"{w}_psi"] = x, y, p
        print(
            w,
            json.dumps(
                {k: v for k, v in out.items() if k.startswith("t_") or k == "nodes"}
            ),
        )
    OUT.write_text(json.dumps(res, indent=1))
    np.savez_compressed(PSI, **psi)
