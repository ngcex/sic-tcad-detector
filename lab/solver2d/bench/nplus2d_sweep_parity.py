"""Production-path parity WITH n+: _sweep_one(solver="auto") vs ("devsim"), and
run_field engine vs forced devsim path.

usage:
  uv run python lab/solver2d/bench/nplus2d_sweep_parity.py sweep START STOP
      -> appends to lab/solver2d/notes/nplus2d_sweep_parity.csv
  uv run python lab/solver2d/bench/nplus2d_sweep_parity.py field
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path
from unittest import mock

import numpy as np

D = 8.823059e13
# (hw, epi um, N_D_bulk, V, n+ um, N_D_nplus, nplus_generation)
POINTS = [
    (25, 10, D, 5, 1, 1e19, False),
    (25, 10, 5e13, 100, 1, 1e19, False),
    (25, 5, 5e14, 100, 5, 1e18, False),
    (25, 20, 5e14, 20, 1, 1e18, False),
    (25, 20, 5e13, 100, 5, 1e19, True),
    (50, 7.5, D, 60, 3, 1e18, True),
    (50, 15, 5e14, 100, 2, 1e19, False),
    (100, 5, 5e13, 40, 1, 1e18, True),
    (150, 10, D, 100, 1, 1e19, False),
    (150, 20, 5e13, 100, 5, 1e18, True),
]
FIELDS = [
    "hw",
    "epi_um",
    "nd",
    "V",
    "np_um",
    "nd_np",
    "np_gen",
    "solver_a",
    "center_d",
    "center_a",
    "edge_d",
    "edge_a",
    "ratio_d",
    "ratio_a",
    "valid_d",
    "valid_a",
    "fail_d",
    "fail_a",
    "t_d",
    "t_a",
]


def sweep(a, b):
    from etna.core.optimization import _sweep_one

    out = Path("lab/solver2d/notes/nplus2d_sweep_parity.csv")
    new = not out.exists()
    with out.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(FIELDS)
        for hw, epi, nd, V, npu, ndn, gen in POINTS[a:b]:
            dkw = dict(
                nplus_thickness_cm=npu * 1e-4, N_D_nplus=ndn, nplus_generation=gen
            )
            task = (0, 1, float(hw), epi * 1e-4, nd, float(V), 10)
            t = time.perf_counter()
            ra = _sweep_one(task + ("auto", dkw))
            ta = time.perf_counter() - t
            t = time.perf_counter()
            rd = _sweep_one(task + ("devsim", dkw))
            td = time.perf_counter() - t
            row = [
                hw,
                epi,
                nd,
                V,
                npu,
                ndn,
                gen,
                ra["solver"],
                rd["center_cce"],
                ra["center_cce"],
                rd["edge_cce"],
                ra["edge_cce"],
                rd["edge_center_ratio"],
                ra["edge_center_ratio"],
                rd["is_valid"],
                ra["is_valid"],
                rd["failure_reason"],
                ra["failure_reason"],
                round(td, 2),
                round(ta, 3),
            ]
            w.writerow(row)
            f.flush()
            print(
                "ROW", row, f"dC={ra['center_cce'] - rd['center_cce']:.2e}", flush=True
            )


def field():
    import etna
    from etna.api import simulation as sim

    for hw, npu, ndn, bias in [
        (25, 1.0, 1e19, -20.0),
        (25, 1.0, 1e19, -100.0),
        (50, 5.0, 1e18, -60.0),
        (150, 2.0, 1e19, -100.0),
    ]:
        cfg = etna.DeviceConfig(
            half_width_um=float(hw), nplus_thickness_um=npu, N_D_nplus=ndn
        )
        t = time.perf_counter()
        re = sim.run_field(cfg, bias_V=bias)
        te = time.perf_counter() - t
        with mock.patch.object(
            sim, "_field_engine_unsupported_reason", lambda c, b: "forced"
        ):
            t = time.perf_counter()
            rd = sim.run_field(cfg, bias_V=bias)
            td = time.perf_counter() - t
        me, md = re.mesh, rd.mesh
        same_xy = np.array_equal(me.x_coords, md.x_coords) and np.array_equal(
            me.y_coords, md.y_coords
        )
        same_nd = np.array_equal(
            me.node_values["NetDoping"], md.node_values["NetDoping"]
        )
        dpsi = np.max(np.abs(me.node_values["Potential"] - md.node_values["Potential"]))
        Ee, Ed = me.node_values["ElectricField"], md.node_values["ElectricField"]
        dE = np.max(np.abs(Ee - Ed) / np.maximum(np.abs(Ed), 1e-3 * np.abs(Ed).max()))
        print(
            f"FIELD hw={hw} np={npu}um {ndn:.0e} bias={bias}: solver {re.metadata['solver']}/"
            f"{rd.metadata['solver']} t {te:.3f}/{td:.2f} s xy_same {same_xy} nd_same {same_nd} "
            f"dpsi {dpsi:.2e} V relE {dE:.2e} contact {me.contacts[1]['position']:.6g}/"
            f"{md.contacts[1]['position']:.6g}",
            flush=True,
        )


if __name__ == "__main__":
    if sys.argv[1] == "sweep":
        sweep(int(sys.argv[2]), int(sys.argv[3]))
    else:
        field()
