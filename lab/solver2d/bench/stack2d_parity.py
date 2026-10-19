"""2D engine vs devsim on the recalibrated stack: p+ anode 0.3-1 um and graded doping
(N_D_junction, L_transition) away from device2d's old defaults, with and without n+.

Modes:
  parity START STOP   engine internals (gate bypassed) vs devsim: dark psi/n/p, W, CCE
                      -> appends to lab/solver2d/notes/stack2d_parity.csv
  wall HW EPI_UM V0 V1 SUB_UM NDJ NDB L_UM
                      no n+: devsim 2D ramp (0.5 V steps) last converged bias and the
                      engine's dark cathode-cell field there (field-rule transfer check)
  sweep               _sweep_one(auto) vs (devsim) with the stack keys in the 9th element
  field               run_field engine vs forced devsim on stack configs

usage: uv run python lab/solver2d/bench/stack2d_parity.py MODE ...
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import nplus2d_parity as base  # noqa: E402  (engine/ref/W helpers of the n+ parity run)

CJ, CB, CL = 2.974271e15, 6.335936e13, 0.987  # NPLUS_CALIBRATION doping, L in um
LJ, LB = 2.93e15, 8.82e13  # legacy DeviceConfig doping
# (hw, epi um, N_D_bulk, V, n+ um or None, N_D_nplus, nplus_generation, p+ um, N_D_junction, L um)
POINTS = [
    (25, 10, CB, 20, 1, 1e19, False, 0.3, CJ, CL),
    (25, 10, CB, 60, 1, 1e19, False, 0.3, CJ, CL),
    (25, 10, CB, 100, 1, 1e19, False, 0.3, CJ, CL),
    (50, 10, CB, 40, 1, 1e19, False, 0.3, CJ, CL),
    (25, 10, CB, 30, 1, 1e19, True, 0.3, CJ, CL),
    (150, 10, CB, 80, 1, 1e19, False, 0.3, CJ, CL),
    (50, 10, CB, 10, None, 1e19, False, 0.3, CJ, CL),
    (25, 10, CB, 40, None, 1e19, False, 0.3, CJ, CL),
    (25, 10, LB, 20, None, 1e19, False, 1.0, LJ, CL),
    (25, 10, LB, 40, None, 1e19, False, 1.0, LJ, CL),
    (25, 5, 5e14, 40, 2, 1e19, False, 0.5, 1e16, 0.3),
    (25, 20, 5e13, 10, 5, 1e18, False, 0.5, 1e15, 3.0),
    (100, 15, 1e14, 60, None, 1e19, False, 0.3, 5e15, 2.0),
    (25, 7.5, 2e14, 30, None, 1e19, False, 0.7, 8e15, 1.5),
    (50, 20, 5e13, 70, 3, 3e18, True, 1.0, 1e15, 0.3),
    (25, 10, LB, 50, 1, 1e19, False, 1.0, LJ, CL),
]
FIELDS = [
    "hw",
    "epi_um",
    "nd",
    "V",
    "np_um",
    "nd_np",
    "np_gen",
    "sub_um",
    "ndj",
    "L_um",
    "W_d_um",
    "dW_um",
    "dpsi_V",
    "dlogn",
    "dlogp",
    "center_d",
    "max_dCCE",
    "cathode_field",
    "t_d",
    "gummel_max",
]


def kwargs(epi, nd, npu, ndn, gen, sub, ndj, L):
    kw = dict(
        epi_thickness_cm=epi * 1e-4,
        N_D_bulk=nd,
        substrate_thickness_cm=sub * 1e-4,
        N_D_junction=ndj,
        L_transition=L * 1e-4,
    )
    if npu is not None:
        kw.update(nplus_thickness_cm=npu * 1e-4, N_D_nplus=ndn, nplus_generation=gen)
    return kw


def parity(a, b):
    from etna.core.engine2d.device import Device2D
    from etna.core.engine2d.grid import grid_from_devsim

    out = Path("lab/solver2d/notes/stack2d_parity.csv")
    new = not out.exists()
    with out.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(FIELDS)
        for hw, epi, nd, V, npu, ndn, gen, sub, ndj, L in POINTS[a:b]:
            kw = kwargs(epi, nd, npu, ndn, gen, sub, ndj, L)
            grid, info, de, cce_e, te, git = base.engine(kw, hw, float(V))
            g2, i2 = grid_from_devsim(half_width_um=hw, **kw)
            d2 = Device2D(g2, i2)
            cfield = d2.dark(float(V))
            dd, cce_d, td = base.ref(kw, hw, float(V))
            inv = np.empty_like(grid.perm)
            inv[grid.perm] = np.arange(grid.N)
            dk = {k: v[inv] for k, v in dd.items()}
            assert np.allclose(dk["y"], np.repeat(grid.ys, grid.nx))
            col = slice(0, None, grid.nx)
            Wd = base.W_col(grid.ys, dk["n"][col], dk["donors"][col], info)
            We = base.W_col(grid.ys, de["n"][col], grid.donors, info)
            dCCE = cce_e - np.clip(cce_d, 0, 1)
            row = [
                hw,
                epi,
                nd,
                V,
                npu,
                ndn,
                gen,
                sub,
                ndj,
                L,
                Wd * 1e4,
                (We - Wd) * 1e4,
                float(np.max(np.abs(de["u"] - dk["u"]))),
                base.logdiff(de["n"], dk["n"]),
                base.logdiff(de["p"], dk["p"]),
                cce_d[0],
                float(np.max(np.abs(dCCE))),
                cfield,
                round(td, 2),
                git,
            ]
            w.writerow(row)
            f.flush()
            print("ROW", row, flush=True)


def wall(hw, epi, v0, v1, sub, ndj, ndb, L):
    from etna.core.charge_collection_2d import create_2d_dd_device
    from etna.core.devsim_reset import reset_devsim_fully
    from etna.core.drift_diffusion import ramp_bias
    from etna.core.engine2d.device import Device2D
    from etna.core.engine2d.grid import grid_from_devsim

    kw = kwargs(epi, ndb, None, None, False, sub, ndj, L)
    t = time.perf_counter()
    V = v0
    failed = False
    try:
        info = create_2d_dd_device(half_width_um=hw, V_bias=v0, **kw)
        while V < v1 - 1e-9:
            ramp_bias(info, V_target=V + 0.5, contact="cathode", V_step=0.5)
            V += 0.5
    except Exception:
        failed = True
    reset_devsim_fully(preserve_solver=True)
    g, gi = grid_from_devsim(half_width_um=hw, **kw)
    fld = lambda v: Device2D(g, gi).dark(v)  # noqa: E731
    msg = (
        f"WALL hw={hw:g} epi={epi:g} sub={sub:g} ndj={ndj:.4g} ndb={ndb:.4g} L={L:g}: "
        f"last converged {V:.1f} V field {fld(V):.4g}"
    )
    msg += (
        f"; failed at {V + 0.5:.1f} V field {fld(V + 0.5):.4g}"
        if failed
        else " (reached V1)"
    )
    print(msg + f" ({time.perf_counter() - t:.0f}s)", flush=True)


SWEEP_POINTS = [  # (hw, epi um, N_D_bulk, V, stack keys)
    (
        25,
        10,
        CB,
        20,
        dict(
            nplus_thickness_cm=1e-4,
            N_D_nplus=1e19,
            nplus_generation=False,
            substrate_thickness_cm=0.3e-4,
            N_D_junction=CJ,
            L_transition=CL * 1e-4,
        ),
    ),
    (
        50,
        10,
        CB,
        60,
        dict(
            nplus_thickness_cm=1e-4,
            N_D_nplus=1e19,
            nplus_generation=False,
            substrate_thickness_cm=0.3e-4,
            N_D_junction=CJ,
            L_transition=CL * 1e-4,
        ),
    ),
    (
        25,
        10,
        CB,
        40,
        dict(substrate_thickness_cm=0.3e-4, N_D_junction=CJ, L_transition=CL * 1e-4),
    ),
    (
        25,
        10,
        LB,
        20,
        dict(substrate_thickness_cm=1e-4, N_D_junction=LJ, L_transition=CL * 1e-4),
    ),
    (
        25,
        5,
        5e14,
        60,
        dict(
            nplus_thickness_cm=2e-4,
            N_D_nplus=1e18,
            nplus_generation=True,
            substrate_thickness_cm=0.5e-4,
            N_D_junction=1e16,
            L_transition=0.3e-4,
        ),
    ),
    (
        100,
        20,
        5e13,
        100,
        dict(
            nplus_thickness_cm=1e-4,
            N_D_nplus=1e19,
            nplus_generation=False,
            substrate_thickness_cm=0.3e-4,
            N_D_junction=1e15,
            L_transition=3e-4,
        ),
    ),
]


def sweep():
    from etna.core.optimization import _sweep_one

    for hw, epi, nd, V, dkw in SWEEP_POINTS:
        task = (0, 1, float(hw), epi * 1e-4, nd, float(V), 10)
        t = time.perf_counter()
        ra = _sweep_one(task + ("auto", dkw))
        ta = time.perf_counter() - t
        t = time.perf_counter()
        rd = _sweep_one(task + ("devsim", dkw))
        td = time.perf_counter() - t
        print(
            f"SWEEP hw={hw} epi={epi} nd={nd:.3g} V={V} {dkw}: solver {ra['solver']} "
            f"dC {ra['center_cce'] - rd['center_cce']:.2e} dE {ra['edge_cce'] - rd['edge_cce']:.2e} "
            f"dratio {ra['edge_center_ratio'] - rd['edge_center_ratio']:.2e} "
            f"valid {ra['is_valid']}/{rd['is_valid']} fd {ra['is_fully_depleted']}/"
            f"{rd['is_fully_depleted']} Vfd {ra['V_fd']:.4g}/{rd['V_fd']:.4g} "
            f"fail {ra['failure_reason']!r}/{rd['failure_reason']!r} t {ta:.2f}/{td:.1f} s",
            flush=True,
        )


def field():
    import etna
    from etna.api import simulation as sim

    cases = [
        (
            25,
            -20.0,
            dict(
                nplus_thickness_um=1.0,
                substrate_thickness_um=0.3,
                N_D_junction=CJ,
                N_D_bulk=CB,
                L_transition_um=CL,
            ),
        ),
        (
            50,
            -50.0,
            dict(
                nplus_thickness_um=1.0,
                substrate_thickness_um=0.3,
                N_D_junction=CJ,
                N_D_bulk=CB,
                L_transition_um=CL,
            ),
        ),
        (
            25,
            -100.0,
            dict(
                nplus_thickness_um=1.0,
                substrate_thickness_um=0.3,
                N_D_junction=CJ,
                N_D_bulk=CB,
                L_transition_um=CL,
            ),
        ),
        (
            25,
            -30.0,
            dict(
                substrate_thickness_um=0.3,
                N_D_junction=CJ,
                N_D_bulk=CB,
                L_transition_um=CL,
            ),
        ),
        (
            50,
            -20.0,
            dict(
                substrate_thickness_um=1.0,
                N_D_junction=LJ,
                N_D_bulk=LB,
                L_transition_um=CL,
            ),
        ),
    ]
    for hw, bias, kw in cases:
        cfg = etna.DeviceConfig(half_width_um=float(hw), **kw)
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
            f"FIELD hw={hw} bias={bias} {kw}: solver {re.metadata['solver']}/"
            f"{rd.metadata['solver']} t {te:.3f}/{td:.2f} s xy_same {same_xy} "
            f"nd_same {same_nd} dpsi {dpsi:.2e} V relE {dE:.2e}",
            flush=True,
        )


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "parity":
        parity(int(sys.argv[2]), int(sys.argv[3]))
    elif mode == "wall":
        a = sys.argv[2:]
        wall(
            float(a[0]),
            float(a[1]),
            float(a[2]),
            float(a[3]),
            float(a[4]),
            float(a[5]),
            float(a[6]),
            float(a[7]),
        )
    elif mode == "sweep":
        sweep()
    else:
        field()
