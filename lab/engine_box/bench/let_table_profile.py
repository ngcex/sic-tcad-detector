"""Profile one CCE(LET) table point (build_cce_let_table inner loop).

Splits the wall time into device build (create_2d_dd_device: mesh, Poisson,
equilibrium, DD setup, bias ramp), ion-track generation and the transient
(simulate_single_particle), and counts devsim.solve calls / Newton iterations per
solve type.

    uv run python lab/engine_box/bench/let_table_profile.py [half_width_um] [LET]
"""

import cProfile
import collections
import pstats
import sys
import time

import devsim

import etna
from etna.core import single_particle as sp
from etna.core.charge_collection_2d import create_2d_dd_device
from etna.core.devsim_reset import delete_device_and_mesh

etna.set_devsim_quiet(True)

HW = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0
LET = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0

_orig_solve = devsim.solve
stats = collections.defaultdict(lambda: {"calls": 0, "iters": 0, "time": 0.0})
phase = {"name": "build"}


def counting_solve(*args, **kw):
    kw = dict(kw)
    kw["info"] = True
    t0 = time.perf_counter()
    try:
        res = _orig_solve(*args, **kw)
    finally:
        s = stats[(phase["name"], kw.get("type"))]
        s["calls"] += 1
        s["time"] += time.perf_counter() - t0
    if isinstance(res, dict) and "iterations" in res:
        s["iters"] += len(res["iterations"])
    return res


devsim.solve = counting_solve

prof = cProfile.Profile()
t0 = time.perf_counter()
prof.enable()
info = create_2d_dd_device(half_width_um=HW, V_bias=50.0)
t_build = time.perf_counter() - t0
n_nodes = len(
    devsim.get_node_model_values(
        device=info["device_name"], region=info["region_name"], name="x"
    )
)
phase["name"] = "generation"
t1 = time.perf_counter()
gen, q_gen = sp.ion_track_generation_2d(info, LET)
t_gen = time.perf_counter() - t1
phase["name"] = "transient"
t2 = time.perf_counter()
res = sp.simulate_single_particle(info, gen)
t_tr = time.perf_counter() - t2
prof.disable()
delete_device_and_mesh(info["device_name"])

print(f"hw={HW} um LET={LET} nodes={n_nodes}")
print(f"build {t_build:.1f} s | generation {t_gen:.2f} s | transient {t_tr:.1f} s")
print(
    f"transient steps={len(res['times'])} t_end={res['times'][-1]:.3e} s "
    f"converged={res['converged']} CCE={res['Q_collected']/q_gen:.9f}"
)
for (ph, ty), s in sorted(stats.items()):
    print(
        f"  {ph:10s} {ty:15s} calls={s['calls']:4d} newton_iters={s['iters']:5d} "
        f"time={s['time']:.1f} s ({s['time']/max(s['calls'],1):.3f} s/call)"
    )
pstats.Stats(prof).sort_stats("tottime").print_stats(12)
