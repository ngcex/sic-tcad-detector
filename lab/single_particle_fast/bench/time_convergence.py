"""Time convergence at one LET (15.8 keV/um, 50 V): the adaptive TR-BDF2 rtol sweep and
the devsim-grid twin with the step fraction and the injection step halved repeatedly
(first-order BDF1: the differences halve towards the converged value).

    uv run python lab/single_particle_fast/bench/time_convergence.py [hw]  > notes/time.out
"""

import sys

import numpy as np

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.single_particle import analyze_current_pulse  # noqa: E402
from etna.core.track.planar import (  # noqa: E402
    build_planar,
    ion_track_generation,
    transient_adaptive,
    transient_devsim_grid,
)

hw = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
LET = 15.8
st = build_planar(hw, 50.0)
gen, qg = ion_track_generation(st, LET)


def show(tag, r):
    pi = analyze_current_pulse(
        r["times"], r["currents"], 0.0, currents_total=r["currents_total"]
    )
    t = r["times"]
    hist = np.histogram(np.log10(t), bins=np.arange(-16, -6))[0]
    bal = r.get("cce_balance", float("nan"))
    print(f"{tag:28s} CCE {r['Q_collected']/qg:.10f} (balance {bal:.10f}) t95 {pi['t_collection']*1e9:.6f} ns "
          f"Ipk {pi['I_peak']:.5e} steps {r['n_steps']} lu {r['lu_n']} {r['wall_s']:.2f}s "
          f"t_end {t[-1]:.2e} conv {r['converged']} hist {hist.tolist()}", flush=True)  # fmt: skip


for rtol in (1e-3, 1e-4, 1e-5, 1e-6):
    show(f"adaptive rtol {rtol:g}", transient_adaptive(st, gen, rtol=rtol))
show("adaptive stop 1e-11", transient_adaptive(st, gen, rtol=1e-5, stop_tol=1e-11))
for f, dti in ((0.1, 1e-12), (0.05, 0.5e-12), (0.025, 0.25e-12), (0.0125, 0.125e-12)):
    show(f"grid frac {f} inj {dti:g}", transient_devsim_grid(
        st, gen, step_fraction=f, dt_inject=dti, dt_min=dti / 10, rtol_converged=1e-9))  # fmt: skip
