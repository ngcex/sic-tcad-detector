"""Mesh convergence at one LET (15.8 keV/um, 50 V): devsim's mesher with every mesh-line
spacing divided by the given factors (etna.core.track.planar._mesh_line_scale), solved
by the adaptive scheme on the 10 um cut (or "full").

    uv run python lab/single_particle_fast/bench/mesh_convergence.py HW WIDTH|full F1 F2 ...
"""

import sys

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.single_particle import analyze_current_pulse  # noqa: E402
from etna.core.track.planar import (
    build_planar,
    ion_track_generation,
    transient_adaptive,
)  # noqa: E402

hw = float(sys.argv[1])
W = None if sys.argv[2] == "full" else float(sys.argv[2])
for f in [float(a) for a in sys.argv[3:]]:
    st = build_planar(hw, 50.0, mesh_refine=f, lateral_width_um=W)
    gen, qg = ion_track_generation(st, 15.8)
    r = transient_adaptive(st, gen)
    pi = analyze_current_pulse(
        r["times"], r["currents"], 0.0, currents_total=r["currents_total"]
    )
    print(f"refine {f} nx {st.mesh.nx} ny {st.mesh.ny}: CCE {r['Q_collected']/qg:.10f} "
          f"bal {r['cce_balance']:.10f} t95 {pi['t_collection']*1e9:.7f} {r['wall_s']:.1f}s "
          f"dark {st.dark_check}", flush=True)  # fmt: skip
