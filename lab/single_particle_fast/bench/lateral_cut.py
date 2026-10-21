"""Lateral cut: the adaptive CCE / t_collection at LET 15.8, 50 V with the device cut at
the given widths (um beyond the track) against the full half-width ("full").

    uv run python lab/single_particle_fast/bench/lateral_cut.py HW W1 W2 ... full
"""

import sys

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.single_particle import analyze_current_pulse  # noqa: E402
from etna.core.track.planar import (  # noqa: E402
    build_planar,
    generated_charge,
    ion_track_generation,
    transient_adaptive,
)

hw = float(sys.argv[1])
widths = [None if a == "full" else float(a) for a in sys.argv[2:]]
for W in widths:
    st = build_planar(hw, 50.0, lateral_width_um=W)
    gen, qg = ion_track_generation(st, 15.8)
    qf = generated_charge(st, 15.8)
    r = transient_adaptive(st, gen)
    pi = analyze_current_pulse(
        r["times"], r["currents"], 0.0, currents_total=r["currents_total"]
    )
    print(f"W {W} nx {st.mesh.nx} N {st.N}: CCE {r['Q_collected']/qg:.10f} bal {r['cce_balance']:.10f} "
          f"abs {r['Q_abs_sum']/qg:.10f} t95 {pi['t_collection']*1e9:.7f} Ipk_tot {pi['I_peak']:.6e} "
          f"qg/qfull-1 {qg/qf-1:.1e} {r['wall_s']:.1f}s steps {r['n_steps']} build {st.build_s:.1f}s",
          flush=True)  # fmt: skip
