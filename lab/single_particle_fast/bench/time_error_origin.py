"""Origin of the devsim-grid time error in the CCE (hw 2, 50 V, LET 15.8): for a LINEAR
transient the backward-Euler (and TR-BDF2) collected charge does not depend on the step
sequence (A sum_k h_k x_k = x_end - x_0), so a grid-dependent CCE needs a nonlinearity.
Candidate: SRH of the excess pairs in the depleted epi, U ~ dn dp / (tau_p dn + tau_n dp)
(homogeneous of degree 1 in the excess, so LET-independent, but not linear). Check: the
same twin / adaptive pair with the SRH lifetimes x 1e6 (SRH negligible).

    uv run python lab/single_particle_fast/bench/time_error_origin.py
"""

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.track.planar import (  # noqa: E402
    build_planar,
    ion_track_generation,
    transient_adaptive,
    transient_devsim_grid,
)

st = build_planar(2.0, 50.0)
gen, qg = ion_track_generation(st, 15.811388300841896)
for scale in (1.0, 1e6):
    st.system.taun = st.system.taun * scale
    st.system.taup = st.system.taup * scale
    g = transient_devsim_grid(st, gen, rtol_converged=1e-10)
    a = transient_adaptive(st, gen)
    cg, ca = g["Q_signed"] / qg, a["Q_collected"] / qg
    print(f"tau x {scale:g}: CCE grid {cg:.10f} adaptive {ca:.10f} adaptive-grid "
          f"{ca - cg:+.2e} (1 - CCE adaptive {1 - ca:.2e})", flush=True)  # fmt: skip
    st.system.taun = st.system.taun / scale
    st.system.taup = st.system.taup / scale
