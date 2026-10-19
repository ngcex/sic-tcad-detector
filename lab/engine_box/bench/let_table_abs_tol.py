"""Is the LET dependence of the CCE(LET) table physics or Newton tolerance?

The 2D single-particle transient at the table's LETs is a low-injection problem
(peak excess density 2.6e9 cm^-3 at 0.5 keV/um, 2.6e12 at 500 keV/um, against
N_D_bulk 3.4e13): if the discrete equations are linear in the excess carriers, the
whole transient scales with the LET and CCE is LET-independent. The transient Newton
stops on TRANSIENT_ABS_ERROR = 1e10 cm^-3 (absolute carrier update), which is ABOVE the
excess density at low LET: the absolute tolerance does not scale with the LET and
breaks that invariance. Running each LET with the absolute tolerance scaled by
LET / 500 (same relative accuracy as the 500 keV/um run) tests it: if the low-LET
CCE then matches the 500 keV/um CCE, the LET dependence of the table is numerical.

    uv run python lab/engine_box/bench/let_table_abs_tol.py HW LET:ABS[:REL] [...]

Result (hw 50, 2026-10-18): the hypothesis is FALSE. The absolute tolerance is not
binding (LET 0.5 with 1e7, LET 15.8 with 3.2e8 and LET 500 with 1e8 are bit-identical
to 1e10) and TRANSIENT_REL_ERROR 1e-8 -> 1e-11 moves LET 0.5 by 3.4e-9 (0.9999623609 ->
0.9999623643). The LET dependence comes from the dark-current subtraction
(|I| - |I_dark| with I_dark of the opposite sign at the anode; build_cce_let_table Notes).

Builds ONE device and restores its dark state per run (build_cce_let_table pattern).
"""

import sys
import time

import etna
from etna.core import single_particle as sp
from etna.core.charge_collection_2d import create_2d_dd_device
from etna.core.devsim_reset import delete_device_and_mesh

etna.set_devsim_quiet(True)

hw = float(sys.argv[1])
runs = [tuple(float(v) for v in a.split(":")) for a in sys.argv[2:]]
info = create_2d_dd_device(half_width_um=hw, V_bias=50.0)
dark = sp._snapshot_dark_state(info)
default_abs, default_rel = sp.TRANSIENT_ABS_ERROR, sp.TRANSIENT_REL_ERROR
try:
    for run in runs:
        L, abs_err = run[:2]
        rel_err = run[2] if len(run) > 2 else default_rel
        sp._restore_dark_state(info, dark)
        sp.TRANSIENT_ABS_ERROR = abs_err
        sp.TRANSIENT_REL_ERROR = rel_err
        t0 = time.perf_counter()
        gen, qg = sp.ion_track_generation_2d(info, L)
        r = sp.simulate_single_particle(info, gen)
        print(
            f"LET {L:g} abs_error {abs_err:.1e} rel_error {rel_err:.1e}: CCE_raw={r['Q_collected'] / qg:.10f} "
            f"steps={len(r['times'])} converged={r['converged']} "
            f"({time.perf_counter() - t0:.1f} s)",
            flush=True,
        )
finally:
    sp.TRANSIENT_ABS_ERROR, sp.TRANSIENT_REL_ERROR = default_abs, default_rel
    delete_device_and_mesh(info["device_name"])
