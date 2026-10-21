"""Parity of the fast single-ion transient with the devsim references (data/*.npz from
devsim_reference.py): same device (devsim's mesh imported verbatim), generation, time
grid (scheme devsim_grid), plus the adaptive TR-BDF2 run (full width and the 10 um cut).

    uv run python lab/single_particle_fast/bench/parity.py [hw ...]  > notes/parity.out

Per case: Q_generated ratio; series agreement (max |I_fast - (I_devsim - I_dark)| / peak
over the common steps, conduction and displacement); CCE of devsim's estimator
(sum (|I| - |I_dark|) dt), of the sign-aware devsim sum |sum (I - I_dark) dt|, of the
fast twin (I_dark = 0 by construction); t_collection; the adaptive CCE / t_collection
and the wall times.
"""

import glob
import os
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

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
EXACT_LETS = (0.5, 5.0, 10.0, 15.811388300841896, 500.0)


def t95(r):
    return analyze_current_pulse(
        r["times"], r["currents"], 0.0, currents_total=r["currents_total"]
    )["t_collection"]


if __name__ == "__main__":
    only = {float(a) for a in sys.argv[1:]}
    setups = {}
    for f in sorted(glob.glob(os.path.join(DATA, "hw*_let*.npz"))):
        name = os.path.basename(f)[:-4]
        hw = float(name.split("_")[0][2:])
        # exact LET of the reference (the file name keeps 6 digits)
        LET = float(name.split("let")[1])
        LET = min(EXACT_LETS, key=lambda v: abs(v - LET) / v)
        V = float(name.split("_V")[1].split("_")[0]) if "_V" in name else 50.0
        if only and hw not in only:
            continue
        ref = np.load(f)
        if (hw, V) not in setups:
            setups[hw, V] = build_planar(hw, V)
            print(f"# hw {hw:g}, {V:g} V: fast build {setups[hw, V].build_s:.2f} s, "
                  f"{setups[hw, V].N} nodes, devsim build {float(ref['build_s']):.1f} s")  # fmt: skip
        st = setups[hw, V]
        gen, qg = ion_track_generation(st, LET)
        qd = float(ref["Q_gen"])
        g = transient_devsim_grid(st, gen)
        Id = float(ref["I_dark"])
        n = min(len(g["times"]), len(ref["times"]))
        sig = ref["currents"][:n] - Id
        peak = np.max(np.abs(sig))
        dI = np.max(np.abs(g["currents"][:n] - sig)) / peak
        dD = np.max(np.abs(g["currents_displacement"][:n] - ref["disp"][:n])) / np.max(
            np.abs(ref["disp"][:n])
        )
        dts = np.diff(ref["times"], prepend=0.0)
        cce_dev = float(ref["Q_collected"]) / qd
        cce_dev_signed = abs(float(np.sum((ref["currents"] - Id) * dts))) / qd
        cce_twin = g["Q_collected"] / qg
        a = transient_adaptive(st, gen)
        cce_a = a["Q_collected"] / qg
        print(f"hw {hw:g} {V:g} V LET {LET:.6g}: Q_gen ratio-1 {qg / qd - 1:.1e} | steps devsim "
              f"{len(ref['times'])} twin {len(g['times'])} | series dI/peak {dI:.1e} "
              f"dDisp/peak {dD:.1e}")  # fmt: skip
        print(f"  CCE devsim {cce_dev:.10f} devsim-signed {cce_dev_signed:.10f} twin "
              f"{cce_twin:.10f} | twin-devsim {cce_twin - cce_dev:+.2e} twin-signed "
              f"{cce_twin - cce_dev_signed:+.2e} | I_dark*t_end/Q_gen "
              f"{2 * abs(Id) * ref['times'][-1] / qd:.1e}")  # fmt: skip
        print(f"  t95 devsim {float(ref['t_collection']) * 1e9:.7f} twin {t95(g) * 1e9:.7f} "
              f"ns | adaptive CCE {cce_a:.10f} (balance {a['cce_balance']:.10f}, ledger "
              f"{a['ledger_residual']:.1e}) t95 {t95(a) * 1e9:.7f} ns, adaptive-twin "
              f"{cce_a - cce_twin:+.2e}")  # fmt: skip
        print(f"  wall: devsim {float(ref['run_s']):.1f} s, twin {g['wall_s']:.2f} s "
              f"({g['lu_n']} LU), adaptive {a['wall_s']:.2f} s ({a['n_steps']} steps, "
              f"{a['lu_n']} LU)", flush=True)  # fmt: skip
        if hw > 10:
            stc = build_planar(hw, 50.0, lateral_width_um=10.0)
            genc, qgc = ion_track_generation(stc, LET)
            c = transient_adaptive(stc, genc)
            print(f"  10 um cut ({stc.N} nodes, build {stc.build_s:.2f} s): adaptive CCE "
                  f"{c['Q_collected'] / qgc:.10f} (cut - full {c['Q_collected'] / qgc - cce_a:+.1e}) "
                  f"t95 {t95(c) * 1e9:.7f} ns, {c['wall_s']:.2f} s", flush=True)  # fmt: skip
