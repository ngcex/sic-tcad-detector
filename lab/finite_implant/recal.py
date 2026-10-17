"""Self-consistent recalibration WITH the fringe term, then the -30 V prediction.

Measured C of the 2x2 mm die = C_plane(V) A + C'_fringe(V) P. The paper's W labels are
eps A / C_measured. Fixed point:
  1. given C'(0), C'(-10) from a 2D run at the current (N_D_junction, N_D_bulk), solve the
     1D plane (fi2d reference builder, same y-lines; L_transition held) for (NJ, NB) with
     C_plane(V) = eps/W_label(V) - C'(V) P/A at V = 0, -10 (W_label 1.7, 9.5 um);
  2. rerun the 2D finite-implant device at the new doping (0, -10, -30 V);
  3. repeat until C' stops moving; report W_meas(-30) = eps A/(C_plane A + C' P) vs 9.73.

usage: uv run python lab/finite_implant/recal.py TAG [build key=value ...]  (e.g. Qf=-5e11)
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fi2d  # noqa: E402
from etna.core.device import DEFAULT_STACK as S  # noqa: E402
from etna.core.devsim_reset import reset_devsim_fully  # noqa: E402

EPS = 9.7 * 8.854e-14
P, A = 0.8, 0.04
WL = {0.0: 1.7e-4, 10.0: 9.5e-4}
_k = [0]


def plane_C(NJ, NB, volts=(0.0, 10.0, 30.0), dv=0.05, **pkw):
    reset_devsim_fully()
    _k[0] += 1
    info = fi2d.build(name=f"r{_k[0]}", x_edge_um=1, x_out_um=0, NJ=NJ, NB=NB, **pkw)
    fi2d.solve()
    out, vnow = {}, 0.0
    for V in volts:
        fi2d.ramp(info, vnow, V)
        vnow = V
        q = {}
        for s in (-1, 1):
            fi2d.ramp(info, V, V + s * dv, step=dv)
            q[s] = fi2d.charges(info)[0]
            fi2d.ramp(info, V + s * dv, V, step=dv)
        out[V] = abs(q[1] - q[-1]) / (2 * dv) / info["xmax"]
    return out


def solve_plane(Cf, x0, **pkw):
    tgt = {V: EPS / WL[V] - Cf[V] * P / A for V in WL}

    def res(z):
        c = plane_C(np.exp(z[0]), np.exp(z[1]), volts=(0.0, 10.0), **pkw)
        return [c[V] / tgt[V] - 1 for V in (0.0, 10.0)]

    z = fsolve(res, np.log(x0), xtol=1e-7, epsfcn=1e-6)
    return np.exp(z)


def run2d(tag, NJ, NB, extra):
    out = HERE / "notes" / f"{tag}.json"
    args = [
        sys.executable,
        str(HERE / "run_one.py"),
        str(out),
        f"NJ={NJ}",
        f"NB={NB}",
        "volts=0,10,30",
        *extra,
    ]
    with open(HERE / "notes" / f"{tag}.log", "w") as log:
        subprocess.check_call(args, stdout=log, stderr=subprocess.STDOUT)
    r = json.load(open(out))
    return {abs(v): abs(c) for v, c in zip(r["V"], r["C_e"]["0.05"])}, r["xe_cm"]


def main():
    tag, extra = sys.argv[1], sys.argv[2:]
    # plane reference must share the y-lines of the 2D run
    pkw = {k: float(v) for k, v in (a.split("=") for a in extra) if k in ("surf_ps_um", "epi_um")}
    NJ, NB = S["N_D_junction"], S["N_D_bulk"]
    hist = []
    Cf_prev = None
    for it in range(5):
        C2, xe = run2d(f"{tag}_it{it}", NJ, NB, extra)
        cp = plane_C(NJ, NB, **pkw)
        Cf = {V: C2[V] - cp[V] * xe for V in C2}
        Wm = {V: EPS * A / (cp[V] * A + Cf[V] * P) * 1e4 for V in C2}
        row = dict(
            it=it,
            NJ=NJ,
            NB=NB,
            Cf=Cf,
            Wplane={V: EPS / cp[V] * 1e4 for V in cp},
            Wmeas=Wm,
            f={V: Cf[V] * P / (cp[V] * A) for V in Cf},
        )
        hist.append(row)
        print(json.dumps(row), flush=True)
        if (
            Cf_prev
            and max(abs(Cf[V] / Cf_prev[V] - 1) for V in Cf) < 2e-3
            and abs(Wm[0.0] - 1.7) < 2e-3
            and abs(Wm[10.0] - 9.5) < 2e-3
        ):
            break
        Cf_prev = Cf
        NJ, NB = solve_plane(Cf, (NJ, NB), **pkw)
    json.dump(hist, open(HERE / "notes" / f"{tag}_recal.json", "w"), indent=1)


if __name__ == "__main__":
    main()
