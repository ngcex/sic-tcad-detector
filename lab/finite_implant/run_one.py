"""One configuration per process: C(V) = dQ/dV of the finite-implant device (or the reference).

usage: uv run python lab/finite_implant/run_one.py OUT.json key=value ...
  keys = fi2d.build kwargs (x_edge_um, x_out_um, t_ox_um, fp_um, Qf, edge_ps_um, far_ps_um,
  epi_um) plus dv (V, default 0.05) and volts (comma list of reverse cathode biases, V > 0)
"""

import json
import sys
import time

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import fi2d  # noqa: E402

VOLTS = (0, 2, 5, 10, 15, 20, 30, 40)


def main():
    out = sys.argv[1]
    kw = {}
    for a in sys.argv[2:]:
        k, v = a.split("=")
        kw[k] = v
    dvs = [float(x) for x in kw.pop("dv", "0.05").split(",")]
    volts = [float(x) for x in kw.pop("volts", ",".join(map(str, VOLTS))).split(",")]
    bkw = {k: float(v) for k, v in kw.items()}
    t0 = time.time()
    info = fi2d.build(**bkw)
    fi2d.solve()
    res = dict(cfg=bkw, nodes=info["nodes"], xe_cm=info["xe"], xmax_cm=info["xmax"], V=[], C_e={}, C_h={}, lateral={})
    vnow = 0.0
    for V in volts:
        fi2d.ramp(info, vnow, V)
        vnow = V
        res["V"].append(-V)
        for dv in dvs:
            qs = {}
            for s in (-1, 1):
                fi2d.ramp(info, vnow, V + s * dv, step=dv)
                qs[s] = fi2d.charges(info)
                fi2d.ramp(info, V + s * dv, V, step=dv)
            Ce = (qs[1][0] - qs[-1][0]) / (2 * dv)
            Ch = (qs[1][1] - qs[-1][1]) / (2 * dv)
            res["C_e"].setdefault(str(dv), []).append(Ce)
            res["C_h"].setdefault(str(dv), []).append(Ch)
        if info["has_ox"]:
            f = fi2d.fields(info)
            lat = {}
            for yd in (0.0, 1.0, 3.0, 5.0, 8.0):
                yv = info["yj"] + yd * 1e-4 if yd > 0 else 0.0
                ys = np.unique(f["y"])
                yy = ys[np.argmin(np.abs(ys - yv))]
                sel = np.abs(f["y"] - yy) < 1e-12
                xs, nr = f["x"][sel], (f["n"] / f["D"])[sel]
                o = np.argsort(xs)
                xs, nr = xs[o], nr[o]
                # rightmost x beyond the edge where n/N_D < 0.5 (depleted)
                dep = xs[(nr < 0.5) & (xs > info["xe"])]
                lat[f"{yd:g}"] = (
                    float((dep.max() - info["xe"]) * 1e4) if dep.size else 0.0
                )
            res["lateral"][str(-V)] = lat
        print(
            f"V={-V:+.0f}  C_e={res['C_e'][str(dvs[0])][-1]:.6e} F/cm  "
            f"lat={res['lateral'].get(str(-V))}  t={time.time()-t0:.0f}s",
            flush=True,
        )
    res["t"] = time.time() - t0
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
