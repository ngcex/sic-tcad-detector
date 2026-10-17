"""Fringe capacitance per unit edge length from run_one.py outputs.

C'_fringe(V) = C_2D(V) - C_plane(V) * x_edge      (F/cm, C_plane per area from ref.json,
                                                    same y-lines)
dW_eq(V)     = C'_fringe / C_plane                  (um: equivalent extra active width)
f(V)         = C'_fringe * P / (C_plane * A)        (P = 0.8 cm, A = 0.04 cm^2: 2x2 mm die)
W_meas(V)    = eps A / (C_plane A + C'_fringe P)    (what eps*A/C would read)

usage: uv run python lab/finite_implant/analyse.py notes/ref.json notes/X.json ...
"""

import json
import sys

EPS = 9.7 * 8.854e-14
P, A = 0.8, 0.04
NEED = 0.0256  # fraction of the plane C missing at -30 V (9.979 -> 9.73 um)


def load(p):
    return json.load(open(p))


def main():
    ref = load(sys.argv[1])
    dv0 = "0.05"
    cpl = {v: abs(c) / ref["xmax_cm"] for v, c in zip(ref["V"], ref["C_e"][dv0])}
    rows = []
    for path in sys.argv[2:]:
        r = load(path)
        tag = path.rsplit("/", 1)[-1].replace(".json", "")
        cfg = r["cfg"]
        xe = r["xmax_cm"] if cfg.get("null") else r["xe_cm"]
        print(f"\n== {tag}  {cfg}  nodes {r['nodes']}  t {r['t']:.0f} s")
        print(
            "   V    W_plane(um)  C'_fringe(F/cm)  dW_eq(um)  f(2x2mm,%)  W_meas(um)"
            "  [dv 0.1: dW_eq]  [holes: dW_eq]  lateral depl. beyond edge (um) at y-y_j=0/1/3/5/8"
        )
        for i, V in enumerate(r["V"]):
            c = cpl[V]
            Cf = abs(r["C_e"][dv0][i]) - c * xe
            dw = Cf / c * 1e4
            f = Cf * P / (c * A)
            wm = EPS * A / (c * A + Cf * P) * 1e4
            extra = ""
            if "0.1" in r["C_e"]:
                extra += f"  {(abs(r['C_e']['0.1'][i]) - c * xe) / c * 1e4:8.3f}"
            extra += f"  {(abs(r['C_h'][dv0][i]) - c * xe) / c * 1e4:8.3f}"
            lat = r["lateral"].get(str(V), {})
            ls = "/".join(f"{lat[k]:.1f}" for k in lat) if lat else ""
            print(
                f"{V:5.0f}  {EPS / c * 1e4:9.4f}  {Cf:14.4e}  {dw:9.3f}  {100 * f:9.3f}"
                f"  {wm:9.4f}{extra}   {ls}"
            )
            rows.append(dict(tag=tag, V=V, Cf=Cf, dW=dw, f=f, Wm=wm))
    print(
        f"\nTarget at -30 V: f = {100 * NEED:.2f} % (dW_eq = {NEED * A / P * 1e4:.1f} um) to bring"
        " eps*A/C from 9.979 to 9.73 um with the plane held."
    )
    return rows


if __name__ == "__main__":
    main()
