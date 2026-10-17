"""Model-independent piecewise log-linear donor profile fitted to the digitised curve."""

import json

import numpy as np

from fitlib import KNOTS_UM, fit_knots, load_photons, report

D = load_photons("lab/cv_gap/notes/photons_w.csv", "lab/cv_gap/notes/photons_w_markers.csv")
out = {}
for lam, Qf in [(1.0, None), (0.3, None), (3.0, None), (1.0, 0.0), (1.0, 1e12)]:
    p, Wm, ok = fit_knots(D, lam=lam, Qf=Qf)
    tag = f"knots lam {lam} Qf {Qf}"
    r = report(D, tag, {}, Wm, ok, npar=len(p))
    r["logN_knots"] = [float(v) for v in p]
    out[tag] = r
    print("   depth(um): " + " ".join(f"{k:5.2f}" for k in KNOTS_UM))
    print("   log10 N  : " + " ".join(f"{v:5.2f}" for v in p))
with open("lab/cv_gap/notes/knots.json", "w") as f:
    json.dump(out, f, indent=1)
