"""Robustness fits: fringe at its base level, combinations, more starts. Writes extra.json."""

import json
import sys

import numpy as np

from fitlib import fit_family, load_photons

D = load_photons("lab/cv_gap/notes/photons_w.csv", "lab/cv_gap/notes/photons_w_markers.csv")
print(len(D["V"]), "points")
out = {}
for fam in [
    "buffer+fringe",
    "buffer+Nnp",
    "buffer+Nnp+fringe0",
    "buffer+Nnp+fringe",
    "buffer+Nnp+area",
]:
    r, Wm = fit_family(D, fam)
    out[fam] = r
# best model, per-bias residual table
for fam in ("buffer+Nnp", "buffer+Nnp+fringe0"):
    print(f"\n{fam}: W_model / residual (data - model) at key biases")
    for k in out[fam]["W_key"]:
        print(
            f"  {float(k):6.1f} V  W {out[fam]['W_key'][k]:.3f}  resid {out[fam]['resid_key'][k]:+.3f}"
        )
with open("lab/cv_gap/notes/extra.json", "w") as f:
    json.dump(out, f, indent=1)
