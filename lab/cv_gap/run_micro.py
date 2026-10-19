"""Stability check: the same fits on the Microdosimetry Fig. 3b digitisation."""

import json

from cvmodel import model_W
from etna.core.device import DEFAULT_STACK as S
from fitlib import fit_family, load_micro, report

D = load_micro("lab/cv_gap/notes/micro_fig3b_band.csv", "lab/cv_gap/notes/micro_fig3b_markers.csv")
print(len(D["V"]), "points")
out = {}
Wm, ok = model_W(D["V"], S["N_D_junction"], S["N_D_bulk"], S["L_transition"])
out["production"] = report(
    D,
    "production (no fit)",
    dict(NJ=S["N_D_junction"], NB=S["N_D_bulk"], L_um=0.987),
    Wm,
    ok,
)
for fam in ["plain", "buffer", "fringe", "area", "buffer+Nnp"]:
    r, _ = fit_family(D, fam)
    out[fam] = r
with open("lab/cv_gap/notes/micro_fits.json", "w") as f:
    json.dump(out, f, indent=1)
