"""Fit the 1D candidate models to the digitised Photons Fig. 4b eps*A/C curve.

usage (repo root): PYTHONPATH=lab/cv_gap uv run python lab/cv_gap/fit_models.py lab/cv_gap/notes/photons_w.csv lab/cv_gap/notes/photons_w_markers.csv lab/cv_gap/notes/fits.json
"""

import json
import sys
import time

import numpy as np
from scipy.optimize import least_squares

from cvmodel import model_W

# ---------------- data table with uncertainties
band = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
mk = np.genfromtxt(sys.argv[2], delimiter=",", names=True)
rows = []  # (V, W, sigma_W, sigma_V, source)
# knee: isolated diamonds, snapped to the 250 mV grid (their centroids sit within 0.04 V of it)
for V, Wv in zip(mk["V"], mk["W_um"]):
    Vn = round(V * 4) / 4
    if -3.3 < V < -0.1:
        rows.append((Vn, Wv, 0.03, 0.04, "marker"))
rows.append((0.0, 1.68, 0.02, 0.0, "label"))
rows.append((-10.0, 9.52, 0.01, 0.0, "label"))
rows.append((-30.0, 9.73, 0.01, 0.0, "label"))
for V, Wv, sd, wt in zip(band["V"], band["W_um"], band["band_sd_um"], band["weight"]):
    if not np.isfinite(Wv) or wt < 10:
        continue  # frame edge / annotation dots
    if V in (-10.0, -30.0, -9.75, -10.25, -10.5, -9.5, -29.75, -30.25) or V > -3.5:
        continue
    sig = 0.02 if V <= -5.0 else 0.06
    rows.append((V, Wv, sig, 0.02, "band"))
rows.sort()
V = np.array([r[0] for r in rows])
Wd = np.array([r[1] for r in rows])
sW = np.array([r[2] for r in rows])
sV = np.array([r[3] for r in rows])
src = [r[4] for r in rows]
slope = np.gradient(Wd, V)
sig = np.sqrt(sW**2 + (slope * sV) ** 2)
print(
    f"{len(V)} data points: {sum(s=='marker' for s in src)} markers, {sum(s=='band' for s in src)} band, 3 labels"
)
KEY = [0, -1, -2, -3, -5, -7, -10, -15, -20, -30, -40, -50, -58]


def at(Varr, Warr, v):
    i = np.argmin(np.abs(Varr - v))
    return Warr[i] if abs(Varr[i] - v) < 0.13 else np.nan


def chi2(Wm):
    return float((((Wd - Wm) / sig) ** 2).sum())


def report(name, pars, Wm, ok, extra=""):
    c2 = chi2(Wm)
    res = Wd - Wm
    knee = np.array([s == "marker" for s in src]) | (V == 0)
    plat = V <= -5
    r = {
        "name": name,
        "pars": pars,
        "chi2": c2,
        "dof": int(len(V) - len(pars)),
        "rms_knee_um": float(np.sqrt((res[knee] ** 2).mean())),
        "rms_plateau_um": float(np.sqrt((res[plat] ** 2).mean())),
        "W_key": {str(v): float(at(V, Wm, v)) for v in KEY},
        "converged": ok,
    }
    print(
        f"{name:34s} chi2 {c2:9.1f} (dof {r['dof']})  rms knee {r['rms_knee_um']:.3f}  plateau {r['rms_plateau_um']:.4f} um  "
        + " ".join(f"{k}={v:.3g}" for k, v in pars.items())
        + extra
    )
    return r


def unpack(p, family):
    d = dict(NJ=10 ** p[0], NB=10 ** p[1], L=10 ** p[2] * 1e-4)
    kw = {}
    if family == "buffer":
        kw["Lb_um"] = 10 ** p[3]
    elif family == "area":
        kw["area_scale"] = p[3]
    elif family == "fringe":
        kw["Qf"] = p[3] * 1e12
    elif family == "buffer+area":
        kw["Lb_um"], kw["area_scale"] = 10 ** p[3], p[4]
    elif family == "buffer+fringe":
        kw["Lb_um"], kw["Qf"] = 10 ** p[3], p[4] * 1e12
    elif family == "buffer+Nnp":
        kw["Lb_um"], kw["N_np"] = 10 ** p[3], 10 ** p[4]
    return d, kw


def resid(p, family, epi):
    d, kw = unpack(p, family)
    Wm, ok = model_W(V, d["NJ"], d["NB"], d["L"], epi_um=epi, **kw)
    if not ok or not np.all(np.isfinite(Wm)):
        return np.full_like(Wd, 1e3)
    return (Wd - Wm) / sig


BOUNDS = {
    "plain": ([14, 12, -1.3], [17.5, 15, 0.7]),
    "buffer": ([14, 12, -1.3, -1.3], [17.5, 15, 0.7, 0.8]),
    "area": ([14, 12, -1.3, 0.7], [17.5, 15, 0.7, 1.4]),
    "fringe": ([14, 12, -1.3, 0.0], [17.5, 15, 0.7, 2.0]),
    "buffer+area": ([14, 12, -1.3, -1.3, 0.7], [17.5, 15, 0.7, 0.8, 1.4]),
    "buffer+fringe": ([14, 12, -1.3, -1.3, 0.0], [17.5, 15, 0.7, 0.8, 2.0]),
    "buffer+Nnp": ([14, 12, -1.3, -1.3, 16.0], [17.5, 15, 0.7, 0.8, 19.5]),
}
STARTS = {
    "plain": [[15.47, 13.92, -0.006], [16.0, 13.6, -0.5], [15.0, 13.5, 0.0]],
    "buffer": [
        [15.47, 13.92, -0.006, -0.28],
        [16.0, 13.6, -0.5, 0.0],
        [15.5, 13.7, -0.3, 0.3],
    ],
    "area": [[15.47, 13.92, -0.006, 1.0], [16.0, 13.6, -0.5, 1.03]],
    "fringe": [[15.47, 13.92, -0.006, 1.0], [16.0, 13.6, -0.5, 0.5]],
    "buffer+area": [[16.0, 13.6, -0.5, 0.0, 1.0]],
    "buffer+fringe": [[16.0, 13.6, -0.5, 0.0, 0.5]],
    "buffer+Nnp": [[16.0, 13.6, -0.5, 0.0, 19.0]],
}


def fit(family, epi=10.0, label=None):
    best = None
    for p0 in STARTS[family]:
        t = time.time()
        r = least_squares(
            resid,
            p0,
            bounds=BOUNDS[family],
            args=(family, epi),
            xtol=1e-10,
            ftol=1e-10,
            max_nfev=400,
        )
        c = float((r.fun**2).sum())
        if best is None or c < best[0]:
            best = (c, r.x)
    d, kw = unpack(best[1], family)
    Wm, ok = model_W(V, d["NJ"], d["NB"], d["L"], epi_um=epi, **kw)
    pars = dict(NJ=d["NJ"], NB=d["NB"], L_um=d["L"] * 1e4, epi_um=epi, **kw)
    return report(label or f"{family} fit (epi {epi} um)", pars, Wm, ok)


out = {
    "data": [
        dict(V=float(v), W=float(w), sigma=float(s), source=t)
        for v, w, s, t in zip(V, Wd, sig, src)
    ]
}
from etna.core.device import DEFAULT_STACK as S

Wm, ok = model_W(V, S["N_D_junction"], S["N_D_bulk"], S["L_transition"])
out["production"] = report(
    "production DEFAULT_STACK (no fit)",
    dict(NJ=S["N_D_junction"], NB=S["N_D_bulk"], L_um=S["L_transition"] * 1e4),
    Wm,
    ok,
)
# the four 3-point candidates as they stand (LOG numbers), no refit
for name, kw in [
    (
        "3pt buffer 529 nm (LOG)",
        dict(NJ=2.9888e15, NB=6.4218e13, L=0.987e-4, Lb_um=0.529),
    ),
    (
        "3pt epi 9.75 um (LOG)",
        dict(NJ=2.9844e15, NB=6.9058e13, L=0.987e-4, epi_um=9.75),
    ),
    (
        "3pt fringe Qf 2e12 (LOG recal)",
        dict(NJ=2.99621e15, NB=5.00892e13, L=0.987e-4, Qf=2e12),
    ),
]:
    L = kw.pop("L")
    Wm, ok = model_W(V, kw.pop("NJ"), kw.pop("NB"), L, **kw)
    out[name] = report(name, dict(L_um=L * 1e4, **kw), Wm, ok)
# area +2.5 %: production profile recalibrated is what the LOG implies; here the plain
# production profile with area_scale chosen to hit 9.73 at -30 V
Wm, ok = model_W(
    V, S["N_D_junction"], S["N_D_bulk"], S["L_transition"], area_scale=9.979 / 9.73
)
out["3pt area +2.56 %"] = report(
    "3pt area x1.0256 (prod profile)", dict(area_scale=9.979 / 9.73), Wm, ok
)

print("\n== refits on the full curve ==")
out["fit_plain"] = fit("plain")
out["fit_buffer"] = fit("buffer")
out["fit_area"] = fit("area")
out["fit_fringe"] = fit("fringe")
for epi in (9.0, 9.25, 9.5, 9.75, 10.25, 10.5, 11.0, 12.0):
    out[f"fit_plain_epi{epi}"] = fit("plain", epi=epi)
print("\n== combinations ==")
out["fit_buffer_area"] = fit("buffer+area")
out["fit_buffer_fringe"] = fit("buffer+fringe")
out["fit_buffer_Nnp"] = fit("buffer+Nnp")
for epi in (9.5, 10.5, 11.0):
    out[f"fit_buffer_epi{epi}"] = fit("buffer", epi=epi)
with open(sys.argv[3], "w") as f:
    json.dump(out, f, indent=1, default=float)
