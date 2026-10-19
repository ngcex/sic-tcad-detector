"""Shared data table + fit helpers for the C-V gap study."""

import numpy as np
from scipy.optimize import least_squares

from cvmodel import (
    N_NP,
    batch_from_profiles,
    fringe_Cprime,
    model_W,
    node_positions,
    sweep_W,
    SUB,
    NP,
    P_OVER_A,
)
from etna.core.batched.material import EPS

KEY = [0, -1, -2, -3, -5, -7, -10, -15, -20, -30, -40, -50, -58]


def load_photons(band_csv, markers_csv):
    band = np.genfromtxt(band_csv, delimiter=",", names=True)
    mk = np.genfromtxt(markers_csv, delimiter=",", names=True)
    rows = []
    for V, Wv in zip(mk["V"], mk["W_um"]):
        if -3.3 < V < -0.1:
            rows.append((round(V * 4) / 4, Wv, 0.03, 0.04, "marker"))
    rows += [
        (0.0, 1.68, 0.02, 0.0, "label"),
        (-10.0, 9.52, 0.01, 0.0, "label"),
        (-30.0, 9.73, 0.01, 0.0, "label"),
    ]
    for V, Wv, wt in zip(band["V"], band["W_um"], band["weight"]):
        if not np.isfinite(Wv) or wt < 10:
            continue
        if V in (-10.0, -30.0, -9.75, -10.25, -10.5, -9.5, -29.75, -30.25) or V > -3.5:
            continue
        rows.append((V, Wv, 0.02 if V <= -5.0 else 0.06, 0.02, "band"))
    return _table(rows)


def load_micro(band_csv, markers_csv):
    """Same data set re-plotted in the Microdosimetry paper (52.6 px/um, 9.75 px/V)."""
    band = np.genfromtxt(band_csv, delimiter=",", names=True)
    mk = np.genfromtxt(markers_csv, delimiter=",", names=True)
    rows = []
    for V, Wv in zip(mk["V"], mk["W_um"]):
        if -3.1 < V < -0.1:
            rows.append((round(V * 4) / 4, Wv, 0.03, 0.05, "marker"))
    for V, Wv, wt in zip(band["V"], band["W_um"], band["W_weight"]):
        if not np.isfinite(Wv) or wt < 5 or V > -3.5 or V < -59.9:
            continue
        rows.append((V, Wv, 0.03 if V <= -5.0 else 0.08, 0.03, "band"))
    return _table(rows)


def _table(rows):
    rows.sort()
    V = np.array([r[0] for r in rows])
    Wd = np.array([r[1] for r in rows])
    sW = np.array([r[2] for r in rows])
    sV = np.array([r[3] for r in rows])
    src = np.array([r[4] for r in rows])
    slope = np.gradient(Wd, V)
    sig = np.sqrt(sW**2 + (slope * sV) ** 2)
    return dict(V=V, W=Wd, sig=sig, src=src)


def at(Varr, Warr, v):
    i = np.argmin(np.abs(Varr - v))
    return float(Warr[i]) if abs(Varr[i] - v) < 0.13 else float("nan")


def report(D, name, pars, Wm, ok, npar=None):
    V, Wd, sig, src = D["V"], D["W"], D["sig"], D["src"]
    res = Wd - Wm
    c2 = float(((res / sig) ** 2).sum())
    knee = (src == "marker") | (V == 0)
    plat = V <= -5
    r = dict(
        name=name,
        pars={k: float(v) for k, v in pars.items()},
        chi2=c2,
        dof=int(len(V) - (npar if npar is not None else len(pars))),
        rms_knee_um=float(np.sqrt((res[knee] ** 2).mean())),
        rms_plateau_um=float(np.sqrt((res[plat] ** 2).mean())),
        W_key={str(v): at(V, Wm, v) for v in KEY},
        resid_key={str(v): at(V, res, v) for v in KEY},
        converged=bool(ok),
    )
    print(
        f"{name:40s} chi2 {c2:9.1f} (dof {r['dof']})  rms knee {r['rms_knee_um']:.3f}  plateau {r['rms_plateau_um']:.4f} um  "
        + " ".join(f"{k}={v:.3g}" for k, v in pars.items()),
        flush=True,
    )
    return r


# ---- parametric families: p -> (NJ, NB, L_cm, kwargs)
FAMILIES = {
    # name: (extra parameter names, lower bounds, upper bounds, default starts)
    "plain": ([], [], [], []),
    "buffer": (["lLb"], [-1.3], [0.8], [[-0.28], [0.0], [0.3]]),
    "area": (["a"], [0.7], [1.4], [[1.0], [1.03]]),
    "fringe": (["Qf12"], [0.0], [2.0], [[1.0], [0.5]]),
    "buffer+Nnp": (
        ["lLb", "lNnp"],
        [-1.3, 15.0],
        [0.8, 19.5],
        [[0.0, 19.0], [0.1, 16.5], [0.3, 17.5]],
    ),
    "buffer+fringe": (
        ["lLb", "Qf12"],
        [-1.3, 0.0],
        [0.8, 2.0],
        [[0.0, 0.5], [0.3, 1.0], [-0.3, 1.5]],
    ),
    "buffer+Nnp+fringe": (
        ["lLb", "lNnp", "Qf12"],
        [-1.3, 15.0, 0.0],
        [0.8, 19.5, 2.0],
        [[0.1, 16.5, 0.5], [0.3, 17.5, 1.0], [0.0, 19.0, 0.0]],
    ),
    "buffer+Nnp+fringe0": (
        ["lLb", "lNnp"],
        [-1.3, 15.0],
        [0.8, 19.5],
        [[0.1, 16.5], [0.3, 17.5]],
    ),  # Qf = 0 fixed
    "buffer+Nnp+area": (
        ["lLb", "lNnp", "a"],
        [-1.3, 15.0, 0.7],
        [0.8, 19.5, 1.4],
        [[0.1, 16.5, 1.0], [0.3, 17.5, 0.95]],
    ),
}
BASE_STARTS = [[15.47, 13.92, -0.006], [16.0, 13.6, -0.5], [15.9, 13.5, -0.35]]


def unpack(p, family):
    NJ, NB, L = 10 ** p[0], 10 ** p[1], 10 ** p[2] * 1e-4
    names = FAMILIES[family][0]
    kw = {}
    for nm, val in zip(names, p[3:]):
        if nm == "lLb":
            kw["Lb_um"] = 10**val
        elif nm == "a":
            kw["area_scale"] = val
        elif nm == "Qf12":
            kw["Qf"] = val * 1e12
        elif nm == "lNnp":
            kw["N_np"] = 10**val
    if family.endswith("fringe0"):
        kw["Qf"] = 0.0
    return NJ, NB, L, kw


def fit_family(D, family, epi=10.0, label=None):
    names, lo, hi, starts = FAMILIES[family]
    V, Wd, sig = D["V"], D["W"], D["sig"]

    def resid(p):
        NJ, NB, L, kw = unpack(p, family)
        Wm, ok = model_W(V, NJ, NB, L, epi_um=epi, **kw)
        if not ok or not np.all(np.isfinite(Wm)):
            return np.full_like(Wd, 1e3)
        return (Wd - Wm) / sig

    best = None
    extra_starts = starts if starts else [[]]
    for b0 in BASE_STARTS:
        for e0 in extra_starts:
            r = least_squares(
                resid,
                b0 + e0,
                bounds=([14, 12, -1.3] + lo, [17.5, 15, 0.7] + hi),
                xtol=1e-10,
                ftol=1e-10,
                max_nfev=300,
            )
            c = float((r.fun**2).sum())
            if best is None or c < best[0]:
                best = (c, r.x)
    NJ, NB, L, kw = unpack(best[1], family)
    Wm, ok = model_W(V, NJ, NB, L, epi_um=epi, **kw)
    pars = dict(NJ=NJ, NB=NB, L_um=L * 1e4, epi_um=epi, **kw)
    return (
        report(D, label or f"{family} (epi {epi})", pars, Wm, ok, npar=3 + len(names)),
        Wm,
    )


# ---- model-independent piecewise log-linear profile (knots in um from the junction)
KNOTS_UM = np.array(
    [
        0.0,
        0.25,
        0.5,
        0.75,
        1.0,
        1.5,
        2.0,
        3.0,
        4.0,
        5.0,
        6.0,
        7.0,
        8.0,
        8.5,
        9.0,
        9.25,
        9.5,
        9.75,
        10.0,
    ]
)


def knot_profile(x, logN_knots, epi_um=10.0):
    xu = (x - SUB) * 1e4
    logN = np.interp(np.clip(xu, 0, KNOTS_UM[-1]), KNOTS_UM, logN_knots)
    d = 10**logN * (x >= SUB)
    # beyond the last knot: hold the last value (the C-V never probes past ~9.9 um)
    return d


def fit_knots(D, lam=1.0, Qf=None, area_scale=1.0, p0=None):
    V, Wd, sig = D["V"], D["W"], D["sig"]
    x = node_positions(10e-4, SUB, NP)

    def model(p):
        pb = batch_from_profiles(x, knot_profile(x, p)[None, :])
        Wc, conv = sweep_W(pb, V)
        C = EPS / Wc[:, 0]
        if Qf is not None:
            C = C + fringe_Cprime(V, Qf) * P_OVER_A
        return EPS / (C * area_scale) * 1e4, conv.all()

    def resid(p):
        Wm, ok = model(p)
        if not ok or not np.all(np.isfinite(Wm)):
            return np.concatenate([np.full_like(Wd, 1e3), np.zeros(len(p) - 2)])
        smooth = (
            lam * np.diff(p, 2) / np.diff(KNOTS_UM)[1:] ** 0.5
        )  # curvature penalty, weaker where knots are dense
        return np.concatenate([(Wd - Wm) / sig, smooth])

    if p0 is None:
        p0 = np.interp(
            KNOTS_UM,
            [0, 0.5, 1.5, 6, 8.5, 9.5, 10],
            [16.0, 15.3, 13.6, 13.6, 14.0, 15.5, 16.5],
        )
    r = least_squares(
        resid, p0, bounds=(12.5, 19.5), xtol=1e-9, ftol=1e-9, max_nfev=200
    )
    Wm, ok = model(r.x)
    return r.x, Wm, ok
