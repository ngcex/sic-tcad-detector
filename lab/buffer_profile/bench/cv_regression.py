"""Measured-C-V preset vs the digitised Photons Fig. 4b curve, through run_cv (batched).

Gate: plateau (-5..-59 V, 250 mV grid) rms <= 0.02 um; knee markers (-3.25..-0.25 V)
within the digitisation uncertainty sigma = 0.03 um (+) 0.04 V x |dW/dV|; labels at
0 / -10 / -30 V. Also a mesh-resolution check of the transition region (production mesh
vs a 4x refined mesh, batched engine on both).

usage: uv run python lab/buffer_profile/bench/cv_regression.py
"""

from __future__ import annotations

import numpy as np

import etna
from etna.api.simulation import run_cv
from etna.core.batched.material import EPS

ROOT = "lab/cv_gap/notes/"


def digitised():
    band = np.genfromtxt(ROOT + "photons_w.csv", delimiter=",", names=True)
    mk = np.genfromtxt(ROOT + "photons_w_markers.csv", delimiter=",", names=True)
    plateau = [
        (float(v), float(w))
        for v, w, wt in zip(band["V"], band["W_um"], band["weight"])
        if np.isfinite(w)
        and wt >= 10
        and v <= -5.0
        and v not in (-10.0, -30.0, -9.75, -10.25, -10.5, -9.5, -29.75, -30.25)
    ]
    knee = [
        (round(float(v) * 4) / 4, float(w))
        for v, w in zip(mk["V"], mk["W_um"])
        if -3.3 < v < -0.1
    ]
    labels = [(0.0, 1.68), (-10.0, 9.52), (-30.0, 9.73)]
    return plateau, knee, labels


def model_W(config, V):
    V = np.asarray(V, float)
    order = np.argsort(-V)  # 0 first, increasingly negative
    Vs = V[order]
    # run_cv solves the grid it is given; use the sorted grid, map back
    r = run_cv(config, v_start=float(Vs[0]), v_stop=float(Vs[-1]), n_points=len(Vs))
    assert r.metadata["solver"] == "batched", r.metadata
    assert np.allclose(r.x, Vs)
    W = EPS * config.area_cm2 / r.y * 1e4
    out = np.empty_like(V)
    out[order] = W
    return out, r


def main():
    cfg = etna.measured_cv_config(area_cm2=0.04)
    plateau, knee, labels = digitised()
    Vp = np.array([v for v, _ in plateau])
    Wp = np.array([w for _, w in plateau])
    Vk = np.array([v for v, _ in knee])
    Wk = np.array([w for _, w in knee])
    # one uniform 250 mV grid 0..-60 covers everything
    V = np.round(np.arange(0.0, -60.001, -0.25), 2)
    W, r = model_W(cfg, V)
    at = lambda v: W[np.argmin(np.abs(V - v))]
    res_p = Wp - np.array([at(v) for v in Vp])
    rms_p = float(np.sqrt(np.mean(res_p**2)))
    slope = np.gradient(W, V)
    sig_k = np.sqrt(
        0.03**2 + (0.04 * np.array([slope[np.argmin(np.abs(V - v))] for v in Vk])) ** 2
    )
    res_k = Wk - np.array([at(v) for v in Vk])
    print(
        f"RES plateau N={len(Vp)} rms {rms_p:.4f} um max|res| {np.abs(res_p).max():.4f} um (gate <= 0.02)"
    )
    print(
        f"RES knee N={len(Vk)} rms {np.sqrt(np.mean(res_k**2)):.3f} um; max |res|/sigma {np.max(np.abs(res_k)/sig_k):.2f}; within 1 sigma {int(np.sum(np.abs(res_k) <= sig_k))}/{len(Vk)}, within 2 sigma {int(np.sum(np.abs(res_k) <= 2*sig_k))}/{len(Vk)}"
    )
    for v, w, rr, s in zip(Vk, Wk, res_k, sig_k):
        print(
            f"    knee V {v:6.2f}  data {w:.3f}  model {w - rr:.3f}  res {rr:+.3f}  sigma {s:.3f}"
        )
    for v, w in labels:
        print(
            f"RES label V {v:6.1f} data {w:.3f} model {at(v):.4f} res {w - at(v):+.4f}"
        )
    chi2 = float(
        np.sum((res_p / 0.02) ** 2)
        + np.sum((res_k / sig_k) ** 2)
        + sum(((w - at(v)) / (0.02 if v == 0 else 0.01)) ** 2 for v, w in labels)
    )
    print(
        f"RES chi2 {chi2:.1f} over {len(Vp) + len(Vk) + 3} points (lab: 198 / 230 with its band sigmas)"
    )
    print(
        f"RES reach_through_V {r.metadata['reach_through_V']}  W(-2) {at(-2):.3f} W(-5) {at(-5):.3f} W(-58) {at(-58):.3f}"
    )

    # ---- mesh resolution of the transition region: production mesh vs 4x refined
    from etna.core.batched import make_batch, node_positions, reverse_sweep
    from etna.core.device import MEASURED_CV_STACK as M

    sub, epi, t_np = M["substrate_thickness_cm"], 10e-4, M["nplus_thickness_cm"]
    x = node_positions(epi, sub, t_np)
    xr = np.sort(
        np.concatenate(
            [
                x,
                (x[:-1] + x[1:]) / 2,
                (3 * x[:-1] + x[1:]) / 4,
                (x[:-1] + 3 * x[1:]) / 4,
            ]
        )
    )
    Vt = np.array([0.0, -2.0, -5.0, -10.0, -30.0, -60.0])
    for lbl, xx in (("production", x), ("4x refined", xr)):
        pb = make_batch(
            xx,
            sub,
            1e19,
            M["N_D_junction"],
            M["N_D_bulk"],
            M["L_transition"],
            300.0,
            nplus=(sub + epi, M["N_D_nplus"]),
            buffer=(sub + epi, M["buffer_thickness_cm"], M["buffer_end_doping"]),
        )
        rr = reverse_sweep(pb, Vt, sub, epi)
        Wc = EPS / rr.C[:, 0] * 1e4
        print(
            f"RES mesh {lbl:11s} N={len(xx)} eps/C um: "
            + " ".join(f"{v:g}V={w:.5f}" for v, w in zip(Vt, Wc))
        )
        if lbl == "production":
            W0 = Wc
        else:
            print(
                f"RES mesh refined-production max |dW| {np.max(np.abs(Wc - W0)):.2e} um, rel {np.max(np.abs(Wc - W0) / W0):.2e}"
            )
    m = (x > sub + epi - M["buffer_thickness_cm"]) & (x < sub + epi)
    print(
        f"RES mesh nodes in the transition {int(m.sum())}, spacing {np.diff(x[m]).min()*1e4:.4f}..{np.diff(x[m]).max()*1e4:.4f} um"
    )


if __name__ == "__main__":
    main()
