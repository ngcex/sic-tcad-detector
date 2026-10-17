"""Robustness of the N_D preference of the (d) columnar and (b-trap) fits (hypotheses.py)
to the generation model: shape bragg/box, tilt 60/65/70 deg, range 17/18/19 um, finer
N_D grid. Reports chi2_min(N_D) per variant (free scale, sigma 1e-3) and the N_D range
within Delta chi2 <= 9 of the variant's best (3 sigma for one interesting parameter).

usage: uv run python lab/microdosimeter/bench/hypotheses_robust.py > lab/microdosimeter/notes/hypotheses_robust.out
"""

from __future__ import annotations

import sys
import time

import numpy as np

sys.path.insert(0, __file__.rsplit("/bench/", 1)[0])
sys.path.insert(0, __file__.rsplit("/hypotheses_robust.py", 1)[0])
import mdlab  # noqa: E402
from hypotheses import (  # noqa: E402
    DATA,
    PEAK,
    SIG,
    V_DATA,
    V_SAT_N,
    V_SAT_P,
    dark_fields,
    hecht_curve,
)

from etna.core.batched.material import PARAMS  # noqa: E402
from etna.core.sic_material import mobility_caughey_thomas_T  # noqa: E402

N_FINE = np.array([0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0575, 4.0, 6.0]) * 1e14
SUB = 0.3e-4
taus = np.geomspace(1e-11, 1e-5, 37)
a_grid = np.linspace(0.0, 0.8, 81)
Ec_grid = np.geomspace(1e2, 1e6, 81)


def fit_variant(shape, theta, R_um, td=0.0):
    x = mdlab.mesh(SUB)
    pb, _, fld, _ = dark_fields(SUB, N_FINE, V_DATA)
    dx = np.diff(x)
    xc = 0.5 * (x[1:] + x[:-1])
    G, _ = mdlab.tilted_generation(
        x,
        theta_deg=theta,
        t_dead_cm=td,
        R_cm=R_um * 1e-4,
        shape=shape,
        epi_end=SUB + mdlab.EPI_CM,
    )
    G[x < SUB] = 0.0
    Gc = 0.5 * (G[1:] + G[:-1]) * dx
    Gtot = Gc.sum()
    ccedd, _ = mdlab.cce_sweep_custom(
        pb, SUB, V_DATA, G / G.max() * PEAK, N_D_mob=N_FINE, hole_mob_from="N_D"
    )
    out = []
    for b, N in enumerate(N_FINE):
        mu_h = mobility_caughey_thomas_T(N, 300.0, "hole", PARAMS)
        mu_e = mobility_caughey_thomas_T(N, 300.0, "electron", PARAMS)
        Qh_all, Qe_all, cce_c = (
            [],
            [],
            np.empty((len(V_DATA), len(a_grid), len(Ec_grid))),
        )
        for k, V in enumerate(V_DATA):
            E, W, n = fld[k]
            Eb, Wb = E[b], W[b]
            i_j = int(np.searchsorted(xc, SUB))
            i_w = max(int(np.searchsorted(xc, SUB + Wb)) - 1, i_j)
            dep = np.zeros(len(xc), bool)
            dep[i_j : i_w + 1] = True
            G_dep, G_und = Gc[dep].sum(), Gc[~dep & (xc >= SUB)].sum()
            eta = (
                1.0
                if G_und <= 0
                else float(np.clip((ccedd[k, b] * Gtot - G_dep) / G_und, 0, 1))
            )
            vh = mu_h * Eb / (1 + mu_h * Eb / V_SAT_P) + 1e-30
            ve = mu_e * Eb / (1 + mu_e * Eb / V_SAT_N) + 1e-30
            Qh, Qe = hecht_curve(dx, vh, ve, Gc, Wb, i_j, i_w, taus, taus, eta, G_und)
            Qh_all.append(Qh)
            Qe_all.append(Qe)
            g, e = Gc[dep], Eb[dep]
            f = 1.0 / (1.0 + e[None, :] / Ec_grid[:, None])
            dep_coll = g.sum() - a_grid[:, None] * (f * g[None, :]).sum(1)[None, :]
            cce_c[k] = (dep_coll + eta * G_und * (1.0 - a_grid)[:, None]) / Gtot
        Qh_all, Qe_all = np.array(Qh_all), np.array(Qe_all)
        res = {}
        for name, cube in (
            ("trap", (Qh_all[:, :, None] + Qe_all[:, None, :]) / Gtot),
            ("col", cce_c),
            ("dd", ccedd[:, b][:, None, None]),
        ):
            M = cube.reshape(len(V_DATA), -1).T
            s = (M @ DATA) / np.einsum("ij,ij->i", M, M)
            c2 = (((s[:, None] * M - DATA) / SIG) ** 2).sum(1)
            ib = int(np.argmin(c2))
            if name == "trap":
                ih, ie = divmod(ib, len(taus))
                par = f"tau_h {taus[ih]:.1e} tau_e {taus[ie]:.1e}"
            elif name == "col":
                ia, iE = divmod(ib, len(Ec_grid))
                par = f"a {a_grid[ia]:.3f} E_c {Ec_grid[iE]:.2e}"
            else:
                par = "tau x1"
            res[name] = (float(c2[ib]), par, float(s[ib]))
        out.append((N, res))
    return out


def main():
    t0 = time.time()
    variants = [
        ("bragg", 65.0, 18.0, 0.0),
        ("box", 65.0, 18.0, 0.0),
        ("bragg", 60.0, 18.0, 0.0),
        ("bragg", 70.0, 18.0, 0.0),
        ("bragg", 65.0, 17.0, 0.0),
        ("bragg", 65.0, 19.0, 0.0),
        ("bragg", 65.0, 18.0, 0.3e-4),
    ]
    print(
        "chi2_min per N_D (free scale, sigma 1e-3); columns: (d) columnar | (b) trapping | DD only"
    )
    for shape, theta, R, td in variants:
        out = fit_variant(shape, theta, R, td)
        print(
            f"\n== {shape}, theta {theta:.0f}, R {R:.0f} um, t_dead {td * 1e4:.1f} um =="
        )
        for N, res in out:
            print(
                f"  N_D {N:.3e}: col {res['col'][0]:8.1f} ({res['col'][1]}, s {res['col'][2]:.4f}) | "
                f"trap {res['trap'][0]:8.1f} ({res['trap'][1]}) | dd {res['dd'][0]:9.1f}"
            )
        for name in ("col", "trap"):
            c2 = np.array([r[name][0] for _, r in out])
            ok = N_FINE[c2 <= c2.min() + 9.0]
            print(
                f"  {name}: best N_D {N_FINE[np.argmin(c2)]:.2e} (chi2 {c2.min():.1f}); "
                f"within Dchi2<=9: {ok.min():.2e} .. {ok.max():.2e}"
            )
    print(f"\ntotal {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
