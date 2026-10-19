"""Phase-4 recalibration of the n+ model: C-V (graded doping) then N_t.

C-V: three depletion widths (W = eps A / C, A = 4 mm2) at 0 / -10 / -30 V:
1.7 / 9.5 / 9.73 um, three unknowns (N_D_junction, N_D_bulk, L_transition):
an EXACT SOLVE (3 equations, 3 unknowns), not a fit -- no goodness of fit exists.
Solved with the batched engine (W parity vs devsim 6e-15 with n+) by least_squares on
log-parameters from several starting points (uniqueness check), then verified on devsim.

N_t: leakage is linear in N_t; one devsim dark_current_sweep at the reference N_t, rescale
to the bound "<18 pA up to -60 V at 4 mm2 (300 K)" used as equality, verify.

usage: uv run python lab/solver2d/bench/nplus_calibrate.py [cv|nt|all|cal-c] [--no-nplus]
(cal-c: the 2026-10-17 recalibration on eps/C of the small-signal C, calibrate_on_C)
"""

from __future__ import annotations

import sys
import uuid

import numpy as np
import scipy.optimize

EPI = 10e-4
SUB = 0.3e-4  # p+ anode layer (paper Fig. 1: 0.3 um)
NPLUS = 1e-4  # n+ layer thickness (truncated substrate; results independent of it)
N_NPLUS = 1e19
V_T = np.array([0.0, -10.0, -30.0])
W_T = np.array([1.7e-4, 9.5e-4, 9.73e-4])
AREA = 0.04  # cm2 (4 mm2)
I_BOUND = 18e-12  # A at -60 V, 300 K
E_T = 0.67  # eV, Z1/2 (radiation_damage.DamageParameters.E_Z12, Burin et al. Table I)


def engine_W(params, nplus=True, sub=SUB, observable="edge"):
    """W at V_T: observable "edge" = n/N_D = 0.5 depletion edge (the pre-2026-10-17
    calibration), "C" = eps/C of the small-signal C = dQ/dV (what W = eps A / C of a
    measured C is; the calibration since 2026-10-17)."""
    from etna.core.batched import make_batch, node_positions, reverse_sweep
    from etna.core.batched.material import EPS

    NJ, NB, L = params
    x = node_positions(EPI, sub, NPLUS if nplus else None)
    pb = make_batch(
        x, sub, 1e19, NJ, NB, L, 300.0, nplus=(sub + EPI, N_NPLUS) if nplus else None, buffer=None
    )
    r = reverse_sweep(pb, V_T, sub, EPI)
    if not r.converged.all():
        return np.full(3, np.nan)
    return r.W[:, 0] if observable == "edge" else EPS / r.C[:, 0]


def solve_cv(nplus=True, sub=SUB):
    def res(z):
        W = engine_W(np.exp(z), nplus, sub)
        return (W - W_T) / W_T

    lo = np.log([1e14, 1e12, 0.05e-4])
    hi = np.log([1e17, 1e15, 5e-4])
    starts = [
        (2.9e15, 8.5e13, 1e-4),
        (2.93e15, 8.82e13, 0.987e-4),
        (1e16, 5e13, 0.5e-4),
        (1e15, 2e14, 2e-4),
        (5e15, 1e14, 0.3e-4),
    ]
    sols = []
    for s in starts:
        r = scipy.optimize.least_squares(
            res,
            np.log(s),
            bounds=(lo, hi),
            xtol=1e-15,
            ftol=1e-15,
            gtol=1e-15,
            max_nfev=2000,
        )
        p = np.exp(r.x)
        sols.append((s, p, res(r.x), r.status))
        print(
            f"RES cv nplus={nplus} sub={sub*1e4:.2f}um start {s} -> NJ {p[0]:.6e} NB {p[1]:.6e} L {p[2]*1e4:.6f} um "
            f"residuals {np.array2string(res(r.x), precision=2)} status {r.status}",
            flush=True,
        )
    return sols


def devsim_W(params, nplus=True, sub=SUB):
    import devsim

    from etna.core.cv_analysis import cv_sweep
    from etna.core.drift_diffusion import create_dd_device

    NJ, NB, L = params
    kw = dict(
        doping_profile="graded",
        N_D_junction=NJ,
        N_D_bulk=NB,
        L_transition=L,
        substrate_thickness_cm=sub,
    )
    if nplus:
        kw.update(nplus_thickness_cm=NPLUS, N_D_nplus=N_NPLUS)
    name = f"cal_{uuid.uuid4().hex[:6]}"
    info = create_dd_device(device_name=name, **kw)
    r = cv_sweep(info, V_T)
    devsim.delete_device(device=name)
    return np.asarray(r["depletion_widths"])


def dark_60(params, N_t, nplus=True, sub=SUB, E_t=E_T, V=(-10.0, -30.0, -60.0)):
    import devsim

    from etna.core.dark_current import (
        create_dark_current_device,
        dark_current_sweep,
        setup_tat_model,
    )

    NJ, NB, L = params
    kw = dict(
        doping_profile="graded",
        N_D_junction=NJ,
        N_D_bulk=NB,
        L_transition=L,
        substrate_thickness_cm=sub,
    )
    if nplus:
        kw.update(nplus_thickness_cm=NPLUS, N_D_nplus=N_NPLUS)
    name = f"dk_{uuid.uuid4().hex[:6]}"
    info = create_dark_current_device(N_t=N_t, device_name=name, **kw)
    setup_tat_model(
        info, E_t=E_t, N_t=N_t
    )  # re-register with the reconciled trap level
    r = dark_current_sweep(info, np.array(V), area=AREA)
    devsim.delete_device(device=name)
    return np.abs(np.asarray(r["I_total"]))


if __name__ == "__main__":
    task = sys.argv[1] if len(sys.argv) > 1 else "all"
    if task == "cal-c":
        pass  # dispatched at the end of the file (calibrate_on_C is defined below)
    nplus = "--no-nplus" not in sys.argv
    if task in ("cv", "all") and task != "cal-c":
        sols = solve_cv(nplus)
        best = min(sols, key=lambda s: np.nansum(s[2] ** 2))[1]
        Wd = devsim_W(best, nplus)
        print(
            f"RES cv devsim check W_um {np.array2string(Wd*1e4, precision=6)} rel {np.array2string((Wd-W_T)/W_T, precision=2)}"
        )
        print(f"RES CV_SOLUTION {best[0]:.6e} {best[1]:.6e} {best[2]:.6e}")
    if task in ("nt", "all"):
        params = [float(v) for v in sys.argv[2:5]] if task == "nt" else best
        I_ref = dark_60(params, 2.8e12, nplus)
        N_t = 2.8e12 * I_BOUND / I_ref[-1]
        I_new = dark_60(params, N_t, nplus)
        print(
            f"RES nt I(N_t=2.8e12) {np.array2string(I_ref, precision=6)} A -> N_t {N_t:.6e} -> I {np.array2string(I_new, precision=6)} A"
        )


# --- 2-point solve (0 V, -10 V) at fixed L: the -30 V target is unreachable ----------------
def solve_2pt(L, nplus=True, sub=SUB, x0=(2.9e15, 8.5e13), observable="edge"):
    """Exact solve of W(0 V) = 1.7 um and W(-10 V) = 9.5 um for (N_D_junction, N_D_bulk) at
    fixed L_transition; returns params and the residuals at 0 / -10 / -30 V."""

    def res(z):
        W = engine_W((np.exp(z[0]), np.exp(z[1]), L), nplus, sub, observable)
        return ((W - W_T) / W_T)[:2]

    r = scipy.optimize.least_squares(
        res, np.log(x0), bounds=(np.log([1e14, 1e12]), np.log([1e17, 1e15])), xtol=1e-15, ftol=1e-15, gtol=1e-15
    )
    p = (float(np.exp(r.x[0])), float(np.exp(r.x[1])), L)
    W = engine_W(p, nplus, sub, observable)
    return p, (W - W_T) / W_T


def devsim_W_C(params, nplus=True, sub=SUB):
    """eps/C at V_T on the devsim path (cv_sweep small-signal C, W_from_C)."""
    import devsim

    from etna.core.cv_analysis import cv_sweep
    from etna.core.drift_diffusion import create_dd_device

    NJ, NB, L = params
    kw = dict(doping_profile="graded", N_D_junction=NJ, N_D_bulk=NB, L_transition=L,
              substrate_thickness_cm=sub)
    if nplus:
        kw.update(nplus_thickness_cm=NPLUS, N_D_nplus=N_NPLUS)
    name = f"calc_{uuid.uuid4().hex[:6]}"
    info = create_dd_device(device_name=name, **kw)
    r = cv_sweep(info, V_T)
    devsim.delete_device(device=name)
    return np.asarray(r["W_from_C"]), np.asarray(r["depletion_widths"])


def calibrate_on_C(L=0.987e-4, N_t_ref=2.8e12):
    """2026-10-17 recalibration on the small-signal observable: (N_D_junction, N_D_bulk)
    exact on eps/C at 0 / -10 V (L held), verified on devsim, then N_t rescaled to
    I(300 K, -60 V, 4 mm2) = 18 pA and verified."""
    p, rel = solve_2pt(L, observable="C", x0=(2.974271e15, 6.335936e13))
    print(f"RES C-cal NJ {p[0]:.6e} NB {p[1]:.6e} L {p[2]*1e4:.4f} um | engine eps/C residuals "
          f"0/-10/-30 V {np.array2string(rel, precision=3)}", flush=True)
    WC, We = devsim_W_C(p)
    print(f"RES C-cal devsim eps/C {np.array2string(WC*1e4, precision=6)} um rel "
          f"{np.array2string((WC-W_T)/W_T, precision=3)}; W_edge {np.array2string(We*1e4, precision=4)} um",
          flush=True)
    I_ref = dark_60(p, N_t_ref)
    N_t = N_t_ref * I_BOUND / I_ref[-1]
    I_new = dark_60(p, N_t)
    print(f"RES C-cal N_t {N_t_ref:.4g} -> I {np.array2string(I_ref, precision=6)} A; N_t {N_t:.6e} -> "
          f"I(-10,-30,-60 V) {np.array2string(I_new, precision=6)} A", flush=True)
    return p, N_t


def cce_engine(params, nplus=True, sub=SUB, V=(-10.0, -50.0)):
    from etna.core.batched import alpha_generation, cce_sweep, make_batch, node_positions

    NJ, NB, L = params
    x = node_positions(EPI, sub, NPLUS if nplus else None)
    pb = make_batch(x, sub, 1e19, NJ, NB, L, 300.0, nplus=(sub + EPI, N_NPLUS) if nplus else None, buffer=None)
    G = alpha_generation(x, sub, 15e-4, 1e18, epi_end=sub + EPI if nplus else None)
    return cce_sweep(pb, np.array(V), sub, G).cce[:, 0]


def scan_L(nplus=True, sub=SUB):
    for L_um in (0.5, 0.75, 0.987, 1.25, 1.5, 2.0):
        p, rel = solve_2pt(L_um * 1e-4, nplus, sub)
        c = cce_engine(p, nplus, sub)
        d = dark_60(p, 2.8e12, nplus, sub)
        print(
            f"RES scanL L {L_um:.3f} um -> NJ {p[0]:.6e} NB {p[1]:.6e} | residuals 0/-10/-30 V "
            f"{np.array2string(rel, precision=2)} | CCE(-10,-50) {np.array2string(c, precision=6)} | "
            f"I(N_t 2.8e12; -10,-30,-60 V) {np.array2string(d, precision=5)} A",
            flush=True,
        )


if __name__ == "__main__" and sys.argv[1:2] == ["cal-c"]:
    calibrate_on_C()
