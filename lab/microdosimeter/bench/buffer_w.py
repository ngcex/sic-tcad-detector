"""Task 3: does a graded n-/n+ buffer reproduce the dosimeter's W creep between -10 and
-30 V (Photons Fig. 4b labels 1.68 / 9.52 / 9.73 um at 0 / -10 / -30 V; production model
1.700 / 9.500 / 9.899 um)?

Buffer: donors rise log-linearly from the local epi value to N_D_nplus over a width Lb,
ending at (x_end + shift): shift = 0 puts the whole transition inside the nominal 10 um epi,
shift = Lb/2 centres it on the interface. For every buffer the two epi parameters
(N_D_junction, N_D_bulk) are re-solved EXACTLY on the 0 V and -10 V targets (L_transition
held at 0.987 um, as the production calibration), then W(-30 V) is reported.

Two W definitions:
  W_edge  production: first n/N_D >= 0.5 node (what the calibration fits and what
          capacitance() inverts, C = eps A / W_edge);
  W_C     eps / (dQ/dV): what eps*A/C of a MEASURED C is. Small-signal: the AC current
          through a carrier-free plane x0 in the depletion region is pure displacement
          current, so C = eps dE(x0)/dV (central difference +-DV). x0 = 0.5 um below the
          junction: past the p+ hole Debye tail (~5e-8 C/cm2 that SHRINKS with bias; the
          junction field peak or the integrated n-side charge mix it in and give nonsense,
          W_C ~ 20 um) and inside W at every bias used (W >= 1.7 um).

usage: uv run python lab/microdosimeter/bench/buffer_w.py
"""

import numpy as np
import scipy.optimize

from etna.core.batched import make_batch, node_positions
from etna.core.batched.extract import depletion_width
from etna.core.batched.material import EPS, Q
from etna.core.batched.poisson import PoissonBatch
from etna.core.batched.sweeps import _Ramp
from etna.core.device import DEFAULT_STACK as S

EPI0, SUB, NP = 10e-4, S["substrate_thickness_cm"], S["nplus_thickness_cm"]
V_T = np.array([0.0, -10.0, -30.0])
W_PAPER = np.array([1.68, 9.52, 9.73])  # um, Photons Fig. 4b labels
W_TGT = np.array([1.7e-4, 9.5e-4, 9.73e-4])  # the production calibration's targets
DV = 0.25
X0 = 0.5e-4


def batch(NJ, NB, Lb_cm=0.0, shift_cm=0.0, epi=EPI0, L=S["L_transition"]):
    x = node_positions(epi, SUB, NP)
    x_end = SUB + epi
    pb = make_batch(x, SUB, 1e19, NJ, NB, L, 300.0, nplus=(x_end, S["N_D_nplus"]))
    if Lb_cm > 0:
        donors = pb.donors.copy()
        x1 = x_end + shift_cm
        x0 = x1 - Lb_cm
        base = NB + (NJ - NB) * np.exp(-(x - SUB) / L)
        t = np.clip((x - x0) / Lb_cm, 0.0, 1.0)
        graded = base * (S["N_D_nplus"] / base) ** t
        m = (x >= x0) & (x <= x1)
        donors[0, m] = np.maximum(donors[0, m], graded[m])
        pb = PoissonBatch(x, donors, pb.donors - pb.netdop, pb.n_i, pb.V_t)
    return pb, epi


def widths(pbe, x0=X0, dv=DV, V_T=V_T):
    """W_edge and W_C at V_T. Central differences everywhere (at 0 V the +dv point is a
    small FORWARD bias, fine for this dark Poisson model)."""
    pb, epi = pbe
    x = pb.x
    volts = sorted({round(v + d, 9) for v in V_T for d in (-dv, 0.0, dv)}, reverse=True)
    ramp = _Ramp(pb, 1e-10)
    xm = 0.5 * (x[1:] + x[:-1])
    out = {}
    for V in volts:
        res = (
            ramp.to(-V)
            if (abs(V) > 1e-12 or ramp.Vc != 0.0)
            else pb.solve(ramp.psi, tol=1e-10)
        )
        E = -np.diff(res.psi[0]) / np.diff(x)
        out[V] = (
            depletion_width(x, res.n, pb.donors, SUB, epi + NP)[0],
            EPS * abs(np.interp(SUB + x0, xm, E)),
        )
    We = np.array([out[round(v, 9)][0] for v in V_T])
    WC = []
    for v in V_T:
        lo, hi = round(v - dv, 9), round(v + dv, 9)
        WC.append(EPS * 2 * dv / abs(out[lo][1] - out[hi][1]))
    return We, np.array(WC), bool(ramp.converged.all())


def resolve(mode, n_tgt=2, **kw):
    """Exact solve of (NJ, NB) [n_tgt=2] or (NJ, NB, L) [n_tgt=3] on the W targets."""

    def res(z):
        p = np.exp(z)
        extra = dict(L=p[2]) if n_tgt == 3 else {}
        We, WC, ok = widths(batch(p[0], p[1], **kw, **extra))
        W = We if mode == "edge" else WC
        return (
            (W[:n_tgt] - W_TGT[:n_tgt]) / W_TGT[:n_tgt] if ok else np.full(n_tgt, 1e3)
        )

    z0 = np.log([S["N_D_junction"], S["N_D_bulk"], S["L_transition"]][:n_tgt])
    lo = np.log([1e14, 1e12, 0.05e-4][:n_tgt])
    hi = np.log([1e17, 1e15, 5e-4][:n_tgt])
    r = scipy.optimize.least_squares(
        res, z0, bounds=(lo, hi), xtol=1e-14, ftol=1e-14, gtol=1e-14
    )
    p = np.exp(r.x)
    return p, res(r.x)


def row(label, p, r, pbe):
    We, WC, ok = widths(pbe)
    L = f" L {p[2]*1e4:.3f}" if len(p) == 3 else ""
    print(
        f"{label:<30} | NJ {p[0]:.4e} NB {p[1]:.4e}{L} res {np.abs(r).max():.0e} | "
        f"edge {We[0]*1e4:.3f} {We[1]*1e4:.3f} {We[2]*1e4:.3f} | C {WC[0]*1e4:.3f} {WC[1]*1e4:.3f} {WC[2]*1e4:.3f} | "
        f"creep e/C {(We[2]-We[1])*1e4:+.3f}/{(WC[2]-WC[1])*1e4:+.3f} conv {ok}"
    )


print(
    f"targets W (um): paper labels {W_PAPER}; production targets {W_TGT*1e4}; paper creep -10->-30 = +0.21"
)
print("== baseline (production DEFAULT_STACK, no buffer) ==")
base = batch(S["N_D_junction"], S["N_D_bulk"])
row("DEFAULT_STACK", [S["N_D_junction"], S["N_D_bulk"]], [0.0], base)
print("   W_C robustness (x0, DV):")
for x0 in (0.3e-4, 0.5e-4, 1.0e-4):
    for dv in (0.1, 0.25):
        We, WC, ok = widths(base, x0, dv)
        print(f"     x0 {x0*1e4:.1f} um DV {dv}: W_C {np.round(WC*1e4, 3)}")
u = batch(3e14, 3e14)
We, WC, _ = widths(u)
print(f"   uniform 3e14 check: W_edge {np.round(We*1e4, 3)}  W_C {np.round(WC*1e4, 3)}")

for mode in ("edge", "C"):
    print(f"== buffer scan, (NJ, NB) re-solved on the 0/-10 V targets with W_{mode} ==")
    for pos in ("inside", "centred"):
        for Lb_nm in (0, 100, 200, 300, 500, 800):
            Lb = Lb_nm * 1e-7
            shift = 0.0 if pos == "inside" else Lb / 2
            if Lb_nm == 0 and pos == "centred":
                continue
            p, r = resolve(mode, Lb_cm=Lb, shift_cm=shift)
            row(
                f"{mode} Lb {Lb_nm:4d} nm {pos if Lb_nm else 'none'}",
                p,
                r,
                batch(*p, Lb, shift),
            )
print("== 3-parameter exact solve (NJ, NB, L) on all three targets ==")
for mode in ("edge", "C"):
    p, r = resolve(mode, n_tgt=3)
    row(f"{mode} 3-par, no buffer", p, r, batch(p[0], p[1], L=p[2]))
print("== epi thickness (nominal 10 um; real epi tolerance is typically several %) ==")
for epi_um in (9.70, 9.75, 9.80):
    for mode in ("edge", "C"):
        p, r = resolve(mode, epi=epi_um * 1e-4)
        row(f"{mode} epi {epi_um} um", p, r, batch(*p, epi=epi_um * 1e-4))


def resolve_lb(mode):
    """Exact solve of (NJ, NB, Lb) on all three targets, buffer inside the nominal epi."""

    def res(z):
        p = np.exp(z)
        We, WC, ok = widths(batch(p[0], p[1], Lb_cm=p[2]))
        W = We if mode == "edge" else WC
        return (W - W_TGT) / W_TGT if ok else np.full(3, 1e3)

    z0 = np.log([S["N_D_junction"], S["N_D_bulk"], 0.5e-4])
    r = scipy.optimize.least_squares(
        res,
        z0,
        bounds=(np.log([1e14, 1e12, 0.01e-4]), np.log([1e17, 1e15, 2e-4])),
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
    )
    return np.exp(r.x), res(r.x)


print(
    "== 3-parameter exact solve (NJ, NB, Lb), L held 0.987 um, buffer inside the epi =="
)
curves = {"production (no buffer)": base}
for mode in ("edge", "C"):
    p, r = resolve_lb(mode)
    print(f"{mode}: Lb = {p[2]*1e7:.0f} nm")
    row(f"{mode} 3-par (NJ,NB,Lb)", p, r, batch(p[0], p[1], Lb_cm=p[2]))
    curves[f"buffer {p[2]*1e7:.0f} nm, fit W_{mode}"] = batch(p[0], p[1], Lb_cm=p[2])
pC, _ = resolve("C")
curves["no buffer, fit W_C 0/-10 V"] = batch(*pC)
VV = np.array([0.0, -2.0, -5.0, -10.0, -15.0, -20.0, -30.0, -40.0, -60.0])
print("== W(V) curves (um), W_edge / W_C ==")
print(f"{'':<32}" + "".join(f"{v:>14.0f}" for v in VV))
for k, pbe in curves.items():
    We, WC, ok = widths(pbe, V_T=VV)
    print(
        f"{k:<32}"
        + "".join(f"{a*1e4:7.3f}/{b*1e4:6.3f}" for a, b in zip(We, WC))
        + f"  conv {ok}"
    )
