"""Which mechanism explains Petringa 2025 Fig. 4b (CCE vs bias, 65-deg 5.5 MeV alpha)?

Lab only. Every hypothesis is judged by chi-square against the digitised curve
(mdlab.FIG4B_*), sigma = 0.001 per point (quote x4 / x0.25 for 0.0005 / 0.002), with a
FREE overall scale (the paper's absolute "100 %" rests on a Si-detector calibration and
8.4 eV/pair; only the shape is trusted). dof = 7 - 1 - n_free_params.

(a)/(b-SRH)/(c): steady-state drift-diffusion (batched engine) on a grid of uniform
    N_D x p+ x dead layer x SRH lifetime scale. The DD model has no trapping in the drift
    region (midgap SRH needs both carriers), so once W passes the track end it is flat.
(b-trap): Hecht-type trapping during drift, post-processed on the engine's dark field
    profile E(x; V): hole/electron trapping times tau_h / tau_e, v = mu E / (1 + mu E / v_sat),
    uniform weighting field 1/W over the depleted region; charge generated beyond W is
    collected with the DD diffusion efficiency and then drifts the full W.
(d): columnar / plasma recombination, loss(E) = a / (1 + E / E_c) on the local field at
    the generation point (Jaffe-like saturating form; a = field-free loss of the dense
    column, E_c = field scale), the same a for the undepleted part of the track.

usage: uv run python lab/microdosimeter/bench/hypotheses.py > lab/microdosimeter/notes/hypotheses.out
"""

from __future__ import annotations

import itertools
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, __file__.rsplit("/bench/", 1)[0])
import mdlab  # noqa: E402
from etna.core.batched.extract import depletion_width  # noqa: E402
from etna.core.batched.material import PARAMS, Q  # noqa: E402
from etna.core.batched.sweeps import _Ramp  # noqa: E402
from etna.core.sic_material import mobility_caughey_thomas_T  # noqa: E402

V_DATA = mdlab.FIG4B_V
DATA = mdlab.FIG4B_CCE
SIG = 0.001
PEAK = 1e18
V_SAT_N, V_SAT_P = 2.2e7, 1.0e7  # cm/s (4H-SiC; low fields here, barely matters)
N_GRID = (
    np.array([0.5, 0.75, 1.0, 1.5, 2.0, 3.0575, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0]) * 1e14
)
TAU_GRID = (0.03, 0.1, 0.3, 1.0, 3.0)
TD_GRID = (0.0, 0.3e-4, 0.6e-4)
SUB_GRID = (0.3e-4, 1.0e-4)


def chi2_free_scale(model, data=DATA, sig=SIG):
    """min over scale s of sum ((s*model - data)/sig)^2; returns (chi2, s)."""
    s = np.dot(model, data) / np.dot(model, model)
    return float(np.sum(((s * model - data) / sig) ** 2)), float(s)


def chi2_fixed(model, data=DATA, sig=SIG):
    """normalised at -50 V (the paper's 100 % point), 6 free points."""
    m = model / model[-1]
    return float(np.sum(((m[:-1] - data[:-1]) / sig) ** 2))


# --------------------------------------------------------------------- engine field data
def dark_fields(sub, Ns, voltages):
    """psi, n, |E| at cell centres, W_edge for a uniform-N_D batch at each bias."""
    pb, x = mdlab.uniform_batch(sub, Ns)
    ramp = _Ramp(pb, 1e-8, 10.0)
    out = []
    for V in voltages:
        res = ramp.to(-V) if abs(V) > 1e-12 else pb.solve(ramp.psi, tol=1e-8)
        E = np.abs(-(res.psi[:, 1:] - res.psi[:, :-1]) / np.diff(x)[None, :])
        W = depletion_width(x, res.n, pb.donors, sub, mdlab.EPI_CM)
        out.append((E, W, res.n))
    return pb, x, out, ramp.converged


@njit(cache=True)
def _hecht_induced(dx, v, W, i_start, i_end, tau, toward_low):
    """Induced charge (units of q) of a carrier starting in cell i_start and drifting to
    the depleted-region boundary with cell velocities v, weighting field 1/W, trapping
    time tau. toward_low: hole (drifts to decreasing x); else electron."""
    q = 0.0
    t = 0.0
    if toward_low:
        i = i_start
        while i >= i_end:
            dt = dx[i] / v[i]
            # average survival over the cell: exp(-t/tau) * (1 - exp(-dt/tau)) * tau / dt
            if dt / tau < 1e-6:
                surv = np.exp(-t / tau)
            else:
                surv = np.exp(-t / tau) * (1.0 - np.exp(-dt / tau)) * tau / dt
            q += dx[i] / W * surv
            t += dt
            i -= 1
    else:
        i = i_start
        while i <= i_end:
            dt = dx[i] / v[i]
            if dt / tau < 1e-6:
                surv = np.exp(-t / tau)
            else:
                surv = np.exp(-t / tau) * (1.0 - np.exp(-dt / tau)) * tau / dt
            q += dx[i] / W * surv
            t += dt
            i += 1
    return q


@njit(cache=True)
def hecht_curve(dx, vh, ve, Gc, W, i_j, i_w, taus_h, taus_e, eta_diff, G_undep):
    """Returns Qh[len(taus_h)], Qe[len(taus_e)] (pairs-weighted induced charge, units of
    q x pairs) for generation Gc per cell in [i_j, i_w] (depleted), plus the undepleted
    part G_undep injected at the edge with efficiency eta_diff (holes only)."""
    Qh = np.zeros(len(taus_h))
    Qe = np.zeros(len(taus_e))
    for k in range(len(taus_h)):
        tot = 0.0
        for i in range(i_j, i_w + 1):
            if Gc[i] > 0:
                tot += Gc[i] * _hecht_induced(dx, vh, W, i, i_j, taus_h[k], True)
        tot += eta_diff * G_undep * _hecht_induced(dx, vh, W, i_w, i_j, taus_h[k], True)
        Qh[k] = tot
    for k in range(len(taus_e)):
        tot = 0.0
        for i in range(i_j, i_w + 1):
            if Gc[i] > 0:
                tot += Gc[i] * _hecht_induced(dx, ve, W, i, i_w, taus_e[k], False)
        Qe[k] = tot
    return Qh, Qe


def main():
    t0 = time.time()
    print("Petringa 2025 Fig. 4b: hypothesis test (free-scale chi2, sigma = 1e-3)")
    print("DATA  " + "  ".join(f"{v:>7.4f}" for v in DATA))
    print()

    # ----------------------------------------------------------------- (a)(b-SRH)(c): DD
    print(
        "== (a)/(b-SRH)/(c): steady-state DD grid, N_D x p+ x t_dead x tau_scale (bragg) =="
    )
    rows = []
    for sub in SUB_GRID:
        x = mdlab.mesh(sub)
        combos = list(itertools.product(N_GRID, TD_GRID, TAU_GRID))
        Gs = []
        for N, td, tau in combos:
            G, _ = mdlab.tilted_generation(x, t_dead_cm=td, epi_end=sub + mdlab.EPI_CM)
            G[x < sub] = 0.0
            Gs.append(G / G.max() * PEAK)
        Ns = np.array([c[0] for c in combos])
        pb, _ = mdlab.uniform_batch(sub, Ns)
        cce, conv = mdlab.cce_sweep_custom(
            pb,
            sub,
            V_DATA,
            np.array(Gs),
            N_D_mob=Ns,
            tau_scale=[c[2] for c in combos],
            hole_mob_from="N_D",
        )
        for j, (N, td, tau) in enumerate(combos):
            c2, s = chi2_free_scale(cce[:, j])
            rows.append(
                dict(
                    sub=sub,
                    N=N,
                    td=td,
                    tau=tau,
                    cce=cce[:, j],
                    chi2=c2,
                    s=s,
                    chi2fix=chi2_fixed(cce[:, j]),
                    conv=bool(conv[:, j].all()),
                )
            )
    rows.sort(key=lambda r: r["chi2"])
    print(
        f"{len(rows)} DD curves in {time.time() - t0:.1f} s; best 12 by chi2 (dof 6 - 0 = 6):"
    )
    for r in rows[:12]:
        print(
            f"  p+{r['sub'] * 1e4:.1f} N{r['N']:.2e} td{r['td'] * 1e4:.1f} tau x{r['tau']:<4} "
            f"chi2 {r['chi2']:9.1f} (scale {r['s']:.4f}; fixed-norm {r['chi2fix']:9.1f}) "
            + " ".join(f"{v * r['s']:.4f}" for v in r["cce"])
            + (" " if r["conv"] else " NOCONV")
        )
    print("  preset (p+0.3, 3.06e14, td 0, tau x1):")
    for r in rows:
        if (
            r["sub"] == 0.3e-4
            and abs(r["N"] - 3.0575e14) < 1
            and r["td"] == 0
            and r["tau"] == 1
        ):
            print(
                f"    chi2 {r['chi2']:.1f} scale {r['s']:.4f} fixed {r['chi2fix']:.1f} "
                + " ".join(f"{v * r['s']:.4f}" for v in r["cce"])
            )
    # (a): where does W reach the track end?
    print(
        "\n  (a) bias at which W_edge reaches the track end (7.61 um below the surface):"
    )
    for sub in SUB_GRID:
        pb, x, fld, _ = dark_fields(sub, N_GRID, np.arange(0, -61, -1.0))
        Wk = np.array([f[1] for f in fld])  # (K, B)
        target = 7.61e-4 - sub
        for b, N in enumerate(N_GRID):
            k = np.flatnonzero(Wk[:, b] >= target)
            vreach = -float(k[0]) if len(k) else None
            print(
                f"    p+{sub * 1e4:.1f} N_D {N:.2e}: W>=track end at {vreach} V; "
                f"W(0/-10/-30/-50) = {Wk[0, b] * 1e4:.2f}/{Wk[10, b] * 1e4:.2f}/{Wk[30, b] * 1e4:.2f}/{Wk[50, b] * 1e4:.2f} um"
            )
    print()

    # ----------------------------------------------------- (b-trap) and (d): post-processing
    print(
        "== (b-trap) Hecht trapping on the engine field; (d) columnar loss a/(1+E/E_c) =="
    )
    taus = np.geomspace(1e-11, 1e-5, 37)
    a_grid = np.linspace(0.0, 0.6, 61)
    Ec_grid = np.geomspace(1e2, 1e6, 81)
    best_trap, best_col, best_trap_e, best_trap_h = [], [], [], []
    for sub in SUB_GRID:
        pb, x, fld, convP = dark_fields(sub, N_GRID, V_DATA)
        dx = np.diff(x)
        xc = 0.5 * (x[1:] + x[:-1])
        epi_end = sub + mdlab.EPI_CM
        # DD CCE (tau x1) for the diffusion efficiency of the undepleted part
        for td in TD_GRID:
            G, _ = mdlab.tilted_generation(x, t_dead_cm=td, epi_end=epi_end)
            G[x < sub] = 0.0
            Gc = 0.5 * (G[1:] + G[:-1]) * dx  # pairs per cell
            Gtot = Gc.sum()
            ccedd, _ = mdlab.cce_sweep_custom(
                pb, sub, V_DATA, G / G.max() * PEAK, N_D_mob=N_GRID, hole_mob_from="N_D"
            )
            for b, N in enumerate(N_GRID):
                mu_h = mobility_caughey_thomas_T(N, 300.0, "hole", PARAMS)
                mu_e = mobility_caughey_thomas_T(N, 300.0, "electron", PARAMS)
                Qh_all, Qe_all, Edep_all, Gdep_all, eta_all, Gund_all = (
                    [],
                    [],
                    [],
                    [],
                    [],
                    [],
                )
                for k, V in enumerate(V_DATA):
                    E, W, n = fld[k]
                    Eb, Wb = E[b], W[b]
                    i_j = int(np.searchsorted(xc, sub))
                    i_w = int(np.searchsorted(xc, sub + Wb)) - 1
                    i_w = max(i_w, i_j)
                    dep = np.zeros(len(xc), bool)
                    dep[i_j : i_w + 1] = True
                    G_dep, G_und = Gc[dep].sum(), Gc[~dep & (xc >= sub)].sum()
                    eta = (
                        1.0
                        if G_und <= 0
                        else float(np.clip((ccedd[k, b] * Gtot - G_dep) / G_und, 0, 1))
                    )
                    vh = mu_h * Eb / (1 + mu_h * Eb / V_SAT_P) + 1e-30
                    ve = mu_e * Eb / (1 + mu_e * Eb / V_SAT_N) + 1e-30
                    Qh, Qe = hecht_curve(
                        dx, vh, ve, Gc, Wb, i_j, i_w, taus, taus, eta, G_und
                    )
                    Qh_all.append(Qh)
                    Qe_all.append(Qe)
                    Edep_all.append(Eb)
                    Gdep_all.append(dep)
                    eta_all.append(eta)
                    Gund_all.append(G_und)
                Qh_all = np.array(Qh_all)
                Qe_all = np.array(Qe_all)  # (K, ntau)
                # (b-trap): CCE(V) = (Qh[tau_h] + Qe[tau_e]) / Gtot
                cce_t = (Qh_all[:, :, None] + Qe_all[:, None, :]) / Gtot  # (K, nh, ne)
                M = cce_t.reshape(len(V_DATA), -1).T  # (nh*ne, K)
                s = (M @ DATA) / np.einsum("ij,ij->i", M, M)
                c2 = (((s[:, None] * M - DATA) / SIG) ** 2).sum(1)
                ib = int(np.argmin(c2))
                ih, ie = divmod(ib, len(taus))
                best_trap.append(
                    (c2[ib], sub, N, td, taus[ih], taus[ie], s[ib], M[ib] * s[ib])
                )
                # hole-only / electron-only trapping (other carrier not trapped)
                ihh = int(np.argmin(c2.reshape(len(taus), len(taus))[:, -1]))
                best_trap_h.append(
                    (
                        c2.reshape(len(taus), len(taus))[ihh, -1],
                        sub,
                        N,
                        td,
                        taus[ihh],
                        s.reshape(len(taus), len(taus))[ihh, -1],
                    )
                )
                iee = int(np.argmin(c2.reshape(len(taus), len(taus))[-1, :]))
                best_trap_e.append(
                    (
                        c2.reshape(len(taus), len(taus))[-1, iee],
                        sub,
                        N,
                        td,
                        taus[iee],
                        s.reshape(len(taus), len(taus))[-1, iee],
                    )
                )
                # (d) columnar: loss a/(1+E/Ec) at the generation point; undepleted part:
                # loss a (E ~ 0) then diffusion efficiency eta
                cce_c = np.empty((len(V_DATA), len(a_grid), len(Ec_grid)))
                for k in range(len(V_DATA)):
                    Eb, dep, eta, G_und = (
                        Edep_all[k],
                        Gdep_all[k],
                        eta_all[k],
                        Gund_all[k],
                    )
                    g = Gc[dep]
                    e = Eb[dep]
                    f = 1.0 / (1.0 + e[None, :] / Ec_grid[:, None])  # (nEc, ncell)
                    dep_coll = (
                        g.sum() - a_grid[:, None] * (f * g[None, :]).sum(1)[None, :]
                    )  # (na, nEc)
                    und_coll = eta * G_und * (1.0 - a_grid)[:, None]
                    cce_c[k] = (dep_coll + und_coll) / Gtot
                M = cce_c.reshape(len(V_DATA), -1).T
                s = (M @ DATA) / np.einsum("ij,ij->i", M, M)
                c2 = (((s[:, None] * M - DATA) / SIG) ** 2).sum(1)
                ib = int(np.argmin(c2))
                ia, iE = divmod(ib, len(Ec_grid))
                best_col.append(
                    (c2[ib], sub, N, td, a_grid[ia], Ec_grid[iE], s[ib], M[ib] * s[ib])
                )
    best_trap.sort(key=lambda r: r[0])
    best_col.sort(key=lambda r: r[0])
    best_trap_h.sort(key=lambda r: r[0])
    best_trap_e.sort(key=lambda r: r[0])
    print(f"post-processing done at {time.time() - t0:.1f} s")
    print(
        "\n(b-trap) best per (p+, N_D, t_dead), 2 params tau_h, tau_e (dof 4); top 12:"
    )
    for c2, sub, N, td, th, te, s, m in best_trap[:12]:
        print(
            f"  p+{sub * 1e4:.1f} N{N:.2e} td{td * 1e4:.1f} tau_h {th:.1e} tau_e {te:.1e} scale {s:.4f} "
            f"chi2 {c2:8.1f}  " + " ".join(f"{v:.4f}" for v in m)
        )
    print("(b-trap) hole-only trapping (tau_e -> inf), 1 param (dof 5); top 6:")
    for c2, sub, N, td, th, s in best_trap_h[:6]:
        print(
            f"  p+{sub * 1e4:.1f} N{N:.2e} td{td * 1e4:.1f} tau_h {th:.1e} scale {s:.4f} chi2 {c2:8.1f}"
        )
    print("(b-trap) electron-only trapping (tau_h -> inf), 1 param (dof 5); top 6:")
    for c2, sub, N, td, te, s in best_trap_e[:6]:
        print(
            f"  p+{sub * 1e4:.1f} N{N:.2e} td{td * 1e4:.1f} tau_e {te:.1e} scale {s:.4f} chi2 {c2:8.1f}"
        )
    print("\n(d) columnar a/(1+E/E_c), 2 params (dof 4); top 12:")
    for c2, sub, N, td, a, Ec, s, m in best_col[:12]:
        print(
            f"  p+{sub * 1e4:.1f} N{N:.2e} td{td * 1e4:.1f} a {a:.3f} E_c {Ec:.2e} V/cm scale {s:.4f} "
            f"chi2 {c2:8.1f}  " + " ".join(f"{v:.4f}" for v in m)
        )
    print("\n(d) best per N_D (p+ 0.3, td 0):")
    for N in N_GRID:
        cand = [r for r in best_col if r[1] == 0.3e-4 and r[2] == N and r[3] == 0.0]
        c2, sub, N, td, a, Ec, s, m = cand[0]
        print(
            f"  N{N:.2e}: a {a:.3f} E_c {Ec:.2e} scale {s:.4f} chi2 {c2:8.1f}  "
            + " ".join(f"{v:.4f}" for v in m)
        )
    print("\n(b-trap) best per N_D (p+ 0.3, td 0):")
    for N in N_GRID:
        cand = [r for r in best_trap if r[1] == 0.3e-4 and r[2] == N and r[3] == 0.0]
        c2, sub, N, td, th, te, s, m = cand[0]
        print(
            f"  N{N:.2e}: tau_h {th:.1e} tau_e {te:.1e} scale {s:.4f} chi2 {c2:8.1f}  "
            + " ".join(f"{v:.4f}" for v in m)
        )
    print(f"\ntotal {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
