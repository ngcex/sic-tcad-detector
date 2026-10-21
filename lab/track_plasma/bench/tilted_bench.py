"""Bench of the tilted-track (2.5D) composition vs Petringa 2025 Fig. 4b (issue #9 follow-up).

    uv run python lab/track_plasma/bench/tilted_bench.py SECTION [--workers 8]

Sections (numbers recorded in lab/track_plasma/notes/LOG.md, sections 10-13):

  slices      the discriminator: transverse-slice loss at the fields of -20 / -50 V
              (E_perp = 1.45e4 / 4.2e4 V/cm), R / field-gradient / mesh / core / physics
              brackets
  scan        PlasmaTables (cached in notes/tables/*.json) + compositions for the parameter
              grid (tau scale, s, R, mobility, Auger, carrier scattering, N_D, T)
  relin       recompute the linear-response rows of the cached tables (LET_LINEAR)
  calibrate   tau-scale calibration on the digitised Fig. 4b (free scale), nuisance
              brackets, identifiability (Jacobian singular values)
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from etna import microdosimeter_config

NOTES = Path(__file__).resolve().parents[1] / "notes"
TABLES = NOTES / "tables"
DATA_V = np.array([0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0])
DATA = np.array([0.9083, 0.9662, 0.9806, 0.9904, 0.9932, 0.9976, 1.0000])


def sec_slices(args):
    from etna.core.track.cross_section import slice_loss
    from etna.core.track.radial import beta_of, gaussian, thermalized, whk

    rad = thermalized(whk(beta_of(0.4)), 10e-7)
    rows = []

    def go(tag, E, L=300.0, radial=rad, **kw):
        r = slice_loss(E, L, radial, **kw)
        rows.append(
            (tag, E, r.loss, r.ledger_residual, r.converged, r.n_steps, r.wall_s)
        )
        print(f"{tag:30s} E={E:8.3g}  loss={r.loss:.4e}  ledger={r.ledger_residual:.0e} "
              f"conv={r.converged} steps={r.n_steps} {r.wall_s:.1f}s", flush=True)  # fmt: skip

    for E in (1.45e4, 4.2e4):
        for R in (5e-4, 10e-4, 20e-4):
            go(f"R = {R * 1e4:.0f} um", E, R_cm=R)
        go("gradient N_space 0.82e14, R 10", E, N_space=0.82e14)
        go("gradient N_space 0.82e14, R 20", E, N_space=0.82e14, R_cm=20e-4)
        go("mesh r_min 0.25 nm, x1.25", E, r_min_cm=0.25e-7, growth=1.25)
        go("WHK s = 1 nm", E, radial=thermalized(whk(beta_of(0.4)), 1e-7))
        go("WHK s = 50 nm", E, radial=thermalized(whk(beta_of(0.4)), 50e-7))
        go("Gaussian 300 nm", E, radial=gaussian(300e-7))
        go("LET 1 (linear)", E, L=1.0)
        go("Auger off", E, C_n=0.0, C_p=0.0)
        go("Auger injection", E, auger_model="injection")
        go("carrier-carrier scattering", E, carrier_scattering=True)
        go("frozen-field v_sat", E, mu_field=True)
        go("tau x1/10", E, tau_n=1e-10, tau_p=60e-9)
        go("tau x1/30", E, tau_n=1e-9 / 30, tau_p=600e-9 / 30)


# ---------------------------------------------------------------- scan / calibrate
BASE = dict(tau_scale=1.0, s_nm=10.0, R_um=10.0, mu_field=False, auger_model="constant",
            carrier_scattering=False, T=300.0, N_space=0.0)  # fmt: skip
TAUS = (1.0, 1 / 3, 1 / 10, 1 / 20, 1 / 30, 1 / 50, 1 / 100)


def _key(p):
    p = {k: v for k, v in p.items() if not (k == "N_space" and v == 0.0)}
    return "_".join(f"{k}={p[k]:.4g}" if isinstance(p[k], float) else f"{k}={p[k]}"
                    for k in sorted(p))  # fmt: skip


def get_table(p, workers):
    from etna.core.track.grid import PARAMS
    from etna.core.sic_material import srh_lifetime
    from etna.core.track.tilted import PlasmaTable, plasma_table

    TABLES.mkdir(parents=True, exist_ok=True)
    f = TABLES / (_key(p) + ".json")
    if f.exists():
        return PlasmaTable.from_dict(json.loads(f.read_text()))
    T = p["T"]
    t0 = time.perf_counter()
    tab = plasma_table(
        junction_um=0.3,
        s_nm=p["s_nm"],
        R_um=p["R_um"],
        workers=workers,
        tau_n=srh_lifetime(T, "electron", PARAMS) * p["tau_scale"],
        tau_p=srh_lifetime(T, "hole", PARAMS) * p["tau_scale"],
        auger_model=p["auger_model"],
        carrier_scattering=p["carrier_scattering"],
        mu_field=p["mu_field"],
        N_space=p["N_space"],
        N_D=1e14,
        T=T,
    )
    print(f"  table {_key(p)}: {time.perf_counter() - t0:.0f} s, "
          f"{tab.diagnostics}", flush=True)  # fmt: skip
    f.write_text(json.dumps(tab.to_dict()))
    return tab


def compose(p, table, N_D=1e14, V=DATA_V):
    from etna.core.track.tilted import tilted_alpha_cce

    cfg = microdosimeter_config(N_D=N_D, T=p["T"])
    return tilted_alpha_cce(cfg, V, table=table, tau_scale=p["tau_scale"],
                            s_nm=p["s_nm"])  # fmt: skip


def chi2_free_scale(shape, sigma=1e-3, pts=slice(None)):
    """min over a of sum ((a m - d) / sigma)^2 on the points ``pts``; (chi2, a)."""
    m, d = shape[pts], DATA[pts]
    a = float(np.dot(m, d) / np.dot(m, m))
    return float(np.sum(((a * m - d) / sigma) ** 2)), a


def variants():
    out = []
    for t in TAUS:
        out.append(dict(BASE, tau_scale=t))
    for t in (1 / 10, 1 / 30):
        for s in (1.0, 50.0):
            out.append(dict(BASE, tau_scale=t, s_nm=s))
        for R in (5.0, 20.0):
            out.append(dict(BASE, tau_scale=t, R_um=R))
        out.append(dict(BASE, tau_scale=t, mu_field=True))
        out.append(dict(BASE, tau_scale=t, auger_model="injection"))
        out.append(dict(BASE, tau_scale=t, carrier_scattering=True))
        out.append(dict(BASE, tau_scale=t, T=295.0))
        # in-plane field gradient of the 1e14 epi at 65 deg: N_D sin^2(theta)
        out.append(dict(BASE, tau_scale=t, N_space=0.82e14, R_um=20.0))
    out.append(dict(BASE, s_nm=1.0))
    out.append(dict(BASE, s_nm=50.0))
    return out


def sec_scan(args):
    rows = []
    for p in variants():
        tab = get_table(p, args.workers)
        for N_D in (0.5e14, 1e14, 1.5e14, 2e14):
            if N_D != 1e14 and not (p["s_nm"] == 10.0 and p["R_um"] == 10.0
                                    and not p["mu_field"] and p["T"] == 300.0
                                    and p["auger_model"] == "constant"
                                    and not p["carrier_scattering"]
                                    and p["N_space"] == 0.0):  # fmt: skip
                continue
            r = compose(p, tab, N_D=N_D)
            sh = r.shape()
            c_all, a_all = chi2_free_scale(sh)
            c_hi, a_hi = chi2_free_scale(sh, pts=slice(3, None))
            rows.append(dict(p, N_D=N_D, cce=r.cce.tolist(), column=r.cce_column.tolist(),
                             loss=r.plasma_loss.tolist(), shape=sh.tolist(),
                             chi2_all=c_all, scale_all=a_all, chi2_tail=c_hi,
                             scale_tail=a_hi, col_ledger=float(r.column_ledger.max()),
                             loss_band=r.loss_by_field().tolist()))  # fmt: skip
            print(f"{_key(p)} N_D={N_D:.2g}: shape {np.round(sh, 4)} chi2(all, sigma 1e-3) "
                  f"{c_all:8.1f} chi2(-20..-50) {c_hi:7.1f} loss(50V) "
                  f"{r.plasma_loss[-1]:.2e}", flush=True)  # fmt: skip
    (NOTES / "tilted_scan.json").write_text(json.dumps(rows, indent=1))


def _is(r, **kw):
    q = dict(BASE, N_D=1e14)
    q.update(kw)
    return all(r.get(k, BASE.get(k)) == v for k, v in q.items())


def sec_calibrate(args):
    """(a) posterior of log(tau_scale) of the base model on the tail (-20..-50 V) and on
    all 7 points (free scale, sigma 5e-4 / 1e-3 / 2e-3); (b) the tail-best tau of every
    systematic variant (its shift from the base interpolated in log tau between the
    1/10 and 1/30 runs); (c) where the low-bias loss comes from (loss by local E_perp);
    (d) identifiability: singular values of the scale-projected shape Jacobian in
    (log tau, log s, N_D) and the s profile."""
    from etna.core.track.tilted import calibrate_tau_scale

    rows = json.loads((NOTES / "tilted_scan.json").read_text())
    base = sorted([r for r in rows if _is(r, tau_scale=r["tau_scale"])],
                  key=lambda r: r["tau_scale"])  # fmt: skip
    taus = np.array([r["tau_scale"] for r in base])
    shapes = np.array([r["shape"] for r in base])
    out = {"base_taus": taus.tolist()}
    print("(a) base model (s 10 nm, R 10 um, N_D 1e14), posterior of tau_scale")
    for name, pts in (("tail", slice(3, None)), ("all", slice(None))):
        for sigma in (5e-4, 1e-3, 2e-3):
            cal = calibrate_tau_scale(taus, shapes, DATA, sigma=sigma, points=pts)
            out[f"{name}_{sigma:g}"] = cal.to_dict()
            print(f"  {name:4s} sigma {sigma:.0e}: tau_scale {cal.tau_scale:.4f} (68 % "
                  f"{cal.ci68[0]:.4f}-{cal.ci68[1]:.4f}, 95 % {cal.ci95[0]:.4f}-"
                  f"{cal.ci95[1]:.4f}) chi2 {cal.chi2:.1f}/dof {cal.dof} scale "
                  f"{cal.free_scale:.4f} tau_p {cal.tau_p_ns:.1f} ns "
                  f"edge={cal.at_edge}", flush=True)  # fmt: skip
    cal = calibrate_tau_scale(taus, shapes, DATA, sigma=1e-3, points=slice(3, None))
    print("  model at the tail mode (all 7 V, free scale of the tail):",
          np.round(cal.model_at_mode, 4), " data", DATA)
    # (b) systematics
    print("(b) tail-best tau_scale per variant (sigma 1e-3)")
    lg = np.log10(taus)
    var_keys = [dict(s_nm=1.0), dict(s_nm=50.0), dict(R_um=5.0), dict(R_um=20.0),
                dict(mu_field=True), dict(auger_model="injection"),
                dict(carrier_scattering=True), dict(T=295.0),
                dict(N_space=0.82e14, R_um=20.0), dict(N_D=0.5e14), dict(N_D=1.5e14),
                dict(N_D=2e14)]  # fmt: skip
    syst = []
    for vk in var_keys:
        if "N_D" in vk:
            vr = sorted([r for r in rows if _is(r, tau_scale=r["tau_scale"], **vk)],
                        key=lambda r: r["tau_scale"])  # fmt: skip
            vt = np.array([r["tau_scale"] for r in vr])
            vs = np.array([r["shape"] for r in vr])
        else:
            vr = {r["tau_scale"]: r for r in rows if _is(r, tau_scale=r["tau_scale"], **vk)}
            d10 = d30 = None
            for t, r in vr.items():
                b = next(x for x in base if abs(x["tau_scale"] - t) < 1e-12)
                d = np.array(r["shape"]) - np.array(b["shape"])
                if abs(t - 0.1) < 1e-9:
                    d10 = d
                elif abs(t - 1 / 30) < 1e-9:
                    d30 = d
            if d10 is None or d30 is None:
                print("  missing", vk)
                continue
            x = np.clip((lg - (-1.0)) / (np.log10(1 / 30) + 1.0), 0.0, 1.0)
            vs = shapes + (1 - x)[:, None] * d10[None, :] + x[:, None] * d30[None, :]
            vt = taus
        c = calibrate_tau_scale(vt, vs, DATA, sigma=1e-3, points=slice(3, None))
        syst.append(dict(variant=vk, tau_scale=c.tau_scale, chi2=c.chi2,
                         ci68=list(c.ci68)))  # fmt: skip
        print(f"  {str(vk):44s} tau_scale {c.tau_scale:.4f} (68 % {c.ci68[0]:.4f}-"
              f"{c.ci68[1]:.4f}) chi2 {c.chi2:.1f}", flush=True)  # fmt: skip
    out["systematics"] = syst
    # (c) loss by field band
    print("(c) plasma loss by local E_perp band (<2e3 | 2e3-1e4 | >1e4 V/cm) and column")
    for r in base:
        if min(abs(r["tau_scale"] - t) for t in (1.0, 0.1, 0.05, 1 / 30)) < 1e-9:
            col = np.array(r["column"])
            print(f"  tau x{r['tau_scale']:.4g}: column/column(50) "
                  f"{np.round(col / col[-1], 4)}")
            for V, band in zip(DATA_V, r["loss_band"]):
                print(f"     {V:4.0f} V: " + " | ".join(f"{b:.2e}" for b in band))
    # (d) identifiability at tau 1/30 (and 1/10): log tau, log s, N_D
    print("(d) identifiability: scale-projected Jacobian singular values (sigma 1e-3)")
    for t0, (tlo, thi) in ((1 / 30, (0.02, 0.05)), (0.1, (0.05, 1 / 3))):
        def shp(**kw):
            r = next(r for r in rows if _is(r, **kw))
            return np.array(r["shape"])

        m = shp(tau_scale=t0)
        J = np.stack([
            (shp(tau_scale=thi) - shp(tau_scale=tlo)) / np.log10(thi / tlo),
            (shp(tau_scale=t0, s_nm=50.0) - shp(tau_scale=t0, s_nm=1.0)) / np.log10(50.0),
            (shp(tau_scale=t0, N_D=1.5e14) - shp(tau_scale=t0, N_D=0.5e14)) / 1.0,
        ], axis=1) / 1e-3  # fmt: skip
        for name, pts in (("tail", slice(3, None)), ("all", slice(None))):
            Jp, mp = J[pts], m[pts]
            P = np.eye(len(mp)) - np.outer(mp, mp) / np.dot(mp, mp)
            sv = np.linalg.svd(P @ Jp, compute_uv=False)
            col = np.linalg.norm(P @ Jp, axis=0)
            print(f"  tau {t0:.4g} {name:4s}: singular values {np.array2string(sv, precision=2)}"
                  f"; column norms (log tau, log10 s, N_D/1e14) "
                  f"{np.array2string(col, precision=2)}")
        for s_nm in (1.0, 10.0, 50.0):
            sh = shp(tau_scale=t0, s_nm=s_nm)
            print(f"     s {s_nm:4.0f} nm: chi2 tail {chi2_free_scale(sh, pts=slice(3, None))[0]:7.1f}"
                  f"  all {chi2_free_scale(sh)[0]:9.1f}")
    (NOTES / "tilted_calibration.json").write_text(json.dumps(out, indent=1, default=float))


def sec_relin(args):
    """Recompute the linear-response row f_lin of every cached table at the current
    ``LET_LINEAR`` (the first scan used LET 1 keV/um, which is not linear at 1e4 V/cm:
    the column self-field matters; tests/test_track_tilted.py)."""
    from etna.core.track.tilted import LET_LINEAR, PlasmaTable, _slice_job
    from etna.core.track.radial import beta_of

    for f in sorted(TABLES.glob("*.json")):
        tab = PlasmaTable.from_dict(json.loads(f.read_text()))
        st = tab.settings
        if st.get("LET_linear") == LET_LINEAR:
            continue
        keep = ("tau_n", "tau_p", "auger_model", "carrier_scattering", "mu_field", "N_D",
                "N_space", "T", "r_min_cm", "growth")  # fmt: skip
        kw = dict(R_cm=st["R_um"] * 1e-4, **{k: st[k] for k in keep if k in st})
        b = beta_of(max(float(tab.E_per_u_MeV[len(tab.s_cm) // 2]), 1e-3))
        old = tab.f_lin.copy()
        tab.f_lin = np.array([_slice_job((float(E), LET_LINEAR, b, st["s_nm"], kw))[0]
                              for E in tab.E_grid])  # fmt: skip
        st["LET_linear"] = LET_LINEAR
        f.write_text(json.dumps(tab.to_dict()))
        print(f"{f.name}: f_lin LET 1 {np.array2string(old, precision=2)}\n"
              f"   -> LET {LET_LINEAR} {np.array2string(tab.f_lin, precision=2)}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("section")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    t0 = time.perf_counter()
    {"slices": sec_slices, "scan": sec_scan, "calibrate": sec_calibrate,
     "relin": sec_relin}[args.section](
        args
    )
    print(f"[{args.section}] {time.perf_counter() - t0:.1f} s")


if __name__ == "__main__":
    main()
