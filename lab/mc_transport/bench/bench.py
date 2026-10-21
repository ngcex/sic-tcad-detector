"""Benchmark of the condensed-history sampler (etna.core.mc_transport), issue #7.

    uv run python lab/mc_transport/bench/bench.py            # ~20 s, prints all tables
    uv run python lab/mc_transport/bench/bench.py --quick    # fewer events

References are NIST data (consistency: SiC is Bragg-composed from the same tables),
the PDG Landau formulas, a numerical Vavilov density (cumulant-checked) and the Bohr
variance. No Geant4 run (deferred by the maintainer). Output copied to
lab/mc_transport/notes/bench.out.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

from etna.core import mc_response as mr
from etna.core import mc_transport as mt
from etna.core import straggling as S
from etna.core.stopping_power import bohr_range_straggling_cm, csda_range_cm

QUICK = "--quick" in sys.argv
N = 200_000 if QUICK else 1_000_000


def mpv_fwhm(x, bins=600):
    lo, hi = np.percentile(x, [0.2, 99.0])
    h, e = np.histogram(x, bins=bins, range=(lo, hi))
    c = 0.5 * (e[1:] + e[:-1])
    k = np.ones(7) / 7
    hs = np.convolve(h.astype(float), k, mode="same")
    i = int(np.argmax(hs))
    half = hs[i] / 2
    left = np.nonzero(hs[:i] < half)[0]
    right = np.nonzero(hs[i:] < half)[0]
    w = c[i + right[0]] - c[left[-1]] if left.size and right.size else np.nan
    return c[i], w


def ref_shape(part, E, t, mean):
    """Most probable value / FWHM of the analytic reference distribution."""
    kap = mt.vavilov_kappa(part, E, t)
    b2, _, _ = mt.kinematics(part, E)
    xi = mt.landau_xi_keV(part, E, t)
    if kap >= 10:
        sd = np.sqrt(mt.bohr_variance_keV2(part, E, t))
        return "Bohr-Gauss", mean, 2.3548 * sd, None
    if kap >= 0.05:
        tmax = xi / kap
        d = np.linspace(
            max(0, mean - 8 * np.sqrt(xi * tmax)), mean + 12 * np.sqrt(xi * tmax), 4001
        )
        p = S.vavilov_pdf_delta(d, mean, xi, kap, b2)
        i = int(np.argmax(p))
        half = p[i] / 2
        lft = d[:i][p[:i] < half][-1]
        rgt = d[i:][p[i:] < half][0]
        return "Vavilov", d[i], rgt - lft, (d, p)
    return (
        "Landau(PDG)",
        mt.landau_mpv_keV(part, E, t),
        mt.landau_fwhm_keV(part, E, t),
        None,
    )


def section(title):
    print("\n" + "=" * 100 + f"\n{title}\n" + "=" * 100)


def main():
    t_start = time.time()
    mt.sample_track_events("proton", 62.0, n_events=10)  # JIT warm-up (cached)

    section(
        "1. Mean deposit vs NIST CSDA integral (delta rays off; consistency), 10 um SiC, normal"
    )
    print(
        f"{'particle':8s} {'E MeV':>8s} {'kappa':>8s} {'mean MC keV':>12s} {'CSDA keV':>10s} "
        f"{'ratio':>7s} {'var MC':>9s} {'var Vav/Bohr':>12s} {'ratio':>6s}"
    )
    cases = [
        ("proton", 1.0),
        ("proton", 10.0),
        ("proton", 62.0),
        ("proton", 150.0),
        ("proton", 250.0),
        ("alpha", 5.486),
        ("c12", 12 * 62.0),
        ("c12", 12 * 290.0),
    ]
    for part, E in cases:
        ev = mt.sample_track_events(part, E, n_events=N, delta_rays=False, seed=11)
        m = ev.edep_keV.mean()
        ref = float(mt.mean_csda_deposit_keV(part, E, 10.0))
        v = ev.edep_keV.var()
        vb = mt.bohr_variance_keV2(part, E, 10.0)
        print(
            f"{part:8s} {E:8.2f} {mt.vavilov_kappa(part, E, 10):8.3g} {m:12.3f} {ref:10.3f} "
            f"{m / ref:7.4f} {v:9.1f} {vb:12.1f} {v / vb:6.3f}"
        )
    print(
        "(1 MeV p and 5.486 MeV alpha lose a large fraction of E in 10 um: the fixed-E "
        "Bohr/Vavilov variance is not the reference there, the slowing adds variance)"
    )

    section(
        "2. Distribution shape vs analytic theory (delta rays off = full loss local)"
    )
    print(
        f"{'case':28s} {'kappa':>7s} {'ref':>12s} {'MPV MC':>8s} {'MPV ref':>8s} "
        f"{'FWHM MC':>8s} {'FWHM ref':>8s}"
    )
    for part, E, t in [
        ("proton", 250.0, 10),
        ("proton", 150.0, 10),
        ("proton", 62.0, 10),
        ("proton", 62.0, 2),
        ("c12", 12 * 290.0, 10),
        ("proton", 20.0, 10),
        ("proton", 10.0, 10),
        ("c12", 12 * 62.0, 10),
        ("alpha", 5.486, 2),
    ]:
        ev = mt.sample_track_events(
            part, E, sv_um=(100, 100, t), n_events=N, delta_rays=False, seed=12
        )
        mpv, fw = mpv_fwhm(ev.edep_keV)
        name, mref, fref, _ = ref_shape(part, E, t, ev.edep_keV.mean())
        print(
            f"{part} {E:7.1f} MeV t={t:>3} um      {mt.vavilov_kappa(part, E, t):7.3g} "
            f"{name:>12s} {mpv:8.2f} {mref:8.2f} {fw:8.2f} {fref:8.2f}"
        )
    print(
        "Landau rows: the PDG notes thin-Si distributions are 'significantly wider than the "
        "Landau width' (Bichsel 1988); the Urban model (Geant4) reproduces that broadening."
    )

    section(
        "3. Alpha 5.486 MeV range in SiC (thick slab, normal): MC vs CSDA and Bohr straggling"
    )
    ev = mt.sample_track_events(
        "alpha", 5.486, sv_um=(100, 100, 40), n_events=N, seed=13
    )
    R = float(csda_range_cm("alpha", "sic", 5.486)) * 1e4
    sb = bohr_range_straggling_cm("alpha", 5.486) * 1e4
    s = ev.stop_path_um
    print(
        f"mean stop depth {np.nanmean(s):.3f} um  (CSDA {R:.3f} um)  ratio {np.nanmean(s) / R:.4f}"
    )
    print(f"range straggling MC {np.nanstd(s):.3f} um  (Bohr integral {sb:.3f} um)")
    print(
        f"deposit {ev.edep_keV.mean():.1f} keV of 5486 (ionisation; nuclear loss "
        f"not deposited: {5486 - ev.edep_keV.mean():.1f} keV)"
    )

    section(
        "4. Adaptive steps: step_fraction convergence (alpha 5.486 MeV, 20 um SiC, stops inside)"
    )
    for sf in (0.4, 0.2, 0.1, 0.05, 0.02, 0.01):
        ev = mt.sample_track_events(
            "alpha",
            5.486,
            sv_um=(100, 100, 20),
            n_events=N // 5,
            step_fraction=sf,
            delta_rays=False,
            fluctuations="none",
        )
        print(
            f"step_fraction {sf:5.2f}: mean stop {np.nanmean(ev.stop_path_um):.4f} um, "
            f"deposit {ev.edep_keV.mean():.3f} keV"
        )

    section("5. Delta-ray escape and equilibrium (62 MeV p, 100x100x10 um SiC, normal)")
    ref = float(mt.mean_csda_deposit_keV("proton", 62.0, 10.0))
    for sur, rf in ((0.0, 1.0), (0.0, 0.7), (20.0, 1.0), (100.0, 1.0), (300.0, 1.0)):
        ev = mt.sample_track_events(
            "proton", 62.0, n_events=N, seed=14, surround_um=sur, delta_range_factor=rf
        )
        hit = ev.chord_um > 0
        per = ev.edep_keV.sum() / hit.sum()
        print(
            f"surround {sur:5.0f} um, range factor {rf}: deposit per crossing {per:8.3f} keV "
            f"= {per / ref:6.4f} x CSDA mean; delta part {ev.edep_delta_keV.sum() / hit.sum():.3f} keV"
        )

    section("6. T_cut plateau (62 MeV p, bare 100x100x10 um SiC, normal, l = t)")
    for tc in (1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 139.0):
        a = mr.adaptive_lineal_energy(
            "proton",
            62.0,
            rel_tol=0.003,
            batch=200_000,
            max_events=N,
            tcut_keV=tc,
            seed=15,
            keep_events=False,
        )
        print(
            f"T_cut {tc:6.1f} keV: y_F {a['y_F']:.4f} +- {a['se_y_F']:.4f}, "
            f"y_D {a['y_D']:.4f} +- {a['se_y_D']:.4f} keV/um (n={a['n_events']})"
        )
    print(
        "(T_cut >= T_max = 139 keV: no explicit delta rays, all loss local: no escape)"
    )

    section(
        "7. Lineal energy: deterministic LET x chord baseline (ETNA before) vs sampled"
    )
    print(
        "columns: base = deterministic LET x chord; local = straggling only (all loss "
        "local, no delta escape);\n bare = straggling + delta rays, SV in vacuum; "
        "embedded = same with 300 um SiC around the SV (alpha: 5 um),\n crossing events only (delta equilibrium)"
    )
    print(
        f"{'case':36s} {'':4s} {'base':>8s} {'local':>8s} {'bare':>8s} {'embedded':>8s}"
    )
    for part, E, irr in [
        ("proton", 62.0, "normal"),
        ("proton", 62.0, "isotropic"),
        ("proton", 150.0, "normal"),
        ("c12", 12 * 290.0, "normal"),
        ("c12", 12 * 62.0, "isotropic"),
        ("alpha", 5.486, "normal"),
    ]:
        lb = mr.default_l_bar_um((100, 100, 10), irr)
        kw = dict(n_events=N // 2, irradiation=irr, seed=16)
        sur = 5.0 if part == "alpha" else 300.0  # the alpha would stop in 300 um
        runs = [
            mt.sample_track_events(part, E, fluctuations="none", delta_rays=False, **kw),
            mt.sample_track_events(part, E, delta_rays=False, **kw),
            mt.sample_track_events(part, E, **kw),
            mt.sample_track_events(part, E, surround_um=sur, **kw),
        ]
        # embedded: SV-crossing events only (delta-only near-miss events, real but
        # removed by any threshold, would dominate the count in a broad beam)
        emb = runs[3]
        runs[3] = emb.edep_keV[emb.chord_um > 0]
        mom = [mr.lineal_moments(r if isinstance(r, np.ndarray) else r.edep_keV, lb)
               for r in runs]
        tag = f"{part} {E:7.1f} MeV {irr:10s} (l={lb:5.2f})"
        print(f"{tag:36s} y_F  " + " ".join(f"{m['y_F']:8.3f}" for m in mom))
        print(f"{'':36s} y_D  " + " ".join(f"{m['y_D']:8.3f}" for m in mom))

    section(
        "8. Throughput (numba float64, parallel; events/s, 62 MeV p, delta rays on)"
    )
    print(
        f"threads: NUMBA_NUM_THREADS={os.environ.get('NUMBA_NUM_THREADS', 'default')}"
    )
    import numba

    print(f"numba threads in use: {numba.get_num_threads()}")
    for n in (10_000, 100_000, 1_000_000, 4_000_000):
        t0 = time.perf_counter()
        mt.sample_track_events("proton", 62.0, n_events=n, seed=17)
        dt = time.perf_counter() - t0
        print(f"n={n:>9d}: {dt * 1e3:8.1f} ms  -> {n / dt / 1e6:6.2f} M events/s")
    numba.set_num_threads(1)
    t0 = time.perf_counter()
    mt.sample_track_events("proton", 62.0, n_events=1_000_000, seed=17)
    dt = time.perf_counter() - t0
    print(f"1 thread, n=1e6: {dt * 1e3:.1f} ms -> {1e6 / dt / 1e6:.2f} M events/s")

    section(
        "9. Nuclear interactions (not modelled): probability per 10 um SiC traversal, C-12"
    )
    n_at = 3.21 / (28.0855 + 12.011) * 6.02214076e23  # SiC molecules / cm^3
    for E, sSi, sC in ((100, 1376, 950), (290, 1273, 857)):
        p = n_at * (sSi + sC) * 1e-27 * 10e-4
        print(
            f"{E} MeV/u: sigma_R(C+Si) {sSi} mb, (C+C) {sC} mb (Tripathi formula, computed, "
            f"not measured) -> P = {p:.2e}"
        )
    print(f"\ntotal {time.time() - t_start:.1f} s")


if __name__ == "__main__":
    main()
