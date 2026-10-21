"""c-axis low-field mobility default (issue #6 follow-up): before/after numbers of
lab/mobility_caxis/notes/LOG.md and docs/USAGE.md section 18.

    uv run python lab/mobility_caxis/bench/caxis_bench.py            # 1D facades, ~1 min
    uv run python lab/mobility_caxis/bench/caxis_bench.py --slow     # + radiation damage,
                                                                     #   FLASH transient

"before" = DeviceConfig(mobility="legacy_perp") (the default until 2026-10-21, the
perpendicular-to-c set), "after" = DeviceConfig() ("c_axis"), "ishikawa" = the opt-in
full-thesis set (holes: Hall mobility parallel to c). Every number is a facade output.
"""

import sys
import time

import numpy as np

from etna import DeviceConfig
from etna.api.models import run_transit_time
from etna.api.simulation import (
    run_cce,
    run_cce_tilted_alpha,
    run_cv,
    run_dark_current,
    run_flash_recombination,
    run_radiation_damage,
    run_temperature_sweep,
    run_transient,
)
from etna.core.physics_models import LOW_FIELD_MOBILITY
from etna.core.sic_material import SiC4H_Parameters, mobility_caughey_thomas_T

SETS = ("legacy_perp", "c_axis", "ishikawa2025")
SHORT = {"legacy_perp": "before", "c_axis": "after", "ishikawa2025": "ishikawa"}


def cfg(name, **kw):
    return DeviceConfig(mobility=name, **kw)


def row(label, vals, fmt="{:.6g}"):
    print(f"  {label:28s}" + "  ".join(fmt.format(v) for v in vals))


def rel(a, b):
    return 100.0 * (np.asarray(b) / np.asarray(a) - 1.0)


def section_mobility():
    print(
        "\n(0) low-field mobility (cm2/Vs): N = epi bulk 3.38e13 / junction 8.19e15 /"
        " p+ and n+ 1e19; T sweep at the bulk"
    )
    N = (3.382506e13, 8.18588e15, 1e19)
    for c in ("electron", "hole"):
        for s in SETS:
            p = SiC4H_Parameters(mobility=s)
            v300 = [mobility_caughey_thomas_T(n, 300.0, c, p) for n in N]
            vT = [mobility_caughey_thomas_T(N[0], T, c, p) for T in (250, 350, 450)]
            row(f"{c[0]} {SHORT[s]:8s} 300K / T250,350,450", v300 + vT, "{:8.1f}")


def section_cv():
    print("\n(1) C-V (electrostatics: must be identical)")
    out = {}
    for s in SETS[:2]:
        r = run_cv(cfg(s), 0.0, -50.0, 6)
        out[s] = np.asarray(r.y)
    d = np.max(np.abs(out["c_axis"] / out["legacy_perp"] - 1))
    print(f"  max |C_after / C_before - 1| = {d:.1e}")


def section_cce():
    print("\n(2) alpha CCE vs bias (run_cce, 18 um alpha, -2 .. -100 V)")
    V = np.array([-2.0, -5.0, -10.0, -20.0, -40.0, -70.0, -100.0])
    res = {}
    for s in SETS:
        r = run_cce(cfg(s), V[0], V[-1], 99)
        x, y = np.asarray(r.x), np.asarray(r.y)
        res[s] = np.interp(-V, -x, y)
        print(f"    [{SHORT[s]}] solver {r.metadata.get('solver')}")
    print("  V                           " + "  ".join(f"{v:9.0f}" for v in V))
    for s in SETS:
        row(f"CCE {SHORT[s]}", res[s], "{:9.6f}")
    row("after - before", res["c_axis"] - res["legacy_perp"], "{:+9.1e}")


def section_decomposition():
    """Which change moves the CCE: bench-only hybrid sets registered at run time."""
    from dataclasses import replace

    L = LOW_FIELD_MOBILITY
    leg, cax = L["legacy_perp"], L["c_axis"]
    hyb = {
        "e_only": replace(cax, name="e_only", mu_max_p=leg.mu_max_p,
                          mu_min_p=leg.mu_min_p),
        "h_only": replace(leg, name="h_only", mu_max_p=cax.mu_max_p,
                          mu_min_p=cax.mu_min_p),
        "e_mumin40": replace(cax, name="e_mumin40", mu_min_n=40.0),
    }  # fmt: skip
    L.update(hyb)
    print(
        "\n(2b) CCE decomposition (bench-only hybrids): e_only = c-axis electrons + "
        "legacy holes; h_only = legacy electrons + c-axis holes; e_mumin40 = c_axis "
        "with the legacy electron mu_min 40 (p+ minority electrons)"
    )
    V = np.array([-2.0, -10.0, -40.0, -100.0])
    base = None
    try:
        for s in ("legacy_perp", "e_only", "h_only", "c_axis", "e_mumin40"):
            r = run_cce(cfg(s), V[0], V[-1], 99)
            y = np.interp(-V, -np.asarray(r.x), np.asarray(r.y))
            base = y if base is None else base
            row(f"{s} minus before", y - base, "{:+9.2e}")
    finally:
        for k in hyb:
            L.pop(k)


def section_transit():
    print("\n(3) full-width drift transit time across the depleted epi (low field)")
    out = {}
    for s in SETS:
        r = run_transit_time(cfg(s), v_start=-10.0, v_stop=-100.0, n_points=4)
        out[s] = (np.asarray(r.metadata["t_e"]), np.asarray(r.metadata["t_h"]))
    print("  V                           " + "  ".join(f"{v:9.0f}" for v in r.x))
    for k, lab in ((0, "t_e ps"), (1, "t_h ps")):
        for s in SETS:
            row(f"{lab} {SHORT[s]}", out[s][k] * 1e12, "{:9.1f}")
        row(f"{lab} after/before", out["c_axis"][k] / out["legacy_perp"][k], "{:9.4f}")


def section_dark():
    print("\n(4) dark current (run_dark_current, 4 mm2, A)")
    out = {}
    for s in SETS[:2]:
        r = run_dark_current(cfg(s), 0.0, -100.0, 11)
        out[s] = (np.asarray(r.x), np.asarray(r.y))
    x = out["c_axis"][0]
    print("  V                           " + "  ".join(f"{v:9.0f}" for v in x[::2]))
    for s in SETS[:2]:
        row(f"I {SHORT[s]}", out[s][1][::2], "{:9.3e}")
    row(
        "after/before - 1",
        rel(out["legacy_perp"][1], out["c_axis"][1])[::2] / 100,
        "{:+9.1e}",
    )


def section_temperature():
    print("\n(5) CCE vs T at -30 V (run_temperature_sweep, dd)")
    for method in ("dd", "hecht"):
        out = {}
        for s in SETS[:2]:
            r = run_temperature_sweep(cfg(s), [250.0, 300.0, 350.0, 450.0], -30.0,
                                      method=method)  # fmt: skip
            out[s] = np.asarray(r.y)
        row(f"{method} before (250/300/350/450)", out["legacy_perp"], "{:9.6f}")
        row(f"{method} after", out["c_axis"], "{:9.6f}")


def section_tilted():
    print("\n(6) tilted-alpha CCE (Petringa Fig. 4b comparison, run_cce_tilted_alpha)")
    import csv

    from etna._paths import data_path

    rows = list(
        csv.reader(open(data_path("validation/micro_fig4b_tilted_alpha_cce.csv")))
    )
    print("   data file header:", rows[0][:4])
    out = {}
    for s in SETS[:2]:
        r = run_cce_tilted_alpha(cfg(s), 0.0, -60.0, 61)
        out[s] = (np.asarray(r.x), np.asarray(r.y))
    x = out["c_axis"][0]
    pick = [0, 1, 2, 5, 10, 20, 40, 60]
    print("  V                           " + "  ".join(f"{x[i]:8.0f}" for i in pick))
    for s in SETS[:2]:
        row(f"CCE {SHORT[s]}", out[s][1][pick], "{:8.5f}")
    return out


def section_validation():
    print("\n(6b) committed validation comparisons (etna.viz.validation)")
    import etna
    from etna.viz.validation import tilted_alpha_shape, w_creep

    for s in SETS[:2]:
        m = tilted_alpha_shape(etna.microdosimeter_config(mobility=s)).metadata
        print(f"  Petringa Micro Fig. 4b [{SHORT[s]}]: rms {m['rms']:.4f}, max|res| "
              f"(>=20 V) {m['max_abs_residual_ge20V']:.4f}, residuals "
              + " ".join(f"{v:+.4f}" for v in m["residual"]))  # fmt: skip
    for s in SETS[:2]:
        m = w_creep(cfg(s)).metadata
        keys = [k for k in m if "resid" in k or "rms" in k]
        print(f"  Photons Fig. 4b W labels [{SHORT[s]}]: "
              + ", ".join(f"{k}={np.round(m[k], 4)}" for k in keys))  # fmt: skip


def section_track():
    print(
        "\n(6c) single-ion track engine (opt-in; MESH_FAST, gaussian 50 nm, 50 V): "
        "CCE, recombined fraction, t95 of the cathode charge (ps)"
    )
    from etna.core.track.phd import MESH_FAST, simulate_track
    from etna.core.track.radial import gaussian

    for let in (10.0, 1000.0):
        for s in SETS[:2]:
            r = simulate_track(cfg(s), LET_keV_um=let, V_bias=50.0,
                               radial=gaussian(50e-7), mesh=MESH_FAST)  # fmt: skip
            t, i = np.asarray(r.times), np.asarray(r.current_cathode)
            q = np.concatenate([[0.0], np.cumsum(0.5 * (i[1:] + i[:-1]) * np.diff(t))])
            t95 = t[np.searchsorted(q / q[-1], 0.95)]
            print(f"  LET {let:6.0f} [{SHORT[s]:6s}] CCE {r.cce:.6f}  recombined "
                  f"{r.recombined / r.N_gen:.3e}  t95 {t95 * 1e12:8.1f} ps  "
                  f"ledger {r.ledger_residual:.1e}")  # fmt: skip


def section_flash():
    print("\n(7) FLASH steady-state recombination CCE (run_flash_recombination, -30 V)")
    out = {}
    for s in SETS[:2]:
        r = run_flash_recombination(cfg(s))
        out[s] = (np.asarray(r.x), np.asarray(r.y))
    x = out["c_axis"][0]
    print("  dose rate Gy/s              " + "  ".join(f"{v:9.1e}" for v in x))
    for s in SETS[:2]:
        row(f"CCE {SHORT[s]}", out[s][1], "{:9.6f}")


def section_slow():
    print("\n(8) radiation damage CCE vs proton fluence (run_radiation_damage, -40 V)")
    out = {}
    for s in SETS[:2]:
        t = time.perf_counter()
        r = run_radiation_damage(cfg(s))
        out[s] = (np.asarray(r.x), np.asarray(r.y))
        print(f"    [{SHORT[s]}] {time.perf_counter() - t:.1f} s")
    x = out["c_axis"][0]
    print("  fluence                     " + "  ".join(f"{v:9.1e}" for v in x))
    for s in SETS[:2]:
        row(f"CCE {SHORT[s]}", out[s][1], "{:9.5f}")
    print("\n(9) FLASH transient CCE (run_transient, -30 V)")
    for s in SETS[:2]:
        t = time.perf_counter()
        r = run_transient(cfg(s))
        print(f"    [{SHORT[s]}] {time.perf_counter() - t:.1f} s")
        row(f"CCE {SHORT[s]}", np.asarray(r.y), "{:9.6f}")


if __name__ == "__main__":
    t0 = time.perf_counter()
    section_mobility()
    section_cv()
    section_cce()
    section_transit()
    section_dark()
    section_temperature()
    section_tilted()
    section_validation()
    section_flash()
    if "--slow" in sys.argv:
        section_slow()
    print(f"\ntotal {time.perf_counter() - t0:.1f} s")
