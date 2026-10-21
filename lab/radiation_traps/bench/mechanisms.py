"""Which defect physics carries the measured C-V donor removal? Candidate mechanisms in
the kinetic trap model vs every published irradiated-C-V number we could read.

uv run python lab/radiation_traps/bench/mechanisms.py   (numbers -> ../notes/MECHANISMS.md)

All numbers are in the si-neq NIEL convention (pinned below: the default is sic-damage
since 2026-10-21; the per-convention refit is bench/convention_refit.py).

Data (only values read in the papers; VEC = read from the PDF's vector figure paths,
DER = derived here):
- Radmanovac et al., arXiv:2510.11304 (252.7 MeV p, 1 MHz), Table 2 g per p/cm^2:
  onsemi 50 um 4e13 (Table 1): 5.60 +- 0.81 (W 10-40 um, fluences 0..2.8e12);
  CNM W2 50 um (Table 1 "(2-10)e13"; N0 7.8e13 VEC Fig. 5 at 39 um): 4.60 +- 0.33
  (W 30-40 um, 0..6.9e12); CNM W4 100 um 4e14 (N0 3.86e14 VEC Fig. 5): 4.48 +- 0.26
  (W 10-40 um, 0..2.1e13, the two highest excluded by the authors). The +- is the
  spread over depletion widths. onsemi: C-V diode-like at 2.8e12, depletion-width jump
  (compensated) at 6.9e12 (VEC Fig. 4b; Table 1 prints 6.9e13 for PIN9, the figures
  6.9e12). W4: not flat at 3.5e13 (Fig. 4a). Depth profiles N_eff(W) (VEC Fig. 4).
- Kroll et al., arXiv:2607.11663 (onsemi 50 um, N0 6.2e13 from C, 24 GeV/c p, 20 kHz):
  diode-like at 1e13, bias-independent from ~5e13 p/cm^2; authors' kappa 0.62.
- Li et al., arXiv:2503.09016 (+ Yang et al., Front. Phys. 10 (2022) 718071 for N0):
  1.5 mm NJU diodes N0 2.7e13, 100 um, 10 kHz, ~80 MeV p, fluences in neq only:
  diode-like at 2e11, nearly flat at 3.5e12 neq.
- Gsponer et al., arXiv:2310.02047: CNM 50 um 1.5e14, 10 kHz, flat at 5e14 neq (lowest).
"""

from __future__ import annotations

import dataclasses
import logging
import time

import numpy as np

from etna.api.device import DeviceConfig
from etna.core.radiation_traps import DEFAULT_TRAP_MODEL, DeepTrap
from etna.core.trap_calibration import (
    RemovalObservation,
    calibrate_removal,
    neff_window,
    parametrized_model,
    removal_observable,
    semi_insulating_fluence,
)
from etna.core.trap_kinetics import KineticTrapSolver, simulate_cv
from etna.core.trap_presets import deep_acceptor_trap, literature_kinetic

logging.disable(logging.WARNING)
T0 = time.time()
RAD = "Radmanovac et al. arXiv:2510.11304 Table 2"
OBS = {
    "onsemi": RemovalObservation(
        "onsemi g", "slope", 5.60, 0.81, 4e13, 50.0, 252.7,
        (0.0, 1.4e11, 6.9e11, 1.4e12, 2.8e12), (10.0, 40.0), 150.0, 5.0, 1e6,
        source=RAD, niel_convention="si-neq",
    ),
    "W2": RemovalObservation(
        "CNM W2 g", "slope", 4.60, 0.33, 7.8e13, 50.0, 252.7,
        (0.0, 1.4e11, 6.9e11, 1.4e12, 2.8e12, 6.9e12), (30.0, 40.0), 200.0, 5.0, 1e6,
        source=RAD, niel_convention="si-neq",
    ),
    "W4": RemovalObservation(
        "CNM W4 g", "slope", 4.48, 0.26, 3.86e14, 100.0, 252.7,
        (0.0, 6.9e11, 2.8e12, 6.9e12, 1.4e13, 2.1e13), (10.0, 40.0), 650.0, 10.0, 1e6,
        source=RAD, niel_convention="si-neq",
    ),
}  # fmt: skip
# flat-C-V brackets (lo = last diode-like, hi = first flat; fluence units of the paper)
FLAT = {
    "onsemi 1 MHz (p 252.7 MeV)": dict(N_D=4e13, E=252.7, f=1e6, lo=2.8e12, hi=6.9e12),
    "W4 1 MHz (p 252.7 MeV)": dict(N_D=3.86e14, E=252.7, f=1e6, lo=3.5e13, hi=np.inf),
    "Kroll 20 kHz (neq, kappa 0.62)": dict(N_D=6.2e13, E=None, f=2e4, lo=0.62e13,
                                           hi=0.62 * 5e13),
    "Li 1.5 mm 10 kHz (neq)": dict(N_D=2.7e13, E=None, f=1e4, lo=2e11, hi=3.5e12),
    "Gsponer 10 kHz (neq)": dict(N_D=1.5e14, E=None, f=1e4, lo=0.0, hi=5e14),
}  # fmt: skip
# onsemi Fig. 4b (VEC): N_eff at 13 / 38 um, pristine PIN6 and PIN4 (1.4e12), PIN5 (2.8e12)
DEPTH = {"PIN4": (1.4e12, (4.23e13 - 3.39e13) / (4.70e13 - 4.03e13)),
         "PIN5": (2.8e12, (4.23e13 - 1.88e13) / (4.70e13 - 3.33e13))}  # fmt: skip


def kappa(E):
    return OBS["onsemi"].kappa if E == 252.7 else 1.0


def flat_table(model):
    out = {}
    for k, d in FLAT.items():
        out[k] = semi_insulating_fluence(
            model, d["N_D"], freq_Hz=d["f"], kappa=kappa(d["E"])
        )
    return out


def flat_verdict(k, F):
    d = FLAT[k]
    return "ok" if d["lo"] < F <= d["hi"] else ("EARLY" if F <= d["lo"] else "LATE")


def depth_ratio(model, rate=1.0):
    cfg = DeviceConfig(epi_thickness_um=50.0, doping_profile="uniform", N_D=4e13)
    F = np.array([0.0, 1.4e12, 2.8e12])
    V = -np.arange(0.0, 81.0, 2.0)
    s = KineticTrapSolver.from_config(cfg, F, energy_MeV=252.7, trap_model=model)
    cv = simulate_cv(s, V, 1e6, rate_V_per_s=rate, dV_max=2.0)
    a = neff_window(V, cv.C_ac, (11.0, 15.0))
    b = neff_window(V, cv.C_ac, (36.0, 40.0))
    return (a[0] - a[1:]) / (b[0] - b[1:])


def with_trap(model, name, **kw):
    tr = tuple(
        dataclasses.replace(t, **kw) if t.name == name else t for t in model.traps
    )
    return dataclasses.replace(model, traps=tr)


# ----------------------------------------------------------------- candidate models
burin_k = dataclasses.replace(DEFAULT_TRAP_MODEL, occupancy="kinetic")
lit = literature_kinetic()
EH4H = with_trap(lit, "EH4", E_below_Ec=1.04, sigma_n=6e-14)
CANDIDATES = {
    "M0 Burin table, kinetic": burin_k,
    "M1 + Z1/Z2 negative-U (2.5 V_C/neq)": literature_kinetic(),
    "M1b + Z1/Z2 negative-U (5.0 V_C/neq)": literature_kinetic(5.0),
    "M2 EH4 as donor (Karsthof CAV +/0)": with_trap(lit, "EH4", kind="donor"),
    "M3a EH4 Hazdra pair (1.04 eV, 6e-14)": EH4H,
    "M3b M3a + EH5 acceptor (1.10 eV, 5e-15; rate ASSUMED = EH4)": dataclasses.replace(
        EH4H, traps=EH4H.traps + (DeepTrap("EH5", "acceptor", 2.4, 1.10, 5e-15, 5e-14),)
    ),
    "M4 EH6/7 as acceptor": with_trap(lit, "EH6/7", kind="acceptor"),
    "M5 + direct donor removal 3.6/neq": dataclasses.replace(lit, donor_removal=3.6),
    "M6 + deep acceptor 3.6/neq": dataclasses.replace(
        lit, traps=lit.traps + (deep_acceptor_trap(3.6),)
    ),
}  # fmt: skip
print("== 1. candidate mechanisms: Radmanovac slopes (per p/cm2), flat-C-V fluences")
print("   measured g: onsemi 5.60+-0.81, W2 4.60+-0.33, W4 4.48+-0.26")
rows = {}
for name, m in CANDIDATES.items():
    g = [removal_observable(OBS[k], m) for k in ("onsemi", "W2", "W4")]
    fl = flat_table(m)
    rows[name] = (g, fl)
    print(f"  {name}")
    print(f"     g onsemi {g[0]:5.2f}  W2 {g[1]:5.2f}  W4 {g[2]:5.2f}  "
          f"neutral removal {m.neutral_removal_rate:5.2f}/neq")  # fmt: skip
    print("     flat at: " + "; ".join(
        f"{k.split(' ')[0]} {v:.2e} {flat_verdict(k, v)}" for k, v in fl.items()))  # fmt: skip
print(f"[{time.time() - T0:.0f} s]")

# ----------------------------------------------------------------- depth shape
print(
    "\n== 2. depth shape: removal at W 11-15 um / at 36-40 um (onsemi, 1.4e12, 2.8e12)"
)
print(
    f"   measured (VEC Fig. 4b, 13 / 38 um): PIN4 {DEPTH['PIN4'][1]:.2f}, PIN5 {DEPTH['PIN5'][1]:.2f}"
)
for name in ("M0 Burin table, kinetic", "M5 + direct donor removal 3.6/neq",
             "M3a EH4 Hazdra pair (1.04 eV, 6e-14)"):  # fmt: skip
    for rate in (1.0, 0.1):
        r = depth_ratio(CANDIDATES[name], rate)
        print(f"  {name[:40]:40s} {rate:4.1f} V/s: {np.array2string(r, precision=2)}")
print(f"[{time.time() - T0:.0f} s]")

# ----------------------------------------------------------------- fits
print("\n== 3. calibration (linear parameters on literature_kinetic): g of 3 devices")
obs3 = [OBS["onsemi"], OBS["W2"], OBS["W4"]]
flat_on = RemovalObservation(
    "onsemi flat (2.8e12, 6.9e12]", "flat_fluence",
    0.5 * (np.log10(2.8e12) + np.log10(6.9e12)),
    0.5 * (np.log10(6.9e12) - np.log10(2.8e12)), 4e13, 50.0, 252.7, freq_Hz=1e6,
    source="Radmanovac Fig. 4b (VEC)", niel_convention="si-neq",
)  # fmt: skip
fits = {
    "d only": (obs3, ("donor_removal",)),
    "EH4 scale only": (obs3, ("EH4",)),
    "d + EH4 (3 g)": (obs3, ("donor_removal", "EH4")),
    "d + EH4 + Z (3 g + onsemi flat)": (
        obs3 + [flat_on],
        ("donor_removal", "EH4", "Z1/2"),
    ),
    "d + Z (3 g + onsemi flat)": (obs3 + [flat_on], ("donor_removal", "Z1/2")),
}
cal_d = None
for name, (obs, params) in fits.items():
    c = calibrate_removal(obs, params, base_model=lit, max_nfev=20)
    if name == "d only":
        cal_d = c
    if name == "EH4 scale only":
        s_eh4 = float(c.scales[0])
    sv = np.array(c.settings["singular_values_sorted"])
    print(f"  {name}: " + ", ".join(
        f"{p} {v:.3g} +- {s:.2g}" for p, v, s in zip(params, c.scales, c.sigma))
        + f"; chi2 {c.chi2:.2f} (dof {c.dof}); model {np.round(c.settings['model_values'], 2)}")  # fmt: skip
    if len(params) > 1:
        print(f"     corr {np.round(c.correlation[np.triu_indices(len(params), 1)], 3)}"
              f"; singular values {np.array2string(sv, precision=3)}"
              f"; weakest direction {np.round(c.settings['null_directions'][0], 3)}")  # fmt: skip
    for w in c.warnings:
        print(f"     warning: {w}")
print(f"[{time.time() - T0:.0f} s]")

# ----------------------------------------------------------------- predictions
d_fit = float(cal_d.scales[0])
cal = dataclasses.replace(lit, donor_removal=d_fit)
cal_da = dataclasses.replace(lit, traps=lit.traps + (deep_acceptor_trap(d_fit),))
print(
    f"\n== 4. PREDICTIONS of the calibrated model (d = {d_fit:.2f} +- {cal_d.sigma[0]:.2f} /neq)"
)
for k, v in flat_table(cal).items():
    print(
        f"  flat C-V {k}: {v:.2e} ({flat_verdict(k, v)}; bracket {FLAT[k]['lo']:.2g}..{FLAT[k]['hi']:.2g})"
    )
print(f"  depth ratio onsemi (1 V/s): {np.array2string(depth_ratio(cal), precision=2)}")
print(f"  depth ratio onsemi, EH4-scale fit (1 / 0.1 V/s): "
      f"{np.array2string(depth_ratio(parametrized_model(lit, ('EH4',), [s_eh4])), precision=2)} / "
      f"{np.array2string(depth_ratio(parametrized_model(lit, ('EH4',), [s_eh4]), 0.1), precision=2)}")
m_eh4 = parametrized_model(lit, ("EH4",), [s_eh4])
for T in (300.0, 325.0, 350.0):
    o = dataclasses.replace(OBS["onsemi"], T=T, fluences=(0.0, 1.4e12, 2.8e12))
    print(f"  T {T:.0f} K: g onsemi: direct removal {removal_observable(o, cal):.2f}, "
          f"deep acceptor {removal_observable(o, cal_da):.2f}, "
          f"EH4-scale fit (s {s_eh4:.2f}) {removal_observable(o, m_eh4):.2f}")  # fmt: skip
print(f"\ntotal {time.time() - T0:.0f} s")
