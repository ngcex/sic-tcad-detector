"""Kinetic (non-steady-state) trap occupancy: verification and literature comparison.

uv run python lab/radiation_traps/bench/kinetics.py   (~1-2 min; numbers in ../notes/KINETICS.md)

Model: etna.core.trap_kinetics (opt-in), calibration etna.core.trap_calibration.
Literature numbers (full text opened; GitHub issue #5 lists the sources):
- Radmanovac et al., arXiv:2510.11304: 252.7 MeV p, onsemi 50 um 4e13, C-V at 1 MHz,
  donor removal g = 5.60 +- 0.81 per p/cm^2 at depletion widths 10-40 um.
- Gsponer et al., arXiv:2310.02047: CNM 50 um 1.5e14, 3x3 mm2, C-V at 10 kHz flat for
  5e14-1e16 neq (geometric ~15.5 pF), forward I < 1 nA at 1 kV forward at 1e15 neq,
  reverse < 10 pA at 1.1 kV (measurement floor).
- Kroll et al., arXiv:2607.11663: onsemi 50 um ~6e13, C bias-independent from ~5e13 p/cm2
  of 24 GeV/c protons (kappa 0.62 Si: ~3.1e13 neq).
- Li et al., arXiv:2503.09016: 100 um 5.2e13, C-V flat at 3.9e13 neq (10 kHz).
"""
# NOTE 2026-10-21: the defaults changed (NIEL convention sic-damage, damage model
# "auto" = calibrated kinetic traps). The numbers logged from this bench predate the flip
# (si-neq); pass niel_convention="si-neq" / niel_table=hardness_table("proton", "si-neq")
# and radiation_traps=False or True explicitly to reproduce them.


import dataclasses
import logging
import time

import numpy as np

from etna.api.device import DeviceConfig
from etna.core.batched.material import EPS
from etna.core.charge_collection import _trap_dark_states
from etna.core.radiation_traps import DEFAULT_TRAP_MODEL, RadiationTrapModel
from etna.core.trap_calibration import scaled_model
from etna.core.trap_kinetics import (
    KineticTrapSolver,
    cv_protocol,
    hold,
    neff_from_cv,
    simulate_cv,
)

logging.disable(logging.WARNING)
T0 = time.time()


def neff_window(V, C, lo=10.0, hi=40.0):
    """Median apparent N_eff over depletion widths lo..hi um (Radmanovac method)."""
    W = EPS / C * 1e4
    ne = neff_from_cv(V, C)
    out = []
    for b in range(C.shape[0]):
        m = (W[b] > lo) & (W[b] < hi)
        out.append(np.median(ne[b, m]) if m.any() else np.nan)
    return np.array(out)


def slope(F, ne):
    """-dN_eff/dPhi by least squares over the finite (uncompensated) points."""
    m = np.isfinite(ne)
    return -np.polyfit(np.asarray(F)[m], ne[m], 1)[0]


# ----------------------------------------------------------------------- 1. limits
print("== 1. kernel limits (default stack, 62 MeV, -40 V)")
F62 = np.array([1e11, 1e12, 3e12, 1e13, 1e14])
s = KineticTrapSolver.from_config(DeviceConfig(), F62)
st = s.steady_state(-40.0)
dark = _trap_dark_states(s.ref, s.fluences_neq, -40.0, s.model)
print(
    f"steady limit vs radiation_traps kernel: max|dpsi| {np.abs(st.psi - dark['res'].psi).max():.1e} V"
)

print(f"[{time.time() - T0:.1f} s]")
# ----------------------------------------------------------------------- 2. tol study
print("\n== 2. adaptive step: C-V error vs tol (20 um, 1e14, 0..-40 V at 1 V/s)")
u = KineticTrapSolver.from_config(
    DeviceConfig(epi_thickness_um=20.0, doping_profile="uniform", N_D=1e14),
    [0.0, 3e12, 6e12, 3e13],
    energy_MeV=None,
)
Vt = -np.arange(0.0, 41.0, 4.0)


def cvrun(tol, norm):
    r = u.run(
        u.initial_state("equilibrium"),
        cv_protocol(Vt, 1.0),
        tol=tol,
        dV_max=1.0,
        norm=norm,
    )
    return np.array([u.capacitance_qs(x)[0] for x in r.states]).T, r


Cref, rr = cvrun(1e-5, "charge")
print(f"reference tol 1e-5: {rr.n_steps} steps")
for norm in ("charge", "occupancy"):
    for tol in (1e-2, 1e-3, 1e-4):
        t0 = time.time()
        C, r = cvrun(tol, norm)
        print(
            f"  {norm:9s} tol {tol:.0e}: {time.time() - t0:5.2f} s, {r.n_steps:5d} steps "
            f"({r.n_rejected} rejected), max|dC/C| {np.max(np.abs(C / Cref - 1)):.1e}"
        )

print(f"[{time.time() - T0:.1f} s]")
# ----------------------------------------------------------------------- 3. removal
print(
    "\n== 3. C-V donor removal, Radmanovac-like device (50 um, 4e13, 252.7 MeV, 1 MHz)"
)
cfgR = DeviceConfig(epi_thickness_um=50.0, doping_profile="uniform", N_D=4e13)
FR = np.array([0.0, 1e12, 2e12, 3e12, 4e12])
VR = -np.arange(0.0, 151.0, 5.0)
acc67 = RadiationTrapModel(
    tuple(
        dataclasses.replace(t, kind="acceptor") if t.name == "EH6/7" else t
        for t in DEFAULT_TRAP_MODEL.traps
    )
)
cases = [
    ("legacy removal (eta 5/neq)", None),
    ("steady state", dict(occupancy="steady")),
    ("kinetic, equilibrium start", dict(occupancy="kinetic", initial="equilibrium")),
    ("kinetic, in-beam start + 10 min at 0 V",
     dict(occupancy="kinetic", initial="in_beam", storage=(hold(600.0, 0.0),))),
    ("kinetic, in-beam start + 1 day at 0 V",
     dict(occupancy="kinetic", initial="in_beam", storage=(hold(86400.0, 0.0),))),
]  # fmt: skip
sR = KineticTrapSolver.from_config(cfgR, FR, energy_MeV=252.7)
kappa = sR.fluences_neq[1] / FR[1]
print(f"kappa(252.7 MeV) = {kappa:.3f}; measured g = 5.60 +- 0.81 /p")
for name, kw in cases:
    if kw is None:
        print(f"  {name:42s}: g = {5.0 * kappa:.2f} /p (scalar, every depth)")
        continue
    cv = simulate_cv(sR, VR, 1e6, rate_V_per_s=1.0, dV_max=5.0, **kw)
    ne = neff_window(VR, cv.C_ac)
    print(
        f"  {name:42s}: g = {slope(FR, ne):5.2f} /p   N_eff {np.array2string(ne, precision=3)}"
    )
sA = KineticTrapSolver.from_config(cfgR, FR, energy_MeV=252.7, trap_model=acc67)
cv = simulate_cv(sA, VR, 1e6, rate_V_per_s=1.0, dV_max=5.0)
print(
    f"  {'ALTERNATIVE table: EH6/7 acceptor, kinetic':42s}: g = {slope(FR, neff_window(VR, cv.C_ac)):5.2f} /p"
)
G_MEAS, G_SIG = 5.60, 0.81


def g_scaled(grp, sc, cfg=cfgR):
    m = scaled_model(DEFAULT_TRAP_MODEL, grp, [sc])
    so = KineticTrapSolver.from_config(cfg, FR[:3], energy_MeV=252.7, trap_model=m)
    cv = simulate_cv(so, VR, 1e6, rate_V_per_s=1.0, dV_max=5.0)
    return slope(FR[:3], neff_window(VR, cv.C_ac))


g1 = g_scaled(("all",), 1.0)  # unscaled, on the 3-fluence subset (linear regime)
sens = {}
for grp in (("all",), ("EH4",), ("Z1/2",)):
    g2 = g_scaled(grp, 2.0)
    sens[grp[0]] = g2 - g1  # dg / ds (g is linear in s below compensation)
    if grp[0] == "Z1/2":
        continue
    s_fit = 1.0 + (G_MEAS - g1) / sens[grp[0]]
    s_sig = G_SIG / abs(sens[grp[0]])
    nrem = sum(
        t.eta
        for t in scaled_model(DEFAULT_TRAP_MODEL, grp, [s_fit]).traps
        if t.kind == "acceptor"
    )
    print(
        f"  CALIBRATED scale {grp[0]!r} to g = 5.60 +- 0.81: s = {s_fit:.2f} +- {s_sig:.2f}"
        f"  -> PREDICTED neutral-bulk removal {nrem * kappa:.1f} /p, compensation of 4e13 at"
        f" {4e13 / (nrem * kappa):.2e} p/cm2 (N0/g_measured = {4e13 / 5.6:.2e})"
    )
print(
    f"  sensitivities dg/ds (/p per unit scale, 3-fluence subset, g(1) = {g1:.2f}): "
    + ", ".join(f"{k} {v:.2f}" for k, v in sens.items())
)
# two-number reading (Z1/2, EH4 scales) of the same measured g: (i) the C-V slope,
# linear in both scales, g = g1 + a4 (s4 - 1) + aZ (sZ - 1); (ii) the compensation
# fluence N0 / g, i.e. neutral removal kappa (5 sZ + 2.4 s4) = g, ASSUMING the measured
# linear N_eff(Phi) extrapolates to compensation (as CNM W4: 8.8e13 p/cm2 ~ N0 / g in
# Radmanovac et al.). Both rows use the same g: sigma by dg propagation (correlated).
a4, aZ = sens["EH4"], sens["Z1/2"]
M = np.array([[a4, aZ], [2.4 * kappa, 5.0 * kappa]])
rhs = np.array([G_MEAS - g1 + a4 + aZ, G_MEAS])
s4, sZ = np.linalg.solve(M, rhs)
ds = np.linalg.solve(M, np.array([1.0, 1.0])) * G_SIG
print(
    f"  2-number reading: s_EH4 = {s4:.2f} +- {abs(ds[0]):.2f}, s_Z1/2 = {sZ:.2f} +- "
    f"{abs(ds[1]):.2f} (assumption-laden, not a fit)"
)
print(
    "\n== 3b. PREDICTED temperature dependence of the kinetic C-V removal (same device)"
)
for T in (300.0, 325.0, 350.0):
    cT = DeviceConfig(epi_thickness_um=50.0, doping_profile="uniform", N_D=4e13, T=T)
    print(
        f"  T {T:5.1f} K: g = {g_scaled(('all',), 1.0, cT):.2f} /p (EH4 emission speeds up)"
    )
print(f"[{time.time() - T0:.1f} s]")
# ----------------------------------------------------------------------- 4. history
print(
    "\n== 4. PREDICTED history dependence (same device, 2e12 p/cm2): up sweep, 1 h at -150 V, down sweep"
)
sH = KineticTrapSolver.from_config(cfgR, [0.0, 2e12], energy_MeV=252.7)
VH = -np.arange(0.0, 151.0, 5.0)
r_up = sH.run(
    sH.initial_state("equilibrium"), cv_protocol(VH, 1.0), tol=1e-3, dV_max=5.0
)
r_hold = sH.run(r_up.final, [hold(3600.0, VH[-1], record=0)], tol=1e-3)
r_dn = sH.run(r_hold.final, cv_protocol(VH[::-1], 1.0), tol=1e-3, dV_max=5.0)


def cvs(r):
    Cq = np.array([sH.capacitance_qs(x)[0] for x in r.states]).T
    return np.array(
        [sH.capacitance_ac(x, 1e6, Cq[:, j]) for j, x in enumerate(r.states)]
    ).T


ne_up = neff_window(VH, cvs(r_up))
ne_dn = neff_window(VH[::-1], cvs(r_dn)[:, ::-1][:, ::-1])
ne_ss = neff_window(VH, simulate_cv(sH, VH, 1e6, occupancy="steady").C_ac)
print(
    f"  N_eff up {ne_up[1]:.3e}, after 1 h hold (down) {ne_dn[1]:.3e}, steady state {ne_ss[1]:.3e}, pristine {ne_up[0]:.3e}"
)

print(f"[{time.time() - T0:.1f} s]")
# ----------------------------------------------------------------------- 5. flattening
print(
    "\n== 5. C-V flattening, Gsponer-like device (50 um, 1.5e14, 3x3 mm2, 10 kHz; neq)"
)
A = 0.09
cfgG = DeviceConfig(epi_thickness_um=50.0, doping_profile="uniform", N_D=1.5e14)
FG = np.array([0.0, 1e13, 2e13, 3e13, 5e14, 1e15, 1e16])
sG = KineticTrapSolver.from_config(cfgG, FG, energy_MeV=None)
VG = -np.array([0.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0])
cvG = simulate_cv(sG, VG, 1e4, rate_V_per_s=2.0, dV_max=5.0)
print(f"  geometric eps A / 50 um = {EPS / 50e-4 * A * 1e12:.2f} pF (paper: ~15.5 pF)")
for b, f in enumerate(FG):
    C = cvG.C_ac[b] * A * 1e12
    print(
        f"  {f:8.1e} neq: C(0) {C[0]:7.2f} pF, C(-500) {C[-1]:6.2f} pF, max/min {C.max() / C.min():.4f}"
    )
R = sG.flatband_resistance()
print(
    "  flat-band epi resistance (Ohm cm2) and ohmic bound of I at 1 kV forward (3x3 mm2):"
)
for f, r in zip(FG, R):
    print(f"    {f:8.1e} neq: R {r:9.2e}, I_max {1000.0 / r * A:9.2e} A")
print("  measured: I_fwd < 1 nA at 1 kV at 1e15 neq; = reverse (< 10 pA) at 1e16 neq")
print(f"  trap SRH generation at -500 V (A, 3x3 mm2): {np.array2string(cvG.J_gen[:, -1] * A, precision=2)}"
      " (measured floor ~1e-11 A: not testable)")  # fmt: skip

print(f"[{time.time() - T0:.1f} s]")
# ----------------------------------------------------------------------- 6. onset
print(
    "\n== 6. flat-C-V onset vs doping (kinetic, 10 kHz, 50 um; max/min C over 0..-100 V < 1.02)"
)
for ND, meas in ((6.2e13, "Kroll: flat from ~3.1e13 neq (5e13 p 24 GeV/c, kappa 0.62)"),
                 (1.5e14, "Gsponer: flat at 5e14 neq (lowest measured)")):  # fmt: skip
    Fs = np.geomspace(4e12, 6.4e13, 5)
    so = KineticTrapSolver.from_config(
        DeviceConfig(epi_thickness_um=50.0, doping_profile="uniform", N_D=ND),
        Fs,
        energy_MeV=None,
    )
    cv = simulate_cv(so, [0.0, -20.0, -50.0, -100.0], 1e4, rate_V_per_s=2.0, dV_max=5.0)
    ratio = cv.C_ac.max(axis=1) / cv.C_ac.min(axis=1)
    on = Fs[np.argmax(ratio < 1.02)] if (ratio < 1.02).any() else np.nan
    print(
        f"  N_D {ND:.2e}: model flat from {on:.2e} neq (N_D / 7.4 = {ND / 7.4:.2e}); {meas}"
    )

print(f"\ntotal {time.time() - T0:.1f} s")
