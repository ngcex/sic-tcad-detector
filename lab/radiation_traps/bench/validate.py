"""Validation numbers of the opt-in charged-trap damage model (etna.core.radiation_traps).

uv run python lab/radiation_traps/bench/validate.py   (~10 s; numbers in ../notes/LOG.md)

1. Local bulk physics (trap_diagnostics) of the default-stack bulk N_D vs Phi_neq.
2. Default stack (10 um, 62 MeV, alpha), CCE vs fluence at -40 V: legacy removal vs traps.
3. Bias dependence beyond compensation (default stack).
4. Literature-shaped devices (qualitative: different geometry / particle / readout):
   Gsponer et al. arXiv:2310.02047 (50 um, N_D ~1.5e14, CCE > 80 % to 1e15 neq) and the
   De Napoli 2022 review (20 um, 1e14: losses 7 % / 22.5 % at 1.3e14 / 7.3e14 cm^-2,
   14 MeV n; normalisation of those fluences to neq not checked).
"""
# NOTE 2026-10-21: the defaults changed (NIEL convention sic-damage, damage model
# "auto" = calibrated kinetic traps). The numbers logged from this bench predate the flip
# (si-neq); pass niel_convention="si-neq" / niel_table=hardness_table("proton", "si-neq")
# and radiation_traps=False or True explicitly to reproduce them.


import logging

import numpy as np

from etna.api.device import DeviceConfig, device_kwargs
from etna.core.charge_collection import cce_vs_fluence
from etna.core.radiation_damage import get_hardness_factor
from etna.core.radiation_traps import trap_diagnostics

logging.disable(logging.WARNING)
KAPPA = get_hardness_factor(62.0)


def sweep(cfg, fluences_p, V, traps):
    kw = device_kwargs(cfg)
    epi = kw.pop("epi_thickness_cm")
    extra = {"radiation_traps": True} if traps else {}
    return cce_vs_fluence(
        np.asarray(fluences_p),
        V_bias=V,
        epi_thickness_cm=epi,
        device_kwargs=kw,
        **extra,
    )


print("== 1. bulk (N_D = 3.38e13)")
for phi in (0, 1e12, 4e12, 1e13, 1e14, 1e15, 1e16):
    d = trap_diagnostics(DeviceConfig().N_D_bulk, phi)
    print(
        f"Phi_neq {phi:8.1e}: n {d['n_neutral']:9.3e}  Ec-EF {d['Ec_minus_EF']:.3f} eV  "
        f"rho {d['resistivity_ohm_cm']:8.2e} Ohm cm  N_eff_dep(ss) "
        f"{d['N_eff_depleted_ss']:9.3e}  (frozen {d['N_eff_depleted_frozen']:9.2e})"
    )
print("emission times (s):", trap_diagnostics(3.38e13, 1e14)["emission_time_s"])

print("\n== 2. default stack, -40 V, 62 MeV (kappa %.3f)" % KAPPA)
F = np.array([0, 1e10, 1e11, 1e12, 3e12, 1e13, 3e13, 1e14, 3e14, 1e15, 5e15])
leg = sweep(DeviceConfig(), F, -40.0, False)
tr = sweep(DeviceConfig(), F, -40.0, True)
for f, a, b in zip(F, leg["cce_values"], tr["cce_values"]):
    print(f"Phi_p {f:8.1e} (neq {f * KAPPA:8.1e}): legacy {a:.6f}  traps {b:.6f}")

print("\n== 3. bias dependence, default stack")
for V in (-40.0, -100.0, -200.0):
    r = sweep(DeviceConfig(), [1e14, 1e15, 5e15], V, True)
    print(f"V {V:6.0f}:", "  ".join(f"{c:.4f}" for c in r["cce_values"]))

print("\n== 4a. Gsponer-like: 50 um uniform 1.5e14 (fluences in neq)")
g = DeviceConfig(epi_thickness_um=50.0, doping_profile="uniform", N_D=1.5e14)
neq = np.array([1e13, 1e14, 1e15, 1e16])
for V in (-100.0, -300.0, -600.0):
    r = sweep(g, neq / KAPPA, V, True)
    print(
        f"V {V:6.0f}:",
        "  ".join(f"{n:.0e}:{c:.3f}" for n, c in zip(neq, r["cce_values"])),
    )

print("\n== 4b. De Napoli-like: 20 um uniform 1e14 (fluences taken as neq)")
dn = DeviceConfig(epi_thickness_um=20.0, doping_profile="uniform", N_D=1e14)
neq = np.array([1.3e14, 7.3e14])
for V in (-100.0, -300.0):
    r = sweep(dn, neq / KAPPA, V, True)
    print(
        f"V {V:6.0f}:",
        "  ".join(f"{n:.1e}:{c:.3f}" for n, c in zip(neq, r["cce_values"])),
    )
