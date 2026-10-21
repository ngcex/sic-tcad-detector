"""CCE vs fluence beyond compensation: steady vs kinetic occupancy, occupancy-weighted
lifetimes, the calibrated slow-removal component, hold / initial-state brackets and the
EH6/7 sigma_n choice.

uv run python lab/radiation_traps/bench/cce_kinetic.py   (numbers -> ../notes/MECHANISMS.md)

Default stack (10 um epi), 62 MeV protons (kappa 1.942), alpha CCE via
run_radiation_damage. All PREDICTED; only the slow-removal rate of the
radmanovac_calibrated presets is CALIBRATED (to 252.7 MeV C-V slopes, MECHANISMS.md).
"""

from __future__ import annotations

import dataclasses
import logging
import time

import numpy as np

from etna.api.device import DeviceConfig
from etna.api.simulation import run_radiation_damage
from etna.core.radiation_damage import RadiationDamageParams, compute_K_tau
from etna.core.radiation_traps import DEFAULT_TRAP_MODEL
from etna.core.trap_presets import literature_kinetic, radmanovac_calibrated

logging.disable(logging.WARNING)
T0 = time.time()
F = np.array([1e11, 1e12, 3e12, 1e13, 1e14, 3e14, 1e15, 5e15])


def cce(model, V=-40.0):
    r = run_radiation_damage(
        DeviceConfig(radiation_traps=model),
        fluences=F,
        V_bias=V,
        proton_energy_MeV=62.0,
    )
    return np.asarray(r.y, dtype=float)


legacy = run_radiation_damage(
    DeviceConfig(), fluences=F, V_bias=-40.0, proton_energy_MeV=62.0
)
print("fluences (p/cm2, 62 MeV):", np.array2string(F, precision=0))
print(
    f"  legacy removal          -40 V: {np.array2string(np.asarray(legacy.y, float), precision=3)}"
)
lit = literature_kinetic()
cal = radmanovac_calibrated()
cal_da = radmanovac_calibrated(True)
MODELS = {
    "burin steady (default traps)": DEFAULT_TRAP_MODEL,
    "burin kinetic, legacy tau": dataclasses.replace(DEFAULT_TRAP_MODEL, occupancy="kinetic"),
    "literature_kinetic, legacy tau": dataclasses.replace(lit, lifetimes="legacy"),
    "literature_kinetic, occ. tau": lit,
    "calibrated direct, occ. tau": cal,
    "calibrated direct, legacy tau": dataclasses.replace(cal, lifetimes="legacy"),
    "calibrated deep acc., occ. tau": cal_da,
}  # fmt: skip
for V in (-40.0, -200.0):
    print(f"\n== CCE at {V:g} V")
    for name, m in MODELS.items():
        print(f"  {name:31s}: {np.array2string(cce(m, V), precision=3)}")
    print(f"[{time.time() - T0:.0f} s]")

print("\n== hold / initial-state bracket (calibrated direct, occ. tau, -40 V)")
for hold_s in (60.0, 3600.0, 86400.0):
    for init in ("equilibrium", "in_beam"):
        m = dataclasses.replace(cal, hold_s=hold_s, initial=init)
        print(
            f"  hold {hold_s:8.0f} s, {init:11s}: {np.array2string(cce(m), precision=3)}"
        )
print(f"[{time.time() - T0:.0f} s]")

print(
    "\n== EH6/7 sigma_n: legacy K_tau (electron, 300 K) and CCE (burin steady, -40 V)"
)
for sn in (1.4e-14, 1e-13, 3e-13, 9e-12):
    K = compute_K_tau(dataclasses.replace(RadiationDamageParams(), sigma_n_EH67=sn))
    tr = tuple(
        dataclasses.replace(t, sigma_n=sn) if t.name == "EH6/7" else t
        for t in DEFAULT_TRAP_MODEL.traps
    )
    m = dataclasses.replace(DEFAULT_TRAP_MODEL, traps=tr)
    print(f"  sigma_n {sn:.1e}: K_n {K:.3e} cm2/s; CCE steady traps (legacy tau, K_n of the "
          f"default params) {np.array2string(cce(m), precision=3)}")  # fmt: skip
for sn in (1.4e-14, 1e-13, 9e-12):
    m = literature_kinetic(eh67_sigma_n=sn)
    print(
        f"  literature_kinetic occ. tau, sigma_n {sn:.1e}: {np.array2string(cce(m), precision=3)}"
    )
print(f"\ntotal {time.time() - T0:.0f} s")
