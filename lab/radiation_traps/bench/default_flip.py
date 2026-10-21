"""Before / after the 2026-10-21 default flip (issues #5, #8): NIEL convention si-neq ->
sic-damage, damage model legacy scalar removal -> calibrated kinetic traps ("auto").

uv run python lab/radiation_traps/bench/default_flip.py   (~10-30 s warm)

Default stack (DeviceConfig()), 62 MeV protons, -40 V. Prints the compensation fluences
and the CCE at 0.1 / 0.5 / 0.9 / 1.5 / 10 x Phi_comp of each model (its own Phi_comp),
plus the wall time of the default run_radiation_damage sweep (cold = first call in the
process, incl. numba cache load; warm = second call).
"""

from __future__ import annotations

import dataclasses
import logging
import time

import numpy as np

from etna.api.device import DeviceConfig
from etna.api.simulation import compensation_fluence, run_radiation_damage

logging.disable(logging.WARNING)
E, V = 62.0, -40.0
FR = np.array([0.1, 0.5, 0.9, 1.5, 10.0])
cfg = DeviceConfig()
CASES = {
    "OLD default (legacy removal, si-neq)": (dataclasses.replace(cfg, radiation_traps=False), "si-neq"),
    "legacy removal, sic-damage": (dataclasses.replace(cfg, radiation_traps=False), "sic-damage"),
    "calibrated traps, si-neq (d 2.87)": (cfg, "si-neq"),
    "NEW default (calibrated traps 'auto', sic-damage)": (cfg, None),
}  # fmt: skip
t0 = time.time()
r = run_radiation_damage(cfg, proton_energy_MeV=E, V_bias=V)
t_cold = time.time() - t0
t0 = time.time()
r = run_radiation_damage(cfg, proton_energy_MeV=E, V_bias=V)
t_warm = time.time() - t0
print(
    f"default sweep ({r.x.size} points, {r.metadata['damage_model']}): cold {t_cold:.2f} s, "
    f"warm {t_warm:.2f} s"
)
for name, (c, conv) in CASES.items():
    phi = compensation_fluence(c, E, niel_convention=conv)
    t0 = time.time()
    res = run_radiation_damage(c, FR * phi, V, E, niel_convention=conv)
    dt = time.time() - t0
    print(
        f"{name:50s} Phi_comp {phi:.3e} p/cm2 ({res.metadata['damage_model']}, kappa "
        f"{res.metadata['kappa']:.4f}); CCE at {FR.tolist()} x Phi_comp: "
        f"{np.array2string(np.asarray(res.y), precision=6)}  [{dt:.2f} s]"
    )
# the same absolute fluences for all models (fractions of the OLD Phi_comp)
phi_old = compensation_fluence(
    CASES["legacy removal, sic-damage"][0], E, niel_convention="si-neq"
)
for name, (c, conv) in CASES.items():
    res = run_radiation_damage(c, FR * phi_old, V, E, niel_convention=conv)
    print(
        f"  at x OLD Phi_comp {phi_old:.3e}: {name:50s} "
        f"{np.array2string(np.asarray(res.y), precision=6)}"
    )
