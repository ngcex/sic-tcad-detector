"""Does any API path use a calibrated preset under another NIEL convention?"""

import dataclasses
import numpy as np
from etna import DeviceConfig, run_radiation_damage
from etna.api.simulation import (
    trap_compensation_fluence,
    compensation_fluence,
    default_fluence_grid,
)
from etna.core.trap_presets import trap_preset

m_sic = trap_preset("radmanovac_calibrated")  # built in the default (sic-damage)
print("preset convention", m_sic.niel_convention, "d", m_sic.donor_removal)
cfg = DeviceConfig(radiation_traps=m_sic)
# A) trap_compensation_fluence with an explicit model and another convention
print(
    "trap_comp si-neq (model sic-damage):",
    trap_compensation_fluence(cfg, 62.0, niel_convention="si-neq"),
)
m_si = trap_preset("radmanovac_calibrated", niel_convention="si-neq")
print(
    "trap_comp si-neq (model si-neq)    :",
    trap_compensation_fluence(
        DeviceConfig(radiation_traps=m_si), 62.0, niel_convention="si-neq"
    ),
)
# B) compensation_fluence (goes through resolve_damage_model -> check)
for fn in (compensation_fluence, default_fluence_grid):
    try:
        fn(cfg, 62.0, niel_convention="si-neq")
        print(fn.__name__, "NO raise")
    except ValueError as e:
        print(fn.__name__, "raises:", str(e)[:80])
try:
    run_radiation_damage(
        cfg,
        fluences=np.array([1e12]),
        V_bias=-40.0,
        proton_energy_MeV=62.0,
        niel_convention="si-neq",
    )
    print("run_radiation_damage NO raise")
except ValueError as e:
    print("run_radiation_damage raises:", str(e)[:80])
