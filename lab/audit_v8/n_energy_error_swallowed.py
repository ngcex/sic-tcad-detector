"""Out-of-table energy in the core cce_vs_fluence: raised or swallowed per point?"""

import numpy as np
from etna import DeviceConfig
from etna.api.device import device_kwargs
from etna.core.charge_collection import cce_vs_fluence

kw = device_kwargs(DeviceConfig())
eth = kw.pop("epi_thickness_cm")
for traps in (None, "auto"):
    r = cce_vs_fluence(
        np.array([0.0, 1e12]),
        V_bias=-40.0,
        epi_thickness_cm=eth,
        energy_MeV=2000.0,
        device_kwargs=kw,
        radiation_traps=traps,
    )
    print(traps, r["cce_values"], r["valid"], r["reason"])
