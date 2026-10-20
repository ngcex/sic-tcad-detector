"""N_t refit with the Hurkx Gamma: I(300 K, -60 V, 4 mm2) = 18 pA per preset (devsim)."""

import dataclasses

import numpy as np

import etna.api.simulation as sim
from etna.api.device import (
    DeviceConfig,
    legacy_config,
    measured_cv_tabulated_config,
    microdosimeter_config,
    nplus_calibrated_config,
)

PRESETS = {
    "default": DeviceConfig,
    "tabulated": measured_cv_tabulated_config,
    "nplus_calibrated": nplus_calibrated_config,
    "legacy": legacy_config,
    "microdosimeter": microdosimeter_config,
}


def I60(cfg):
    r = sim.run_dark_current(cfg, 0.0, -60.0, 31, solver="devsim")
    return float(r.y[-1])


for name, mk in PRESETS.items():
    cfg = mk()
    I0 = I60(cfg)
    Nt = cfg.N_t * 18e-12 / I0
    for _ in range(3):
        I1 = I60(dataclasses.replace(cfg, N_t=Nt))
        if abs(I1 - 18e-12) < 1e-18:
            break
        Nt *= 18e-12 / I1
    print(
        f"RESULT {name}: area {cfg.area_cm2} N_t {cfg.N_t:.6e} I60 {I0*1e12:.4f} pA -> "
        f"N_t {Nt:.6e} I60 {I1*1e12:.6f} pA"
    )
