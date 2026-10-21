"""Issue #13: mesh-refinement check of the default-path CCE (occupancy lifetimes) with a
legacy-lifetime control. Scales the spacing of the build_diode_mesh lines with
PSMIN <= ps <= PSMAX by FACT (devsim.add_1d_mesh_line patched), default stack, 62 MeV,
-40 V, 1e13 / 1e14 / 1e15 p/cm^2.

uv run python lab/audit_v8/m2_mesh_refinement.py FACT [PSMIN [PSMAX]]
  0.25 0 1e-6 : 4x finer 1-nm lines (junction, n-/n+ step)
  0.5         : every spacing halved

Results 2026-10-21 (dCCE refined - base at 1e13 / 1e14 / 1e15):
  0.25 0 1e-6, node mean (before #13): -2.9e-5 / -1.0e-4 / -2.7e-5
  0.25 0 1e-6, epi length mean      : -3.5e-7 / -5.7e-7 / +9.9e-7
  0.25 0 1e-6, legacy lifetimes     : +1.0e-7 / +7.5e-7 / +4.1e-6
  0.5,         epi length mean      : +4.2e-6 / +1.7e-5 / +1.2e-4
  0.5,         legacy lifetimes     : +7.7e-6 / +3.8e-5 / +2.4e-4
"""

import dataclasses
import logging
import sys
import time

import numpy as np

import devsim
import etna.core.device as dev
from etna import DeviceConfig, run_radiation_damage
from etna.core.trap_presets import trap_preset

logging.disable(logging.WARNING)
FACT = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
PSMIN = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0  # scale only lines ps >= PSMIN
PSMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
orig = devsim.add_1d_mesh_line


def scaled(**kw):
    kw = dict(kw)
    if PSMIN <= kw["ps"] <= PSMAX:
        kw["ps"] = kw["ps"] * FACT
    return orig(**kw)


Fp = np.array([1e13, 1e14, 1e15])
cfg = DeviceConfig()
leg = dataclasses.replace(
    cfg,
    radiation_traps=dataclasses.replace(
        trap_preset("radmanovac_calibrated"), lifetimes="legacy"
    ),
)
out = {}
for refine in (False, True):
    dev.devsim.add_1d_mesh_line = scaled if refine else orig
    t = time.time()
    a = run_radiation_damage(cfg, fluences=Fp, V_bias=-40.0, proton_energy_MeV=62.0)
    b = run_radiation_damage(leg, fluences=Fp, V_bias=-40.0, proton_energy_MeV=62.0)
    out[refine] = (np.asarray(a.y), np.asarray(b.y))
    print(
        "refine",
        refine,
        "occupancy",
        np.array2string(a.y, precision=7),
        "legacy-tau",
        np.array2string(b.y, precision=7),
        f"{time.time()-t:.1f}s",
    )
dev.devsim.add_1d_mesh_line = orig
d_occ = out[True][0] - out[False][0]
d_leg = out[True][1] - out[False][1]
print("dCCE occupancy:", d_occ)
print("dCCE legacy   :", d_leg)
print("d(1-CCE)/(1-CCE) occ:", d_occ / (1 - out[False][0]))
