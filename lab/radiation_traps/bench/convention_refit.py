"""Re-fit the Radmanovac-calibrated slow removal d in every NIEL convention (issue #5/#8).

uv run python lab/radiation_traps/bench/convention_refit.py   (~1-3 min)

The data are per PROTON (Radmanovac et al. arXiv:2510.11304 Table 2, 252.7 MeV): g =
5.60 +- 0.81 (onsemi), 4.60 +- 0.33 (CNM W2), 4.48 +- 0.26 (CNM W4). The model reads
them as g = kappa_c(252.7) x [G_lit(kappa_c Phi) + d_c] (per neq of convention c), with
G_lit the removal of the neutron-anchored Burin/Capan levels (convention independent
per neq). So d_c must be RE-FITTED per convention; the naive rescaling
d_c = d_si x kappa_si / kappa_c is shown for comparison (it ignores that G_lit does not
scale). Checks printed:
1. the fitted model g per proton is the same in every convention (to the fit
   tolerance): the per-proton calibration is convention independent;
2. what is NOT invariant: the Burin-level share per proton, the 62 MeV per-proton
   neutral-bulk removal rate kappa_c(62) x (7.4 + d_c) and the compensation fluence
   of the default stack.
"""

from __future__ import annotations

import dataclasses
import logging
import time

import numpy as np

from etna.api.device import DeviceConfig
from etna.core.radiation_damage import (
    NIEL_CONVENTIONS,
    get_hardness_factor,
    hardness_table,
)
from etna.core.trap_calibration import (
    RemovalObservation,
    calibrate_removal,
    removal_observable,
)
from etna.core.trap_presets import (
    CALIBRATED_SLOW_REMOVAL_BY_CONVENTION,
    literature_kinetic,
)

logging.disable(logging.WARNING)
T0 = time.time()
RAD = "Radmanovac et al. arXiv:2510.11304 Table 2"
BASE = [
    RemovalObservation(
        "onsemi g", "slope", 5.60, 0.81, 4e13, 50.0, 252.7,
        (0.0, 1.4e11, 6.9e11, 1.4e12, 2.8e12), (10.0, 40.0), 150.0, 5.0, 1e6, source=RAD,
    ),
    RemovalObservation(
        "CNM W2 g", "slope", 4.60, 0.33, 7.8e13, 50.0, 252.7,
        (0.0, 1.4e11, 6.9e11, 1.4e12, 2.8e12, 6.9e12), (30.0, 40.0), 200.0, 5.0, 1e6,
        source=RAD,
    ),
    RemovalObservation(
        "CNM W4 g", "slope", 4.48, 0.26, 3.86e14, 100.0, 252.7,
        (0.0, 6.9e11, 2.8e12, 6.9e12, 1.4e13, 2.1e13), (10.0, 40.0), 650.0, 10.0, 1e6,
        source=RAD,
    ),
]  # fmt: skip
lit = literature_kinetic()
k = {
    c: get_hardness_factor(252.7, hardness_table("proton", c)) for c in NIEL_CONVENTIONS
}
k62 = {
    c: get_hardness_factor(62.0, hardness_table("proton", c)) for c in NIEL_CONVENTIONS
}
res = {}
for c in NIEL_CONVENTIONS:
    obs = [dataclasses.replace(o, niel_convention=c) for o in BASE]
    g_lit = [removal_observable(o, lit) for o in obs]
    cal = calibrate_removal(obs, ("donor_removal",), base_model=lit)
    d, sd = float(cal.scales[0]), float(cal.sigma[0])
    res[c] = (d, sd, cal.chi2, np.array(cal.settings["model_values"]), np.array(g_lit))
    print(
        f"{c:10s} kappa(252.7) {k[c]:.4f}  d = {d:.3f} +- {sd:.3f} /neq  chi2 {cal.chi2:.3f}"
        f" (dof {cal.dof}); model g/p {np.round(res[c][3], 3)}; Burin-level share g/p "
        f"{np.round(g_lit, 3)}; lin.err {cal.settings['linearity_error']:.2g}"
        f"  [{time.time() - T0:.0f} s]"
    )
d_si = res["si-neq"][0]
print("\nshipped si-neq value:", CALIBRATED_SLOW_REMOVAL_BY_CONVENTION["si-neq"])
for c in NIEL_CONVENTIONS:
    d = res[c][0]
    naive = d_si * k["si-neq"] / k[c]
    rate62 = k62[c] * (lit.neutral_removal_rate + d)
    phi = DeviceConfig().N_D_bulk / rate62
    print(
        f"{c:10s} naive d_si*k_si/k_c = {naive:.3f} (vs fitted {d:.3f}); max |dg/p| vs "
        f"si-neq {np.max(np.abs(res[c][3] - res['si-neq'][3])):.3f}; 62 MeV: kappa "
        f"{k62[c]:.4f}, neutral-bulk removal {rate62:.2f} /p, Phi_comp(N_D_bulk "
        f"{DeviceConfig().N_D_bulk:.3g}) = {phi:.3e} p/cm2"
    )
print(f"\ntotal {time.time() - T0:.0f} s")
