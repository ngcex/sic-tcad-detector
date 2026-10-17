"""1-parameter fit: uniform epi N_D such that the full-depletion (reach-through) voltage is
30 V, with the production definition (etna.api.simulation._reach_through_V +
_refine_reach_through: depletion edge within one Debye length of the epi end, 0.05 V lattice).

usage: uv run python lab/microdosimeter/bench/fit_vfd.py
"""

import numpy as np
import scipy.optimize

from etna import DeviceConfig, run_cv
from etna.core.batched.material import EPS, Q

TARGET = 30.0


def cfg(N_D, sub_um):
    return DeviceConfig(
        doping_profile="uniform", N_D=N_D, substrate_thickness_um=sub_um, area_cm2=0.25
    )


def vfd(N_D, sub_um):
    r = run_cv(cfg(N_D, sub_um), 0.0, -60.0, 121)
    assert r.metadata["solver"] == "batched", r.metadata
    v = r.metadata["reach_through_V"]
    return 100.0 if v is None else v


for sub in (0.3, 1.0):
    f = lambda lg: vfd(10**lg, sub) - TARGET
    lg = scipy.optimize.brentq(f, 14.0, 15.0, xtol=1e-5)
    N = 10**lg
    # analytic one-sided abrupt estimate, V_bi ~ 3 V
    an = 2 * EPS * (TARGET + 3.0) / (Q * (10e-4) ** 2)
    print(
        f"RES p+ {sub} um: N_D = {N:.4e} cm^-3 -> V_fd {vfd(N, sub):.2f} V "
        f"(analytic abrupt 2eps(V+Vbi)/qd^2 = {an:.3e});"
        f" V_fd(1e14) = {vfd(1e14, sub)}, V_fd(0.95N) = {vfd(0.95 * N, sub)}, V_fd(1.05N) = {vfd(1.05 * N, sub)}"
    )
    r = run_cv(cfg(N, sub), 0.0, -60.0, 13)
    print(
        "   C plateau (-60 V) = %.3e F (25 mm2); same device at 4 mm2 would give %.3e F"
        % (r.y[-1], r.y[-1] * 0.04 / 0.25)
    )
