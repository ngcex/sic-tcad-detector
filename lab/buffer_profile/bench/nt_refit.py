"""Refit N_t of the measured-C-V preset to I(300 K, -60 V, 4 mm2) = 18 pA (leakage is
linear in N_t: one solve at the default N_t, rescale, one verification solve).

usage: uv run python lab/buffer_profile/bench/nt_refit.py
"""

from __future__ import annotations

import uuid

import devsim
import numpy as np

from etna.api.device import device_kwargs, measured_cv_config
from etna.core.dark_current import (
    create_dark_current_device,
    dark_current_sweep,
    setup_tat_model,
)
from etna.core.device import DEFAULT_STACK


def I60(N_t, E_t):
    cfg = measured_cv_config()
    info = create_dark_current_device(
        N_t=N_t,
        E_t=E_t,
        device_name=f"ntr_{uuid.uuid4().hex[:6]}",
        **device_kwargs(cfg),
    )
    try:
        setup_tat_model(info, E_t=E_t, N_t=N_t)
        r = dark_current_sweep(info, np.array([-10.0, -30.0, -60.0]), area=0.04)
    finally:
        devsim.delete_device(device=info["device_name"])
    return r


def main():
    E_t = DEFAULT_STACK["E_t"]
    N0 = DEFAULT_STACK["N_t"]
    r0 = I60(N0, E_t)
    I0 = abs(r0["I_total"][-1])
    print(
        f"RES default N_t {N0:.7g}: I(-10/-30/-60) = "
        + " / ".join(f"{abs(i)*1e12:.3f}" for i in r0["I_total"])
        + " pA"
    )
    N1 = N0 * 18e-12 / I0
    r1 = I60(N1, E_t)
    I1 = abs(r1["I_total"][-1])
    N2 = N1 * 18e-12 / I1
    print(
        f"RES refit N_t {N1:.7g}: I(-60) = {I1*1e12:.5f} pA (target 18); second pass N_t {N2:.7g}"
    )
    r2 = I60(float(f"{N2:.7g}"), E_t)
    print(
        f"RES rounded N_t {float(f'{N2:.7g}'):.7g}: I(-10/-30/-60) = "
        + " / ".join(f"{abs(i)*1e12:.4f}" for i in r2["I_total"])
        + " pA; components at -60 V: SRH "
        f"{abs(r2['I_SRH'][-1])*1e12:.3f} TAT {abs(r2['I_TAT'][-1])*1e12:.3f} SRV {abs(r2['I_SRV'][-1])*1e12:.3f} pA"
    )


if __name__ == "__main__":
    main()
