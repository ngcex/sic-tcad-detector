"""cv_sweep solves exactly at each requested bias, also off the 0.5 V lattice (A-01)."""

import devsim
import numpy as np
import pytest

from etna.api.device import DeviceConfig, build_device
from etna.core.cv_analysis import cv_sweep
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.drift_diffusion import ramp_bias
from etna.core.poisson import extract_depletion_width_numerical


@pytest.mark.parametrize("V", [-1.282, -3.3, -7.77])
def test_cv_sweep_width_is_exact_off_the_lattice(V):
    cfg = DeviceConfig()
    reset_devsim_fully()
    res = cv_sweep(build_device(cfg), np.array([0.0, V]), area=cfg.area_cm2)
    assert not res["truncated"]
    reset_devsim_fully()
    info = build_device(cfg)
    ramp_bias(info, -V, contact="cathode", V_step=0.05)
    W_ref = extract_depletion_width_numerical(info)
    devsim.delete_device(device=info["device_name"])
    assert res["depletion_widths"][1] == pytest.approx(W_ref, rel=1e-6)
