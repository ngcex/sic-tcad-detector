"""Dark-current ramp step: the 2 V default reproduces the 0.5 V reference (fixed-point TAT state)."""

import devsim
import numpy as np
import pytest

from etna.api.device import DeviceConfig, device_kwargs
from etna.core import dark_current as dc
from etna.core.devsim_reset import reset_devsim_fully

GRID = np.array([-10.0, -30.0, -60.0])


def _sweep(T, V_step):
    reset_devsim_fully(preserve_solver=True)
    cfg = DeviceConfig(T=T)
    info = dc.create_dark_current_device(N_t=cfg.N_t, E_t=cfg.E_t, **device_kwargs(cfg))
    try:
        res = dc.dark_current_sweep(info, GRID, V_step=V_step)
        # the state is self-consistent: one more refresh + re-solve does not move it
        dc._refresh_and_resolve(info["device_name"], info["region_name"])
        again = dc.extract_dark_current_components(info)["I_total"]
    finally:
        devsim.delete_device(device=info["device_name"])
    return res, again


def test_default_step_is_2V():
    assert dc.DARK_RAMP_STEP_V == 2.0


@pytest.mark.parametrize(
    "T", [300.0, 400.0]
)  # 400 K: largest lag of the single refresh
def test_2V_ramp_matches_half_volt_reference(T):
    ref, _ = _sweep(T, 0.5)
    res, again = _sweep(T, dc.DARK_RAMP_STEP_V)
    assert not res["truncated"] and not ref["truncated"]
    np.testing.assert_allclose(res["I_total"], ref["I_total"], rtol=1e-12, atol=0)
    assert again == pytest.approx(res["I_total"][-1], rel=1e-9)


def test_default_area_is_the_4mm2_dosimeter():
    assert dc.AREA_DOSIMETER_CM2 == 0.04
    import inspect

    for fn in (
        dc.extract_dark_current_components,
        dc.dark_current_sweep,
        dc.dark_current_at_bias,
        dc.sensitivity_sweep,
        dc.dark_current_vs_fluence,
        dc.dark_current_post_anneal,
    ):
        assert inspect.signature(fn).parameters["area"].default == dc.AREA_DOSIMETER_CM2
