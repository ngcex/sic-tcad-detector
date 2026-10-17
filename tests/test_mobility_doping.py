"""Both mobilities follow the config's epi doping (graded: N_D_bulk) in every builder.

Before 2026-10-17 graded devices evaluated mu_n / mu_p at DEFAULT_N_D = 1.07e15 whatever
N_D_bulk was, so a sidebar change of N_D_bulk was silently ignored by the transport
(config-field-ignored bug class). One rule now: etna.core.device.mobility_doping.
"""

from dataclasses import replace

import numpy as np
import pytest

devsim = pytest.importorskip("devsim")

from etna.api import simulation as sim  # noqa: E402
from etna.api.device import DeviceConfig, build_device, device_kwargs  # noqa: E402
from etna.core.charge_collection import cce_vs_bias  # noqa: E402
from etna.core.device import DEFAULT_N_D, DEFAULT_STACK, mobility_doping  # noqa: E402
from etna.core.sic_material import mobility_caughey_thomas_T  # noqa: E402

N_D_BULK_ALT = 5e13


def test_mobility_doping_rule():
    assert mobility_doping(None, "graded", None) == DEFAULT_STACK["N_D_bulk"]
    assert mobility_doping(3e15, "graded", 4e13) == 4e13
    assert mobility_doping(None, "uniform", 4e13) == DEFAULT_N_D
    assert mobility_doping(2e15, "uniform", 4e13) == 2e15


@pytest.mark.parametrize("half_width_um", [None, 20.0])
def test_devsim_builders_follow_N_D_bulk(half_width_um):
    cfg = DeviceConfig(N_D_bulk=N_D_BULK_ALT, half_width_um=half_width_um)
    info = build_device(cfg)
    try:
        assert info["N_D_mobility"] == N_D_BULK_ALT
        for carrier, key in (("electron", "mu_n"), ("hole", "mu_p")):
            mu = mobility_caughey_thomas_T(N_D_BULK_ALT, cfg.T, carrier)
            assert info[key] == pytest.approx(mu, rel=1e-12)
            got = devsim.get_parameter(
                device=info["device_name"], region=info["region_name"], name=key
            )
            assert got == pytest.approx(mu, rel=1e-12)
    finally:
        devsim.delete_device(device=info["device_name"])


def test_engine_batch_follows_N_D_bulk():
    cfgs = [DeviceConfig(), DeviceConfig(N_D_bulk=N_D_BULK_ALT)]
    _, _, nd_mob = sim._engine_batch([device_kwargs(c) for c in cfgs])
    np.testing.assert_array_equal(nd_mob, [DEFAULT_STACK["N_D_bulk"], N_D_BULK_ALT])


def test_engine_devsim_cce_parity_non_default_N_D_bulk():
    """Same mobility in both engines: CCE parity <= 1e-10 at a non-default N_D_bulk."""
    cfg = replace(DeviceConfig(), N_D_bulk=N_D_BULK_ALT)
    eng = sim.run_cce(cfg, v_start=-10.0, v_stop=-30.0, n_points=3)
    assert eng.metadata["solver"] == "batched"
    kw = device_kwargs(cfg)
    ref = cce_vs_bias(
        np.asarray(eng.x), epi_thickness_cm=kw.pop("epi_thickness_cm"), device_kwargs=kw
    )
    np.testing.assert_allclose(eng.y, ref["cce_values"], rtol=0, atol=1e-10)
