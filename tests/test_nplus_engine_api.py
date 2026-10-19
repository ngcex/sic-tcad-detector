"""run_cce / run_cv route the n+ default to the batched engine and match devsim.

Split out of test_nplus_engine.py so --dist loadfile can run it on its own worker.
"""

import numpy as np

from etna.api import simulation as sim
from etna.api.device import DeviceConfig, build_device, device_kwargs
from etna.core.charge_collection import cce_vs_bias
from etna.core.cv_analysis import cv_sweep
from etna.core.devsim_reset import reset_devsim_fully


def test_api_routes_nplus_to_the_engine_and_matches_devsim():
    cfg = DeviceConfig(nplus_thickness_um=1.0)
    assert sim.validated_domain_reason(cfg) is None
    res = sim.run_cce(cfg, v_start=-10.0, v_stop=-100.0, n_points=4)
    assert res.metadata["solver"] == "batched"
    assert res.metadata["cross_checked_to_V"] == 100.0
    kw = device_kwargs(cfg)
    ref = cce_vs_bias(
        np.linspace(-10.0, -100.0, 4),
        epi_thickness_cm=kw.pop("epi_thickness_cm"),
        device_kwargs=kw,
    )
    np.testing.assert_allclose(res.y, ref["cce_values"], atol=1e-6)
    cv = sim.run_cv(cfg, v_start=0.0, v_stop=-20.0, n_points=9)
    assert cv.metadata["solver"] == "batched"
    reset_devsim_fully()
    ref_cv = cv_sweep(build_device(cfg), np.linspace(0.0, -20.0, 9), area=cfg.area_cm2)
    reset_devsim_fully()
    np.testing.assert_allclose(
        cv.metadata["depletion_widths"], ref_cv["depletion_widths"], rtol=1e-10
    )
    # reach-through with n+ and the graded transition (default since 2026-10-18): edge
    # within one Debye length of the transition start, 3.1 V (12.3 V with the abrupt n+)
    assert 2.0 < cv.metadata["reach_through_V"] < 5.0
