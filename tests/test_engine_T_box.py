"""Engine T box 150-450 K: run_cv / run_cce served by the engine, parity with devsim."""

import numpy as np
import pytest

from etna.api import simulation as sim
from etna.api.device import DeviceConfig, legacy_config

BIAS_CV = np.linspace(0.0, -40.0, 9)
BIAS_CCE = np.linspace(-10.0, -40.0, 4)


@pytest.mark.parametrize(
    "cfg",
    [DeviceConfig(T=150.0), DeviceConfig(T=200.0), legacy_config(T=450.0)],
    ids=["default-150K", "default-200K", "legacy-450K"],
)
def test_engine_matches_devsim_outside_old_box(cfg):
    assert sim.validated_domain_reason(cfg) is None
    cv = sim.run_cv(cfg, 0.0, -40.0, 9)
    assert cv.metadata["solver"] == "batched"
    ref = sim._run_cv_devsim(cfg, BIAS_CV)
    n = len(ref.x)
    np.testing.assert_allclose(cv.y[:n], ref.y, rtol=1e-6)
    np.testing.assert_allclose(
        cv.metadata["depletion_widths"][:n],
        ref.metadata["depletion_widths"],
        rtol=0,
        atol=1e-14,
    )
    cce = sim.run_cce(cfg, -10.0, -40.0, 4)
    assert cce.metadata["solver"] == "batched"
    ref = sim._run_cce_devsim(cfg, BIAS_CCE)
    np.testing.assert_allclose(cce.y[: len(ref.x)], ref.y, rtol=0, atol=1e-10)


@pytest.mark.parametrize("T", [140.0, 460.0])
def test_outside_T_box_is_devsim(T):
    reason = sim.validated_domain_reason(DeviceConfig(T=T))
    assert reason is not None and reason.startswith("T=")
