"""Tests for the public build_device() API facade.

WR-03 regression guard: build_device() is the dispatch point between the
1D and 2D construction paths (etna/api/device.py) and was previously
only exercised indirectly through run_cv/run_field. These tests assert its
DD-initialization and dimension-tagging contract directly for both
branches, so a regression that broke device_name uniqueness or the
dd_initialized/dimension bookkeeping would be caught here without needing
a full run_cv/run_field solve.
"""

import pytest

from etna import DeviceConfig
from etna.api.device import build_device


@pytest.mark.slow
class TestBuildDevice:
    """Integration test: build_device() with live devsim devices."""

    def test_build_device_1d(self):
        device_info = build_device(DeviceConfig())

        assert device_info["dd_initialized"] is True
        # 1D branch (create_dd_device) does not set a "dimension" key.
        assert device_info.get("dimension") is None

    def test_build_device_2d(self):
        device_info = build_device(DeviceConfig(half_width_um=50.0))

        assert device_info["dd_initialized"] is True
        assert device_info["dimension"] == 2

    def test_build_device_generates_unique_device_names(self):
        info_a = build_device(DeviceConfig())
        info_b = build_device(DeviceConfig())

        assert info_a["device_name"] != info_b["device_name"]


def test_build_device_forwards_explicit_device_name(monkeypatch):
    """device_name= reaches the 1D and 2D core builders (audit C-13); no devsim."""
    import etna.api.device as api_device

    seen = []

    def fake_builder(device_name, **kwargs):
        seen.append(device_name)
        return {"device_name": device_name}

    monkeypatch.setattr(api_device, "create_dd_device", fake_builder)
    monkeypatch.setattr(api_device, "create_sic_2d_device", fake_builder)
    for fn in ("setup_poisson", "solve_equilibrium", "setup_sic_drift_diffusion"):
        monkeypatch.setattr(api_device, fn, lambda info: None)

    assert build_device(DeviceConfig(), device_name="my_1d")["device_name"] == "my_1d"
    info_2d = build_device(DeviceConfig(half_width_um=25.0), device_name="my_2d")
    assert info_2d["device_name"] == "my_2d"
    assert seen == ["my_1d", "my_2d"]
    # without a name, a unique generated one is used
    assert build_device(DeviceConfig())["device_name"].startswith("device1d_")
