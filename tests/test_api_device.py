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
