"""2D device convergence at reverse bias (split out of test_device2d.py so that
``pytest -n auto --dist loadfile`` can run it on its own worker)."""

import pytest

devsim = pytest.importorskip("devsim")

_device_counter = 0


def _unique_name():
    global _device_counter
    _device_counter += 1
    return f"test_2d_conv_{_device_counter}"


class TestReverseBiasConvergence:
    """Phase 26 / CONS-01 SC#1: 2D device converges at -15, -30, -50 V on both SV sizes.

    Plan 03 replaces xfail with real assertions after calibration completes.
    """

    @pytest.mark.slow
    def test_converges_at_minus_15V_100um(self):
        from etna.core.charge_collection_2d import create_2d_dd_device
        from etna.core.devsim_reset import reset_devsim_fully

        try:
            dev = create_2d_dd_device(
                device_name=_unique_name(), half_width_um=50.0, V_bias=15.0
            )
            assert dev["dd_initialized"] is True
            assert dev["dimension"] == 2
        finally:
            reset_devsim_fully()

    @pytest.mark.slow
    def test_converges_at_minus_30V_100um(self):
        from etna.core.charge_collection_2d import create_2d_dd_device
        from etna.core.devsim_reset import reset_devsim_fully

        try:
            dev = create_2d_dd_device(
                device_name=_unique_name(), half_width_um=50.0, V_bias=30.0
            )
            assert dev["dd_initialized"] is True
            assert dev["dimension"] == 2
        finally:
            reset_devsim_fully()

    @pytest.mark.slow
    def test_converges_at_minus_50V_100um(self):
        from etna.core.charge_collection_2d import create_2d_dd_device
        from etna.core.devsim_reset import reset_devsim_fully

        try:
            dev = create_2d_dd_device(
                device_name=_unique_name(), half_width_um=50.0, V_bias=50.0
            )
            assert dev["dd_initialized"] is True
            assert dev["dimension"] == 2
        finally:
            reset_devsim_fully()

    @pytest.mark.slow
    @pytest.mark.exhaustive  # 300 um variants (~6 s each); the build is seeded (engine start), so the -50 V build no longer ramps through this bias; run with -m exhaustive
    def test_converges_at_minus_15V_300um(self):
        from etna.core.charge_collection_2d import create_2d_dd_device
        from etna.core.devsim_reset import reset_devsim_fully

        try:
            dev = create_2d_dd_device(
                device_name=_unique_name(), half_width_um=150.0, V_bias=15.0
            )
            assert dev["dd_initialized"] is True
            assert dev["dimension"] == 2
        finally:
            reset_devsim_fully()

    @pytest.mark.slow
    @pytest.mark.exhaustive  # 300 um variants (~6 s each); the build is seeded (engine start), so the -50 V build no longer ramps through this bias; run with -m exhaustive
    def test_converges_at_minus_30V_300um(self):
        from etna.core.charge_collection_2d import create_2d_dd_device
        from etna.core.devsim_reset import reset_devsim_fully

        try:
            dev = create_2d_dd_device(
                device_name=_unique_name(), half_width_um=150.0, V_bias=30.0
            )
            assert dev["dd_initialized"] is True
            assert dev["dimension"] == 2
        finally:
            reset_devsim_fully()

    @pytest.mark.slow
    def test_converges_at_minus_50V_300um(self):
        from etna.core.charge_collection_2d import create_2d_dd_device
        from etna.core.devsim_reset import reset_devsim_fully

        try:
            dev = create_2d_dd_device(
                device_name=_unique_name(), half_width_um=150.0, V_bias=50.0
            )
            assert dev["dd_initialized"] is True
            assert dev["dimension"] == 2
        finally:
            reset_devsim_fully()
