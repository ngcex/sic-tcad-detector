"""Tests for src/devsim_reset.py — the single-responsibility full reset utility.

Covers PITFALLS P03/P30 (cylindrical-axis globals leak from
src/alternative_structures.py) and P20 (enumerate devices, never hardcode).

The pure-constant tests (import + _CYLINDRICAL_GLOBALS contents) run without
devsim. The behavioural tests (device deletion, cylindrical-state restore,
solver preservation) require a real devsim session.
"""

import pytest


def test_import_reset_devsim_fully():
    """Test 1: public API imports succeed."""
    from src.devsim_reset import reset_devsim_fully, _CYLINDRICAL_GLOBALS

    assert callable(reset_devsim_fully)
    assert isinstance(_CYLINDRICAL_GLOBALS, tuple)


def test_cylindrical_globals_exact_contents():
    """Test 2: _CYLINDRICAL_GLOBALS is exactly the seven leaking globals."""
    from src.devsim_reset import _CYLINDRICAL_GLOBALS

    expected = (
        "raxis_zero",
        "raxis_variable",
        "node_volume_model",
        "edge_couple_model",
        "element_edge_couple_model",
        "element_node0_volume_model",
        "element_node1_volume_model",
    )
    assert _CYLINDRICAL_GLOBALS == expected
    assert len(_CYLINDRICAL_GLOBALS) == 7


def test_idempotent_on_empty_state():
    """Test 6 (partial): reset is safe to call on an empty session."""
    pytest.importorskip("devsim")
    from src.devsim_reset import reset_devsim_fully

    # Should not raise even with no devices / no cylindrical state.
    reset_devsim_fully()
    reset_devsim_fully(preserve_solver=False)


class TestResetClearsState:
    """Behavioural tests requiring a real devsim session."""

    def test_restores_cartesian_node_volume_model(self):
        """Test 3: after a cylindrical session, node_volume_model is Cartesian again."""
        devsim = pytest.importorskip("devsim")
        from src.devsim_reset import reset_devsim_fully

        # Simulate the leak: set the cylindrical assembly-model global.
        devsim.set_parameter(name="node_volume_model", value="CylindricalNodeVolume")
        reset_devsim_fully()
        assert (
            devsim.get_parameter(name="node_volume_model") == "NodeVolume"
        ), "reset_devsim_fully must restore the Cartesian node_volume_model"

    def test_deletes_all_devices(self):
        """Test 4: after reset, get_device_list() is empty."""
        devsim = pytest.importorskip("devsim")
        from src.devsim_reset import reset_devsim_fully
        from src.device2d import create_sic_2d_device

        create_sic_2d_device(device_name="reset_canary_dev", half_width_um=50)
        assert len(devsim.get_device_list()) >= 1
        reset_devsim_fully()
        assert list(devsim.get_device_list()) == []

    def test_preserves_direct_solver(self):
        """Test 5: with preserve_solver=True the direct_solver is preserved."""
        devsim = pytest.importorskip("devsim")
        from src.devsim_reset import reset_devsim_fully

        devsim.set_parameter(name="direct_solver", value="mkl_pardiso")
        reset_devsim_fully(preserve_solver=True)
        assert devsim.get_parameter(name="direct_solver") == "mkl_pardiso"
