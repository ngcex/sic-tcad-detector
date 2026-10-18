"""Tests for etna/core/devsim_reset.py — the single-responsibility full reset utility.

Covers PITFALLS P03/P30 (cylindrical-axis globals leak from
etna/core/alternative_structures.py) and P20 (enumerate devices, never hardcode).

The pure-constant tests (import + _CYLINDRICAL_GLOBALS contents) run without
devsim. The behavioural tests (device deletion, cylindrical-state restore,
solver preservation) require a real devsim session.
"""

import pytest


@pytest.fixture(autouse=True)
def _restore_solver_state():
    """Put devsim's solver parameters back after each test.

    These tests deliberately set a bogus ``direct_solver`` and call resets that
    can leave it "unknown"; without this restore the next ``devsim.solve`` in
    the same process aborts the interpreter.
    """
    devsim = pytest.importorskip("devsim")
    saved = {}
    for key in ("direct_solver", "solver_callback"):
        try:
            saved[key] = devsim.get_parameter(name=key)
        except Exception:
            pass
    yield
    for key, value in saved.items():
        devsim.set_parameter(name=key, value=value)


def test_import_reset_devsim_fully():
    """Test 1: public API imports succeed."""
    from etna.core.devsim_reset import reset_devsim_fully, _CYLINDRICAL_GLOBALS

    assert callable(reset_devsim_fully)
    assert isinstance(_CYLINDRICAL_GLOBALS, tuple)


def test_cylindrical_globals_exact_contents():
    """Test 2: _CYLINDRICAL_GLOBALS is exactly the seven leaking globals."""
    from etna.core.devsim_reset import _CYLINDRICAL_GLOBALS

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
    from etna.core.devsim_reset import reset_devsim_fully

    # Should not raise even with no devices / no cylindrical state.
    reset_devsim_fully()
    reset_devsim_fully(preserve_solver=False)


class TestResetClearsState:
    """Behavioural tests requiring a real devsim session."""

    def test_restores_cartesian_node_volume_model(self):
        """Test 3: after a cylindrical session, node_volume_model is Cartesian again."""
        devsim = pytest.importorskip("devsim")
        from etna.core.devsim_reset import reset_devsim_fully

        # Simulate the leak: set the cylindrical assembly-model global.
        devsim.set_parameter(name="node_volume_model", value="CylindricalNodeVolume")
        reset_devsim_fully()
        assert (
            devsim.get_parameter(name="node_volume_model") == "NodeVolume"
        ), "reset_devsim_fully must restore the Cartesian node_volume_model"

    def test_deletes_all_devices(self):
        """Test 4: after reset, get_device_list() is empty."""
        devsim = pytest.importorskip("devsim")
        from etna.core.devsim_reset import reset_devsim_fully
        from etna.core.device2d import create_sic_2d_device

        create_sic_2d_device(device_name="reset_canary_dev", half_width_um=50)
        assert len(devsim.get_device_list()) >= 1
        reset_devsim_fully()
        assert list(devsim.get_device_list()) == []

    def test_preserves_direct_solver(self):
        """Test 5: with preserve_solver=True the direct_solver is preserved."""
        devsim = pytest.importorskip("devsim")
        from etna.core.devsim_reset import reset_devsim_fully

        devsim.set_parameter(name="direct_solver", value="mkl_pardiso")
        reset_devsim_fully(preserve_solver=True)
        assert devsim.get_parameter(name="direct_solver") == "mkl_pardiso"


def test_reset_undoes_every_cylindrical_global_without_the_removed_step():
    """devsim.reset_devsim() alone (the former "restore Cartesian defaults" step was
    dead and is removed) restores the five assembly models to devsim's own Cartesian
    values == _CARTESIAN_DEFAULTS and removes raxis_zero / raxis_variable."""
    devsim = pytest.importorskip("devsim")
    from etna.core.devsim_reset import (
        _CARTESIAN_DEFAULTS,
        _CYLINDRICAL_GLOBALS,
        reset_devsim_fully,
    )

    # exactly what alternative_structures._activate_cylindrical_coords sets
    cyl = {
        "raxis_zero": 0.0,
        "raxis_variable": "x",
        "node_volume_model": "CylindricalNodeVolume",
        "edge_couple_model": "CylindricalEdgeCouple",
        "element_edge_couple_model": "ElementCylindricalEdgeCouple",
        "element_node0_volume_model": "ElementCylindricalNodeVolume@en0",
        "element_node1_volume_model": "ElementCylindricalNodeVolume@en1",
    }
    assert set(cyl) == set(_CYLINDRICAL_GLOBALS)
    for name, value in cyl.items():
        devsim.set_parameter(name=name, value=value)
    reset_devsim_fully()
    for name, value in _CARTESIAN_DEFAULTS.items():
        assert devsim.get_parameter(name=name) == value
    for name in ("raxis_zero", "raxis_variable"):
        with pytest.raises(devsim.error):
            devsim.get_parameter(name=name)


def test_solve_still_works_after_reset_tests():
    """The solver parameters are intact, so a real equilibrium solve converges."""
    pytest.importorskip("devsim")
    from etna.core.drift_diffusion import create_dd_device

    info = create_dd_device(device_name="post_reset_canary")
    assert info["dd_initialized"]
