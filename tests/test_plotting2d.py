"""Tests for 2D Poisson solve, 1D-vs-2D validation, and 2D visualization.

Tests cover:
- 2D Poisson equation convergence at equilibrium and under reverse bias
- Quantitative 1D-vs-2D validation (MESH-02): potential and E-field within 1%
- Visualization functions return valid matplotlib figures (MESH-03)

All tests use unique device names to avoid devsim global state conflicts.
"""

import matplotlib

matplotlib.use("Agg")  # headless backend for CI/testing

import matplotlib.pyplot as plt
import numpy as np
import pytest

devsim = pytest.importorskip("devsim")

# Unique device name counter to avoid devsim name collisions
_device_counter = 0


def _unique_name():
    global _device_counter
    _device_counter += 1
    return f"test_p2d_dev_{_device_counter}"


class TestPoisson2DSolve:
    """Test that 2D Poisson equation converges on device2d mesh."""

    def test_equilibrium_converges_2d(self):
        """2D Poisson solve should converge at 0V equilibrium."""
        from src.device2d import create_sic_2d_device
        from src.poisson import setup_poisson, solve_equilibrium

        dev = create_sic_2d_device(
            device_name=_unique_name(), half_width_um=50, doping_profile="graded"
        )
        setup_poisson(dev)
        solve_equilibrium(dev)  # should not raise

        # Verify Potential exists and is non-trivial
        pot = np.array(
            devsim.get_node_model_values(
                device=dev["device_name"],
                region=dev["region_name"],
                name="Potential",
            )
        )
        assert len(pot) > 0
        assert np.any(
            pot != 0.0
        ), "Potential should be non-zero after equilibrium solve"

    def test_potential_range_physical(self):
        """Equilibrium potential should be in a physically reasonable range."""
        from src.device2d import create_sic_2d_device
        from src.poisson import setup_poisson, solve_equilibrium

        dev = create_sic_2d_device(
            device_name=_unique_name(), half_width_um=50, doping_profile="graded"
        )
        setup_poisson(dev)
        solve_equilibrium(dev)

        pot = np.array(
            devsim.get_node_model_values(
                device=dev["device_name"],
                region=dev["region_name"],
                name="Potential",
            )
        )
        assert (
            pot.min() > -5.0
        ), f"Potential min {pot.min():.2f}V too low for equilibrium"
        assert (
            pot.max() < 5.0
        ), f"Potential max {pot.max():.2f}V too high for equilibrium"

    def test_reverse_bias_solves(self):
        """2D Poisson should converge under +5V cathode reverse bias."""
        import devsim.python_packages.simple_physics as simple_physics
        from src.device2d import create_sic_2d_device
        from src.poisson import setup_poisson, solve_equilibrium

        dev = create_sic_2d_device(
            device_name=_unique_name(), half_width_um=50, doping_profile="graded"
        )
        setup_poisson(dev)
        solve_equilibrium(dev)

        # Ramp cathode to +5V in small steps (reverse bias for p+/n-)
        bias_name = simple_physics.GetContactBiasName("cathode")
        for v in [1.0, 2.0, 3.0, 4.0, 5.0]:
            devsim.set_parameter(device=dev["device_name"], name=bias_name, value=v)
            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=1e-10,
                maximum_iterations=40,
            )

        # Verify potential range expanded under reverse bias
        pot = np.array(
            devsim.get_node_model_values(
                device=dev["device_name"],
                region=dev["region_name"],
                name="Potential",
            )
        )
        pot_range = pot.max() - pot.min()
        assert (
            pot_range > 3.0
        ), f"Potential range {pot_range:.2f}V too small under 5V reverse bias"


class TestValidation2Dvs1D:
    """Quantitative 1D-vs-2D validation (MESH-02 requirement)."""

    @pytest.fixture(scope="class")
    def solved_devices(self):
        """Create and solve both 1D and 2D devices with matching parameters."""
        from src.device import create_sic_device
        from src.device2d import create_sic_2d_device
        from src.poisson import setup_poisson, solve_equilibrium

        # Matching parameters for fair comparison
        params = dict(
            epi_thickness_cm=10e-4,
            substrate_thickness_cm=1e-4,
            N_A=1e19,
            T=300,
            doping_profile="graded",
        )

        dev_1d = create_sic_device(device_name=_unique_name(), **params)
        setup_poisson(dev_1d)
        solve_equilibrium(dev_1d)

        dev_2d = create_sic_2d_device(
            device_name=_unique_name(), half_width_um=50, **params
        )
        setup_poisson(dev_2d)
        solve_equilibrium(dev_2d)

        return dev_1d, dev_2d

    def test_potential_matches_1d_within_1pct(self, solved_devices):
        """2D center-column potential should match 1D within 1%."""
        from src.plotting2d import validate_2d_vs_1d

        dev_1d, dev_2d = solved_devices
        result = validate_2d_vs_1d(dev_2d, dev_1d)
        assert (
            result["potential_max_rel_error"] < 0.01
        ), f"Potential error {result['potential_max_rel_error']:.4e} exceeds 1%"

    def test_efield_matches_1d_within_1pct(self, solved_devices):
        """2D center-column E-field should match 1D within 1%."""
        from src.plotting2d import validate_2d_vs_1d

        dev_1d, dev_2d = solved_devices
        result = validate_2d_vs_1d(dev_2d, dev_1d)
        assert (
            result["efield_max_rel_error"] < 0.01
        ), f"E-field error {result['efield_max_rel_error']:.4e} exceeds 1%"


class TestVisualization2D:
    """Test that 2D visualization functions produce valid figures."""

    @pytest.fixture(scope="class")
    def solved_2d_device(self):
        """Create and solve a 2D device for visualization tests."""
        from src.device2d import create_sic_2d_device
        from src.poisson import setup_poisson, solve_equilibrium

        dev = create_sic_2d_device(
            device_name=_unique_name(), half_width_um=50, doping_profile="graded"
        )
        setup_poisson(dev)
        solve_equilibrium(dev)
        return dev

    def test_plot_potential_returns_figure(self, solved_2d_device):
        """plot_potential_2d should return (fig, ax)."""
        from src.plotting2d import plot_potential_2d

        dev = solved_2d_device
        fig, ax = plot_potential_2d(dev["device_name"], dev["region_name"])
        assert isinstance(fig, plt.Figure)
        assert ax is not None
        plt.close(fig)

    def test_plot_efield_returns_figure(self, solved_2d_device):
        """plot_efield_2d should return (fig, ax)."""
        from src.plotting2d import plot_efield_2d

        dev = solved_2d_device
        fig, ax = plot_efield_2d(dev["device_name"], dev["region_name"])
        assert isinstance(fig, plt.Figure)
        assert ax is not None
        plt.close(fig)

    def test_triangulation_has_correct_shape(self, solved_2d_device):
        """Triangulation should have correct node count and triangle shape."""
        from src.plotting2d import get_triangulation

        dev = solved_2d_device
        tri = get_triangulation(dev["device_name"], dev["region_name"])

        # Number of triangulation nodes should match device node count
        assert tri.x.shape[0] == dev["num_nodes"]
        # Triangles should have 3 vertices each
        assert tri.triangles.shape[1] == 3
        # Should have non-trivial number of triangles
        assert tri.triangles.shape[0] > 100
