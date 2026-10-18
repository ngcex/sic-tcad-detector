"""Tests for plotting utilities.

Smoke tests for all 5 public functions in src/plotting.py.
Uses synthetic data only -- no devsim dependency needed.
"""

import os

import matplotlib

matplotlib.use("Agg")  # Non-interactive backend, must be before pyplot import

import matplotlib.pyplot as plt
import numpy as np

from etna.core.plotting import (
    plot_electric_field,
    plot_depletion_width_vs_bias,
    plot_doping_profile,
    save_figure,
)


class TestPlotElectricField:
    """Test plot_electric_field function."""

    def teardown_method(self):
        plt.close("all")

    def test_returns_axes(self):
        x = np.linspace(0, 10e-4, 100)
        E = np.linspace(5e4, 0, 100)
        ax = plot_electric_field(x, E)
        assert isinstance(ax, matplotlib.axes.Axes)

    def test_plots_data(self):
        x = np.linspace(0, 10e-4, 100)
        E = np.linspace(5e4, 0, 100)
        ax = plot_electric_field(x, E, voltage_label="0V")
        assert len(ax.get_lines()) > 0

    def test_uses_provided_ax(self):
        fig, ax = plt.subplots()
        x = np.linspace(0, 10e-4, 100)
        E = np.linspace(5e4, 0, 100)
        returned_ax = plot_electric_field(x, E, ax=ax)
        assert returned_ax is ax

    def test_sets_labels(self):
        x = np.linspace(0, 10e-4, 100)
        E = np.linspace(5e4, 0, 100)
        ax = plot_electric_field(x, E)
        assert ax.get_xlabel() != ""
        assert ax.get_ylabel() != ""


class TestPlotDepletionWidthVsBias:
    """Test plot_depletion_width_vs_bias function."""

    def teardown_method(self):
        plt.close("all")

    def test_returns_axes(self):
        V = np.linspace(0, -30, 10)
        W = np.linspace(1.7e-4, 9e-4, 10)
        ax = plot_depletion_width_vs_bias(V, W)
        assert isinstance(ax, matplotlib.axes.Axes)

    def test_with_analytical(self):
        V = np.linspace(0, -30, 10)
        W_num = np.linspace(1.7e-4, 9e-4, 10)
        W_ana = np.linspace(1.8e-4, 8.5e-4, 10)
        ax_without = plot_depletion_width_vs_bias(V, W_num)
        n_lines_without = len(ax_without.get_lines())
        plt.close("all")
        ax_with = plot_depletion_width_vs_bias(V, W_num, W_analytical=W_ana)
        n_lines_with = len(ax_with.get_lines())
        assert n_lines_with > n_lines_without

    def test_default_experimental_data(self):
        """Without W_experimental, default reference data should be plotted."""
        V = np.linspace(0, -30, 10)
        W = np.linspace(1.7e-4, 9e-4, 10)
        ax = plot_depletion_width_vs_bias(V, W)
        # Should have at least 2 lines: numerical + experimental points
        assert len(ax.get_lines()) >= 2


class TestPlotDopingProfile:
    """Test plot_doping_profile function."""

    def teardown_method(self):
        plt.close("all")

    def test_returns_axes(self):
        x = np.linspace(0, 10e-4, 200)
        # Mixed: first half p-type, second half n-type
        nd = np.concatenate([np.full(100, -1e18), np.full(100, 1e15)])
        ax = plot_doping_profile(x, nd)
        assert isinstance(ax, matplotlib.axes.Axes)

    def test_all_positive_doping(self):
        x = np.linspace(0, 10e-4, 100)
        nd = np.full(100, 1e15)
        ax = plot_doping_profile(x, nd)
        assert isinstance(ax, matplotlib.axes.Axes)

    def test_all_negative_doping(self):
        x = np.linspace(0, 10e-4, 100)
        nd = np.full(100, -1e18)
        ax = plot_doping_profile(x, nd)
        assert isinstance(ax, matplotlib.axes.Axes)


class TestSaveFigure:
    """Test save_figure function."""

    def teardown_method(self):
        plt.close("all")

    def test_creates_png_and_pdf_in_cwd_figures(self, tmp_path, monkeypatch):
        """Default: ./figures relative to the CWD (not the package directory)."""
        import etna.core.plotting as plotting_mod

        monkeypatch.delenv("ETNA_FIGURE_DIR", raising=False)
        monkeypatch.chdir(tmp_path)
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 4, 9])
        paths = save_figure(fig, "test_output")

        fig_dir = tmp_path / "figures"
        assert (fig_dir / "test_output.png").exists()
        assert (fig_dir / "test_output.pdf").exists()
        assert paths == [
            str(fig_dir / "test_output.png"),
            str(fig_dir / "test_output.pdf"),
        ]
        pkg_dir = os.path.dirname(os.path.dirname(plotting_mod.__file__))
        assert not os.path.exists(os.path.join(pkg_dir, "figures", "test_output.png"))

    def test_env_var_overrides_the_directory(self, tmp_path, monkeypatch):
        """ETNA_FIGURE_DIR, read at call time, wins over ./figures."""
        target = tmp_path / "elsewhere"
        monkeypatch.chdir(tmp_path)
        monkeypatch.setenv("ETNA_FIGURE_DIR", str(target))
        fig, _ = plt.subplots()
        save_figure(fig, "env_output")
        assert (target / "env_output.png").exists()
        assert not (tmp_path / "figures").exists()
