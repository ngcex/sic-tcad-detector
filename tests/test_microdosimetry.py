"""Tests for microdosimetry module.

Validates:
- Mean chord length (slab, 3D rectangular parallelepiped, cube)
- Log-spaced y-bins (count, range, geometric mean centers)
- Lineal energy spectrum f(y) and d(y) normalization
- y_D >= y_F for all distributions (Jensen's inequality)
- Monoenergetic vs broad distribution behavior
- Kappa table computation from stopping power CSVs
- Tissue-equivalence correction (constant and energy-dependent)
- Publication-quality plot functions (Agg backend, no exceptions)
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from src.microdosimetry import (
    compute_kappa_table,
    compute_microdosimetric_means,
    lineal_energy_spectrum,
    make_y_bins,
    mean_chord_length,
    plot_yd_spectrum,
    plot_yf_spectrum,
    tissue_equivalence_correction,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def broad_energies():
    """Broad lognormal energy distribution (2000 events)."""
    rng = np.random.default_rng(42)
    return rng.lognormal(mean=2.0, sigma=1.0, size=2000)


@pytest.fixture
def monoenergetic_energies():
    """Monoenergetic distribution (1000 events, all 10 keV)."""
    return np.full(1000, 10.0)


@pytest.fixture
def spectrum_broad(broad_energies):
    """Pre-computed spectrum from broad distribution."""
    l_bar = 20.0  # 10 um slab
    return lineal_energy_spectrum(broad_energies, l_bar)


@pytest.fixture
def spectrum_mono(monoenergetic_energies):
    """Pre-computed spectrum from monoenergetic distribution."""
    l_bar = 20.0
    return lineal_energy_spectrum(monoenergetic_energies, l_bar)


@pytest.fixture
def kappa_table():
    """Pre-computed kappa table from bundled stopping power data."""
    return compute_kappa_table()


# ---------------------------------------------------------------------------
# Mean chord length tests
# ---------------------------------------------------------------------------


class TestMeanChordLength:
    def test_slab(self):
        """l_bar = 2*thickness for slab approximation."""
        assert mean_chord_length(10.0) == pytest.approx(20.0)

    def test_slab_different_thickness(self):
        """Slab approximation scales linearly with thickness."""
        assert mean_chord_length(5.0) == pytest.approx(10.0)

    def test_3d_rectangular(self):
        """l_bar = 4V/S for 100x100x10 um rectangular parallelepiped."""
        # V = 100*100*10 = 100000
        # S = 2*(100*100 + 100*10 + 100*10) = 2*(10000+1000+1000) = 24000
        # l_bar = 4*100000/24000 = 16.667
        l_bar = mean_chord_length(10.0, sv_width_um=100.0, sv_depth_um=100.0)
        expected = 4.0 * 100000.0 / 24000.0
        assert l_bar == pytest.approx(expected, rel=1e-6)

    def test_cube(self):
        """l_bar = 2/3 * side for a cube (known analytical result)."""
        side = 50.0
        l_bar = mean_chord_length(side, sv_width_um=side, sv_depth_um=side)
        # 4V/S = 4*side^3 / (6*side^2) = 2/3 * side
        expected = 2.0 / 3.0 * side
        assert l_bar == pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------------
# Y-bins tests
# ---------------------------------------------------------------------------


class TestMakeYBins:
    def test_bin_count_default(self):
        """Default params: 6 decades * 50/decade = 300 bins."""
        bins = make_y_bins()
        assert len(bins["bin_centers"]) == 300
        assert len(bins["bin_widths"]) == 300
        assert len(bins["bin_edges"]) == 301

    def test_bin_range(self):
        """Edges span from y_min to y_max."""
        bins = make_y_bins(y_min=0.01, y_max=1e4)
        assert bins["bin_edges"][0] == pytest.approx(0.01, rel=1e-6)
        assert bins["bin_edges"][-1] == pytest.approx(1e4, rel=1e-6)

    def test_geometric_mean_centers(self):
        """Bin centers are geometric means of adjacent edges."""
        bins = make_y_bins()
        expected = np.sqrt(bins["bin_edges"][:-1] * bins["bin_edges"][1:])
        np.testing.assert_allclose(bins["bin_centers"], expected, rtol=1e-10)

    def test_custom_bins_per_decade(self):
        """Custom bins_per_decade produces correct count."""
        bins = make_y_bins(y_min=1.0, y_max=100.0, bins_per_decade=20)
        # 2 decades * 20 = 40 bins
        assert len(bins["bin_centers"]) == 40


# ---------------------------------------------------------------------------
# Lineal energy spectrum tests
# ---------------------------------------------------------------------------


class TestLinealEnergySpectrum:
    def test_f_y_normalization(self, spectrum_broad):
        """Integral f(y)*dy should be approximately 1.0 for large sample."""
        integral = np.sum(spectrum_broad["f_y"] * spectrum_broad["bin_widths"])
        assert integral == pytest.approx(1.0, abs=0.05)

    def test_y_values_computation(self):
        """y = E / l_bar for known inputs."""
        energies = np.array([20.0, 40.0, 60.0])
        l_bar = 10.0
        result = lineal_energy_spectrum(energies, l_bar)
        np.testing.assert_allclose(result["y_values"], [2.0, 4.0, 6.0])

    def test_y_D_ge_y_F(self, spectrum_broad):
        """y_D >= y_F (Jensen's inequality) for broad distribution."""
        assert spectrum_broad["y_D"] >= spectrum_broad["y_F"]

    def test_monoenergetic_y_F_approx_y_D(self, spectrum_mono):
        """For delta-like distribution, y_F should approximately equal y_D."""
        # Monoenergetic: all events in same bin, so y_F ~ y_D
        # Allow some tolerance due to bin discretization
        ratio = spectrum_mono["y_D"] / spectrum_mono["y_F"]
        assert ratio == pytest.approx(1.0, abs=0.15)

    def test_broad_y_D_gt_y_F(self, spectrum_broad):
        """For broad distribution, y_D should be strictly greater than y_F."""
        assert spectrum_broad["y_D"] > spectrum_broad["y_F"] * 1.1

    def test_d_y_normalization(self, spectrum_broad):
        """Integral d(y)*dy should be approximately 1.0."""
        integral = np.sum(spectrum_broad["d_y"] * spectrum_broad["bin_widths"])
        assert integral == pytest.approx(1.0, abs=0.05)

    def test_n_events(self, spectrum_broad):
        """n_events should match input array length."""
        assert spectrum_broad["n_events"] == 2000

    def test_empty_bins_handled(self):
        """Spectrum with events in narrow range: empty bins are zero."""
        energies = np.full(100, 5.0)
        result = lineal_energy_spectrum(energies, 10.0)
        # Most bins should be zero
        n_nonzero = np.count_nonzero(result["f_y"])
        assert n_nonzero < 10  # only a few bins hit


# ---------------------------------------------------------------------------
# Compute microdosimetric means (standalone)
# ---------------------------------------------------------------------------


class TestComputeMicrodosimetricMeans:
    def test_agrees_with_spectrum(self, spectrum_broad):
        """Standalone means should match those from lineal_energy_spectrum."""
        means = compute_microdosimetric_means(
            spectrum_broad["bin_centers"],
            spectrum_broad["f_y"],
            spectrum_broad["bin_widths"],
        )
        assert means["y_F"] == pytest.approx(spectrum_broad["y_F"], rel=1e-6)
        assert means["y_D"] == pytest.approx(spectrum_broad["y_D"], rel=1e-6)


# ---------------------------------------------------------------------------
# Kappa table tests
# ---------------------------------------------------------------------------


class TestComputeKappaTable:
    def test_loads_and_computes(self, kappa_table):
        """Kappa table loads CSV files and produces values."""
        assert len(kappa_table["energy_MeV"]) > 20
        assert len(kappa_table["kappa"]) == len(kappa_table["energy_MeV"])

    def test_kappa_range(self, kappa_table):
        """Kappa values should be in physically reasonable range [0.3, 1.0]."""
        assert np.all(kappa_table["kappa"] > 0.3)
        assert np.all(kappa_table["kappa"] < 1.0)


# ---------------------------------------------------------------------------
# Tissue-equivalence correction tests
# ---------------------------------------------------------------------------


class TestBraggKappaMachinery:
    """Phase 27 machinery: Bragg-additivity SiC + source switch (data-blocked)."""

    def test_bragg_composer_weights(self):
        """SiC = 0.7004*Si + 0.2996*C (Bragg additivity), log-log interpolated."""
        from src.microdosimetry import sic_stopping_power_bragg

        e = np.array([1.0, 10.0, 100.0])
        s_si = np.array([100.0, 20.0, 4.0])
        s_c = np.array([130.0, 26.0, 5.2])
        grid, s_sic = sic_stopping_power_bragg(e, s_si, e, s_c)
        expected = 0.7004 * s_si + 0.2996 * s_c
        np.testing.assert_allclose(s_sic, expected, rtol=1e-6)

    def test_bragg_source_refuses_placeholder(self):
        """compute_kappa_table(source='bragg') must REFUSE placeholder data
        (no fabrication) until real PSTAR files are dropped in."""
        from src.microdosimetry import compute_kappa_table

        with pytest.raises(FileNotFoundError, match="PSTAR"):
            compute_kappa_table(source="bragg")

    def test_legacy_source_back_compat(self):
        """Legacy source still returns the (placeholder) flat table + source key."""
        from src.microdosimetry import compute_kappa_table

        res = compute_kappa_table(source="legacy")
        assert res["source"] == "legacy"
        assert len(res["kappa"]) == len(res["energy_MeV"])


class TestTissueEquivalenceCorrection:
    def test_constant_kappa(self):
        """y_tissue = kappa_constant * y_SiC for constant kappa."""
        y_sic = np.array([1.0, 2.0, 5.0, 10.0])
        energies = np.array([100.0, 200.0, 500.0, 1000.0])
        y_tissue = tissue_equivalence_correction(
            y_sic, energies, kappa_table=None, kappa_constant=0.58
        )
        np.testing.assert_allclose(y_tissue, 0.58 * y_sic)

    def test_table_kappa_uses_kinetic_energy(self, kappa_table):
        """AUDIT C-2: kappa(E) lookup uses particle KINETIC energy, scalar or
        per-event, NOT deposited energy."""
        y_sic = np.array([5.0, 5.0, 5.0])
        # Per-event KINETIC energies (MeV) -- correct variable
        Ek = np.array([1.0, 10.0, 100.0])
        y_tissue = tissue_equivalence_correction(
            y_sic, kappa_table=kappa_table, particle_energy_MeV=Ek
        )
        assert y_tissue.shape == y_sic.shape
        # scalar (mono-energetic beam) path also works
        y_beam = tissue_equivalence_correction(
            y_sic, kappa_table=kappa_table, particle_energy_MeV=62.0
        )
        assert y_beam.shape == y_sic.shape

    def test_no_kinetic_energy_falls_back_to_average(self, kappa_table):
        """AUDIT C-2: without particle_energy_MeV, fall back to energy-averaged
        kappa (documented approximation), NOT the deposited-energy lookup."""
        y_sic = np.array([5.0, 5.0])
        y_tissue = tissue_equivalence_correction(y_sic, kappa_table=kappa_table)
        kappa_avg = float(np.mean(kappa_table["kappa"]))
        np.testing.assert_allclose(y_tissue, kappa_avg * y_sic)


# ---------------------------------------------------------------------------
# Plot tests
# ---------------------------------------------------------------------------


class TestPlots:
    def test_plot_yd_spectrum_creates_figure(self, spectrum_broad):
        """plot_yd_spectrum runs without error and returns axes."""
        ax = plot_yd_spectrum(
            spectrum_broad["bin_centers"],
            spectrum_broad["d_y"],
            y_F=spectrum_broad["y_F"],
            y_D=spectrum_broad["y_D"],
            label="Test",
        )
        assert ax is not None
        plt.close("all")

    def test_plot_yf_spectrum_creates_figure(self, spectrum_broad):
        """plot_yf_spectrum runs without error and returns axes."""
        ax = plot_yf_spectrum(
            spectrum_broad["bin_centers"],
            spectrum_broad["f_y"],
            y_F=spectrum_broad["y_F"],
            label="Test",
        )
        assert ax is not None
        plt.close("all")

    def test_plot_yd_on_existing_axes(self, spectrum_broad):
        """plot_yd_spectrum can plot on provided axes."""
        fig, ax = plt.subplots()
        returned_ax = plot_yd_spectrum(
            spectrum_broad["bin_centers"],
            spectrum_broad["d_y"],
            ax=ax,
        )
        assert returned_ax is ax
        plt.close("all")

    def test_plot_yf_on_existing_axes(self, spectrum_broad):
        """plot_yf_spectrum can plot on provided axes."""
        fig, ax = plt.subplots()
        returned_ax = plot_yf_spectrum(
            spectrum_broad["bin_centers"],
            spectrum_broad["f_y"],
            ax=ax,
        )
        assert returned_ax is ax
        plt.close("all")
