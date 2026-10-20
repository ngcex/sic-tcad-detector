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

from etna.core.microdosimetry import (
    compute_kappa_table,
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
    """Real proton kappa(E) table (NIST PSTAR, SiC by Bragg additivity, LET ratio)."""
    return compute_kappa_table(source="bragg")


@pytest.fixture
def legacy_kappa_table():
    """The FABRICATED flat ~0.58 legacy table (explicit opt-in; C-1). Only used to pin
    that the placeholder stays reachable and stays labelled as such."""
    return compute_kappa_table(source="legacy_unsafe")


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
# Kappa table tests
# ---------------------------------------------------------------------------


class TestComputeKappaTable:
    def test_loads_and_computes(self, kappa_table):
        """Real NIST-based kappa table: grid and metadata."""
        assert len(kappa_table["energy_MeV"]) > 100
        assert len(kappa_table["kappa"]) == len(kappa_table["energy_MeV"])
        assert kappa_table["source"] == "bragg"
        assert kappa_table["particle"] == "proton"
        assert kappa_table["ratio"] == "linear"

    def test_real_kappa_range_and_energy_dependence(self, kappa_table):
        """SiC -> water LET ratio from NIST PSTAR: 0.426 @1 MeV -> 0.375 @150 MeV,
        monotonic, NOT flat (C-1 fixed with real data)."""
        E = kappa_table["energy_MeV"]
        k = kappa_table["kappa"]
        sel = (E >= 1.0) & (E <= 250.0)
        assert np.all((k[sel] > 0.37) & (k[sel] < 0.43))
        assert np.all(np.diff(k[sel]) < 0)
        assert k[sel][0] / k[sel][-1] > 1.12

    def test_mass_ratio_option(self):
        """ratio='mass' is the Bragg-Gray dose factor (S/rho ratio), > 1."""
        km = compute_kappa_table(source="bragg", ratio="mass")
        E = km["energy_MeV"]
        sel = (E >= 1.0) & (E <= 250.0)
        assert np.all((km["kappa"][sel] > 1.19) & (km["kappa"][sel] < 1.37))

    def test_legacy_unsafe_is_placeholder(self, legacy_kappa_table):
        """AUDIT C-1: the legacy_unsafe table is the FABRICATED flat placeholder
        (~0.58). Real SiC->water LET ratio is ~0.38-0.43 (protons); this pins the
        known-wrong placeholder range, it is not an endorsement."""
        k = legacy_kappa_table["kappa"]
        assert np.all(k > 0.5) and np.all(k < 0.65)
        assert k.max() / k.min() < 1.05  # flat: the defect

    def test_no_source_raises(self):
        """v6 physics-honesty fix: omitting source= must raise, not silently
        return the fabricated placeholder."""
        with pytest.raises(ValueError, match="source"):
            compute_kappa_table()

    def test_invalid_source_raises(self):
        """Any source other than 'bragg'/'legacy_unsafe' must raise."""
        with pytest.raises(ValueError, match="source"):
            compute_kappa_table(source="legacy")


class TestBraggKappaMachinery:
    """Bragg-additivity SiC + source switch."""

    def test_bragg_composer_weights(self):
        """SiC = 0.7004*Si + 0.2996*C (Bragg additivity), log-log interpolated."""
        from etna.core.microdosimetry import sic_stopping_power_bragg

        e = np.array([1.0, 10.0, 100.0])
        s_si = np.array([100.0, 20.0, 4.0])
        s_c = np.array([130.0, 26.0, 5.2])
        grid, s_sic = sic_stopping_power_bragg(e, s_si, e, s_c)
        expected = 0.7004 * s_si + 0.2996 * s_c
        np.testing.assert_allclose(s_sic, expected, rtol=1e-6)

    def test_bragg_source_uses_real_nist_data(self):
        """compute_kappa_table(source='bragg') == stopping_power.kappa_table (NIST)."""
        from etna.core.stopping_power import kappa_table as sp_kappa

        res = compute_kappa_table(source="bragg", particle="alpha", tissue="muscle_icru")
        ref = sp_kappa("alpha", tissue="muscle_icru")
        np.testing.assert_allclose(res["kappa"], ref["kappa"])

    def test_legacy_unsafe_requires_explicit_opt_in(self, legacy_kappa_table):
        """The fabricated flat-kappa table is only reachable via the explicit
        source='legacy_unsafe' -- the old silent source='legacy' default is gone."""
        assert legacy_kappa_table["source"] == "legacy_unsafe"
        assert len(legacy_kappa_table["kappa"]) == len(legacy_kappa_table["energy_MeV"])


class TestTissueEquivalenceCorrection:
    def test_constant_kappa(self):
        """y_tissue = kappa_constant * y_SiC for a caller-supplied constant kappa."""
        y_sic = np.array([1.0, 2.0, 5.0, 10.0])
        energies = np.array([100.0, 200.0, 500.0, 1000.0])
        y_tissue = tissue_equivalence_correction(
            y_sic, energies, kappa_table=None, kappa_constant=1.2
        )
        np.testing.assert_allclose(y_tissue, 1.2 * y_sic)

    def test_constant_kappa_legacy_unsafe(self):
        """kappa_constant='legacy_unsafe' explicitly opts into the fabricated
        0.58 placeholder (v6 physics-honesty fix: no longer the silent default)."""
        y_sic = np.array([1.0, 2.0, 5.0, 10.0])
        y_tissue = tissue_equivalence_correction(
            y_sic, kappa_table=None, kappa_constant="legacy_unsafe"
        )
        np.testing.assert_allclose(y_tissue, 0.58 * y_sic)

    def test_no_kappa_table_no_constant_raises(self):
        """v6 physics-honesty fix: omitting both kappa_table and kappa_constant
        must raise, not silently fall back to the fabricated 0.58."""
        y_sic = np.array([1.0, 2.0, 5.0, 10.0])
        with pytest.raises(ValueError, match="kappa_constant"):
            tissue_equivalence_correction(y_sic, kappa_table=None)

    def test_table_kappa_uses_kinetic_energy(self, kappa_table):
        """AUDIT C-2: kappa(E) lookup uses particle KINETIC energy, scalar or
        per-event, NOT deposited energy."""
        y_sic = np.array([5.0, 5.0, 5.0])
        Ek = np.array([1.0, 10.0, 100.0])  # per-event KINETIC energies (MeV)
        y_tissue = tissue_equivalence_correction(
            y_sic, kappa_table=kappa_table, particle_energy_MeV=Ek
        )
        assert y_tissue.shape == y_sic.shape
        # real kappa(E): lower kinetic energy -> larger LET ratio
        assert y_tissue[0] > y_tissue[1] > y_tissue[2]
        np.testing.assert_allclose(y_tissue / y_sic, [0.426, 0.391, 0.376], atol=2e-3)
        # scalar (mono-energetic beam) path
        y_beam = tissue_equivalence_correction(
            y_sic, kappa_table=kappa_table, particle_energy_MeV=62.0
        )
        np.testing.assert_allclose(y_beam / y_sic, 0.378, atol=2e-3)

    def test_deposited_energy_is_not_used_for_lookup(self, kappa_table):
        """C-2: passing deposited energies (keV) must not change the result."""
        y_sic = np.array([5.0, 5.0])
        a = tissue_equivalence_correction(
            y_sic, [10.0, 900.0], kappa_table=kappa_table, particle_energy_MeV=62.0
        )
        b = tissue_equivalence_correction(
            y_sic, None, kappa_table=kappa_table, particle_energy_MeV=62.0
        )
        np.testing.assert_allclose(a, b)

    def test_kinetic_energy_outside_table_raises(self, kappa_table):
        with pytest.raises(ValueError, match="kappa table"):
            tissue_equivalence_correction(
                [1.0], kappa_table=kappa_table, particle_energy_MeV=[np.nan]
            )
        with pytest.raises(ValueError, match="kappa table"):
            tissue_equivalence_correction(
                [1.0], kappa_table=kappa_table, particle_energy_MeV=1e6
            )

    def test_particle_mismatch_raises(self, kappa_table):
        with pytest.raises(ValueError, match="particle-specific"):
            tissue_equivalence_correction(
                [1.0], kappa_table=kappa_table, particle_energy_MeV=5.0, particle="alpha"
            )

    def test_no_kinetic_energy_falls_back_to_average(self, kappa_table):
        """AUDIT C-2: without particle_energy_MeV, fall back to the table-mean
        kappa (documented crude approximation), NOT the deposited-energy lookup."""
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
        assert isinstance(ax, plt.Axes) and ax.has_data()
        plt.close("all")

    def test_plot_yf_spectrum_creates_figure(self, spectrum_broad):
        """plot_yf_spectrum runs without error and returns axes."""
        ax = plot_yf_spectrum(
            spectrum_broad["bin_centers"],
            spectrum_broad["f_y"],
            y_F=spectrum_broad["y_F"],
            label="Test",
        )
        assert isinstance(ax, plt.Axes) and ax.has_data()
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
