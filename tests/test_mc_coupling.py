"""Tests for Monte Carlo coupling module.

Validates:
- CSV import with default and custom column mappings
- Unit conversion (mm/MeV -> cm/keV, unknown units raise ValueError)
- ROOT import via mocking (no real ROOT files needed)
- process_mc_ensemble: shape, zero-energy filtering, constant CCE
- pulse_height_distribution: shape, total counts, log-spacing
- events_to_charge_profiles: returns correct structure, calls ion_track_generation_2d
- Integration test with real CCE table (marked slow)
"""

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

from src.mc_coupling import (
    convert_units,
    events_to_charge_profiles,
    list_root_trees,
    load_mc_events_csv,
    load_mc_events_root,
    process_mc_ensemble,
    pulse_height_distribution,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_events_df():
    """DataFrame with 50 events, 5 steps each (250 rows)."""
    rng = np.random.default_rng(42)
    n_events = 50
    n_steps = 5
    rows = []
    for eid in range(n_events):
        for _ in range(n_steps):
            rows.append(
                {
                    "event_id": eid,
                    "x_cm": rng.uniform(-0.005, 0.005),
                    "y_cm": rng.uniform(-0.001, 0.0),
                    "z_cm": rng.uniform(-0.005, 0.005),
                    "edep_keV": rng.uniform(10, 100),
                }
            )
    return pd.DataFrame(rows)


@pytest.fixture
def sample_csv_file(sample_events_df, tmp_path):
    """Write sample events to a CSV with default column names."""
    # CSV uses the default column mapping: event_id, x, y, z, edep
    df = sample_events_df.rename(
        columns={
            "x_cm": "x",
            "y_cm": "y",
            "z_cm": "z",
            "edep_keV": "edep",
        }
    )
    path = tmp_path / "events.csv"
    df.to_csv(path, index=False)
    return path


@pytest.fixture
def mock_cce_interp():
    """Constant CCE = 0.95 for all LET values."""
    return lambda let: 0.95


# ---------------------------------------------------------------------------
# Unit conversion tests
# ---------------------------------------------------------------------------


class TestConvertUnits:
    def test_mm_mev(self):
        """Verify mm -> cm (0.1) and MeV -> keV (1000)."""
        df = pd.DataFrame(
            {
                "event_id": [0],
                "x_cm": [10.0],
                "y_cm": [20.0],
                "z_cm": [30.0],
                "edep_keV": [0.5],
            }
        )
        convert_units(df, pos_unit="mm", energy_unit="MeV")
        assert df["x_cm"].iloc[0] == pytest.approx(1.0)
        assert df["y_cm"].iloc[0] == pytest.approx(2.0)
        assert df["z_cm"].iloc[0] == pytest.approx(3.0)
        assert df["edep_keV"].iloc[0] == pytest.approx(500.0)

    def test_unknown_pos_unit_raises(self):
        df = pd.DataFrame({"x_cm": [1.0], "edep_keV": [1.0]})
        with pytest.raises(ValueError, match="Unknown position unit"):
            convert_units(df, pos_unit="inches", energy_unit="keV")

    def test_unknown_energy_unit_raises(self):
        df = pd.DataFrame({"x_cm": [1.0], "edep_keV": [1.0]})
        with pytest.raises(ValueError, match="Unknown energy unit"):
            convert_units(df, pos_unit="cm", energy_unit="joules")

    @pytest.mark.parametrize(
        "pos_unit,factor",
        [
            ("mm", 0.1),
            ("cm", 1.0),
            ("um", 1e-4),
        ],
    )
    def test_position_factors(self, pos_unit, factor):
        df = pd.DataFrame(
            {"x_cm": [100.0], "y_cm": [100.0], "z_cm": [100.0], "edep_keV": [1.0]}
        )
        convert_units(df, pos_unit=pos_unit, energy_unit="keV")
        assert df["x_cm"].iloc[0] == pytest.approx(100.0 * factor)

    @pytest.mark.parametrize(
        "energy_unit,factor",
        [
            ("MeV", 1e3),
            ("keV", 1.0),
            ("eV", 1e-3),
        ],
    )
    def test_energy_factors(self, energy_unit, factor):
        df = pd.DataFrame({"x_cm": [1.0], "edep_keV": [100.0]})
        convert_units(df, pos_unit="cm", energy_unit=energy_unit)
        assert df["edep_keV"].iloc[0] == pytest.approx(100.0 * factor)


# ---------------------------------------------------------------------------
# CSV import tests
# ---------------------------------------------------------------------------


class TestLoadCSV:
    def test_default_columns(self, sample_csv_file, sample_events_df):
        """Load CSV with default column names, verify shape and columns."""
        df = load_mc_events_csv(sample_csv_file)
        assert set(df.columns) >= {"event_id", "x_cm", "y_cm", "z_cm", "edep_keV"}
        assert len(df) == 250
        assert df["event_id"].nunique() == 50

    def test_custom_column_map(self, tmp_path):
        """Load CSV with non-standard column names via custom mapping."""
        df_src = pd.DataFrame(
            {
                "evt": [0, 0, 1, 1],
                "px": [1.0, 2.0, 3.0, 4.0],
                "py": [0.1, 0.2, 0.3, 0.4],
                "pz": [0.0, 0.0, 0.0, 0.0],
                "energy": [50.0, 60.0, 70.0, 80.0],
            }
        )
        path = tmp_path / "custom.csv"
        df_src.to_csv(path, index=False)

        col_map = {
            "event_id": "evt",
            "x_cm": "px",
            "y_cm": "py",
            "z_cm": "pz",
            "edep_keV": "energy",
        }
        df = load_mc_events_csv(path, column_map=col_map)
        assert "event_id" in df.columns
        assert "x_cm" in df.columns
        assert len(df) == 4

    def test_unit_conversion(self, tmp_path):
        """CSV in mm/MeV: verify values converted to cm/keV."""
        df_src = pd.DataFrame(
            {
                "event_id": [0, 1],
                "x": [10.0, 20.0],  # mm
                "y": [5.0, 10.0],  # mm
                "z": [0.0, 0.0],  # mm
                "edep": [0.5, 1.0],  # MeV
            }
        )
        path = tmp_path / "mm_mev.csv"
        df_src.to_csv(path, index=False)

        df = load_mc_events_csv(path, pos_unit="mm", energy_unit="MeV")
        assert df["x_cm"].iloc[0] == pytest.approx(1.0)  # 10 mm * 0.1
        assert df["edep_keV"].iloc[0] == pytest.approx(500.0)  # 0.5 MeV * 1000


# ---------------------------------------------------------------------------
# process_mc_ensemble tests
# ---------------------------------------------------------------------------


class TestProcessMCEnsemble:
    def test_shape(self, sample_events_df, mock_cce_interp):
        """Output dict has correct keys and array lengths."""
        result = process_mc_ensemble(
            sample_events_df, mock_cce_interp, sv_thickness_um=10.0
        )
        assert set(result.keys()) == {
            "event_energies_keV",
            "event_LET_keV_um",
            "event_CCE",
            "event_collected_keV",
            "n_events",
            "n_zero_energy",
        }
        assert len(result["event_energies_keV"]) == 50
        assert result["n_events"] == 50

    def test_zero_energy_filtering(self, mock_cce_interp):
        """Zero-edep events are filtered and counted."""
        df = pd.DataFrame(
            {
                "event_id": [0, 1, 2, 2],
                "x_cm": [0, 0, 0, 0],
                "y_cm": [0, 0, 0, 0],
                "z_cm": [0, 0, 0, 0],
                "edep_keV": [50.0, 0.0, 30.0, 20.0],
            }
        )
        result = process_mc_ensemble(df, mock_cce_interp)
        assert result["n_zero_energy"] == 1  # event 1 has edep=0
        assert result["n_events"] == 2  # events 0 and 2

    def test_constant_cce(self, sample_events_df, mock_cce_interp):
        """With CCE=0.95, collected = 0.95 * deposited for all events."""
        result = process_mc_ensemble(sample_events_df, mock_cce_interp)
        np.testing.assert_allclose(
            result["event_collected_keV"],
            0.95 * result["event_energies_keV"],
            rtol=1e-10,
        )


# ---------------------------------------------------------------------------
# Pulse height distribution tests
# ---------------------------------------------------------------------------


class TestPulseHeightDistribution:
    def test_shape(self):
        """Bin centers, counts, bin_edges have correct lengths."""
        energies = np.random.default_rng(0).uniform(10, 500, size=200)
        result = pulse_height_distribution(energies, n_bins=50)
        assert len(result["bin_centers_keV"]) == 50
        assert len(result["counts"]) == 50
        assert len(result["bin_edges_keV"]) == 51

    def test_total_counts(self):
        """Sum of counts equals number of positive-energy events."""
        energies = np.concatenate(
            [
                np.random.default_rng(1).uniform(10, 500, size=150),
                np.zeros(10),  # zero-energy events excluded
            ]
        )
        result = pulse_height_distribution(energies, n_bins=100)
        assert result["counts"].sum() == 150

    def test_log_spacing(self):
        """Bin edges are log-spaced (constant ratio)."""
        energies = np.random.default_rng(2).uniform(10, 500, size=100)
        result = pulse_height_distribution(energies, n_bins=20)
        edges = result["bin_edges_keV"]
        ratios = edges[1:] / edges[:-1]
        np.testing.assert_allclose(ratios, ratios[0], rtol=1e-10)


# ---------------------------------------------------------------------------
# events_to_charge_profiles tests
# ---------------------------------------------------------------------------


class TestEventsToChargeProfiles:
    def test_returns_list(self):
        """Returns list of dicts with correct keys, one per event."""
        df = pd.DataFrame(
            {
                "event_id": [0, 0, 1, 2, 2],
                "x_cm": [0, 0, 0, 0, 0],
                "y_cm": [0, 0, 0, 0, 0],
                "z_cm": [0, 0, 0, 0, 0],
                "edep_keV": [50.0, 30.0, 100.0, 40.0, 60.0],
            }
        )
        mock_device = {"device": "test", "region": "test"}

        with patch("src.single_particle.ion_track_generation_2d") as mock_itg:
            mock_itg.return_value = (np.zeros(100), 1e-15)
            profiles = events_to_charge_profiles(df, mock_device, sv_thickness_um=10.0)

        assert len(profiles) == 3  # 3 unique events
        for p in profiles:
            assert set(p.keys()) == {
                "event_id",
                "LET_keV_um",
                "generation",
                "Q_generated_C_per_cm",
            }
            assert isinstance(p["generation"], np.ndarray)

    def test_calls_ion_track_with_correct_let(self):
        """ion_track_generation_2d called with LET = total_edep / sv_thickness."""
        df = pd.DataFrame(
            {
                "event_id": [0, 0],
                "x_cm": [0, 0],
                "y_cm": [0, 0],
                "z_cm": [0, 0],
                "edep_keV": [50.0, 30.0],  # total = 80 keV
            }
        )
        mock_device = {"device": "test"}
        sv_um = 10.0
        expected_let = 80.0 / sv_um  # 8.0 keV/um

        with patch("src.single_particle.ion_track_generation_2d") as mock_itg:
            mock_itg.return_value = (np.zeros(100), 1e-15)
            events_to_charge_profiles(df, mock_device, sv_thickness_um=sv_um)

        assert mock_itg.call_count == 1
        call_args = mock_itg.call_args
        assert call_args[0][1] == pytest.approx(expected_let)


# ---------------------------------------------------------------------------
# ROOT tests (mocked -- no real ROOT files needed)
# ---------------------------------------------------------------------------


class TestROOTImport:
    def test_list_root_trees(self):
        """list_root_trees filters to TTree objects only."""
        mock_tree = MagicMock()
        mock_tree.num_entries = 1000

        mock_histogram = MagicMock(spec=[])  # no num_entries attribute

        mock_file = MagicMock()
        mock_file.keys.return_value = ["Hits;1", "RunInfo;1"]
        mock_file.__getitem__ = lambda self, k: {
            "Hits;1": mock_tree,
            "RunInfo;1": mock_histogram,
        }[k]
        mock_file.__enter__ = lambda self: self
        mock_file.__exit__ = MagicMock(return_value=False)

        mock_uproot = MagicMock()
        mock_uproot.open.return_value = mock_file

        with patch.dict("sys.modules", {"uproot": mock_uproot}):
            trees = list_root_trees("fake.root")

        assert "Hits" in trees
        # RunInfo has no num_entries so should be excluded
        assert "RunInfo" not in trees

    def test_load_root_events(self):
        """Mock uproot.open to verify DataFrame construction and unit conversion."""
        n = 100
        mock_arrays = {
            "EventID": np.arange(n),
            "PosX": np.random.default_rng(0).uniform(0, 10, n),  # mm
            "PosY": np.random.default_rng(1).uniform(-5, 0, n),  # mm
            "PosZ": np.random.default_rng(2).uniform(0, 10, n),  # mm
            "Edep": np.random.default_rng(3).uniform(0.01, 1.0, n),  # MeV
        }

        mock_tree = MagicMock()
        mock_tree.keys.return_value = ["EventID", "PosX", "PosY", "PosZ", "Edep"]
        mock_tree.arrays.return_value = mock_arrays

        mock_file = MagicMock()
        mock_file.__getitem__ = lambda self, k: mock_tree
        mock_file.__enter__ = lambda self: self
        mock_file.__exit__ = MagicMock(return_value=False)

        mock_uproot = MagicMock()
        mock_uproot.open.return_value = mock_file

        with patch.dict("sys.modules", {"uproot": mock_uproot}):
            df = load_mc_events_root("fake.root", tree_name="Hits")

        assert len(df) == n
        assert "event_id" in df.columns
        assert "x_cm" in df.columns
        # Verify unit conversion: mm -> cm (factor 0.1)
        np.testing.assert_allclose(
            df["x_cm"].values,
            mock_arrays["PosX"] * 0.1,
            rtol=1e-10,
        )
        # MeV -> keV (factor 1000)
        np.testing.assert_allclose(
            df["edep_keV"].values,
            mock_arrays["Edep"] * 1000,
            rtol=1e-10,
        )


# ---------------------------------------------------------------------------
# Integration test (slow)
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestFullPipeline:
    def test_with_cce_table(self):
        """Full pipeline: synthetic events -> process_mc_ensemble -> PHD."""
        import json
        import time

        cce_path = Path("data/cce_let_table_100um.json")
        if not cce_path.exists():
            pytest.skip("CCE table not found at data/cce_let_table_100um.json")

        from src.single_particle import load_cce_let_table

        cce_interp, metadata = load_cce_let_table(cce_path)

        # Synthetic 1000-event dataset
        rng = np.random.default_rng(123)
        n_events = 1000
        sv_um = 10.0
        # LET range 1-100 keV/um -> edep = LET * sv_um = 10-1000 keV
        event_edeps = rng.uniform(10, 1000, size=n_events)

        df = pd.DataFrame(
            {
                "event_id": np.arange(n_events),
                "x_cm": np.zeros(n_events),
                "y_cm": np.zeros(n_events),
                "z_cm": np.zeros(n_events),
                "edep_keV": event_edeps,
            }
        )

        t0 = time.perf_counter()
        result = process_mc_ensemble(df, cce_interp, sv_thickness_um=sv_um)
        phd = pulse_height_distribution(result["event_collected_keV"])
        elapsed = time.perf_counter() - t0

        # (a) All 1000 events processed
        assert result["n_events"] == 1000

        # (b) CCE values in [0, 1.05]
        assert np.all(result["event_CCE"] >= 0)
        assert np.all(result["event_CCE"] <= 1.05)

        # (c) PHD has non-zero counts
        assert phd["counts"].sum() > 0

        # (d) Completes in < 1 second
        assert elapsed < 1.0, f"Pipeline took {elapsed:.2f}s (limit: 1s)"
