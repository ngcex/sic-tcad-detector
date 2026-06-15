"""Unit tests for optimization module scoring logic and noise estimation.

Tests pure-computation functions (estimate_noise_floor, score_structures)
that do not require devsim device creation. Integration tests for
microdosimetric_sweep and get_dark_current_2d are in the notebook.
"""

import numpy as np
import pandas as pd
import pytest

from src.optimization import (
    estimate_noise_floor,
    score_structures,
    _rank_sweep_results,
)


class TestSweepValidityGate:
    """Mj-3 full-depletion + CCE-floor ranking gate for the sweep results."""

    def test_deceptive_undepleted_ranks_below_valid(self):
        """A high-ratio but NOT-fully-depleted config must rank BELOW a
        lower-ratio but valid (fully-depleted, high-CCE) config. This is the
        exact deceptive-uniformity trap the audit flagged."""
        df = pd.DataFrame(
            [
                # deceptive: both-low CCE, perfect ratio, but undepleted
                {
                    "edge_center_ratio": 0.99,
                    "center_cce": 0.30,
                    "edge_cce": 0.297,
                    "is_fully_depleted": False,
                    "passes_cce_floor": False,
                    "is_valid": False,
                },
                # genuine: fully depleted, high CCE, slightly lower ratio
                {
                    "edge_center_ratio": 0.96,
                    "center_cce": 0.99,
                    "edge_cce": 0.95,
                    "is_fully_depleted": True,
                    "passes_cce_floor": True,
                    "is_valid": True,
                },
            ]
        )
        ranked = _rank_sweep_results(df)
        assert ranked.iloc[0]["is_valid"]
        assert ranked.iloc[0]["center_cce"] == 0.99
        assert not ranked.iloc[1]["is_valid"]

    def test_valid_configs_ordered_by_uniformity(self):
        """Within the valid group, higher edge_center_ratio ranks first."""
        df = pd.DataFrame(
            [
                {"edge_center_ratio": 0.93, "is_valid": True},
                {"edge_center_ratio": 0.98, "is_valid": True},
                {"edge_center_ratio": 0.50, "is_valid": False},
            ]
        )
        ranked = _rank_sweep_results(df)
        assert ranked.iloc[0]["edge_center_ratio"] == 0.98
        assert ranked.iloc[1]["edge_center_ratio"] == 0.93
        assert not ranked.iloc[2]["is_valid"]


# ---------------------------------------------------------------------------
# estimate_noise_floor tests
# ---------------------------------------------------------------------------


class TestEstimateNoiseFloor:
    """Tests for noise floor estimation from dark current."""

    def test_values_sanity(self):
        """Known I_dark=1e-11 A produces positive Q_min, y_min in sane range."""
        result = estimate_noise_floor(
            I_dark_A=1e-11,
            t_shaping_s=1e-6,
            sv_thickness_um=10.0,
            sv_width_um=100.0,
        )
        assert result["Q_min_fC"] > 0, "Q_min must be positive"
        assert result["y_min_keV_um"] > 0, "y_min must be positive"
        # Sanity: y_min should be well below 1 keV/um for ~10 pA dark current
        assert (
            result["y_min_keV_um"] < 1.0
        ), f"y_min={result['y_min_keV_um']:.4f} keV/um too high for 10 pA"
        assert result["E_min_keV"] > 0, "E_min must be positive"
        assert result["l_bar_um"] > 0, "l_bar must be positive"

    def test_scaling_with_dark_current(self):
        """Doubling I_dark increases sigma_shot by sqrt(2)."""
        r1 = estimate_noise_floor(I_dark_A=1e-11)
        r2 = estimate_noise_floor(I_dark_A=2e-11)

        ratio = r2["sigma_shot_C"] / r1["sigma_shot_C"]
        expected = np.sqrt(2.0)
        assert (
            abs(ratio - expected) < 0.01
        ), f"sigma_shot ratio={ratio:.4f}, expected sqrt(2)={expected:.4f}"

    def test_zero_dark_current(self):
        """Zero dark current gives zero noise floor."""
        result = estimate_noise_floor(I_dark_A=0.0)
        assert result["sigma_shot_C"] == 0.0
        assert result["Q_min_fC"] == 0.0
        assert result["y_min_keV_um"] == 0.0

    def test_slab_vs_3d_chord_length(self):
        """Slab approximation gives different l_bar than 3D geometry."""
        r_slab = estimate_noise_floor(
            I_dark_A=1e-11, sv_thickness_um=10.0, sv_width_um=100.0, sv_depth_um=None
        )
        r_3d = estimate_noise_floor(
            I_dark_A=1e-11, sv_thickness_um=10.0, sv_width_um=100.0, sv_depth_um=100.0
        )
        # Slab l_bar = 2*t = 20 um; 3D l_bar = 4V/S < 20 um
        assert r_slab["l_bar_um"] == 20.0, "Slab l_bar should be 2*thickness"
        assert (
            r_3d["l_bar_um"] < r_slab["l_bar_um"]
        ), "3D chord length should be smaller than slab approximation"


# ---------------------------------------------------------------------------
# score_structures tests
# ---------------------------------------------------------------------------


class TestScoreStructures:
    """Tests for multi-criteria structure scoring."""

    @pytest.fixture
    def three_structures(self):
        """Three structures with known metrics for deterministic ranking."""
        return {
            "planar": {
                "cce_uniformity": 0.85,
                "noise_floor": 0.05,
                "spectral_resolution": 0.3,
                "fabrication_complexity": 1,
            },
            "guard_ring": {
                "cce_uniformity": 0.95,
                "noise_floor": 0.04,
                "spectral_resolution": 0.25,
                "fabrication_complexity": 2,
            },
            "3d_electrode": {
                "cce_uniformity": 0.99,
                "noise_floor": 0.03,
                "spectral_resolution": 0.20,
                "fabrication_complexity": 4,
            },
        }

    def test_ranking_order(self, three_structures):
        """Verify ranking order matches expected for default weights."""
        df = score_structures(three_structures)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert "weighted_score" in df.columns

        # 3d_electrode ranks highest with default weights: best in 3 of 4
        # metrics (cce, noise, spectral_res), only loses on fabrication
        top = df.iloc[0]["structure"]
        assert top == "3d_electrode", (
            f"Expected 3d_electrode at top, got {top}. "
            f"Scores: {df[['structure', 'weighted_score']].to_dict()}"
        )
        # Planar should be last (worst in 3 of 4 metrics)
        last = df.iloc[-1]["structure"]
        assert last == "planar", f"Expected planar last, got {last}"

    def test_normalization_range(self, three_structures):
        """All normalized scores must be in [0, 1]."""
        df = score_structures(three_structures)
        norm_cols = [c for c in df.columns if c.endswith("_norm")]
        assert len(norm_cols) == 4, f"Expected 4 normalized columns, got {norm_cols}"

        for col in norm_cols:
            vals = df[col].values
            assert np.all(vals >= 0.0), f"{col} has values < 0: {vals}"
            assert np.all(vals <= 1.0), f"{col} has values > 1: {vals}"

    def test_equal_inputs_equal_scores(self):
        """All structures with identical metrics get equal scores."""
        equal_metrics = {
            "A": {
                "cce_uniformity": 0.9,
                "noise_floor": 0.05,
                "spectral_resolution": 0.3,
                "fabrication_complexity": 2,
            },
            "B": {
                "cce_uniformity": 0.9,
                "noise_floor": 0.05,
                "spectral_resolution": 0.3,
                "fabrication_complexity": 2,
            },
            "C": {
                "cce_uniformity": 0.9,
                "noise_floor": 0.05,
                "spectral_resolution": 0.3,
                "fabrication_complexity": 2,
            },
        }
        df = score_structures(equal_metrics)
        scores = df["weighted_score"].values
        assert np.allclose(
            scores, scores[0]
        ), f"Equal inputs should give equal scores: {scores}"

    def test_custom_weights(self, three_structures):
        """Custom weights change the ranking."""
        # Heavily weight fabrication -> planar should win
        weights = {
            "cce_uniformity": 0.10,
            "noise_floor": 0.10,
            "spectral_resolution": 0.10,
            "fabrication_complexity": 0.70,
        }
        df = score_structures(three_structures, weights=weights)
        top = df.iloc[0]["structure"]
        assert (
            top == "planar"
        ), f"With high fabrication weight, planar should win, got {top}"

    def test_output_columns(self, three_structures):
        """Output DataFrame has expected column structure."""
        df = score_structures(three_structures)
        expected_cols = {
            "structure",
            "cce_uniformity",
            "noise_floor",
            "spectral_resolution",
            "fabrication_complexity",
            "cce_uniformity_norm",
            "noise_floor_norm",
            "spectral_resolution_norm",
            "fabrication_complexity_norm",
            "weighted_score",
        }
        assert (
            set(df.columns) == expected_cols
        ), f"Missing columns: {expected_cols - set(df.columns)}"
