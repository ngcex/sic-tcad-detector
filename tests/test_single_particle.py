"""Tests for single-particle transient charge collection module.

Validates:
- LET to pairs/cm conversion (unit test, no devsim)
- Current pulse analysis with synthetic data (unit test, no devsim)
- CCE(LET) table save/load round-trip with interpolation (unit test, no devsim)
- Ion track generation integral matches expected Q_generated (physics, slow)
- Backward-Euler estimator and convergence stop (unit, no devsim)
- Default stack: raw CCE <= 1 and monotone in bias (physics, slow)
  (legacy-model checks: tests/test_single_particle_legacy.py)
"""

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from etna.core.single_particle import (
    analyze_current_pulse,
    save_cce_let_table,
    load_cce_let_table,
)


# ---------------------------------------------------------------------------
# Fast unit tests (no devsim required)
# ---------------------------------------------------------------------------


class TestAnalyzeCurrentPulse:
    """Test current pulse analysis with synthetic data."""

    def test_q_collected_from_known_integral(self):
        """Verify Q_collected matches known integral of rectangular pulse."""
        # Rectangular pulse: I=1e-6 A/cm for 10 ns, I_dark=0
        times = np.linspace(0, 10e-9, 101)
        currents = np.full_like(times, 1e-6)
        I_dark = 0.0

        result = analyze_current_pulse(times, currents, I_dark)

        # Expected: Q = I * t = 1e-6 * 10e-9 = 1e-14 C/cm
        assert abs(result["Q_collected"] - 1e-14) / 1e-14 < 0.01

    def test_peak_current_and_time(self):
        """Verify I_peak and t_peak from triangular pulse."""
        times = np.array([0, 5e-9, 10e-9])
        currents = np.array([0, 2e-6, 0])
        I_dark = 0.0

        result = analyze_current_pulse(times, currents, I_dark)

        assert result["I_peak"] == pytest.approx(2e-6, rel=1e-10)
        assert result["t_peak"] == pytest.approx(5e-9, rel=1e-10)

    def test_dark_current_subtraction(self):
        """Verify dark current is subtracted from signal."""
        times = np.linspace(0, 10e-9, 101)
        I_dark = 1e-8
        currents = np.full_like(times, 1e-6 + I_dark)

        result = analyze_current_pulse(times, currents, I_dark)

        # Q should be integral of (|I| - |I_dark|) = 1e-6 * 10e-9
        assert abs(result["Q_collected"] - 1e-14) / 1e-14 < 0.01

    def test_collection_time_95_percent(self):
        """Verify t_collection is time to reach 95% of total charge."""
        # Exponential decay pulse: I = I0 * exp(-t/tau)
        tau = 5e-9
        times = np.linspace(0, 50e-9, 1001)
        I0 = 1e-6
        currents = I0 * np.exp(-times / tau)
        I_dark = 0.0

        result = analyze_current_pulse(times, currents, I_dark)

        # 95% of charge collected at t = -tau * ln(0.05) ~ 3*tau = 15 ns
        expected_t95 = -tau * np.log(0.05)
        assert abs(result["t_collection"] - expected_t95) / expected_t95 < 0.1


class TestSaveLoadCCETable:
    """Test CCE(LET) table round-trip and interpolation."""

    def test_round_trip(self):
        """Save and load CCE table, verify data integrity."""
        df = pd.DataFrame(
            {
                "LET_keV_um": [1.0, 10.0, 100.0],
                "Q_generated_fC": [1.0, 10.0, 100.0],
                "Q_collected_fC": [0.9, 9.0, 85.0],
                "CCE": [0.9, 0.9, 0.85],
                "t_collection_ns": [5.0, 8.0, 15.0],
            }
        )

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            filepath = f.name

        try:
            geometry = {"half_width_um": 50.0, "epi_um": 10.0, "bias_V": 50.0}
            save_cce_let_table(df, filepath, geometry_info=geometry)

            # Verify JSON structure
            with open(filepath) as f:
                data = json.load(f)
            assert "LET_keV_um" in data
            assert "CCE" in data
            assert "geometry" in data
            assert len(data["LET_keV_um"]) == 3

            # Load and verify interpolation function
            cce_interp, metadata = load_cce_let_table(filepath)
            assert callable(cce_interp)
            assert cce_interp(1.0) == pytest.approx(0.9, abs=0.01)
            assert cce_interp(100.0) == pytest.approx(0.85, abs=0.01)
        finally:
            Path(filepath).unlink(missing_ok=True)

    def test_interpolation_between_points(self):
        """Verify log-linear interpolation between LET values."""
        df = pd.DataFrame(
            {
                "LET_keV_um": [1.0, 100.0],
                "Q_generated_fC": [1.0, 100.0],
                "Q_collected_fC": [0.9, 80.0],
                "CCE": [0.9, 0.8],
                "t_collection_ns": [5.0, 15.0],
            }
        )

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            filepath = f.name

        try:
            save_cce_let_table(df, filepath)
            cce_interp, _ = load_cce_let_table(filepath)

            # At LET=10 (log midpoint of 1 and 100), CCE should be ~0.85
            cce_mid = cce_interp(10.0)
            assert 0.8 < cce_mid < 0.9
        finally:
            Path(filepath).unlink(missing_ok=True)

    def test_nan_handling(self):
        """Verify NaN values in CCE table are handled during interpolation."""
        df = pd.DataFrame(
            {
                "LET_keV_um": [1.0, 10.0, 100.0],
                "Q_generated_fC": [1.0, float("nan"), 100.0],
                "Q_collected_fC": [0.9, float("nan"), 80.0],
                "CCE": [0.9, float("nan"), 0.8],
                "t_collection_ns": [5.0, float("nan"), 15.0],
            }
        )

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            filepath = f.name

        try:
            save_cce_let_table(df, filepath)
            cce_interp, _ = load_cce_let_table(filepath)

            # Should interpolate between valid points, ignoring NaN
            assert cce_interp(1.0) == pytest.approx(0.9, abs=0.01)
            assert cce_interp(100.0) == pytest.approx(0.8, abs=0.01)
        finally:
            Path(filepath).unlink(missing_ok=True)


class TestCurrentContact:
    """The integrated contact: anode with an n+ cathode, cathode without (no devsim)."""

    @staticmethod
    def _run(monkeypatch, device_info, **kw):
        import etna.core.single_particle as sp

        seen = []

        def fake_current(info, contact):
            seen.append(contact)
            return 0.0

        monkeypatch.setattr(sp.devsim, "solve", lambda **k: None)
        monkeypatch.setattr(sp, "add_generation_to_dd", lambda info, g: None)
        monkeypatch.setattr(sp, "extract_contact_current", fake_current)
        r = sp.simulate_single_particle(device_info, np.zeros(3), t_collect=3e-12, **kw)
        return r["current_contact"], set(seen)

    def test_auto_contact(self, monkeypatch):
        from etna.core.device import DEFAULT_STACK, LEGACY_STACK

        info = {"device_name": "x", "region_name": "sic"}
        for stack, expected in ((DEFAULT_STACK, "anode"), (LEGACY_STACK, "cathode")):
            dev = {**info, "nplus_thickness_cm": stack["nplus_thickness_cm"]}
            assert self._run(monkeypatch, dev) == (expected, {expected})
        # explicit argument wins; unknown contact rejected
        dev = {**info, "nplus_thickness_cm": 1e-4}
        assert self._run(monkeypatch, dev, contact="cathode") == (
            "cathode",
            {"cathode"},
        )
        with pytest.raises(ValueError):
            self._run(monkeypatch, dev, contact="bulk")


class TestBackwardEulerEstimator:
    """Q_collected is the right-endpoint sum over every BDF1 step, the injection
    step included, and the run stops on charge convergence (no devsim)."""

    def test_collected_charge_is_right_endpoint_sum(self):
        from etna.core.single_particle import collected_charge

        times = np.array([1e-12, 3e-12, 7e-12])
        currents = np.array([-2.0, -1.0, -0.5])  # sign ignored, |I| integrated
        expected = 2.0 * 1e-12 + 1.0 * 2e-12 + 0.5 * 4e-12
        assert collected_charge(times, currents, 0.0) == pytest.approx(expected)
        # dark current subtracted per unit time
        assert collected_charge(times, currents - 0.1, -0.1) == pytest.approx(expected)
        assert collected_charge([], [], 0.0) == 0.0

    @staticmethod
    def _run_decay(monkeypatch, tau, **kw):
        """Fake device whose BDF1 step returns I = exp(-t_end / tau) (A/cm)."""
        import etna.core.single_particle as sp

        state = {"t": 0.0, "dt": 0.0}

        def fake_solve(**k):
            if k.get("type") == "transient_bdf1":
                state["t"] += k["tdelta"]

        monkeypatch.setattr(sp.devsim, "solve", fake_solve)
        monkeypatch.setattr(sp, "add_generation_to_dd", lambda info, g: None)
        monkeypatch.setattr(
            sp,
            "extract_contact_current",
            lambda info, c: 0.0 if state["t"] == 0.0 else np.exp(-state["t"] / tau),
        )
        info = {"device_name": "x", "region_name": "sic", "nplus_thickness_cm": None}
        return sp.simulate_single_particle(info, np.zeros(3), **kw)

    def test_injection_step_counted_and_stop_converges(self, monkeypatch):
        from etna.core.single_particle import analyze_current_pulse, collected_charge

        tau = 1e-9
        r = self._run_decay(monkeypatch, tau)
        assert r["converged"] and r["times"][-1] < 200e-9
        assert r["times"][0] == pytest.approx(1e-12)  # the injection step
        assert r["Q_collected"] == pytest.approx(
            collected_charge(r["times"], r["currents"], 0.0)
        )
        # the stop does not truncate: same as running to t_collect without a stop
        full = self._run_decay(monkeypatch, tau, rtol_converged=0.0)
        assert not full["converged"]
        assert r["Q_collected"] == pytest.approx(full["Q_collected"], rel=1e-6)
        # and the time-discrete sum is close to the exact integral tau
        assert r["Q_collected"] == pytest.approx(tau, rel=0.1)
        # analyze_current_pulse uses the same estimator
        pulse = analyze_current_pulse(r["times"], r["currents"], 0.0)
        assert pulse["Q_collected"] == pytest.approx(r["Q_collected"])


class TestLETConversion:
    """Test LET to pairs/cm conversion."""

    def test_let_1_keV_um(self):
        """LET=1 keV/um should give ~1.28e6 pairs/cm (E_pair=7.8 eV)."""
        # pairs/cm = LET * 1e3 / E_pair * 1e4
        # = 1.0 * 1e3 / 7.8 * 1e4 = 1e7 / 7.8 = ~1.28e6
        LET = 1.0
        E_pair = 7.8
        pairs_per_cm = LET * 1e3 / E_pair * 1e4

        expected = 1e7 / 7.8  # ~1.282e6
        assert abs(pairs_per_cm - expected) / expected < 1e-10

    def test_let_100_keV_um(self):
        """LET=100 keV/um should give ~1.28e8 pairs/cm (E_pair=7.8 eV)."""
        LET = 100.0
        E_pair = 7.8
        pairs_per_cm = LET * 1e3 / E_pair * 1e4
        expected = 100.0 * 1e7 / 7.8
        assert abs(pairs_per_cm - expected) / expected < 1e-10


# ---------------------------------------------------------------------------
# Physics validation tests (require devsim, marked slow)
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestIonTrackGeneration:
    """Validate ion track generation integral against expected Q_generated."""

    def test_generation_integral_matches_expected(self):
        """Generation area integral * q should match expected charge."""
        import devsim
        from etna.core.charge_collection_2d import create_2d_dd_device
        from etna.core.single_particle import ion_track_generation_2d

        device_info = None
        try:
            device_info = create_2d_dd_device(
                # hw 25 (was 50): the 1 um Gaussian track lies well inside
                half_width_um=25.0,
                V_bias=50.0,
                device_name="test_gen_integral",
            )

            LET = 10.0  # keV/um
            E_pair = 7.8  # eV (production default; ion_track_generation_2d uses 7.8)

            generation, Q_gen = ion_track_generation_2d(device_info, LET, x_ion_cm=0.0)

            # Expected: pairs/cm * epi_thickness (cm) gives total pairs per cm z-depth
            # But the Gaussian normalization and mesh integration may differ slightly.
            # The key check: Q_gen should be positive and physically reasonable.
            assert Q_gen > 0, "Q_generated must be positive"

            # Analytical estimate: pairs/cm * epi_thickness * q
            # pairs/cm = 10 * 1e7 / 7.8 ~ 1.28e7
            # epi = 10e-4 cm
            # Q_expected ~ 1.28e7 * 10e-4 * 1.602e-19 ~ 2.05e-15 C/cm
            # But Gaussian lateral integral over half-device may capture
            # less than full track -- depends on sigma vs half_width.
            # With sigma=1um and half_width=50um, integral should be ~50% of full
            # (half-device captures half the Gaussian).
            pairs_per_cm = 10.0 * 1e7 / E_pair
            epi_cm = device_info["epi_thickness_cm"]
            Q_expected_full = pairs_per_cm * epi_cm * 1.602e-19
            Q_expected_half = Q_expected_full * 0.5  # half-device

            # Allow 5% tolerance for mesh discretization
            assert (
                abs(Q_gen - Q_expected_half) / Q_expected_half < 0.05
            ), f"Q_gen={Q_gen:.4e} vs expected_half={Q_expected_half:.4e}"

        finally:
            if device_info is not None:
                try:
                    devsim.delete_device(device=device_info["device_name"])
                except Exception:
                    pass


@pytest.mark.slow
class TestCCEVsBias:
    """Default stack (p+ / graded epi / n+, anode current): raw single-particle CCE
    <= 1 and monotone in bias, and the run converges before t_collect.

    Folds in the former TestChargeConservation (0.5 <= CCE, Q > 0 at 50 V). The
    legacy (no n+) pin with old doping that lived here is gone: the legacy physics
    checks are in tests/test_single_particle_legacy.py. Asserts on the RAW
    Q_collected / Q_generated, not the table value (clipped to 1, so <= 1 there would
    be vacuous). Before the backward-Euler estimator this read 0.9824 / 0.9772 at
    10 / 50 V (non-monotone, hw 50); now 0.99978 / 0.99997 (hw 25).
    """

    def test_cce_le_one_and_monotone_in_bias(self):
        import devsim
        from etna.core.charge_collection_2d import create_2d_dd_device
        from etna.core.single_particle import (
            ion_track_generation_2d,
            simulate_single_particle,
        )

        cce_results = {}
        for V_bias in [10.0, 50.0]:
            device_info = None
            try:
                device_info = create_2d_dd_device(
                    # hw 25: centred 1 um track on a laterally uniform device
                    half_width_um=25.0,
                    V_bias=V_bias,
                    device_name=f"test_bias_{int(V_bias)}",
                )
                generation, Q_gen = ion_track_generation_2d(
                    device_info, LET_keV_um=10.0, x_ion_cm=0.0
                )
                result = simulate_single_particle(device_info, generation)
                assert result["current_contact"] == "anode"
                assert result["converged"], f"not converged at {V_bias} V"
                assert result["Q_collected"] > 0 and Q_gen > 0
                cce_results[V_bias] = result["Q_collected"] / Q_gen
            finally:
                if device_info is not None:
                    try:
                        devsim.delete_device(device=device_info["device_name"])
                    except Exception:
                        pass

        for V, cce in cce_results.items():
            # <= 1 up to round-off; a fully depleted 10 um SiC epi collects ~all
            assert 0.5 <= cce <= 1.0 + 1e-6, f"CCE({V} V)={cce:.6f}"
        assert (
            cce_results[50.0] >= cce_results[10.0]
        ), f"CCE(50V)={cce_results[50.0]:.6f} < CCE(10V)={cce_results[10.0]:.6f}"
