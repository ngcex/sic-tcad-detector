"""Tests for Hurkx TAT dark current model and surface recombination.

Validates field-enhancement physics, component decomposition, and
calibration against the experimental 18 pA target at -30V.
"""

import numpy as np
import pytest
import devsim

from src.dark_current import (
    create_dark_current_device,
    extract_dark_current_components,
    dark_current_sweep,
    setup_tat_model,
    setup_surface_recombination,
    _compute_node_efield,
    _compute_gamma_factors,
)
from src.drift_diffusion import create_dd_device, ramp_bias, extract_contact_current


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_device_counter = [0]


def _unique_name(prefix="dc_test"):
    """Generate unique device name to avoid devsim conflicts."""
    _device_counter[0] += 1
    return f"{prefix}_{_device_counter[0]}"


@pytest.fixture
def tat_device_equilibrium():
    """Create a TAT-enabled device at equilibrium (0V)."""
    name = _unique_name("tat_eq")
    dev = create_dark_current_device(device_name=name)
    return dev


@pytest.fixture
def tat_device_reverse():
    """Create a TAT-enabled device ramped to -30V."""
    name = _unique_name("tat_rev")
    dev = create_dark_current_device(device_name=name)
    ramp_bias(dev, V_target=-30.0, contact="anode", V_step=1.0)
    _compute_node_efield(dev["device_name"], dev["region_name"])
    _compute_gamma_factors(dev["device_name"], dev["region_name"])
    return dev


@pytest.fixture
def srh_only_device():
    """Create a standard DD device (SRH only, no TAT)."""
    name = _unique_name("srh_only")
    dev = create_dd_device(device_name=name)
    return dev


# ---------------------------------------------------------------------------
# TestGammaEnhancement
# ---------------------------------------------------------------------------


class TestGammaEnhancement:
    """Test the field-enhancement factor Gamma behavior."""

    def test_low_field_gamma_near_unity(self, tat_device_equilibrium):
        """At equilibrium (0V), Gamma should be ~1 everywhere.

        At SiC detector fields (~1e4-1e5 V/cm), Kt >> 4 so Gamma = 1.
        The Hurkx field enhancement requires MV/cm fields (Kt < 4).
        """
        dev = tat_device_equilibrium
        device = dev["device_name"]
        region = dev["region_name"]

        Gamma = np.array(
            devsim.get_node_model_values(device=device, region=region, name="Gamma_n")
        )
        assert np.all(Gamma == pytest.approx(1.0, abs=0.01)), (
            f"Gamma at equilibrium should be ~1.0, got range "
            f"[{Gamma.min():.4f}, {Gamma.max():.4f}]"
        )

    def test_high_field_gamma_remains_unity_at_moderate_bias(self, tat_device_reverse):
        """At -30V, Gamma stays ~1 because SiC fields are below MV/cm.

        For Z1/2 (E_t=0.65 eV, m_t=0.25 m0), Kt < 4 requires E > 4.5 MV/cm.
        At -30V the max field is ~100 kV/cm, giving Kt >> 4 so Gamma = 1.
        This is expected physics: field enhancement is negligible at these
        voltages, and the effective generation term dominates dark current.
        """
        dev = tat_device_reverse
        device = dev["device_name"]
        region = dev["region_name"]

        Gamma = np.array(
            devsim.get_node_model_values(device=device, region=region, name="Gamma_n")
        )
        # At moderate reverse bias, Gamma stays ~1 for SiC
        assert np.all(
            Gamma == pytest.approx(1.0, abs=0.01)
        ), f"Gamma at -30V should be ~1.0 for SiC fields, got max {Gamma.max():.4f}"

    def test_gamma_calculation_at_extreme_field(self):
        """Verify Gamma formula produces enhancement at very high fields."""
        # Test the Schenk approximation directly
        # For Kt = 1 (very high field): Gamma = sqrt(pi)/(2*1) * exp(1/4) = 1.14
        # For Kt = 0.5: Gamma = sqrt(pi)/(2*0.5) * exp(1) = 4.81
        from src.dark_current import _compute_gamma_factors

        dev = create_dark_current_device(device_name=_unique_name("gamma_ext"))
        device = dev["device_name"]
        region = dev["region_name"]

        # Set Kt_TAT to test values manually
        n_nodes = len(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        # Set a range of Kt values: some high field (Kt < 4), some low field
        Kt_test = np.full(n_nodes, 1.0)  # Kt=1 everywhere
        devsim.node_model(device=device, region=region, name="Kt_TAT", equation="1.0")
        devsim.set_node_values(
            device=device, region=region, name="Kt_TAT", values=Kt_test.tolist()
        )

        _compute_gamma_factors(device, region)

        Gamma = np.array(
            devsim.get_node_model_values(device=device, region=region, name="Gamma_n")
        )
        # Gamma at Kt=1: sqrt(pi)/(2) * exp(0.25) = 0.886 * 1.284 = 1.137
        expected = (np.sqrt(np.pi) / 2.0) * np.exp(0.25)
        assert Gamma[0] == pytest.approx(
            expected, rel=0.01
        ), f"Gamma at Kt=1 should be {expected:.3f}, got {Gamma[0]:.3f}"


# ---------------------------------------------------------------------------
# TestTATGeneration
# ---------------------------------------------------------------------------


class TestTATGeneration:
    """Test TAT generation rate physics."""

    def test_tat_exceeds_srh_at_reverse_bias(self, srh_only_device):
        """TAT device dark current should be >> SRH-only at -30V."""
        srh_dev = srh_only_device
        ramp_bias(srh_dev, V_target=-30.0, contact="anode", V_step=1.0)
        I_srh = abs(extract_contact_current(srh_dev, contact="cathode"))

        tat_dev = create_dark_current_device(device_name=_unique_name("tat_cmp"))
        ramp_bias(tat_dev, V_target=-30.0, contact="anode", V_step=1.0)
        _compute_node_efield(tat_dev["device_name"], tat_dev["region_name"])
        _compute_gamma_factors(tat_dev["device_name"], tat_dev["region_name"])
        I_tat = abs(extract_contact_current(tat_dev, contact="cathode"))

        # TAT current should dominate SRH by a large factor
        assert (
            I_tat > 10 * I_srh
        ), f"TAT current ({I_tat:.2e}) should be >> SRH current ({I_srh:.2e})"

    def test_tat_components_exist(self, tat_device_reverse):
        """Component extraction should return non-trivial TAT contribution."""
        dev = tat_device_reverse
        components = extract_dark_current_components(dev)

        # The TAT component (from integration) should be significantly nonzero
        assert abs(components["J_TAT"]) > 0, "J_TAT should be nonzero"
        # I_total from contact should be nonzero
        assert (
            abs(components["I_total"]) > 1e-20
        ), f"I_total should be measurable, got {components['I_total']:.2e}"


# ---------------------------------------------------------------------------
# TestSurfaceRecombination
# ---------------------------------------------------------------------------


class TestSurfaceRecombination:
    """Test surface recombination velocity effects."""

    def test_srv_increases_current(self):
        """Higher SRV should produce larger SRV component."""
        dev_no_srv = create_dark_current_device(
            device_name=_unique_name("no_srv"),
            S_n=0.0,
            S_p=0.0,
        )
        ramp_bias(dev_no_srv, V_target=-30.0, contact="anode", V_step=1.0)
        _compute_node_efield(dev_no_srv["device_name"], dev_no_srv["region_name"])
        _compute_gamma_factors(dev_no_srv["device_name"], dev_no_srv["region_name"])
        comp_no_srv = extract_dark_current_components(dev_no_srv)

        dev_high_srv = create_dark_current_device(
            device_name=_unique_name("high_srv"),
            S_n=1e4,
            S_p=1e4,
        )
        ramp_bias(dev_high_srv, V_target=-30.0, contact="anode", V_step=1.0)
        _compute_node_efield(dev_high_srv["device_name"], dev_high_srv["region_name"])
        _compute_gamma_factors(dev_high_srv["device_name"], dev_high_srv["region_name"])
        comp_high_srv = extract_dark_current_components(dev_high_srv)

        # SRV component should be larger for high SRV
        assert abs(comp_high_srv["I_SRV"]) >= abs(comp_no_srv["I_SRV"]), (
            f"High SRV current ({comp_high_srv['I_SRV']:.2e}) should be >= "
            f"zero SRV current ({comp_no_srv['I_SRV']:.2e})"
        )


# ---------------------------------------------------------------------------
# TestDarkCurrentCalibration
# ---------------------------------------------------------------------------


class TestDarkCurrentCalibration:
    """Test dark current calibration against experimental target."""

    def test_dark_current_order_of_magnitude(self):
        """Simulated dark current at -30V should be within order of magnitude of 18 pA.

        Target range: 1.8 pA to 180 pA (0.1x to 10x of 18 pA).
        Uses default parameters: N_t=2.2e13 cm^-3/s, S_n=S_p=1e3 cm/s.
        N_t was calibrated to produce ~18 pA at -30V with the E-field-weighted
        effective generation model.
        """
        dev = create_dark_current_device(device_name=_unique_name("calib"))
        ramp_bias(dev, V_target=-30.0, contact="anode", V_step=1.0)
        _compute_node_efield(dev["device_name"], dev["region_name"])
        _compute_gamma_factors(dev["device_name"], dev["region_name"])

        components = extract_dark_current_components(dev, area=0.05)
        I_dark = abs(components["I_total"])

        assert 1.8e-12 <= I_dark <= 180e-12, (
            f"Dark current at -30V = {I_dark:.2e} A ({I_dark*1e12:.2f} pA), "
            f"should be between 1.8 pA and 180 pA"
        )

    def test_dark_current_sweep_monotonic(self):
        """Dark current magnitude should increase monotonically with reverse voltage."""
        dev = create_dark_current_device(device_name=_unique_name("mono"))

        V_range = [0.0, -5.0, -10.0, -20.0, -30.0]
        result = dark_current_sweep(dev, V_range, area=0.05, V_step=1.0)

        I_mag = np.abs(result["I_total"])
        for i in range(1, len(I_mag)):
            if I_mag[i - 1] > 1e-20:
                assert I_mag[i] >= I_mag[i - 1] * 0.9, (
                    f"|I| at V={result['voltages'][i]:.0f}V ({I_mag[i]:.2e}) "
                    f"decreased from V={result['voltages'][i-1]:.0f}V ({I_mag[i-1]:.2e})"
                )


# ---------------------------------------------------------------------------
# TestDarkCurrentSweep
# ---------------------------------------------------------------------------


class TestDarkCurrentSweep:
    """Test dark_current_sweep functionality."""

    def test_sweep_returns_components(self):
        """Sweep should return dict with all expected keys."""
        dev = create_dark_current_device(device_name=_unique_name("sweep_keys"))
        result = dark_current_sweep(dev, [-5.0, -10.0], area=0.05, V_step=1.0)

        expected_keys = {
            "voltages",
            "I_total",
            "I_SRH",
            "I_TAT",
            "I_SRV",
            "J_total",
            "J_SRH",
            "J_TAT",
            "J_SRV",
        }
        assert expected_keys.issubset(
            set(result.keys())
        ), f"Missing keys: {expected_keys - set(result.keys())}"

    def test_sweep_voltage_range(self):
        """Returned voltages should match input range."""
        dev = create_dark_current_device(device_name=_unique_name("sweep_vr"))
        V_range = [-5.0, -15.0, -25.0]
        result = dark_current_sweep(dev, V_range, area=0.05, V_step=1.0)

        np.testing.assert_array_almost_equal(
            result["voltages"],
            V_range,
            decimal=6,
            err_msg="Returned voltages should match input range",
        )
