"""Microdosimeter preset (Petringa 2025, Phys. Med. 140) and the tilted-alpha CCE helper.

Numbers from lab/microdosimeter (notes/LOG.md): doping (section 9), tilted generation,
CCE(V). All on the batched engine (~1 s total).
"""

import numpy as np
import pytest

import etna
from etna.api.device import DEFAULT_STACK, MICRODOSIMETER_N_D, DeviceConfig
from etna.api.simulation import _ENGINE_BOX, OutsideModelValidity
from etna.core.generation_profiles import (
    ALPHA_RANGE_SIC_CM,
    alpha_path_stopping,
    tilted_alpha_generation,
)


def test_preset_fields_and_default_stack_untouched():
    c = etna.microdosimeter_config()
    assert c.area_cm2 == 0.25 and c.epi_thickness_um == 10.0
    assert c.doping_profile == "uniform" and c.N_D == MICRODOSIMETER_N_D == 1.0e14
    assert c.substrate_thickness_um == 0.3 and c.nplus_thickness_um == 1.0
    assert etna.microdosimeter_config(1.0).substrate_thickness_um == 1.0
    assert etna.microdosimeter_config(T=310.0).T == 310.0
    # the paper's own 1/C^2 order of magnitude, inside the engine's validated box
    assert (
        _ENGINE_BOX["N_D_uniform"][0]
        <= MICRODOSIMETER_N_D
        <= _ENGINE_BOX["N_D_uniform"][1]
    )
    d = DeviceConfig()  # the dosimeter default is a different device
    # 4 mm2 dosimeter area (0.04 cm2) since 2026-10-18 (was 1e-4 cm2)
    assert d.area_cm2 == 0.04 and d.doping_profile == "graded"
    assert d.N_D_bulk == DEFAULT_STACK["N_D_bulk"]


@pytest.mark.parametrize("p_plus", [0.3, 1.0])
def test_preset_full_depletion_near_10V(p_plus):
    """Petringa 2025 JINST / SiC_Photons (same epi spec, 4 mm2): "full depletion at about
    10 V", W ~9.5 um at -10 V; the Phys. Med. sentence "approximately -30 V" (which the
    preset followed as 3.0575e14 until 2026-10-17) is the over-depleted operating point.
    """
    r = etna.run_cv(etna.microdosimeter_config(p_plus), 0.0, -60.0, 121)
    assert r.metadata["solver"] == "batched"
    assert r.metadata["reach_through_V"] == pytest.approx(9.0, abs=0.05)
    W = np.array(r.metadata["W_from_C"]) * 1e4
    assert W[20] > 9.5  # -10 V
    # 25 mm2 plateau, small-signal C = dQ/dV; 4 mm2-equivalent 3.44e-11 F vs the paper's
    # Fig. 3b 3.46e-11 (that figure is the 4 mm2 dosimeter's curve)
    assert r.y[-1] == pytest.approx(2.1485e-10, rel=2e-3)


def test_tilted_generation_geometry():
    x = np.linspace(0, 12e-4, 4001)
    G, lost = tilted_alpha_generation(x)
    pairs = np.trapezoid(G, x)
    assert pairs == pytest.approx(5.5e6 / 7.8, rel=1e-3)  # whole track in SiC
    assert lost == 0.0
    # ends at 18 cos65 = 7.61 um (straggling sigma 0.15 um projected)
    assert np.trapezoid(G[x > 8e-4], x[x > 8e-4]) < 1e-4 * pairs
    assert 6.5e-4 < x[np.argmax(G)] < 7.3e-4  # Bragg peak before the end (NIST: 6.99 um)
    G2, lost2 = tilted_alpha_generation(x, dead_layer_cm=0.3e-4)
    assert 0.0 < lost2 < 0.05
    assert np.trapezoid(G2, x) == pytest.approx(pairs * (1 - lost2), rel=2e-3)
    s = np.linspace(0, 23e-4, 20001)
    f = alpha_path_stopping(s)
    assert 2.4 < f.max() / f[0] < 2.6  # NIST ASTAR (SiC) peak/entrance: 2.49
    with pytest.raises(ValueError):
        tilted_alpha_generation(x, theta_deg=90.0)


def test_tilted_cce_normal_incidence_bragg_matches_run_cce():
    """theta=0, NIST Bragg curve at the run_cce range, no p+ generation ~ run_cce
    (differs by the 0.3 um p+ origin shift and the 5.5 vs 5.486 MeV curve)."""
    c = etna.microdosimeter_config()
    a = etna.run_cce(c, -10.0, -40.0, 7)
    b = etna.run_cce_tilted_alpha(
        c,
        -10.0,
        -40.0,
        7,
        theta_deg=0.0,
        alpha_range_um=ALPHA_RANGE_SIC_CM * 1e4,
        shape="bragg",
        generation_in_p_plus=False,
    )
    assert a.metadata["solver"] == b.metadata["solver"] == "batched"
    np.testing.assert_allclose(b.y, a.y, atol=2e-4)


def test_tilted_cce_vs_bias_pins_lab_numbers():
    r = etna.run_cce_tilted_alpha(etna.microdosimeter_config(), 0.0, -60.0, 13)
    V = np.round(r.x, 6)
    n = r.y / r.y[V == -50.0][0]
    at = lambda v: n[V == v][0]
    # default: epi-only generation (p+ = dead layer), NIST bragg shape (since
    # 2026-10-20; the qualitative curve gave 0.9945 / 0.9689 / 0.9845), no metal. N_D 1e14:
    # W(0 V) = 5.4 um, the depletion edge passes the track end (7.6 um) by -3 V, so the
    # steady-state DD curve is flat to 0.7 % from 0 V. Fig. 4b (0.908 at 0 V, +1 % from
    # -20 to -50 V) is NOT reproduced by the DD model at any doping (LOG.md section 9).
    # c-axis mobility default since 2026-10-21: 0.9920 at 0 V (legacy_perp 0.9928; the
    # holes diffusing from the undepleted tail are 10 % slower)
    assert at(0.0) == pytest.approx(0.9920, abs=5e-4)
    assert at(-10.0) == pytest.approx(0.9998, abs=3e-4)
    assert at(-20.0) > 0.999
    assert np.all(np.diff(r.y) >= -1e-6)
    assert r.y[V == -50.0][0] == pytest.approx(0.9999, abs=2e-4)
    # relative to E0/E_pair: the 2.7 % deposited in the 0.3 um p+ is lost
    assert r.metadata["cce_vs_E0"][V == -50.0][0] == pytest.approx(0.9731, abs=5e-4)
    assert r.metadata["projected_end_um"] == pytest.approx(7.607, abs=1e-3)
    assert r.metadata["reach_through_V"] == pytest.approx(9.0, abs=0.05)
    p = etna.run_cce_tilted_alpha(
        etna.microdosimeter_config(), -50.0, -50.0, 1, generation_in_p_plus=True
    )
    assert p.y[0] == pytest.approx(0.9863, abs=5e-4)  # p+ pairs: electron L_n ~ 0.3 um


def test_old_30V_fit_doping_is_still_simulable_but_not_the_preset():
    """The pre-2026-10-17 value (3.0575e14, V_fd 30 V) remains a valid override; its
    CCE(V) shape (0.967 at 0 V, flat from -12 V) is what LOG.md section 3 compared."""
    c = etna.microdosimeter_config(N_D=3.0575e14)
    r = etna.run_cce_tilted_alpha(c, 0.0, -50.0, 11)
    assert r.metadata["reach_through_V"] == pytest.approx(30.0, abs=0.05)
    n = r.y / r.y[-1]
    # 0.9668 with the qualitative curve; 0.9622 with the perpendicular-to-c mobility
    # (legacy_perp, before 2026-10-21)
    assert n[0] == pytest.approx(0.9583, abs=5e-4)
    assert n[2] > 0.999  # -10 V


def test_tilted_cce_refuses_outside_engine_domain():
    with pytest.raises(OutsideModelValidity):
        etna.run_cce_tilted_alpha(etna.microdosimeter_config(N_D=5e13), 0.0, -10.0, 3)
    with pytest.raises(OutsideModelValidity):
        etna.run_cce_tilted_alpha(etna.microdosimeter_config(), 1.0, -10.0, 3)
    with pytest.raises(NotImplementedError):
        etna.run_cce_tilted_alpha(etna.microdosimeter_config(half_width_um=50.0))


def test_run_cce_cold_start_at_minus_50V_stays_on_the_engine():
    """A cold 0 -> -50 V Poisson jump does not converge on low-doped uniform devices:
    run_cce used to fall back to devsim silently. The sweeps now walk 0 V in <= 10 V
    Poisson sub-steps (as run_cce_tilted_alpha)."""
    from etna.api import simulation as sim

    c = etna.microdosimeter_config()
    r = etna.run_cce(c, -50.0, -60.0, 3)
    assert r.metadata["solver"] == "batched"
    assert "fallback_reason" not in r.metadata
    assert r.metadata["reach_through_V"] == pytest.approx(9.0, abs=0.05)
    d = sim._run_cce_devsim(c, np.array([-50.0, -55.0, -60.0]))
    np.testing.assert_allclose(r.y, d.y, atol=1e-6)
    cv = etna.run_cv(c, -50.0, -60.0, 3)
    assert cv.metadata["solver"] == "batched"
