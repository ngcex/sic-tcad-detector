"""Microdosimeter preset (Petringa 2025, Phys. Med. 140) and the tilted-alpha CCE helper.

Numbers from lab/microdosimeter (notes/LOG.md): V_fd fit, tilted generation, CCE(V).
All on the batched engine (~1 s total).
"""

import numpy as np
import pytest

import etna
from etna.api.device import DEFAULT_STACK, MICRODOSIMETER_N_D, DeviceConfig
from etna.api.simulation import OutsideModelValidity
from etna.core.generation_profiles import alpha_path_stopping, tilted_alpha_generation


def test_preset_fields_and_default_stack_untouched():
    c = etna.microdosimeter_config()
    assert c.area_cm2 == 0.25 and c.epi_thickness_um == 10.0
    assert c.doping_profile == "uniform" and c.N_D == MICRODOSIMETER_N_D == 3.0575e14
    assert c.substrate_thickness_um == 0.3 and c.nplus_thickness_um == 1.0
    assert etna.microdosimeter_config(1.0).substrate_thickness_um == 1.0
    assert etna.microdosimeter_config(T=310.0).T == 310.0
    d = DeviceConfig()  # the dosimeter default is a different device
    assert d.area_cm2 == 1e-4 and d.doping_profile == "graded"
    assert d.N_D_bulk == DEFAULT_STACK["N_D_bulk"]


@pytest.mark.parametrize("p_plus", [0.3, 1.0])
def test_preset_full_depletion_is_30V(p_plus):
    r = etna.run_cv(etna.microdosimeter_config(p_plus), 0.0, -60.0, 121)
    assert r.metadata["solver"] == "batched"
    assert r.metadata["reach_through_V"] == pytest.approx(30.0, abs=0.05)
    # 25 mm2 plateau, small-signal C = dQ/dV (2026-10-17; it was eps*A/W_edge =
    # 2.160e-10, -0.50 %); 4 mm2-equivalent 3.44e-11 F vs the paper's Fig. 3b 3.46e-11
    assert r.y[-1] == pytest.approx(2.1493e-10, rel=2e-3)


def test_tilted_generation_geometry():
    x = np.linspace(0, 12e-4, 4001)
    G, lost = tilted_alpha_generation(x)
    pairs = np.trapezoid(G, x)
    assert pairs == pytest.approx(5.5e6 / 7.8, rel=1e-3)  # whole track in SiC
    assert lost == 0.0
    # ends at 18 cos65 = 7.61 um (straggling sigma 0.15 um projected)
    assert np.trapezoid(G[x > 8e-4], x[x > 8e-4]) < 1e-4 * pairs
    assert 6.0e-4 < x[np.argmax(G)] < 7.0e-4  # Bragg peak before the end
    G2, lost2 = tilted_alpha_generation(x, dead_layer_cm=0.3e-4)
    assert 0.0 < lost2 < 0.05
    assert np.trapezoid(G2, x) == pytest.approx(pairs * (1 - lost2), rel=2e-3)
    s = np.linspace(0, 23e-4, 20001)
    f = alpha_path_stopping(s)
    assert 1.9 < f.max() / f[0] < 2.3  # qualitative peak/entrance after straggling
    with pytest.raises(ValueError):
        tilted_alpha_generation(x, theta_deg=90.0)


def test_tilted_cce_normal_incidence_box_matches_run_cce():
    """theta=0, box, 15 um, no p+ generation ~ run_cce (differs only by the 0.3 um p+
    origin shift of the profile)."""
    c = etna.microdosimeter_config()
    a = etna.run_cce(c, -10.0, -40.0, 7)
    b = etna.run_cce_tilted_alpha(
        c,
        -10.0,
        -40.0,
        7,
        theta_deg=0.0,
        alpha_range_um=15.0,
        shape="box",
        generation_in_p_plus=False,
    )
    assert a.metadata["solver"] == b.metadata["solver"] == "batched"
    np.testing.assert_allclose(b.y, a.y, atol=2e-4)


def test_tilted_cce_vs_bias_pins_lab_numbers():
    r = etna.run_cce_tilted_alpha(etna.microdosimeter_config(), 0.0, -60.0, 13)
    V = np.round(r.x, 6)
    n = r.y / r.y[V == -50.0][0]
    at = lambda v: n[V == v][0]
    # default: epi-only generation (p+ = dead layer), bragg shape, no metal
    # hole mobility at the epi doping since 2026-10-17 (was N_A, audit-v7 A-03): 0.9488 -> 0.9668 (moves AWAY from the digitised 0.9083; not tuned)
    assert at(0.0) == pytest.approx(0.9668, abs=5e-4)
    assert at(-10.0) == pytest.approx(0.9993, abs=3e-4)
    # documented model limitation: flat from -20 V (data: 0.9904 -> 1.000 to -50 V)
    assert at(-20.0) > 0.999
    assert np.all(np.diff(r.y) >= -1e-6)
    assert r.y[V == -50.0][0] == pytest.approx(0.9999, abs=2e-4)
    # relative to E0/E_pair: the 3.1 % deposited in the 0.3 um p+ is lost
    assert r.metadata["cce_vs_E0"][V == -50.0][0] == pytest.approx(0.9688, abs=5e-4)
    assert r.metadata["projected_end_um"] == pytest.approx(7.607, abs=1e-3)
    assert r.metadata["reach_through_V"] == pytest.approx(30.0, abs=0.05)
    p = etna.run_cce_tilted_alpha(
        etna.microdosimeter_config(), -50.0, -50.0, 1, generation_in_p_plus=True
    )
    assert p.y[0] == pytest.approx(0.9845, abs=5e-4)  # p+ pairs: electron L_n ~ 0.3 um


def test_tilted_cce_refuses_outside_engine_domain():
    with pytest.raises(OutsideModelValidity):
        etna.run_cce_tilted_alpha(etna.microdosimeter_config(N_D=5e13), 0.0, -10.0, 3)
    with pytest.raises(OutsideModelValidity):
        etna.run_cce_tilted_alpha(etna.microdosimeter_config(), 1.0, -10.0, 3)
    with pytest.raises(NotImplementedError):
        etna.run_cce_tilted_alpha(etna.microdosimeter_config(half_width_um=50.0))
