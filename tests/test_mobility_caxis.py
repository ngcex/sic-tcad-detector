"""c-axis low-field mobility default (2026-10-21, issue #6 follow-up).

The solvers' default set is "c_axis" (parallel to c, the vertical current of a (0001)
diode); "legacy_perp" (the perpendicular-to-c set, the previous default) stays
selectable through DeviceConfig(mobility=...) and reproduces the earlier numbers bit for
bit. Provenance: etna.core.physics_models.LOW_FIELD_MOBILITY; numbers:
lab/mobility_caxis/notes/LOG.md.
"""

import numpy as np
import pytest

from etna import DeviceConfig
from etna.api import simulation as sim
from etna.api.device import device_kwargs
from etna.api.models import run_transit_time
from etna.core.charge_collection import cce_vs_bias
from etna.core.physics_models import (
    HIGH_FIELD_MOBILITY,
    LOW_FIELD_MOBILITY,
    PhysicsModels,
)
from etna.core.sic_material import (
    DEFAULT_MOBILITY,
    SiC4H_Parameters,
    mobility_caughey_thomas,
    mobility_caughey_thomas_T,
)

L = LOW_FIELD_MOBILITY
N = np.geomspace(1e13, 1e20, 22)
LEGACY = DeviceConfig(mobility="legacy_perp")


# ------------------------------------------------------------------ the sets


def test_default_is_c_axis_everywhere():
    assert DEFAULT_MOBILITY == "c_axis"
    assert SiC4H_Parameters().mobility == "c_axis"
    assert DeviceConfig().mobility == "c_axis"
    assert L["c_axis"].direction == "c_axis"
    assert L["legacy_perp"].direction == "perp"
    assert "mobility" not in device_kwargs(DeviceConfig())  # default keyword set kept
    assert device_kwargs(LEGACY)["mobility"] == "legacy_perp"


@pytest.mark.parametrize("bad", ["etna_default", "C_axis", "", None])
def test_unknown_set_is_refused(bad):
    with pytest.raises(ValueError, match="legacy_perp"):
        DeviceConfig(mobility=bad)


@pytest.mark.parametrize("carrier", ["electron", "hole"])
def test_legacy_perp_is_the_previous_formula(carrier):
    p = SiC4H_Parameters(mobility="legacy_perp")
    s = "n" if carrier == "electron" else "p"
    mu_max, mu_min = getattr(p, f"mu_{s}_max"), getattr(p, f"mu_{s}_min")
    N_ref, alpha, gamma = (getattr(p, f"{k}_{s}") for k in ("N_ref", "alpha", "gamma"))
    for T in (150.0, 300.0, 450.0):
        old = mu_min + (mu_max * (T / 300.0) ** gamma - mu_min) / (
            1.0 + (N / N_ref) ** alpha
        )
        np.testing.assert_array_equal(mobility_caughey_thomas_T(N, T, carrier, p), old)
        np.testing.assert_array_equal(L["legacy_perp"].mobility(N, T, carrier), old)
    np.testing.assert_array_equal(
        mobility_caughey_thomas_T(N, 300.0, carrier, p),
        mobility_caughey_thomas(N, carrier),
    )


def test_c_axis_electrons_are_ishikawa_eq_2_16_at_300K():
    """At 300 K the c_axis electrons equal the thesis's drift fit (Eq. 2.16, the
    'arora' form of ishikawa2025) to round-off; away from 300 K the T laws differ by
    design (c_axis keeps Schaffer's <0001> -2.40 on mu_max)."""
    np.testing.assert_allclose(
        L["c_axis"].mobility(N, 300.0, "electron"),
        L["ishikawa2025"].mobility(N, 300.0, "electron"),
        rtol=1e-14,
    )
    # Eq. 2.16 at 300 K written out
    ref = (1190.0 - 10.0) / (1.0 + (N / 1.0e17) ** 0.54) + 10.0
    np.testing.assert_allclose(
        L["c_axis"].mobility(N, 300.0, "electron"), ref, rtol=1e-14
    )
    # the full thesis T law at 450 K, written out
    t = 1.5
    ref450 = (1180.0 * t**-3.4) / (1.0 + (N / 1e17) ** 0.54 * t**-2.2) + 10.0 * t**0.8
    np.testing.assert_allclose(
        L["ishikawa2025"].mobility(N, 450.0, "electron"), ref450, rtol=1e-14
    )


def test_c_axis_against_the_independent_cross_checks():
    """Opened primaries (lab/mobility_caxis/notes/LOG.md): vertical-SBD drift mobility
    1070 +- 290 at N_D 9.7e14 (Ishikawa 2024); the thesis's Hall fit parallel to c (Eq.
    2.11; its pulsed-conductance values of Table 3.1, 1190-1200 at 2e15-5e15, use the
    Hall density and agree with it) is the drift value times a Hall factor 1.05-1.30;
    the thesis's perpendicular DRIFT fit (Eq. 2.17) equals the legacy set (Schaffer
    Hall) within 3.5 % up to 5e15 (4.2 % at 1e16); mu_par / mu_perp 1.16 (Kimoto JJAP
    2019) to 1.2 (Ayalew eq. 3.81) at low doping."""
    mu = lambda n: mobility_caughey_thomas_T(n, 300.0, "electron")
    assert abs(mu(9.7e14) - 1070.0) < 290.0
    n = np.geomspace(1e14, 1e18, 13)
    hall_par = 20.0 + 1240.0 / (1.0 + (n / 2.0e17) ** 0.64)
    assert np.all((hall_par / mu(n) > 1.05) & (hall_par / mu(n) < 1.31))
    legacy = SiC4H_Parameters(mobility="legacy_perp")
    lo = n <= 5e15
    drift_perp = 30.0 + 950.0 / (1.0 + (n / 1.2e17) ** 0.55)
    leg = mobility_caughey_thomas_T(n, 300.0, "electron", legacy)
    np.testing.assert_allclose(drift_perp[lo], leg[lo], rtol=0.035)
    for nn in (1e14, 1e15):
        r = mu(nn) / mobility_caughey_thomas_T(nn, 300.0, "electron", legacy)
        assert 1.16 <= r <= 1.25
    # holes: unchanged in c_axis (the primaries disagree on them); x 0.90 (thesis drift
    # ratio) only in the opt-in c_axis_holes_0p9, whose electrons are c_axis's
    for T in (250.0, 300.0, 450.0):
        np.testing.assert_array_equal(
            L["c_axis"].mobility(N, T, "hole"), L["legacy_perp"].mobility(N, T, "hole")
        )
        np.testing.assert_allclose(
            L["c_axis_holes_0p9"].mobility(N, T, "hole"),
            0.9 * L["legacy_perp"].mobility(N, T, "hole"),
            rtol=1e-12,
        )
        np.testing.assert_array_equal(
            L["c_axis_holes_0p9"].mobility(N, T, "electron"),
            L["c_axis"].mobility(N, T, "electron"),
        )


@pytest.mark.filterwarnings("ignore:low-field mobility:UserWarning")
@pytest.mark.parametrize(
    "name", ["legacy_perp", "c_axis", "c_axis_holes_0p9", "ishikawa2025"]
)
@pytest.mark.parametrize("carrier", ["electron", "hole"])
def test_devsim_expression_is_the_numpy_evaluation(name, carrier):
    """The devsim string (device.mobility_ct_expression) evaluates to the numpy value
    (same constants, same order of operations)."""
    from etna.core.device import mobility_ct_expression

    p = SiC4H_Parameters(mobility=name)
    for T in (250.0, 300.0, 410.0):
        expr = mobility_ct_expression("NN", T, carrier, p)
        # eval of a string ETNA itself generates (floats, pow, NN): safe in this test
        assert set(expr) <= set("0123456789.e+-*/() ,NNpow")
        got = np.array([eval(expr, {"pow": pow, "NN": float(n)}) for n in N])
        np.testing.assert_allclose(
            got, mobility_caughey_thomas_T(N, T, carrier, p), rtol=2e-15
        )


# ---------------------------------------------------------- the solvers


def test_cv_does_not_depend_on_the_mobility():
    a = sim.run_cv(DeviceConfig(), 0.0, -40.0, 5)
    b = sim.run_cv(LEGACY, 0.0, -40.0, 5)
    np.testing.assert_array_equal(a.y, b.y)


def test_legacy_selection_reaches_engine_and_devsim():
    """DeviceConfig(mobility="legacy_perp") is honoured by the batched engine and by the
    devsim path (they agree), and differs from the default."""
    V = np.linspace(-2.0, -10.0, 3)
    res = {}
    for cfg in (DeviceConfig(), LEGACY):
        r = sim.run_cce(cfg, v_start=V[0], v_stop=V[-1], n_points=3)
        assert r.metadata["solver"] == "batched"
        kw = device_kwargs(cfg)
        ref = cce_vs_bias(
            V, epi_thickness_cm=kw.pop("epi_thickness_cm"), device_kwargs=kw
        )
        np.testing.assert_allclose(r.y, ref["cce_values"], atol=1e-6)
        res[cfg.mobility] = np.asarray(r.y)
    # only the electrons change: the CCE moves by < 1e-5 (the losses are hole
    # diffusion out of the undepleted tail, unchanged)
    d = res["c_axis"] - res["legacy_perp"]
    assert np.all((d > 0.0) & (d < 1e-5))


def test_transit_times_scale_with_the_mobility():
    a = run_transit_time(DeviceConfig(), v_start=-40.0, v_stop=-100.0, n_points=2)
    b = run_transit_time(LEGACY, v_start=-40.0, v_stop=-100.0, n_points=2)
    np.testing.assert_array_equal(a.metadata["t_h"], b.metadata["t_h"])
    r_e = np.asarray(a.metadata["t_e"]) / np.asarray(b.metadata["t_e"])
    assert np.all((r_e > 0.80) & (r_e < 0.86))


# ------------------------------------------------------- PhysicsModels


def test_anisotropy_only_rotates_a_perpendicular_set():
    with pytest.raises(ValueError, match="already parallel to c"):
        PhysicsModels(anisotropy="burin2024")
    with pytest.raises(ValueError, match="already parallel to c"):
        PhysicsModels(low_field_mobility="c_axis", anisotropy="ayalew2004")
    m = PhysicsModels(low_field_mobility="legacy_perp", anisotropy="burin2024")
    assert m.vertical_is_parallel()


def test_high_field_vsat_follows_the_transport_direction():
    hf = HIGH_FIELD_MOBILITY["lades2000"]
    m = PhysicsModels(high_field_mobility="lades2000")
    assert m.high_field_params(300.0)[0] == hf.vsat_par_n  # default device: c-axis
    assert m.high_field_params(300.0, "legacy_perp")[0] == hf.vsat_perp_n
    assert (
        m.describe("legacy_perp")["low_field_mobility"]["device_set"] == "legacy_perp"
    )
    ish = PhysicsModels(high_field_mobility="ishikawa2025")
    assert ish.high_field_params(300.0)[:2] == (1.7e7, 0.76)


def test_legacy_perp_reproduces_the_pre_change_pins():
    """The numbers pinned before 2026-10-21 (tests/test_default_stack.py history) come
    back bit for bit (to their printed digits) with mobility="legacy_perp"."""
    cce = sim.run_cce(LEGACY, -10.0, -50.0, 3)
    np.testing.assert_allclose(cce.y, [0.9994291, 0.9997970, 0.9998749], atol=1e-7)
    fl = sim.run_flash_recombination(LEGACY, dose_rates_Gy_s=[20.0])
    assert fl.y[0] == pytest.approx(0.999809, abs=1e-6)
    fl_new = sim.run_flash_recombination(DeviceConfig(), dose_rates_Gy_s=[20.0])
    assert abs(fl_new.y[0] - fl.y[0]) < 1e-5  # electrons only: within 1e-5


def test_ishikawa2025_warns_outside_its_fitted_T_range():
    """The thesis T laws were fitted at 300-600 K; below, the electron t^-3.4 grows fast
    (2159 cm2/Vs at 250 K vs 1819 with the kept -2.40): a warning, not an error. The
    array-T path (track engine) is the scalar expression, broadcast."""
    import warnings

    from etna.core.physics_models import MobilityExtrapolationWarning

    s = L["ishikawa2025"]
    with pytest.warns(MobilityExtrapolationWarning, match="300-600 K"):
        s.mobility(1e14, 250.0, "electron")
    with warnings.catch_warnings():
        warnings.simplefilter("error", MobilityExtrapolationWarning)
        L["c_axis"].mobility(1e14, 150.0, "electron")  # no stated range: no warning
        T = np.array([300.0, 350.0, 450.0, 600.0])
        for c in ("electron", "hole"):
            np.testing.assert_allclose(
                s.mobility(np.full(4, 1e16), T, c),
                [s.mobility(1e16, t, c) for t in T],
                rtol=1e-14,
            )
