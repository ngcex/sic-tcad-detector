"""Opt-in carrier-dependent incomplete ionization (audit-v7 A-02) in every solver.

``DeviceConfig(incomplete_ionization=True)``: N_D+(n) / N_A-(p) (Burin et al.
arXiv:2410.06798 Eq. 82, Kagamihara 2004 N, Achatz 2008 Al; etna.core.incomplete_ionization)
inside Poisson of devsim (1D / 2D), the batched 1D engine (C-V, CCE) and engine2d (dark
field, lateral CCE). Gates: engine vs devsim with the model ON at the existing parity
levels (C: devsim's +-1 mV central difference, <= 3.2e-8 as the _ENGINE_BOX T gate; W, psi,
CCE <= 1e-9), the model changes the physics (V_bi, C(0)), the legacy model is untouched,
and the paths whose engine does not carry the model run on devsim.
"""

import uuid

import devsim
import numpy as np
import pytest

from etna.api import simulation as sim
from etna.api.device import (
    DeviceConfig,
    device_kwargs,
    legacy_config,
    measured_cv_tabulated_config,
)
from etna.core.batched import make_batch, node_positions
from etna.core.batched.material import Q
from etna.core.batched.sweeps import _Ramp
from etna.core.device import apply_damaged_params, create_sic_device
from etna.core.drift_diffusion import create_dd_device
from etna.core.incomplete_ionization import (
    DEFAULT_IONIZATION_MODEL,
    IonizationModel,
    neutral_ionized_acceptor_concentration,
)

ON = DeviceConfig(incomplete_ionization=True)


def _values(info, name):
    return np.array(
        devsim.get_node_model_values(
            device=info["device_name"], region=info["region_name"], name=name
        )
    )


# ------------------------------------------------------------------- config plumbing
def test_config_default_is_legacy_and_emits_no_keyword():
    assert DeviceConfig().incomplete_ionization is False
    assert "incomplete_ionization" not in device_kwargs(DeviceConfig())
    assert device_kwargs(ON)["incomplete_ionization"] is DEFAULT_IONIZATION_MODEL
    custom = IonizationModel(g_A=2.0)
    kw = device_kwargs(DeviceConfig(incomplete_ionization=custom))
    assert kw["incomplete_ionization"] is custom
    with pytest.raises(ValueError, match="incomplete_ionization"):
        DeviceConfig(incomplete_ionization="yes")


def test_unwired_engines_route_to_devsim():
    from etna.core.charge_collection import _fluence_engine_reason
    from etna.core.cv_analysis import _cv_engine_reason
    from etna.core.dark_current import dark_engine_reason
    from etna.core.flash_recombination import flash_engine_reason

    kw = device_kwargs(ON)
    V = np.array([-10.0])
    for why in (
        dark_engine_reason(kw, V),
        flash_engine_reason(kw, -10.0),
        _fluence_engine_reason(kw, -10.0),
        _cv_engine_reason(kw, V),
    ):
        assert "incomplete_ionization" in why
    # the legacy keyword set is unaffected
    assert dark_engine_reason(device_kwargs(DeviceConfig()), V) is None
    # the C-V / CCE engine and engine2d take it
    assert sim._engine_unsupported_reason(ON, V) is None
    assert (
        sim._field_engine_unsupported_reason(
            DeviceConfig(half_width_um=25.0, incomplete_ionization=True), -20.0
        )
        is None
    )


def test_legacy_batch_is_bit_identical():
    sub, epi, t_np = 0.3e-4, 10e-4, 1e-4
    x = node_positions(epi, sub, t_np)
    args = (x, sub, 1e19, 8.2e15, 3.4e13, 0.46e-4)
    kw = dict(nplus=(sub + epi, 1e19), buffer=None)
    a = make_batch(*args, **kw)
    b = make_batch(*args, incomplete_ionization=False, **kw)
    assert a.ionization is None and b.ionization is None
    np.testing.assert_array_equal(a.netdop, b.netdop)
    ra = a.solve(a.neutral_guess(), V_cathode=10.0)
    c = make_batch(*args, incomplete_ionization=True, **kw)
    rc = c.solve(c.neutral_guess(), V_cathode=10.0)
    assert np.abs(ra.psi - rc.psi).max() > 1e-3  # the model moves the potential


# ------------------------------------------------------------------- devsim physics
def test_devsim_depleted_dopants_fully_ionized_neutral_ones_not():
    name = f"ii_dd_{uuid.uuid4().hex[:6]}"
    info = create_dd_device(device_name=name, **device_kwargs(ON))
    try:
        x, nd = _values(info, "x"), _values(info, "Donors")
        dp, am = _values(info, "IonizedDonors"), _values(info, "IonizedAcceptors")
        na_tot = info["N_A"]
        j = info["junction_pos"]
        # inside the 0 V depletion (W(0) ~ 1.6 um): donors fully ionized
        k = np.argmin(np.abs(x - (j + 0.5e-4)))
        assert dp[k] / nd[k] == pytest.approx(1.0, abs=1e-9)
        # neutral n+ (1e19): about half ionized (cubic site deeper); neutral p+ Al ~2.9 %
        k = np.argmin(np.abs(x - (info["epi_end"] + 0.5e-4)))
        assert 0.45 < dp[k] / nd[k] < 0.56
        k = np.argmin(np.abs(x - 0.1e-4))
        assert am[k] / na_tot == pytest.approx(0.02913, rel=2e-3)
        assert info["N_A_ionized"] == pytest.approx(
            neutral_ionized_acceptor_concentration(1e19, 300.0), rel=1e-15
        )
        psi = _values(info, "Potential")
        vbi_on = psi[-1] - psi[0]
    finally:
        devsim.delete_device(device=name)
    name = f"ii_dd_{uuid.uuid4().hex[:6]}"
    info = create_dd_device(device_name=name, **device_kwargs(DeviceConfig()))
    try:
        psi = _values(info, "Potential")
        vbi_off = psi[-1] - psi[0]
    finally:
        devsim.delete_device(device=name)
    # V_bi drops by ~56 mV (p+: kT ln(0.029/0.132) = -39 mV; n+ contact ~ -17 mV)
    assert vbi_off == pytest.approx(3.19832, abs=1e-4)
    assert vbi_on == pytest.approx(3.14189, abs=1e-4)


def test_damage_override_rebuilds_the_ionization_fields():
    name = f"ii_dmg_{uuid.uuid4().hex[:6]}"
    info = create_sic_device(device_name=name, **device_kwargs(ON))
    try:
        x, nd = _values(info, "x"), _values(info, "Donors")
        epi = x >= info["junction_pos"]
        apply_damaged_params(
            info, {"tau_n": 1e-9, "tau_p": 1e-9, "N_D_profile": 0.5 * nd[epi]}
        )
        nd2 = _values(info, "Donors")
        np.testing.assert_allclose(_values(info, "IonD0"), 0.5 * nd2, rtol=1e-15)
        np.testing.assert_allclose(_values(info, "IonD1"), 0.5 * nd2, rtol=1e-15)
        net = _values(info, "NetDoping")
        k = np.argmin(np.abs(x - (info["epi_end"] + 0.5e-4)))  # neutral n+, halved
        assert 0.5 < net[k] / nd2[k] < 0.75
    finally:
        devsim.delete_device(device=name)


# ------------------------------------------------- batched C: independent derivation
def test_engine_capacitance_includes_the_bound_hole_charge():
    """Analytic engine C == central difference of q*sum(v (p - N_A-)) (neutral Al = a
    bound hole), and differs from the free-hole-only difference: the term matters."""
    sub, epi, t_np = 0.3e-4, 10e-4, 1e-4
    x = node_positions(epi, sub, t_np)
    pb = make_batch(
        x, sub, 1e19, 8.2e15, 3.4e13, 0.46e-4,
        nplus=(sub + epi, 1e19), buffer=None, incomplete_ionization=True,
    )  # fmt: skip
    v = pb.v.copy()
    v[0] = v[-1] = 0.0
    ramp = _Ramp(pb, 1e-12)
    dv = 0.05
    hole = lambda r: r.p[0] - pb.ionized_doping(r)[1][0]
    for Vc in (0.5, 10.0, 30.0):
        r_lo = ramp.to(Vc - dv)
        r_mid = ramp.to(Vc)
        C = pb.small_signal_capacitance(r_mid)[0]
        r_hi = ramp.to(Vc + dv)
        dQ = Q * np.sum(v * (hole(r_lo) - hole(r_hi))) / (2 * dv)
        assert C == pytest.approx(dQ, rel=2e-5)
        dQ_free = Q * np.sum(v * (r_lo.p[0] - r_hi.p[0])) / (2 * dv)
        assert abs(dQ_free / dQ - 1) > 1e-3


# ------------------------------------------------------------- 1D engine vs devsim
PARITY_1D = {
    "default": (ON, -100.0, 21),
    "tabulated": (measured_cv_tabulated_config(incomplete_ionization=True), -100.0, 21),
    "T150": (DeviceConfig(incomplete_ionization=True, T=150.0), -100.0, 11),
    "T450": (DeviceConfig(incomplete_ionization=True, T=450.0), -100.0, 11),
    "legacy_no_nplus": (legacy_config(incomplete_ionization=True), -40.0, 9),
}


@pytest.mark.slow
@pytest.mark.parametrize("name", PARITY_1D)
def test_1d_engine_matches_devsim(name):
    cfg, v_stop, n = PARITY_1D[name]
    e = sim.run_cv(cfg, 0.0, v_stop, n)
    d = sim._run_cv_devsim(cfg, np.linspace(0.0, v_stop, n))
    assert e.metadata["solver"] == "batched" and d.metadata["solver"] == "devsim"
    assert np.max(np.abs(e.y / d.y - 1)) < 3.2e-8  # devsim C: +-1 mV central difference
    np.testing.assert_allclose(
        e.metadata["depletion_widths"],
        d.metadata["depletion_widths"],
        rtol=0,
        atol=1e-13,
    )
    Vc = np.linspace(-10.0, v_stop, 4)
    ce = sim.run_cce(cfg, Vc[0], Vc[-1], 4)
    cd = sim._run_cce_devsim(cfg, Vc)
    assert ce.metadata["solver"] == "batched"
    assert np.max(np.abs(ce.y - cd.y)) < 1e-9


@pytest.mark.slow
def test_model_moves_the_zero_bias_capacitance_not_the_plateau():
    off = sim.run_cv(DeviceConfig(), 0.0, -30.0, 4)
    on = sim.run_cv(ON, 0.0, -30.0, 4)
    rel = on.y / off.y - 1
    # 0 V: V_bi -56 mV and the frozen-out neutral-side Debye tails, C +2.2 %
    assert rel[0] == pytest.approx(0.022, abs=0.003)
    assert np.all(np.abs(rel[1:]) < 2e-3)  # depleted epi: fully ionized either way


# ---------------------------------------------------------------- engine2d vs devsim
@pytest.mark.slow
@pytest.mark.parametrize("bias", [-20.0, -50.0])
def test_engine2d_field_matches_devsim(monkeypatch, bias):
    cfg = DeviceConfig(half_width_um=25.0, incomplete_ionization=True)
    eng = sim.run_field(cfg, bias)
    monkeypatch.setattr(
        sim, "_field_engine_unsupported_reason", lambda c, v: "forced devsim"
    )
    ref = sim.run_field(cfg, bias)
    assert eng.metadata["solver"] == "engine2d" and ref.metadata["solver"] == "devsim"
    me, md = eng.mesh, ref.mesh
    np.testing.assert_array_equal(
        me.node_values["NetDoping"], md.node_values["NetDoping"]
    )
    assert (
        np.abs(me.node_values["Potential"] - md.node_values["Potential"]).max() < 1e-9
    )


@pytest.mark.slow
def test_engine2d_lateral_cce_matches_devsim():
    from etna.core.charge_collection_2d import cce_lateral_scan, create_2d_dd_device
    from etna.core.engine2d import lateral_cce_scan_2d

    kw = device_kwargs(ON)
    kw.pop("N_D")
    e = lateral_cce_scan_2d(25.0, 20.0, n_lateral_points=3, **kw)
    name = f"ii2d_{uuid.uuid4().hex[:8]}"
    info = create_2d_dd_device(half_width_um=25.0, V_bias=20.0, device_name=name, **kw)
    try:
        d = cce_lateral_scan(info, n_points=3)
    finally:
        devsim.delete_device(device=name)
    np.testing.assert_allclose(e["cce_values"], d["cce_values"], rtol=0, atol=1e-10)
