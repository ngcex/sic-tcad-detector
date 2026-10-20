"""Engine-seeded bias (drift_diffusion.seeded_bias) == the bias ramp it replaces.

State: the seeded devsim state equals the ramped one to <= 1e-12 (potential in V,
carriers normwise relative); downstream quantities (1D CCE, 2D CCE, 2D devsim dark
current -- its TAT model is added after the seeded bias) to <= 1e-10 relative. Every
parity test also asserts the seed was actually used (a fallback to the ramp would pass
any parity trivially). Devices with physics the seed engine does not carry (Auger,
TAT / Hurkx, carrier-dependent ionization) and unsupported biases must be ramped: those
cases assert the reason and the bit-for-bit ramp state. Auger and TAT / Hurkx are
seeded where the call site allows them (seeded_bias(allow=...), the strict solve
carries the model): the dark Auger device of transient_cce_vs_dose_rate and the devsim
dark-current sweep (dark_current._bias_anode), with the same state / downstream gates.
"""

import uuid

import devsim
import numpy as np
import pytest

from etna.core import drift_diffusion as dd
from etna.core.device import LEGACY_STACK, MEASURED_CV_TABULATED_STACK
from etna.core.devsim_reset import delete_device_and_mesh

STATE_TOL = 1e-12
DOWNSTREAM_TOL = 1e-10
_FIT = ("N_t", "E_t")  # dark-current fit keys, not device keywords
LEGACY = {k: v for k, v in LEGACY_STACK.items() if k not in _FIT}
TABLE = {k: v for k, v in MEASURED_CV_TABULATED_STACK.items() if k not in _FIT}


def _state(info):
    return {
        m: np.array(
            devsim.get_node_model_values(
                device=info["device_name"], region=info["region_name"], name=m
            )
        )
        for m in ("Potential", "Electrons", "Holes")
    }


def _assert_same_state(got, ref, tol=STATE_TOL):
    assert np.max(np.abs(got["Potential"] - ref["Potential"])) <= tol
    for m in ("Electrons", "Holes"):
        scale = np.max(np.abs(ref[m]))
        assert np.max(np.abs(got[m] - ref[m])) <= tol * scale, m


def _build_1d(auger=False, **kw):
    info = dd.create_dd_device(device_name=f"seed1d_{uuid.uuid4().hex[:8]}", **kw)
    if auger:
        from etna.core.flash_recombination import add_auger_recombination

        add_auger_recombination(info)
    return info


def _build_tat():
    from etna.api.device import DeviceConfig, device_kwargs
    from etna.core.dark_current import create_dark_current_device

    cfg = DeviceConfig()
    info = create_dark_current_device(N_t=cfg.N_t, E_t=cfg.E_t, **device_kwargs(cfg))
    return info


def _build_2d(**kw):
    from etna.core.device2d import create_sic_2d_device
    from etna.core.poisson import setup_poisson, solve_equilibrium

    info = create_sic_2d_device(
        device_name=f"seed2d_{uuid.uuid4().hex[:8]}", half_width_um=5.0, **kw
    )
    setup_poisson(info)
    solve_equilibrium(info)
    dd.setup_sic_drift_diffusion(info)
    info["dd_initialized"] = True
    return info


def _ramped_and_seeded(build, V, contact, ramp_kw):
    info = build()
    try:
        dd.ramp_bias(info, V, contact, **ramp_kw)
        assert not info.get("ramp_relaxed_retries") and not info.get(
            "ramp_aggressive_retries"
        )
        ref = _state(info)
    finally:
        delete_device_and_mesh(info["device_name"])
    info = build()
    try:
        dd.seeded_bias(info, V, contact, **ramp_kw)
        got = _state(info)
        bias = devsim.get_parameter(
            device=info["device_name"],
            name=dd.simple_physics.GetContactBiasName(contact),
        )
    finally:
        delete_device_and_mesh(info["device_name"])
    return info, got, ref, bias


ADAPT = dict(V_step=0.5, V_step_max=2.0)


TWO_D = dict(V_step=0.5, V_step_max=5.0)


@pytest.mark.parametrize(
    "build, V, contact, ramp_kw",
    [
        (_build_1d, 50.0, "cathode", ADAPT),
        (_build_1d, -60.0, "anode", dict(V_step=0.5)),
        (lambda: _build_1d(T=400.0), -30.0, "anode", dict(V_step=0.5)),
        (lambda: _build_1d(**LEGACY), 40.0, "cathode", ADAPT),
        (lambda: _build_1d(**TABLE), 40.0, "cathode", ADAPT),
        (_build_2d, 30.0, "cathode", TWO_D),
        (lambda: _build_2d(**LEGACY), 30.0, "cathode", TWO_D),
    ],
    ids=["1d-cathode", "1d-anode", "1d-400K", "1d-legacy", "1d-table", "2d", "2d-legacy"],
)
def test_seeded_state_equals_ramped_state(build, V, contact, ramp_kw):
    info, got, ref, bias = _ramped_and_seeded(build, V, contact, ramp_kw)
    assert info["bias_seeded"] is True, info["bias_seed_reason"]
    assert bias == V
    _assert_same_state(got, ref)


@pytest.mark.parametrize(
    "build, V, contact, ramp_kw, reason",
    [
        (_build_1d, -1.0, "cathode", ADAPT, "forward bias"),
        (lambda: _build_1d(auger=True), -10.0, "anode", dict(V_step=0.5), "Auger"),
        (_build_tat, -10.0, "anode", dict(V_step=0.5), "TAT / Hurkx"),
        (
            lambda: _build_1d(incomplete_ionization=True),
            10.0,
            "cathode",
            ADAPT,
            "incomplete) ionization",
        ),
        (lambda: _build_2d(**LEGACY), 60.0, "cathode", TWO_D, "past the devsim wall"),
    ],
    ids=["forward", "auger", "tat-hurkx", "incomplete-ionization", "2d-legacy-wall"],
)
def test_ineligible_devices_are_ramped(build, V, contact, ramp_kw, reason):
    """Physics the seed engine does not carry (or a bias it does not represent): the
    ramp, bit for bit -- never a seed from an inconsistent model."""
    info, got, ref, bias = _ramped_and_seeded(build, V, contact, ramp_kw)
    assert info["bias_seeded"] is False
    assert reason in info["bias_seed_reason"]
    assert bias == V
    _assert_same_state(got, ref, tol=0.0)


def test_generation_present_is_ramped():
    from etna.core.charge_collection import add_generation_to_dd

    info = _build_1d()
    try:
        x = np.array(
            devsim.get_node_model_values(
                device=info["device_name"], region=info["region_name"], name="x"
            )
        )
        add_generation_to_dd(info, np.full_like(x, 1e10))
        dd.seeded_bias(info, -5.0, "anode")
        assert info["bias_seeded"] is False
        assert "generation present" in info["bias_seed_reason"]
    finally:
        delete_device_and_mesh(info["device_name"])


def test_failed_seed_solve_restores_and_ramps(monkeypatch):
    """A seed whose strict solve fails must leave the old ramp result: solution and
    bias restored before the ramp (else ramp_bias would see delta 0 and keep it)."""
    V = 30.0
    info = _build_1d()
    try:
        dd.ramp_bias(info, V, "cathode", **ADAPT)
        ref = _state(info)
    finally:
        delete_device_and_mesh(info["device_name"])

    real = dd.solve_tier
    calls = []

    def first_fails(tier):
        calls.append(tier)
        if len(calls) == 1:  # the seed's strict solve: leave a wrecked state behind
            d, r = info["device_name"], info["region_name"]
            psi = np.array(devsim.get_node_model_values(device=d, region=r, name="Potential"))
            devsim.set_node_values(device=d, region=r, name="Potential", values=(psi + 1.0).tolist())
            return None
        return real(tier)

    monkeypatch.setattr(dd, "solve_tier", first_fails)
    info = _build_1d()
    try:
        dd.seeded_bias(info, V, "cathode", **ADAPT)
        assert info["bias_seeded"] is False
        assert "strict solve" in info["bias_seed_reason"]
        assert len(calls) > 2  # the ramp ran after the failed seed
        _assert_same_state(_state(info), ref, tol=0.0)
    finally:
        delete_device_and_mesh(info["device_name"])


def test_seed_switch_off_ramps(monkeypatch):
    monkeypatch.setattr(dd, "SEED_BIAS", False)
    info = _build_1d()
    try:
        dd.seeded_bias(info, 10.0, "cathode", **ADAPT)
        assert info["bias_seeded"] is False
        assert info["bias_seed_reason"] == "seed disabled"
    finally:
        delete_device_and_mesh(info["device_name"])


# ----------------------------------------------------------------- downstream
def _seed_log(monkeypatch):
    """Record every seed attempt's outcome (None = seeded) of the run."""
    log = []
    real = dd._seed_bias

    def spy(*a, **kw):
        out = real(*a, **kw)
        log.append(out)
        return out

    monkeypatch.setattr(dd, "_seed_bias", spy)
    return log


def _both(monkeypatch, fn):
    """fn() with the ramp (SEED_BIAS False), then seeded; asserts every seed attempt of
    the seeded run succeeded."""
    monkeypatch.setattr(dd, "SEED_BIAS", False)
    ref = np.asarray(fn(), dtype=float)
    monkeypatch.setattr(dd, "SEED_BIAS", True)
    log = _seed_log(monkeypatch)
    got = np.asarray(fn(), dtype=float)
    assert log and all(r is None for r in log), log
    return got, ref


def _assert_close(got, ref, tol=DOWNSTREAM_TOL):
    np.testing.assert_allclose(got, ref, rtol=tol, atol=0)


def test_cce_vs_bias_parity(monkeypatch):
    from etna.core.charge_collection import cce_vs_bias

    def run():
        r = cce_vs_bias(np.array([-10.0, -40.0]))
        return np.concatenate([r["cce_values"], r["I_collected"]])

    _assert_close(*_both(monkeypatch, run))


def test_2d_cce_and_dark_parity(monkeypatch):
    from etna.core.charge_collection_2d import cce_lateral_scan, create_2d_dd_device
    from etna.core.optimization import get_dark_current_2d

    def run():
        info = create_2d_dd_device(half_width_um=5.0, V_bias=30.0)
        try:
            cce = cce_lateral_scan(info, n_points=2)["cce_values"]
        finally:
            delete_device_and_mesh(info["device_name"])
        dark = get_dark_current_2d(half_width_um=5.0, V_bias=30.0, solver="devsim")
        return np.concatenate([np.ravel(cce), [dark]])

    _assert_close(*_both(monkeypatch, run))


# ------------------------------------------- Auger / TAT allowed in the strict solve
def _tat_state_at(V, monkeypatch, seed, **kw):
    """TAT / Hurkx device at anode bias V via dark_current_at_bias (frozen-field fixed
    point), seeded or ramped; returns (state, components, bias_seeded)."""
    from etna.api.device import DeviceConfig, device_kwargs
    from etna.core.dark_current import create_dark_current_device, dark_current_at_bias

    monkeypatch.setattr(dd, "SEED_BIAS", seed)
    cfg = DeviceConfig()
    dev_kw = {**device_kwargs(cfg), **kw}
    info = create_dark_current_device(N_t=cfg.N_t, E_t=cfg.E_t, **dev_kw)
    try:
        comp = dark_current_at_bias(info, V)
        return _state(info), comp, info["bias_seeded"], info["bias_seed_reason"]
    finally:
        delete_device_and_mesh(info["device_name"])


def test_auger_allowed_seed_equals_ramp():
    """transient_cce_vs_dose_rate's dark Auger device: seed + strict solve WITH Auger
    == the 0.5 V ramp with Auger. (Auger moves the minority carriers of the dark state
    pointwise by up to ~8 % in the n+, so the strict solve must carry it.)"""
    V = -30.0
    info = _build_1d(auger=True)
    try:
        dd.ramp_bias(info, V, "anode")
        assert not info.get("ramp_relaxed_retries")
        assert not info.get("ramp_aggressive_retries")
        ref = _state(info)
    finally:
        delete_device_and_mesh(info["device_name"])
    info = _build_1d(auger=True)
    try:
        dd.seeded_bias(info, V, "anode", allow=("UAuger",))
        assert info["bias_seeded"] is True, info["bias_seed_reason"]
        got = _state(info)
    finally:
        delete_device_and_mesh(info["device_name"])
    _assert_same_state(got, ref)


@pytest.mark.parametrize("V", [-10.0, -60.0])
def test_tat_seeded_state_equals_ramped_state(V, monkeypatch):
    got, cg, seeded, reason = _tat_state_at(V, monkeypatch, True)
    ref, cr, ref_seeded, _ = _tat_state_at(V, monkeypatch, False)
    assert seeded is True, reason
    assert ref_seeded is False
    _assert_same_state(got, ref)
    _assert_close(cg["I_total"], cr["I_total"])


def test_tat_legacy_past_the_wall_is_ramped(monkeypatch):
    """Contact-on-epi cathode past FIELD_ENGINE_MAX: the dark sweep keeps its own ramp
    (and its convergence verdict), bit for bit."""
    got, cg, seeded, reason = _tat_state_at(-60.0, monkeypatch, True, **LEGACY)
    ref, cr, _, _ = _tat_state_at(-60.0, monkeypatch, False, **LEGACY)
    assert seeded is False and "past the devsim wall" in reason
    _assert_same_state(got, ref, tol=0.0)
    assert cg["I_total"] == cr["I_total"]


def test_tat_failed_seed_solve_restores_and_ramps(monkeypatch):
    """Seed strict solve fails on the TAT device: restored, then the dark sweep's own
    field-refreshing ramp (_ramp_anode) -- the ramped result bit for bit."""
    import etna.core.dark_current as dc

    ref, cr, _, _ = _tat_state_at(-20.0, monkeypatch, False)
    monkeypatch.setattr(dd, "solve_tier", lambda tier: None)
    ramps = []
    real = dc._ramp_anode
    monkeypatch.setattr(
        dc, "_ramp_anode", lambda *a, **k: (ramps.append(a[1]), real(*a, **k))
    )
    got, cg, seeded, reason = _tat_state_at(-20.0, monkeypatch, True)
    assert seeded is False and "strict solve" in reason
    assert ramps == [-20.0]
    _assert_same_state(got, ref, tol=0.0)
    assert cg["I_total"] == cr["I_total"]


def test_transient_dose_rate_parity(monkeypatch):
    from etna.core.transient import transient_cce_vs_dose_rate

    seen = []

    def run():
        df = transient_cce_vs_dose_rate(dose_rates=[20.0, 230.0])
        seen.append(df.attrs["bias_seeded"])
        return df["transient_cce_raw"].to_numpy()

    _assert_close(*_both(monkeypatch, run))
    assert seen == [False, True]


def test_run_dark_current_devsim_parity(monkeypatch):
    from etna.api.device import DeviceConfig
    from etna.api.simulation import run_dark_current

    seen = []

    def run():
        r = run_dark_current(DeviceConfig(), n_points=6, solver="devsim")
        assert not r.metadata["truncated"]
        seen.append(r.metadata["bias_seeded"])
        return r.y

    _assert_close(*_both(monkeypatch, run))
    assert not seen[0].any()
    # every bias point but the first (0 V: no bias change) was seeded
    assert seen[1].tolist() == [False] + [True] * 5
