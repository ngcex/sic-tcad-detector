"""Steady-state FLASH CCE (cce_vs_dose_rate) on the batched engine vs devsim, and the
engine's Auger n^3 term vs devsim's Auger model.

Facade parity (gate 1e-9) on the default stack, 150/450 K, legacy -30 V and notebook
05b's grid corners, with Auger active on both paths since 2026-10-20 (before, the devsim
path's add_generation_to_dd overwrote the Auger-augmented generation models and the engine
reproduced that with C_n = C_p = 0). test_auger_changes_cce fails on that old behaviour:
the "with Auger" CCE was bit-equal to the SRH-only reference.
"""

import itertools

import devsim
import numpy as np
import pytest

import etna.api.simulation as sim
from etna.api.device import (
    DeviceConfig,
    device_kwargs,
    legacy_config,
    microdosimeter_config,
)
from etna.core.batched.dark import generation_cce
from etna.core.charge_collection import add_generation_to_dd, compute_cce_from_dd
from etna.core.dark_current import dark_reference
from etna.core.device import DEFAULT_STACK
from etna.core.devsim_reset import delete_device_and_mesh
from etna.core.drift_diffusion import create_dd_device, ramp_bias
from etna.core.flash_recombination import (
    add_auger_recombination,
    cce_vs_dose_rate,
    flash_engine_reason,
)
from etna.core.generation_profiles import (
    confine_generation_to_epi,
    proton_generation_profile,
)
from etna.core.sic_material import srh_lifetime

GATE = 1e-9
DR = np.array([20.0, 50.0, 100.0, 150.0, 200.0, 230.0])


def _scaled(epi, nd, vb):
    nj = DEFAULT_STACK["N_D_junction"] * nd / DEFAULT_STACK["N_D_bulk"]
    return dict(epi_thickness_cm=epi, N_D_bulk=nd, N_D_junction=nj, V_bias=vb)


def _legacy(fn=legacy_config):
    kw = device_kwargs(fn())
    return dict(
        epi_thickness_cm=kw.pop("epi_thickness_cm"), device_kwargs=kw, V_bias=-30.0
    )


CASES = [
    dict(),
    dict(device_kwargs=dict(T=150.0)),
    dict(device_kwargs=dict(T=450.0)),
    _legacy(),
    _legacy(microdosimeter_config),
] + [
    _scaled(*c) for c in itertools.product([5e-4, 20e-4], [5e13, 5e14], [-10.0, -50.0])
]


@pytest.mark.slow
@pytest.mark.parametrize("case", CASES)
def test_engine_matches_devsim(case):
    e = cce_vs_dose_rate(DR, solver="engine", **case)
    d = cce_vs_dose_rate(DR, solver="devsim", **case)
    assert e["solver"] == "batched" and d["solver"] == "devsim"
    assert np.abs(e["cce_values"] - d["cce_values"]).max() <= GATE
    # the (numerical) dose-rate dependence itself is reproduced
    assert np.abs(np.diff(e["cce_values"]) - np.diff(d["cce_values"])).max() <= GATE
    assert abs(e["cce_no_auger_ref"] - d["cce_no_auger_ref"]) <= GATE


# notebook 05b corner where Auger matters most: 5 um epi, N_D_bulk 5e14, -10 V
# (minority-hole recombination in the 1e19 n+; measured offset 5.604e-5 on 2026-10-20)
AUGER_CORNER = _scaled(5e-4, 5e14, -10.0)


@pytest.mark.parametrize("solver", ["engine", "devsim"])
def test_auger_changes_cce(solver):
    """Regression (2026-10-20): Auger is genuinely in the swept device. The old code gave
    cce_values[0] == cce_no_auger_ref bit for bit (Auger dropped by add_generation_to_dd
    on devsim, C_n = C_p = 0 on the engine)."""
    r = cce_vs_dose_rate([20.0, 230.0], solver=solver, **AUGER_CORNER)
    assert r["solver"] == ("batched" if solver == "engine" else "devsim")
    assert (r["C_n"], r["C_p"]) == (5e-31, 2e-31)
    # the old behaviour fails here: the difference was exactly 0
    assert r["cce_no_auger_ref"] - r["cce_values"][0] > 4e-5
    # per dose rate vs the same sweep without Auger: a static ~5.6e-5 offset
    srh = cce_vs_dose_rate([20.0, 230.0], solver=solver, C_n=0.0, C_p=0.0, **AUGER_CORNER)
    drop = srh["cce_values"] - r["cce_values"]
    assert np.all(drop > 4e-5) and np.all(drop < 7e-5), drop
    assert abs(drop[1] - drop[0]) < 1e-7  # Auger adds no dose-rate dependence (~1e-8)


@pytest.mark.parametrize("solver", ["engine", "devsim"])
def test_auger_off_reproduces_srh_reference(solver):
    """C_n = C_p = 0: the sweep equals the SRH-only reference at the lowest rate."""
    r = cce_vs_dose_rate([20.0], solver=solver, C_n=0.0, C_p=0.0, **AUGER_CORNER)
    assert abs(r["cce_values"][0] - r["cce_no_auger_ref"]) <= GATE


def test_set_generation_keeps_auger_active():
    """devsim level: after set_generation the generation models still carry UAuger,
    after a bare add_generation_to_dd they do not (the original bug)."""
    from etna.core.flash_recombination import set_generation

    kw = dict(doping_profile="graded", L_transition=DEFAULT_STACK["L_transition"])
    c = dict(AUGER_CORNER)
    c.pop("V_bias")
    kw.update(c)
    name = "flash_setgen_ref"
    info = create_dd_device(device_name=name, **kw)
    region = info["region_name"]
    try:
        add_auger_recombination(info)
        get = lambda m: np.array(
            devsim.get_node_model_values(device=name, region=region, name=m)
        )
        # out of equilibrium: 1e10 cm^-3 excess carriers everywhere, so UAuger != 0
        for var in ("Electrons", "Holes"):
            devsim.set_node_values(
                device=name, region=region, name=var, values=list(get(var) + 1e10)
            )
        q = devsim.get_parameter(device=name, region=region, name="ElectronCharge")
        g = np.zeros(len(get("x")))
        set_generation(info, g)
        u_a = get("UAuger")
        assert u_a.max() > 1e17  # 1e19 n+: C_n N_D^2 * 1e10 = 5e17 cm^-3 s^-1
        with_auger = -q * (get("USRH") + u_a)
        np.testing.assert_allclose(get("ElectronGeneration"), with_auger, rtol=1e-12)
        add_generation_to_dd(info, g)  # bare: SRH-only again
        np.testing.assert_allclose(get("ElectronGeneration"), -q * get("USRH"), rtol=1e-12)
    finally:
        delete_device_and_mesh(name)


def test_engine_auger_term_matches_devsim_with_auger_active():
    """The n^3 term: devsim with Auger really in the equations (generation first, then
    add_auger_recombination, which keeps RadGenRate) vs the engine with C_n/C_p."""
    kw = dict(doping_profile="graded", L_transition=DEFAULT_STACK["L_transition"])
    c = _scaled(5e-4, 5e14, -10.0)
    V = c.pop("V_bias")
    kw.update(c)
    info = create_dd_device(device_name="flash_auger_ref", **kw)
    try:
        ramp_bias(info, -V, contact="cathode", V_step=0.5, V_step_max=2.0)
        x = np.array(
            devsim.get_node_model_values(
                device="flash_auger_ref", region=info["region_name"], name="x"
            )
        )
        g = proton_generation_profile(x, 62, 20.0)
        g[x < info["junction_pos"]] = 0.0
        g = confine_generation_to_epi(info, x, g)
        add_generation_to_dd(info, g)
        add_auger_recombination(info)
        devsim.solve(
            type="dc", absolute_error=1e10, relative_error=1e-10, maximum_iterations=60
        )
        ref_auger = compute_cce_from_dd(info, g)
    finally:
        delete_device_and_mesh("flash_auger_ref")
    r = dark_reference(kw)
    i, p = r["info"], r["info"]["params"]
    args = (r["x"], r["donors"], r["acceptors"], i["T"], i["N_A"], i["n_i"], V, g[None])
    tau = dict(
        tau_n=srh_lifetime(i["T"], "electron", p), tau_p=srh_lifetime(i["T"], "hole", p)
    )
    with_auger, ok = generation_cce(*args, C_n=p.C_n, C_p=p.C_p, **tau)
    without, ok0 = generation_cce(*args, **tau)
    assert ok.all() and ok0.all()
    assert abs(with_auger[0] - ref_auger) <= GATE
    # discriminating: Auger moves the CCE far more than the gate at this corner
    assert without[0] - with_auger[0] > 1e-5


def test_facade_provenance():
    r = sim.run_flash_recombination(DeviceConfig())
    md = r.metadata
    assert md["solver"] == "batched" and md["validated_domain"] is True
    assert md["cross_checked_to_V"] == 50.0 and "fallback_reason" not in md
    d = sim.run_flash_recombination(DeviceConfig(), solver="devsim")
    assert d.metadata["solver"] == "devsim"
    np.testing.assert_allclose(r.y, d.y, rtol=0, atol=GATE)


def test_scope():
    kw = device_kwargs(legacy_config())
    assert flash_engine_reason(kw, -40.0) is None
    assert "legacy" in flash_engine_reason(kw, -50.0)
    assert "T=500" in flash_engine_reason(device_kwargs(DeviceConfig(T=500.0)), -30.0)
    with pytest.raises(ValueError, match="engine unavailable"):
        sim.run_flash_recombination(DeviceConfig(T=500.0), solver="engine")


def test_engine_path_forwards_device_kwargs(monkeypatch):
    import etna.core.dark_current as dc

    got = {}
    real = dc.dark_reference

    def spy(kw):
        got.update(kw)
        return real(kw)

    monkeypatch.setattr(dc, "dark_reference", spy)
    r = cce_vs_dose_rate([20.0], device_kwargs={"T": 350.0, "N_A": 5e18})
    assert r["solver"] == "batched"
    assert (got["T"], got["N_A"]) == (350.0, 5e18)
    assert got["N_D_junction"] == DEFAULT_STACK["N_D_junction"]
