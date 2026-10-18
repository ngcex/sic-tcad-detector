"""The default model is the p+ / graded n- epi (graded n-/n+ transition) / n+ stack
(etna.core.device.DEFAULT_STACK; since 2026-10-18 the measured-C-V fit).

* one source of truth: DeviceConfig, the core graded-doping fallbacks, cce_vs_bias,
  device2d and the 2D engine gate all read DEFAULT_STACK; N_t / E_t live in
  SiC4H_Parameters and are read back into it;
* every facade hands the config's stack (p+, doping, n+) to the device builder -- a facade
  that dropped the n+ keys would compute the legacy model under an n+ label;
* the default and the legacy preset reach -100 V in the dark sweep (E_t 0.67 made the
  no-n+ devsim sweep stop at -89.5 V; with n+ it must not).
"""


import numpy as np
import pytest

import etna
import etna.api.simulation as sim
import etna.core.device as core_device
import etna.core.device2d as core_device2d
from etna.api.device import (
    DeviceConfig,
    device_kwargs,
    legacy_config,
    measured_cv_config,
    nplus_calibrated_config,
)
from etna.core.device import DEFAULT_STACK, LEGACY_STACK, NPLUS_CALIBRATION
from etna.core.engine2d import (
    engine2d_dark_unsupported_reason,
    engine2d_unsupported_reason,
)
from etna.core.sic_material import SiC4H_Parameters


def test_one_source_of_truth():
    c = DeviceConfig()
    # 2026-10-18 flip: the default is the measured-C-V fit (was NPLUS_CALIBRATION:
    # 2.958237e15 / 8.252732e13 / 0.987 um, abrupt n+, N_t 2.805336e12)
    assert core_device.MEASURED_CV_STACK is DEFAULT_STACK
    assert NPLUS_CALIBRATION is not DEFAULT_STACK
    assert c == measured_cv_config()
    assert c.substrate_thickness_um == 0.3 and c.nplus_thickness_um == 1.0
    assert c.N_D_nplus == DEFAULT_STACK["N_D_nplus"] == 1e19
    assert c.N_D_junction == DEFAULT_STACK["N_D_junction"] == 8.185880e15
    assert c.N_D_bulk == DEFAULT_STACK["N_D_bulk"] == 3.382506e13
    assert c.L_transition_um == 0.463863
    assert c.buffer_thickness_um == 1.136634
    assert c.buffer_end_doping == DEFAULT_STACK["buffer_end_doping"] == 4.516421e16
    assert c.N_t == SiC4H_Parameters().N_t == DEFAULT_STACK["N_t"] == 2.849024e12
    assert c.E_t == SiC4H_Parameters().E_t == DEFAULT_STACK["E_t"] == 0.67
    assert core_device2d._N_D_JUNCTION_DEFAULT == DEFAULT_STACK["N_D_junction"]
    assert core_device2d._N_D_BULK_DEFAULT == DEFAULT_STACK["N_D_bulk"]
    assert core_device2d._L_TRANSITION_DEFAULT == DEFAULT_STACK["L_transition"]
    assert (
        core_device2d._BUFFER_THICKNESS_DEFAULT_CM
        == DEFAULT_STACK["buffer_thickness_cm"]
    )
    assert (
        core_device2d._BUFFER_END_DOPING_DEFAULT == DEFAULT_STACK["buffer_end_doping"]
    )
    import inspect

    for fn in (core_device.create_sic_device, core_device2d.create_sic_2d_device):
        p = inspect.signature(fn).parameters
        assert p["buffer_thickness_cm"].default == DEFAULT_STACK["buffer_thickness_cm"]
        assert p["buffer_end_doping"].default == DEFAULT_STACK["buffer_end_doping"]
    # the previous default stays reproducible as a named preset
    p = nplus_calibrated_config()
    assert (p.N_D_junction, p.N_D_bulk, p.L_transition_um) == (
        2.958237e15,
        8.252732e13,
        0.987,
    )
    assert p.buffer_thickness_um is None and p.buffer_end_doping is None
    assert (p.N_t, p.E_t) == (2.805336e12, 0.67)
    # device_kwargs always carries the transition keys (None = abrupt epi end)
    kp = device_kwargs(p)
    assert kp["buffer_thickness_cm"] is None and kp["buffer_end_doping"] is None
    # N_t / E_t are not builder keywords (the engine gates reject unknown keys)
    assert "N_t" not in device_kwargs(c) and "E_t" not in device_kwargs(c)


def test_default_config_inside_every_engine_domain():
    c = DeviceConfig()
    assert sim.validated_domain_reason(c) is None
    kw = device_kwargs(c)
    kw2 = {k: v for k, v in kw.items() if k not in ("epi_thickness_cm", "N_D_bulk")}
    assert (
        engine2d_unsupported_reason(
            50.0, 60.0, kw["epi_thickness_cm"], kw["N_D_bulk"], 10, **kw2
        )
        is None
    )
    assert engine2d_dark_unsupported_reason(50.0, 50.0, **kw) is None
    legacy = device_kwargs(legacy_config())
    # explicit None: the builders default to the stack's n+ layer since 2026-10-17
    assert legacy["nplus_thickness_cm"] is None


def test_legacy_preset_matches_the_old_default():
    c = legacy_config()
    assert (c.substrate_thickness_um, c.nplus_thickness_um) == (1.0, None)
    assert (c.N_D_junction, c.N_D_bulk, c.L_transition_um) == (2.93e15, 8.82e13, 0.987)
    assert (
        (c.N_t, c.E_t) == (LEGACY_STACK["N_t"], LEGACY_STACK["E_t"]) == (2.8e12, 0.65)
    )
    assert legacy_config(T=310.0).T == 310.0


def test_legacy_preset_reproduces_the_pre_flip_numbers():
    """Values of DeviceConfig() before the n+ flip (commit f60be71, scratch baseline)."""
    c = legacy_config()
    cv = sim.run_cv(c, 0.0, -30.0, 31)
    W = np.asarray(cv.metadata["depletion_widths"]) * 1e4
    assert W[0] == pytest.approx(1.6970381228506364, rel=1e-9)
    assert W[10] == pytest.approx(9.51526809913126, rel=1e-9)
    assert W[30] == pytest.approx(9.984257844169887, rel=1e-9)
    assert cv.metadata["reach_through_V"] == 17.0
    cce = sim.run_cce(c, -10.0, -50.0, 5)
    # hole mobility at the epi doping since 2026-10-17 (was N_A, audit-v7 A-03); pre-flip values were
    # 0.9632660166 0.9961156665 0.9977188375 0.9982665005 0.9984736379; both mobilities at
    # N_D_bulk (was DEFAULT_N_D 1.07e15) since 2026-10-17: 0.9636306340 0.9962744523
    # 0.9978247277 0.9983467119 0.9985384182 before; local per-node mobility (p+ at its
    # own 1e19, A-03) since 2026-10-18: 0.9636418771662469 0.9962794936625294
    # 0.9978281262365173 0.9983493034607709 0.998540520801036 before
    np.testing.assert_allclose(
        cce.y,
        [
            0.9636400757290315,
            0.9962782971599875,
            0.9978271896050328,
            0.9983485188030706,
            0.9985398378691338,
        ],
        rtol=1e-9,
    )
    d = sim.run_dark_current(c, 0.0, -30.0, 4)  # default area: 4 mm2
    np.testing.assert_allclose(
        d.y * 1e12,
        # hole mobility at the epi doping (2026-10-17, A-03); pre-flip:
        # 1.5958223017 5.9621516945 10.4615603668 14.5641749059; mobilities at N_D_bulk
        # (2026-10-17): 1.5988415042 5.9622111020 10.4615631238 14.5641756901 before;
        # local per-node mobility (p+ electrons at their own 1e19, A-03; 2026-10-18):
        # 1.5989592488465492 5.962238289889233 10.461588518070496 14.564200940400156 before
        [1.5940349754412595, 5.957760064773854, 10.457105099357713, 14.5597135908548],
        rtol=1e-9,
    )
    assert d.metadata["N_t"] == 2.8e12 and d.metadata["E_t"] == 0.65


def test_default_numbers_are_the_measured_cv_record():
    """Default since 2026-10-18 (measured-C-V fit, lab/cv_gap + lab/buffer_profile):
    eps*A/C 1.655 / 9.519 / 9.750 um at 0 / -10 / -30 V (digitised labels 1.68 / 9.52 /
    9.73; the whole curve is a fit input), depletion edge 1.649 / 9.480 / 9.736 um,
    reach-through 3.1 V (edge within one L_D of the transition start), CCE 0.9994734 /
    0.9998111 / 0.9998832 at -10 / -30 / -50 V (local per-node mobility since 2026-10-18,
    A-03; scalar epi mobility before: 0.9994722 / 0.9998104 / 0.9998829), dark 6.296 /
    15.763 / 18.000 pA (4 mm2; N_t refitted to the 18 pA bound; 17.99985 pA at -60 V
    with the local mobility)."""
    c = DeviceConfig()
    cv = sim.run_cv(c, 0.0, -30.0, 31)
    WC = np.asarray(cv.metadata["W_from_C"]) * 1e4
    np.testing.assert_allclose(
        WC[[0, 10, 30]], [1.654732, 9.518526, 9.75029], atol=1e-6
    )
    W = np.asarray(cv.metadata["depletion_widths"]) * 1e4
    np.testing.assert_allclose(
        W[[0, 10, 30]], [1.649339, 9.480291, 9.736051], atol=1e-6
    )
    assert cv.metadata["reach_through_V"] == pytest.approx(3.1)
    cce = sim.run_cce(c, -10.0, -50.0, 3)
    np.testing.assert_allclose(cce.y, [0.9994734, 0.9998111, 0.9998832], atol=1e-7)
    d = sim.run_dark_current(c, 0.0, -60.0, 7)  # default area: 4 mm2
    np.testing.assert_allclose(
        d.y[[1, 3, 6]] * 1e12, [6.296, 15.763, 18.000], atol=1e-3
    )


def test_previous_default_numbers_are_the_calibration_record():
    """Previous default (nplus_calibrated_config; 2026-10-17 record, on the small-signal
    C; lab/solver2d/notes/LOG.md): eps*A/C 1.700 / 9.500 (inputs) / 9.979 um (out of fit,
    +2.56 %), depletion edge 1.689 / 9.267 / 9.896 um, dark 6.077 / 14.720 / 18.000 pA
    (4 mm2). Previous record (edge calibration): edge 1.700 / 9.500 / 9.899,
    CCE 0.999021 / 0.999706 / 0.999819, dark 6.132 / 14.873 / 18.000, reach-through 10.3 V.
    """
    c = nplus_calibrated_config()
    cv = sim.run_cv(c, 0.0, -30.0, 31)
    WC = np.asarray(cv.metadata["W_from_C"]) * 1e4
    np.testing.assert_allclose(WC[[0, 10, 30]], [1.7, 9.5, 9.978827], rtol=0, atol=1e-6)
    W = np.asarray(cv.metadata["depletion_widths"]) * 1e4
    np.testing.assert_allclose(
        W[[0, 10, 30]], [1.688703, 9.266588, 9.896195], rtol=0, atol=1e-6
    )
    # Debye-length reach-through rule on the depletion edge (a definition; paper: ~10 V)
    assert cv.metadata["reach_through_V"] == pytest.approx(12.3)
    cce = sim.run_cce(c, -10.0, -50.0, 3)
    # hole mobility at the epi doping since 2026-10-17 (was N_A, audit-v7 A-03): 0.998899 / 0.999707 / 0.999820 before;
    # both mobilities at N_D_bulk (was DEFAULT_N_D 1.07e15) since 2026-10-17: 0.999303 / 0.999812 / 0.999884 before;
    # local per-node mobility (A-03) since 2026-10-18: 0.999316 / 0.999816 / 0.999887 before
    np.testing.assert_allclose(cce.y, [0.9993150, 0.9998155, 0.9998862], atol=1e-7)
    d = sim.run_dark_current(c, 0.0, -60.0, 7)  # default area: 4 mm2
    np.testing.assert_allclose(
        d.y[[1, 3, 6]] * 1e12, [6.077, 14.720, 18.000], atol=1e-3
    )


def test_dark_sweep_reaches_minus_100V():
    d = sim.run_dark_current(DeviceConfig(), 0.0, -100.0, 20)
    assert len(d.x) == 20 and d.x[-1] == pytest.approx(-100.0)
    assert np.all(np.isfinite(d.y))


def test_legacy_dark_sweep_truncation_is_honest_and_the_plateau_smooth():
    """Legacy (cathode on the n- epi) dark sweep: no fixed '-60 V wall', but a fragile one.

    With 0.5 V ramp steps (before 2026-10-18) the TAT sweep stopped converging at
    -57.9 V; with the 2 V steps of DARK_RAMP_STEP_V it converges to -100 V on macOS/
    Accelerate BLAS (agrees with the 0.5 V points to 1.9e-14 where both exist), but on the
    ubuntu CI runner (OpenBLAS/UMFPACK) it stops at -78.9 V with 'final re-solve failed'.
    The cathode-on-epi model has no n+ layer to stop the contact, so how far Newton gets
    past -60 V depends on the linear-algebra build: that is a property of the legacy
    model (documented wall), so the test pins what must hold on EVERY platform: the sweep
    gets well past the old -57.9 V artefact, the truncation record is self-consistent (never
    silent), and the converged points are a smooth plateau. The truncation record itself is
    pinned by test_dark_current.py::test_sweep_reports_truncation."""
    d = sim.run_dark_current(legacy_config(), 0.0, -100.0, 20)
    reached = d.metadata["reached_v_stop"]
    assert reached <= -70.0  # past the 0.5 V-step artefact (-57.9 V) on every platform
    assert d.metadata["truncated"] == (reached != -100.0)  # never silent
    assert len(d.y) == len(d.x) and np.all(np.isfinite(d.y))
    plateau = d.y[d.x <= -60.0]
    assert len(plateau) >= 3
    assert np.all(np.diff(plateau) >= 0) and plateau[-1] / plateau[0] - 1 < 1e-3
    # ... but past reach-through the contact sits in depleted epi: flagged per point
    rt = d.metadata["reach_through_V"]
    assert rt == pytest.approx(run_cv_legacy_reach_through(), abs=1e-9)
    np.testing.assert_array_equal(d.metadata["past_reach_through"], np.abs(d.x) > rt)
    assert (
        d.metadata["past_reach_through"][-1] and not d.metadata["past_reach_through"][0]
    )
    assert "outside model validity" in d.metadata["validity_note"]
    # the default (n+) model carries no such record
    n = sim.run_dark_current(DeviceConfig(), 0.0, -10.0, 3)
    assert (
        n.metadata["validity_note"] is None and n.metadata["past_reach_through"] is None
    )


def test_legacy_dark_sweep_before_reach_through_is_all_valid():
    """A legacy sweep ending before reach-through (0..-10 V < 17 V) still knows where
    reach-through lies (fixed lattice, not the caller's grid): every point valid."""
    d = sim.run_dark_current(legacy_config(), 0.0, -10.0, 3)
    assert d.metadata["reach_through_V"] == pytest.approx(17.0)
    np.testing.assert_array_equal(d.metadata["past_reach_through"], [False] * 3)


def run_cv_legacy_reach_through():
    """The legacy reach-through as run_cv reports it (same engine rule, 0.05 V)."""
    return sim.run_cv(legacy_config(), 0.0, -40.0, 21).metadata["reach_through_V"]


class _Spy:
    def __init__(self, monkeypatch):
        self.calls = []
        real1 = core_device.create_sic_device
        real2 = core_device2d.create_sic_2d_device

        def sic(*a, **k):
            self.calls.append(("1d", k))
            return real1(*a, **k)

        def sic2d(*a, **k):
            self.calls.append(("2d", k))
            return real2(*a, **k)

        monkeypatch.setattr(core_device, "create_sic_device", sic)
        monkeypatch.setattr("etna.core.drift_diffusion.create_sic_device", sic)
        monkeypatch.setattr(core_device2d, "create_sic_2d_device", sic2d)
        monkeypatch.setattr(
            "etna.core.charge_collection_2d.create_sic_2d_device", sic2d
        )

        self.batches = []
        real_batch = sim._engine_batch

        def batch(kws):
            self.batches.extend(kws)
            return real_batch(kws)

        monkeypatch.setattr(sim, "_engine_batch", batch)

    def assert_stack(self, cfg):
        assert self.calls or self.batches
        for _, k in self.calls:
            assert k.get("nplus_thickness_cm") == pytest.approx(
                cfg.nplus_thickness_um * 1e-4
            )
            assert k.get("N_D_nplus") == cfg.N_D_nplus
            assert k.get("substrate_thickness_cm") == pytest.approx(
                cfg.substrate_thickness_um * 1e-4
            )
            assert k.get("N_D_junction") == cfg.N_D_junction
            assert k.get("N_D_bulk") == cfg.N_D_bulk
            assert k.get("L_transition") == pytest.approx(cfg.L_transition_um * 1e-4)
        for k in self.batches:
            assert k["nplus_thickness_cm"] == pytest.approx(1e-4)
            assert k["substrate_thickness_cm"] == pytest.approx(0.3e-4)


FACADES = {
    "run_cv": lambda c: sim.run_cv(c, 0.0, -20.0, 3),
    "run_cce": lambda c: sim.run_cce(c, -10.0, -20.0, 2),
    "run_field_1d": lambda c: sim.run_field(c, -20.0),
    "run_dark_current": lambda c: sim.run_dark_current(c, 0.0, -20.0, 3),
    "run_radiation_damage": lambda c: sim.run_radiation_damage(c, np.array([1e13])),
    "run_temperature_sweep": lambda c: sim.run_temperature_sweep(c, np.array([300.0])),
    "run_flash_recombination": lambda c: sim.run_flash_recombination(
        c, np.array([20.0])
    ),
    "run_transient": lambda c: sim.run_transient(c, np.array([20.0])),
}


@pytest.mark.parametrize("name", list(FACADES))
def test_facade_hands_the_stack_to_the_builder(monkeypatch, name):
    spy = _Spy(monkeypatch)
    cfg = DeviceConfig()
    res = FACADES[name](cfg)
    assert len(res.y) > 0 or name == "run_field_1d"
    spy.assert_stack(cfg)


def test_run_field_2d_hands_the_stack_to_the_builder(monkeypatch):
    spy = _Spy(monkeypatch)
    cfg = DeviceConfig(half_width_um=25.0)
    res = sim.run_field(cfg, -20.0)
    assert res.metadata["solver"] == "engine2d"
    assert spy.calls and all(kind == "2d" for kind, _ in spy.calls)
    spy.assert_stack(cfg)


def test_microdosimetry_builds_no_device(monkeypatch, tmp_path):
    spy = _Spy(monkeypatch)
    csv = tmp_path / "ev.csv"
    csv.write_text("event_id,x,y,z,edep\n0,0,0,0,10\n1,0,0,0,20\n2,0,0,0,15\n")
    sim.run_microdosimetry(DeviceConfig(), str(csv))
    assert spy.calls == [] and spy.batches == []


def test_transient_uses_the_anode_with_nplus():
    r = sim.run_transient(DeviceConfig(), np.array([20.0, 100.0]))
    assert r.metadata["current_contact"] == "anode"
    # 2026-10-18 (local mobility): 0.999725 / 0.999869 at 20 / 100 Gy/s; balanced
    # steady-state CCE at 20 Gy/s 0.999809; direct anode round-off ~1e-4
    np.testing.assert_allclose(r.y, [0.9997250, 0.9998694], atol=1e-6)
    r = sim.run_transient(legacy_config(), np.array([20.0]))
    assert r.metadata["current_contact"] == "cathode"
    # pre-flip value 0.9975336675676356; hole mobility at the epi doping since 2026-10-17 (was N_A, audit-v7 A-03)
    # 0.9976403048899449; both mobilities at N_D_bulk (was DEFAULT_N_D) since 2026-10-17
    # 0.9976437322922994; step-end stamps/BDF1 estimator (A-04) 0.9976437319364355; local
    # per-node mobility (A-03) since 2026-10-18
    assert r.y[0] == pytest.approx(0.9976427859179097, rel=1e-9)


def test_run_all_default_every_section_ok():
    rep = etna.run_all(DeviceConfig())
    got = {s.name: (s.status, s.solver) for s in rep.sections}
    assert all(st == "ok" for st, _ in got.values()), got
    assert got["cv"][1] == "batched" and got["cce"][1] == "batched"
    rep2 = etna.run_all(DeviceConfig(half_width_um=50.0), sections=("field", "sweep"))
    got2 = {s.name: (s.status, s.solver) for s in rep2.sections}
    assert got2 == {"field": ("ok", "engine2d"), "sweep": ("ok", "engine2d")}, got2
    table = rep2.sections[1].result.metadata["table"]
    assert {"nplus_thickness_um", "N_D_nplus", "nplus_generation"} <= set(table.columns)
    assert (table["nplus_thickness_um"] == 1.0).all() and not table[
        "nplus_generation"
    ].any()
    # legacy 2D sweep row (no n+): NaN / False stack columns. One 20 V engine row instead of
    # the run_all legacy sweep, whose 60 V rows fall back to devsim near the wall (~35 s).
    from etna.core.optimization import microdosimetric_sweep

    t3 = microdosimetric_sweep(
        half_widths_um=(25,),
        epi_thicknesses_cm=(10e-4,),
        N_D_bulks=(8.82e13,),
        V_biases=(20,),
        n_lateral_points=2,
        nplus_thickness_cm=None,
        substrate_thickness_cm=1e-4,
        N_D_junction=2.93e15,
        L_transition=0.987e-4,
    )
    assert t3["solver"].tolist() == ["engine2d"]
    assert t3["nplus_thickness_um"].isna().all() and not t3["nplus_generation"].any()


def test_uniform_profile_does_not_inherit_the_fitted_transition():
    """The transition was fitted with the graded epi: a uniform device stays abrupt."""
    import etna
    from app.components.device_sidebar import assemble_config
    from etna.api.device import device_kwargs

    c = etna.DeviceConfig(doping_profile="uniform", N_D=1e14)
    assert c.buffer_thickness_um is None and c.buffer_end_doping is None
    assert device_kwargs(c)["buffer_thickness_cm"] is None
    # an explicit non-default transition is kept
    c2 = etna.DeviceConfig(
        doping_profile="uniform",
        N_D=1e14,
        buffer_thickness_um=0.5,
        buffer_end_doping=1e16,
    )
    assert c2.buffer_thickness_um == 0.5
    # the graded default keeps it
    assert etna.DeviceConfig().buffer_thickness_um == pytest.approx(
        1.136634, rel=1e-12
    )
    # the sidebar path (uniform chosen there) goes through the same rule
    v = dict(
        epi_thickness_um=10.0,
        substrate_thickness_um=0.3,
        half_width_um=None,
        N_A=1e19,
        doping_profile="uniform",
        N_D=1e14,
        N_D_junction=2.9e15,
        N_D_bulk=8e13,
        L_transition_um=0.9,
        T=300.0,
        area_cm2=0.04,
    )
    assert assemble_config(v).buffer_thickness_um is None


def test_core_builders_uniform_profile_is_abrupt():
    import uuid

    import devsim

    from etna.core.device import create_sic_device

    info = create_sic_device(
        device_name=f"u_{uuid.uuid4().hex[:6]}", doping_profile="uniform", N_D=1e14
    )
    try:
        assert info.get("buffer_thickness_cm") is None
    finally:
        devsim.delete_device(device=info["device_name"])


def test_bare_builder_is_the_documented_uniform_test_diode():
    """Decision 2026-10-18 (WP-A): a bare create_sic_device() stays a uniform 1.07e15
    epi on the default geometry (p+ 0.3 um, n+ 1 um, abrupt end), NOT DeviceConfig();
    the facades pass the full stack (device_kwargs carries every stack key)."""
    import uuid

    import devsim

    from etna.api.device import device_kwargs
    from etna.core.device import DEFAULT_N_D, create_sic_device

    info = create_sic_device(device_name=f"bare_{uuid.uuid4().hex[:6]}")
    try:
        assert info["doping_profile"] == "uniform" and info["N_D"] == DEFAULT_N_D
        assert info["substrate_thickness_cm"] == DEFAULT_STACK["substrate_thickness_cm"]
        assert info["nplus_thickness_cm"] == DEFAULT_STACK["nplus_thickness_cm"]
        assert info["buffer_thickness_cm"] is None
    finally:
        devsim.delete_device(device=info["device_name"])
    kw = device_kwargs(DeviceConfig())
    assert kw["doping_profile"] == "graded"
    for key in (
        "substrate_thickness_cm",
        "nplus_thickness_cm",
        "N_D_junction",
        "N_D_bulk",
        "L_transition",
        "buffer_thickness_cm",
        "buffer_end_doping",
    ):
        assert key in kw, key
