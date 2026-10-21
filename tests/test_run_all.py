"""etna.run_all: the one-shot orchestrator (sections, failure isolation, repeatability)."""

import json

import numpy as np
import pytest

import etna
from etna.api import run_all as ra


def test_default_config_runs_every_section():
    rep = etna.run_all(etna.DeviceConfig())
    assert [s.name for s in rep.sections] == list(ra.DEFAULT_SECTIONS)
    assert rep.ok, rep.rows()
    assert all(s.solver or s.name == "flash" for s in rep.sections)
    json.loads(rep.to_json())  # serialisable


def test_twice_in_one_process_gives_identical_numbers():
    cfg = etna.DeviceConfig()
    a, b = etna.run_all(cfg), etna.run_all(cfg)
    for sa, sb in zip(a.sections, b.sections):
        np.testing.assert_array_equal(sa.result.y, sb.result.y, err_msg=sa.name)


def test_section_failure_is_isolated(monkeypatch):
    from etna.api import simulation as s

    def boom(*a, **k):
        raise RuntimeError("synthetic")

    monkeypatch.setattr(s, "run_cce", boom)
    rep = etna.run_all(etna.DeviceConfig(), ("cv", "cce", "field"))
    status = {x.name: x.status for x in rep.sections}
    assert status == {"cv": "ok", "cce": "error", "field": "ok"}
    assert "RuntimeError: synthetic" in rep.get("cce").reason
    assert not rep.ok


def test_2d_config_skips_1d_only_sections():
    rep = etna.run_all(etna.DeviceConfig(half_width_um=25.0), ("cv", "field"))
    assert rep.get("cv").status == "skipped" and "1D" in rep.get("cv").reason
    assert rep.get("field").status == "ok" and rep.get("field").solver == "engine2d"


def test_unknown_section_rejected():
    with pytest.raises(ValueError, match="unknown section"):
        etna.run_all(etna.DeviceConfig(), ("nope",))


def test_leaked_device_is_reported(monkeypatch):
    deleted = []
    calls = []

    def alive():  # first call = snapshot before the run (none), then one leak
        calls.append(1)
        return [] if len(calls) == 1 or deleted else ["leaky"]

    monkeypatch.setattr(ra, "_devices_alive", alive)
    monkeypatch.setattr(ra, "_delete_device", lambda name: deleted.append(name))
    rep = etna.run_all(etna.DeviceConfig(), ("cv",))
    assert rep.get("cv").status == "error" and "leaked" in rep.get("cv").reason
    assert deleted == ["leaky"]


def test_sweep_is_opt_in_and_reports_rows():
    assert "sweep" in ra.SECTION_NAMES and "sweep" not in ra.DEFAULT_SECTIONS
    rep = etna.run_all(etna.DeviceConfig(), ("sweep",))
    sec = rep.get("sweep")
    assert sec.status == "ok" and sec.solver == "engine2d", sec.reason
    assert "valid" in sec.summary and "n+ 1 um" in sec.summary  # default stack
    assert len(sec.result.metadata["table"]) == 6
    json.loads(rep.to_json())


def test_sweep_outside_engine_box_is_skipped_not_run():
    rep = etna.run_all(etna.DeviceConfig(epi_thickness_um=40.0), ("sweep",))
    assert rep.get("sweep").status == "skipped"
    assert "epi_thickness_cm" in rep.get("sweep").reason


def test_sweep_runs_on_the_engine_when_config_has_nplus():
    rep = etna.run_all(etna.DeviceConfig(nplus_thickness_um=1.0), ("sweep",))
    sec = rep.get("sweep")
    assert sec.status == "ok", sec
    df = sec.result.metadata["table"]
    assert (
        set(df["solver"]) == {"engine2d"}
        and "n+ 1 um" in sec.result.metadata["summary"]
    )


def test_sweep_skipped_when_nplus_outside_engine_box():
    rep = etna.run_all(etna.DeviceConfig(nplus_thickness_um=10.0), ("sweep",))
    assert rep.get("sweep").status == "skipped" and "nplus" in rep.get("sweep").reason


def test_run_field_2d_with_nplus_uses_engine():
    res = etna.run_field(
        etna.DeviceConfig(half_width_um=25.0, nplus_thickness_um=1.0), bias_V=-20.0
    )
    assert res.metadata["solver"] == "engine2d", res.metadata


def test_2d_field_summary_reports_the_map():
    rep = etna.run_all(etna.DeviceConfig(half_width_um=25.0), ("field",))
    assert (
        "2D map" in rep.get("field").summary and "max |E|" in rep.get("field").summary
    )


def test_1d_field_summary_excludes_the_nplus_interface_spike():
    rep = etna.run_all(etna.DeviceConfig(), ("field",))
    res, summ = rep.get("field").result, rep.get("field").summary
    assert "excl. spikes at the abrupt doping steps" in summ
    assert float(summ.split("peak |E| ")[1].split(" ")[0]) < 0.7 * float(res.y.max())


def test_truncated_sweep_is_flagged_in_the_summary(monkeypatch):
    from etna.api import simulation as s

    def cut(cfg, *a, **k):
        return etna.SimResult(
            config=cfg,
            sim_type="dark_current",
            x=np.array([0.0, -10.0, -20.0]),
            y=np.array([1e-12, 2e-12, 3e-12]),
            metadata={"truncated": True, "requested_v_stop": -100.0},
        )

    monkeypatch.setattr(s, "run_dark_current", cut)
    rep = etna.run_all(etna.DeviceConfig(), sections=("dark",))
    summary = rep.get("dark").summary
    assert "TRUNCATED: reached -20 V of requested -100 V" in summary, summary
    full = etna.SimResult(
        config=etna.DeviceConfig(),
        sim_type="dark_current",
        x=np.array([0.0, -100.0]),
        y=np.array([1e-12, 2e-12]),
        metadata={"truncated": False, "requested_v_stop": -100.0},
    )
    assert "TRUNCATED" not in ra._summarize(full)


@pytest.mark.parametrize(
    "factory",
    [etna.DeviceConfig, etna.measured_cv_tabulated_config, etna.legacy_config],
    ids=["default", "tabulated", "legacy"],
)
def test_damage_default_grid_adapts_to_the_compensation_fluence(factory):
    """LEGACY damage model (radiation_traps=False): the default fluence grid starts
    below the device's compensation fluence (the old fixed 1e13..1e16 grid was entirely
    beyond it on the default and tabulated stacks)."""
    from etna.api import simulation as s

    cfg = factory(radiation_traps=False)
    phi_c = s.compensation_fluence(cfg, 5.6)
    grid = s.default_fluence_grid(cfg, 5.6)
    assert grid.size == s.DEFAULT_FLUENCE_GRID_POINTS == 8
    np.testing.assert_allclose(grid[[0, -1]], [1e-2 * phi_c, 3.0 * phi_c], rtol=1e-12)
    rep = etna.run_all(cfg, ("damage",))
    sec = rep.get("damage")
    res = sec.result
    np.testing.assert_allclose(res.x, grid, rtol=0, atol=0)
    valid = np.asarray(res.metadata["valid"])
    # every point below the compensation fluence is solved, the last one (3x) is not
    assert valid[res.x < phi_c].all() and (res.x < phi_c).sum() == 6
    assert not valid[-1] and np.isnan(res.y[-1])
    assert "full donor compensation" in res.metadata["reason"][-1]
    assert sec.status == "ok" and sec.solver == "batched", (sec.status, sec.solver)
    assert sec.summary.startswith(f"{int(valid.sum())} of 8 points valid; "), sec.summary
    assert f"full donor compensation at {phi_c:.3g} p/cm^2" in sec.summary, sec.summary


def test_damage_explicit_grid_is_used_unchanged():
    from etna.api import simulation as s

    # default 5.6 MeV: compensation 1.01e12 p/cm2 with the SR-NIEL kappa (2026-10-20)
    grid = np.array([1e10, 3e11, 5e14])
    # re-pinned 2026-10-21: the default (calibrated traps) solves beyond compensation
    res = s.run_radiation_damage(etna.DeviceConfig(), fluences=grid)
    np.testing.assert_array_equal(res.x, grid)
    assert list(res.metadata["valid"]) == [True, True, True]
    res = s.run_radiation_damage(etna.DeviceConfig(radiation_traps=False), fluences=grid)
    assert list(res.metadata["valid"]) == [True, True, False]


@pytest.mark.parametrize(
    "factory",
    [etna.DeviceConfig, etna.measured_cv_tabulated_config, etna.legacy_config],
    ids=["default", "tabulated", "legacy"],
)
def test_damage_default_grid_is_the_calibrated_trap_grid(factory):
    """Default model ("auto" -> calibrated kinetic traps, sic-damage): 10 points from
    0.01x its compensation fluence to 1e16 neq, every one solved."""
    from etna.api import simulation as s

    cfg = factory()
    phi_c = s.compensation_fluence(cfg, 5.6)
    grid = s.default_fluence_grid(cfg, 5.6)
    assert grid.size == s.TRAP_FLUENCE_GRID_POINTS == 10
    assert grid[0] == pytest.approx(1e-2 * phi_c, rel=1e-12)
    assert grid[-1] * s._kappa(5.6) == pytest.approx(1e16, rel=1e-9)
    rep = etna.run_all(cfg, ("damage",))
    sec = rep.get("damage")
    res = sec.result
    assert res.metadata["damage_model"] == "charged_traps"
    assert np.asarray(res.metadata["valid"]).all(), res.metadata["reason"]
    assert sec.status == "ok", (sec.status, sec.reason)
    assert f"charged traps: bulk compensated at ~{phi_c:.3g} p/cm^2" in sec.summary


def test_all_nan_section_is_a_labelled_warning(monkeypatch):
    from etna.api import simulation as s

    def beyond(cfg, *a, **k):
        return etna.SimResult(
            config=cfg,
            sim_type="damage",
            x=np.array([1e15, 1e16]),
            y=np.array([np.nan, np.nan]),
            metadata={
                "valid": np.array([False, False]),
                "reason": ["full donor compensation", "full donor compensation"],
                "compensation_fluence": 1.35e13,
                "solver": "none (no point solved)",
            },
        )

    monkeypatch.setattr(s, "run_radiation_damage", beyond)
    rep = etna.run_all(etna.DeviceConfig(), ("damage",))
    sec = rep.get("damage")
    assert sec.status == "warning" and not rep.ok
    assert sec.reason.startswith("no valid point (0 of 2)") and "compensation" in sec.reason
    assert "all NaN" in sec.summary and "1.35e+13" in sec.summary


def test_ok_ignores_skipped_sections():
    rep = etna.run_all(etna.DeviceConfig(half_width_um=25.0), ("cv", "field"))
    assert rep.get("cv").status == "skipped" and rep.ok
