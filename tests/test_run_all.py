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
    import devsim

    deleted = []
    calls = []

    def alive():  # first call = snapshot before the run (none), then one leak
        calls.append(1)
        return [] if len(calls) == 1 or deleted else ["leaky"]

    monkeypatch.setattr(ra, "_devices_alive", alive)
    monkeypatch.setattr(devsim, "delete_device", lambda device: deleted.append(device))
    rep = etna.run_all(etna.DeviceConfig(), ("cv",))
    assert rep.get("cv").status == "error" and "leaked" in rep.get("cv").reason
    assert deleted == ["leaky"]


def test_config_key_distinguishes_configs():
    a, b = etna.DeviceConfig(), etna.DeviceConfig(epi_thickness_um=12.0)
    assert ra.config_key(a) == ra.config_key(etna.DeviceConfig())
    assert ra.config_key(a) != ra.config_key(b)


def test_sweep_is_opt_in_and_reports_rows():
    assert "sweep" in ra.SECTION_NAMES and "sweep" not in ra.DEFAULT_SECTIONS
    rep = etna.run_all(etna.DeviceConfig(), ("sweep",))
    sec = rep.get("sweep")
    assert sec.status == "ok" and sec.solver == "engine2d", sec.reason
    assert "valid" in sec.summary and "no n+" in sec.summary
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
    assert set(df["solver"]) == {"engine2d"} and "n+ 1 um" in sec.result.metadata["summary"]


def test_sweep_skipped_when_nplus_outside_engine_box():
    rep = etna.run_all(etna.DeviceConfig(nplus_thickness_um=10.0), ("sweep",))
    assert rep.get("sweep").status == "skipped" and "nplus" in rep.get("sweep").reason


def test_run_field_2d_with_nplus_uses_engine():
    res = etna.run_field(etna.DeviceConfig(half_width_um=25.0, nplus_thickness_um=1.0), bias_V=-20.0)
    assert res.metadata["solver"] == "engine2d", res.metadata
