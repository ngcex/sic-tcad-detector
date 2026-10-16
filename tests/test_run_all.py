"""etna.run_all: the one-shot orchestrator (sections, failure isolation, repeatability)."""

import json

import numpy as np
import pytest

import etna
from etna.api import run_all as ra


def test_default_config_runs_every_section():
    rep = etna.run_all(etna.DeviceConfig())
    assert [s.name for s in rep.sections] == list(ra.SECTION_NAMES)
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
    monkeypatch.setattr(ra, "_devices_alive", lambda: ["leaky"] if not deleted else [])
    monkeypatch.setattr(devsim, "delete_device", lambda device: deleted.append(device))
    rep = etna.run_all(etna.DeviceConfig(), ("cv",))
    assert rep.get("cv").status == "error" and "leaked" in rep.get("cv").reason
    assert deleted == ["leaky"]


def test_config_key_distinguishes_configs():
    a, b = etna.DeviceConfig(), etna.DeviceConfig(epi_thickness_um=12.0)
    assert ra.config_key(a) == ra.config_key(etna.DeviceConfig())
    assert ra.config_key(a) != ra.config_key(b)
