"""Allpix field maps keep their true depth scale (issue #19).

Allpix (master ElectricFieldReader / DopingProfileReader / WeightingPotentialReader,
opened 2026-10-21) maps a mesh file's z grid onto (sensor_max_z - depth, sensor_max_z),
depth = field_depth / doping_depth / potential_depth, DEFAULT THE SENSOR THICKNESS. The
emulation below applies that rule to the exported file and the returned config snippet:
before the fix the snippet set no depth key, so in a 10 um sensor the 11.3 um map was
compressed (potential integral ~47 V instead of 53.2 V at -50 V) and in a 350 um sensor
stretched 31x.
"""

import re

import numpy as np
import pytest

from etna import DeviceConfig, run_field
from etna.core import mc_exchange as mx


def _depth_key(snippet, section, key):
    body = snippet.split(f"[{section}]", 1)[1].split("[", 1)[0]
    m = re.search(rf"^{key}\s*=\s*([0-9.eE+-]+)um\s*$", body, re.M)
    return float(m.group(1)) if m else None


def _allpix_potential_drop(path, snippet, sensor_um):
    """Line integral of E_z as Allpix would place the map in a sensor_um sensor."""
    depth = _depth_key(snippet, "ElectricFieldReader", "field_depth")
    depth = sensor_um if depth is None else depth
    if depth - sensor_um > 1e-9:
        return None  # Allpix rejects the configuration
    Ez = mx.read_init(path)["values"][0, 0, :, 2]  # V/cm
    return float(np.sum(Ez) * depth / Ez.size * 1e-4)


@pytest.fixture(scope="module")
def field_50V():
    return run_field(DeviceConfig(), bias_V=-50.0)


def test_potential_drop_survives_allpix_mapping(field_50V, tmp_path):
    out = mx.export_allpix_fields(field_50V, tmp_path, nz=400)
    psi = np.asarray(field_50V.mesh.node_values["Potential"])
    o = np.argsort(np.asarray(field_50V.mesh.x_coords))
    drop = float(psi[o][-1] - psi[o][0])
    assert drop == pytest.approx(53.2, abs=0.1)  # 50 V + built-in
    th = out["thickness_um"]
    snip = out["config_snippet"]
    for sec, key in (
        ("ElectricFieldReader", "field_depth"),
        ("DopingProfileReader", "doping_depth"),
        ("WeightingPotentialReader", "potential_depth"),
    ):
        assert _depth_key(snip, sec, key) == pytest.approx(th, rel=1e-5)
    for sensor in (th, 12.0, 350.0):
        v = _allpix_potential_drop(out["paths"]["efield"], snip, sensor)
        assert v == pytest.approx(drop, rel=2e-3), sensor


def test_sensor_thickness_mismatch_warns(field_50V, tmp_path):
    with pytest.warns(UserWarning, match="Allpix rejects"):
        out = mx.export_allpix_fields(
            field_50V, tmp_path, nz=50, sensor_thickness_um=10.0
        )
    # the 10 um sensor cannot hold the 11.3 um map: Allpix refuses, nothing is stretched
    assert (
        _allpix_potential_drop(out["paths"]["efield"], out["config_snippet"], 10.0)
        is None
    )
    with pytest.warns(UserWarning, match="no field below"):
        mx.export_allpix_fields(field_50V, tmp_path, nz=50, sensor_thickness_um=350.0)
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        out = mx.export_allpix_fields(
            field_50V,
            tmp_path,
            nz=50,
            depth_range_um=(0.3, 10.3),
            sensor_thickness_um=10.0,
        )
    assert out["thickness_um"] == pytest.approx(10.0)
    assert (
        _depth_key(out["config_snippet"], "ElectricFieldReader", "field_depth") == 10.0
    )
