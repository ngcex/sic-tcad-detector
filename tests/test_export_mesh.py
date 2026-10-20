"""Mesh / profile export: carrier columns, VTU writer, readback by independent parsers.

The legacy-VTK / CSV basics are in tests/test_export.py; here: Electrons / Holes, the
XML UnstructuredGrid (VTU) file parsed with xml.etree and with meshio (dev extra), and
the engine2d-vs-devsim carrier parity behind the documented provenance note.
"""

from __future__ import annotations

import csv
import dataclasses
import io
import xml.etree.ElementTree as ET

import numpy as np
import pytest

import etna
from etna import DeviceConfig
from etna.api.results import MeshData, SimResult
from etna.export import VALIDATION_LABEL, profile_to_csv, profile_to_vtk
from etna.export.mesh_export import VTK_LINE, VTK_QUAD, VTK_VERTEX, profile_to_vtu


def _fake_1d(with_carriers=True):
    x_cm = np.array([5.0, 0.0, 0.3, 11.3]) * 1e-4  # unsorted on purpose
    nv = {
        "NetDoping": np.array([8e13, -1e18, 3e15, 1e19]),
        "Potential": np.array([18.0, -1.5, -1.4, 31.6]),
        "ElectricField": np.array([np.nan, 0.0, 1e5, 0.0]),
    }
    if with_carriers:
        nv["Electrons"] = np.array([1e-13, 1e2, 1e5, 1e19])
        nv["Holes"] = np.array([1e-12, 1e18, 1e11, 1e-36])
    return SimResult(
        DeviceConfig(),
        "field",
        x_cm * 1e4,
        nv["ElectricField"],
        dict(bias_V=-20.0, solver="devsim", validated_domain=True, domain_note=None),
        mesh=MeshData(x_cm, None, nv, [], []),
    )


def _vtu(text):
    root = ET.fromstring(text)  # well-formed XML or this raises
    assert root.tag == "VTKFile" and root.get("type") == "UnstructuredGrid"
    piece = root.find("UnstructuredGrid/Piece")
    arrays = {}
    for da in root.iter("DataArray"):
        name = da.get("Name") or "Points"
        dtype = float if da.get("type").startswith("Float") else int
        arrays[name] = np.array((da.text or "").split(), dtype=dtype)
    return root, piece, arrays


def test_csv_has_carrier_columns_when_present():
    text = profile_to_csv(_fake_1d())
    rows = [r for r in csv.reader(io.StringIO(text)) if r and not r[0].startswith("#")]
    assert rows[0] == [
        "x_um",
        "NetDoping_cm-3",
        "Potential_V",
        "ElectricField_V_per_cm",
        "Electrons_cm-3",
        "Holes_cm-3",
    ]
    assert float(rows[1][4]) == 1e-13 and float(rows[4][5]) == 1e-36
    base = [
        r
        for r in csv.reader(io.StringIO(profile_to_csv(_fake_1d(False))))
        if r and not r[0].startswith("#")
    ]
    assert len(base[0]) == 4  # older results: base columns only


def test_legacy_vtk_has_carriers_sorted_by_depth():
    vtk = profile_to_vtk(_fake_1d())
    assert "SCALARS Electrons_cm-3 double 1" in vtk
    i = vtk.splitlines().index("SCALARS Holes_cm-3 double 1")
    vals = np.array(vtk.splitlines()[i + 2].split(), dtype=float)
    np.testing.assert_array_equal(vals, [1e18, 1e11, 1e-12, 1e-36])  # x sorted


def test_vtu_1d_lines_and_values():
    text = profile_to_vtu(_fake_1d())
    assert VALIDATION_LABEL in text
    root, piece, a = _vtu(text)
    assert piece.get("NumberOfPoints") == "4" and piece.get("NumberOfCells") == "3"
    np.testing.assert_allclose(a["Points"].reshape(-1, 3)[:, 0], [0, 0.3, 5, 11.3])
    np.testing.assert_array_equal(a["types"], [VTK_LINE] * 3)
    np.testing.assert_array_equal(a["connectivity"], [0, 1, 1, 2, 2, 3])
    np.testing.assert_array_equal(a["offsets"], [2, 4, 6])
    np.testing.assert_allclose(a["Potential_V"], [-1.5, -1.4, 18.0, 31.6])
    assert np.isnan(a["ElectricField_V_per_cm"][2])
    np.testing.assert_array_equal(a["Electrons_cm-3"], [1e2, 1e5, 1e-13, 1e19])
    assert a["bias_V"][0] == -20.0 and a["validated_domain"][0] == 1


def test_vtu_non_grid_2d_uses_vertices():
    x = np.array([0.0, 1.0, 0.5]) * 1e-4
    y = np.array([0.0, 0.0, 1.0]) * 1e-4
    nv = {
        k: np.array([1.0, 2.0, 3.0])
        for k in ("NetDoping", "Potential", "ElectricField")
    }
    res = SimResult(
        DeviceConfig(half_width_um=1.0),
        "field",
        np.array([]),
        np.array([]),
        dict(bias_V=-1.0, solver="devsim"),
        MeshData(x, y, nv, [], []),
    )
    _, piece, a = _vtu(profile_to_vtu(res))
    assert piece.get("NumberOfCells") == "3"
    np.testing.assert_array_equal(a["types"], [VTK_VERTEX] * 3)
    assert "Electrons_cm-3" not in a


def test_meshio_reads_vtu_and_vtk(tmp_path):
    meshio = pytest.importorskip("meshio")
    res = _fake_1d()
    (tmp_path / "p.vtu").write_text(profile_to_vtu(res))
    (tmp_path / "p.vtk").write_text(profile_to_vtk(res))
    m = meshio.read(tmp_path / "p.vtu")
    np.testing.assert_allclose(m.points[:, 0], [0, 0.3, 5, 11.3])
    assert m.cells[0].type == "line" and len(m.cells[0].data) == 3
    np.testing.assert_array_equal(
        m.point_data["Holes_cm-3"], [1e18, 1e11, 1e-12, 1e-36]
    )
    v = meshio.read(tmp_path / "p.vtk")
    np.testing.assert_allclose(
        np.ravel(v.point_data["Potential_V"]), [-1.5, -1.4, 18.0, 31.6]
    )


@pytest.fixture(scope="module")
def field2d():
    return etna.run_field(DeviceConfig(half_width_um=25.0), bias_V=-30.0)


def test_run_field_2d_engine_carries_carriers(field2d, tmp_path):
    assert field2d.metadata["solver"] == "engine2d"
    m = field2d.mesh
    for k in ("Electrons", "Holes"):
        v = m.node_values[k]
        assert v.shape == m.x_coords.shape and np.all(v >= 0)
    text = profile_to_vtu(field2d)
    _, piece, a = _vtu(text)
    nxu, nyu = np.unique(m.x_coords).size, np.unique(m.y_coords).size
    assert int(piece.get("NumberOfPoints")) == m.x_coords.size
    assert int(piece.get("NumberOfCells")) == (nxu - 1) * (nyu - 1)
    assert set(a["types"].tolist()) == {VTK_QUAD}
    # every point value sits at its own node
    pts = a["Points"].reshape(-1, 3)
    lookup = {
        (round(p * 1e4, 9), round(q * 1e4, 9)): i
        for i, (p, q) in enumerate(zip(m.x_coords, m.y_coords))
    }
    idx = np.array([lookup[(round(p, 9), round(q, 9))] for p, q, _ in pts])
    np.testing.assert_array_equal(a["Electrons_cm-3"], m.node_values["Electrons"][idx])
    np.testing.assert_array_equal(a["NetDoping_cm-3"], m.node_values["NetDoping"][idx])
    meshio = pytest.importorskip("meshio")
    (tmp_path / "f.vtu").write_text(text)
    mm = meshio.read(tmp_path / "f.vtu")
    assert mm.cells[0].type == "quad" and mm.points.shape[0] == m.x_coords.size


def test_run_field_1d_devsim_carriers():
    res = etna.run_field(DeviceConfig(), bias_V=-10.0)
    nv = res.mesh.node_values
    assert res.metadata["solver"] == "devsim"
    assert nv["Electrons"].shape == res.mesh.x_coords.shape
    # p+ anode is hole-rich, n+ end electron-rich (sanity of the read-out)
    o = np.argsort(res.mesh.x_coords)
    assert nv["Holes"][o][0] > 1e17 and nv["Electrons"][o][-1] > 1e18


def test_engine2d_carriers_match_devsim():
    """Backs the provenance note in etna/export/mesh_export.py (measured 2026-10-20)."""
    import devsim

    from etna.api.device import build_device, device_kwargs
    from etna.core.devsim_reset import delete_device_and_mesh, reset_devsim_fully
    from etna.core.drift_diffusion import ramp_bias
    from etna.core.engine2d.scan import dark_state_2d

    cfg = dataclasses.replace(DeviceConfig(), half_width_um=25.0)
    st = dark_state_2d(cfg.half_width_um, 20.0, **device_kwargs(cfg))
    reset_devsim_fully()
    di = build_device(cfg)
    try:
        ramp_bias(di, V_target=20.0, contact="cathode", V_step=0.5, V_step_max=2.0)
        d, r = di["device_name"], di["region_name"]

        def g(name):
            return np.array(devsim.get_node_model_values(device=d, region=r, name=name))

        np.testing.assert_allclose(g("x"), st["x"])
        np.testing.assert_allclose(g("y"), st["y"])
        for k in ("Electrons", "Holes"):
            a, b = g(k), st[k]
            big = a > 1e10
            assert big.any()
            np.testing.assert_allclose(b[big], a[big], rtol=1e-9)
            assert np.all(a[~big] < 1e10) and np.all(b[~big] < 1e10)
            # depleted region: physically zero in both
            low = a < 1.0
            assert np.all(b[low] < 1.0)
    finally:
        try:
            delete_device_and_mesh(di["device_name"])
        except Exception:
            reset_devsim_fully()
