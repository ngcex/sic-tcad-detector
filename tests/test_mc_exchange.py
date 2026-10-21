"""Allpix Squared exchange formats (etna.core.mc_exchange), issue #7."""

import numpy as np
import pytest

from etna import DeviceConfig, run_field
from etna.api.results import peak_field_1d
from etna.core import mc_exchange as mx
from etna.core import mc_transport as mt
from etna.core.mc_coupling import process_mc_ensemble


@pytest.fixture(scope="module")
def tracks():
    return mt.sample_track_events("proton", 62.0, n_events=500, max_hits=32, seed=4)


def test_deposition_csv_round_trip_into_cce_chain(tracks, tmp_path):
    p = tmp_path / "dep.csv"
    snip = mx.write_deposition_csv(p, tracks, unit_length="mm", unit_energy="MeV")
    assert 'model = "csv"' in snip and 'unit_length = "mm"' in snip
    text = p.read_text()
    assert text.endswith("\n")
    assert text.count("Event: ") == tracks.n_events
    assert "\nEvent: 0\n" in text
    df = mx.read_deposition_csv(p)  # Allpix defaults mm / MeV / ns
    per = df.groupby("event_id")["edep_keV"].sum()
    assert np.allclose(per.values, tracks.edep_keV[per.index.values], rtol=1e-7)
    # parent ids exist before use, primary parent 0
    for _, g in df.groupby("event_id"):
        seen = {0}
        for trk, par in zip(g["track_id"], g["parent_id"]):
            assert par in seen
            seen.add(trk)
    # z: beam from the implant side -> z' in [-t/2, t/2], primary hits move to -z
    assert df["z_cm"].abs().max() <= 5e-4 + 1e-12
    r = process_mc_ensemble(df, lambda L: np.ones_like(L), sv_thickness_um=10.0)
    pos = tracks.edep_keV > 0
    assert np.allclose(
        np.sort(r["event_energies_keV"]), np.sort(tracks.edep_keV[pos]), rtol=1e-7
    )


def test_read_deposition_variants(tmp_path):
    p = tmp_path / "x.csv"
    p.write_text(
        "# comment\n\nEvent: 0\n211, 3.2e+01, 1.0e-02, -2.5, 4.4, -0.25, MyDetector, 1, 0\n"
        "Event: 1\n11, 1.0e-03, 1.0, 2.0, 3.0, other, 2, 1\n"
    )
    with pytest.raises(Exception):
        mx.read_deposition_csv(p)
    p.write_text(
        "Event: 0\n211, 1.0e-02, -2.5, 4.4, -0.25, MyDetector\nEvent: 1\n11, 1.0, 1.0, 2.0, 3.0, other\n"
    )
    df = mx.read_deposition_csv(
        p, assign_timestamps=False, create_mcparticles=False, detector="MyDetector"
    )
    assert len(df) == 1 and df["edep_keV"].iloc[0] == pytest.approx(10.0)
    assert df["x_cm"].iloc[0] == pytest.approx(-0.25)


def test_tracks_to_dataframe_carries_kinetic_energy(tracks):
    df = mx.tracks_to_dataframe(tracks)
    r = process_mc_ensemble(df, lambda L: np.ones_like(L), sv_thickness_um=10.0)
    assert np.allclose(r["event_kinetic_MeV"], 62.0)  # bare SV: entry energy = beam
    sur = mt.sample_track_events("proton", 62.0, n_events=200, surround_um=100.0, seed=1)
    e_in = sur.E_sv_entry_MeV[sur.chord_um > 0]
    assert np.all(e_in < 62.0) and np.all(e_in > 61.0)


def test_init_round_trip_asymmetric_grid(tmp_path):
    nx, ny, nz = 3, 4, 5
    i, j, k = np.indices((nx, ny, nz))
    v = np.stack([100 * i + 10 * j + k, -i, 1.5 * k], axis=-1).astype(float)
    p = tmp_path / "f.init"
    mx.write_init(p, v, (30.0, 40.0, 10.0), "V/cm", header="test")
    lines = p.read_text().splitlines()
    assert lines[1].split()[0] == "V/cm"
    hdr = lines[4].split()
    assert hdr[:3] == ["10", "30", "40"]  # thickness, size_x, size_y (um)
    assert hdr[7:10] == ["3", "4", "5"]  # bins x, y, z
    assert lines[5].split()[:3] == ["1", "1", "1"]  # 1-based, z fastest
    assert lines[6].split()[:3] == ["1", "1", "2"]
    r = mx.read_init(p)
    assert r["shape"] == (nx, ny, nz) and r["size_um"] == (30.0, 40.0, 10.0)
    assert np.array_equal(r["values"], v)
    s = tmp_path / "s.init"
    mx.write_init(s, v[..., 0], (30, 40, 10), "")
    assert s.read_text().splitlines()[1].startswith("internal")
    assert np.array_equal(mx.read_init(s)["values"][..., 0], v[..., 0])
    with pytest.raises(ValueError):
        mx.write_init(p, np.zeros((2, 2)), (1, 1, 1), "V/cm")


def test_export_allpix_fields_1d(tmp_path):
    res = run_field(DeviceConfig(), bias_V=-30.0)
    out = mx.export_allpix_fields(res, tmp_path, nz=400)
    ef = mx.read_init(out["paths"]["efield"])
    th = out["thickness_um"]
    Ez = ef["values"][0, 0, :, 2]
    psi = np.asarray(res.mesh.node_values["Potential"])
    depth = np.asarray(res.mesh.x_coords)
    o = np.argsort(depth)
    # line integral of E_z over the sensor = psi(back) - psi(implant side)
    assert np.sum(Ez) * th / 400 * 1e-4 == pytest.approx(
        psi[o][-1] - psi[o][0], rel=0.01
    )
    # field points to +z (implant / p+ side) where depleted; |E| peak close to run_field
    assert Ez.max() > 0 and abs(Ez.min()) < 1e-3 * Ez.max()
    # away from the abrupt doping steps the exported field equals run_field's |E|
    zdep = th - (np.arange(400) + 0.5) / 400 * th  # ETNA depth of the Allpix cells
    for d_um in (1.0, 2.0, 3.0):
        k = int(np.argmin(np.abs(zdep - d_um)))
        assert Ez[k] == pytest.approx(np.interp(zdep[k], res.x, res.y), rel=0.01)
    peak, _, _ = peak_field_1d(res.x, res.y, res.config)  # skips nm spikes at steps
    assert Ez.max() == pytest.approx(peak, rel=0.15)
    dop = mx.read_init(out["paths"]["doping"])["values"][0, 0, :, 0]
    assert dop[-1] < 0 < dop[0]  # p+ (negative) at the implant side, n at the back
    wp = mx.read_init(out["paths"]["weighting"])["values"][0, 0, :, 0]
    assert wp[0] < wp[-1] and wp[-1] < 1.0
    assert "ElectricFieldReader" in out["config_snippet"]


def test_export_allpix_fields_2d_matches_1d(tmp_path):
    # the default 2D device is laterally uniform: centre column = 1D field
    r2 = run_field(DeviceConfig(half_width_um=50.0), bias_V=-30.0)
    out = mx.export_allpix_fields(r2, tmp_path, nz=100, nx=21)
    assert out["size_um"][0] == pytest.approx(100.0)  # half-cell mirrored to full plane
    assert "PIXEL_FULL" in out["config_snippet"]
    v = mx.read_init(out["paths"]["efield"])["values"]
    assert v.shape == (21, 1, 100, 3)
    th = out["thickness_um"]
    zdep = th - (np.arange(100) + 0.5) / 100 * th
    r1 = run_field(DeviceConfig(), bias_V=-30.0)
    k = int(np.argmin(np.abs(zdep - 2.0)))
    assert v[10, 0, k, 2] == pytest.approx(np.interp(zdep[k], r1.x, r1.y), rel=0.02)
    assert np.max(np.abs(v[..., 0])) < 1e-3 * np.max(v[..., 2])
