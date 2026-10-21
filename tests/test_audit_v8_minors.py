"""Regression tests for the minor findings of .planning/AUDIT_v8_MC_IO_APP.md (5-17).

Each test failed on main at edf4483 (the finding) and passes after the fix.
"""

from __future__ import annotations

import warnings
import xml.etree.ElementTree as ET

import numpy as np
import pytest

import etna
import etna.measurements as M
from etna import DeviceConfig, SimResult
from etna.core import mc_exchange as mx
from etna.core import mc_transport as mt
from etna.core.mc_coupling import process_mc_ensemble
from etna.core.microdosimetry import (
    compute_kappa_table,
    lineal_energy_spectrum,
    tissue_equivalence_correction,
)
from tests._measurement_synth import cv_csv, it_csv


# ---- 5: zero-deposit events do not scale y_F ----------------------------------------
def test_zero_deposit_events_do_not_bias_y_F():
    e = np.random.default_rng(1).gamma(4.0, 5.0, 20000)
    a = lineal_energy_spectrum(e, 10.0)
    b = lineal_energy_spectrum(np.concatenate([e, np.zeros(e.size)]), 10.0)
    assert b["y_F"] == pytest.approx(a["y_F"], rel=1e-12)  # was 0.5x
    assert b["y_D"] == pytest.approx(a["y_D"], rel=1e-12)
    assert np.sum(b["f_y"] * b["bin_widths"]) == pytest.approx(1.0, rel=1e-12)
    assert (b["n_events"], b["n_zero"], b["n_in_range"]) == (40000, 20000, 20000)


# ---- 6: delta-only events with a surround get a finite kinetic energy ----------------
def test_delta_only_events_have_finite_ekin_and_flow_through_the_chain():
    ev = mt.sample_track_events(
        "proton", 20.0, sv_um=(20, 20, 10), surround_um=50.0, n_events=5000, seed=3
    )
    df = mx.tracks_to_dataframe(ev)
    dep = df["edep_keV"] > 0
    assert df["delta_only"].sum() > 0  # the case of the finding exists
    assert np.all(np.isfinite(df.loc[dep, "ekin_MeV"]))
    r = process_mc_ensemble(df, lambda L: np.ones_like(L), sv_thickness_um=10.0)
    ek = r["event_kinetic_MeV"]
    assert np.all(np.isfinite(ek))
    y = r["event_energies_keV"] / 10.0
    kt = compute_kappa_table(source="bragg", particle="proton")
    out = tissue_equivalence_correction(
        y, kappa_table=kt, particle_energy_MeV=ek, particle="proton"
    )
    assert np.all(np.isfinite(out))  # raised "must be finite" before


# ---- 7/8: a partially filled sigma column is refused, not silently dropped ----------
def test_partial_sigma_columns_are_refused():
    V = -np.arange(0.0, 6.0)
    C = 1e-10 / np.sqrt(1 + np.abs(V))
    txt = cv_csv(V, C, sigma_rel=0.01).splitlines()
    txt[-2] = txt[-2].rsplit(",", 1)[0] + ","  # one empty sigma_C cell
    with pytest.raises(M.MeasurementFileError, match="sigma_C_pF.*empty cell"):
        M.load_cv_csv("\n".join(txt) + "\n")
    M.load_cv_csv(cv_csv(V, C, sigma_rel=0.01))  # complete column: fine
    T = np.linspace(280, 340, 5)
    I = 1e-12 * np.exp(-1.0 / (8.617e-5 * T)) / np.exp(-1.0 / (8.617e-5 * 300))
    txt = it_csv(T, I, sigma_A=0.02 * I).splitlines()
    txt[-1] = txt[-1].rsplit(",", 1)[0] + ","
    with pytest.raises(M.MeasurementFileError, match="sigma_I_pA.*empty cell"):
        M.load_it_csv("\n".join(txt) + "\n")
    M.load_it_csv(it_csv(T, I, sigma_A=0.02 * I))


# ---- 9: non-UTF-8 files are a MeasurementFileError, not a UnicodeDecodeError -------
def test_latin1_file_is_a_measurement_error(tmp_path):
    raw = cv_csv(-np.arange(0.0, 4.0), np.full(4, 1e-10)).replace(
        "# notes:", "# T = 25 °C, µm scale\n# notes:"
    )
    data = raw.encode("latin-1")
    with pytest.raises(M.MeasurementFileError, match="not UTF-8"):
        M.load_cv_csv(data)
    p = tmp_path / "latin1.csv"
    p.write_bytes(data)
    with pytest.raises(M.MeasurementFileError, match="line 1:"):
        M.load_cv_csv(p)
    M.load_cv_csv(raw.encode("utf-8"))  # same content in UTF-8 loads


# ---- 11 / 15: CSV headers say what was used --------------------------------------------
def test_dark_current_csv_records_the_override_N_t():
    from app.components.results import to_csv_bytes

    cfg = DeviceConfig()
    res = SimResult(
        cfg,
        "dark_current_vs_T",
        np.array([300.0]),
        np.array([-1e-11]),
        {"V_bias": -20.0, "N_t": 5e12},
    )
    text = to_csv_bytes(res).decode()
    assert "N_t=5000000000000.0" in text and "# N_t_override: 5e+12" in text
    assert f"N_t={cfg.N_t}" not in text


def test_neutron_damage_csv_column_is_per_neutron():
    from app.components.results import to_csv_bytes

    res = SimResult(
        DeviceConfig(),
        "damage",
        np.array([1e12, 1e13]),
        np.array([0.99, 0.9]),
        {"V_bias": -20.0, "energy_MeV": 1.0, "particle": "neutron"},
    )
    assert "fluence_n_per_cm2,CCE" in to_csv_bytes(res).decode()


# ---- 12: proton -> neutron switch flags the cached result as stale ---------------------
def test_particle_switch_marks_damage_result_stale(monkeypatch):
    from streamlit.testing.v1 import AppTest

    def fake(cfg, **kw):
        return SimResult(
            cfg,
            "damage",
            np.geomspace(1e11, 1e13, 3),
            np.array([0.99, 0.95, 0.9]),
            {
                "V_bias": -20.0,
                "energy_MeV": 62.0,
                "particle": kw.get("particle", "proton"),
            },
        )

    def page():
        from app.workflows.radiation_damage import render

        render()

    monkeypatch.setattr(etna, "run_radiation_damage", fake)
    at = AppTest.from_function(page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    at.button[0].click()
    at.run()
    assert at.exception == []
    at.radio(key="rad_particle").set_value("neutron").run()
    assert any("particle" in (w.value or "") for w in at.warning)


# ---- 16: spec-sheet W(-10 V) is the solved value, not an interpolation --------------
def test_spec_sheet_W_at_minus_10V_is_exact():
    from etna.export.spec_sheet import build_spec

    cfg = etna.nplus_calibrated_config()
    rep = etna.run_all(cfg, sections=["cv"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = build_spec(cfg, rep)
    assert spec["cv"]["W_um"]["-10 V"] == pytest.approx(9.5, abs=5e-5)  # was 9.4723
    assert spec["cv"]["W_um"]["-30 V"] == pytest.approx(9.97883, abs=5e-5)


# ---- 17: VTU comment can never be closed or invalidated by the metadata -------------
@pytest.mark.parametrize("note", ["a ---> b", "x--", "--->", "a - -- b ----"])
def test_vtu_comment_sanitiser(note):
    from etna.export.mesh_export import _xml_comment_safe

    s = _xml_comment_safe(note)
    assert "--" not in s and not s.endswith("-")
    ET.fromstring(f"<a><!-- {s} --></a>")
