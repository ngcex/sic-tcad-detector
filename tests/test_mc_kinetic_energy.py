"""Audit C-2: MC events carry the primary KINETIC energy for the kappa(E) lookup."""

import numpy as np
import pandas as pd
import pytest

from etna.core.mc_coupling import load_mc_events_csv, process_mc_ensemble
from etna.core.microdosimetry import compute_kappa_table, tissue_equivalence_correction


def _cce_one(let):
    return 1.0


def _events():
    # event 0: primary (track 1) + a secondary electron (track 2) listed first
    # event 1: zero deposit (filtered); event 2: primary only
    return pd.DataFrame(
        {
            "event_id": [0, 0, 0, 1, 2, 2],
            "x": [0.0] * 6,
            "y": [0.0] * 6,
            "z": [0.0] * 6,
            "edep": [5.0, 30.0, 20.0, 0.0, 40.0, 10.0],
            "ekin": [0.02, 62.0, 61.99, 10.0, 10.0, 9.99],
            "track_id": [2, 1, 1, 1, 1, 1],
        }
    )


def test_csv_loader_maps_ekin_and_track_id(tmp_path):
    p = tmp_path / "ev.csv"
    _events().to_csv(p, index=False)
    df = load_mc_events_csv(p, energy_unit="keV")
    assert "ekin_MeV" in df.columns and "track_id" in df.columns
    # ekin is never unit-converted (MeV by name), even when edep is converted
    df2 = load_mc_events_csv(p, energy_unit="MeV")
    np.testing.assert_allclose(df2["ekin_MeV"], df["ekin_MeV"])
    np.testing.assert_allclose(df2["edep_keV"], 1e3 * df["edep_keV"])


def test_event_kinetic_energy_is_primary_first_step_and_masked(tmp_path):
    p = tmp_path / "ev.csv"
    _events().to_csv(p, index=False)
    res = process_mc_ensemble(load_mc_events_csv(p), _cce_one)
    # event 1 filtered (zero deposit); event 0 -> primary's first step 62.0, not 0.02
    np.testing.assert_allclose(res["event_energies_keV"], [55.0, 50.0])
    np.testing.assert_allclose(res["event_kinetic_MeV"], [62.0, 10.0])


def test_without_ekin_column_returns_none():
    df = _events().rename(columns={"edep": "edep_keV"}).drop(columns=["ekin"])
    res = process_mc_ensemble(df, _cce_one)
    assert res["event_kinetic_MeV"] is None


def test_end_to_end_kinetic_kappa(tmp_path):
    p = tmp_path / "ev.csv"
    _events().to_csv(p, index=False)
    res = process_mc_ensemble(load_mc_events_csv(p), _cce_one)
    y = res["event_collected_keV"] / 20.0
    kt = compute_kappa_table(source="bragg", particle="proton")
    yt = tissue_equivalence_correction(
        y,
        kappa_table=kt,
        particle_energy_MeV=res["event_kinetic_MeV"],
        particle="proton",
    )
    k = yt / y
    assert k[0] == pytest.approx(0.378, abs=2e-3)  # 62 MeV
    assert k[1] == pytest.approx(0.391, abs=2e-3)  # 10 MeV
