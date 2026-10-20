"""The embedded SR-NIEL tables (etna.core.niel_tables) are a verbatim transcription of
the raw data files in data/niel/ (provenance: data/niel/README.md)."""

import csv
from pathlib import Path

import pytest

from etna.core.niel_tables import (
    NIEL_NEUTRON_1MEV_SI_MEV_CM2_G,
    SR_NIEL_PROTON_SI,
    SR_NIEL_PROTON_SIC,
)

DATA = Path(__file__).parent.parent / "data" / "niel"


def _csv_rows(name):
    with open(DATA / name) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    rows = list(csv.DictReader(lines))
    return [
        (float(r["energy_MeV"]), float(r["niel_MeV_cm2_per_g"]))
        for r in rows
        if float(r["energy_MeV"]) >= 1.0
    ]


@pytest.mark.parametrize(
    "table, name",
    [
        (SR_NIEL_PROTON_SI, "sr_niel_proton_si_Ed21.csv"),
        (SR_NIEL_PROTON_SIC, "sr_niel_proton_sic_Ed21-21.csv"),
    ],
)
def test_embedded_table_equals_the_raw_csv(table, name):
    assert list(table) == _csv_rows(name)


def test_csv_headers_record_the_sr_niel_version_and_thresholds():
    head = (DATA / "sr_niel_proton_sic_Ed35-21.csv").read_text().splitlines()[:2]
    assert "SR-NIEL" in head[0] and "ver. 11.1" in head[0]
    assert "35.00; 21.00" in head[1]


def test_astm_silicon_normalisation():
    # 95 MeV mb; 100 MeV mb = 2.144 keV cm^2/g (Fretwurst et al., physics/0211118)
    assert NIEL_NEUTRON_1MEV_SI_MEV_CM2_G == pytest.approx(0.95 * 2.144e-3, rel=1e-3)
