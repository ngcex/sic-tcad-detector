"""Neutron NIEL (NJOY HEATR / ENDF-B-VIII.0) tables, the constants derived from them and
the opt-in NIEL conventions of the radiation-damage chain (data/niel/README.md)."""

import hashlib
import re

import numpy as np
import pytest

from etna.core import neutron_niel as nn
from etna.core import niel_tables as nt
from etna.core.radiation_damage import (
    NIEL_HARDNESS_PROTON_SIC,
    get_hardness_factor,
    hardness_table,
    niel_ratio_sic_si,
)

ED21 = nn.ED_SETS["Ed21"]


@pytest.fixture(scope="module")
def curves():
    return {m: nn.material_niel_from_isotopes(m, ED21) for m in ("si", "sic")}


def test_endf_float_format():
    assert nn._endf_float(" 1.234567+5") == pytest.approx(1.234567e5)
    assert nn._endf_float("-2.5-3") == pytest.approx(-2.5e-3)
    assert nn._endf_float(" 6012.00000") == pytest.approx(6012.0)
    assert nn._endf_float("") == 0.0


def test_parse_mf3_section_roundtrip():
    # hand-built MF3 TAB1 of 3 points, one lin-lin range (ENDF 80-column layout)
    def rec(body, mat=1425, mf=3, mt=444):
        return f"{body:66s}{mat:4d}{mf:2d}{mt:3d}{1:5d}"

    lines = [
        rec(" 1.402800+4 2.773700+1          0          0          0          0"),
        rec(" 0.000000+0 0.000000+0          0          0          1          3"),
        rec("          3          2"),
        rec(" 1.000000-5 1.000000+2 1.000000+6 9.900000+4 2.000000+7 3.000000+5"),
    ]
    E, s = nn.parse_mf3_section("\n".join(lines), 1425, 444)
    assert list(E) == pytest.approx([1e-5, 1e6, 2e7])
    assert list(s) == pytest.approx([100.0, 9.9e4, 3e5])


def test_lethargy_average_is_exact_for_linear_curves():
    E = np.array([1.0, 2.0, 4.0])
    y = 3.0 * E + 1.0
    edges = np.array([1.0, 1.5, 4.0])
    got = nn.lethargy_average(E, y, edges)
    for (lo, hi), g in zip(zip(edges[:-1], edges[1:]), got):
        exact = (3.0 * (hi - lo) + np.log(hi / lo)) / np.log(hi / lo)
        assert g == pytest.approx(exact, rel=1e-12)


@pytest.mark.parametrize("tag", list(nn.ED_SETS))
def test_averaged_tables_match_the_pointwise_data(tag):
    eds = nn.ED_SETS[tag]
    for material in ("si", "sic"):
        path = nn.DATA_DIR / nn.niel_csv_name(material, tag)
        stored = nn._load_numeric(path)
        E, y = nn.material_niel_from_isotopes(material, eds)
        avg = nn.lethargy_average(E, y, nn.NIEL_GRID_MEV)
        assert np.allclose(stored[:, 3], avg, rtol=2e-6, atol=0)
        assert np.allclose(stored[:, 1], nn.NIEL_GRID_MEV[:-1], rtol=1e-6)
        # the header lists the sha256[:16] of every pointwise input it used
        head = path.read_text().splitlines()[3]
        for name, h in re.findall(r"(neutron_damage_\S+?)=([0-9a-f]{16})", head):
            text = nn.read_text(nn.DATA_DIR / name)
            assert hashlib.sha256(text.encode()).hexdigest()[:16] == h


def test_pointwise_files_record_the_endf_hash_and_njoy_version():
    for iso, (fname, mat, _f, sha) in nn.ENDF_FILES.items():
        head = nn.read_text(nn.DATA_DIR / nn.damage_csv_name(iso, 21.0)).splitlines()
        assert fname in head[1] and f"MAT {mat}" in head[1] and sha in head[1]
        assert "njoy 2016.78" in head[2] and "ed = 21 eV" in head[2]


def test_embedded_constants_are_recomputed_from_the_files(curves):
    (E1, si), (E2, sic) = curves["si"], curves["sic"]
    w = nn.watt_spectrum()
    assert np.interp(1.0, E1, si) == pytest.approx(
        nt.NIEL_NEUTRON_1MEV_SI_NJOY_MEV_CM2_G, rel=2e-5
    )
    assert np.interp(1.0, E2, sic) == pytest.approx(
        nt.NIEL_NEUTRON_1MEV_SIC_MEV_CM2_G, rel=2e-5
    )
    f_si = nn.spectrum_average(E1, si, w, 1e-3, 20.0)
    f_sic = nn.spectrum_average(E2, sic, w, 1e-3, 20.0)
    assert f_si == pytest.approx(nt.NIEL_NEUTRON_FISSION_SI_MEV_CM2_G, rel=2e-5)
    assert f_sic == pytest.approx(nt.NIEL_NEUTRON_FISSION_SIC_MEV_CM2_G, rel=2e-5)
    assert nt.R_N_SIC_SI_FISSION == pytest.approx(1.0687, abs=1e-4)


def test_silicon_reproduces_the_astm_normalisation(curves):
    # ASTM E722: D(1 MeV) = 95 MeV mb is a defined value; HEATR gives 0.998 of it
    E, si = curves["si"]
    ratio = np.interp(1.0, E, si) / nt.NIEL_NEUTRON_1MEV_SI_MEV_CM2_G
    assert ratio == pytest.approx(1.0, abs=0.005)


def test_silicon_matches_the_griffin_damage_function(curves):
    """Griffin et al. SAND92-0094 (RD48 compilation, data/niel/): fission-spectrum
    average within 1 %, lethargy bins 0.1-19.9 MeV within 6 %."""
    g = nn._load_numeric(nn.DATA_DIR / "si_neutron_damage_griffin_rd48.csv")
    Eg, Dg = g[:, 0], g[:, 1] * nt.NIEL_NEUTRON_1MEV_SI_MEV_CM2_G
    E, si = curves["si"]
    w = nn.watt_spectrum()
    a = nn.spectrum_average(E, si, w, 1e-3, 19.9)
    b = nn.spectrum_average(Eg, Dg, w, 1e-3, 19.9)
    assert a / b == pytest.approx(1.0, abs=0.01)
    edges = np.geomspace(0.1, 19.9, 12)
    r = nn.lethargy_average(E, si, edges) / nn.lethargy_average(Eg, Dg, edges)
    assert np.all(np.abs(r - 1) < 0.06), r


def test_displacement_threshold_does_not_matter_for_fast_neutrons():
    E21, s21 = nn.material_niel_from_isotopes("sic", nn.ED_SETS["Ed21"])
    E35, s35 = nn.material_niel_from_isotopes("sic", nn.ED_SETS["Ed35-21"])
    for e in (0.1, 1.0, 14.0):
        assert np.interp(e, E35, s35) == pytest.approx(np.interp(e, E21, s21), rel=1e-4)


def test_sic_heavier_damage_per_gram_than_si_for_fission_neutrons():
    # r_n > r_p at every proton energy: the si-neq kappa (the default until 2026-10-21)
    # overstates the SiC damage
    for e in (1.0, 5.6, 62.0, 150.0, 1000.0):
        assert niel_ratio_sic_si(e) < nt.R_N_SIC_SI_FISSION


def test_hardness_table_default_is_sic_damage():
    # re-pinned 2026-10-21 (issue #8 decision): default si-neq -> sic-damage
    from etna.core.radiation_damage import DEFAULT_NIEL_CONVENTION, table_convention

    assert DEFAULT_NIEL_CONVENTION == "sic-damage"
    assert hardness_table() is hardness_table("proton", "sic-damage")
    assert get_hardness_factor(62.0) == pytest.approx(1.552, abs=1e-3)
    # the previous default stays selectable, unchanged
    assert hardness_table("proton", "si-neq") is NIEL_HARDNESS_PROTON_SIC
    assert get_hardness_factor(62.0, NIEL_HARDNESS_PROTON_SIC) == pytest.approx(
        1.942, abs=1e-3
    )
    assert table_convention(None) == "sic-damage"
    assert table_convention(NIEL_HARDNESS_PROTON_SIC) == "si-neq"
    assert table_convention(hardness_table("neutron", "sic-1mev")) == "sic-1mev"
    assert table_convention({1.0: 1.0, 10.0: 1.0}) is None
    with pytest.raises(ValueError, match="convention"):
        hardness_table("proton", "si")


@pytest.mark.parametrize("E", [1.0, 5.6, 62.0, 250.0])
def test_proton_sic_damage_convention(E):
    k = get_hardness_factor(E, hardness_table("proton", "sic-damage"))
    si = hardness_table("proton", "si-neq")
    expect = get_hardness_factor(E, si) * niel_ratio_sic_si(E) / nt.R_N_SIC_SI_FISSION
    assert k == pytest.approx(expect, rel=1e-3)


def test_proton_sic_1mev_convention():
    k = get_hardness_factor(62.0, hardness_table("proton", "sic-1mev"))
    assert k == pytest.approx(
        nt.SR_NIEL_PROTON_SIC[[e for e, _ in nt.SR_NIEL_PROTON_SIC].index(62.0)][1]
        / nt.NIEL_NEUTRON_1MEV_SIC_MEV_CM2_G,
        rel=1e-9,
    )


def test_neutron_hardness_tables():
    si = hardness_table("neutron", "si-neq")
    # 14 MeV Si hardness: Griffin 1.79-1.82 (RD48 table); the HEATR/ENDF-B-VIII.0 bin
    # gives 1.89, +4-5 % (the 10-20 MeV systematic quoted in data/niel/README.md)
    assert get_hardness_factor(14.0, si) == pytest.approx(1.80, rel=0.06)
    dmg = hardness_table("neutron", "sic-damage")
    # fission-spectrum average of the sic-damage kappa equals the Si one by construction
    E = np.array(sorted(dmg))
    w = nn.watt_spectrum()
    k_si = nn.spectrum_average(E, np.array([si[e] for e in E]), w, 1e-3, 18.0)
    k_sic = nn.spectrum_average(E, np.array([dmg[e] for e in E]), w, 1e-3, 18.0)
    assert k_sic == pytest.approx(k_si, rel=0.01)
    with pytest.raises(ValueError):
        get_hardness_factor(25.0, dmg)  # beyond the ENDF range: no extrapolation


def test_hardness_table_rejects_unknown_names():
    with pytest.raises(ValueError, match="particle"):
        hardness_table("pion")
    with pytest.raises(ValueError, match="convention"):
        hardness_table("proton", "gaas")


def test_api_compensation_fluence_scales_with_the_convention():
    from etna import DeviceConfig
    from etna.api.simulation import compensation_fluence

    # legacy scalar model (pinned: 3.48e12 -> 4.36e12 p/cm^2 at 62 MeV)
    cfg = DeviceConfig(radiation_traps=False)
    base = compensation_fluence(cfg, 62.0, niel_convention="si-neq")
    sic = compensation_fluence(cfg, 62.0)  # default: sic-damage
    assert base == pytest.approx(3.484e12, rel=1e-3)
    assert sic == pytest.approx(4.358e12, rel=1e-3)
    assert sic / base == pytest.approx(
        nt.R_N_SIC_SI_FISSION / niel_ratio_sic_si(62.0), rel=1e-3
    )
    si = hardness_table("proton", "si-neq")
    n = compensation_fluence(cfg, 1.0, particle="neutron", niel_convention="si-neq")
    k = get_hardness_factor(1.0, hardness_table("neutron", "si-neq"))
    assert n == pytest.approx(base * get_hardness_factor(62.0, si) / k, rel=1e-9)


def test_damage_figures_label_neutron_runs():
    from app.components.results import build_damage_figure
    from etna import DeviceConfig, SimResult

    r = SimResult(
        config=DeviceConfig(),
        sim_type="damage",
        x=np.geomspace(1e11, 1e12, 3),
        y=np.array([0.99, 0.98, 0.97]),
        metadata={"V_bias": -40.0, "energy_MeV": 14.0, "particle": "neutron"},
        mesh=None,
    )
    fig = build_damage_figure(r)
    assert "Neutron" in fig.layout.xaxis.title.text
