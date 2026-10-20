"""Validation of the bundled NIST PSTAR/ASTAR tables and the quantities derived from them.

Anchor values are read from the downloaded NIST files themselves (data/nist_star/raw,
retrieved 2026-10-20) or from the published Si microdosimetry literature (Si -> muscle
lineal-energy factor ~0.58 for protons). No number here is hand-typed physics data
other than those anchors.
"""

from pathlib import Path

import numpy as np
import pytest

from etna.core import stopping_power as sp

RAW = sp.DATA_DIR / "raw"


@pytest.mark.parametrize("particle", ["proton", "alpha"])
@pytest.mark.parametrize("material", sorted(sp.NIST_MATERIALS))
def test_csv_is_exact_reformat_of_raw_nist_response(particle, material):
    raw = RAW / f"{sp.PROGRAMS[particle]}_{sp.NIST_MATERIALS[material]}.html"
    if not raw.exists():
        pytest.skip(f"raw NIST response {raw.name} not present")
    name, table = sp.parse_nist_star_html(raw.read_text(encoding="latin-1"))
    t = sp.load_star(particle, material)
    assert len(t["energy_MeV"]) == len(table)
    for j, col in enumerate(sp.COLUMNS):
        np.testing.assert_allclose(t[col], table[:, j], rtol=1e-5)
    # header names the NIST material that was requested
    first = Path(sp._csv_path(particle, material)).read_text().splitlines()[0]
    assert name in first


def test_raw_material_names():
    expect = {
        "014": "SILICON",
        "006": "CARBON (amorphous",
        "906": "CARBON (Graphite",
        "276": "WATER, LIQUID",
        "202": "MUSCLE, STRIATED (ICRU)",
    }
    for prog in ("pstar", "astar"):
        for matno, nm in expect.items():
            name, _ = sp.parse_nist_star_html(
                (RAW / f"{prog}_{matno}.html").read_text(encoding="latin-1")
            )
            assert name.startswith(nm), (prog, matno, name)


def test_alpha_5p5MeV_range_in_silicon_nist_column():
    """ASTAR Si: 5.5 MeV alpha CSDA range 6.487e-3 g/cm^2 = 27.84 um (~28 um anchor)."""
    t = sp.load_star("alpha", "si")
    i = int(np.argmin(abs(t["energy_MeV"] - 5.5)))
    assert t["energy_MeV"][i] == pytest.approx(5.5)
    assert t["csda_range_g_cm2"][i] == pytest.approx(6.487e-3, rel=1e-4)
    um = t["csda_range_g_cm2"][i] / sp.density("si") * 1e4
    assert 27.5 < um < 28.2


@pytest.mark.parametrize(
    "particle,E", [("alpha", 5.5), ("proton", 10.0), ("proton", 62.0)]
)
def test_integrated_csda_reproduces_nist_silicon(particle, E):
    """Our 1/S integration (used for SiC, which NIST lacks) reproduces NIST's own range."""
    t = sp.load_star(particle, "si")
    nist = np.exp(
        np.interp(np.log(E), np.log(t["energy_MeV"]), np.log(t["csda_range_g_cm2"]))
    )
    e, r = sp._csda_table(particle, "si")  # the integration path used for SiC
    ours = np.interp(E, e, r)
    assert ours == pytest.approx(nist, rel=2e-3)


def test_proton_62MeV_range_in_water():
    """PSTAR water: 62 MeV CSDA range 3.281 g/cm^2 (log-log between the 60/65 MeV rows; the old 31.0 mm table entry was wrong)."""
    assert sp.csda_range_cm("proton", "water", 62.0) * 10 == pytest.approx(
        32.81, abs=0.01
    )


def test_si_to_muscle_linear_ratio_matches_literature_factor():
    """Si microdosimetry literature: ~0.58 Si -> striated muscle for protons (linear ratio)."""
    for E in (10.0, 62.0, 150.0):
        r = sp.linear_stopping("proton", "muscle_icru", E) / sp.linear_stopping(
            "proton", "si", E
        )
        assert 0.54 < r < 0.59  # NIST: 0.584 @10, 0.556 @62, 0.550 @150 MeV


def test_sic_bragg_additivity_between_constituents():
    E = np.array([1.0, 10.0, 100.0])
    s_si = sp.mass_stopping("proton", "si", E)
    s_c = sp.mass_stopping("proton", "c_amorphous", E)
    s = sp.sic_mass_stopping("proton", E)
    np.testing.assert_allclose(s, 0.7004 * s_si + 0.2996 * s_c, rtol=1e-12)
    assert np.all((s > np.minimum(s_si, s_c)) & (s < np.maximum(s_si, s_c)))


def test_carbon_phase_choice_is_a_small_effect_above_1MeV():
    E = np.geomspace(1.0, 250.0, 30)
    a = sp.sic_mass_stopping("proton", E, carbon="c_amorphous")
    g = sp.sic_mass_stopping("proton", E, carbon="c_graphite")
    assert np.max(np.abs(a / g - 1)) < 0.01


class TestKappa:
    def test_linear_proton_kappa_values_and_trend(self):
        k = sp.kappa_table("proton")
        E = np.array([1.0, 10.0, 62.0, 100.0, 150.0, 250.0])
        v = np.exp(np.interp(np.log(E), np.log(k["energy_MeV"]), np.log(k["kappa"])))
        assert np.all(np.diff(v) < 0)  # monotonic decrease with kinetic energy
        np.testing.assert_allclose(v[:4], [0.426, 0.391, 0.378, 0.376], atol=2e-3)
        # NOT flat: ~13% variation 1 -> 250 MeV (the legacy table was flat to ~2%)
        assert v[0] / v[-1] > 1.12

    def test_mass_ratio_is_density_times_linear(self):
        lin = sp.kappa_table("proton", ratio="linear")["kappa"]
        mass = sp.kappa_table("proton", ratio="mass")["kappa"]
        np.testing.assert_allclose(mass, lin * sp.density("sic") / 1.0, rtol=1e-12)
        E = sp.kappa_table("proton")["energy_MeV"]
        sel = (E >= 1.0) & (E <= 250.0)
        assert np.all((mass[sel] > 1.19) & (mass[sel] < 1.37))

    def test_alpha_kappa_differs_from_proton(self):
        kp = sp.kappa_table("proton")
        ka = sp.kappa_table("alpha")
        E = 5.0
        vp = np.interp(E, kp["energy_MeV"], kp["kappa"])
        va = np.interp(E, ka["energy_MeV"], ka["kappa"])
        assert abs(va - vp) > 0.01
        assert ka["particle"] == "alpha"

    def test_bad_arguments(self):
        with pytest.raises(ValueError):
            sp.kappa_table("electron")
        with pytest.raises(ValueError):
            sp.kappa_table("proton", tissue="bone")
        with pytest.raises(ValueError):
            sp.kappa_table("proton", ratio="dose")
        with pytest.raises(ValueError):
            sp.mass_stopping("proton", "si", 2e4)  # beyond the NIST grid


class TestAlphaBraggCurve:
    def test_am241_alpha_numbers_in_sic(self):
        s = sp.alpha_bragg_summary()
        # SiC range: between the Si-mass-range and C-mass-range scaled values, ~18.4 um
        assert s["csda_range_um"] == pytest.approx(18.43, abs=0.1)
        assert 0.1 < s["bohr_straggling_um"] < 0.25
        assert 2.3 < s["peak_to_entrance"] < 2.7
        # peak close to the end of range (within ~2 um)
        assert s["csda_range_um"] - 2.0 < s["peak_depth_um"] < s["csda_range_um"]

    def test_curve_integrates_to_deposited_energy(self):
        s, f, R, _ = sp.alpha_bragg_curve_sic(5.486)
        # electronic part only: total minus the (small) nuclear share
        dep = np.trapezoid(f, s)
        assert 5.40 < dep < 5.486

    def test_residual_energy_monotonic_and_zero_beyond_range(self):
        R = sp.csda_range_cm("alpha", "sic", 5.486)
        d = np.linspace(0, 1.2 * R, 50)
        E = sp.residual_energy_MeV("alpha", "sic", 5.486, d)
        assert E[0] == pytest.approx(5.486, rel=1e-6)
        assert np.all(np.diff(E) <= 0)
        assert np.all(E[d > R * 1.001] == 0.0)
