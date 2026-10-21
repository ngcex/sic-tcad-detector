"""Independent check of the NJOY HEATR neutron damage cross sections (elastic part).

Recomputes the elastic damage-energy cross section of Si-28 and C-12 directly from the
ENDF/B-VIII.0 file (MF3/MT2 cross section + MF4/MT2 centre-of-mass angular distribution,
Legendre and tabulated parts) with the same Lindhard/Robinson partition HEATR uses,
and compares it with the stored HEATR MT444 (``data/niel/neutron_damage_*``) where
elastic scattering is the only open channel that matters (below the first inelastic
level: Si-28 1.779 MeV, C-12 4.439 MeV). Then varies the stopping LATTICE of the
partition (NJOY: recoil in its own element; SiC: recoil in a Si/C mixture) to bound
the compound effect on the SiC/Si neutron ratio r_n.

Needs the ENDF files and the RECONR PENDF tapes of the build (not in the repo; see
scripts/build_neutron_niel.py; <work> is its --work directory):

    uv run python lab/neutron_niel/elastic_crosscheck.py <endf-dir> <work>
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from numpy.polynomial import legendre as L

from etna.core.neutron_niel import (
    _endf_float,
    isotope_damage,
    parse_mf3_section,
    watt_spectrum,
)

ED = 21.0


def lindhard_damage(E_R, Zr, Ar, Zl, Al, ed=ED):
    """NJOY HEATR ``df``: Robinson's fit of the Lindhard partition (eV in, eV out)."""
    E_R = np.asarray(E_R, float)
    el = 30.724 * Zr * Zl * np.sqrt(Zr ** (2 / 3) + Zl ** (2 / 3)) * (Ar + Al) / Al
    fl = (
        0.0793
        * Zr ** (2 / 3)
        * np.sqrt(Zl)
        * (Ar + Al) ** 1.5
        / ((Zr ** (2 / 3) + Zl ** (2 / 3)) ** 0.75 * Ar**1.5 * np.sqrt(Al))
    )
    ep = E_R / el
    d = E_R / (1 + fl * (ep + 3.4008 * ep ** (1 / 6) + 0.40244 * ep**0.75))
    return np.where(E_R < ed, 0.0, d)


def _records(text, mat, mf, mt):
    tag = f"{mat:4d}{mf:2d}{mt:3d}"
    return [ln for ln in text.splitlines() if ln[66:75] == tag]


def _floats(lines, n, start):
    """n ENDF floats starting at line index ``start``; returns (values, next line)."""
    vals, i = [], start
    while len(vals) < n:
        ln = lines[i]
        for k in range(6):
            if len(vals) < n:
                vals.append(_endf_float(ln[11 * k : 11 * k + 11]))
        i += 1
    return vals, i


def parse_mf4_elastic(text, mat):
    """MF4/MT2 (LTT 1/3, LCT 2): Legendre list [(E, a_l...)] and tabulated [(E, mu, p)]."""
    ls = _records(text, mat, 4, 2)
    head = ls[0]
    ltt = int(head[33:44])
    lct = int(ls[1][33:44])
    assert lct == 2, "expect centre-of-mass distributions"
    i = 2
    leg, tab = [], []
    if ltt in (1, 3):
        nr, ne = int(ls[i][44:55]), int(ls[i][55:66])
        i += 1 + (2 * nr + 5) // 6
        for _ in range(ne):
            E = _endf_float(ls[i][11:22])
            nl = int(ls[i][44:55])
            i += 1
            a, i = _floats(ls, nl, i)
            leg.append((E, np.asarray(a)))
    if ltt in (2, 3):
        nr, ne = int(ls[i][44:55]), int(ls[i][55:66])
        i += 1 + (2 * nr + 5) // 6
        for _ in range(ne):
            E = _endf_float(ls[i][11:22])
            nr2, np2 = int(ls[i][44:55]), int(ls[i][55:66])
            i += 1 + (2 * nr2 + 5) // 6
            v, i = _floats(ls, 2 * np2, i)
            v = np.asarray(v).reshape(np2, 2)
            tab.append((E, v[:, 0], v[:, 1]))
    return leg, tab


MU = np.linspace(-1, 1, 2001)


def _pdf_at(E, leg, tab):
    """p(mu) on MU at incident energy E (lin-lin in E between tabulated energies)."""

    def leg_pdf(a):
        c = np.concatenate([[1.0], a])
        c = c * (2 * np.arange(len(c)) + 1) / 2
        return L.legval(MU, c)

    def tab_pdf(mu, p):
        return np.interp(MU, mu, p)

    entries = [(e, ("l", a)) for e, a in leg] + [(e, ("t", (m, p))) for e, m, p in tab]
    entries.sort(key=lambda t: t[0])
    Es = np.array([e for e, _ in entries])
    j = np.searchsorted(Es, E)
    j = min(max(j, 1), len(Es) - 1)
    out = []
    for e, (kind, d) in (entries[j - 1], entries[j]):
        out.append(leg_pdf(d) if kind == "l" else tab_pdf(*d))
    w = (E - Es[j - 1]) / (Es[j] - Es[j - 1])
    w = min(max(w, 0.0), 1.0)
    pdf = (1 - w) * out[0] + w * out[1]
    pdf = np.clip(pdf, 0, None)
    return pdf / np.trapezoid(pdf, MU)


def elastic_damage(text, pendf, mat, Z, energies_eV, lattice):
    """sigma_D,el(E) [eV b] for the given lattice (Zl, Al).

    Cross section from the RECONR PENDF (resonances reconstructed), angular
    distribution from the ENDF file itself."""
    awr = _endf_float(_records(text, mat, 3, 2)[0][11:22])
    Es, sig = parse_mf3_section(pendf, mat, 2)
    leg, tab = parse_mf4_elastic(text, mat)
    out = []
    for E in energies_eV:
        pdf = _pdf_at(E, leg, tab)
        E_R = E * 2 * awr / (awr + 1) ** 2 * (1 - MU)
        d = lindhard_damage(E_R, Z, awr, lattice[0], lattice[1])
        out.append(np.interp(E, Es, sig) * np.trapezoid(pdf * d, MU))
    return np.asarray(out)


def main(endf_dir, work):
    d, w_ = Path(endf_dir), Path(work)
    si = (d / "n-014_Si_028.endf").read_text()
    c = (d / "n-006_C_012.endf").read_text()
    pe = {"Si28": (w_ / "Si28_Ed21" / "tape22").read_text(),
          "C12": (w_ / "C12_Ed21" / "tape22").read_text()}
    awr_si, awr_c = 27.737, 11.89365
    A_mix = 0.5 * (awr_si + awr_c)
    lattices = {
        "own element (NJOY)": None,
        "SiC mean lattice (Z 10, A_mix)": (10.0, A_mix),
        "other element": "swap",
    }
    print("Elastic damage (eV b), python vs NJOY HEATR MT444 (total), Ed 21 eV")
    for name, text, mat, Z, awr, E_list in (
        ("Si28", si, 1425, 14.0, awr_si, [0.2e6, 0.5e6, 0.8e6, 1.0e6, 1.2e6, 1.5e6]),
        ("C12", c, 625, 6.0, awr_c, [0.2e6, 0.5e6, 1.0e6, 2.0e6, 3.0e6, 4.0e6]),
    ):
        Eh, sh = isotope_damage(name, ED)
        mine = elastic_damage(text, pe[name], mat, Z, E_list, (Z, awr))
        for E, m in zip(E_list, mine):
            nj = np.interp(E * 1e-6, Eh, sh) * 1e30  # MeV cm^2 -> eV b
            print(
                f"  {name} {E / 1e6:5.2f} MeV  python {m:10.1f}  NJOY {nj:10.1f}  "
                f"ratio {m / nj:.4f}"
            )

    # fission-averaged elastic damage under the lattice variants (0.01-20 MeV)
    Eg = np.geomspace(1e4, 2e7, 400)
    w = watt_spectrum()(Eg * 1e-6)
    print(
        "\nWatt-fission-averaged ELASTIC damage, lattice variants (relative to NJOY's):"
    )
    res = {}
    for name, text, mat, Z, awr, other in (
        ("Si28", si, 1425, 14.0, awr_si, (6.0, awr_c)),
        ("C12", c, 625, 6.0, awr_c, (14.0, awr_si)),
    ):
        base = elastic_damage(text, pe[name], mat, Z, Eg, (Z, awr))
        bavg = np.trapezoid(base * w, Eg)
        for lname, lat in lattices.items():
            lat = (Z, awr) if lat is None else (other if lat == "swap" else lat)
            v = elastic_damage(text, pe[name], mat, Z, Eg, lat)
            r = np.trapezoid(v * w, Eg) / bavg
            res[(name, lname)] = r
            print(f"  {name:4s} in {lname:32s}: {r:.4f}")
    # SiC/Si elastic damage ratio shift: only the SiC side changes lattice
    for lname in lattices:
        if lname.startswith("own"):
            continue
        print(
            f"  -> r_n factor if both SiC sublattices use '{lname}': "
            f"Si x{res[('Si28', lname)]:.4f}, C x{res[('C12', lname)]:.4f}"
        )


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
