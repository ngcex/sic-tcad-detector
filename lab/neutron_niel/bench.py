"""Neutron NIEL / SiC normalisation numbers from the stored tables (seconds, no NJOY).

    uv run python lab/neutron_niel/bench.py

Prints the checks of data/niel/README.md (Si vs ASTM E722 and vs the Griffin damage
function), the SiC values, the spectrum sensitivity of r_n, and what the opt-in
"sic-damage" convention changes in kappa and in the compensation fluence. The NJOY
reruns (elastic_crosscheck.py, library_spread.py) need NJOY + ENDF files.
"""

import time

import numpy as np

from etna.core import neutron_niel as nn
from etna.core import niel_tables as nt
from etna.core.radiation_damage import (
    get_hardness_factor,
    hardness_table,
    niel_ratio_sic_si,
)


def main():
    t0 = time.perf_counter()
    ed = nn.ED_SETS["Ed21"]
    E1, si = nn.material_niel_from_isotopes("si", ed)
    E2, sic = nn.material_niel_from_isotopes("sic", ed)
    astm = nt.NIEL_NEUTRON_1MEV_SI_MEV_CM2_G
    w = nn.watt_spectrum()

    print("== Si checks")
    print(
        f"  NIEL_n,Si(1 MeV) HEATR pointwise / ASTM 95 MeV mb : "
        f"{np.interp(1.0, E1, si) / astm:.4f}"
    )
    g = nn._load_numeric(nn.DATA_DIR / "si_neutron_damage_griffin_rd48.csv")
    Eg, Dg = g[:, 0], g[:, 1] * astm
    fa = nn.spectrum_average(E1, si, w, 1e-3, 19.9)
    fb = nn.spectrum_average(Eg, Dg, w, 1e-3, 19.9)
    print(
        f"  Watt-fission hardness kappa_Si: HEATR {fa / astm:.4f}, Griffin "
        f"{fb / astm:.4f} (ratio {fa / fb:.4f})"
    )
    edges = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 19.9])
    r = nn.lethargy_average(E1, si, edges) / nn.lethargy_average(Eg, Dg, edges)
    print(
        "  HEATR/Griffin per lethargy bin: "
        + ", ".join(
            f"{lo:g}-{hi:g} {v:.3f}" for lo, hi, v in zip(edges[:-1], edges[1:], r)
        )
    )

    print("== SiC (4H, 1:1, natural abundance), E_d 21/21 eV")
    p1 = np.interp(1.0, E2, sic)
    per_atom_mbmev = p1 * (nn.M_SI + nn.M_C) / nn.AVOGADRO / 2 * 1e27
    print(
        f"  NIEL_n,SiC(1 MeV) = {p1:.4e} MeV cm2/g ({per_atom_mbmev:.1f} MeV mb per "
        f"atom); SiC/Si at 1 MeV = {p1 / np.interp(1.0, E1, si):.4f}"
    )
    print("  r_n = <SiC>/<Si> for spectrum variants:")
    for lab, lo, hi in (
        ("Watt U-235, 1 keV-20 MeV (default)", 1e-3, 20.0),
        ("Watt, 0.1-20 MeV", 0.1, 20.0),
        ("Watt, 1-20 MeV", 1.0, 20.0),
    ):
        rn = nn.spectrum_average(E2, sic, w, lo, hi) / nn.spectrum_average(
            E1, si, w, lo, hi
        )
        print(f"    {lab:36s} {rn:.4f}")
    flat = lambda E: 1.0 / E  # noqa: E731  (1/E epithermal slowing-down shape)
    rn_ep = nn.spectrum_average(E2, sic, flat, 1e-3, 0.1) / nn.spectrum_average(
        E1, si, flat, 1e-3, 0.1
    )
    print(f"    1/E epithermal part only, 1-100 keV     {rn_ep:.4f}")
    # reactor-like mix: Watt fission + a 1/E slowing-down tail (1 keV-1 MeV) carrying a
    # fraction f of the >1 keV flux (illustrative shapes, NOT a measured TRIGA spectrum)
    x = np.geomspace(1e-3, 20.0, 40001)
    ww = w(x) / np.trapezoid(w(x), x)
    tail = np.where(x < 1.0, 1.0 / x, 0.0)
    tail /= np.trapezoid(tail, x)
    for f in (0.3, 0.5, 0.7):
        phi = (1 - f) * ww + f * tail
        num = np.trapezoid(np.interp(x, E2, sic) * phi, x)
        den = np.trapezoid(np.interp(x, E1, si) * phi, x)
        print(f"    Watt + 1/E tail, tail flux fraction {f:.1f}     {num / den:.4f}")
    print(f"  embedded R_N_SIC_SI_FISSION = {nt.R_N_SIC_SI_FISSION:.4f}")

    print("== kappa per proton (neq/p): si-neq (default) -> sic-damage")
    dmg = hardness_table("proton", "sic-damage")
    for E in (5.6, 30.0, 62.0, 150.0, 250.0):
        k0, k1 = get_hardness_factor(E), get_hardness_factor(E, dmg)
        print(
            f"  {E:6.1f} MeV: {k0:.3f} -> {k1:.3f} (r_p {niel_ratio_sic_si(E):.3f}, "
            f"x{k1 / k0:.3f}; compensation fluence x{k0 / k1:.3f})"
        )
    nk = hardness_table("neutron", "si-neq")
    print(
        "== neutron kappa (si-neq, bin-averaged): "
        + ", ".join(
            f"{E:g} MeV {get_hardness_factor(E, nk):.3f}" for E in (0.5, 1.0, 2.5, 14.0)
        )
    )
    print(f"[{time.perf_counter() - t0:.2f} s]")


if __name__ == "__main__":
    main()
