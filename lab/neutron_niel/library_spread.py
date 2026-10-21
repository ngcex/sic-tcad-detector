"""Evaluated-library spread of the neutron NIEL numbers: ENDF/B-VII.1 vs the stored
ENDF/B-VIII.0 tables (data/niel/). VII.1 has natural carbon (MAT 600) instead of
C-12/C-13 and independent Si evaluations of an older vintage.

Needs NJOY and the VII.1 files (NNDC, ENDF-B-VII.1-neutrons.zip):

    uv run python lab/neutron_niel/library_spread.py <njoy> <endf-b71-dir> <work>
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from etna.core import neutron_niel as nn

B71 = {  # iso -> (file, MAT, atom fraction)
    "Si28": ("n-014_Si_028.endf", 1425, 0.92223),
    "Si29": ("n-014_Si_029.endf", 1428, 0.04685),
    "Si30": ("n-014_Si_030.endf", 1431, 0.03092),
    "Cnat": ("n-006_C_000.endf", 600, 1.0),
}


def summary(E_si, si, E_sic, sic):
    w = nn.watt_spectrum()
    f_si = nn.spectrum_average(E_si, si, w, 1e-3, 20.0)
    f_sic = nn.spectrum_average(E_sic, sic, w, 1e-3, 20.0)
    return {
        "Si(1 MeV)": float(np.interp(1.0, E_si, si)),
        "SiC(1 MeV)": float(np.interp(1.0, E_sic, sic)),
        "Si fission": f_si,
        "SiC fission": f_sic,
        "r_n fission": f_sic / f_si,
        "r_n 1 MeV": float(np.interp(1.0, E_sic, sic) / np.interp(1.0, E_si, si)),
    }


def main(njoy, endf_dir, work):
    curves = {}
    for iso, (fname, mat, f) in B71.items():
        E, s, _ = nn.run_njoy(
            Path(njoy), Path(endf_dir) / fname, mat, 21.0, Path(work) / f"b71_{iso}"
        )
        curves[iso] = (E * 1e-6, s * 1e-30, f)  # MeV, MeV cm^2
    grid = np.unique(np.concatenate([c[0] for c in curves.values()]))

    def el(prefix):
        return sum(
            f * np.interp(grid, E, s, left=0, right=0)
            for iso, (E, s, f) in curves.items()
            if iso.startswith(prefix)
        )

    si = el("Si") * nn.AVOGADRO / nn.M_SI
    sic = (el("Si") + el("C")) * nn.AVOGADRO / (nn.M_SI + nn.M_C)
    b71 = summary(grid, si, grid, sic)
    E1, si8 = nn.material_niel_from_isotopes("si", nn.ED_SETS["Ed21"])
    E2, sic8 = nn.material_niel_from_isotopes("sic", nn.ED_SETS["Ed21"])
    b80 = summary(E1, si8, E2, sic8)
    print(f"{'quantity':14s} {'B-VIII.0':>11s} {'B-VII.1':>11s}  VII.1/VIII.0")
    for k in b80:
        print(f"{k:14s} {b80[k]:11.4e} {b71[k]:11.4e}  {b71[k] / b80[k]:.4f}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
