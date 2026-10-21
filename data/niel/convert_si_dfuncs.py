"""Convert the Si neutron damage function compilation (Vasilescu & Lindstroem) to CSV.

Source file: ``neutrons.xls`` of the ROSE/RD48 "Displacement damage in silicon, on-line
compilation" by A. Vasilescu and G. Lindstroem (archived copy of
http://sesam.desy.de/members/gunnar/Si-dfuncs.html, retrieved 2026-10-21 from
https://web.archive.org/web/20091008104308/http://sesam.desy.de/members/gunnar/Si-dfuncs.html).
The page states for neutrons: "P.J. Griffin et al., SAND92-0094 ... E = 1.025E-10 -
1.995E+01 MeV". Only that first (Griffin) column block is written; the values are the
compilation's "D/(95MeVmb)" verbatim.

    uv run --with xlrd python data/niel/convert_si_dfuncs.py <neutrons.xls>
"""

import hashlib
import sys
from pathlib import Path

import xlrd


def main(src):
    src = Path(src)
    sha = hashlib.sha256(src.read_bytes()).hexdigest()
    sh = xlrd.open_workbook(str(src)).sheet_by_index(0)
    rows = []
    for r in range(sh.nrows):
        e, d = sh.cell_value(r, 1), sh.cell_value(r, 2)
        if isinstance(e, float) and isinstance(d, float):
            rows.append((e, d))
    out = Path(__file__).resolve().parent / "si_neutron_damage_griffin_rd48.csv"
    lines = [
        "# Si neutron displacement damage function D(E)/(95 MeV mb), Griffin et al.",
        "# SAND92-0094 block of the Vasilescu & Lindstroem RD48 compilation (neutrons.xls)",
        "# archived: https://web.archive.org/web/20091008104308/http://sesam.desy.de/"
        "members/gunnar/Si-dfuncs.html",
        f"# neutrons.xls sha256={sha}; converted by convert_si_dfuncs.py (verbatim cells)",
        "energy_MeV,D_over_95MeVmb",
    ]
    lines += [f"{e!r},{d!r}" for e, d in rows]
    out.write_text("\n".join(lines) + "\n")
    print(out.name, len(rows), "rows")


if __name__ == "__main__":
    main(sys.argv[1])
