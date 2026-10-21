"""Compute neutron NIEL (displacement damage energy) for Si and SiC with NJOY2016 HEATR.

Regenerates ``data/niel/neutron_damage_*.csv`` (pointwise damage-energy production
cross sections per isotope, ENDF MF3/MT444 written by HEATR) and
``data/niel/neutron_niel_*.csv`` (natural-abundance Si and SiC NIEL, lethargy-averaged
on a fixed grid). Provenance and the method are in ``data/niel/README.md``; nothing
here is typed by hand.

Needs (not ETNA dependencies, nothing in the package calls this at run time):

- NJOY2016 (open source, https://github.com/njoy/NJOY2016). The numbers in the repo were
  made with conda-forge ``njoy2016`` 2016.78::

      conda create -p <scratch>/njoyenv -c conda-forge --override-channels njoy2016

- ENDF/B-VIII.0 neutron files for Si-28/29/30 and C-12/13 (NNDC,
  https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_neutrons.zip; the SHA-256 of
  each file used is checked against ``ENDF_FILES`` below).

Usage::

    uv run python scripts/build_neutron_niel.py --njoy <path/to/njoy> \
        --endf-dir <dir with the n-0*.endf files> --work <scratch dir>
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path

import numpy as np

from etna.core.neutron_niel import (
    ENDF_FILES,
    ED_SETS,
    NJOY_VERSION,
    NIEL_GRID_MEV,
    damage_csv_name,
    material_niel_from_isotopes,
    niel_csv_name,
    read_text,
    run_njoy,
    lethargy_average,
)

OUT = Path(__file__).resolve().parent.parent / "data" / "niel"


def njoy_version_line(work: Path) -> str:
    for line in (work / "output.log").read_text().splitlines():
        if "njoy" in line.lower() and ("2016" in line or "version" in line.lower()):
            # "njoy 2016.78  03Feb25   <run date/time>": keep the version only, so the
            # output files are byte-reproducible
            return " ".join(line.split()[:3])
    return "unknown"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--njoy", type=Path, required=True)
    ap.add_argument("--endf-dir", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    args = ap.parse_args(argv)

    for iso, (fname, mat, _abund, sha) in ENDF_FILES.items():
        got = hashlib.sha256((args.endf_dir / fname).read_bytes()).hexdigest()
        if got != sha:
            raise SystemExit(f"{fname}: sha256 {got} != expected {sha}")

    done = {}
    for tag, eds in ED_SETS.items():
        for iso, (fname, mat, _abund, sha) in ENDF_FILES.items():
            ed = eds["Si"] if iso.startswith("Si") else eds["C"]
            key = (iso, ed)
            if key in done:
                continue
            work = args.work / f"{iso}_Ed{ed:g}"
            E, sig, inp = run_njoy(args.njoy, args.endf_dir / fname, mat, ed, work)
            ver = njoy_version_line(work)
            out = OUT / damage_csv_name(iso, ed)
            lines = [
                f"# NJOY HEATR damage-energy production cross section (MF3/MT444), {iso}",
                f"# ENDF/B-VIII.0 file {fname} (MAT {mat}) sha256={sha}",
                f"# NJOY: {NJOY_VERSION} ({ver}); RECONR err 0.001, 0 K; HEATR ed = {ed:g} eV",
                "# NJOY input: " + " | ".join(inp.strip().splitlines()),
                "# partition: NJOY HEATR Lindhard/Robinson (recoil in its own element)",
                "energy_eV,damage_eV_barn",
            ]
            lines += [f"{e:.7e},{s:.7e}" for e, s in zip(E, sig)]
            # gzip with mtime 0: byte-reproducible for the same NJOY output
            out.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0))
            done[key] = out
            print(f"wrote {out.name}: {len(E)} points, {ver}")

    for tag, eds in ED_SETS.items():
        for material in ("si", "sic"):
            Eg = NIEL_GRID_MEV
            pts = material_niel_from_isotopes(material, eds, OUT)
            avg = lethargy_average(pts[0], pts[1], Eg)
            out = OUT / niel_csv_name(material, tag)
            hashes = " ".join(
                f"{p.name}={hashlib.sha256(read_text(p).encode()).hexdigest()[:16]}"
                for (iso, ed), p in sorted(done.items())
                if (material == "sic" or iso.startswith("Si"))
                and ed == (eds["Si"] if iso.startswith("Si") else eds["C"])
            )
            lines = [
                f"# Neutron NIEL in {material.upper()} (natural abundance), NJOY HEATR on ENDF/B-VIII.0",
                f"# displacement thresholds [eV]: Si {eds['Si']:g}"
                + (f", C {eds['C']:g}" if material == "sic" else ""),
                "# value = lethargy average of the pointwise NIEL over [E_lo, E_hi] (MeV);",
                "#   energy_MeV = geometric bin centre. Pointwise inputs (sha256[:16]): "
                + hashes,
                "energy_MeV,E_lo_MeV,E_hi_MeV,niel_MeV_cm2_per_g",
            ]
            for i in range(len(Eg) - 1):
                c = float(np.sqrt(Eg[i] * Eg[i + 1]))
                lines.append(f"{c:.6e},{Eg[i]:.6e},{Eg[i + 1]:.6e},{avg[i]:.6e}")
            out.write_text("\n".join(lines) + "\n")
            print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
