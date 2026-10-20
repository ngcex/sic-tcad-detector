# data/srim/ -- superseded

This directory used to hold three PLACEHOLDER proton stopping-power CSVs
(`water_proton.csv`, `si_proton.csv`, `c_proton.csv`) that blocked
`compute_kappa_table(source="bragg")`. They were removed on 2026-10-20, when real NIST
PSTAR/ASTAR tables were downloaded: see **`data/nist_star/README.md`** (provenance,
SHA-256 of the raw responses, Bragg-additivity caveats, anchor values).

No SRIM output is vendored. A SRIM-2013 (or ICRU 73/90) 4H-SiC table would still be a
useful cross-check of the Bragg-additivity composition; if one is added, record the
SRIM version and parameters here.

Note on the expected values previously written here ("kappa ~1.24 at 1 MeV -> ~1.13 at
100 MeV"): they were not taken from data and match neither convention of the real
tables (mass ratio 1.37 -> 1.21, LET ratio 0.43 -> 0.38).
