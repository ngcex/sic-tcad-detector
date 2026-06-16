# Real proton stopping-power data for Phase 27 (κ tissue-equivalence)

The files in this directory are **PLACEHOLDERS**. `compute_kappa_table(source="bragg")`
refuses to run against them (raises `FileNotFoundError`) until they are replaced
with real tabulated data — fabricating these is the exact defect physics audit
C-1 flagged.

## What to drop in (3 CSV files, same format)

Each file: header `energy_MeV,stopping_power_MeV_cm2_per_g` followed by rows of
proton **mass** stopping power (MeV·cm²/g) vs proton kinetic energy (MeV),
covering at least 0.1–250 MeV.

| File               | Source         | What to download                          |
| ------------------ | -------------- | ----------------------------------------- |
| `water_proton.csv` | NIST **PSTAR** | Liquid water, total mass stopping power   |
| `si_proton.csv`    | NIST **PSTAR** | Silicon (Z=14), total mass stopping power |
| `c_proton.csv`     | NIST **PSTAR** | Carbon (Z=6), total mass stopping power   |

NIST PSTAR: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html
(this environment has no network access — download elsewhere and copy the files in).

SiC is then composed automatically by Bragg additivity in
`sic_stopping_power_bragg()`:
`(S/ρ)_SiC = 0.7004·(S/ρ)_Si + 0.2996·(S/ρ)_C` (4H-SiC mass fractions).

## Expected sanity check once real data is in place

`compute_kappa_table(source="bragg")` should give **κ ≈ 1.24 at 1 MeV decreasing
to ≈ 1.13 at 100 MeV** (a ~10% monotonic decrease), i.e. κ > 1 everywhere —
NOT the legacy flat ~0.58. After dropping the data, update
`tests/test_microdosimetry.py` κ-range assertions to ~1.0–1.4 and switch the
microdosimetry notebooks to `source="bragg"`.

## Alternative source

SRIM (`SRIM-2013`) can produce the SiC stopping power directly; if used, write a
single `sic_proton.csv` and load it via the legacy path instead of Bragg additivity.
Document the SRIM version and parameters here.
