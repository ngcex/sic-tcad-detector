---
id: A-11
title: Dead / test-only code and unused locals in the core slice
severity: low
category: dead-code
status: fixed 2026-10-18 (wp-c2 8136768: unused import + kT_eV locals removed; f256d10: mobility_caughey_thomas kept as the 300 K reference, depletion_width_from_capacitance kept as public utility, M_c/B/E_D_* marked informational, NC_300/NV_300 pinned to compute_ni(300) within 1 %; 8c8dec4: effective_dos kept for A-02; skip_init via A-05)
location: etna/core/sic_material.py (mobility_caughey_thomas ~156-187, effective_dos ~308-328, fields M_c, E_D_hex, E_D_cub, B, NC_300, NV_300); etna/core/device.py:23 (unused import), :140 (kT_eV); etna/core/device2d.py:318 (kT_eV); etna/core/cv_analysis.py:~66-84 (depletion_width_from_capacitance); etna/core/transient.py:232 (skip_init)
labels: [audit-v7, slice-A, dead-code, cleanup]
source: raw-A-core-solvers.md F11
related: A-05 (skip_init), A-10 (M_c), A-02 (effective_dos would gain a consumer)
---

## Summary

Several symbols are referenced only by their own tests, or not at all. Several locals are
computed and then discarded. Scope: "production" here means etna/ and app/. Test
references are counted separately.

## Reproduction / Evidence

Command I ran:
`grep -rnw --include='*.py' --include='*.ipynb' <sym> etna app tests scripts notebooks`,
excluding etna/core/sic_material.py (the defining file).

| symbol                           | production refs                                             | test refs                   |
| -------------------------------- | ----------------------------------------------------------- | --------------------------- |
| mobility_caughey_thomas          | 1 (device.py:22, unused import only)                        | 14 (tests/test_material.py) |
| effective_dos                    | 0                                                           | 6 (tests/test_material.py)  |
| depletion_width_from_capacitance | 0 (only its own def in cv_analysis.py)                      | 4 (tests/test_cv.py)        |
| M_c                              | 0                                                           | 1 (tests/test_material.py)  |
| E_D_hex                          | 0                                                           | 0                           |
| E_D_cub                          | 0                                                           | 0                           |
| NC_300                           | 0                                                           | 5 (tests/test_material.py)  |
| NV_300                           | 0                                                           | 1 (tests/test_material.py)  |
| `.B` field                       | 0 (`grep -rn '\.B\b' etna app --include='*.py'` -> no hits) | -                           |

Other checks:

- `grep -n kT_eV etna/core/device.py etna/core/device2d.py` finds only the assignments
  (device.py:140, device2d.py:318) and no reads.
- For `skip_init`, see A-05.

## Impact

This is maintenance noise. `NC_300`/`NV_300` are literature constants that can drift
from what `compute_ni(300)` actually computes without anyone noticing.

## Proposed fix

- Delete `mobility_caughey_thomas`, and port its tests to
  `mobility_caughey_thomas_T(N, 300, ...)`.
- Drop the unused import at device.py:22 and the `kT_eV` locals.
- Remove `E_D_hex`, `E_D_cub`, `B` and `M_c`, or keep the ones intended as documentation,
  marked as informational.
- Keep `NC_300`/`NV_300` only if a test asserts `NC_300 ≈ compute_ni(300)[1]` and
  `NV_300 ≈ compute_ni(300)[2]` within 5 %.
- Keep `effective_dos` if A-02 is adopted. Otherwise delete it.
- Remove `depletion_width_from_capacitance`, or use it in `plot_cv_evolution`.
- Resolve `skip_init` through A-05.

## Test plan

- `ruff check --select F401,F841 etna/core` is clean.
- The full test suite passes after the tests are ported.
- Add the new consistency test for NC_300/NV_300 if those fields are kept.

## Effort

S
