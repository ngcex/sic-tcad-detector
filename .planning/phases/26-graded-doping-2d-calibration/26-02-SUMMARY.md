---
phase: 26-graded-doping-2d-calibration
plan: 02
subsystem: infra
tags:
  [devsim, tcad, 2d, depletion-width, global-state, cylindrical, reset, poisson]

requires:
  - phase: 26-01
    provides: xfail test skeletons (TestResetStateLeak, TestGradedDopingSmoothness, TestCalibrationCV) + H1 diagnosis
provides:
  - "src/devsim_reset.py::reset_devsim_fully(preserve_solver=True) — clears all 7 cylindrical-axis globals, restores Cartesian assembly-model defaults, preserves solver settings"
  - "src/poisson.py::extract_depletion_width_2d_center(device_info, center_x_tol_cm=1e-6) — center-column W for 2D devices (y as depth)"
  - "src/optimization.py refactored to import reset_devsim_fully (no inline devsim.reset_devsim())"
  - "TestGradedDopingSmoothness + TestResetStateLeak now pass with real assertions (xfail removed)"
affects: [26-03, 26-04, graded-doping-calibration, phase-31-anisotropic]

tech-stack:
  added: []
  patterns:
    [
      "single-responsibility devsim_reset module replaces scattered inline resets",
      "2D depletion-width extraction mirrors 1D logic along y with center-column (|x|<tol) filter",
    ]

key-files:
  created:
    - src/devsim_reset.py
    - tests/test_devsim_reset.py
  modified:
    - src/optimization.py
    - src/poisson.py
    - tests/test_device2d.py

key-decisions:
  - "reset_devsim_fully restores Cartesian DEFAULTS for the 5 assembly-model globals (not empty string) and clears raxis_zero/raxis_variable to empty — empty assembly-model names would break the next solve."
  - "extract_depletion_width_2d_center is appended to poisson.py (not a new module) to sit beside the 1D twin extract_depletion_width_numerical, reusing module-level numpy/devsim imports."
  - "Reset-canary scalar is total cathode contact current (electrons+holes) via drift_diffusion.extract_contact_current — a reproducible value that shifts if cylindrical assembly weights leak."

patterns-established:
  - "Cross-module devsim reset: import reset_devsim_fully from src.devsim_reset instead of inlining devsim.reset_devsim()."
  - "2D field extraction: filter center column |x| < 1e-6 (x=0 symmetry plane), sort by y, apply depletion-edge logic along y."

requirements-completed: [CONS-01]

duration: ~95min
completed: 2026-06-15
---

# Phase 26 Plan 02: Reset hygiene + 2D depletion-width extractor Summary

**`reset_devsim_fully()` closes the seven-global cylindrical-axis leak (P03/P30) and `extract_depletion_width_2d_center()` measures center-column W along the y depth axis for 2D devices — the two infrastructure pieces Plan 03's calibration loop and Plan 04's regression canary depend on.**

## Performance

- **Duration:** ~95 min
- **Started:** 2026-06-15T12:?? (worktree base 387a995)
- **Completed:** 2026-06-15T13:38Z
- **Tasks:** 3
- **Files created:** 2 · **Files modified:** 3

## Accomplishments

- New `src/devsim_reset.py` with `reset_devsim_fully(preserve_solver=True)`:
  - Module constant `_CYLINDRICAL_GLOBALS` = the seven leaking globals (verbatim):
    `"raxis_zero"`, `"raxis_variable"`, `"node_volume_model"`, `"edge_couple_model"`,
    `"element_edge_couple_model"`, `"element_node0_volume_model"`,
    `"element_node1_volume_model"`.
  - Five-step workflow: snapshot solver -> enumerate-and-delete devices (P20) ->
    restore Cartesian defaults / clear `raxis_*` -> `devsim.reset_devsim()` ->
    restore solver. Idempotent on empty state.
- `src/optimization.py` refactored: imports `reset_devsim_fully`, removed the inline
  solver save/restore block and the `devsim.reset_devsim()` call. The Mj-3 sweep
  validity gate (`full_depletion_voltage_graded` + validity columns) was preserved untouched.
- New `src/poisson.py::extract_depletion_width_2d_center(device_info, center_x_tol_cm=1e-6)`:
  coordinate convention **x = lateral, y = depth, x ≈ 0 is the center column**. Filters
  `|x| < center_x_tol_cm`, sorts by y, finds the junction-side depletion edge where
  electrons recover to 50% of local Donors, returns W (cm) clamped to `[0, epi_thickness]`.
  Raises `ValueError` for non-2D devices or an empty center column.
- Previously-xfail tests now pass with real assertions: `TestGradedDopingSmoothness`
  (P27 no->50%-jump in center-column Donors) and `TestResetStateLeak` (cylindrical-leak
  canary, CONS-01 success criterion #3). New `TestExtractDepletionWidth2DCenter` class
  (import / 1D-rejection / equilibrium band / reverse-bias expansion) and
  `tests/test_devsim_reset.py` (6 tests).

## Task Commits

1. **Task 1 (RED): failing tests for reset_devsim_fully** - `828fca0` (test)
2. **Task 1 (GREEN): reset_devsim_fully + optimization refactor** - `71d5512` (feat)
3. **Task 2 (RED): 2D W extractor + smoothness tests** - `8b8ec7d` (test)
4. **Task 2 (GREEN): extract_depletion_width_2d_center** - `2bd4177` (feat)
5. **Task 3: TestResetStateLeak cylindrical-leak canary** - `318afca` (test)

## Files Created/Modified

- `src/devsim_reset.py` (created) - `reset_devsim_fully` + `_CYLINDRICAL_GLOBALS` + `_CARTESIAN_DEFAULTS`
- `tests/test_devsim_reset.py` (created) - 6 tests for the reset utility
- `src/optimization.py` (modified) - imports reset utility; inline reset removed; Mj-3 gate preserved
- `src/poisson.py` (modified) - appended `extract_depletion_width_2d_center` after the 1D extractor
- `tests/test_device2d.py` (modified) - wired `TestGradedDopingSmoothness`, `TestResetStateLeak`; added `TestExtractDepletionWidth2DCenter`

## Test Status (Phase 26 surface)

| Class / test                                                         | Status after Plan 02     | Owner   |
| -------------------------------------------------------------------- | ------------------------ | ------- |
| `TestGradedDopingSmoothness::test_graded_doping_smoothness_no_kinks` | PASS (real assertion)    | Plan 02 |
| `TestResetStateLeak::test_reset_after_alt_structures`                | PASS (slow, real canary) | Plan 02 |
| `TestExtractDepletionWidth2DCenter` (4 tests)                        | PASS                     | Plan 02 |
| `tests/test_devsim_reset.py` (6 tests)                               | PASS                     | Plan 02 |
| `TestReverseBiasConvergence` (6 tests)                               | still xfail              | Plan 03 |
| `TestCalibrationCV::test_2d_vs_1d_cv_centerline`                     | still xfail              | Plan 03 |

Plan 03 consumes both new helpers: `reset_devsim_fully` between calibration trials and
`extract_depletion_width_2d_center` to compute the cost function (center-column W vs the
1D-twin C-V targets). Plan 04's cylindrical-leak regression also reuses `reset_devsim_fully`.

## Decisions Made

- Restore Cartesian default assembly-model names (not empty strings) — empty model names
  would leave devsim unable to assemble the next device. Only `raxis_zero`/`raxis_variable`
  are cleared to empty.
- Append the 2D extractor to `poisson.py` rather than create a new module, so it sits beside
  the 1D twin and shares the existing numpy/devsim imports (no new imports, per plan).

## Deviations from Plan

None - plan executed exactly as written. (No Rule 1-4 deviations; the Mj-3 sweep gate flagged
in the project-specific note was preserved as instructed.)

## Issues Encountered

- **devsim multi-DD-build process abort (environmental, not a code defect):** Running the
  _entire_ `tests/test_device2d.py` file (all Poisson-only classes + every slow DD-building
  class + the reset canary) in a single interpreter triggers `Fatal Python error: Aborted`
  inside devsim after ~26 tests pass — a devsim process-global/resource exhaustion when many
  full DD device builds stack in one process. **Every plan verification command passes
  individually:** `tests/test_optimization.py` (11), `TestGradedDopingSmoothness` (1),
  `TestResetStateLeak -m slow` (1), `TestDevice2DCreation`+`TestDoping2D` (7),
  `TestExtractDepletionWidth2DCenter`+`TestGradedDopingSmoothness -m slow` (2),
  `tests/test_devsim_reset.py` (6). This matches the project's documented quick-vs-slow,
  class-by-class run guidance (RESEARCH.md sampling section). No logic in this plan is the
  cause; the original Poisson-only classes also abort the same way when stacked with the new
  DD classes. Recommendation for Plan 03/04 and CI: run devsim DD classes class-by-class (or
  with `pytest-forked`/`--forked`), not as one monolithic file invocation.

## Known Stubs

None — both helpers are fully wired to real devsim node data; no placeholder/empty-data paths.

## Next Phase Readiness

- `reset_devsim_fully` and `extract_depletion_width_2d_center` are ready for Plan 03's
  `calibrate_graded_doping_2d` loop and the cost function.
- `TestReverseBiasConvergence` and `TestCalibrationCV` remain xfail for Plan 03 to wire after
  calibration produces the new 2D defaults.

## Self-Check: PASSED

- Files verified present: `src/devsim_reset.py`, `tests/test_devsim_reset.py`, `26-02-SUMMARY.md`
- Commits verified present: `828fca0`, `71d5512`, `8b8ec7d`, `2bd4177`, `318afca`

---

_Phase: 26-graded-doping-2d-calibration_
_Completed: 2026-06-15_
