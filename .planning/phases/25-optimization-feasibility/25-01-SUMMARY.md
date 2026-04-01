---
phase: 25-optimization-feasibility
plan: 01
subsystem: optimization
tags: [parametric-sweep, noise-floor, structure-scoring, cce-uniformity, pandas]

# Dependency graph
requires:
  - phase: 20-charge-collection-2d
    provides: create_2d_dd_device, cce_lateral_scan for parametric sweep
  - phase: 14-dark-current
    provides: setup_tat_model, extract_contact_current for noise estimation
  - phase: 23-microdosimetry
    provides: mean_chord_length for y_min computation
provides:
  - microdosimetric_sweep: parametric grid search over SV geometry/doping/bias
  - estimate_noise_floor: shot-noise-limited detection thresholds from dark current
  - score_structures: multi-criteria weighted ranking of detector structures
  - get_dark_current_2d: TAT dark current extraction for 2D devices
affects: [25-02-feasibility-report, optimization-notebooks]

# Tech tracking
tech-stack:
  added: []
  patterns:
    [
      multi-criteria-scoring,
      parametric-sweep-with-cleanup,
      shot-noise-estimation,
    ]

key-files:
  created:
    - src/optimization.py
    - tests/test_optimization.py
  modified:
    - src/charge_collection_2d.py

key-decisions:
  - "3 public functions + 1 helper in single optimization.py module"
  - "Shot noise only for noise floor (readout electronics noise excluded, documented)"
  - "Min-max normalization with higher/lower-is-better classification for scoring"

patterns-established:
  - "Parametric sweep with try/finally devsim.delete_device cleanup per iteration"
  - "max_configs truncation with warning for large parameter grids"
  - "Pure computation tests only for devsim-dependent modules (integration tested in notebooks)"

requirements-completed: [FEAS-01, FEAS-02, FEAS-03]

# Metrics
duration: 13min
completed: 2026-04-01
---

# Phase 25 Plan 01: Optimization Module Summary

**Parametric sweep, shot-noise floor estimation, and multi-criteria structure scoring for SiC microdosimeter design space exploration**

## Performance

- **Duration:** 13 min
- **Started:** 2026-04-01T10:38:45Z
- **Completed:** 2026-04-01T10:52:18Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments

- Extended create_2d_dd_device with \*\*device_kwargs for parametric sweep parameterization
- Created optimization.py with 3 public functions and 1 helper integrating charge collection, dark current, and microdosimetry modules
- 9 unit tests covering noise floor estimation and structure scoring logic

## Task Commits

Each task was committed atomically:

1. **Task 1: Extend create_2d_dd_device to accept device parameter kwargs** - `7d3e990` (feat)
2. **Task 2: Create optimization.py with sweep, noise floor, and scoring functions** - `6e21f0b` (feat)

## Files Created/Modified

- `src/optimization.py` - Parametric sweep, noise floor estimation, structure scoring, dark current helper
- `src/charge_collection_2d.py` - Added \*\*device_kwargs forwarding to create_sic_2d_device
- `tests/test_optimization.py` - 9 unit tests for scoring normalization, ranking, and noise estimation

## Decisions Made

- Kept all optimization functions in a single module (not split by concern) for discoverability
- estimate_noise_floor is pure computation (no device creation) for fast unit testing
- Shot noise only for noise floor -- documented that real systems will have higher noise from readout electronics
- Min-max normalization with 0.5 fallback for equal-value metrics in score_structures

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed test_ranking_order expected winner**

- **Found during:** Task 2 (test creation)
- **Issue:** Test expected guard_ring as top-ranked structure, but with default weights (0.30 CCE, 0.20 noise, 0.20 spectral, 0.30 fab), 3d_electrode wins (best in 3 of 4 metrics)
- **Fix:** Corrected test to expect 3d_electrode as winner, added planar-last assertion
- **Files modified:** tests/test_optimization.py
- **Verification:** All 9 tests pass

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Test expectation corrected to match actual scoring math. No scope change.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- optimization.py ready for feasibility report notebook (Plan 25-02)
- All 4 exports available: microdosimetric_sweep, estimate_noise_floor, score_structures, get_dark_current_2d
- Full test suite passes (no regression from charge_collection_2d change)

---

_Phase: 25-optimization-feasibility_
_Completed: 2026-04-01_
