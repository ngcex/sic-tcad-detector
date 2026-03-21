---
phase: 05-parametric-studies-and-publication
plan: 01
subsystem: simulation
tags: [parametric-sweep, cce, doping, itertools, json-caching]

# Dependency graph
requires:
  - phase: 04-flash-plasma-recombination
    provides: "cce_vs_dose_rate function, Auger recombination model"
provides:
  - "parametric_cce_sweep function for multi-dimensional CCE sweeps"
  - "Doping-parametrized cce_vs_dose_rate (N_D_junction, N_D_bulk, L_transition kwargs)"
  - "JSON save/load helpers for parametric results caching"
affects: [05-02-publication-figures]

# Tech tracking
tech-stack:
  added: [itertools, json, ast]
  patterns:
    [tuple-keyed-dict, proportional-doping-scaling, graceful-degradation-sweep]

key-files:
  created: []
  modified:
    - src/flash_recombination.py
    - tests/test_flash_recombination.py

key-decisions:
  - "N_D_junction scaled proportionally with N_D_bulk to preserve graded profile shape"
  - "ast.literal_eval for safe tuple-key reconstruction from JSON strings"
  - "Failed sweep combinations stored as None with warning logged (graceful degradation)"

patterns-established:
  - "Proportional doping scaling: N_D_junction = base * (N_D_bulk / base_bulk)"
  - "Tuple-keyed results dict with string serialization for JSON persistence"

requirements-completed: [FLASH-04]

# Metrics
duration: 3min
completed: 2026-03-21
---

# Phase 5 Plan 1: Parametric CCE Sweep Infrastructure Summary

**Parametric cce_vs_dose_rate with doping kwargs and multi-dimensional sweep wrapper using itertools.product over epi/doping/bias combinations with JSON caching**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-21T21:00:05Z
- **Completed:** 2026-03-21T21:03:00Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Extended cce_vs_dose_rate with N_D_junction, N_D_bulk, L_transition keyword arguments (backward-compatible defaults)
- Created parametric_cce_sweep that iterates over all epi/doping/bias combinations with proportional junction doping scaling
- Added save_parametric_results and load_parametric_results for JSON serialization with tuple key handling
- All 8 tests pass including 3 new tests (doping params, save/load roundtrip, sweep structure)

## Task Commits

Each task was committed atomically:

1. **Task 1: Parametrize cce_vs_dose_rate and create parametric_cce_sweep** - `f16cb45` (feat)
2. **Task 2: Add tests for parametric sweep and doping parametrization** - `5756151` (test)

## Files Created/Modified

- `src/flash_recombination.py` - Added doping kwargs to cce_vs_dose_rate, parametric_cce_sweep, save/load helpers
- `tests/test_flash_recombination.py` - Added 3 new test classes for doping params, JSON roundtrip, sweep structure

## Decisions Made

- N_D_junction scaled proportionally with N_D_bulk (ratio preserved from base values) to maintain graded profile shape across sweep
- ast.literal_eval used for safe tuple-key reconstruction from JSON (no arbitrary code execution)
- Failed sweep combinations stored as None with exc_info logging for debugging

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- parametric_cce_sweep ready for 05-02 publication figure generation
- JSON caching enables re-running analysis without re-solving (hours of compute time saved)

---

_Phase: 05-parametric-studies-and-publication_
_Completed: 2026-03-21_
