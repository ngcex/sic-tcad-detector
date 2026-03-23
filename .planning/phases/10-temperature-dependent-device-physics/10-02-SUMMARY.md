---
phase: 10-temperature-dependent-device-physics
plan: 02
subsystem: device-physics
tags:
  [
    temperature,
    4h-sic,
    intrinsic-concentration,
    mobility,
    srh-lifetime,
    devsim,
    regression,
  ]

# Dependency graph
requires:
  - phase: 10-01
    provides: "T-dependent material functions (intrinsic_concentration, mobility_caughey_thomas_T, srh_lifetime)"
provides:
  - "T-dependent device parameter setting in create_sic_device"
  - "T-dependent V_bi calculation in poisson.py"
  - "T-dependent Hecht CCE defaults in charge_collection.py"
  - "n_i and E_g keys in device_info dict for downstream use"
affects: [10-03, 11-dark-current, 12-transient]

# Tech tracking
tech-stack:
  added: []
  patterns:
    [
      "sentinel-based default detection (_UNSET) for backward-compatible T-dependent defaults",
    ]

key-files:
  created: []
  modified:
    - src/device.py
    - src/poisson.py
    - src/charge_collection.py
    - tests/test_material.py
    - tests/test_poisson.py
    - tests/test_charge_collection.py

key-decisions:
  - "Used _UNSET sentinel (not None) in hecht_cce to distinguish 'caller passed no value' from 'caller passed None', preserving full backward compatibility"
  - "For hecht_cce T-dependent defaults, used mu_max*(T/300)^gamma (low-doping limit) rather than Caughey-Thomas at some doping, since Hecht assumes bulk mobility"

patterns-established:
  - "device_info dict includes n_i and E_g keys for T-aware downstream consumers"
  - "Functions accepting T parameter default to T=300 for backward compatibility"

requirements-completed: [TEMP-06, TEMP-07]

# Metrics
duration: 5min
completed: 2026-03-23
---

# Phase 10 Plan 02: T-Dependent Pipeline Integration Summary

**Temperature threaded through device.py/poisson.py/charge_collection.py using Plan 01's material functions, with 8 regression tests proving bit-for-bit v1.0 parity at T=300K**

## Performance

- **Duration:** 5 min
- **Started:** 2026-03-23T15:58:09Z
- **Completed:** 2026-03-23T16:03:16Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments

- Replaced all hardcoded 300K values in device.py with T-dependent function calls (intrinsic_concentration, mobility_caughey_thomas_T, srh_lifetime)
- poisson.py extract_depletion_width now uses T-dependent n_i for V_bi calculation
- charge_collection.py hecht_cce and hecht_cce_partial_depletion accept T parameter with sentinel-based backward compatibility
- 8 new regression/physics tests confirm v1.0 parity at T=300K and correct T-dependence at 280K/350K
- All 186 tests pass (178 existing + 8 new), zero regression

## Task Commits

Each task was committed atomically:

1. **Task 1: Wire T-dependent functions into device.py, poisson.py, charge_collection.py** - `4e2c664` (feat)
2. **Task 2: Add regression tests confirming v1.0 parity at T=300K** - `32d9e4a` (test)

## Files Created/Modified

- `src/device.py` - Uses intrinsic_concentration(T), mobility_caughey_thomas_T(N,T), srh_lifetime(T); returns n_i and E_g in device_info
- `src/poisson.py` - extract_depletion_width uses T-dependent n_i for V_bi
- `src/charge_collection.py` - hecht_cce and hecht_cce_partial_depletion accept T/params with \_UNSET sentinel defaults
- `tests/test_material.py` - TestRegressionT300K (2 tests), TestTemperaturePhysics (3 tests)
- `tests/test_poisson.py` - TestVbiRegression (1 test)
- `tests/test_charge_collection.py` - TestHechtCCE300KRegression (2 tests)

## Decisions Made

- Used `_UNSET = object()` sentinel pattern in charge_collection.py instead of None defaults, because None could be a legitimate caller value and we need to distinguish "not provided" from "explicitly None"
- For hecht_cce T-dependent mobility defaults, used the low-doping limit `mu_max*(T/300)^gamma` rather than calling mobility_caughey_thomas_T with some doping level, since Hecht equation assumes bulk/intrinsic mobility

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All simulation modules now support T parameter (280-350K range verified)
- device_info dict includes n_i and E_g for downstream T-aware analysis
- Plan 03 (temperature sweep notebook and analysis) can proceed immediately

---

_Phase: 10-temperature-dependent-device-physics_
_Completed: 2026-03-23_

## Self-Check: PASSED

- [x] src/device.py exists
- [x] src/poisson.py exists
- [x] src/charge_collection.py exists
- [x] tests/test_material.py exists
- [x] tests/test_poisson.py exists
- [x] tests/test_charge_collection.py exists
- [x] 10-02-SUMMARY.md exists
- [x] Commit 4e2c664 (Task 1) found
- [x] Commit 32d9e4a (Task 2) found
