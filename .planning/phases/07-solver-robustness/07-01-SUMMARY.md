---
phase: 07-solver-robustness
plan: 01
subsystem: solver
tags:
  [devsim, transient, time_node_model, equation-registration, regression-test]

# Dependency graph
requires:
  - phase: 03-charge-collection
    provides: "add_generation_to_dd with devsim.equation re-registration"
  - phase: 04-flash-recombination
    provides: "add_auger_recombination with devsim.equation re-registration"
provides:
  - "Correct time_node_model preservation in all equation re-registrations"
  - "Regression test for transient solve capability after generation and Auger setup"
  - "Unambiguous ROADMAP SC-3 documentation of flat CCE null result"
affects: [future-transient-simulations]

# Tech tracking
tech-stack:
  added: []
  patterns: ["time_node_model preservation on devsim.equation re-registration"]

key-files:
  created: []
  modified:
    - src/charge_collection.py
    - src/flash_recombination.py
    - tests/test_flash_recombination.py
    - .planning/ROADMAP.md

key-decisions:
  - "transient_dc with tdelta=0 required before BDF1 step to initialize devsim time data storage"

patterns-established:
  - "Always include time_node_model when re-registering devsim continuity equations to preserve transient capability"

requirements-completed: []

# Metrics
duration: 3min
completed: 2026-03-21
---

# Phase 7 Plan 1: Solver Robustness Summary

**Fixed latent transient-solve bug by preserving time_node_model (NCharge/PCharge) in all devsim.equation() re-registrations, with BDF1 regression test**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-21T22:01:55Z
- **Completed:** 2026-03-21T22:04:39Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments

- Fixed 4 devsim.equation() re-registration calls to include time_node_model, preventing silent loss of charge storage terms during transient solves
- Added TestTransientCapabilityPreserved regression test that creates a full DD device, adds generation + Auger, and verifies a BDF1 transient step succeeds
- Updated ROADMAP Phase 4 SC-3 to unambiguously state flat CCE as accepted null result with physical justification

## Task Commits

Each task was committed atomically:

1. **Task 1: Add time_node_model to equation re-registrations and add regression test** - `e99ec8a` (fix)
2. **Task 2: Update ROADMAP Phase 4 SC-3 wording** - `daab947` (docs)

## Files Created/Modified

- `src/charge_collection.py` - Added time_node_model="NCharge" and "PCharge" to add_generation_to_dd re-registrations
- `src/flash_recombination.py` - Added time_node_model="NCharge" and "PCharge" to add_auger_recombination re-registrations
- `tests/test_flash_recombination.py` - Added TestTransientCapabilityPreserved with BDF1 transient regression test
- `.planning/ROADMAP.md` - Updated SC-3 wording for flat CCE null result

## Decisions Made

- Used transient_dc with tdelta=0 as initialization step before BDF1 transient solve, since devsim requires time data storage to be initialized before BDF1 can execute

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Added transient_dc initialization before BDF1 step**

- **Found during:** Task 1 (regression test)
- **Issue:** devsim requires transient_dc solve with tdelta=0 to initialize time data storage before a BDF1 transient step can run; without it, devsim crashes with "UNEXPECTED missing time data"
- **Fix:** Added transient_dc initialization step in the regression test before the BDF1 solve
- **Files modified:** tests/test_flash_recombination.py
- **Verification:** Test passes successfully with transient_dc init + BDF1 step
- **Committed in:** e99ec8a (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Auto-fix necessary for correct test implementation. No scope creep.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Transient solve capability is now preserved through generation and Auger setup
- All 147 tests pass with no regressions

---

_Phase: 07-solver-robustness_
_Completed: 2026-03-21_
