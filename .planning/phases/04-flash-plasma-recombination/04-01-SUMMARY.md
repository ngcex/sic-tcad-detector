---
phase: 04-flash-plasma-recombination
plan: 01
subsystem: simulation
tags: [auger, recombination, flash, devsim, continuation-solver, high-injection]

requires:
  - phase: 03-charge-collection
    provides: "DD solver with SRH recombination, CCE computation, generation injection"
provides:
  - "Auger recombination model (UAuger) with full Jacobian derivatives"
  - "Generation-rate continuation solver for high-injection convergence"
  - "Combined SRH+Auger recombination in DD continuity equations"
affects: [04-02, dose-rate-dependent-cce, flash-plasma-recombination]

tech-stack:
  added: []
  patterns:
    [
      "generation-rate continuation ramping with bisection fallback",
      "Auger node model with analytical Jacobian derivatives",
    ]

key-files:
  created:
    - src/flash_recombination.py
    - tests/test_flash_recombination.py
  modified: []

key-decisions:
  - "Auger added after bias ramp but before generation for Jacobian stability"
  - "Continuation solver uses 5-step linear ramp with up to 3 bisection retries per step"
  - "RadGenRate existence checked at Auger setup time for flexible call ordering"

patterns-established:
  - "Continuation solver pattern: ramp physics parameter in small increments for Newton convergence at extreme conditions"
  - "Auger model pattern: devsim expression + CreateNodeModelDerivative for all solution variables"

requirements-completed: [FLASH-01, FLASH-02]

duration: 2min
completed: 2026-03-21
---

# Phase 4 Plan 1: Auger Recombination Model Summary

**Auger recombination (R_Auger = (C_n*n + C_p*p)*(n*p - n_i^2)) with generation-rate continuation solver for FLASH dose rates up to 230 Gy/s**

## Performance

- **Duration:** 2 min
- **Started:** 2026-03-21T20:15:37Z
- **Completed:** 2026-03-21T20:18:06Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Auger recombination model with correct expression and all Jacobian derivatives for Newton solver convergence
- Combined SRH+Auger recombination active in electron and hole continuity equations
- Generation-rate continuation solver with bisection fallback handles 230 Gy/s without divergence
- 5 integration tests verify model correctness: creation, equilibrium, high-injection CCE reduction, convergence, coefficients

## Task Commits

Each task was committed atomically:

1. **Task 1: Create flash_recombination module with Auger model and continuation solver** - `acecf24` (feat)
2. **Task 2: Integration tests for Auger model and high-injection convergence** - `ff6ce9e` (test)

## Files Created/Modified

- `src/flash_recombination.py` - Auger recombination model setup and generation-rate continuation solver
- `tests/test_flash_recombination.py` - 5 integration tests for Auger model and convergence

## Decisions Made

- Auger model is added after bias ramping but before generation injection, ensuring Jacobian includes Auger terms from the start
- Continuation solver defaults to 5 steps (20% increments) with up to 3 bisection retries per failed step
- RadGenRate model existence is checked at Auger setup time; zero placeholder created if not yet present, allowing flexible call ordering (Auger before or after generation)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Auger model and continuation solver ready for 04-02 (CCE vs dose-rate sweep)
- Both add_auger_recombination and solve_with_continuation are importable and tested
- Combined SRH+Auger generation terms properly update DD continuity equations

---

_Phase: 04-flash-plasma-recombination_
_Completed: 2026-03-21_
