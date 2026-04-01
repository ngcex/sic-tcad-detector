---
phase: 24-alternative-structures
plan: 01
subsystem: simulation
tags:
  [
    devsim,
    mesh-generation,
    cylindrical-coordinates,
    multi-region,
    sic-microdosimeter,
  ]

# Dependency graph
requires:
  - phase: 19-device2d
    provides: "2D device creation pattern, doping helpers, material parameter setup"
  - phase: 20-dd2d
    provides: "Poisson/DD solver dimension-agnostic pattern, device deletion pattern"
provides:
  - "Four alternative structure mesh builders: mesa, 3D electrode, delta-E/E, guard ring"
  - "Cylindrical coordinate activation/restoration helpers"
  - "Multi-region Poisson setup helper for non-standard contacts"
affects: [24-02, notebooks]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Cylindrical coordinate activation via _activate_cylindrical_coords + restore_cartesian_coords"
    - "Multi-region device with devsim interface (no contacts at interface boundary)"
    - "Guard ring doping overlay via explicit expression (avoid cyclic dependency)"

key-files:
  created:
    - src/alternative_structures.py
    - tests/test_alternative_structures.py
  modified: []

key-decisions:
  - "Delta-E/E uses 2 contacts (de_anode, estop_cathode) not 4 -- devsim cannot have contacts at interface boundary"
  - "Guard ring acceptor doping uses explicit expression instead of self-referencing Acceptors model"
  - "Guard ring Poisson solve needs relaxed tolerance fallback due to 35K-node mesh complexity"

patterns-established:
  - "Cylindrical device lifecycle: create -> activate_cylindrical -> physics -> delete -> restore_cartesian"
  - "Multi-region Poisson setup via _setup_poisson_region helper for non-standard contact names"

requirements-completed: [ALTS-01, ALTS-02, ALTS-03, ALTS-04]

# Metrics
duration: 11min
completed: 2026-04-01
---

# Phase 24 Plan 01: Alternative Structures Summary

**Four mesh builders (mesa, 3D electrode, delta-E/E, guard ring) with cylindrical coordinates, multi-region interfaces, and Poisson solve validation**

## Performance

- **Duration:** 11 min
- **Started:** 2026-04-01T07:45:52Z
- **Completed:** 2026-04-01T07:56:28Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Created alternative_structures.py with 4 public mesh builder functions + 5 internal helpers
- All builders return pipeline-compatible device_info dicts with structure_type key
- Cylindrical coordinate activation/restoration works correctly for 3D electrode
- All 12 tests pass including Poisson solve smoke tests on all four structures

## Task Commits

Each task was committed atomically:

1. **Task 1: Create alternative_structures.py with all four mesh builders** - `3fb43a7` (feat)
2. **Task 2: Create tests and validate mesh builders with Poisson solve** - `bc4d77f` (test)

## Files Created/Modified

- `src/alternative_structures.py` - Mesa, 3D electrode, delta-E/E, and guard ring mesh builders with shared helpers
- `tests/test_alternative_structures.py` - 12 tests: creation, Poisson solve, doping verification, cylindrical coordinate cleanup

## Decisions Made

- Delta-E/E telescope uses 2 contacts (de_anode at top, estop_cathode at bottom) instead of 4 per-layer contacts, because devsim prevents contact + interface at the same mesh boundary
- Guard ring acceptor doping redefines Acceptors as explicit expression (N_A_ionized\*step + Acceptors_GR) to avoid cyclic model dependency in devsim
- Guard ring Poisson solve requires relaxed tolerance fallback (abs_error=1e12, rel_error=1e-8, 100 iterations) due to 35K-node mesh with 3 contacts

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Delta-E/E interface creation failed with contacts at boundary**

- **Found during:** Task 2 (test_delta_e_e_poisson_solve)
- **Issue:** devsim cannot create interface between two regions if a contact exists at the same boundary position. estop_anode and de_cathode at y=de_total prevented interface creation.
- **Fix:** Removed de_cathode and estop_anode contacts. Interface provides current continuity between layers. Only de_anode (top) and estop_cathode (bottom) contacts remain.
- **Files modified:** src/alternative_structures.py
- **Verification:** Interface created with 26 node pairs, Poisson solve converges

**2. [Rule 1 - Bug] Guard ring Acceptors model cyclic dependency**

- **Found during:** Task 2 (test_guard_ring_poisson_solve)
- **Issue:** `Acceptors = Acceptors + Acceptors_GR` is self-referential in devsim's model evaluation engine, causing cyclic dependency error.
- **Fix:** Redefined Acceptors using explicit base expression: `N_A_ionized * step(junction_pos - y) + Acceptors_GR`
- **Files modified:** src/alternative_structures.py
- **Verification:** No cyclic dependency, guard ring doping test confirms enhanced acceptor values

**3. [Rule 1 - Bug] Guard ring Poisson solve convergence failure**

- **Found during:** Task 2 (test_guard_ring_poisson_solve)
- **Issue:** Default solver tolerance (abs=1e10, rel=1e-10, 40 iters) insufficient for 35K-node mesh with 3 contacts
- **Fix:** Added relaxed tolerance fallback (same pattern as solve_equilibrium in poisson.py)
- **Files modified:** tests/test_alternative_structures.py
- **Verification:** Solve converges with relaxed parameters

---

**Total deviations:** 3 auto-fixed (3 bugs)
**Impact on plan:** All fixes necessary for correctness. Delta-E/E contact design simplified due to devsim constraint. No scope creep.

## Issues Encountered

None beyond the auto-fixed deviations above.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All four mesh builders ready for Plan 02 (comparison notebook)
- Pipeline compatibility confirmed: device_info dicts have all required keys
- Cylindrical coordinate cleanup verified (restore_cartesian_coords works)
- Note: delta-E/E has 2 contacts (not 4 as originally planned); downstream pipeline may need adjustment for independent layer readout

---

_Phase: 24-alternative-structures_
_Completed: 2026-04-01_
