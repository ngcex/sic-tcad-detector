---
phase: 03-charge-collection-efficiency
plan: 01
subsystem: physics
tags:
  [
    hecht-equation,
    alpha-particle,
    proton-bragg-peak,
    cce,
    generation-rate,
    numpy,
    scipy,
  ]

# Dependency graph
requires:
  - phase: 02-electrical-characterization
    provides: "SiC4H_Parameters (mu, tau), DD solver infrastructure"
provides:
  - "alpha_generation_profile: Am-241 spatial generation profile with erfc roll-off"
  - "proton_generation_profile: 30-150 MeV flat entrance-dose profiles"
  - "dose_rate_to_generation: Gy/s to carrier pair rate conversion"
  - "hecht_cce: two-carrier analytical CCE benchmark"
  - "compute_cce_from_current: CCE extraction from DD simulation"
  - "hecht_cce_partial_depletion: extended Hecht with diffusion collection"
affects: [03-02-PLAN, 03-03-PLAN, 04-flash-dynamics]

# Tech tracking
tech-stack:
  added: []
  patterns: [erfc-smoothed-generation-profile, log-log-range-interpolation]

key-files:
  created:
    - src/generation_profiles.py
    - src/charge_collection.py
    - tests/test_generation_profiles.py
    - tests/test_charge_collection.py
  modified: []

key-decisions:
  - "Alpha profile uses erfc roll-off at 0.8*range with sigma=0.1*range for smooth DD solver compatibility"
  - "Proton profiles flat within detector for all therapeutic energies (range >> 10um)"
  - "Partial depletion Hecht uses average diffusion collection probability in neutral region"

patterns-established:
  - "Generation profiles return spatial rates (cm^-1 or cm^-3 s^-1) with CGS units"
  - "CCE functions accept scalar or array voltage, return clipped [0,1] values"

requirements-completed: [CCE-04, CCE-02]

# Metrics
duration: 3min
completed: 2026-03-21
---

# Phase 3 Plan 01: Generation Profiles & Hecht Equation Summary

**Alpha/proton generation profiles and two-carrier Hecht equation with partial-depletion extension for CCE validation**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-21T18:13:58Z
- **Completed:** 2026-03-21T18:17:34Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments

- Generation profile module with alpha (Am-241 erfc-smoothed) and proton (NIST PSTAR range-scaled) models
- Hecht equation producing CCE = 1.0 at V=40V for 10um SiC, confirming literature consistency
- Partial-depletion extension combining drift and diffusion collection for comparison with DD solver
- Full test coverage: 35 tests covering normalization, smoothness, physics limits, vectorization

## Task Commits

Each task was committed atomically:

1. **Task 1: Generation profile models (alpha + proton Bragg peak)** - `6142ef6` (feat)
2. **Task 2: Hecht equation and CCE utilities** - `2289007` (feat)

## Files Created/Modified

- `src/generation_profiles.py` - Alpha and proton generation rate profiles with dose rate conversion
- `src/charge_collection.py` - Hecht equation, CCE from current ratio, partial depletion extension
- `tests/test_generation_profiles.py` - 16 tests for profile normalization, smoothness, range scaling
- `tests/test_charge_collection.py` - 19 tests for Hecht physics limits, vectorization, partial depletion

## Decisions Made

- Alpha profile uses erfc roll-off (not step function) at 0.8\*range to avoid numerical ringing in DD solver
- Proton profiles are flat within the thin detector for all therapeutic energies (30-150 MeV ranges far exceed 10um)
- Partial depletion model uses average exponential diffusion collection probability in neutral region
- Log-log interpolation for proton ranges between tabulated NIST energies

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed partial depletion test comparison logic**

- **Found during:** Task 2
- **Issue:** Test compared partial-depletion CCE against standard Hecht with d=d_epi, but SiC's long mu\*tau gives CCE~1.0 even at 5V with d=10um, making the comparison invalid
- **Fix:** Changed comparison to drift-fraction-only CCE (W/d_epi \* Hecht(V,W)), which correctly shows diffusion adds to collection
- **Files modified:** tests/test_charge_collection.py
- **Verification:** All 19 tests pass
- **Committed in:** 2289007

---

**Total deviations:** 1 auto-fixed (1 bug fix)
**Impact on plan:** Test logic correction, no scope change.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Generation profiles ready for integration into DD solver (Plan 02)
- Hecht equation provides analytical benchmark for CCE validation (Plan 03)
- No devsim dependency in either module (pure numpy/scipy)

---

_Phase: 03-charge-collection-efficiency_
_Completed: 2026-03-21_
