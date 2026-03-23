---
phase: 10-temperature-dependent-device-physics
plan: 01
subsystem: material-physics
tags:
  [
    varshni,
    caughey-thomas,
    srh,
    temperature,
    4h-sic,
    bandgap,
    mobility,
    intrinsic-concentration,
  ]

# Dependency graph
requires: []
provides:
  - "T-dependent bandgap via Varshni equation"
  - "Calibrated T-dependent intrinsic carrier concentration anchored at n_i_300=5e-9"
  - "T-dependent Caughey-Thomas mobility with phonon scattering exponents"
  - "T-dependent effective density of states with T^(3/2) scaling"
  - "T-dependent SRH lifetime with power-law scaling"
affects: [10-02, 10-03, 11-dark-current, 12-transient]

# Tech tracking
tech-stack:
  added: []
  patterns:
    [
      "calibration-factor pattern for anchoring first-principles to literature values",
    ]

key-files:
  created: []
  modified:
    - src/sic_material.py
    - tests/test_material.py

key-decisions:
  - "Calibrated E_g_0 to 3.2965625 eV so Varshni gives exactly 3.26 eV at 300K (original 3.265 gave 3.228)"
  - "Used calibration factor n_i(T) = n_i_300 * compute_ni(T)/compute_ni(300) to anchor n_i at 5e-9"
  - "T-scaling applied only to mu_max in Caughey-Thomas; mu_min kept constant"

patterns-established:
  - "T-dependent functions accept optional params=None defaulting to SiC4H_Parameters()"
  - "All T-dependent functions return exact 300K values when T=300 for regression safety"

requirements-completed: [TEMP-01, TEMP-02, TEMP-03, TEMP-04, TEMP-05]

# Metrics
duration: 3min
completed: 2026-03-23
---

# Phase 10 Plan 01: T-Dependent Material Properties Summary

**Five T-dependent material functions (bandgap, n_i, mobility, DOS, SRH lifetime) with Varshni/Caughey-Thomas/power-law models anchored to exact 300K values**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-23T15:52:27Z
- **Completed:** 2026-03-23T15:55:25Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Added 5 new T-dependent public functions to sic_material.py (bandgap, intrinsic_concentration, mobility_caughey_thomas_T, effective_dos, srh_lifetime)
- Added 3 new dataclass fields (gamma_n, gamma_p, alpha_tau) for T-scaling exponents from Ayalew thesis
- Calibrated E_g_0 so Varshni equation reproduces exactly 3.26 eV at 300K
- 24 new unit tests validating T-dependence at 280K, 300K, 350K against literature values
- All 45 tests pass (21 existing + 24 new), zero regression

## Task Commits

Each task was committed atomically:

1. **Task 1: Add T-dependent material property functions** - `5881766` (feat)
2. **Task 2: Add unit tests for T-dependent functions** - `d55bd69` (test)

## Files Created/Modified

- `src/sic_material.py` - Added 5 T-dependent functions + 3 dataclass fields, calibrated E_g_0
- `tests/test_material.py` - Added 5 test classes (24 tests) for T-dependent validation

## Decisions Made

- Calibrated E_g_0 from 3.265 to 3.2965625 eV: the original Ioffe value produced 3.228 eV at 300K via Varshni, not 3.26. Recalculated E_g_0 = E_g(300) + alpha\*300^2/(300+beta) to get exact match. This only affects the new bandgap() function; compute_ni() still uses its hardcoded 3.265.
- Kept compute_ni() completely unchanged for backward compatibility -- new intrinsic_concentration() wraps it with the calibration factor approach.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] E_g_0 Varshni calibration mismatch**

- **Found during:** Task 1 (bandgap function implementation)
- **Issue:** E_g_0=3.265 with alpha=6.5e-4, beta=1300 produces 3.228 eV at 300K, not 3.26 as required
- **Fix:** Recalculated E_g_0 = 3.2965625 eV to produce exactly 3.26 eV at 300K
- **Files modified:** src/sic_material.py (E_g_0 field in SiC4H_Parameters)
- **Verification:** bandgap(300) returns exactly 3.26; all existing tests still pass (compute_ni uses hardcoded value)
- **Committed in:** 5881766 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Essential for meeting the exact-300K regression requirement. No scope creep.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All 5 T-dependent material functions ready for Plan 02 to thread temperature through the simulation pipeline
- Functions follow consistent API (T, params=None) for easy integration
- Exact 300K regression guaranteed by calibration factor approach

---

_Phase: 10-temperature-dependent-device-physics_
_Completed: 2026-03-23_

## Self-Check: PASSED

- [x] src/sic_material.py exists
- [x] tests/test_material.py exists
- [x] 10-01-SUMMARY.md exists
- [x] Commit 5881766 (Task 1) found
- [x] Commit d55bd69 (Task 2) found
