---
phase: 16-carrier-removal-cv-evolution
plan: 01
subsystem: physics-simulation
tags: [carrier-removal, cv-analysis, phi-crit, devsim, radiation-damage]

# Dependency graph
requires:
  - phase: 14-cce-vs-fluence
    provides: staged device creation pattern, compute_damaged_params, apply_damaged_params
  - phase: 08-cv-depletion
    provides: cv_sweep, junction_capacitance
provides:
  - compute_phi_crit() for critical fluence detection from graded doping profile
  - cv_at_fluence() for C-V curves at arbitrary proton fluence
  - plot_cv_evolution() for overlaying C-V curves across fluence levels
affects: [16-02, carrier-removal-cv-evolution]

# Tech tracking
tech-stack:
  added: []
  patterns: [fluence-as-temperature C-V sweep, Phi_crit threshold warning/abort]

key-files:
  created: []
  modified:
    - src/radiation_damage.py
    - src/cv_analysis.py
    - tests/test_radiation_damage.py
    - tests/test_cv.py

key-decisions:
  - "Phi_crit computed from min(N_D > 0) in profile, not from mean or bulk value"
  - "cv_at_fluence returns None (not raises) when fluence >= Phi_crit, matching prior fluence-sweep patterns"
  - "Fluence threshold warning at 90% of Phi_crit (configurable)"

patterns-established:
  - "Phi_crit check before solver: always compute and check before attempting device creation at high fluence"

requirements-completed: [CRMV-01, CRMV-02]

# Metrics
duration: 16min
completed: 2026-03-25
---

# Phase 16 Plan 01: Carrier Removal C-V Evolution Summary

**compute_phi_crit() for critical fluence detection and cv_at_fluence() for C-V curves under irradiation with Phi_crit safety checks**

## Performance

- **Duration:** 16 min
- **Started:** 2026-03-25T08:52:09Z
- **Completed:** 2026-03-25T09:08:00Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments

- compute_phi_crit() correctly returns ~4.86e13 protons/cm^2 for the Petringa graded profile at 62 MeV
- cv_at_fluence() produces valid C-V curves at fluence=0 and progressively flatter curves at moderate fluence
- cv_at_fluence() returns None with error log when fluence >= Phi_crit and warns at >= 90%
- All 63 tests pass including full regression suite (no regressions in existing tests)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add compute_phi_crit, cv_at_fluence, plot_cv_evolution** - `0deb473` (feat)
2. **Task 2: Integration tests for compute_phi_crit and cv_at_fluence** - `5ddcf80` (test)

## Files Created/Modified

- `src/radiation_damage.py` - Added compute_phi_crit() after compute_damaged_params()
- `src/cv_analysis.py` - Added cv_at_fluence() and plot_cv_evolution() after cv_sweep(); added imports for uuid, matplotlib, device, radiation_damage, poisson, drift_diffusion, sic_material
- `tests/test_radiation_damage.py` - Added TestComputePhiCrit class (4 tests)
- `tests/test_cv.py` - Added TestCvAtFluence class (4 integration tests with devsim)

## Decisions Made

- Phi_crit computed from min(N_D > 0) in the profile, not mean or bulk value -- ensures the weakest point is the reference for compensation
- cv_at_fluence returns None (not raises) when fluence >= Phi_crit, consistent with prior fluence-sweep functions that handle per-point failures gracefully
- 90% Phi_crit warning threshold is configurable via phi_crit_threshold parameter

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- cv_at_fluence() and plot_cv_evolution() ready for 16-02 publication notebook
- Phi_crit detection provides solver safety for any future fluence sweep near compensation
- C-V flattening physics validated: damaged device shows reduced capacitance spread vs pristine

---

_Phase: 16-carrier-removal-cv-evolution_
_Completed: 2026-03-25_
