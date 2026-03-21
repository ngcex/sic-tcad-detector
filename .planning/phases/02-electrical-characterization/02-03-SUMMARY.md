---
phase: 02-electrical-characterization
plan: 03
subsystem: validation
tags: [jupyter, notebook, iv-curve, cv-curve, devsim, graded-doping, validation]

# Dependency graph
requires:
  - phase: 02-electrical-characterization (plans 01-02)
    provides: drift-diffusion solver, I-V/C-V sweep, validation framework, plotting
provides:
  - Phase 2 validation notebook with end-to-end I-V and C-V characterization
  - Publication-quality figures (doping profile, C-V comparison, forward/reverse I-V)
affects: [phase-2-verification, phase-3-cce]

# Tech tracking
tech-stack:
  added: [nbformat, jupyter]
  patterns: [notebook-based validation workflow]

key-files:
  created:
    - notebooks/02_electrical_characterization.ipynb
    - figures/phase2_doping_profile.png
    - figures/phase2_cv_comparison.png
    - figures/phase2_forward_iv.png
    - figures/phase2_reverse_iv.png
    - figures/phase2_iv_validation.png

key-decisions:
  - "Graded doping defaults (N_D_junction=1e15, N_D_bulk=5e13, L_transition=200um) produce incorrect equilibrium W -- needs recalibration"
  - "Notebook artifacts complete but physics results do not meet validation targets -- gap closure required"

patterns-established:
  - "Validation notebook pattern: device creation, simulation, experimental overlay, metrics table"

requirements-completed: [ELEC-01, ELEC-02, VAL-01]

# Metrics
duration: 5min
completed: 2026-03-21
---

# Phase 2 Plan 03: Validation Notebook Summary

**End-to-end I-V and C-V validation notebook with graded doping, but graded doping calibration produces incorrect equilibrium depletion width cascading into unrealistic I-V -- gap closure needed**

## Performance

- **Duration:** ~5 min (task 1) + checkpoint review
- **Started:** 2026-03-21T08:46:31Z
- **Completed:** 2026-03-21T09:05:17Z
- **Tasks:** 2 (1 auto, 1 checkpoint)
- **Files modified:** 11 (1 notebook + 10 figure files)

## Accomplishments

- Created 9-cell validation notebook covering device creation, doping visualization, C-V simulation, C-V metrics, forward I-V, reverse I-V, I-V metrics, and summary assessment
- Generated publication-quality figures saved to figures/ directory (PDF + PNG pairs)
- Notebook executes end-to-end without errors -- pipeline is functional

## Critical Issues Found During Verification

The user executed the notebook and found significant physics/calibration problems:

### C-V Results

- **W(0V) = 10.0 um** instead of expected 1.7 um (488% error)
- **R-squared = -0.66** (negative, indicating model is worse than mean)
- Root cause: Epi layer is fully depleted at 0V with current graded doping parameters
- Graded doping defaults (N_D_junction=1e15, N_D_bulk=5e13, L_transition=200um) do not produce correct equilibrium depletion width

### I-V Results

- **Dark current = 6.7e-49 A** (unrealistically small, physically impossible)
- **Rectification ratio = 10.6** vs target 1e5
- **Series resistance = 0 Ohm** vs target ~3 kOhm
- Root cause: Cascading from incorrect C-V -- fully depleted epi leads to wrong carrier distributions

### Root Cause Analysis

The graded doping profile parameters are not calibrated to produce the correct equilibrium (0V) depletion width. When the epi is fully depleted at 0V, all downstream results (C-V agreement, I-V dark current, rectification) are wrong. This requires recalibrating the graded doping profile, likely adjusting N_D levels and/or transition length to match the 1.7 um equilibrium target before re-running the full characterization.

## Task Commits

Each task was committed atomically:

1. **Task 1: Phase 2 validation notebook** - `b6e0c13` (feat)
2. **Task 2: Human verification checkpoint** - No commit (checkpoint reviewed, issues found)

## Files Created/Modified

- `notebooks/02_electrical_characterization.ipynb` - 9-cell validation notebook with I-V, C-V, doping profile
- `figures/phase2_doping_profile.{png,pdf}` - Graded N_D(x) profile visualization
- `figures/phase2_cv_comparison.{png,pdf}` - C-V with experimental overlay
- `figures/phase2_forward_iv.{png,pdf}` - Forward bias I-V curve
- `figures/phase2_reverse_iv.{png,pdf}` - Reverse bias dark current
- `figures/phase2_iv_validation.{png,pdf}` - I-V validation metrics

## Decisions Made

- Proceed via gap closure rather than blocking here -- notebook infrastructure is sound, physics calibration needs work
- Requirements ELEC-01, ELEC-02, VAL-01 were already marked complete in plan 02-02 (the simulation and validation code exists); the notebook revealed that the results do not yet meet targets

## Deviations from Plan

None in notebook creation -- plan executed as written for Task 1.

The checkpoint (Task 2) was NOT approved due to physics calibration issues. The user directed completion via gap closure to address:

1. Graded doping parameter recalibration (N_D values and transition length)
2. Re-validation of C-V and I-V after recalibration

## Issues Encountered

- Graded doping defaults from plan 02-01 (N_D_junction=1e15, N_D_bulk=5e13, L_transition=200um) fully deplete the epi at 0V, producing W(0V)=10um instead of 1.7um
- This cascades into physically impossible I-V results (dark current ~1e-49 A)
- These are physics/calibration issues, not code bugs -- the simulation infrastructure works correctly

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**NOT ready for Phase 3.** Phase 2 verification must first close the graded doping calibration gap:

- Recalibrate graded doping profile to match W(0V)=1.7um equilibrium target
- Re-run I-V and C-V validation to confirm targets are met
- Only then can Phase 3 (CCE) build on correct electrical characterization

## Self-Check: PASSED

- FOUND: notebooks/02_electrical_characterization.ipynb
- FOUND: figures/phase2_cv_comparison.png
- FOUND: figures/phase2_doping_profile.png
- FOUND: 02-03-SUMMARY.md
- FOUND: commit b6e0c13

---

_Phase: 02-electrical-characterization_
_Completed: 2026-03-21_
