---
phase: 02-electrical-characterization
plan: 02
subsystem: physics
tags:
  [
    4H-SiC,
    I-V-sweep,
    C-V-analysis,
    validation,
    depletion-width,
    capacitance,
    Mott-Schottky,
    experimental-comparison,
  ]

# Dependency graph
requires:
  - phase: 02-electrical-characterization (plan 01)
    provides: drift-diffusion solver with SRH recombination, contact current extraction, graded doping
provides:
  - I-V sweep function with incremental voltage ramping and convergence fallback
  - C-V analysis from depletion width (parallel-plate capacitance model)
  - Validation framework with Petringa experimental targets and agreement metrics
  - I-V and C-V plotting functions with experimental data overlay
affects:
  [02-electrical-characterization (plan 03), 03-charge-collection-efficiency]

# Tech tracking
tech-stack:
  added: []
  patterns:
    [
      incremental voltage ramping with fallback convergence for I-V sweep,
      parallel-plate C = eps*A/W junction capacitance model,
      R-squared and RMSE agreement metrics without sklearn,
      Mott-Schottky 1/C^2 vs V analysis for doping extraction,
    ]

key-files:
  created:
    - src/cv_analysis.py
    - src/validation.py
    - tests/test_cv.py
    - tests/test_validation.py
  modified:
    - src/drift_diffusion.py
    - src/plotting.py

key-decisions:
  - "iv_sweep ramps incrementally from current bias state (not from 0V each time) for convergence stability"
  - "C-V uses parallel-plate approximation C=eps*A/W (exact for abrupt junction, good approximation for graded)"
  - "Validation pass/fail thresholds allow 2 orders of magnitude tolerance (simulation vs measurement)"

patterns-established:
  - "I-V sweep: iv_sweep() with V_step_forward=0.1V, V_step_reverse=0.5V, fallback to relaxed tolerances"
  - "C-V from W(V): compute_cv_from_depletion() for analytical, cv_sweep() for numerical"
  - "Agreement metrics: compute_agreement_metrics() returns r_squared, rmse, max_deviation, relative errors"
  - "Plotting: plot_iv_comparison/plot_cv_comparison with experimental overlay and metric annotations"

requirements-completed: [ELEC-01, ELEC-02, VAL-01]

# Metrics
duration: 4min
completed: 2026-03-21
---

# Phase 2 Plan 02: I-V Sweep, C-V Analysis, and Validation Summary

**I-V sweep with convergence ramping, C-V from depletion width, and Petringa experimental validation metrics with 18 new tests**

## Performance

- **Duration:** 4 min
- **Started:** 2026-03-21T08:42:48Z
- **Completed:** 2026-03-21T08:46:31Z
- **Tasks:** 2
- **Files created/modified:** 6

## Accomplishments

- I-V sweep function with incremental voltage ramping and dual-tolerance convergence fallback
- C-V analysis module computing junction capacitance from depletion width (parallel-plate model) with Mott-Schottky 1/C^2 support
- Validation framework with Petringa experimental targets (dark current, rectification ratio, W(V) at three bias points) and comprehensive agreement metrics
- Four publication-quality plotting functions for I-V and C-V with experimental data overlay and annotation
- 18 new unit tests (8 C-V + 10 validation), 96 total passing

## Task Commits

Each task was committed atomically:

1. **Task 1: I-V sweep, C-V analysis, and validation modules** - `5c60f31` (feat)
2. **Task 2: I-V and C-V plotting functions** - `79f6d5e` (feat)

## Files Created/Modified

- `src/drift_diffusion.py` - Added iv_sweep() with incremental ramping and fallback convergence (349 lines total)
- `src/cv_analysis.py` - New: junction_capacitance, depletion_width_from_capacitance, compute_cv_from_depletion, cv_sweep (177 lines)
- `src/validation.py` - New: EXPERIMENTAL_IV/CV targets, compute_agreement_metrics, validate_iv, validate_cv (207 lines)
- `src/plotting.py` - Added plot_iv_curve, plot_iv_comparison, plot_cv_curve, plot_cv_comparison (489 lines total)
- `tests/test_cv.py` - New: 8 tests for capacitance computation and round-trip conversion (95 lines)
- `tests/test_validation.py` - New: 10 tests for agreement metrics and experimental data constants (107 lines)

## Decisions Made

- iv_sweep tracks current bias state and ramps incrementally to each target voltage rather than resetting to 0V, which improves convergence for large sweeps
- Validation thresholds use 2 orders of magnitude tolerance (100x) for dark current and rectification ratio, reflecting expected simulation-vs-measurement discrepancy
- C-V uses parallel-plate approximation which is exact for abrupt junctions and a good approximation for the graded doping profile

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - all modules implemented cleanly, all tests passed on first run.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- I-V sweep ready for full forward+reverse characterization with dd_device
- C-V sweep ready for numerical W(V) extraction and comparison with Petringa data
- Validation framework ready for automated pass/fail checking in Plan 03 integration
- All plotting functions available for publication-quality figure generation
- 96 tests passing across full suite (78 Phase 1/2.01 + 18 new)

## Self-Check: PASSED

All 6 files verified on disk. Both task commits (5c60f31, 79f6d5e) verified in git history. cv_analysis.py: 177 lines (>40 min). validation.py: 207 lines (>60 min). test_cv.py: 95 lines (>30 min). test_validation.py: 107 lines (>40 min). All required exports verified. 96 tests passing.

---

_Phase: 02-electrical-characterization_
_Completed: 2026-03-21_
