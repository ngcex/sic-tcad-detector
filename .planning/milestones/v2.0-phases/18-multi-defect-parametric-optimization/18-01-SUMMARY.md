---
phase: 18-multi-defect-parametric-optimization
plan: 01
subsystem: radiation-damage
tags: [parametric-sweep, uncertainty-envelope, single-defect, device-geometry]

# Dependency graph
requires:
  - phase: 14-radiation-damage-cce
    provides: cce_vs_fluence, cce_vs_bias_at_fluence, compute_damaged_params
  - phase: 15-dark-current-radiation
    provides: dark_current_vs_fluence
  - phase: 16-cv-radiation
    provides: cv_at_fluence
provides:
  - Parameterized N_D_junction/N_D_bulk/L_transition in cce_vs_fluence, cce_vs_bias_at_fluence, dark_current_vs_fluence, cv_at_fluence
  - make_single_defect_params() for K_tau-matched single-defect model
  - cce_uncertainty_envelope() for min/max CCE from per-defect eta scatter
  - radiation_hardness_sweep() for ranked device geometry optimization
affects: [18-02, 18-03, notebooks]

# Tech tracking
tech-stack:
  added: []
  patterns:
    [
      parametric-geometry-kwargs,
      lazy-import-for-devsim-isolation,
      K_tau-matching-construction,
    ]

key-files:
  created: []
  modified:
    - src/charge_collection.py
    - src/dark_current.py
    - src/cv_analysis.py
    - src/radiation_damage.py
    - tests/test_radiation_damage.py
    - tests/test_charge_collection.py

key-decisions:
  - "Geometry kwargs placed after epi_thickness_cm in signatures (before alpha_range_cm/energy_MeV)"
  - "Near-zero eta (1e-10) not zero for disabled defects to satisfy __post_init__ validation"
  - "cce_uncertainty_envelope and radiation_hardness_sweep use lazy imports to keep radiation_damage.py devsim-free"

patterns-established:
  - "Parameterized geometry: N_D_junction/N_D_bulk/L_transition kwargs with backward-compatible defaults"
  - "K_tau-matching: single-defect sigma derived from K_tau / (eta_eff * v_th)"

requirements-completed: [PARM-01, PARM-02, PARM-03]

# Metrics
duration: 19min
completed: 2026-03-26
---

# Phase 18 Plan 01: Parametric Optimization Infrastructure Summary

**Parameterized device geometry in sweep functions and added single-defect constructor, CCE uncertainty envelope, and radiation hardness sweep infrastructure**

## Performance

- **Duration:** 19 min
- **Started:** 2026-03-26T00:05:06Z
- **Completed:** 2026-03-26T00:24:00Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments

- Added N_D_junction, N_D_bulk, L_transition kwargs to all four fluence sweep functions (backward compatible)
- Implemented make_single_defect_params() with verified K_tau matching for both carriers
- Implemented cce_uncertainty_envelope() generating 8-combination eta scatter bounds
- Implemented radiation_hardness_sweep() returning ranked DataFrame of design parameter sweeps
- 10 new tests (4 unit + 6 integration) all passing

## Task Commits

Each task was committed atomically:

1. **Task 1: Parameterize device geometry in fluence sweep functions** - `2d8171c` (feat)
2. **Task 2: Add single-defect constructor, uncertainty envelope, and parametric sweep functions with tests** - `a084f8c` (feat)

## Files Created/Modified

- `src/charge_collection.py` - Added N_D_junction/N_D_bulk/L_transition to cce_vs_fluence, cce_vs_bias_at_fluence
- `src/dark_current.py` - Added N_D_junction/N_D_bulk/L_transition to dark_current_vs_fluence
- `src/cv_analysis.py` - Added N_D_junction/N_D_bulk/L_transition to cv_at_fluence
- `src/radiation_damage.py` - Added make_single_defect_params, cce_uncertainty_envelope, radiation_hardness_sweep
- `tests/test_radiation_damage.py` - Added TestMakeSingleDefectParams (4 tests)
- `tests/test_charge_collection.py` - Added TestParametricFunctions (2 tests)

## Decisions Made

- Geometry kwargs placed after epi_thickness_cm in all function signatures for consistent ordering
- Near-zero eta values (1e-10) used for disabled defects in single-defect model to satisfy RadiationDamageParams validation
- cce_uncertainty_envelope and radiation_hardness_sweep placed in radiation_damage.py with lazy imports to maintain devsim-free module-level imports (AST test passes)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All parametric functions ready for Phase 18 Plan 02 (uncertainty band notebooks)
- Device geometry parameterization enables Phase 18 Plan 03 (design optimization study)
- All existing tests pass unchanged (backward compatibility confirmed)

---

_Phase: 18-multi-defect-parametric-optimization_
_Completed: 2026-03-26_
