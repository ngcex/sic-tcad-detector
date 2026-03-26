---
phase: 18-multi-defect-parametric-optimization
plan: 02
subsystem: radiation-damage
tags: [notebook, multi-defect, uncertainty-bands, parametric-sweep, heatmap]

# Dependency graph
requires:
  - phase: 18-multi-defect-parametric-optimization
    provides: make_single_defect_params, cce_uncertainty_envelope, radiation_hardness_sweep, parameterized geometry
provides:
  - Notebook 12: single-vs-multi-defect CCE/dark current/C-V comparison with uncertainty bands
  - Notebook 13: parametric radiation hardness optimization with ranked table and heatmaps
affects: [18-03, publications]

# Tech tracking
tech-stack:
  added: []
  patterns:
    [
      nbformat-notebook-creation,
      fill-between-uncertainty-envelope,
      pivot-heatmap-visualization,
    ]

key-files:
  created:
    - notebooks/12_multi_defect_comparison.ipynb
    - notebooks/13_parametric_optimization.ipynb
  modified: []

key-decisions:
  - "Notebook 12 uses serif fonts, semilog-x, and 150 DPI matching project conventions"
  - "Notebook 13 designed for offline execution (~15-20 min sweep) with no auto-execution"

patterns-established:
  - "Multi-panel heatmap: pivot_table + imshow with annotated cell values for parameter sweeps"

requirements-completed: [PARM-01, PARM-02, PARM-03]

# Metrics
duration: 3min
completed: 2026-03-26
---

# Phase 18 Plan 02: Multi-Defect Comparison and Parametric Optimization Notebooks Summary

**Publication-quality notebooks comparing single-vs-three-defect models with uncertainty bands and 64-point parametric radiation hardness sweep with heatmap visualization**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-26T00:30:26Z
- **Completed:** 2026-03-26T00:33:21Z
- **Tasks:** 2
- **Files created:** 2

## Accomplishments

- Created notebook 12 with 10 cells: single-vs-three-defect CCE/dark current/C-V overlays, model equivalence discussion, per-defect uncertainty bands at 4 bias voltages, and summary table
- Created notebook 13 with 10 cells: 64-point parametric sweep (4 epi x 4 doping x 4 bias), ranked table, 2x2 heatmap, effect-of-parameter line plots, and design recommendations
- All code cells parse as valid Python AST
- Both notebooks follow project conventions (serif fonts, figure saving, labeled axes)

## Task Commits

Each task was committed atomically:

1. **Task 1: Multi-defect comparison and uncertainty bands notebook** - `558ddd9` (feat)
2. **Task 2: Parametric radiation hardness optimization notebook** - `e544d68` (feat)

## Files Created/Modified

- `notebooks/12_multi_defect_comparison.ipynb` - Single-vs-three-defect model comparison with CCE/dark current/C-V overlays, uncertainty bands, and summary table
- `notebooks/13_parametric_optimization.ipynb` - Parametric sweep of epi thickness, doping, and bias voltage with ranked table and heatmap visualization

## Decisions Made

- Notebook 12 uses serif fonts, semilog-x for fluence plots, and 150 DPI consistent with all prior notebooks
- Notebook 13 designed for user-initiated execution due to ~15-20 min runtime for 64-point sweep
- Heatmap uses pivot_table + imshow with annotated cell values for clear parameter comparison

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Both notebooks ready for execution; notebook 13 requires ~15-20 min for the parametric sweep
- All figures will be saved to figures/ directory on execution
- Phase 18 Plan 03 (final verification/documentation) can proceed

---

## Self-Check: PASSED

All files and commits verified.

_Phase: 18-multi-defect-parametric-optimization_
_Completed: 2026-03-26_
