---
phase: 14-cce-vs-fluence
plan: 02
subsystem: simulation
tags:
  [
    cce,
    fluence,
    radiation-damage,
    proton-irradiation,
    sensitivity-analysis,
    notebook,
    publication,
  ]

# Dependency graph
requires:
  - phase: 14-cce-vs-fluence-plan-01
    provides: "cce_vs_fluence() and cce_vs_bias_at_fluence() sweep functions"
provides:
  - "Publication-quality notebook 10 with CCE degradation curves under proton irradiation"
  - "4 figures: CCE vs fluence, multi-bias overlay, CCE vs bias at fixed fluence, sensitivity analysis"
  - "Summary table comparing linear and logarithmic lifetime models"
affects: [15-design-optimization, 16-validation, 18-documentation]

# Tech tracking
tech-stack:
  added: []
  patterns: [sensitivity-envelope-from-scaled-damage-constants]

key-files:
  created:
    - notebooks/10_cce_vs_fluence.ipynb
    - figures/10_cce_vs_fluence.png
    - figures/10_cce_vs_fluence_multibias.png
    - figures/10_cce_vs_bias_damaged.png
    - figures/10_cce_sensitivity.png
  modified: []

key-decisions:
  - "Reduced bias-sweep max fluence from 1e14 to 5e13 to stay within solver convergence regime (cce_vs_bias_at_fluence creates single device, no per-point error handling)"
  - "Used 8-point fluence grid for sensitivity analysis (vs 12 for main sweeps) to manage runtime"

patterns-established:
  - "Sensitivity envelope: scale all eta values (Z12, EH67, EH4, removal) uniformly by 0.5x-2.0x, fill_between for uncertainty bands"

requirements-completed: [NBKV-02]

# Metrics
duration: 15min
completed: 2026-03-24
---

# Phase 14 Plan 02: CCE vs Fluence Publication Notebook Summary

**Publication-quality notebook with CCE degradation curves at 5 bias voltages, bias recovery at 4 fluence levels, and linear vs logarithmic model comparison with 0.5x-2.0x uncertainty bands**

## Performance

- **Duration:** 15 min
- **Started:** 2026-03-24T22:05:35Z
- **Completed:** 2026-03-24T22:21:34Z
- **Tasks:** 1
- **Files created:** 5

## Accomplishments

- Created notebook 10 with 9 cells covering CCE degradation under 62 MeV proton irradiation
- Generated 4 publication-quality figures saved to figures/ directory
- CCE vs fluence at -40V shows monotonic degradation from ~1.0 to NaN (solver divergence above ~5e13 p/cm^2)
- Multi-bias plot (CCED-02) demonstrates higher bias partially compensating radiation damage
- CCE vs bias at fixed fluence (CCED-03) shows CCE recovery with increasing reverse bias
- Sensitivity analysis (NBKV-02) with side-by-side linear/logarithmic models and 0.5x-2.0x uncertainty envelope

## Task Commits

Each task was committed atomically:

1. **Task 1: Create CCE vs fluence publication notebook** - `3ffb1a1` (feat)

## Files Created/Modified

- `notebooks/10_cce_vs_fluence.ipynb` - Publication notebook with 9 cells: intro, imports, fluence grid, 4 figures, summary table, discussion
- `figures/10_cce_vs_fluence.png` - CCE vs fluence at -40V reference bias
- `figures/10_cce_vs_fluence_multibias.png` - CCE vs fluence at 5 bias voltages (-20V to -100V)
- `figures/10_cce_vs_bias_damaged.png` - CCE vs bias at pristine, 1e12, 1e13, 5e13 p/cm^2
- `figures/10_cce_sensitivity.png` - Two-panel linear vs logarithmic model comparison with uncertainty bands

## Decisions Made

- Reduced max fluence for bias-sweep figure from 1e14 to 5e13 p/cm^2: the `cce_vs_bias_at_fluence()` function creates a single device for the entire bias sweep, so if the device setup fails at full compensation the entire call crashes (unlike `cce_vs_fluence()` which catches per-point errors)
- Used 8-point fluence grid for sensitivity analysis to keep total runtime manageable (~150+ device solves total)
- Added try/except wrapping around `cce_vs_bias_at_fluence()` calls in case solver diverges at boundary fluences

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Reduced bias-sweep fluence from 1e14 to 5e13**

- **Found during:** Task 1 (Notebook execution)
- **Issue:** `cce_vs_bias_at_fluence(fluence=1e14)` crashes because full doping compensation causes Poisson equilibrium solver to diverge during device setup (not caught by per-point error handling)
- **Fix:** Changed highest fluence from 1e14 to 5e13 p/cm^2 and added try/except wrapper around each `cce_vs_bias_at_fluence()` call
- **Files modified:** notebooks/10_cce_vs_fluence.ipynb
- **Verification:** Notebook executes end-to-end without errors; all 4 figures generated
- **Committed in:** 3ffb1a1 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug - solver convergence boundary)
**Impact on plan:** Physics correctly demonstrated within solver convergence regime. The 5e13 fluence point shows significant damage and is near the convergence boundary. Full-compensation handling deferred to Phase 16.

## Issues Encountered

- Newton solver diverges at 1e14 p/cm^2 for `cce_vs_bias_at_fluence()` -- known limitation (full doping compensation). Handled by reducing fluence to 5e13.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 14 complete: both fluence sweep infrastructure (Plan 01) and publication notebook (Plan 02) delivered
- All 4 CCED requirements met, NBKV-02 delivered
- Solver convergence limit near full compensation remains documented for Phase 16
- Notebook and figures ready for Phase 18 documentation integration

---

_Phase: 14-cce-vs-fluence_
_Completed: 2026-03-24_
