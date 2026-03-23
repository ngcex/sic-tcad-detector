---
phase: 11-dark-current-modeling
plan: 02
subsystem: physics
tags:
  [
    dark-current,
    sensitivity-analysis,
    visualization,
    jupyter,
    decomposition,
    parameter-sweep,
  ]

# Dependency graph
requires:
  - phase: 11-dark-current-modeling
    plan: 01
    provides: TAT dark current model, component extraction, device creation
provides:
  - sensitivity_sweep for parameter studies (epi, doping, SRV, N_t)
  - plot_dark_current_decomposition for I-V visualization
  - Notebook 07 with dark current analysis workflow
affects: [12-transient-pulse, validation, documentation]

# Tech tracking
tech-stack:
  added: [pandas]
  patterns:
    - Parameter sweep with device creation/cleanup per point
    - Publication-quality decomposition plots (semilogy, multi-component)
    - Sensitivity DataFrame output for downstream analysis

key-files:
  created:
    - notebooks/07_dark_current.ipynb
  modified:
    - src/dark_current.py

key-decisions:
  - "Used pandas DataFrame for sensitivity_sweep output (consistent with temperature_sweep patterns)"
  - "Device cleanup via devsim.delete_device in finally block prevents resource leaks during sweeps"
  - "Notebook uses semilogx for SRV subplot since S=0 is a valid data point (shifted to 1 on log axis)"

patterns-established:
  - "sensitivity_sweep pattern: create device, ramp, extract, cleanup per parameter point"
  - "plot_dark_current_decomposition: reusable plotting function for I-V decomposition"

requirements-completed: [DARK-04, DARK-05, NOTE-02]

# Metrics
duration: 9min
completed: 2026-03-23
---

# Phase 11 Plan 02: Dark Current Visualization and Sensitivity Summary

**Sensitivity sweep utility and notebook 07 with dark current decomposition plots and epi/doping/SRV parameter studies**

## Performance

- **Duration:** 9 min
- **Started:** 2026-03-23T21:12:22Z
- **Completed:** 2026-03-23T21:21:42Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Added sensitivity_sweep function for parameter studies (epi thickness, N_D, N_t, S_n/S_p) with per-point device creation and cleanup
- Added plot_dark_current_decomposition for publication-quality log-scale I-V decomposition plots
- Created notebook 07 (17 cells, 6 sections) covering full dark current analysis: setup, I-V decomposition, calibration check, sensitivity studies, and summary
- Zero regression: 203 existing tests pass

## Task Commits

1. **Task 1: Add sensitivity_sweep and plot_dark_current_decomposition** - `db59fac` (feat)
2. **Task 2: Create Jupyter notebook 07** - `49479ca` (feat)

## Files Created/Modified

- `src/dark_current.py` - Added sensitivity_sweep (parameter studies) and plot_dark_current_decomposition (visualization)
- `notebooks/07_dark_current.ipynb` - 17-cell notebook with dark current analysis workflow

## Decisions Made

1. **pandas DataFrame output for sensitivity_sweep:** Consistent with temperature_sweep module patterns (10-02). Provides clean tabular output for display and further analysis.

2. **Device cleanup in finally block:** Each sweep point creates a fresh devsim device. The finally-block cleanup pattern (matching charge_collection.py and temperature_sweep.py) prevents resource leaks even on solver convergence failures.

3. **Semilogx for SRV subplot:** S=0 is a physically meaningful data point (no surface recombination). Shifted to 1 on the log-x axis for display while keeping the S=0 data point visible.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 11 (dark current modeling) is fully complete
- Dark current module provides: device creation, I-V sweep, component decomposition, sensitivity analysis, and visualization
- Ready for Phase 12 (transient pulse response) which will build on the established device models

---

_Phase: 11-dark-current-modeling_
_Completed: 2026-03-23_
