---
phase: 15-dark-current-vs-fluence
plan: 02
subsystem: simulation
tags:
  [
    dark-current,
    fluence-sweep,
    notebook,
    publication,
    delta-J,
    component-decomposition,
  ]

# Dependency graph
requires:
  - phase: 15-dark-current-vs-fluence
    provides: dark_current_vs_fluence() and plot_dark_current_vs_fluence() functions
provides:
  - Publication-quality dark current vs fluence notebook (05_dark_current_vs_fluence.ipynb)
  - Component decomposition visualization (SRH, TAT, SRV vs fluence)
  - Delta-J analysis plot (radiation-induced increase above baseline)
  - Multi-bias comparison (V = -10, -30, -50 V)
affects: [dark current analysis, publication figures]

# Tech tracking
tech-stack:
  added: []
  patterns:
    [delta-J notebook structure with baseline annotation and multi-bias overlay]

key-files:
  created:
    - notebooks/05_dark_current_vs_fluence.ipynb
  modified: []

key-decisions:
  - "Notebook numbered 05 per plan specification despite existing notebooks going up to 10"
  - "Multi-bias comparison uses V = -10, -30, -50 V with reduced fluence grid (10 points) for faster computation"

patterns-established:
  - "Dark current notebook structure: sweep -> component plot -> delta-J -> multi-bias -> discussion -> summary table"

requirements-completed: [DCRR-02]

# Metrics
duration: 2min
completed: 2026-03-25
---

# Phase 15 Plan 02: Dark Current vs Fluence Publication Notebook Summary

**Publication notebook with component-decomposed dark current vs fluence, delta-J analysis, and multi-bias comparison for 4H-SiC detector**

## Performance

- **Duration:** 2 min
- **Started:** 2026-03-25T07:54:59Z
- **Completed:** 2026-03-25T07:56:39Z
- **Tasks:** 1
- **Files created:** 1

## Accomplishments

- 8-cell publication notebook covering introduction, fluence sweep, component-decomposed plot, delta-J analysis, multi-bias comparison, discussion, and summary table
- Delta-J analysis separates radiation-induced increase from pristine baseline with log-log visualization
- Multi-bias comparison at V = -10, -30, -50 V shows voltage dependence of dark current evolution
- Discussion section explains counterintuitive SiC high-fluence behavior vs silicon monotonic increase
- Summary table provides full component decomposition (I_SRH, I_TAT, I_SRV, delta_I) at all fluence points

## Task Commits

Each task was committed atomically:

1. **Task 1: Create dark current vs fluence publication notebook** - `c3ad5fc` (feat)

## Files Created/Modified

- `notebooks/05_dark_current_vs_fluence.ipynb` - Publication-quality dark current vs fluence notebook with 8 cells: intro, imports, sweep, component plot, delta-J, multi-bias, discussion, summary table

## Decisions Made

- Notebook numbered 05 per plan specification (existing notebooks go up to 10, but plan explicitly names 05_dark_current_vs_fluence.ipynb)
- Multi-bias comparison uses reduced fluence grid (10 irradiated points vs 12) for faster computation while maintaining curve quality

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 15 (dark current vs fluence) fully complete: both infrastructure (15-01) and publication notebook (15-02) delivered
- Dark current vs fluence analysis ready for integration with broader radiation damage study
- Figures saved to figures/ directory upon notebook execution

---

_Phase: 15-dark-current-vs-fluence_
_Completed: 2026-03-25_
