---
phase: 16-carrier-removal-cv-evolution
plan: 02
subsystem: analysis
tags:
  [
    notebook,
    dark-current,
    cv-analysis,
    carrier-removal,
    mott-schottky,
    publication,
  ]

# Dependency graph
requires:
  - phase: 16-carrier-removal-cv-evolution-01
    provides: "cv_at_fluence, compute_phi_crit, plot_cv_evolution functions"
  - phase: 15-dark-current-fluence
    provides: "dark_current_vs_fluence and plot_dark_current_vs_fluence functions"
provides:
  - "Publication notebook 11 combining dark current + C-V evolution under proton irradiation"
  - "Four individual figure panels and combined 2x2 publication figure"
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns:
    [
      "4-panel combined figure layout for multi-observable radiation damage visualization",
    ]

key-files:
  created:
    - notebooks/11_dark_current_cv_evolution.ipynb
    - figures/11_dark_current_cv_evolution.png
    - figures/11a_dark_current_vs_fluence.png
    - figures/11b_cv_evolution.png
    - figures/11c_doping_profiles.png
  modified: []

key-decisions:
  - "Phi_crit ~4.86e13 protons/cm^2 annotated on dark current and C-V plots"
  - "Fluence levels [0, 1e12, 5e12, 1e13, 5e13] chosen to show C-V progression without hitting Phi_crit"

patterns-established:
  - "Multi-observable notebook pattern: combine related radiation damage observables in single notebook with shared parameters"

requirements-completed: [NBKV-03]

# Metrics
duration: 8min
completed: 2026-03-25
---

# Phase 16 Plan 02: Combined Dark Current + C-V Evolution Notebook Summary

**Publication notebook combining dark current vs fluence with C-V evolution, Mott-Schottky analysis, and position-dependent doping profiles under 62 MeV proton irradiation**

## Performance

- **Duration:** ~8 min (including human verification)
- **Started:** 2026-03-25T09:07:00Z
- **Completed:** 2026-03-25T09:15:05Z
- **Tasks:** 2 (1 auto + 1 human-verify)
- **Files created:** 5

## Accomplishments

- Created notebook 11 with 4-panel publication figure combining dark current and C-V observables
- Dark current panel reuses Phase 15 functions directly (no reimplementation)
- C-V panel uses Plan 01's cv_at_fluence() with progressive flattening visible across fluence levels
- Mott-Schottky (1/C^2 vs V) panel shows slope decrease confirming doping reduction
- Position-dependent doping profiles show bulk epi compensating first, junction region last
- Phi_crit (~4.86e13 protons/cm^2) annotated on relevant plots
- Human-verified as publication quality

## Task Commits

Each task was committed atomically:

1. **Task 1: Create combined dark current + C-V evolution notebook** - `3e5c79a` (feat)
2. **Task 2: Visual verification of publication notebook** - human-approved checkpoint

## Files Created/Modified

- `notebooks/11_dark_current_cv_evolution.ipynb` - Combined dark current + C-V evolution publication notebook
- `figures/11_dark_current_cv_evolution.png` - Combined 2x2 publication figure (300 DPI)
- `figures/11a_dark_current_vs_fluence.png` - Dark current vs fluence with component decomposition
- `figures/11b_cv_evolution.png` - C-V curves and Mott-Schottky plots at multiple fluences
- `figures/11c_doping_profiles.png` - Position-dependent doping profiles at selected fluences

## Decisions Made

- Phi_crit ~4.86e13 protons/cm^2 for Petringa device at 62 MeV, computed from minimum of graded N_D profile
- Fluence levels [0, 1e12, 5e12, 1e13, 5e13] selected to show C-V progression without exceeding Phi_crit

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 16 complete: carrier removal C-V evolution fully implemented and visualized
- All notebook requirements (NBKV-03) satisfied
- Ready for subsequent phases

---

_Phase: 16-carrier-removal-cv-evolution_
_Completed: 2026-03-25_

## Self-Check: PASSED

All files and commits verified present.
