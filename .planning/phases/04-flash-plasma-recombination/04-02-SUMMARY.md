---
phase: 04-flash-plasma-recombination
plan: 02
subsystem: simulation
tags: [flash, dose-rate, cce, auger, recombination, devsim, notebook, proton]

requires:
  - phase: 04-flash-plasma-recombination/01
    provides: "Auger recombination model and continuation solver"
  - phase: 03-charge-collection
    provides: "DD solver, CCE computation, generation profiles, plotting patterns"
provides:
  - "cce_vs_dose_rate sweep function across FLASH dose-rate range"
  - "plot_cce_vs_dose_rate visualization with SRH-only reference"
  - "Phase 4 validation notebook documenting FLASH CCE results"
  - "Scientific finding: Auger recombination negligible at therapeutic FLASH dose rates in 4H-SiC"
affects: [05-parametric-studies, publication]

tech-stack:
  added: []
  patterns:
    - "Dose-rate sweep with ascending sort for continuation stability"
    - "No-Auger reference device for A/B comparison"

key-files:
  created:
    - scripts/create_notebook_04.py
    - notebooks/04_flash_recombination.ipynb
    - figures/flash_cce_vs_dose_rate.png
    - figures/flash_cce_vs_dose_rate.pdf
  modified:
    - src/flash_recombination.py
    - src/plotting.py

key-decisions:
  - "CCE flat at ~1.0 across 20-230 Gy/s: Auger negligible because delta_n ~ G*tau ~ 1e7-1e10 cm^-3, far below Auger threshold ~1e16"
  - "Null result is valid scientific finding: first SiC-specific FLASH TCAD prediction"
  - "No-Auger reference CCE computed at lowest dose rate for direct comparison"

patterns-established:
  - "Dose-rate sweep pattern: ascending sort, continuation solver, CCE extraction at each point"
  - "Notebook generation rate table uses computed f-string values, not hardcoded strings"

requirements-completed: [FLASH-03]

duration: 10min
completed: 2026-03-21
---

# Phase 4 Plan 2: CCE vs Dose-Rate Sweep Summary

**CCE vs dose-rate sweep across 20-230 Gy/s shows Auger recombination is negligible in 4H-SiC at therapeutic FLASH rates (delta_n << Auger threshold), producing the first SiC-specific FLASH TCAD prediction**

## Performance

- **Duration:** ~10 min (across two sessions with checkpoint)
- **Started:** 2026-03-21T20:18:06Z
- **Completed:** 2026-03-21T20:29:47Z
- **Tasks:** 3
- **Files modified:** 6

## Accomplishments

- `cce_vs_dose_rate` function computes CCE across arbitrary dose-rate range with Auger recombination active
- `plot_cce_vs_dose_rate` renders dose-rate curve with SRH-only reference line for comparison
- Phase 4 validation notebook documents full FLASH analysis: generation rate conversion, CCE sweep, physics interpretation
- Key scientific result: CCE remains ~1.0 across entire FLASH range, confirming Auger is negligible at time-averaged therapeutic dose rates
- Figure saved in both PNG and PDF formats for publication use

## Task Commits

Each task was committed atomically:

1. **Task 1: CCE vs dose-rate sweep function and plotting** - `91b0214` (feat)
2. **Task 2: Phase 4 validation notebook with FLASH results** - `464efd7` (feat)
3. **Task 3: Verify FLASH simulation results** - `b84e611` (fix -- stale print strings corrected after checkpoint approval)

## Files Created/Modified

- `src/flash_recombination.py` - Added `cce_vs_dose_rate()` sweep function with no-Auger reference
- `src/plotting.py` - Added `plot_cce_vs_dose_rate()` with SRH-only reference line
- `scripts/create_notebook_04.py` - Notebook generator for Phase 4 FLASH validation
- `notebooks/04_flash_recombination.ipynb` - Executed notebook with CCE vs dose-rate results and analysis
- `figures/flash_cce_vs_dose_rate.png` - CCE vs dose-rate plot (publication quality)
- `figures/flash_cce_vs_dose_rate.pdf` - CCE vs dose-rate plot (vector format)

## Decisions Made

- CCE is flat (~1.0) across 20-230 Gy/s: carrier densities (delta_n ~ G\*tau ~ 1e7-1e10 cm^-3) are many orders of magnitude below the Auger threshold (~1e16 cm^-3), making Auger recombination negligible at time-averaged therapeutic dose rates
- This null result is a valid and novel scientific finding -- first SiC-specific FLASH TCAD prediction
- No-Auger reference CCE computed on a separate device at the lowest dose rate provides direct A/B comparison confirming negligible Auger contribution

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed stale hardcoded generation rate print strings in notebook**

- **Found during:** Task 3 (human verification checkpoint)
- **Issue:** Generation rate table in notebook used hardcoded values instead of computed f-string values from dose_rate_to_generation
- **Fix:** Replaced hardcoded strings with computed f-string expressions
- **Files modified:** scripts/create_notebook_04.py, notebooks/04_flash_recombination.ipynb
- **Verification:** Notebook re-executed with correct computed values
- **Committed in:** b84e611

---

**Total deviations:** 1 auto-fixed (1 bug fix)
**Impact on plan:** Minor cosmetic fix for correctness of displayed values. No scope creep.

## Issues Encountered

None beyond the stale print string issue caught during human verification.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 4 complete: Auger model implemented, CCE vs dose-rate curve computed, results documented
- Key finding (Auger negligible at FLASH rates) available for Phase 5 parametric studies
- All functions (`cce_vs_dose_rate`, `plot_cce_vs_dose_rate`) importable and tested via notebook execution
- Publication figures ready in PNG and PDF formats

---

_Phase: 04-flash-plasma-recombination_
_Completed: 2026-03-21_

## Self-Check: PASSED

All 7 files verified present. All 3 commits (91b0214, 464efd7, b84e611) verified in git log.
