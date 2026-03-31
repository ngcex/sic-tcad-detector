---
phase: 22-monte-carlo-coupling
plan: 02
subsystem: simulation
tags:
  [
    monte-carlo,
    notebook,
    pulse-height-distribution,
    cce,
    visualization,
    microdosimetry,
  ]

requires:
  - phase: 22-monte-carlo-coupling
    plan: 01
    provides: "mc_coupling.py module: load_mc_events_csv, process_mc_ensemble, pulse_height_distribution"
  - phase: 21-single-particle-transient
    provides: "load_cce_let_table for CCE(LET) interpolation, cce_let_table_100um.json data"
provides:
  - "Publication-quality notebook 17 demonstrating full MC coupling pipeline"
  - "Synthetic mixed-field MC event dataset (2000 events, bimodal proton + heavy-ion)"
  - "4 figures: pulse height distribution, CCE vs LET overlay, deposited vs collected scatter, CCE comparison"
  - "Reproducible notebook generator script (scripts/create_notebook_17.py)"
affects: [phase-23, microdosimetric-spectra, publications]

tech-stack:
  added: []
  patterns:
    - "Notebook generator pattern: scripts/create_notebook_17.py generates + executes via nbconvert"
    - "Demo CCE curve for visual comparison when lookup table is flat near full depletion"
    - "Bimodal synthetic event generation: 70% proton-like (50 keV) + 30% heavy-ion-like (500 keV)"

key-files:
  created:
    - "scripts/create_notebook_17.py"
    - "notebooks/17_mc_coupling.ipynb"
    - "data/synthetic_mc_events.csv"
  modified: []

key-decisions:
  - "Added demo partial-depletion CCE curve to show CCE variation, since real CCE table is ~1.0 at full depletion"
  - "2000 synthetic events with bimodal distribution to mimic mixed proton + heavy-ion field"
  - "554 events show CCE < 0.95 with demo curve, mean CCE 0.9350 -- demonstrates pipeline sensitivity"

patterns-established:
  - "Demo CCE curve pattern: when real CCE is flat, add illustrative partial-depletion curve for visual comparison"

requirements-completed: [MCCP-04]

duration: ~15min
completed: 2026-03-31
---

# Phase 22 Plan 02: MC Coupling Notebook Summary

**Publication-quality notebook 17 with 4 figures demonstrating full MC coupling pipeline: synthetic event generation, CSV import, batch CCE lookup, and pulse height distribution for 2000 mixed-field events at 530k+ events/sec throughput**

## Performance

- **Duration:** ~15 min (across two sessions with checkpoint)
- **Started:** 2026-03-31T16:07:00Z
- **Completed:** 2026-03-31T16:27:41Z
- **Tasks:** 2 (+ 1 revision after checkpoint feedback)
- **Files created:** 3

## Accomplishments

- Notebook 17 with 20 cells and 4 publication-quality figures demonstrating the full MC coupling workflow
- Synthetic mixed-field dataset: 2000 events with bimodal proton/heavy-ion energy distribution saved to CSV
- Pipeline processes 530k+ events/sec using pre-computed CCE(LET) lookup table
- Demo CCE comparison showing real (full depletion ~1.0) vs partial-depletion curve with visible CCE variation

## Task Commits

Each task was committed atomically:

1. **Task 1: Create notebook 17 with MC coupling pipeline demonstration** - `636df24` (feat)
2. **Task 1-rev: Revise notebook 17 per checkpoint feedback** - `3cffeb2` (fix)
3. **Task 2: Human verification of notebook 17** - checkpoint approved, no commit

**Plan metadata:** (this commit)

## Files Created/Modified

- `scripts/create_notebook_17.py` - Generator script for notebook 17 (nbformat + nbconvert execution)
- `notebooks/17_mc_coupling.ipynb` - Executed notebook with 20 cells, 4 figures, summary statistics
- `data/synthetic_mc_events.csv` - 2000 synthetic MC events with bimodal energy distribution

## Decisions Made

- Added demo partial-depletion CCE curve alongside the real CCE lookup table. The real table gives CCE ~1.0 at 50V full depletion, so all events cluster at CCE=1. The demo curve (Gaussian dip centered at 5 keV/um) creates visible variation for illustrative purposes.
- Used 2000 events (70% proton-like at ~50 keV, 30% heavy-ion-like at ~500 keV) to produce a realistic bimodal pulse height distribution typical of mixed-field microdosimetry.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] Added demo CCE curve for visual differentiation**

- **Found during:** Task 1 revision (checkpoint feedback)
- **Issue:** Real CCE table at 50V full depletion gives CCE ~1.0 for all LET values, making CCE vs LET and deposited vs collected plots visually uninformative (all points on y=x line)
- **Fix:** Added a demo partial-depletion CCE curve (Gaussian dip at LET=5 keV/um, min CCE ~0.85) and a 4th figure comparing real vs demo CCE results
- **Files modified:** scripts/create_notebook_17.py, notebooks/17_mc_coupling.ipynb
- **Verification:** 554 events show CCE < 0.95 with demo curve, mean CCE 0.9350
- **Committed in:** 3cffeb2

---

**Total deviations:** 1 auto-fixed (1 missing critical visualization)
**Impact on plan:** Enhancement makes figures scientifically informative. No scope creep.

## Issues Encountered

None - pipeline executed correctly on first run; revision was driven by checkpoint feedback for visual quality.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- MC coupling pipeline fully demonstrated and validated in notebook form
- Synthetic event dataset available for downstream testing
- Ready for Phase 23: converting collected energies to lineal energy spectra (microdosimetric y\*d(y))
- SiC-specific kappa tissue-equivalence factor still needed before Phase 23 (documented blocker)

---

_Phase: 22-monte-carlo-coupling_
_Completed: 2026-03-31_

## Self-Check: PASSED

All files and commits verified:

- scripts/create_notebook_17.py: FOUND
- notebooks/17_mc_coupling.ipynb: FOUND
- data/synthetic_mc_events.csv: FOUND
- Commit 636df24: FOUND
- Commit 3cffeb2: FOUND
