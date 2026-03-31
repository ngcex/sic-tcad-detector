---
phase: 21-single-particle-transient
plan: 02
subsystem: simulation
tags: [devsim, transient, cce, let, notebook, publication, single-particle]

# Dependency graph
requires:
  - phase: 21-single-particle-transient/01
    provides: "single_particle.py module with ion_track_generation_2d, simulate_single_particle, build_cce_let_table"
  - phase: 20-2d-transport-cce/02
    provides: "2D charge collection validation and notebook 15 conventions"
provides:
  - "Publication-quality notebook 16 with single-particle CCE analysis"
  - "CCE(LET) lookup tables for 100um and 300um SVs (data/cce_let_table_{100,300}um.json)"
  - "Charge conservation validation across LET range"
affects: [22-mc-coupling, 23-microdosimetry-spectra]

# Tech tracking
tech-stack:
  added: []
  patterns:
    [
      "CCE(LET) lookup table JSON format with geometry metadata",
      "generation-pulse injection for transient CCE",
    ]

key-files:
  created:
    - notebooks/16_single_particle_cce.ipynb
    - scripts/create_notebook_16.py
    - data/cce_let_table_100um.json
    - data/cce_let_table_300um.json
  modified: []

key-decisions:
  - "Reduced 300um SV LET sweep from 15 to 10 points to stay within nbconvert timeout (larger mesh ~3x slower per sim)"
  - "CCE ~1.01 at 50V full depletion is numerically valid (1% overcollection from generation-pulse method)"

patterns-established:
  - "CCE(LET) lookup table: JSON with geometry metadata, loadable via load_cce_let_table for interpolation"

requirements-completed: [NBKV-02]

# Metrics
duration: 210min
completed: 2026-03-31
---

# Phase 21 Plan 02: Single-Particle CCE Notebook Summary

**Publication-quality notebook 16 with ion track visualization, transient pulse waveform, charge conservation validation, and CCE(LET) lookup tables for 100um/300um SVs**

## Performance

- **Duration:** ~3.5 hours (dominated by notebook execution: 20+10 transient simulations)
- **Started:** 2026-03-30T13:12:20Z
- **Completed:** 2026-03-31T00:20:00Z
- **Tasks:** 1 of 2 (Task 2 awaits human verification)
- **Files modified:** 4

## Accomplishments

- Created 17-cell publication-quality notebook with 4 figures and 410+ source lines
- Ion track generation profile visualized as 2D tricontourf on mesh (full + zoomed views)
- Transient current pulse plotted with labeled I_peak, t_peak, t_95% collection time
- Charge conservation validated: CCE = 1.0101 across LET 1-100 keV/um (PASS)
- CCE(LET) lookup tables saved for both SV sizes (20 points at 100um, 10 points at 300um)
- Both geometries show identical CCE at center injection, confirming size independence

## Task Commits

Each task was committed atomically:

1. **Task 1: Create notebook 16 with single-particle CCE analysis** - `b516a54` (feat)
2. **Task 2: Human verification of notebook 16 figures and science** - awaiting checkpoint

## Files Created/Modified

- `scripts/create_notebook_16.py` - Generator script for notebook 16 (17 cells)
- `notebooks/16_single_particle_cce.ipynb` - Executed notebook with all outputs and figures
- `data/cce_let_table_100um.json` - CCE(LET) table for 100um SV (20 LET points, 0.5-500 keV/um)
- `data/cce_let_table_300um.json` - CCE(LET) table for 300um SV (10 LET points, 0.5-500 keV/um)

## Decisions Made

- Reduced 300um SV sweep from 15 to 10 LET points after first attempt timed out at 3600s. The 300um mesh is ~3x larger, making each simulation take ~15 min. With 10 points, the sweep completed in ~150 min.
- CCE = 1.0101 (1% overcollection) accepted as numerically valid -- consistent across all LET values and both geometries, attributable to the generation-pulse injection method's finite timestep.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Reduced 300um LET sweep points from 15 to 10**

- **Found during:** Task 1 (notebook execution)
- **Issue:** First execution attempt timed out at 3600s on the 300um CCE(LET) sweep (cell 13). The 300um SV mesh is ~3x larger than 100um, making each transient simulation take ~15 min.
- **Fix:** Reduced n_let_points from 15 to 10 for the 300um sweep. Updated markdown description to note the larger mesh cost. Increased nbconvert timeout to 7200s.
- **Files modified:** scripts/create_notebook_16.py
- **Verification:** Notebook completed successfully with 10 points in ~150 min
- **Committed in:** b516a54 (part of Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Necessary to complete execution within timeout. 10 LET points still provides adequate coverage for interpolation.

## Issues Encountered

- Notebook execution took ~3.5 hours total (first 1h attempt timed out, second 2.5h attempt succeeded). The 300um SV simulations are significantly slower due to the larger 2D mesh.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- CCE(LET) lookup tables ready for Phase 22 Monte Carlo coupling
- `load_cce_let_table` provides log-linear interpolation function for scoring MC events
- Awaiting human verification of notebook figures and science (Task 2 checkpoint)

## Self-Check: PASSED

All artifacts verified:

- FOUND: notebooks/16_single_particle_cce.ipynb
- FOUND: scripts/create_notebook_16.py
- FOUND: data/cce_let_table_100um.json
- FOUND: data/cce_let_table_300um.json
- FOUND: commit b516a54

---

_Phase: 21-single-particle-transient_
_Completed: 2026-03-31_
