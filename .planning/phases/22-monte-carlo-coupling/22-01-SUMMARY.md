---
phase: 22-monte-carlo-coupling
plan: 01
subsystem: simulation
tags: [monte-carlo, csv, root, uproot, cce, pulse-height, geant4, let]

requires:
  - phase: 21-single-particle-transient
    provides: "ion_track_generation_2d for mesh charge profiles, load_cce_let_table for CCE interpolation"
provides:
  - "CSV and ROOT MC event importers with configurable column/branch mapping"
  - "Unit conversion (mm/MeV -> cm/keV) for Geant4 data"
  - "Per-event mesh charge profile generation via ion_track_generation_2d (full path)"
  - "Batch CCE(LET) lookup for 1000+ events in <1 second (fast path)"
  - "Log-spaced pulse height distribution histogram"
affects: [22-02, phase-23, notebooks, microdosimetric-spectra]

tech-stack:
  added: []
  patterns:
    [
      "Lazy uproot import: ROOT functions import uproot inside function body so CSV-only use works without uproot",
      "Standardized event DataFrame: all importers produce event_id/x_cm/y_cm/z_cm/edep_keV columns",
      "column_map convention: keys=standard names, values=source column names (consistent across CSV and ROOT)",
      "Two processing paths: full (mesh profiles per event) vs fast (CCE lookup table batch)",
    ]

key-files:
  created:
    - "src/mc_coupling.py"
    - "tests/test_mc_coupling.py"
  modified: []

key-decisions:
  - "column_map keys are standard names, values are source names -- consistent across CSV and ROOT loaders"
  - "uproot imported lazily inside ROOT functions so module works without uproot for CSV-only workflows"
  - "process_mc_ensemble filters zero-energy events and reports count rather than raising errors"

patterns-established:
  - "Standardized event DataFrame format: event_id, x_cm, y_cm, z_cm, edep_keV columns for all MC import sources"
  - "Lazy uproot import pattern for optional ROOT file support"
  - "Two-path processing: events_to_charge_profiles (mesh-resolved, slow) vs process_mc_ensemble (CCE lookup, fast)"

requirements-completed: [MCCP-01, MCCP-02, MCCP-03, MCCP-04]

duration: 4min
completed: 2026-03-31
---

# Plan 22-01: MC Coupling Module Summary

**CSV and ROOT event importers with unit conversion, per-event mesh charge profiles via ion_track_generation_2d, batch CCE(LET) lookup processing 1000 events in 0.34s, and log-spaced pulse height distribution -- all 23 tests pass**

## Performance

- **Duration:** ~4 min
- **Started:** 2026-03-31T16:03:11Z
- **Completed:** 2026-03-31T16:07:15Z
- **Tasks:** 2
- **Files created:** 2

## Accomplishments

- 7 public functions: convert_units, load_mc_events_csv, list_root_trees, load_mc_events_root, events_to_charge_profiles, process_mc_ensemble, pulse_height_distribution
- Full path (mesh-resolved charge profiles) and fast path (CCE lookup table) processing modes
- All 23 tests pass: 22 fast unit tests (0.68s) + 1 slow integration test (0.34s with 1000 events)
- Configurable column/branch mapping supports unknown Geant4 TTree naming conventions

## Task Commits

1. **Task 1: Create mc_coupling.py module** - `5b884c9` (feat)
2. **Task 2: Create test_mc_coupling.py** - `311fea0` (test)

## Files Created/Modified

- `src/mc_coupling.py` - MC event import (CSV + ROOT), unit conversion, mesh charge mapping, batch CCE lookup, pulse height distribution (514 lines)
- `tests/test_mc_coupling.py` - 23 tests: unit conversion (parametrized), CSV import, ensemble processing, PHD, charge profiles (mocked), ROOT import (mocked), integration pipeline (467 lines)

## Decisions Made

- column_map convention: keys are standard column names (event_id, x_cm, etc.), values are source column names in CSV/ROOT file. This is consistent across both loaders and matches the plan's specified convention.
- uproot imported lazily inside ROOT functions (list_root_trees, load_mc_events_root) so the module loads without uproot for CSV-only use cases.
- process_mc_ensemble filters zero-energy events silently (with logging) rather than raising errors, since zero-energy events are expected in MC output (particles that miss the SV).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed mock target for events_to_charge_profiles tests**

- **Found during:** Task 2 (test execution)
- **Issue:** `patch("src.mc_coupling.ion_track_generation_2d")` failed because the function is imported locally inside events_to_charge_profiles, not at module level
- **Fix:** Changed mock target to `patch("src.single_particle.ion_track_generation_2d")` which patches the actual function before local import resolves it
- **Files modified:** tests/test_mc_coupling.py
- **Verification:** Both charge profile tests pass
- **Committed in:** 311fea0

**2. [Rule 1 - Bug] Fixed ROOT test mocking strategy**

- **Found during:** Task 2 (test execution)
- **Issue:** `patch("src.mc_coupling.uproot")` failed because uproot is lazily imported, not a module-level attribute
- **Fix:** Used `patch.dict("sys.modules", {"uproot": mock_uproot})` to inject mock at the sys.modules level before lazy import
- **Files modified:** tests/test_mc_coupling.py
- **Verification:** Both ROOT tests pass with mocked uproot
- **Committed in:** 311fea0

---

**Total deviations:** 2 auto-fixed (both test mocking issues, Rule 1)
**Impact on plan:** Both fixes necessary for correct test mocking with lazy imports. No scope creep.

## Issues Encountered

None beyond the test mocking fixes documented above.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- mc_coupling.py ready for notebook visualization (Plan 22-02)
- All 7 functions importable and tested
- Full pipeline validated: synthetic events -> process_mc_ensemble -> pulse_height_distribution completes 1000 events in <1 second
- ROOT import tested via mocking; real ROOT file integration will need sample file from INFN-LNS group (existing blocker in STATE.md)

---

_Phase: 22-monte-carlo-coupling_
_Completed: 2026-03-31_
