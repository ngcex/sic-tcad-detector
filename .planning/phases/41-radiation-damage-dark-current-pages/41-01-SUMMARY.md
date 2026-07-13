---
phase: 41-radiation-damage-dark-current-pages
plan: 01
subsystem: ui
tags: [streamlit, plotly, devsim, parametric-sweep, csv-export]

# Dependency graph
requires:
  - phase: 37-core-api-cce-facades-parametric-sweep
    provides: petringa.ParametricSweep, petringa.run_dark_current, petringa.run_radiation_damage facades
provides:
  - build_damage_figure(result) -> go.Figure (CCE vs proton fluence, log-x, NaN-tolerant)
  - build_dark_current_figure(result) -> go.Figure (4-trace decomposition vs temperature, log-y, abs()+zero-guarded)
  - to_csv_bytes branches for sim_type="damage" and sim_type="dark_current"
  - Spike-confirmed safe widget defaults for radiation damage V_bias and dark current
    ParametricSweep(param="T") sim_kwargs, recorded in 41-01-SPIKE-NOTES.md
affects: [41-02-radiation-damage-page, 41-03-dark-current-page]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Pure single-SimResult-in / go.Figure-out builders extended to sim_type=damage/dark_current, matching the Phase 39/40 build_cv_figure/build_cce_figure/build_field_figures precedent exactly"
    - "to_csv_bytes elif-chain dispatch extended with two new branches following the exact existing df + extra_header_lines convention"
    - "Live-devsim spike-before-build: confirm facade/ParametricSweep convergence at chosen widget defaults BEFORE building UI pages around them, recorded in a phase-scoped SPIKE-NOTES.md for downstream plans to consume verbatim"

key-files:
  created:
    - .planning/phases/41-radiation-damage-dark-current-pages/41-01-SPIKE-NOTES.md
  modified:
    - app/components/results.py
    - tests/test_app_csv_export.py

key-decisions:
  - "ParametricSweep(param='T', sim_fn=run_dark_current) with sim_kwargs={v_start:-20.0, v_stop:-20.0, n_points:1} confirmed working on first attempt, no n_points=2 fallback needed"
  - "run_radiation_damage(V_bias=-20.0) still produces one NaN at fluence~=3.98e13, disconfirming RESEARCH.md Assumption A3 (shallower V_bias does NOT eliminate the NaN) — accepted via NaN-tolerant rendering per plan design, not chased with a different bias"
  - "test_unknown_sim_type_raises_value_error fixture changed from sim_type='damage' to sim_type='not_a_real_sim_type' to fix the Pitfall-4 regression now that 'damage' is a real to_csv_bytes branch"

patterns-established:
  - "Dark current decomposition trace order/colors fixed: Total (#1A1A1A), SRH (bulk) (#1F6FEB), TAT (effective) (#D32F2F), SRV (surface) (#2E7D32) — Phase 42 batch sweep page should reuse this palette per UI-SPEC"

requirements-completed: [FEAT-01, FEAT-02]

# Metrics
duration: 25min
completed: 2026-07-13
---

# Phase 41 Plan 01: Shared Foundation — Spike + Results Builders Summary

**Live-devsim spike confirms `ParametricSweep(param="T")` + `run_dark_current` tolerates a single-point bias request at `n_points=1`; `results.py` gains `build_damage_figure`/`build_dark_current_figure` pure Plotly builders and two new `to_csv_bytes` branches for Wave 2 to consume.**

## Performance

- **Duration:** 25 min
- **Started:** 2026-07-13T17:36:00Z
- **Completed:** 2026-07-13T18:01:33Z
- **Tasks:** 3
- **Files modified:** 2 (+ 1 new spike-notes doc)

## Accomplishments

- Confirmed via live devsim that `ParametricSweep(base_config=DeviceConfig(), param="T", values=[...], sim_fn=run_dark_current, sim_kwargs={"v_start": -20.0, "v_stop": -20.0, "n_points": 1})` returns a valid `list[SimResult]` of the expected length with no exception — the single largest unverified assumption in RESEARCH.md is now empirically confirmed
- Confirmed `run_radiation_damage(DeviceConfig(), V_bias=-20.0)` still produces a NaN at fluence≈3.98e13 (disconfirms an assumption in RESEARCH.md), and documented that this is an accepted, NaN-tolerant-rendering outcome rather than something Wave 2 should try to "fix" via a different bias
- Added `build_damage_figure` and `build_dark_current_figure` as pure Plotly builders to `app/components/results.py`, following the module's existing pure-function convention exactly (no `st.*` calls, importable without a Streamlit runtime)
- Extended `to_csv_bytes` with `damage` and `dark_current` branches following the module's exact existing `df` + `extra_header_lines` dispatch convention
- Fixed the Pitfall-4 regression in `tests/test_app_csv_export.py::test_unknown_sim_type_raises_value_error`, which had been using `sim_type="damage"` as its "unimplemented type" fixture — now uses `sim_type="not_a_real_sim_type"`

## Task Commits

1. **Task 1: Live-devsim spike — confirm ParametricSweep(param="T") + run_dark_current tolerates a minimal single-point bias request** - `4f2d6f7` (docs)
2. **Task 2: Extend app/components/results.py with build_damage_figure, build_dark_current_figure, and two new to_csv_bytes branches** - `3bd4ce0` (feat)
3. **Task 3: Add CSV export tests for damage/dark_current + fix the Pitfall-4 regression in test_unknown_sim_type_raises_value_error** - `d54b33c` (test)

**Plan metadata:** (this commit, following SUMMARY)

## Files Created/Modified

- `.planning/phases/41-radiation-damage-dark-current-pages/41-01-SPIKE-NOTES.md` - Confirmed safe widget defaults for both Wave 2 pages, sourced from live-devsim runs
- `app/components/results.py` - Added `build_damage_figure`, `build_dark_current_figure`, and two `to_csv_bytes` elif branches
- `tests/test_app_csv_export.py` - Added `test_damage_csv_columns_and_header`, `test_dark_current_csv_columns_and_header`, `test_build_dark_current_figure_guards_zero_and_negative`; fixed `test_unknown_sim_type_raises_value_error`

## Decisions Made

- Kept the plan's Task 1 contingency (`n_points=2` fallback) unused since `n_points=1` with `v_start == v_stop` worked cleanly on the first attempt — recorded explicitly in spike notes so Wave 2 doesn't re-test this.
- Did NOT attempt to find a `V_bias` value that eliminates the radiation-damage NaN (out of scope per the plan and confirmed via advisor consultation) — the NaN-tolerant rendering already specified in `build_damage_figure`'s design is the correct, already-planned handling.

## Deviations from Plan

None - plan executed exactly as written. The plan's own contingency path (Task 1's `n_points=2` fallback) was available but not needed.

## Issues Encountered

None. All three tasks' acceptance criteria passed on first verification without requiring fixes.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

Wave 2 (41-02 radiation damage page, 41-03 dark current page) can now proceed:

- `app/components/results.py::build_damage_figure` / `build_dark_current_figure` / extended `to_csv_bytes` are implemented, tested, and importable
- `.planning/phases/41-radiation-damage-dark-current-pages/41-01-SPIKE-NOTES.md` provides the confirmed-safe kwargs for both pages' facade/ParametricSweep calls — Wave 2 plans should read this file and use its values verbatim rather than re-deriving or re-testing
- No blockers. The one open physical caveat (radiation damage NaN at fluence≈3.98e13 even at `V_bias=-20.0`) is documented and has an already-designed UI handling path (NaN-tolerant Plotly rendering + UI-SPEC's `st.info` partial-failure banner), not a blocker for 41-02.

---

_Phase: 41-radiation-damage-dark-current-pages_
_Completed: 2026-07-13_
