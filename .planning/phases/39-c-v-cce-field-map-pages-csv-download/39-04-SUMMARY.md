---
phase: 39-c-v-cce-field-map-pages-csv-download
plan: 04
subsystem: ui
tags: [streamlit, plotly, appTest, devsim, csv-export]

# Dependency graph
requires:
  - phase: 39-01
    provides: proven petringa.run_field module-attribute mockability seam (monkeypatch.setattr under AppTest.from_function)
  - phase: 39-02
    provides: app/components/results.py (build_field_figures, to_csv_bytes)
provides:
  - Working Field Map page (app/workflows/field_map.py): 1D pre-check guard, Run button, cached SimResult, dual Plotly charts (E-field + potential), CSV download, optional net-doping expander
  - tests/test_app_field_page.py: AppTest coverage for happy path, 2D pre-check guard, empty-state guard
affects: [phase-40-geometry-viewer]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "1D-only PRE-CHECK before calling a facade that silently returns empty arrays for unsupported configs (never wrap in try/except when the failure mode is silent, not an exception)"
    - "Single run_field() call feeds two independent Plotly figures (E-field, potential) via build_field_figures returning a 2-tuple"

key-files:
  created: [tests/test_app_field_page.py]
  modified: [app/workflows/field_map.py]

key-decisions:
  - "Mirrored the exact page/test structure already established by the sibling 39-03 plan (app/workflows/cv.py, tests/test_app_cv_page.py) for consistency across the three results pages"
  - "Added optional st.expander('Net doping vs depth') per RESEARCH A4 bonus, plotting metadata['net_doping'] on a log y-axis, opt-in and not covered by acceptance criteria"

patterns-established:
  - "Facade referenced as petringa.run_field (module attribute), never `from petringa import run_field`, to preserve monkeypatch mockability"

requirements-completed: [UI-05, UI-06]

# Metrics
duration: 15min
completed: 2026-07-11
---

# Phase 39 Plan 04: Field Map Page Summary

**Field map page runs run_field on click, caches the SimResult, renders both E-field-vs-depth and potential-vs-depth Plotly charts from one solve, and offers CSV download — guarded by a 1D-only pre-check that fires before run_field is ever called (since run_field silently returns empty arrays for 2D instead of raising).**

## Performance

- **Duration:** ~15 min
- **Completed:** 2026-07-11
- **Tasks:** 1 (TDD: test first, then implementation)
- **Files modified:** 2

## Accomplishments

- Replaced the placeholder `app/workflows/field_map.py` with a full Run -> cache -> render -> download page
- 1D-only guard implemented as a pre-check (`if cfg.half_width_um is not None: st.warning(...); st.stop()`) placed BEFORE any call to `petringa.run_field`, correctly handling the documented trap that `run_field` does not raise for 2D configs
- Both required UI-05 charts (E-field vs depth, potential vs depth) rendered from a single cached `run_field` result via `build_field_figures`
- CSV download wired via `to_csv_bytes`, satisfying the download half of UI-06 for this page
- Added `tests/test_app_field_page.py` with three AppTest cases (happy path with cache verification, 2D pre-check guard, empty-state guard), all mocking `petringa.run_field` as a module attribute per the proven 39-01 seam

## Task Commits

Each task was committed atomically:

1. **Task 1: Field map page + AppTest — 1D pre-check, Run→cache→render(E-field + potential)→download** - `1051b41` (feat)

**Plan metadata:** committed as part of this final documentation commit.

_Note: This task followed TDD in spirit (test file written and run to confirm RED before writing the page to make it GREEN), but both files were committed together in a single `feat` commit per plan convention (no separate `test(...)` commit was required by the plan) — see TDD Gate Compliance note below._

## Files Created/Modified

- `app/workflows/field_map.py` - Field map page: empty-state guard, 1D-only pre-check, Run button, cached SimResult rendering, dual Plotly charts, CSV download, optional net-doping expander
- `tests/test_app_field_page.py` - AppTest suite mocking `petringa.run_field`: run-caches-result, 2D-guard, empty-state-guard

## Decisions Made

- Mirrored the sibling plan's (39-03, `cv.py`) exact page and test structure for consistency across all three results pages built in this phase
- Implemented the optional net-doping-vs-depth expander (RESEARCH A4) since it was explicitly called out in the plan's action steps, even though not required by acceptance criteria

## Deviations from Plan

None - plan executed exactly as written. The task's `<action>` block explicitly asked for the fake `SimResult` shape, guard strings, and rendering order used verbatim; all were implemented as specified.

## TDD Gate Compliance

The plan's task had `tdd="true"` and prescribed writing the AppTest first. The test file was written first and confirmed to fail (RED — `IndexError: list index out of range` on `at.button[0]`, since the placeholder page had no button) before the page implementation was written to make it pass (GREEN). Both the RED confirmation and GREEN implementation were done in this session, but only a single `feat(39-04): ...` commit was created (no separate `test(...)` commit exists in git log for this plan) — the plan's task-level commit protocol treats this as one atomic task commit rather than a full RED/GREEN/REFACTOR multi-commit sequence, consistent with how the sibling 39-03 plan's commits are structured. This is a documentation note only; the plan's own acceptance criteria (verify commands, string/pattern checks) are all satisfied.

## Issues Encountered

None. The RED failure was the expected `IndexError` (no "Run simulation" button existed on the placeholder page), confirming the test suite genuinely exercises new behavior rather than passing vacuously.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Field map page is fully functional against the mocked facade; ready for the 2D geometry viewer to be added in Phase 40 on this same page (per the docstring note)
- No blockers. This plan's files (`app/workflows/field_map.py`, `tests/test_app_field_page.py`) are disjoint from the sibling 39-03 plan's files (`app/workflows/cv.py`, `app/workflows/cce.py`, `tests/test_app_cv_page.py`, `tests/test_app_cce_page.py`) — no merge conflicts encountered
- Full page + regression suites green: `uv run pytest tests/test_app_field_page.py -x` (3 passed), `uv run pytest tests/test_app_pages.py -x` (4 passed)

---

_Phase: 39-c-v-cce-field-map-pages-csv-download_
_Completed: 2026-07-11_

## Self-Check: PASSED

- FOUND: app/workflows/field_map.py
- FOUND: tests/test_app_field_page.py
- FOUND: .planning/phases/39-c-v-cce-field-map-pages-csv-download/39-04-SUMMARY.md
- FOUND: commit 1051b41
