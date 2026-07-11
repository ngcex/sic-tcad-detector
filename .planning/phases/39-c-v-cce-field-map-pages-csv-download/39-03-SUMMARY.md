---
phase: 39-c-v-cce-field-map-pages-csv-download
plan: 03
subsystem: ui
tags: [streamlit, plotly, appTest, devsim, csv-export]

# Dependency graph
requires:
  - phase: 39-01
    provides: proven petringa.run_cv / run_cce module-attribute mockability seam (monkeypatch.setattr under AppTest.from_function)
  - phase: 39-02
    provides: app/components/results.py (build_cv_figure, build_mott_schottky_figure, build_cce_figure, to_csv_bytes)
provides:
  - Working C-V page (app/workflows/cv.py): 1D-only pre-check guard, Run button, cached SimResult, C-V + Mott-Schottky Plotly charts, CSV download, optional depletion-width expander
  - Working CCE page (app/workflows/cce.py): 1D-only pre-check guard, Run button, cached SimResult, CCE-vs-bias Plotly chart, CSV download, I_collected/I_generated caption
  - tests/test_app_cv_page.py: AppTest coverage for happy path, 2D pre-check guard, empty-state guard
  - tests/test_app_cce_page.py: AppTest coverage for happy path, 2D pre-check guard, empty-state guard
affects: [phase-40-geometry-viewer]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "1D-only PRE-CHECK before calling a facade that raises NotImplementedError for 2D configs (never wrap the call in try/except; the pre-check keeps the failure path exception-free and testable via at.warning)"
    - "Run button writes SimResult into st.session_state[<page>_result]; a subsequent rerun (e.g. clicking Download) reads the cache instead of re-invoking the facade — no re-solve on rerun"

key-files:
  created: [tests/test_app_cv_page.py, tests/test_app_cce_page.py]
  modified: [app/workflows/cv.py, app/workflows/cce.py]

key-decisions:
  - "Mirrored home.py's docstring/empty-state-guard idiom verbatim across both pages for UI consistency"
  - "CCE page surfaces metadata[I_collected]/metadata[I_generated] as a small st.caption rather than a required plot, per RESEARCH Extra Plots guidance"
  - "C-V page implements the optional depletion-width-vs-V expander (RESEARCH A4 bonus), converting cm to um for display"

patterns-established:
  - "Facades referenced as petringa.run_cv / petringa.run_cce (module attribute), never `from petringa import run_cv`, to preserve monkeypatch mockability"

requirements-completed: [UI-03, UI-04, UI-06]

# Metrics
duration: 20min
completed: 2026-07-11
---

# Phase 39 Plan 03: C-V + CCE Pages Summary

**C-V and CCE pages both run their respective devsim facade on Run click, cache the SimResult in st.session_state, render Plotly charts from the cache, and offer a CSV download — each guarded by a 1D-only pre-check that fires before the facade is called, since run_cv/run_cce raise NotImplementedError for 2D configs.**

## Performance

- **Duration:** ~20 min
- **Completed:** 2026-07-11
- **Tasks:** 2 (TDD: test first, then implementation, per task)
- **Files modified:** 4

## Accomplishments

- Replaced the placeholder `app/workflows/cv.py` with a full Run -> cache -> render -> download page: C-V curve + Mott-Schottky Plotly figures, CSV download, optional depletion-width expander (UI-03)
- Replaced the placeholder `app/workflows/cce.py` with a full Run -> cache -> render -> download page: CCE-vs-bias Plotly figure, CSV download, I_collected/I_generated caption (UI-04)
- 1D-only guard implemented as a pre-check (`if cfg.half_width_um is not None: st.warning(...); st.stop()`) placed BEFORE `petringa.run_cv`/`petringa.run_cce` on both pages — no try/except around the facade call
- Added `tests/test_app_cv_page.py` and `tests/test_app_cce_page.py`, each with three AppTest cases (happy path with cache verification, 2D pre-check guard, empty-state guard), mocking the respective facade as a module attribute per the proven 39-01 seam
- Download half of UI-06 delivered for both pages via `to_csv_bytes`

## Task Commits

Each task was committed atomically:

1. **Task 1: C-V page + AppTest — 1D guard, Run→cache→render→download** - `d8fb564` (feat)
2. **Task 2: CCE page + AppTest — 1D guard, Run→cache→render→download** - `304424a` (feat)

**Plan metadata:** committed as part of this final documentation commit.

_Note: Both tasks followed TDD in spirit (test file written first and confirmed RED before writing the page to make it GREEN), but each task's test+implementation files were committed together in a single `feat` commit per plan convention (no separate `test(...)` commit), matching the pattern already established by the sibling 39-04 plan._

## Files Created/Modified

- `app/workflows/cv.py` - C-V page: empty-state guard, 1D-only pre-check, Run button, cached SimResult rendering (C-V + Mott-Schottky Plotly figures), CSV download, optional depletion-width expander
- `tests/test_app_cv_page.py` - AppTest suite mocking `petringa.run_cv`: run-caches-result, 2D-guard, empty-state-guard
- `app/workflows/cce.py` - CCE page: empty-state guard, 1D-only pre-check, Run button, cached SimResult rendering (CCE-vs-bias Plotly figure), CSV download, I_collected/I_generated caption
- `tests/test_app_cce_page.py` - AppTest suite mocking `petringa.run_cce`: run-caches-result, 2D-guard, empty-state-guard

## Decisions Made

- Mirrored the sibling 39-04 plan's page/test structure for consistency across all three results pages built in this phase
- CCE metadata surfaced as a caption (not a required plot), matching the plan's explicit guidance to keep I_collected/I_generated as supplementary text

## Deviations from Plan

None - plan executed exactly as written. Both tasks' `<action>` blocks specified the fake `SimResult` shapes, guard strings, and rendering order verbatim; all were implemented as specified.

## TDD Gate Compliance

Both tasks had `tdd="true"` and prescribed writing the AppTest first. For each task the test file was written first and confirmed to fail (RED — `IndexError: list index out of range` on `at.button[0]`, since the placeholder pages had no "Run simulation" button) before the page implementation was written to make it pass (GREEN). Only a single `feat(39-03): ...` commit exists per task in the git log (no separate `test(...)` commit) — the plan's task-level commit protocol treats this as one atomic task commit rather than a full RED/GREEN/REFACTOR multi-commit sequence, matching the sibling 39-04 plan's convention. This is a documentation note only; the plan's own acceptance criteria (verify commands, string/pattern checks) are all satisfied and were checked via grep/pytest before each commit.

## Issues Encountered

None. The RED failures on both tasks were the expected `IndexError` (no "Run simulation" button existed on the placeholder pages), confirming the test suites genuinely exercise new behavior rather than passing vacuously.

A concurrent sibling agent (plan 39-04) was executing against `app/workflows/field_map.py` and `tests/test_app_field_page.py` at the same time. Those files are disjoint from this plan's files, and no file conflicts occurred during either task's commit. STATE.md/ROADMAP.md/REQUIREMENTS.md were re-read fresh (after the sibling's own doc commits landed) before this plan's doc updates, per the concurrency instructions in the execution prompt.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- C-V and CCE pages are fully functional against mocked facades; ready for the 2D geometry viewer to be added in Phase 40 on these same pages
- No blockers. This plan's files (`app/workflows/cv.py`, `app/workflows/cce.py`, `tests/test_app_cv_page.py`, `tests/test_app_cce_page.py`) are disjoint from the sibling 39-04 plan's files (`app/workflows/field_map.py`, `tests/test_app_field_page.py`) — no merge conflicts encountered
- Full page + regression suites green: `uv run pytest tests/test_app_cv_page.py tests/test_app_cce_page.py -x` (6 passed), `uv run pytest tests/test_app_pages.py -x` (4 passed), `uv run python -c "import app.workflows.cv, app.workflows.cce"` exits 0
- Bare full-suite `uv run pytest tests/ -q` was not used as the acceptance gate — per STATE.md's documented devsim resource-exhaustion constraint (pre-existing, unrelated to this plan), per-file verification is the durable convention for this project

---

_Phase: 39-c-v-cce-field-map-pages-csv-download_
_Completed: 2026-07-11_

## Self-Check: PASSED

- FOUND: app/workflows/cv.py
- FOUND: app/workflows/cce.py
- FOUND: tests/test_app_cv_page.py
- FOUND: tests/test_app_cce_page.py
- FOUND: .planning/phases/39-c-v-cce-field-map-pages-csv-download/39-03-SUMMARY.md
- FOUND: commit d8fb564
- FOUND: commit 304424a
