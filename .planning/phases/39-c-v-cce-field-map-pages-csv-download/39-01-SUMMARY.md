---
phase: 39-c-v-cce-field-map-pages-csv-download
plan: 01
subsystem: testing
tags: [plotly, streamlit, apptest, monkeypatch, uv]

# Dependency graph
requires:
  - phase: 38
    provides: app/workflows/ page structure, AppTest.from_function pattern (test_empty_state_guard), streamlit session_state conventions
provides:
  - plotly materialized into the project venv (was declared but not installed)
  - Empirical proof that etna.run_* module-attribute references are mockable via monkeypatch under AppTest.from_function
  - Mandated page import structure for all Phase 39 page plans: `import etna; etna.run_*(cfg)`
affects: [39-02, 39-03, 39-04, "any future phase adding Streamlit pages that call etna run_* facades"]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "etna.run_* facades must be referenced as module attributes (import etna; etna.run_cv(cfg)) in page code, never `from etna import run_cv`, so page tests can monkeypatch.setattr(etna, 'run_cv', fake) to avoid real devsim solves"

key-files:
  created: [tests/test_app_run_mockability.py]
  modified: []

key-decisions:
  - "Confirmed RESEARCH A6 empirically: etna.run_cv as a module attribute is interceptable by monkeypatch under AppTest.from_function, because AppTest.from_function runs the wrapper body in the same process as the test, so the monkeypatch-patched module object is the same object the wrapper's `import etna` resolves to."

patterns-established:
  - "Spike-test pattern for facade mockability: define a fake_run_* returning a hand-built SimResult, monkeypatch.setattr(etna, 'run_*', fake), call the facade inside a self-contained AppTest.from_function wrapper body, assert at.exception == [] and inspect at.session_state for the stashed result."

requirements-completed: []  # NOTE: plan frontmatter lists [UI-03, UI-04, UI-05], but this
  # plan's tasks/acceptance-criteria are entirely plotly materialization + mockability
  # spike (no page code shipped). UI-03/04/05 are genuinely delivered by 39-02/39-03/39-04,
  # which also carry those IDs in their own frontmatter. Marking them complete here would
  # be false — see Deviations section below.

# Metrics
duration: 12min
completed: 2026-07-11
---

# Phase 39 Plan 01: Wave 0 Spike — plotly materialization + run_cv mockability Summary

**Materialized the already-declared plotly>=5.0 dependency via `uv sync` and empirically proved that `etna.run_cv` referenced as a module attribute (`etna.run_cv(cfg)`) is interceptable by `monkeypatch.setattr` under `AppTest.from_function`, dictating the mandatory page import structure for all Phase 39 page plans.**

## Performance

- **Duration:** 12 min
- **Started:** 2026-07-11T11:39:00Z (approx)
- **Completed:** 2026-07-11T11:51:09Z
- **Tasks:** 2
- **Files modified:** 1 created

## Accomplishments

- `plotly` and `plotly.graph_objects` now import cleanly in the project venv (`uv sync` materialized the pre-declared `plotly>=5.0` pin) — unblocks every downstream AppTest and page that imports plotly, per RESEARCH Pitfall 2
- Empirically confirmed RESEARCH A6: `etna.run_cv` used as a module attribute (`import etna; etna.run_cv(cfg)`) is mockable via `monkeypatch.setattr(etna, "run_cv", fake_run_cv)` inside an `AppTest.from_function` wrapper — the fake (2-point array) is proven invoked instead of a real devsim solve (which defaults to 40 points)
- Decided and documented the mandatory page import style for 39-03/39-04: `import etna; etna.run_*(cfg)`, never `from etna import run_*`

## Task Commits

Each task was committed atomically:

1. **Task 1: Materialize plotly into the venv via uv** — no source files modified (dependency materialization only); verified inline, no separate commit (folded into Task 2's commit per plan's `files_modified` scope, which lists only the test file)
2. **Task 2: Spike — prove etna.run_cv is mockable under AppTest.from_function** - `c8d1c7d` (test)

**Plan metadata:** (this commit, following SUMMARY write)

_Note: Task 1 has no `files_modified` entry in the plan frontmatter (dependency sync only) — its verification is folded into the Task 2 commit message._

## Files Created/Modified

- `tests/test_app_run_mockability.py` - Spike test proving `etna.run_cv` is mockable as a module attribute under `AppTest.from_function`; establishes the monkeypatch seam pattern for all future page tests

## Decisions Made

- Confirmed RESEARCH A6 (module-attribute reference style for `etna.run_*` facades) via a passing empirical spike, rather than accepting it purely on paper. This is now a hard constraint for 39-03/39-04 page implementations.
- Task 1 (uv sync) produced no diff to commit on its own since `pyproject.toml` was correctly left untouched (dependency was already pinned) and `uv.lock`/venv state is not itself a tracked deliverable file per the plan's `files_modified` scope — verified via `uv run python -c "import plotly.graph_objects"` exiting 0 and `git diff --stat pyproject.toml` being empty.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] Corrected plan-frontmatter/requirements mismatch — did not mark UI-03/04/05 complete**

- **Found during:** Post-execution state-update step (requirements.mark-complete)
- **Issue:** This plan's frontmatter lists `requirements: [UI-03, UI-04, UI-05]`, but the plan's own tasks and acceptance criteria are entirely about materializing `plotly` and proving a monkeypatch mock seam — no C-V/CCE/field-map page code was written or shipped by this plan. Running `requirements.mark-complete` mechanically against the frontmatter would have checked off UI-03/UI-04/UI-05 as `[x]` Complete in `REQUIREMENTS.md` (both the checkbox list and the traceability table), which is factually false: those user-facing requirements ("User can click Run simulation on the C-V/CCE/field page...") are genuinely delivered by later plans 39-02/39-03/39-04, which independently carry the same requirement IDs in their own frontmatter (`39-02: [UI-03, UI-04, UI-05, UI-06]`, `39-03: [UI-03, UI-04, UI-06]`, `39-04: [UI-05, UI-06]`).
- **Fix:** Reverted `REQUIREMENTS.md` to its pre-mark-complete state (`git diff` on the file is empty) — UI-03/04/05 remain `[ ]` Pending in both the checklist and traceability table, to be marked complete only when 39-02/39-03/39-04 actually ship the pages. Set this SUMMARY's `requirements-completed` frontmatter field to `[]` instead of copying the plan's `requirements` field verbatim, with an inline note explaining why.
- **Files modified:** `.planning/REQUIREMENTS.md` (reverted), `.planning/phases/39-c-v-cce-field-map-pages-csv-download/39-01-SUMMARY.md` (requirements-completed field)
- **Verification:** `git diff .planning/REQUIREMENTS.md` is empty (file byte-identical to pre-execution state)
- **Committed in:** plan metadata commit (this SUMMARY's commit)

---

**Total deviations:** 1 auto-fixed (1 missing-critical correction to state integrity)
**Impact on plan:** No code/test impact. This is a planning-doc quality issue flagged for the next planner: 39-01's frontmatter `requirements` field should not have listed UI-03/04/05 since this plan is a pure Wave-0 prerequisite spike with no user-facing deliverable — those IDs correctly belong only to 39-02/39-03/39-04.

## Issues Encountered

`uv sync` uninstalled 84 packages during materialization (jupyter, pytest, and their transitive deps) that were present in the environment's Python but not declared in `pyproject.toml`'s tracked dependency groups. This did not break anything: `uv run pytest` and `uv run python -c "import plotly"` both continued to work correctly afterward (uv manages an ephemeral/synced run environment per invocation), and both plan-level verification commands pass. No action was needed — flagged here only as an observed side effect of running `uv sync`, not a defect requiring a fix.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- plotly is confirmed importable in the venv; downstream page plans (39-02, 39-03, 39-04) can safely import `plotly.graph_objects` in AppTest-covered code without import failures.
- The mock seam for `etna.run_*` facades is proven and documented: all Phase 39 page implementations MUST use `import etna; etna.run_*(cfg)` (module-attribute call style), never `from etna import run_*`, to remain testable via `monkeypatch.setattr(etna, "run_*", fake)`.
- No blockers for Wave 1+ plans.

---

_Phase: 39-c-v-cce-field-map-pages-csv-download_
_Completed: 2026-07-11_

## Self-Check: PASSED

- FOUND: `tests/test_app_run_mockability.py`
- FOUND commit: `c8d1c7d`
- Re-ran all task-level acceptance criteria: all PASS
- Re-ran plan-level `<verification>`: `uv run python -c "import plotly.graph_objects"` exit 0; `uv run pytest tests/test_app_run_mockability.py -x` exit 0, 1 passed
