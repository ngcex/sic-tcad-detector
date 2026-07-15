---
phase: 42-microdosimetry-page-batch-sweep-page
plan: 03
subsystem: ui
tags: [streamlit, parametric-sweep, apptest, plotly, batch-sweep]

# Dependency graph
requires:
  - phase: 42-microdosimetry-page-batch-sweep-page (plan 01)
    provides: build_sweep_overlay_figure (4-arg), sweep_results_to_csv_bytes, convergence-safe default combo (run_cce + epi_thickness_um=[10,15,20])
provides:
  - Batch Sweep Streamlit page (app/workflows/batch_sweep.py) driving the real etna.ParametricSweep(...).run()
  - Curated SWEEPABLE_FIELDS (8 numeric 1D-safe fields) and SIM_FACADES (run_cce/run_cv/run_temperature_sweep) selectboxes
  - Value-list float()-per-token parse (no code eval); partial-failure-tolerant overlay; bulk "Download all results as CSV" export
  - AppTest suite (6 cases) monkeypatching the facade only, exercising real ParametricSweep orchestration
affects:
  [
    phase-42 microdosimetry page (parallel plan 42-02),
    future sweep-driven pages,
  ]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Curated selectbox constrains ParametricSweep.param to vetted numeric DeviceConfig fields (no free-text attribute injection)"
    - "Widget-key vs run-snapshot-key disjoint namespaces (sweep_param/sweep_values widgets vs sweep_run_param/sweep_run_values snapshot) to avoid StreamlitAPIException on streamlit 1.58"
    - "Per-swept-value skip-empty aggregation + n_ok/n_requested partial-failure banner (generalized from dark_current.py's per-temperature pattern)"

key-files:
  created:
    - tests/test_app_batch_sweep_page.py
  modified:
    - app/workflows/batch_sweep.py

key-decisions:
  - "Referenced etna.ParametricSweep and facades as module attributes via getattr so tests monkeypatch the facade (run_cce) while the real .run() executes"
  - "float() per token inside try/except ValueError for the value list; never eval/exec (T-42-03-T2 mitigation)"
  - "Hardcoded the 42-01-SPIKE-confirmed convergence-safe defaults: run_cce + epi_thickness_um + '10, 15, 20'"

patterns-established:
  - "Batch sweep page = general-case ParametricSweep page; dark_current.py is its hard-coded param=T special case"

requirements-completed: [FEAT-04]

# Metrics
duration: 12 min
completed: 2026-07-14
---

# Phase 42 Plan 03: Batch Sweep Page Summary

**General-case parametric-sweep Streamlit page: curated param/facade selectboxes drive the real etna.ParametricSweep(...).run(), rendering one overlay trace per swept value with a partial-failure-tolerant chart and a bulk "Download all results as CSV" export.**

## Performance

- **Duration:** 12 min
- **Started:** 2026-07-14T15:24:00Z
- **Completed:** 2026-07-14T15:37:00Z
- **Tasks:** 2
- **Files modified:** 2 (1 created, 1 modified)

## Accomplishments

- Replaced the batch_sweep.py placeholder with a full parametric-sweep page driving the real `ParametricSweep.run()` (FEAT-04)
- Closed both security surfaces: curated `SWEEPABLE_FIELDS`/`SIM_FACADES` selectboxes (no attribute injection) + `float()`-per-token value parsing (no code eval)
- Partial per-value convergence failures degrade to a partial overlay + warning; all-fail shows an error; RuntimeError from the sweep shows an error — never crashes
- Bulk "Download all results as CSV" export via `sweep_results_to_csv_bytes`
- Six AppTest cases pass, exercising genuine sweep wiring (facade monkeypatched, `ParametricSweep` never monkeypatched)

## Task Commits

1. **Task 1: Implement the batch sweep page** - `b41f7a2` (feat)
2. **Task 2: AppTest coverage for the batch sweep page** - `2e75c46` (test)

**Plan metadata:** _(this SUMMARY commit)_

_Note: Task 2 is a `tdd="true"` task, but Task 1 delivered the implementation first by design, so the tests went green immediately — no manufactured RED phase._

## Files Created/Modified

- `app/workflows/batch_sweep.py` - Batch Sweep page: empty-state + 1D guards, curated selectboxes, value-list parse, real `ParametricSweep(...).run()`, per-value skip-empty aggregation, partial-failure banner, overlay chart, bulk CSV download
- `tests/test_app_batch_sweep_page.py` - Six AppTest cases: empty-state, 1D-only guard, real-sweep-caches-≥3-results, bad-value-list error, partial-value-failure warning, facade RuntimeError error

## Decisions Made

- Referenced `etna.ParametricSweep` and facades as module attributes (`getattr(etna, ...)`) so tests monkeypatch the facade (`run_cce`) while the real `.run()` orchestration executes — asserts genuine sweep wiring, not a hand-rolled loop.
- Value list parsed with `float()` per token inside `try/except ValueError`; no `eval`/`exec` anywhere in source (T-42-03-T2).
- Run snapshot written to renamed keys `sweep_run_param`/`sweep_run_values` (disjoint from the `sweep_param`/`sweep_values` widget keys) to avoid the streamlit 1.58 `StreamlitAPIException`.
- Hardcoded the 42-01-SPIKE-confirmed default combo (`run_cce` + `epi_thickness_um` + `"10, 15, 20"`).

## Deviations from Plan

None to the production code — plan executed exactly as written for both tasks. One environment blocker was auto-resolved (below), which did not alter any plan-specified code.

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Worktree venv missing dev deps + stale base — pytest/streamlit resolved from miniconda without plotly**

- **Found during:** Task 2 (running the AppTest suite)
- **Issue:** Two coupled blockers in the freshly-created worktree. (a) The worktree branch was created from an old Phase-39-era base (`5e0b686`) via the known Claude Code worktree bug (#2015), so the Wave 1 (42-01) foundation — `build_sweep_overlay_figure`, `sweep_results_to_csv_bytes` — was absent and every import would fail. (b) The worktree `.venv` had only the base (production) deps, not the `dev` optional group, so `uv run pytest` fell back to miniconda's pytest, which ran the AppTest ScriptRunner against miniconda's streamlit where `plotly` is not installed → `ModuleNotFoundError: No module named 'plotly'` (the identical pre-existing `test_app_dark_current_page.py` failed the same way in this worktree but passes in the main repo).
- **Fix:** (a) Per the Pattern A worktree startup contract, hard-reset the per-agent branch to main's Wave-1 tip (`git reset --hard 1a3051b`) — the worktree branch had zero unique commits (merge-base == its own HEAD), so no work was lost; this brought in the 42-01 foundation. (b) Ran the suite via `uv run --extra dev pytest`, which installs the already-declared `pytest>=7.0` dev dependency into the worktree venv so streamlit/plotly resolve from the venv. No new or unlisted package was installed.
- **Files modified:** None (environment/venv only; no source changes)
- **Verification:** `uv run --extra dev pytest tests/test_app_batch_sweep_page.py -q` → 6 passed; combined run with `test_app_dark_current_page.py` + `test_app_csv_export.py` → 20 passed (no regression). NOTE for the verifier/orchestrator: the FIRST test invocation in this worktree MUST include `--extra dev` (e.g. `uv run --extra dev pytest ...`) — a bare `uv run pytest` re-syncs the venv to the base (production-only) group and falls back to the system (miniconda) interpreter where `plotly` is absent. Once `--extra dev` has run once, the dev deps persist and the plan's literal bare command (`uv run pytest tests/test_app_batch_sweep_page.py -q`) also passes (confirmed → 6 passed). In the main repo (dev deps already installed) the bare command passes directly.
- **Committed in:** N/A (no code change — reset predates the task commits; venv state is not tracked)

---

**Total deviations:** 1 auto-fixed (1 blocking, environment-only).
**Impact on plan:** No scope creep and no production-code change. The reset restored the required Wave 1 base; the `--extra dev` invocation is the correct way to run the declared test dependency in this venv.

## Issues Encountered

None during planned work beyond the environment blocker documented above.

## Threat Surface Notes

No new security-relevant surface beyond the plan's threat register. The two mitigations the plan mandated are implemented and verified: curated selectbox (`SWEEPABLE_FIELDS`) keeps `ParametricSweep.param` a vetted field (T-42-03-T1), and `float()`-per-token parsing keeps the value list off any code evaluator (T-42-03-T2). No `eval`/`exec` in source; `run_field` excluded from `SIM_FACADES` (T-42-03-D2).

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Batch Sweep page complete and tested; FEAT-04 satisfied.
- Wave 2 sibling (plan 42-02, microdosimetry page) runs in a separate worktree touching disjoint files — no conflict expected on merge.
- Note for the orchestrator: this worktree was reset to `1a3051b` (main's Wave-1 tip) at startup to recover from the #2015 old-base bug; the merge should be clean.

## Self-Check: PASSED

- `app/workflows/batch_sweep.py` exists (modified) — FOUND
- `tests/test_app_batch_sweep_page.py` exists (created) — FOUND
- Commit `b41f7a2` (feat 42-03) — FOUND
- Commit `2e75c46` (test 42-03) — FOUND
- `uv run --extra dev pytest tests/test_app_batch_sweep_page.py -q` — 6 passed

---

_Phase: 42-microdosimetry-page-batch-sweep-page_
_Completed: 2026-07-14_
