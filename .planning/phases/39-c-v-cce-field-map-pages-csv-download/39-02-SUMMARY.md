---
phase: 39-c-v-cce-field-map-pages-csv-download
plan: 02
subsystem: ui
tags: [plotly, pandas, streamlit, csv-export, results-rendering]

# Dependency graph
requires:
  - phase: 39-01
    provides: plotly materialized in venv via uv sync; proved petringa.run_* module-attribute mockability under AppTest
provides:
  - app/components/results.py — five pure Plotly go.Figure builders (build_cv_figure, build_mott_schottky_figure, build_cce_figure, build_field_figures) plus to_csv_bytes(result) CSV serializer, all st.*-free
  - Confirmed shape ground truth: I_collected is a bias-aligned numpy array, I_generated is a scalar total generated current (read directly from petringa/core/charge_collection.py cce_vs_bias return block)
affects: [39-03, 39-04]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Pure-seam module pattern (mirrors app/components/device_sidebar.py): all chart/CSV logic lives in st.*-free functions so it is unit-testable without a Streamlit runtime or devsim"
    - "CSV metadata header: leading '# '-commented lines (software_version, generated timestamp, full asdict(config) device provenance) prepended to a pandas.to_csv(index=False) table body, re-readable via pandas.read_csv(comment='#')"

key-files:
  created:
    - app/components/results.py
    - tests/test_app_csv_export.py
  modified: []

key-decisions:
  - "CCE metadata: I_collected is a bias-aligned array -> becomes the I_collected_A_per_cm2 column; I_generated is a scalar total generated current -> moved to a '# I_generated_A_per_cm2:' header line instead of a misleading broadcast column. Confirmed by reading cce_vs_bias's return block in petringa/core/charge_collection.py (I_collected=I_sorted array, I_generated=Q*np.trapezoid(...) scalar) before finalizing the column split, per plan instruction — matches RESEARCH's [ASSUMED] A4 exactly, no deviation needed."
  - "uv sync only materializes pyproject.toml's base [project.dependencies] group; pytest lives in [project.optional-dependencies].dev and was not present in .venv (uv run pytest silently fell back to a system/conda pytest lacking the project's installed packages, causing a plotly ModuleNotFoundError on collection). Fixed by uv sync --extra dev to materialize pytest into the project's own venv."

patterns-established:
  - "Pattern 1: Shared results.py module built once — all three Phase 39 result pages (C-V, CCE, field-map) import the same builders/serializer rather than duplicating chart or CSV logic per page."

requirements-completed: [] # 39-02 builds the shared substrate consumed by UI-03..UI-06; those requirements describe end-to-end page behavior (a user clicking "Run simulation" and seeing plots, or downloading a CSV from a page) and are only satisfied once 39-03/39-04 build the actual result pages. Frontmatter listed [UI-03, UI-04, UI-05, UI-06] as "contributes to," not "completes."

# Metrics
duration: 25min
completed: 2026-07-11
---

# Phase 39 Plan 02: Shared results.py (Plotly builders + CSV serializer) Summary

**Built app/components/results.py with five pure Plotly go.Figure builders and a to_csv_bytes(result) serializer with a commented metadata header, dispatching on sim_type across cv/cce/field — fully covered by a pure unit test with no Streamlit or devsim dependency.**

## Performance

- **Duration:** 25 min
- **Started:** 2026-07-11T00:00:00Z (approx)
- **Completed:** 2026-07-11 (see git log)
- **Tasks:** 2 completed
- **Files modified:** 2 created (app/components/results.py, tests/test_app_csv_export.py)

## Accomplishments

- Five pure `go.Figure` builders (`build_cv_figure`, `build_mott_schottky_figure`, `build_cce_figure`, `build_field_figures`) that transform a `SimResult` into Plotly figures with no `st.*` calls, mirroring `petringa/core/plotting.py`'s matplotlib titles/axis labels/reference lines exactly (C-V, Mott-Schottky, CCE-vs-bias with y=1.0 reference line, E-field/potential vs depth).
- `to_csv_bytes(result: SimResult) -> bytes` dispatches on `result.sim_type` to produce exact-schema CSV columns per type (cv/cce/field) with a leading `#`-commented metadata header (software version, ISO-8601 UTC timestamp, full `asdict(config)` device provenance, and a CCE-only `I_generated_A_per_cm2` header line).
- Confirmed via direct source read (not devsim execution) that `I_collected` from `cce_vs_bias` is a bias-aligned numpy array while `I_generated` is a scalar — validated the plan's pre-adopted assumption A4 and avoided a misleading broadcast column.
- `tests/test_app_csv_export.py`: 4 pure unit tests (cv, cce, field, unknown-sim_type ValueError) — all pass without Streamlit or devsim.

## Task Commits

Each task was committed atomically:

1. **Task 1: Plotly figure builders in results.py** - `aeb2126` (feat)
2. **Task 2: to_csv_bytes serializer + pure unit test** - `53236ef` (test)

**Plan metadata:** (this commit)

_Note: Task 1's commit also included the `to_csv_bytes` function body since both were authored in a single file write before task-by-task verification began — see Deviations below._

## Files Created/Modified

- `app/components/results.py` - Five pure Plotly builders + `to_csv_bytes` CSV serializer (no `st.*` calls; imports cleanly without a Streamlit runtime)
- `tests/test_app_csv_export.py` - Pure unit test covering all three `sim_type` CSV schemas, metadata header contents, and the unknown-sim_type `ValueError` path

## Decisions Made

- CCE column split (I_collected as array column, I_generated as scalar header line) confirmed by reading `petringa/core/charge_collection.py`'s `cce_vs_bias` return block directly, per the plan's mandatory read-the-source-first instruction — no devsim run was needed or performed.
- `uv sync --extra dev` was run to materialize `pytest` (declared under `[project.optional-dependencies].dev` in `pyproject.toml`, not the base dependency group) into the project's own `.venv`, fixing a `uv run pytest` fallback to a system pytest that could not see `plotly` in the project venv, causing a collection-time `ModuleNotFoundError`. No `pyproject.toml` changes were needed — the dependency was already correctly declared; only the venv sync state was stale (consistent with the same observation logged in 39-01-SUMMARY.md).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Synced the `dev` dependency group so `pytest` runs from the project venv**

- **Found during:** Task 2 (running the plan's verification command)
- **Issue:** `uv run pytest tests/test_app_csv_export.py -x` collected via a system/conda `pytest` (found on `PATH`, not in `.venv/bin`) that could not import `plotly` from the project's own venv, failing test collection with `ModuleNotFoundError: No module named 'plotly'`. `pytest` was declared in `pyproject.toml`'s `[project.optional-dependencies].dev` group but had never been synced into `.venv` (base `uv sync` only materializes `[project.dependencies]`).
- **Fix:** Ran `uv sync --extra dev`, which installed `pytest==9.1.1` (plus `jupyter` and transitive deps) into `.venv`. No `pyproject.toml` edits were needed — the dependency pin was already correct.
- **Files modified:** None tracked (venv-only change, `pyproject.toml`/`uv.lock` diff is empty — same non-tracked-deliverable pattern noted in 39-01-SUMMARY.md).
- **Verification:** `uv run pytest tests/test_app_csv_export.py -x -v` now runs from `/Users/ngcex/projects/physics/petringa/.venv/bin/python3` and all 4 tests pass.
- **Committed in:** N/A (no file change to commit; environment-only fix)

**2. [Process note, not a Rule 1-4 deviation] Task 1 and Task 2's edits to `results.py` landed in one file-write, split across two commits**

- **Found during:** Staging Task 2's commit
- **Issue:** The entire `app/components/results.py` (builders + `to_csv_bytes`) was authored in a single `Write` call before running Task 1's acceptance-criteria gate. Task 1's commit (`aeb2126`) therefore already contains the `to_csv_bytes` function body that Task 2's action described adding.
- **Fix:** No code fix needed — content is correct and matches both tasks' specs exactly. Task 2's commit (`53236ef`) captures the new test file, which is the file genuinely unique to Task 2's `<files>` list.
- **Files modified:** None (documentation-only note).
- **Verification:** `git diff --stat HEAD` before Task 2's commit was empty for `results.py`, confirming no re-edit occurred; all Task 1 and Task 2 acceptance criteria independently verified pass against the final committed state.
- **Committed in:** N/A (process observation only)

**3. [Documentation correction, not a Rule 1-4 code deviation] Did not mark UI-03/UI-04/UI-05/UI-06 complete in REQUIREMENTS.md**

- **Found during:** Close-out (requirements-tracking step)
- **Issue:** The plan's frontmatter lists `requirements: [UI-03, UI-04, UI-05, UI-06]`, and the standard close-out step marks all listed requirement IDs complete verbatim. Running `requirements.mark-complete` on all four flipped them to `[x]` Complete in `.planning/REQUIREMENTS.md`. But the requirement text describes end-to-end page behavior ("User can click 'Run simulation' on the C-V page... and see an interactive Plotly... plot", "User can download simulation results as a CSV file from any result page") that does not exist after this plan — 39-02 only builds the shared `app/components/results.py` substrate; no Streamlit page imports or calls it yet. 39-03 (frontmatter `[UI-03, UI-04, UI-06]`) and 39-04 (frontmatter `[UI-05, UI-06]`) build the actual pages. Notably, 39-01 also listed `[UI-03, UI-04, UI-05]` in its own frontmatter and correctly did NOT mark them complete (`mark-complete` returned `already_complete: []` for 39-01's run) — this phase's own precedent is "don't mark until delivered."
- **Fix:** Reverted `.planning/REQUIREMENTS.md` — both the `- [x]`/`- [ ]` checkboxes (lines ~30-33) and the traceability table rows (UI-03..UI-06) — back to `Pending`/`[ ]`. Set this Summary's `requirements-completed:` frontmatter to `[]` with an explanatory comment, and added a "Requirements note" under Next Phase Readiness so 39-03/39-04's executors mark these complete themselves once the actual pages exist.
- **Files modified:** `.planning/REQUIREMENTS.md`, this SUMMARY's frontmatter.
- **Verification:** `grep -n "UI-0[3-6]" .planning/REQUIREMENTS.md` shows all four as `[ ]` / `Pending` after the revert.
- **Committed in:** Plan metadata commit (this close-out), alongside SUMMARY.md/STATE.md/ROADMAP.md/REQUIREMENTS.md.

---

**Total deviations:** 1 auto-fixed (1 blocking - environment/venv sync), 2 process notes (no code impact — one commit-split note, one requirements-tracking correction).
**Impact on plan:** The dev-dependency sync was necessary to run the plan's own verification command correctly and is fully reversible/re-derivable from `pyproject.toml` (no drift risk). The commit-split note has zero functional impact. The requirements correction prevents REQUIREMENTS.md from asserting user-facing functionality that does not exist yet — necessary for traceability accuracy, not scope creep.

## Issues Encountered

None beyond the venv/pytest sync documented above.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- `app/components/results.py` exports the exact five function signatures (`build_cv_figure`, `build_mott_schottky_figure`, `build_cce_figure`, `build_field_figures`, `to_csv_bytes`) that plans 39-03 and 39-04 are specified to consume verbatim.
- CCE CSV/header shape ground truth (array vs scalar) is now empirically confirmed from source, removing ambiguity for any future CCE-related work.
- No blockers for 39-03 (page plans) or 39-04.
- **Requirements note:** UI-03/UI-04/UI-05/UI-06 remain `[ ]` Pending in REQUIREMENTS.md after this plan. 39-02 provides the shared substrate (Plotly builders + CSV serializer) those requirements depend on, but the requirement text describes end-to-end page behavior (clicking "Run simulation" and seeing a plot, downloading a CSV from a page) that does not exist until 39-03 (UI-03, UI-04, UI-06) and 39-04 (UI-05, UI-06) build the actual Streamlit pages. Do not mark them complete until those plans execute.

---

_Phase: 39-c-v-cce-field-map-pages-csv-download_
_Completed: 2026-07-11_

## Self-Check: PASSED

- FOUND: app/components/results.py
- FOUND: tests/test_app_csv_export.py
- FOUND: .planning/phases/39-c-v-cce-field-map-pages-csv-download/39-02-SUMMARY.md
- FOUND commit aeb2126 (feat: Plotly figure builders)
- FOUND commit 53236ef (test: to_csv_bytes unit test)
- Re-ran plan-level `<verification>`: `uv run pytest tests/test_app_csv_export.py -x` → 4 passed; `uv run python -c "from app.components.results import build_cv_figure, build_cce_figure, build_field_figures, build_mott_schottky_figure, to_csv_bytes"` → exit 0; `grep -c "import streamlit" app/components/results.py tests/test_app_csv_export.py` → 0 for both
