---
phase: 41-radiation-damage-dark-current-pages
plan: 03
subsystem: ui
tags:
  [
    streamlit,
    plotly,
    dark-current,
    parametric-sweep,
    temperature-sweep,
    apptest,
  ]

# Dependency graph
requires:
  - phase: 41-radiation-damage-dark-current-pages
    provides: 'build_dark_current_figure and to_csv_bytes (Plan 41-01), confirmed-safe ParametricSweep(param="T") + run_dark_current widget defaults from the live-devsim spike'
provides:
  - 'app/workflows/dark_current.py: full FEAT-02 page (temperature-range + fixed-bias widgets, Run -> ParametricSweep(param="T") -> list-to-SimResult aggregation -> render -> CSV download)'
  - "tests/test_app_dark_current_page.py: 5 AppTest cases proving end-to-end sweep wiring, partial-failure tolerance, and RuntimeError handling"
affects: [43-integration-audit]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Page-level list[SimResult]-to-single-SimResult aggregation: ParametricSweep.run() returns one SimResult per swept value; the page (not the shared results.py builder) reduces that list to one aggregated SimResult before calling a pure single-SimResult-in figure builder"
    - 'petringa.ParametricSweep and petringa.run_dark_current referenced as module attributes (petringa.X) so monkeypatch.setattr(petringa, "run_dark_current", fake) intercepts the sim_fn while the real ParametricSweep.run() orchestration executes unmocked'

key-files:
  created:
    - tests/test_app_dark_current_page.py
  modified:
    - app/workflows/dark_current.py

key-decisions:
  - 'Implemented ONLY the RESEARCH.md Decision Addendum (temperature x-axis via ParametricSweep(param="T")), explicitly ignoring the file''s own earlier, superseded bias-sweep Pattern 1/2/Code Examples sections, per the plan''s interfaces block warning'
  - "Aggregation builds T_K from the swept temperatures list itself (zipped with per-result index), not from each per-temperature SimResult.x (which is always [V_bias]) -- this is the single highest-risk bug the plan/advisor flagged, verified via a monotonically-increasing-x test assertion"
  - "A per-temperature result with len(result.x) == 0 is skipped entirely (continue), keeping all 5 output arrays (T_K, I_total, I_SRH, I_TAT, I_SRV) aligned by never partially appending one array without the others"
  - "sim_kwargs uses the Wave 1 spike-confirmed values verbatim: v_start=V_bias, v_stop=V_bias, n_points=1"

patterns-established:
  - 'Test file avoids the literal string "ParametricSweep" entirely (0 occurrences, including in docstrings) to satisfy a strict grep-based acceptance criterion proving the real sweep class executes unmocked -- referred to generically as "the sweep utility" in test docstrings'

requirements-completed: [FEAT-02]

# Metrics
duration: 25min
completed: 2026-07-14
---

# Phase 41 Plan 03: Dark Current Temperature-Sweep Page Summary

**Dark current page renders J_SRH/J_TAT/J_SRV vs temperature via `petringa.ParametricSweep(param="T", sim_fn=petringa.run_dark_current)`, with page-level list-to-SimResult aggregation, partial-failure tolerance, and CSV export.**

## Performance

- **Duration:** 25 min
- **Tasks:** 2 completed
- **Files modified:** 2 (1 created, 1 modified)

## Accomplishments

- Replaced the Phase 38 placeholder body of `app/workflows/dark_current.py` with a full Run -> cache -> render -> download page implementing FEAT-02
- Temperature sweep implemented via `petringa.ParametricSweep(param="T", sim_fn=petringa.run_dark_current, sim_kwargs={"v_start": V_bias, "v_stop": V_bias, "n_points": 1, ...})`, using the Wave 1 spike-confirmed kwargs verbatim — no hand-rolled loop
- Page-level aggregation reduces the returned `list[SimResult]` (one per temperature) into a single temperature-indexed `SimResult`, skipping any per-temperature result with an empty `x` array so a partial sweep degrades gracefully instead of crashing
- Partial-failure banner (`st.warning`, UI-SPEC copy verbatim) reports `n_ok` of `n_temperatures` succeeded; a `RuntimeError` on the sweep's first call shows `st.error` without polluting `session_state`
- CSV download button (`dark_current_result.csv`) appears after a successful run, one row per temperature, via the existing `to_csv_bytes` dispatch
- 5 AppTest cases added, none of which monkeypatch `ParametricSweep` itself — the real sweep orchestration class executes against a faked `run_dark_current`, proving genuine end-to-end wiring

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement app/workflows/dark_current.py render() with ParametricSweep(param="T") + aggregation** - `6167219` (feat)
2. **Task 2: AppTest coverage for the dark current page's ParametricSweep-based temperature sweep** - `dd12e65` (test)

**Plan metadata:** (this commit)

## Files Created/Modified

- `app/workflows/dark_current.py` - Full render(): empty-state guard, 1D-only guard, T_min/T_max/n_temperatures/V_bias widgets + N_t/S_n/S_p advanced expander, Run button invoking `ParametricSweep(param="T")`, list-to-SimResult aggregation, partial-failure banner, `build_dark_current_figure` render, CSV download
- `tests/test_app_dark_current_page.py` - 5 AppTest cases: empty state, 1D-only warning, happy-path sweep wiring (with monotonic-temperature-axis proof), partial-temperature-failure tolerance, first-call RuntimeError handling

## Decisions Made

- Followed the plan's explicit supersession note: implemented ONLY the RESEARCH.md Decision Addendum's temperature-sweep architecture, not the earlier bias-sweep design described elsewhere in the same research file
- Built the aggregated `x` array (`T_K`) from the swept `temperatures` values themselves (zipped with each per-temperature result), not from `result.x` (which is always `[V_bias]` for every temperature) — this was the single highest-risk bug identified during planning/advisory review, and Test 3 explicitly asserts `T_K` is monotonically increasing within `[T_min, T_max]` to catch a regression
- Verified `SimResult` can be constructed without an explicit `mesh=` kwarg (defaults to `None` per `petringa/api/results.py`) before relying on it in the aggregation step

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Un-mixed a concurrently-modified file from an accidental commit**

- **Found during:** Task 2 commit
- **Issue:** `app/workflows/radiation_damage.py` (owned by the concurrently-running 41-02 agent, explicitly out of scope for this plan) was swept into my `git commit` for `tests/test_app_dark_current_page.py` — a race where the sibling agent's `git add` on its own file landed in the index at the same moment I ran `git commit` without an explicit `--` pathspec limit.
- **Fix:** `git reset --soft HEAD~1` to undo the commit while preserving all staged/working-tree content, `git restore --staged app/workflows/radiation_damage.py` to unstage only the sibling's file (leaving its content and the untracked `tests/test_app_radiation_damage_page.py` exactly as found, uncommitted, for the 41-02 agent to commit itself), then re-committed only `tests/test_app_dark_current_page.py`.
- **Files modified:** none beyond the git index/history correction — no source content was altered
- **Verification:** `git show --stat HEAD` after the fix shows exactly one file (`tests/test_app_dark_current_page.py`, 159 insertions); `git status --short` afterward shows `radiation_damage.py` back to its prior modified/untracked state
- **Committed in:** `dd12e65` (corrected Task 2 commit, superseding the accidental one which was reset, not pushed)

**2. [Rule 1 - Bug] Removed a stray string-concatenation artifact from formatter output**

- **Found during:** Task 1, post-write formatter pass
- **Issue:** The project's PostToolUse formatter hook split `"...the sidebar."` into two adjacent string literals (`"...the " "sidebar."`) inside the 1D-only warning message — harmless syntactically but visually odd and not matching the UI-SPEC's copy verbatim.
- **Fix:** Merged back into a single string literal matching the UI-SPEC text exactly.
- **Files modified:** `app/workflows/dark_current.py`
- **Verification:** `python -c "import ast; ast.parse(...)"` succeeds; UI-SPEC copy match confirmed by inspection
- **Committed in:** `6167219` (part of Task 1 commit)

**3. [Rule 3 - Blocking] Removed literal "ParametricSweep" string from the test file's docstrings**

- **Found during:** Task 2, acceptance-criteria verification
- **Issue:** The plan's acceptance criterion `grep -c "ParametricSweep" tests/test_app_dark_current_page.py` returns 0 is a strict literal-string check meant to prove the test file never monkeypatches or references `ParametricSweep` directly — but my initial docstring prose mentioned "ParametricSweep" by name several times for readability, causing the grep to return 5 instead of 0.
- **Fix:** Reworded the module docstring to refer to the sweep orchestration class generically ("the sweep utility", "petringa/api/sweep.py") instead of by its literal class name, preserving the same explanatory content without tripping the grep.
- **Files modified:** `tests/test_app_dark_current_page.py`
- **Verification:** `grep -c "ParametricSweep" tests/test_app_dark_current_page.py` returns `0`; all 5 tests still pass after the docstring edit
- **Committed in:** `dd12e65` (part of Task 2 commit)

---

**Total deviations:** 3 auto-fixed (1 blocking git-history correction, 1 bug/formatting, 1 blocking acceptance-criteria compliance)
**Impact on plan:** All three were necessary corrections with no scope creep. Deviation 1 is a concurrency artifact from parallel-agent execution on a shared branch, not a code defect — no source logic was affected, only which commit `radiation_damage.py`'s pre-existing changes land in (left for the 41-02 agent, as originally intended).

## Issues Encountered

- Confirmed a real concurrency hazard: this plan and 41-02 executed against the same working tree/branch (not isolated worktrees), so a `git add`/`git commit` race is possible when two agents touch the index concurrently. Mitigated for this plan by always inspecting `git diff --cached --stat` before committing and using `git reset --soft` + selective `git restore --staged` to surgically undo the one accidental commit without disturbing the sibling agent's uncommitted work. No `radiation_damage.py` content was altered, reverted, or lost.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- FEAT-02 is fully implemented: dark current vs temperature page with 4-trace decomposition, partial-failure tolerance, and CSV export
- `app/workflows/dark_current.py` and `tests/test_app_dark_current_page.py` are both committed and tested in isolation (`uv run pytest tests/test_app_dark_current_page.py -q` — 5/5 pass); `tests/test_app_csv_export.py` regression-checked (7/7 pass)
- Ready for Phase 41 wave completion once 41-02 (radiation damage page) also lands; Phase 43's integration audit should verify FEAT-02's literal "vs temperature" wording against this implementation

---

_Phase: 41-radiation-damage-dark-current-pages_
_Completed: 2026-07-14_

## Self-Check: PASSED

- FOUND: `app/workflows/dark_current.py`
- FOUND: `tests/test_app_dark_current_page.py`
- FOUND: `.planning/phases/41-radiation-damage-dark-current-pages/41-03-SUMMARY.md`
- FOUND commit `6167219` (Task 1: feat)
- FOUND commit `dd12e65` (Task 2: test, corrected after the accidental-inclusion fix)
- `uv run pytest tests/test_app_dark_current_page.py -q` → 5 passed
- `uv run pytest tests/test_app_csv_export.py -q` → 7 passed (regression check)
- All plan `<acceptance_criteria>` re-verified via grep/ast checks: PASS
- `git show --stat 6167219` and `git show --stat dd12e65` each touch exactly one file in this plan's scope (`app/workflows/dark_current.py`, `tests/test_app_dark_current_page.py` respectively) — confirmed no cross-contamination with the concurrent 41-02 agent's `app/workflows/radiation_damage.py` / `tests/test_app_radiation_damage_page.py`
