---
phase: 42-microdosimetry-page-batch-sweep-page
plan: 02
subsystem: ui
tags:
  [
    streamlit,
    microdosimetry,
    file-upload,
    tempfile,
    plotly,
    csv-export,
    apptest,
  ]

# Dependency graph
requires:
  - phase: 42-microdosimetry-page-batch-sweep-page
    provides: 42-01's build_microdosimetry_figure (log-x spectrum builder) + microdosimetry branch of to_csv_bytes in app/components/results.py
  - phase: 41-radiation-damage-dark-current-pages
    provides: app/workflows/radiation_damage.py Run->cache->render->download page skeleton + module-attribute mockability seam
provides:
  - app/workflows/microdosimetry.py — full upload -> tempfile bridge -> run_microdosimetry -> spectrum + y_F/y_D caption + CSV download page (FEAT-03)
  - tests/test_app_microdosimetry_page.py — AppTest coverage driving st.file_uploader via streamlit 1.58 .upload(filename, content, mime_type)
affects: [microdosimetry-page, app-main-navigation]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "st.file_uploader -> tempfile.NamedTemporaryFile bridge: in-memory UploadedFile bytes written to a SERVER-generated temp path (no traversal), removed in finally (single-request lifetime)"
    - "AppTest file-upload driving: at.file_uploader[0].upload(filename, content, mime_type) — THREE positional args in streamlit 1.58 (not a single tuple)"

key-files:
  created:
    - tests/test_app_microdosimetry_page.py
  modified:
    - app/workflows/microdosimetry.py

key-decisions:
  - "FileUploader.upload takes three POSITIONAL args (filename, content, mime_type) — the plan/PATTERNS tuple form ((name, bytes, mime)) was wrong; verified against the installed streamlit.testing.v1 signature and adapted per the plan's explicit 'confirm and adapt' instruction"
  - "Dropped the radiation_damage.py 1D half_width_um guard, NaN-tolerance handling, and kappa banner — run_microdosimetry is a config-independent pure pipeline with no partial-convergence mode"
  - "y_F/y_D readout uses st.caption (not st.metric) per UI-SPEC Typography Discretion; docstring reworded to avoid the literal 'st.metric' substring that tripped the naive inline verify"

patterns-established:
  - "Pattern: server-generated tempfile bridge for the app's only file-upload surface — os.remove in finally guarantees no disk accumulation and no path-traversal surface"

requirements-completed: [FEAT-03]

# Metrics
duration: 15 min
completed: 2026-07-14
---

# Phase 42 Plan 02: Microdosimetry Page Summary

**Full microdosimetry Streamlit page — CSV upload bridged to run_microdosimetry via a server-generated tempfile (os.remove in finally), y·d(y) log-x spectrum, y_F/y_D st.caption readout, and single-result CSV download — plus a 4-case AppTest suite driving the file uploader with the synthetic MC fixture.**

## Performance

- **Duration:** 15 min
- **Started:** 2026-07-14T17:25:00Z (approx)
- **Completed:** 2026-07-14
- **Tasks:** 2
- **Files modified:** 2 (1 modified, 1 created)

## Accomplishments

- Replaced the `app/workflows/microdosimetry.py` placeholder with the app's first (and only) file-upload page: `st.file_uploader(type=["csv"])` -> `tempfile.NamedTemporaryFile` bridge -> `etna.run_microdosimetry(cfg, mc_csv_path=...)` -> cache in `st.session_state["microdosimetry_result"]` -> spectrum + `st.caption` y_F/y_D readout + `st.download_button` via `to_csv_bytes`.
- Wired the security-critical tempfile lifecycle: server-generated path (no user-supplied path, no traversal), `os.remove` in a `finally` block guarded by `os.path.exists` (single-request temp lifetime, no disk accumulation) — mitigations T-42-02-T1/T2/D1/D2 all present.
- Added `tests/test_app_microdosimetry_page.py` (4 tests): empty-state guard, no-file warning, upload->Run->cache happy path (driven with `data/synthetic_mc_events.csv` bytes), and malformed-CSV error-not-crash.
- Full plan verification green: microdosimetry page (4) + frozen `test_api_microdosimetry.py` contract + 42-01 `test_app_csv_export.py` (14 tests total pass).

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement the microdosimetry page** - `dc29016` (feat)
2. **Task 2: AppTest coverage for the microdosimetry page** - `1bde76f` (test)

**Plan metadata:** (final docs commit — see below)

_Task 2 is a `tdd="true"` task, but by plan design Task 1 builds the page first, so the tests are test-after and passed immediately on first run (GREEN, no RED phase) — the exact precedent set by 42-01. This is expected, not a fail-fast trip._

## Files Created/Modified

- `app/workflows/microdosimetry.py` - Fleshed out from placeholder into the full upload/run/render/download page with empty-state, no-file, and malformed-CSV guards; module-attribute `etna.run_microdosimetry` seam; no 1D/NaN/kappa cruft.
- `tests/test_app_microdosimetry_page.py` - New AppTest suite mirroring `test_app_radiation_damage_page.py`; monkeypatches `etna.run_microdosimetry` as a module attribute; drives `at.file_uploader[0].upload(...)` with the synthetic fixture bytes.

## Decisions Made

- **FileUploader.upload signature:** Confirmed via `inspect.signature` that streamlit 1.58's `FileUploader.upload(filename, content, mime_type)` takes THREE positional args — the PLAN/PATTERNS sketch passing a single `(name, bytes, mime)` tuple was wrong. Adapted to `at.file_uploader[0].upload("events.csv", data, "text/csv")` per the plan's explicit "confirm the real API and adapt" clause. Also verified the delivered object is a real `UploadedFile` supporting `.getvalue()`, so the page's bridge works unchanged.
- **Dropped radiation_damage.py cruft:** No 1D `half_width_um` guard, no `np.isnan` NaN-tolerance handling, no kappa data-blocked banner — all radiation-damage-specific and inapplicable to the config-independent pure `run_microdosimetry` pipeline (per PATTERNS DROP list).
- **Docstring wording:** Reworded "NOT st.metric" to "rather than metric tiles" because the inline verify does a naive `'st.metric' not in src` substring check that the docstring text otherwise tripped — behavior (st.caption, no metric call) is unchanged.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Merged main to obtain the declared 42-01 dependency**

- **Found during:** Task 1 (page import verification)
- **Issue:** The orchestrator cut this worktree from a stale base (`5e0b686`) that predated the 42-01 merge on `main` (`1a3051b`). `app/components/results.py` in the worktree lacked `build_microdosimetry_figure`, so the page's import failed. The plan's frontmatter declares `depends_on: [42-01]` and the orchestrator confirmed 42-01 was merged to main.
- **Fix:** Fast-forwarded the worktree-agent branch with `git merge main` (main HEAD was a strict descendant of the merge-base; main touches `results.py`/`test_app_csv_export.py`/docs but NOT `app/workflows/microdosimetry.py`, so the uncommitted Task-1 change was preserved and HEAD stayed on the worktree-agent branch). No rebase (would rewrite history / trip branch guards).
- **Files modified:** none of the task files (merge only pulled in 42-01 + prior phases' already-merged work)
- **Verification:** Task 1 inline verify prints OK; full test suite green post-merge.
- **Committed in:** merge commit (fast-forward, not a task commit)

**2. [Rule 3 - Blocking] Corrected the test runner invocation to `uv run --extra dev pytest`**

- **Found during:** Task 2 (running the AppTest suite)
- **Issue:** `uv run pytest` provisioned pytest into an ephemeral tool environment lacking the project runtime deps, so imports failed with `ModuleNotFoundError: No module named 'plotly'` (plotly is a declared runtime dep; pytest lives in the `dev` optional-dependencies group). This is not a missing/unknown package — both are already in `pyproject.toml`/`uv.lock`.
- **Fix:** Ran tests with `uv run --extra dev pytest ...`, which installs pytest into the project venv alongside plotly. No package was installed, substituted, or renamed.
- **Files modified:** none
- **Verification:** `uv run --extra dev pytest tests/test_app_microdosimetry_page.py -q` -> 4 passed.
- **Committed in:** n/a (invocation change only)

---

**Total deviations:** 2 auto-fixed (both Rule 3 - blocking, both mechanical/environmental)
**Impact on plan:** Both are environment/dependency plumbing, not code-scope changes. No source deviated from the plan's specified behavior; no scope creep.

## Issues Encountered

- A PostToolUse formatter hook reformatted both files after each write (cosmetic — line wrapping); no functional impact, all verifications green.

## User Setup Required

None - no external service configuration required.

## Threat Flags

None — the file-upload surface introduced here is exactly the surface the plan's `<threat_model>` anticipated (T-42-02-T1/T2/D1/D2), all mitigated in the implementation (server-generated tempfile path, `os.remove` in finally, `pd.read_csv`-only parsing behind `type=["csv"]`, try/except on malformed input).

## Next Phase Readiness

- FEAT-03 satisfied: the microdosimetry page is complete, tested, and ready for st.Page registration in `app/main.py` (if not already wired).
- The tempfile-upload bridge pattern is established for any future file-ingesting page.
- Note for the orchestrator/merge step: this worktree was branched from a base predating the 42-01 merge and required a `git merge main` to obtain its declared dependency — Wave 2's sibling plan 42-03 (also `depends_on: [42-01]`) was very likely cut from the same stale base and may need the same fast-forward.

## Self-Check: PASSED

- Files exist: app/workflows/microdosimetry.py, tests/test_app_microdosimetry_page.py — both FOUND
- Commits exist: dc29016 (feat), 1bde76f (test) — both FOUND
- `uv run --extra dev pytest tests/test_app_microdosimetry_page.py tests/test_api_microdosimetry.py tests/test_app_csv_export.py -q` -> 14 passed
- Inline page verify prints OK (no st.metric, os.remove+finally present, module-attribute seam, no half_width_um guard)

---

_Phase: 42-microdosimetry-page-batch-sweep-page_
_Completed: 2026-07-14_
