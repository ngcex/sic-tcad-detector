---
phase: 38-streamlit-shell-device-config-page
plan: 03
subsystem: ui
tags: [streamlit, st.navigation, widget-state, pytest, apptest]

requires:
  - phase: 38-streamlit-shell-device-config-page (plans 01, 02)
    provides: render_device_sidebar(), st.navigation shell, app/pages workflow modules
provides:
  - app/workflows/ package (renamed from app/pages/) — removes Streamlit magic-multipage collision
  - Explicit key= on all 12 device-config sidebar widgets
  - Regression coverage for first-rerun widget-value persistence (test_first_edit_survives_rerun)
  - Structural guard against reintroducing a sibling app/pages directory (test_no_magic_pages_directory)
affects:
  [
    39-cv-cce-fieldmap-pages,
    40-geometry-viewer,
    41-radiation-dark-current,
    42-microdosimetry-batch-sweep,
  ]

tech-stack:
  added: []
  patterns:
    - "Streamlit entry-script sibling directories must never be named literally 'pages' when st.navigation is used explicitly — triggers legacy magic-multipage auto-detection collision"
    - "Explicit key= on every interactive widget as defense-in-depth against implicit widget-ID instability"

key-files:
  created: []
  modified:
    - app/main.py
    - app/components/device_sidebar.py
    - app/workflows/__init__.py
    - app/workflows/home.py
    - app/workflows/cv.py
    - app/workflows/cce.py
    - app/workflows/field_map.py
    - app/workflows/radiation_damage.py
    - app/workflows/dark_current.py
    - app/workflows/microdosimetry.py
    - app/workflows/batch_sweep.py
    - tests/test_app_pages.py

key-decisions:
  - "Root-caused UI-07 first-edit-loss defect to app/pages/ colliding with Streamlit's legacy magic-multipage auto-detection alongside explicit st.navigation, not to missing widget keys alone (keys alone, without the rename, do not fix the real app — verified empirically by plan-checker before this plan was finalized)."
  - "Used git mv to preserve file history across the directory rename."
  - "Chose a structural guard test (app/pages must not exist) over re-asserting raw Streamlit proto.id, since once explicit keys exist the proto.id is no longer the discriminating signal — the directory absence is."

patterns-established:
  - "RED-then-GREEN gap-closure plans: write the regression test first, prove it fails against unmodified code (with an automated pytest-summary-line gate, not a raw substring grep), only then implement the fix."

requirements-completed: [UI-02, UI-07]

duration: ~15min (automated tasks) + interactive browser verification
completed: 2026-07-10
---

# Phase 38: Streamlit Shell + Device Config Page — Gap Closure (38-03) Summary

**Closed the UI-07 persistence gap: a user's first sidebar edit after cold app boot no longer reverts to its default on the next Streamlit rerun.**

## Performance

- **Tasks:** 3 completed (2 automated + 1 human-verify checkpoint)
- **Files modified:** 12 (9 renamed app/pages→app/workflows + app/main.py + device_sidebar.py + test_app_pages.py)

## Accomplishments

- Root-caused and fixed the Streamlit magic-multipage / `st.navigation` collision by renaming `app/pages/` → `app/workflows/`
- Added explicit `key=` to all 12 device-config sidebar widgets as defense in depth
- Added a real regression test (`test_first_edit_survives_rerun`) that reproduces the exact `38-VERIFICATION.md` repro against the real `app/main.py`, plus a structural guard (`test_no_magic_pages_directory`) protecting the root-cause fix from silent regression
- Verified end-to-end in a real browser (Playwright): edited Epi thickness to 7.5 as the very first interaction after cold load, navigated Home → C-V Analysis → Home via sidebar links, confirmed the config summary still showed `epi_thickness_um: 7.5`

## Task Commits

1. **Task 1: Write failing first-rerun persistence + magic-pages guard tests (RED)** - `82ccfc7` (test)
2. **Task 2: Rename app/pages→app/workflows, fix imports, add 12 explicit sidebar keys (GREEN)** - `f7c76e9` (feat)
3. **Task 3: Human-verify first-edit persistence across navigation in a real browser** - checkpoint, verified via Playwright (see below)

**Plan metadata:** this commit (docs: complete plan)

## Files Created/Modified

- `app/workflows/*.py` (9 files, renamed from `app/pages/` via `git mv`) - workflow page modules, unchanged render logic
- `app/main.py` - 8 imports changed from `app.pages.X` to `app.workflows.X`; `url_path` values on `st.Page(...)` registrations unchanged
- `app/components/device_sidebar.py` - explicit unique `key=` added to all 12 widgets (`cfg_dimensionality`, `cfg_doping_profile`, `cfg_epi_thickness_um`, `cfg_substrate_thickness_um`, `cfg_half_width_um`, `cfg_N_A`, `cfg_N_D`, `cfg_N_D_junction`, `cfg_N_D_bulk`, `cfg_L_transition_um`, `cfg_T`, `cfg_area_cm2`)
- `tests/test_app_pages.py` - added `test_first_edit_survives_rerun` and `test_no_magic_pages_directory`; updated stale `app.pages.cv` references to `app.workflows.cv`

## Decisions Made

- Root cause is the `app/pages/` directory name colliding with Streamlit's legacy magic-multipage detection alongside the explicit `st.navigation` API — not primarily the missing `key=` attributes (confirmed empirically: `key=` alone without the rename does not fix the real app; see 38-VERIFICATION.md and plan-checker findings recorded in `38-03-PLAN.md`).
- Widget keys still added as defense-in-depth per VERIFICATION.md's missing-item #2, scoped to _all_ 12 widgets (not just mode-toggle-gated ones, correcting the under-scoped WR-01 finding from `38-REVIEW.md`).

## Deviations from Plan

None — plan executed exactly as written. The plan itself was revised once before execution (see `3249f5f`, applied to `38-03-PLAN.md` prior to this SUMMARY) to fix a self-contradicting RED-gate grep command found by `gsd-plan-checker`; no deviation occurred during execution of the corrected plan.

## Issues Encountered

- The RED-gate verify command in the original plan draft (before this SUMMARY's tasks ran) grepped case-insensitively for "error", which always matches pytest's own `AssertionError` traceback text — this was caught and fixed by `gsd-plan-checker` at planning time, not during execution.
- During interactive browser verification, an initial `browser_navigate` (full `goto`) to `http://localhost:8501/` was mistakenly used to simulate "returning to Home" — this triggers a full page reload that resets `st.session_state` entirely and is not equivalent to a user clicking the sidebar Home link. Re-ran the verification using an actual sidebar-link click for both the outbound and return navigation, which correctly exercises Streamlit's client-side `st.navigation` routing and confirmed persistence.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

UI-07's `st.session_state["device_config"]` persistence contract is now solid across the first post-boot edit and subsequent navigation — Phases 39-43 can consume it without inheriting this defect. All 8 app tests pass (`tests/test_app_pages.py`, `tests/test_app_session.py`, `tests/test_app_device_sidebar.py`). Real-browser verification (Playwright) confirmed the fix reaches the actual user-facing app, not just the AppTest headless harness.

---

_Phase: 38-streamlit-shell-device-config-page_
_Completed: 2026-07-10_
