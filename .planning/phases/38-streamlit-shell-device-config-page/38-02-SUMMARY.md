---
phase: 38-streamlit-shell-device-config-page
plan: 02
subsystem: ui
tags: [streamlit, st.navigation, apptest, multi-page]

# Dependency graph
requires:
  - phase: 38-streamlit-shell-device-config-page
    plan: 01
    provides: "assemble_config seam + render_device_sidebar() writing st.session_state['device_config']"
provides:
  - "app/main.py: st.navigation entry script (streamlit run app/main.py)"
  - "Home page (orientation + config summary) + 7 workflow placeholder pages, all with the empty-state guard"
  - ".streamlit/config.toml theme tokens"
  - "tests/test_app_pages.py: UI-01 empty-state guard test + AppTest nav/sidebar smoke"
affects:
  [
    39-cv-cce-field-map-pages,
    40-geometry-viewer,
    41-radiation-damage-dark-current,
    42-microdosimetry-batch-sweep,
  ]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Every page module exposes render() with zero module-level st.* side effects, so pages are safe both as st.Page callables and as AppTest.from_function targets"
    - "Explicit url_path= on every st.Page registration — required because all 8 page modules use the same function name (render), which would otherwise collide under Streamlit's inferred URL-pathname"
    - "AppTest.from_function callables must be self-contained (import the page's render function inside the wrapper body), since from_function executes the callable body in isolation without the caller's module-level imports"

key-files:
  created:
    - app/main.py
    - app/pages/__init__.py
    - app/pages/home.py
    - app/pages/cv.py
    - app/pages/cce.py
    - app/pages/field_map.py
    - app/pages/radiation_damage.py
    - app/pages/dark_current.py
    - app/pages/microdosimetry.py
    - app/pages/batch_sweep.py
    - .streamlit/config.toml
    - tests/test_app_pages.py

key-decisions:
  - "Added explicit url_path= to every st.Page(...) registration (deviation from the plan's literal code sketch, which omitted it) — AppTest boot failed with StreamlitAPIException 'Multiple Pages specified with URL pathname render' because every page module exposes a function literally named render, and Streamlit infers the URL pathname from the callable name when no title-derived override applies cleanly across duplicates"
  - 'test_nav_sidebar_smoke asserts only exception==[] and device_config populated, not a nav-page-count number — AppTest''s public API (streamlit.testing.v1.AppTest) exposes no element type for st.navigation''s sidebar page list, so a numeric page-count assertion would have no real accessor to check against; used the plan''s documented fallback ("otherwise assert the boot succeeded and device_config is set as the minimum smoke")'

patterns-established:
  - "Placeholder-page template (title, guard, st.json summary, phase caption) is now locked; Phases 39-42 replace the guard-then-caption body with real facade calls on the same render() functions"

requirements-completed: [UI-01, UI-07]

# Metrics
duration: 35min
completed: 2026-07-10
---

# Phase 38 Plan 02: Streamlit Shell + Pages Summary

**`app/main.py` as the `st.navigation` entry script wiring `render_device_sidebar()` before `pg.run()`, Home + 7 workflow placeholder pages each with the empty-state guard, theme tokens in `.streamlit/config.toml`, and AppTest-based UI-01 test coverage.**

## Performance

- **Duration:** ~35 min
- **Completed:** 2026-07-10
- **Tasks:** 3 completed
- **Files modified:** 12 (10 created page/entry files, 1 theme config, 1 test file)

## Accomplishments

- `.streamlit/config.toml` sets the exact theme tokens from 38-UI-SPEC (`primaryColor="#1F6FEB"`, `backgroundColor="#FFFFFF"`, `secondaryBackgroundColor="#F0F2F6"`, `textColor="#1A1A1A"`)
- `app/pages/home.py` renders the orientation line + config summary (guarded)
- 7 workflow placeholder pages (`cv`, `cce`, `field_map`, `radiation_damage`, `dark_current`, `microdosimetry`, `batch_sweep`) each implement the empty-state guard verbatim and a `st.caption` naming the phase where real behavior lands (39 / 41 / 42 per the page_registration_contract)
- `app/main.py` is the `st.navigation` entry script: `sys.path` fix first, `set_page_config` as the first Streamlit call, `render_device_sidebar()` called before `pg.run()` (config persists across nav — UI-07), all 8 pages registered via callable `st.Page(...)` with Home as `default=True`
- `tests/test_app_pages.py` adds `test_empty_state_guard` (AppTest.from_function proving no crash + guard text fires on empty session_state) and `test_nav_sidebar_smoke` (AppTest.from_file boots the full app, confirms `device_config` is set)
- All 6 `tests/test_app_*.py` tests pass (3 from 38-01 + 2 new + 1 session persistence)
- `AppTest.from_file("app/main.py").run()` completes with no exception

## Task Commits

Each task was committed atomically:

1. **Task 1: Theme config.toml + Home page + 7 placeholder pages with empty-state guard** - `66b9388` (feat)
2. **Task 2: st.navigation entry script (app/main.py) wiring sidebar + all 8 pages** - `bab05f0` (feat)
3. **Task 3: UI-01 empty-state guard test + AppTest nav/sidebar smoke** - `8e1821d` (test)

## Files Created/Modified

- `.streamlit/config.toml` - theme tokens
- `app/pages/__init__.py` - empty package marker
- `app/pages/home.py` - Home landing page (orientation + config summary)
- `app/pages/cv.py`, `cce.py`, `field_map.py`, `radiation_damage.py`, `dark_current.py`, `microdosimetry.py`, `batch_sweep.py` - workflow placeholders, each with the empty-state guard + phase caption
- `app/main.py` - `st.navigation` entry script
- `tests/test_app_pages.py` - UI-01 test coverage

## Decisions Made

- Added `url_path=` explicitly to every `st.Page(...)` registration. The plan's code sketch (from 38-RESEARCH Pattern 1 and the plan's own action text) did not mention it, but every page module exposes a function literally named `render` (per the `verified_import_convention` in the plan), and `AppTest.from_file("app/main.py").run()` failed at boot with `StreamlitAPIException: Multiple Pages specified with URL pathname render` — Streamlit infers the URL pathname from the callable's `__name__` when not otherwise disambiguated, and 8 identically-named callables collided. Explicit `url_path=` per page (e.g. `"home"`, `"cv"`, `"field-map"`) resolves this without changing any page module's function name (which the plan's verified import convention requires to stay `render`).
- `test_empty_state_guard` wraps `from app.pages.cv import render; render()` inside a local function passed to `AppTest.from_function`, rather than passing `app.pages.cv.render` directly. `AppTest.from_function`'s documented contract requires the callable's body to be fully self-contained ("must include any necessary imports") — `render()` alone references the module-level `st` import from `app/pages/cv.py`, which is not visible when Streamlit re-executes just the function body in isolation; without the wrapper, the test failed with `NameError: name 'st' is not defined`, confirmed empirically before settling on the wrapper. The wrapper still calls the real `app.pages.cv.render` (no duplicated guard logic), so it faithfully tests the shipped page.
- `test_nav_sidebar_smoke` does not assert a specific page count. Explored `at.get("page_link")` and the full `AppTest` public attribute list — there is no element type in `streamlit.testing.v1.AppTest` that exposes `st.navigation`'s registered page list or count. Per the plan's own fallback clause ("otherwise assert the boot succeeded and device_config is set as the minimum smoke"), the test asserts `at.exception == []` and `device_config` populated in `session_state`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking issue] Added explicit `url_path=` to every `st.Page(...)` registration**

- **Found during:** Task 2 verification (`AppTest.from_file("app/main.py").run()`)
- **Issue:** The app failed to boot with `StreamlitAPIException: Multiple Pages specified with URL pathname render. URL pathnames must be unique.` — all 8 page modules export a function named `render` (required by the plan's `verified_import_convention`), and Streamlit's default URL-pathname inference collided across all of them.
- **Fix:** Added an explicit, distinct `url_path=` string to each of the 8 `st.Page(...)` calls in `app/main.py` (e.g. `"home"`, `"cv"`, `"field-map"`, `"radiation-damage"`).
- **Files modified:** `app/main.py`
- **Verification:** `AppTest.from_file("app/main.py").run()` now completes with no exception; `device_config` is set in `session_state` after boot.
- **Committed in:** `bab05f0` (Task 2 commit)

**2. [Rule 3 - Blocking issue] Wrapped the page render call for `AppTest.from_function` instead of passing the bare callable**

- **Found during:** Task 3 implementation, exploratory verification before writing the final test
- **Issue:** `AppTest.from_function(render)` where `render` was imported at module scope in the test file raised `NameError: name 'st' is not defined` — `AppTest.from_function` re-executes only the callable's own body in isolation and requires the body to carry its own imports; it does not inherit the test module's `import streamlit as st`.
- **Fix:** Defined a small local wrapper function (`def _run_cv_page(): from app.pages.cv import render; render()`) and passed the wrapper to `AppTest.from_function`. This satisfies the "self-contained" requirement while still exercising the real, shipped `app.pages.cv.render` guard logic (no duplicated assertion text).
- **Files modified:** `tests/test_app_pages.py`
- **Verification:** `test_empty_state_guard` passes; `at.exception == []` and the guard's `st.info` text is present.
- **Committed in:** `8e1821d` (Task 3 commit)

---

**Total deviations:** 2 auto-fixed (both Rule 3 — blocking issues discovered during verification, fixed inline, no scope creep)
**Impact on plan:** No architectural change. Both fixes are narrow, verified corrections to make the plan's own literal code sketches boot/run correctly under the real `AppTest` API contract; page module names, the guard contract, and the sidebar contract from 38-01 are unchanged.

## Issues Encountered

None beyond the two documented deviations above, both resolved without blocking.

## Next Phase Readiness

- `streamlit run app/main.py` boots (verified via `AppTest.from_file` headless equivalent — no exception); `st.navigation` lists Home + all 7 workflow pages.
- The device-config sidebar renders on every page (called before `pg.run()` in the entry script), so `device_config` persists across navigation.
- Every page (Home + 7 workflow) implements the empty-state guard; `test_empty_state_guard` proves no crash on empty `session_state`.
- Phases 39-43 can now add real facade calls inside the existing `render()` functions in `app/pages/*.py` — the guard, title, and page registration scaffolding are already in place and require no restructuring.
- No blockers.

---

_Phase: 38-streamlit-shell-device-config-page_
_Completed: 2026-07-10_

## Self-Check: PASSED

All created files verified present on disk: `app/main.py`, `app/pages/__init__.py`, `app/pages/home.py`, `app/pages/cv.py`, `app/pages/cce.py`, `app/pages/field_map.py`, `app/pages/radiation_damage.py`, `app/pages/dark_current.py`, `app/pages/microdosimetry.py`, `app/pages/batch_sweep.py`, `.streamlit/config.toml`, `tests/test_app_pages.py`. All commit hashes verified present in `git log --oneline`: `66b9388`, `bab05f0`, `8e1821d`.
