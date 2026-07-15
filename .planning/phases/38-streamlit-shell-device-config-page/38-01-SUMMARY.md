---
phase: 38-streamlit-shell-device-config-page
plan: 01
subsystem: ui
tags: [streamlit, pytest, dataclass, tdd]

# Dependency graph
requires:
  - phase: 37-core-api-cce-facades-parametricsweep
    provides: "etna.DeviceConfig public API (11-field dataclass) and stable facade contract"
provides:
  - "Pure assemble_config(values: dict) -> DeviceConfig seam, Streamlit-free, unit-tested"
  - "render_device_sidebar() Streamlit renderer writing st.session_state['device_config']"
  - "pytest pythonpath config so app.* imports resolve under pytest"
  - "UI-02 unit coverage (config assembly, doping-mode mapping, dimensionality mapping)"
  - "UI-07 persistence-key contract test (session_state round-trip)"
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
    - "Streamlit-free pure seam pattern: assemble_config(values: dict) -> DeviceConfig has zero st.* calls, imported directly by unit tests; render_device_sidebar() is the only function touching st.*"
    - "Single non-widget session_state key (device_config) as the cross-page persistence contract, distinct from ephemeral widget-keyed state"
    - "Reactive (non-form) mode selectors: dimensionality/doping-profile render outside any submit-gated container so gated fields toggle immediately"

key-files:
  created:
    - app/__init__.py
    - app/components/__init__.py
    - app/components/device_sidebar.py
    - tests/test_app_device_sidebar.py
    - tests/test_app_session.py
  modified:
    - pyproject.toml

key-decisions:
  - "assemble_config only enforces the doping_profile -> N_D consistency rule (graded forces N_D=None); it does NOT reset the graded triplet or half_width_um when the other mode is active — resetting hidden/gated fields to defaults is render_device_sidebar's job, keeping the pure seam a straightforward field-mapper that never silently drops caller-supplied values"
  - "Persistence test (test_config_persistence_key) asserts against an independently-constructed expected DeviceConfig built directly from the same input values, not against assemble_config's own output — this avoids a vacuous self-round-trip that would pass against the Task-2 stub"

patterns-established:
  - "Pure logic / Streamlit-renderer split: every future app/components/*.py should isolate assembly logic testable without a running Streamlit server, mirroring assemble_config vs render_device_sidebar"

requirements-completed: [UI-02, UI-07]

# Metrics
duration: 25min
completed: 2026-07-10
---

# Phase 38 Plan 01: Streamlit Config Core Summary

**Pure `assemble_config(values: dict) -> DeviceConfig` seam plus `render_device_sidebar()` Streamlit renderer storing device config under a single `st.session_state["device_config"]` key, with pytest pythonpath wiring and 4 passing unit tests locking UI-02/UI-07.**

## Performance

- **Duration:** ~25 min
- **Completed:** 2026-07-10T09:56Z
- **Tasks:** 3 completed
- **Files modified:** 6 (2 created packages, 1 new component, 2 new test files, 1 modified pyproject.toml)

## Accomplishments

- `pyproject.toml` carries `[tool.pytest.ini_options] pythonpath = ["."]` so `app.*` imports resolve under pytest, verified to work even from outside the project root and independent of the pre-existing `pytest.ini` (which only carries the `slow` marker registration)
- `app/components/device_sidebar.py` implements a pure `assemble_config(values: dict) -> DeviceConfig` (zero Streamlit dependency, unit-tested directly) and `render_device_sidebar()` (reactive mode selectors, grouped geometry/doping/operating widgets with `format="%.3e"` on all doping/area fields, writes the single `st.session_state["device_config"]` key)
- 4 unit tests pass: `test_assemble_config_all_fields`, `test_doping_mode_mapping`, `test_dimensionality_mapping` (UI-02), `test_config_persistence_key` (UI-07)
- TDD RED/GREEN gate sequence followed and verified in git log: `test(38-01)` commit (RED, all 4 tests failing on genuine assertion errors against the default-returning stub) followed by `feat(38-01)` commit (GREEN, all 4 passing)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add pytest pythonpath config + app package scaffold** - `0c821e1` (chore)
2. **Task 2: Stub assemble_config + write failing unit tests (RED)** - `e69bcf9` (test)
3. **Task 3: Implement assemble_config + render_device_sidebar (GREEN)** - `f51c162` (feat)

_TDD plan: RED (Task 2) -> GREEN (Task 3), no REFACTOR commit needed — GREEN implementation required no cleanup pass._

## Files Created/Modified

- `pyproject.toml` - added `[tool.pytest.ini_options] pythonpath = ["."]`
- `app/__init__.py` - empty package marker
- `app/components/__init__.py` - empty package marker
- `app/components/device_sidebar.py` - `assemble_config()` pure seam + `render_device_sidebar()` Streamlit renderer
- `tests/test_app_device_sidebar.py` - UI-02 unit coverage (3 tests: all-fields assembly, doping-mode mapping, dimensionality mapping)
- `tests/test_app_session.py` - UI-07 persistence-key contract test

## Decisions Made

- `assemble_config` applies only the `doping_profile` -> `N_D` consistency rule; it passes through the graded triplet and `half_width_um` unchanged regardless of mode, so it never silently overwrites caller-supplied values for fields that happen to be hidden in the UI. Rationale: keeps the pure seam a faithful field-mapper; the responsibility for "what value does a hidden/gated field carry" belongs to the renderer that decides visibility, not the assembly function that just marshals whatever dict it receives. This also keeps `test_assemble_config_all_fields` (uniform mode, all 11 fields non-default) meaningful — it would fail if assemble_config reset the graded triplet to defaults in uniform mode.
- `test_config_persistence_key` asserts the round-tripped config against an independently-built `DeviceConfig(**NON_DEFAULT_VALUES)`, not against `assemble_config`'s own return value compared to itself. A self-referential round-trip assertion would pass vacuously even against the Task-2 default-returning stub (default equals itself after a dict round-trip), defeating the RED gate. Verified empirically: with this anchor, the test failed genuinely (`epi_thickness_um` mismatch, 10.0 vs 5.0) against the stub.

## Deviations from Plan

**1. [Rule 1 - Bug] Reworded "st.form" mentions in comments/docstrings to avoid tripping the literal `grep -c "st.form"` acceptance check**

- **Found during:** Task 3 verification
- **Issue:** The plan's acceptance criterion `grep -c "st.form" app/components/device_sidebar.py` must return 0, but explanatory comments/docstrings describing that the sidebar is "NOT wrapped in st.form" contained the literal substring `st.form`, so the initial implementation returned a count of 3 despite having zero actual `st.form(...)` calls.
- **Fix:** Reworded the comments/docstrings to describe the same fact ("not wrapped in a submit-gated container") without using the literal string `st.form`.
- **Files modified:** app/components/device_sidebar.py
- **Verification:** `grep -c "st.form" app/components/device_sidebar.py` now returns 0; tests still pass (4/4).
- **Committed in:** f51c162 (Task 3 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - cosmetic/verification-script fix, no behavior change)
**Impact on plan:** No scope creep — text-only edit to satisfy a literal grep-based acceptance check; no functional change to the sidebar.

## Issues Encountered

- Confirmed the pre-existing `tests/pytest.ini` (registers only the `slow` marker) does not shadow or conflict with the new `[tool.pytest.ini_options]` block in `pyproject.toml` — `pythonpath = ["."]` resolves `app.*` imports correctly in both the presence and temporary absence of `pytest.ini`, and even when pytest is invoked from a directory outside the project root. No action needed; documented here for downstream phases relying on the same import convention.

## Next Phase Readiness

- Phases 39-43 can now read `st.session_state["device_config"]` as their input contract — the key is written by `render_device_sidebar()` and its shape is locked by the 4 unit tests in this plan.
- `app/main.py` (entry script wiring `st.navigation`, page registration, and calling `render_device_sidebar()` before `pg.run()`) is NOT part of this plan — it is out of scope for 38-01 and expected in a subsequent plan of Phase 38 or the start of Phase 39's page work. No blocker: `render_device_sidebar()` is a complete, tested, importable function ready to be wired into an entry script.
- No blockers.

---

_Phase: 38-streamlit-shell-device-config-page_
_Completed: 2026-07-10_

## Self-Check: PASSED

All created files verified present on disk: `app/__init__.py`, `app/components/__init__.py`, `app/components/device_sidebar.py`, `tests/test_app_device_sidebar.py`, `tests/test_app_session.py`. All commit hashes verified present in `git log --oneline --all`: `0c821e1`, `e69bcf9`, `f51c162`, `554ad48`.
