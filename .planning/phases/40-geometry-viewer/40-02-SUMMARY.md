---
phase: 40-geometry-viewer
plan: 02
subsystem: app-ui
tags: [streamlit, plotly, geometry-viewer, field-map, selectbox, viz, apptest]
requires:
  - "app.components.geometry_viewer.build_geometry_figure (from 40-01)"
  - "app.components.geometry_viewer.QUANTITIES (from 40-01)"
  - "petringa.run_field (module-attribute monkeypatch seam)"
  - "petringa.MeshData (top-level re-export)"
provides:
  - "app/workflows/field_map.py: 2D-through routing, dimensionality-branched render, quantity selectbox feeding build_geometry_figure"
  - "tests/test_app_field_page.py: 2D-route + selectbox-present + selectbox-no-resolve (VIZ-03) coverage"
affects:
  - "Phase 40 is now user-observable: VIZ-01/VIZ-02/VIZ-03 reachable on the field-map page"
tech-stack:
  added: []
  patterns:
    - "dimensionality-branched render (result.mesh.y_coords is None) — 1D line charts + CSV + supplemental bar vs 2D heatmap only"
    - "persistent-key st.selectbox reads cached session_state result; run_field call-counter proves no re-solve on quantity change"
key-files:
  created: []
  modified:
    - app/workflows/field_map.py
    - tests/test_app_field_page.py
decisions:
  - "A1 1D geometry bar SUPPLEMENTS the existing line charts (does not replace)"
  - "A2 2D branch SKIPS the CSV download (to_csv_bytes has no 2D branch; 2D CSV out of scope)"
  - "A6 dropdown default = 'Electric field' (first QUANTITIES key, index=0)"
metrics:
  duration: ~10min
  tasks: 2
  files: 2
  completed: 2026-07-13
---

# Phase 40 Plan 02: Field-Map Wiring Summary

Wired the plan 40-01 geometry viewer into the field-map page and removed the Phase-40 blocker: the 2D `st.stop()` pre-check is deleted (both 1D and 2D configs now flow through `run_field`), the render branches on `result.mesh.y_coords is None`, and a quantity `st.selectbox` feeds `build_geometry_figure(result.mesh, quantity)` — re-rendering from the cached result with no devsim re-solve on a quantity change (VIZ-01/VIZ-02/VIZ-03).

## What Was Built

- **`app/workflows/field_map.py`** (modified):
  - DELETED the 2D stop-guard (`if cfg.half_width_um is not None: st.warning("...1D-only..."); st.stop()`). 2D now routes through `petringa.run_field`.
  - REWROTE the module docstring to describe 2D-through routing, the `result.mesh.y_coords is None` render branch, and the preserved `try/except RuntimeError` convergence guard.
  - Added `from app.components.geometry_viewer import build_geometry_figure, QUANTITIES`.
  - Render is now guarded by `if result is not None and result.mesh is not None:` and branches:
    - **1D** (`y_coords is None`): keeps `build_field_figures` line charts, the `to_csv_bytes` download button, and the net-doping expander UNCHANGED, THEN appends the supplemental geometry bar (A1).
    - **2D** (`y_coords is not None`): SKIPS line charts and CSV (empty `result.x`/`result.y`; `to_csv_bytes` has no 2D branch), renders only the geometry heatmap (A2).
  - BOTH branches: `quantity = st.selectbox("Quantity", list(QUANTITIES.keys()), index=0, key="geo_quantity")` then `st.plotly_chart(build_geometry_figure(result.mesh, quantity))`. The persistent key + cached-result reread satisfies VIZ-03 "without re-running" (Run button returns False on a selectbox rerun).
  - PRESERVED: `try/except RuntimeError` around `st.session_state["field_result"] = petringa.run_field(cfg)` and the module-attribute `petringa.run_field` reference (no `from petringa import run_field`). No `@st.cache_data` added (DeviceConfig unhashable).

- **`tests/test_app_field_page.py`** (rewritten):
  - Updated docstring: Phase 40 routes 2D through, fakes return a populated mesh, VIZ-03 verified via a call counter.
  - `_fake_run_field` (1D) now returns a POPULATED mesh (`y_coords=None`, three node_values keys) so the viewer branch renders.
  - Added `_fake_run_field_2d`: empty `x`/`y`, empty `metadata`, populated 2D mesh (irregular 5-node scatter, `y_coords` non-None).
  - DELETED `test_2d_config_warns_and_skips` (asserted the removed `st.stop()`).
  - Added `test_2d_config_routes_through_and_caches`: `half_width_um=50.0`, asserts `at.exception == []`, `"field_result" in at.session_state`, and no "1D-only" warning.
  - Added `test_quantity_selectbox_present`: default `at.selectbox[0].value == "Electric field"`; all three quantities selectable via the verified `.select()` accessor with `at.exception == []`.
  - Added `test_selectbox_change_does_not_resolve` (VIZ-03 core): call counter stays at 1 across a `at.selectbox[0].select("Net doping").run()`, and the cached `field_result` is the same object (no re-solve).
  - Kept the three legacy tests: `test_run_caches_field_result`, `test_empty_state_guard`, `test_solver_convergence_failure_shows_error_not_crash`.

## Verification

- `uv run pytest tests/test_app_field_page.py -x` → 6 passed.
- `grep -c '1D-only' app/workflows/field_map.py` → 0 (2D guard removed).
- `grep -c 'test_2d_config_warns_and_skips' tests/test_app_field_page.py` → 0 (obsolete test removed).
- `grep -c 'mesh=None' tests/test_app_field_page.py` → 0 (both fakes populate mesh).
- Field-map page contains `build_geometry_figure(result.mesh, quantity)` and `st.selectbox("Quantity", ..., index=0, key="geo_quantity")`.
- Task 1 AST/grep verify command → `OK`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] `uv run pytest` failed with `ModuleNotFoundError: No module named 'plotly'`**

- **Found during:** Task 2 verification (first `uv run pytest` after `uv venv` was freshly created by the first `uv run`).
- **Issue:** The project venv did not yet have the `dev` optional-dependency group materialized, so plotly (a runtime dep pulled via the dev sync path) was absent and the page import failed under test collection — the exact pitfall documented in STATE.md:106 and the 40-01 summary.
- **Fix:** Ran `uv sync --extra dev` to materialize the extras into the venv; the plan's verify command `uv run pytest tests/test_app_field_page.py -x` then passed (6 green). No source change; `uv.lock` unchanged (extras already locked, only not materialized).
- **Files modified:** none (source); venv materialization only.
- **Commit:** n/a (env step, no tracked change).

## Known Stubs

None — the page is wired against the real `MeshData` interface and the 40-01 builder. 2D end-to-end browser verification remains best-effort (a non-converging 2D `ramp_bias` is upstream physics, out of scope for this UI phase, per 40-RESEARCH §Environment); the `try/except RuntimeError` + `result.mesh is not None` guard handle that gracefully.

## Threat Flags

None — no new security surface. The quantity selectbox is a fixed hardcoded option set (`list(QUANTITIES.keys())`) and `build_geometry_figure` does `QUANTITIES[quantity]` (KeyError on any non-member), so no arbitrary key reaches `node_values` (T-40-03 mitigated by construction). No file upload, no free text, no 2D CSV write.

## Self-Check: PASSED

- FOUND: app/workflows/field_map.py
- FOUND: tests/test_app_field_page.py
- FOUND commit 03b4b64 (feat: route 2D through, branch render, selectbox)
- FOUND commit 5a0ebc1 (test: 2D route + selectbox + no-resolve)
