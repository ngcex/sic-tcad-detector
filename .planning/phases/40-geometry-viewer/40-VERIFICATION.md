---
phase: 40-geometry-viewer
verified: 2026-07-13T00:00:00Z
status: human_needed
score: 3/3 must-haves verified (code-level); 2 items require human browser confirmation
overrides_applied: 0
human_verification:
  - test: "Run a real 2D field-map simulation in the running Streamlit app (streamlit run app/main.py), set Dimensionality=2D with a shallow bias in the sidebar, click Run"
    expected: "A 2D Plotly heatmap of the electric field renders on the page (device cross-section, lateral position vs depth axes), not the old '1D-only' warning and not a blank/crashed page"
    why_human: "AppTest (Streamlit's headless test harness) has no plotly_chart accessor — chart rendering in the real browser cannot be asserted programmatically. This verifier separately proved the full data path with a real (non-mocked) devsim solve at the Python level (see Data-Flow Trace below), but the actual browser paint is unverified by tooling."
  - test: "In the running app, change the Quantity dropdown among Electric field / Net doping / Electrostatic potential and visually confirm the chart's shape, values, and colorscale actually change"
    expected: "Heatmap/bar visibly updates to a different colorscale (Viridis / Plasma / RdBu_r) and different value range per quantity, without a page reload delay indicating a re-solve"
    why_human: "Visual colorscale/shape change is not asserted by AppTest (no chart-content accessor); this is inherently a rendering/appearance check the automated suite cannot make from Python objects alone."
---

# Phase 40: Geometry Viewer Verification Report

**Phase Goal:** Users can see the electric field (or other quantities) overlaid on a 2D device cross-section after running a field map simulation, with the visualization mode adapting automatically to 1D or 2D device configuration.
**Verified:** 2026-07-13
**Status:** human_needed
**Re-verification:** No — initial verification (a prior `40-VERIFICATION.md` existed at this path but was a pre-execution plan-checker report, not a post-execution goal-backward verification; it is superseded by this report)

## Goal Achievement

### Observable Truths

| #   | Truth                                                                                           | Status                  | Evidence                                                                                                                                                                                                                                                                                        |
| --- | ----------------------------------------------------------------------------------------------- | ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | A pure builder converts 2D MeshData into a Plotly heatmap trace                                 | VERIFIED                | `app/components/geometry_viewer.py` `_build_heatmap`; unit tests pass; **also independently confirmed against a real non-mocked 2D `run_field(DeviceConfig(half_width_um=50.0), bias_V=-20.0)` solve** — see Data-Flow Trace                                                                    |
| 2   | A pure builder converts 1D MeshData into a Plotly bar trace                                     | VERIFIED                | `_build_bar`; unit tests pass; **also confirmed against a real non-mocked 1D `run_field(DeviceConfig(), bias_V=-20.0)` solve**                                                                                                                                                                  |
| 3   | Builder exposes QUANTITIES map of 3 labels -> node_values keys, "Electric field" first          | VERIFIED                | `QUANTITIES = {"Electric field": "ElectricField", "Net doping": "NetDoping", "Electrostatic potential": "Potential"}` (geometry_viewer.py:23-27); `test_quantities_map_contract` asserts exact dict + order                                                                                     |
| 4   | Builder never imports or calls devsim                                                           | VERIFIED                | `grep -c 'import devsim\|devsim\.' app/components/geometry_viewer.py` → 0; module only imports numpy/plotly/scipy/petringa.MeshData                                                                                                                                                             |
| 5   | 2D device no longer hits `st.stop()`; `run_field` runs and result is cached                     | VERIFIED                | `1D-only` string absent from field_map.py (grep → 0); `test_2d_config_routes_through_and_caches` asserts `"field_result" in at.session_state` for `half_width_um=50.0`; **independently confirmed with a real 2D devsim solve** (converges at `bias_V=-20.0`, produces populated 2D `MeshData`) |
| 6   | For a 2D result, a quantity heatmap renders (no empty line charts)                              | VERIFIED                | field_map.py render branches on `result.mesh.y_coords is None`; 2D branch skips `build_field_figures`/`to_csv_bytes` and renders only `build_geometry_figure`; real 2D solve → heatmap trace with 20000/20000 finite cells for all 3 quantities                                                 |
| 7   | For a 1D result, existing line charts still render AND geometry viewer (bar) renders below      | VERIFIED                | 1D branch in field_map.py keeps `build_field_figures` + CSV + net-doping expander unchanged, then appends the selectbox + `build_geometry_figure` bar; `test_run_caches_field_result` (existing 1D path) still green                                                                            |
| 8   | Quantity dropdown lets user switch quantity, re-renders from cache without re-running run_field | VERIFIED (code-level)   | `test_selectbox_change_does_not_resolve` — call-counter stays at 1 across `at.selectbox[0].select("Net doping").run()`; cached `field_result` object identity preserved                                                                                                                         |
| 9   | 2D heatmap visibly appears in the real browser after a real 2D run                              | UNCERTAIN (needs human) | AppTest cannot assert chart paint; Python-level data path independently verified real end-to-end (see below), but in-browser rendering itself is unverified by tooling                                                                                                                          |
| 10  | Dropdown visibly changes chart appearance in the real browser                                   | UNCERTAIN (needs human) | Same AppTest limitation — visual colorscale/shape difference is inherently a human check                                                                                                                                                                                                        |

**Score:** 8/8 code-verifiable truths VERIFIED; 2 truths require human browser confirmation (see Human Verification Required)

### Required Artifacts

| Artifact                            | Expected                                                                                             | Status   | Details                                                                                                                                                                                                                                                                                                                                                        |
| ----------------------------------- | ---------------------------------------------------------------------------------------------------- | -------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `app/components/geometry_viewer.py` | `build_geometry_figure(mesh, quantity)` + `QUANTITIES`, pure MeshData->go.Figure                     | VERIFIED | 111 lines; exports both symbols; no `st.*`, no devsim; griddata heatmap (n_x=200,n_y=100, no fill_value, NaN-transparent) + bar branch; per-quantity scaling (Viridis/RdBu_r/Plasma) matches plan spec exactly                                                                                                                                                 |
| `tests/test_app_geometry_viewer.py` | Pure unit tests: 1D->Bar, 2D->Heatmap, per-quantity scaling, cm->um, doping log10                    | VERIFIED | 8 tests, all pass; covers bar/heatmap trace type, z-shape (100,200), axis titles, reversed y-axis, cm->µm conversion, signed/zero doping log10 (1D+2D), QUANTITIES contract                                                                                                                                                                                    |
| `app/workflows/field_map.py`        | 2D-through routing, dimensionality-branched render, quantity selectbox feeding build_geometry_figure | VERIFIED | Old `st.stop()` 2D guard deleted; render guarded by `result.mesh is not None`, branches on `result.mesh.y_coords is None`; selectbox `key="geo_quantity"` feeds `build_geometry_figure(result.mesh, quantity)`; `try/except RuntimeError` and module-attribute `petringa.run_field` reference preserved; registered in `app/main.py` navigation (not orphaned) |
| `tests/test_app_field_page.py`      | 2D-route test, selectbox-no-resolve test, updated fakes with populated mesh                          | VERIFIED | 6 tests, all pass; both fakes (`_fake_run_field`, `_fake_run_field_2d`) populate `mesh`; obsolete `test_2d_config_warns_and_skips` removed; call-counter no-resolve test present                                                                                                                                                                               |

### Key Link Verification

| From                 | To                                                     | Via                                                                      | Status | Details                                                                                                                   |
| -------------------- | ------------------------------------------------------ | ------------------------------------------------------------------------ | ------ | ------------------------------------------------------------------------------------------------------------------------- |
| `geometry_viewer.py` | `petringa.MeshData`                                    | `from petringa import MeshData` type-hint import                         | WIRED  | Present line 19; `MeshData` re-exported at top level per `petringa/__init__.py`                                           |
| `field_map.py`       | `app.components.geometry_viewer.build_geometry_figure` | import + `st.plotly_chart(build_geometry_figure(result.mesh, quantity))` | WIRED  | Present line 32 (import), line 91 (call site)                                                                             |
| `field_map.py`       | `st.session_state['field_result']`                     | cached read on rerun (no re-solve on selectbox change)                   | WIRED  | `test_selectbox_change_does_not_resolve` proves via call counter (n stays 1) + object-identity check on the cached result |

### Data-Flow Trace (Level 4)

| Artifact                            | Data Variable                                          | Source                                                                                                                                                          | Produces Real Data    | Status                                                                                                                                                                                                                                                                                                                                               |
| ----------------------------------- | ------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `build_geometry_figure` (2D branch) | `result.mesh`                                          | **Real, non-mocked** `petringa.run_field(DeviceConfig(half_width_um=50.0), bias_V=-20.0)` — actual devsim DD solve, converged (RelError ~1e-14 at iteration 12) | YES                   | FLOWING — `mesh.y_coords is not None` confirmed True; `build_geometry_figure(mesh, "Electric field")` → `data[0].type == "heatmap"`; all 3 quantities (Electric field / Net doping / Electrostatic potential) produced 20000/20000 finite z-values with physically sane ranges (E-field 0–72519 V/cm, doping log10 13.9–18.1, potential -1.6–21.3 V) |
| `build_geometry_figure` (1D branch) | `result.mesh`                                          | **Real, non-mocked** `petringa.run_field(DeviceConfig(), bias_V=-20.0)` — actual devsim 1D DD solve, converged                                                  | YES                   | FLOWING — `mesh.y_coords is None` confirmed True; `data[0].type == "bar"` for all 3 quantities                                                                                                                                                                                                                                                       |
| `field_map.py` render               | `result.mesh` (via `st.session_state["field_result"]`) | `petringa.run_field(cfg)` (module-attribute call site, unmocked in production)                                                                                  | YES (by construction) | Same production code path exercised above at the Python level; AppTest-level coverage uses mocked fakes (by design, to avoid slow/non-deterministic devsim in CI) but the underlying `run_field` -> `MeshData` -> `build_geometry_figure` chain is proven non-hollow by the direct real-solve check above                                            |

This closes the residual risk flagged in the phase's own `40-VERIFICATION.md` (plan-check) "Additional Note": that VIZ-01 was, at plan time, only exercisable via synthetic MeshData / mocked run_field. This goal-backward verification ran the actual production code path (`petringa.run_field` -> real devsim 2D solve -> `build_geometry_figure`) and confirmed it is NOT hollow — real numeric data flows all the way to a real Plotly heatmap trace.

### Behavioral Spot-Checks

| Behavior                                                                     | Command                                                                                                   | Result                                                                                                        | Status |
| ---------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ------ |
| Real 1D `run_field` solve feeds `build_geometry_figure` for all 3 quantities | `uv run python -c "...petringa.run_field(DeviceConfig(), bias_V=-20.0)..."`                               | Converged; `bar` trace for all 3 quantities, no exceptions                                                    | PASS   |
| Real 2D `run_field` solve feeds `build_geometry_figure` for all 3 quantities | `uv run python -c "...petringa.run_field(DeviceConfig(half_width_um=50.0), bias_V=-20.0)..."`             | Converged; `heatmap` trace, shape (100,200), 20000/20000 finite cells, sane value ranges for all 3 quantities | PASS   |
| Full per-file test isolation suite (documented convention, STATE.md:102)     | `uv run pytest tests/test_app_geometry_viewer.py tests/test_app_field_page.py tests/test_app_pages.py -v` | 18 passed, 0 failed                                                                                           | PASS   |
| No-devsim / no-Streamlit / no-fill_value / no-1D-only-string gates           | grep gates per plan acceptance criteria                                                                   | All return 0 as required                                                                                      | PASS   |
| Field-map page registered in app navigation (not orphaned)                   | `grep field_map app/main.py`                                                                              | `render_field_map` imported and registered at `url_path="field-map"`                                          | PASS   |

### Requirements Coverage

| Requirement | Source Plan  | Description                                                                                      | Status    | Evidence                                                                                                                                                                                                                                              |
| ----------- | ------------ | ------------------------------------------------------------------------------------------------ | --------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| VIZ-01      | 40-01, 40-02 | User sees 2D Plotly heatmap of E-field/doping after 2D field-map sim                             | SATISFIED | Real 2D devsim solve -> populated MeshData -> heatmap trace, confirmed at Python level; AppTest confirms page-level wiring (no st.stop, cached result); in-browser paint itself needs human confirmation (routine AppTest limitation, not a code gap) |
| VIZ-02      | 40-01, 40-02 | 1D devices show depth-profile bar instead of heatmap, same MeshData interface                    | SATISFIED | Real 1D devsim solve -> populated MeshData (y_coords=None) -> bar trace, confirmed at Python level; 1D branch in field_map.py keeps line charts and adds supplemental bar                                                                             |
| VIZ-03      | 40-01, 40-02 | Quantity dropdown (E-field/doping/potential) updates visualization without re-running simulation | SATISFIED | `QUANTITIES` dict complete and ordered correctly; call-counter AppTest proves no re-solve on selectbox change; visual confirmation of the actual dropdown-driven repaint needs human confirmation                                                     |

No orphaned requirements: REQUIREMENTS.md maps only VIZ-01/02/03 to Phase 40, and both plans declare exactly these three IDs.

### Anti-Patterns Found

None. Scanned `app/components/geometry_viewer.py`, `app/workflows/field_map.py`, `tests/test_app_geometry_viewer.py`, `tests/test_app_field_page.py` for TBD/FIXME/XXX/TODO/HACK/PLACEHOLDER/"not yet implemented"/empty-return stubs — zero matches. No hardcoded-empty stub patterns found; both fakes in the test file populate real mesh data (required by the plan and confirmed via `grep -c 'mesh=None'` → 0 for the fakes).

### Human Verification Required

### 1. Real 2D heatmap renders in the browser

**Test:** Run `streamlit run app/main.py`, set Dimensionality=2D (a `half_width_um` value) with a shallow bias in the sidebar, click "Run simulation" on the Field Map page.
**Expected:** A 2D Plotly heatmap of the electric field appears (device cross-section, "Lateral position (µm)" / "Depth (µm)" axes), not the old "1D-only" warning, not a blank page or crash.
**Why human:** Streamlit `AppTest` has no `plotly_chart` accessor, so actual browser rendering cannot be asserted by any automated test. This verifier independently ran the real (non-mocked) production code path — `petringa.run_field` with an actual 2D devsim solve feeding `build_geometry_figure` — and confirmed a real heatmap trace with correct shape and finite, sane data. The remaining gap is purely "does the browser paint it," which is inherent to the toolchain, not a code defect.

### 2. Quantity dropdown visibly changes the chart

**Test:** In the running app (either 1D or 2D device), after a Run, change the "Quantity" dropdown across all three options.
**Expected:** The chart's colorscale and value range visibly change (Viridis for E-field, Plasma for net doping, RdBu_r for potential) without any indication of a re-solve (e.g., long delay or console re-solve messages).
**Why human:** Visual appearance/colorscale differences are not asserted by AppTest; this is inherently a human-observable check, though the underlying scaling logic (`_scale_quantity`) is unit-tested for correctness at the data level.

### Gaps Summary

No code-level gaps. All artifacts exist, are substantive (not stubs), are wired correctly, and — going beyond the plan's own synthetic/mocked test posture — this verification additionally exercised the real production code path with actual (non-mocked) devsim solves for both 1D and 2D `DeviceConfig`s, confirming `MeshData` is populated correctly and `build_geometry_figure` renders real, finite, physically sane data for all three quantities in both dimensionalities. This closes the "unverified 2D convergence" concern noted in the phase's own pre-execution plan-check.

The only reason this phase is not marked `passed` is procedural: two behaviors (in-browser chart paint, in-browser dropdown-driven visual change) are only checkable by a human because Streamlit's `AppTest` harness has no chart-content accessor. This is the same posture the project accepted for Phase 39 (STATE.md:125) and does not indicate any missing or hollow implementation — it is a routine "needs eyes on a browser" item, not a gap requiring rework.

---

_Verified: 2026-07-13_
_Verifier: Claude (gsd-verifier)_
