---
phase: 40-geometry-viewer
plan: 01
subsystem: app-ui
tags: [plotly, scipy, griddata, geometry-viewer, pure-builder, viz]
requires:
  - petringa.MeshData (top-level re-export of petringa.api.results.MeshData)
  - scipy.interpolate.griddata
  - plotly.graph_objects
provides:
  - "app/components/geometry_viewer.py: build_geometry_figure(mesh, quantity) -> go.Figure"
  - "app/components/geometry_viewer.py: QUANTITIES label->key map (Electric field first)"
affects:
  - 40-02 (field_map page wiring will import build_geometry_figure + QUANTITIES)
tech-stack:
  added:
    - "scipy.interpolate.griddata (first use in app/ layer)"
    - "plotly go.Heatmap + go.Bar (first heatmap/bar in repo)"
  patterns:
    - "pure MeshData -> go.Figure builder (mirrors results.py purity discipline)"
    - "branch on mesh.y_coords is None: 1D bar vs 2D griddata heatmap"
key-files:
  created:
    - app/components/geometry_viewer.py
    - tests/test_app_geometry_viewer.py
  modified: []
decisions:
  - "A3 NaN-transparent: no fill_value on griddata (avoids painting fake zeros)"
  - "A4 per-quantity scaling: E-field linear/Viridis, potential linear/RdBu_r, doping log10(|.|)/Plasma"
  - "A5 grid resolution n_x=200, n_y=100"
  - "A6 dropdown default = Electric field (first key in QUANTITIES)"
  - "A7 x=lateral, y=depth; 2D y-axis autorange=reversed (surface at top)"
metrics:
  duration: ~15min
  tasks: 2
  files: 2
  completed: 2026-07-13
---

# Phase 40 Plan 01: Geometry Viewer Builder Summary

Pure `MeshData -> go.Figure` builder (`build_geometry_figure` + `QUANTITIES`) that renders a scipy-griddata Plotly heatmap for 2D meshes and a depth-profile bar for 1D meshes from the same interface, with per-quantity scaling and cm→µm conversion, plus its pure unit-test suite. No Streamlit, no devsim.

## What Was Built

- **`app/components/geometry_viewer.py`** (pure, no Streamlit, reads only `MeshData`):
  - `QUANTITIES` maps the three dropdown labels to `node_values` keys in insertion order, with `"Electric field"` first (dropdown default): `{"Electric field": "ElectricField", "Net doping": "NetDoping", "Electrostatic potential": "Potential"}`.
  - `build_geometry_figure(mesh, quantity)` looks up the key, applies `_scale_quantity`, then branches on `mesh.y_coords is None`.
  - `_scale_quantity`: E-field linear/Viridis, Potential linear/RdBu_r, NetDoping `log10(|z|)` with a `|z| < 1.0 → 1.0` floor / Plasma (mirrors `plotting2d.py` lines 259-261).
  - `_build_bar` (1D): `go.Bar` with x = `x_coords * 1e4` (µm depth), axis title "Depth (µm)".
  - `_build_heatmap` (2D): regular grid `xi[200] × yi[100]`, `griddata(..., method="linear")` (no `fill_value` — NaN outside the convex hull renders transparent), `go.Heatmap`, "Lateral position (µm)" / "Depth (µm)" titles, `update_yaxes(autorange="reversed")`.
- **`tests/test_app_geometry_viewer.py`** — 8 pure unit tests over synthetic `_mesh_1d()` / `_mesh_2d()` (irregular 5-node scatter) fixtures: 1D→Bar, 2D→Heatmap with z shape (100, 200), reversed depth axis + µm axis titles, cm→µm bar x-values, doping log10 on signed/zero without raising (1D + 2D), and the `QUANTITIES` contract.

## Verification

- `uv run python -c "from app.components.geometry_viewer import build_geometry_figure, QUANTITIES"` → import OK; `QUANTITIES` order/mapping asserted.
- `uv run pytest tests/test_app_geometry_viewer.py -x` → 8 passed.
- Acceptance greps: `griddata`/`MeshData` imports present; `method="linear"` present; `fill_value` count 0; `autorange="reversed"` present; `* 1e4` count 3 (≥2); no-devsim gate 0; no-`st.` gate 0; test file `run_field`/`devsim` count 0; heatmap + bar assertions present.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] `uv run pytest` fell back to a non-project pytest (ModuleNotFoundError: plotly)**

- **Found during:** Task 2 verification.
- **Issue:** Bare `uv run pytest` resolved a non-venv `pytest` on PATH that could not see project packages (plotly), causing collection to fail — the exact pitfall documented in STATE.md decision line 106.
- **Fix:** Ran `uv sync --extra dev` to materialize the `dev` optional-dependency group (pytest) into the project venv; after that `uv run pytest tests/test_app_geometry_viewer.py -x` (the plan's exact verify command) passes. No source change required; this is an environment materialization step, not a code fix.
- **Files modified:** none (source); venv/lock only.
- **Commit:** n/a (env step, no tracked source change).

### Docstring wording (grep-gate safe)

Per the plan's own guidance and to satisfy the `grep -c 'st\.'` == 0 and no-devsim gates, the module docstring/comments avoid the literal `st.` token (phrased "no Streamlit calls") and avoid a `devsim`-adjacent period ("reads only MeshData"). This is compliance with the stated acceptance gates, not a behavioral deviation.

## Known Stubs

None — the builder is fully wired against the real `MeshData` interface. (Page wiring / `st.selectbox` / `st.plotly_chart` are intentionally out of scope for this interface-first plan and land in 40-02.)

## Self-Check: PASSED

- FOUND: app/components/geometry_viewer.py
- FOUND: tests/test_app_geometry_viewer.py
- FOUND commit c9a8cee (feat: builder)
- FOUND commit a751389 (test: unit tests)
