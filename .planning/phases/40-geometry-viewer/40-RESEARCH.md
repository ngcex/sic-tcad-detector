# Phase 40: Geometry Viewer — Research

**Researched:** 2026-07-13
**Domain:** Streamlit + Plotly 2D visualization of irregular (triangular) devsim mesh node data via scipy griddata; 1D depth-profile bar chart. Extends the existing Phase 39 field-map page.
**Confidence:** HIGH (all contracts read verbatim from source; griddata→go.Heatmap pattern and AppTest selectbox accessor empirically verified this session)

## Summary

Phase 40 adds a **geometry viewer** to the existing field-map workflow. The heavy lifting is already done: `run_field(config)` (etna/api/simulation.py:144) populates a `MeshData` on `SimResult.mesh` for **both** 1D and 2D devices — irregular devsim node coordinates plus `node_values` for `NetDoping`, `Potential`, and `ElectricField` (already NaN-cleaned and node-aligned by the API layer). This phase is pure presentation: read `SimResult.mesh` from the cached result and render it — never call devsim.

The rendering approach is a **locked decision** (STATE.md:92): _"2D Plotly heatmap via `scipy.interpolate.griddata` onto regular grid; 1D bar for 1D devices."_ The existing matplotlib module `etna/core/plotting2d.py` is the authoritative **spec** for axis labels, colormaps, coordinate convention, and the log10 doping treatment — but it is matplotlib (`tricontourf`), not importable into `st.plotly_chart`, so the Plotly builder is written from scratch mirroring its labels/units. `griddata((x,y), z, (Xi,Yi)) → go.Heatmap(x=xi, y=yi, z=Zi)` is verified working this session (NaN outside the convex hull renders as transparent gaps; `update_yaxes(autorange="reversed")` reproduces matplotlib's `invert_yaxis()`).

**Three integration seams dominate the plan** and are NOT stated in the success criteria: **(1)** the current field-map page has a hard 2D guard (`if cfg.half_width_um is not None: st.warning(...); st.stop()`) that short-circuits 2D **before** `run_field` runs — Phase 40's job is to route 2D _through_, so that guard must be removed/replaced. **(2)** For 2D, `run_field` returns **empty `result.x`/`result.y`** and routes all data to `mesh`, so the existing `build_field_figures(result)` line charts operate on empty arrays for 2D and must be skipped on the 2D branch. **(3)** VIZ-03's "without re-running the simulation" is a **server-side Streamlit rerun** that reads the cached `SimResult.mesh` from `session_state` — the devsim solve does not re-run, but the griddata interpolation _does_ re-execute per dropdown change (fast, acceptable). It is NOT client-side.

**Primary recommendation:** Add a pure `app/components/geometry_viewer.py` with one builder `build_geometry_figure(mesh: MeshData, quantity: str) -> go.Figure` that branches on `mesh.y_coords is None` (1D→bar, 2D→griddata heatmap), applies quantity-appropriate scaling (E-field linear/viridis, potential diverging/RdBu, doping log10/plasma), and always converts its own coords cm→µm (`*1e4`). In `field_map.py`, replace the 2D `st.stop()` guard with a branch that calls `run_field` for both dims, renders the existing line charts only for 1D, and renders the geometry viewer (with an `st.selectbox` quantity dropdown feeding the builder) for both dims from `result.mesh`. Test the builder as a pure function against a hand-built synthetic `MeshData` (irregular points) — no devsim, no Streamlit — exactly like `to_csv_bytes`.

## Architectural Responsibility Map

| Capability                      | Primary Tier                      | Secondary Tier                | Rationale                                                             |
| ------------------------------- | --------------------------------- | ----------------------------- | --------------------------------------------------------------------- |
| Physics solve + mesh extraction | etna API (`run_field`)        | devsim                        | Already done; `mesh` populated post-build via `get_node_model_values` |
| Irregular→grid interpolation    | app UI (`geometry_viewer.py`)     | scipy `griddata`              | Locked decision (STATE.md:92); pure transform on `MeshData`           |
| Heatmap / bar construction      | app UI (`geometry_viewer.py`)     | Plotly `go.Figure`            | No Plotly code exists; `plotting2d.py` is matplotlib spec only        |
| Quantity selection              | Streamlit widget (`field_map.py`) | session_state (cached result) | `st.selectbox` triggers rerun; reads cached `mesh`, no re-solve       |
| 1D/2D dispatch                  | app UI (`build_geometry_figure`)  | —                             | Branch on `mesh.y_coords is None`, per VIZ-01/VIZ-02 same-interface   |
| Result caching                  | Streamlit session_state           | —                             | Inherited from Phase 39 (`field_result` key); DeviceConfig unhashable |

<user_constraints>

## User Constraints (from CONTEXT.md)

**No CONTEXT.md exists** for Phase 40 — discuss-phase was not run (same convention as Phase 39, "facciamo la plan direttamente"). The planner must make default decisions. This research supplies concrete `[ASSUMED]` recommendations (bar-chart vs line-chart for 1D, per-quantity scaling, NaN fill, dropdown default) that the planner should adopt as defaults; all such items are in the Assumptions Log for later user confirmation.

### Locked Decisions (carried from STATE.md — treat as CONTEXT.md Decisions)

- **Geometry viewer: 2D Plotly heatmap via `scipy.interpolate.griddata` onto a regular grid; 1D bar for 1D devices** (STATE.md:92). Do NOT re-litigate go.Contour / go.Scatter / tricontour alternatives — griddata→`go.Heatmap` is the chosen path.
- **`MeshData` is populated post-build from devsim node extraction; the geometry viewer never calls devsim directly** (STATE.md:91). The viewer reads only `SimResult.mesh`.
- **x=lateral, y=depth coordinate convention for all 2D modules** (STATE.md:79). Heatmap x-axis = lateral position, y-axis = depth (inverted, surface at top).
- **No physics changes allowed in any v5.0 phase — refactor/presentation only** (STATE.md:49).

### Claude's Discretion

- 1D visualization form (bar vs line), per-quantity color scaling, griddata resolution and NaN handling, dropdown default quantity and label strings, whether the viewer replaces or supplements the existing 1D E-field/potential line charts.

### Deferred Ideas (OUT OF SCOPE)

- Interactive lateral-slice extraction from a 2D mesh (picking a depth profile at a chosen lateral x). Not in VIZ-01/02/03.
- Region/contact overlay boxes on the heatmap (MeshData carries `regions`/`contacts`, but no requirement asks to draw them). Bonus-only if trivial.
- Any second devsim call or re-solve.
  </user_constraints>

<phase_requirements>

## Phase Requirements

| ID                    | Description                                                                                                              | Research Support                                                                                                                                                               |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| VIZ-01                | 2D Plotly heatmap of E-field (or doping) on device cross-section after a 2D field-map run (`half_width_um` set)          | `run_field` 2D path (simulation.py:216-293) populates `mesh.x_coords`/`y_coords` (cm) + `node_values`; griddata→`go.Heatmap` verified §3; existing 2D guard must be removed §2 |
| VIZ-02                | For 1D devices, a depth-profile **bar** chart instead of a heatmap, using the same `MeshData` interface                  | `mesh` is populated for 1D too (`y_coords=None`, all 3 node_values present, simulation.py:283-293); single builder branches on `mesh.y_coords is None` §4                      |
| VIZ-03                | Quantity selector dropdown (E-field, net doping, potential); visualization updates **without re-running** the simulation | `st.selectbox` rerun reads cached `SimResult.mesh` from `session_state["field_result"]`; devsim not re-run (griddata re-runs, fast) §5                                         |
| </phase_requirements> |

---

## 1. `MeshData` Contract (verbatim) — what the viewer consumes

`etna/api/results.py:20-33`:

```python
@dataclass
class MeshData:
    x_coords: np.ndarray            # node x coordinates (cm)   <-- CM, not µm
    y_coords: "np.ndarray | None"   # node y coordinates (cm), None for 1D
    node_values: dict[str, np.ndarray]  # "NetDoping", "Potential", "ElectricField"
    regions: list[dict]             # [{name, x_min, x_max, y_min, y_max}]  (bounds in cm)
    contacts: list[dict]            # [{name, position}]
```

`SimResult.mesh` is `MeshData | None`, populated by `run_field` for BOTH dims.

### What `run_field` actually puts in `mesh` (simulation.py:283-293, verbatim-verified)

For **both** 1D and 2D:

```python
mesh = MeshData(
    x_coords=x_coords,          # cm (node "x" from get_node_model_values)
    y_coords=y_coords,          # cm for 2D, None for 1D
    node_values={
        "NetDoping": net_doping,     # cm^-3, node array (can be signed: p+ negative, n- positive)
        "Potential": potential,      # V, node array
        "ElectricField": field_nodes # V/cm magnitude, node-aligned, NaN already replaced with 0.0
    },
    regions=[{name, x_min, x_max, y_min, y_max}],   # y_min/y_max None for 1D
    contacts=[{"name":"anode","position":0.0}, {"name":"cathode","position":total_length}],
)
```

**Key facts the builder must honor:**

- `x_coords` / `y_coords` are in **cm** → multiply by `1e4` for µm axis labels. (Distinct from `result.x`, which is already µm — the viewer uses `mesh`, so it owns the cm→µm conversion.)
- `ElectricField` is already a **node-aligned magnitude with NaN→0.0** done in the API layer (simulation.py:259-262). No further cleaning needed for the source values; griddata gaps outside the convex hull are a separate concern (§3).
- The three `node_values` keys are exactly `"NetDoping"`, `"Potential"`, `"ElectricField"` — the dropdown maps friendly labels to these keys.
- 2D convergence for `run_field` is genuinely **unverified** in this session (see §Environment) — do NOT run a real 2D solve to test. Test the builder against a synthetic `MeshData`.

## 2. The existing field-map page — what must change (integration seam #1)

`app/workflows/field_map.py` (read this session) currently:

```python
if cfg.half_width_um is not None:
    st.warning("These workflows are 1D-only. 2D field visualization arrives in "
               "Phase 40 (geometry viewer). Set Dimensionality to 1D ...")
    st.stop()                              # <-- BLOCKS 2D before run_field
...
if st.button("Run simulation"):
    try:
        st.session_state["field_result"] = etna.run_field(cfg)
    except RuntimeError as e:
        st.error(...)                      # <-- keep: 2D ramp may not converge
result = st.session_state.get("field_result")
if result is not None:
    efield_fig, potential_fig = build_field_figures(result)   # <-- EMPTY arrays for 2D
    st.plotly_chart(efield_fig); st.plotly_chart(potential_fig)
    st.download_button("Download CSV", data=to_csv_bytes(result), ...)
    with st.expander("Net doping vs depth"): ...
```

**Required changes:**

1. **Remove the `st.stop()` 2D guard.** Phase 40's entire purpose is to route 2D through. (The warning string was even authored as a placeholder for "Phase 40" — this is that phase.)
2. **Keep the `try/except RuntimeError`** around `run_field` — 2D `ramp_bias` may not converge (STATE.md blockers document `ramp_bias` non-convergence at deep bias). The viewer must render only when `result.mesh is not None`.
3. **Branch the render on dimensionality.** `build_field_figures(result)` reads `result.x`/`result.y`/`metadata["potential"]` which are **empty for 2D** (simulation.py:301-311 returns `x=np.array([])`, `y=np.array([])`). So:
   - **1D branch** (`result.mesh.y_coords is None`): existing line charts still valid (they use `result.x`/`result.y`, populated for 1D). Add the geometry-viewer bar chart + dropdown per VIZ-02/VIZ-03. **Decision needed** (A1): does the bar chart _replace_ or _supplement_ the existing E-field/potential line charts? Recommend **supplement** (line charts are richer; bar satisfies the literal VIZ-02 wording). Planner to confirm.
   - **2D branch**: **skip** `build_field_figures` (empty arrays), render only the geometry heatmap + dropdown. CSV download for 2D is out of scope for this phase (no VIZ requirement asks for it; `to_csv_bytes` has no 2D branch and would `raise ValueError` — guard it or skip the download button for 2D). **Decision needed** (A2).

## 3. Rendering pattern — griddata → go.Heatmap (2D) — VERIFIED

Locked path (STATE.md:92). Verified working this session:

```python
# Source: verified locally 2026-07-13 (scipy 1.18.0, plotly 6.8.0)
import numpy as np
from scipy.interpolate import griddata
import plotly.graph_objects as go

# mesh.x_coords / mesh.y_coords are CM -> µm
x_um = mesh.x_coords * 1e4          # lateral (µm)
y_um = mesh.y_coords * 1e4          # depth (µm)
z = mesh.node_values[key]           # "ElectricField" | "Potential" | "NetDoping"

# regular grid spanning the mesh bounding box
xi = np.linspace(x_um.min(), x_um.max(), 200)   # n_x
yi = np.linspace(y_um.min(), y_um.max(), 100)   # n_y (fewer: epi is thin)
Xi, Yi = np.meshgrid(xi, yi)
Zi = griddata((x_um, y_um), z, (Xi, Yi), method="linear")   # NaN outside convex hull

fig = go.Figure(data=go.Heatmap(x=xi, y=yi, z=Zi, colorbar_title=..., colorscale=...))
fig.update_layout(title=..., xaxis_title="Lateral position (µm)", yaxis_title="Depth (µm)")
fig.update_yaxes(autorange="reversed")   # matches plotting2d.py invert_yaxis() — surface at top
```

- **NaN outside the convex hull:** Plotly `go.Heatmap` renders `NaN` cells as **transparent gaps** (verified — no error). For a near-rectangular device mesh the hull ≈ the bounding box so gaps are minor. **Decision (A3):** leave as NaN (transparent) OR pass `fill_value=0.0` to `griddata` (fills exterior with 0). Recommend **NaN/transparent** — a `0` fill would paint fake zero-field regions. Planner to confirm.
- The API's internal 2D E-field computation uses the _same_ `LinearNDInterpolator` grid pattern at n_x=100, n_y=200 (simulation.py:241-259). `griddata(method="linear")` is the functionally-equivalent one-call form of `LinearNDInterpolator`; either is acceptable. `griddata` is the STATE.md-named choice — use it.

### 1D bar chart (VIZ-02)

```python
# 1D: mesh.y_coords is None; x_coords is depth (cm) for the 1D device
depth_um = mesh.x_coords * 1e4
z = mesh.node_values[key]
fig = go.Figure(data=go.Bar(x=depth_um, y=z))   # bar per literal VIZ-02 wording
fig.update_layout(xaxis_title="Depth (µm)", yaxis_title=<quantity label+units>)
```

Note: a 1D devsim mesh may have many closely-spaced nodes → a dense bar chart looks near-solid. That is acceptable and matches the requirement's face-value "bar chart" wording (VIZ-02 + ROADMAP criterion 2). Do not silently substitute a line chart for the _geometry viewer's_ 1D output (though the page may still also show the existing line charts — see §2 A1).

## 4. Single builder, branch on `mesh.y_coords is None` (VIZ-01 + VIZ-02 "same interface")

VIZ-01/VIZ-02 require the 2D heatmap and 1D bar to use the **same `MeshData` interface**. Implement one pure function:

```python
# app/components/geometry_viewer.py  (PURE — no st.* calls)
QUANTITIES = {                       # friendly label -> node_values key
    "Electric field": "ElectricField",
    "Net doping": "NetDoping",
    "Electrostatic potential": "Potential",
}

def build_geometry_figure(mesh: MeshData, quantity: str) -> go.Figure:
    key = QUANTITIES[quantity]
    z = mesh.node_values[key]
    z, colorbar_label, colorscale = _scale_quantity(key, z)   # see per-quantity scaling
    if mesh.y_coords is None:
        return _build_bar(mesh, z, quantity, colorbar_label)   # 1D
    return _build_heatmap(mesh, z, quantity, colorbar_label, colorscale)  # 2D
```

Keep it pure so a unit test builds a synthetic `MeshData` and asserts on the returned `go.Figure` (trace type `Heatmap` vs `Bar`, axis titles, z shape) with **no devsim and no Streamlit**.

### Per-quantity scaling (mirror plotting2d.py — `[ASSUMED]` A4)

| node_values key | plotting2d.py precedent                         | Scaling                            | Colorscale | Colorbar label             |
| --------------- | ----------------------------------------------- | ---------------------------------- | ---------- | -------------------------- |
| `ElectricField` | plot_efield_2d (viridis)                        | linear                             | Viridis    | \|Electric Field\| (V/cm)  |
| `Potential`     | plot_potential_2d (RdBu_r)                      | linear (diverging)                 | RdBu_r     | Potential (V)              |
| `NetDoping`     | plot_doping_2d — `log10(\|NetDoping\|)`, plasma | `np.log10(np.abs(z))` (guard z==0) | Plasma     | log₁₀\|Net Doping\| (cm⁻³) |

`NetDoping` can be signed and spans many orders of magnitude → the log10(|·|) treatment (from plotting2d.py:255+ `plot_doping_2d`) is required for a readable map. Guard `z==0` before log10 (`np.where(z==0, tiny, z)` or mask). E-field and potential plot raw.

## 5. Quantity dropdown + "no re-run" mechanic (VIZ-03) — integration seam #3

**Accurate mechanic:** `st.selectbox` change → **server-side Streamlit rerun** → the page reads the already-cached `SimResult` from `st.session_state["field_result"]` (the `run_field` devsim solve is NOT re-executed because the Run button returns `False` on this rerun) → `build_geometry_figure(result.mesh, selected_quantity)` re-runs the griddata interpolation (fast, ~ms) → new heatmap. This is a server-side re-render of cached mesh data, **not** a client-side update.

```python
# in field_map.py, after result is retrieved from session_state and result.mesh is not None:
from app.components.geometry_viewer import build_geometry_figure, QUANTITIES
quantity = st.selectbox("Quantity", list(QUANTITIES.keys()))   # default = first ("Electric field")
st.plotly_chart(build_geometry_figure(result.mesh, quantity))
```

The `st.selectbox` widget key persists across reruns automatically; no manual session_state bookkeeping needed for the selection. The cached `field_result` (already established by Phase 39) is what makes "without re-running" true — do **not** add `@st.cache_data` (DeviceConfig is unhashable — verified in Phase 39; STATE.md caching decisions).

## 6. Coordinate & label spec (from plotting2d.py, authoritative)

| Element                | Value                                                       | Source                                                                               |
| ---------------------- | ----------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| x-axis (2D)            | "Lateral position (µm)"                                     | plotting2d.py:140,229                                                                |
| y-axis (2D)            | "Depth (µm)", **inverted** (surface at top)                 | plotting2d.py:141-143 `invert_yaxis()` → Plotly `update_yaxes(autorange="reversed")` |
| x-axis (1D bar)        | "Depth (µm)"                                                | 1D x_coords is depth                                                                 |
| Potential title / cbar | "2D Potential Map" / "Potential (V)"                        | plotting2d.py:138,142                                                                |
| E-field title / cbar   | "2D Electric Field Magnitude" / "\|Electric Field\| (V/cm)" | plotting2d.py:227,231                                                                |
| Doping                 | `log10(\|NetDoping\|)`, plasma                              | plotting2d.py:255+                                                                   |

Coordinate convention **x=lateral, y=depth** is locked (STATE.md:79). Mesh `x_coords`=lateral, `y_coords`=depth for 2D.

---

## Architecture Patterns

### Data-flow diagram

```
DeviceConfig (sidebar, session_state["device_config"])
      │
      ▼
[Run simulation] button ──► etna.run_field(cfg)  ──► devsim solve (ONCE)
      │                            │
      │                            ▼
      │                     SimResult (cached in session_state["field_result"])
      │                       ├─ x, y, metadata   (populated for 1D; EMPTY for 2D)
      │                       └─ mesh: MeshData    (populated for BOTH dims)
      │                              ├─ x_coords, y_coords (cm)
      │                              └─ node_values{ElectricField, Potential, NetDoping}
      ▼
  every rerun (incl. dropdown change — NO re-solve):
      result = session_state["field_result"]
      │
      ├── 1D (mesh.y_coords is None) ──► build_field_figures (line, from result.x/y)  [existing]
      │                                └► st.selectbox ─► build_geometry_figure(mesh, q) ─► go.Bar
      │
      └── 2D (mesh.y_coords not None) ─► skip line charts (result.x/y empty)
                                        └► st.selectbox ─► build_geometry_figure(mesh, q)
                                                           └► griddata → go.Heatmap
```

### Recommended file layout (mirror Phase 39 split)

```
app/
├── components/
│   ├── results.py            # existing: build_field_figures, to_csv_bytes (PURE)
│   └── geometry_viewer.py    # NEW: build_geometry_figure(mesh, quantity) + QUANTITIES (PURE)
└── workflows/
    └── field_map.py          # MODIFIED: remove 2D stop-guard, branch render, add selectbox
```

Keep `geometry_viewer.py` **pure** (no `st.*`) so it is unit-testable exactly like `results.py`. The `st.selectbox` / `st.plotly_chart` calls live in `field_map.py`.

### Anti-patterns to avoid

- **Calling devsim from the viewer.** Forbidden (STATE.md:91). Read only `result.mesh`.
- **Re-running `run_field` on dropdown change.** Read the cached result; the builder re-runs, the solve does not.
- **Using `result.x`/`result.y` for 2D.** They are empty arrays for 2D (simulation.py:301-311). Use `result.mesh`.
- **Mixing cm and µm.** `mesh.x_coords`/`y_coords` are **cm**; convert `*1e4` in the viewer. (Contrast: `result.x` is already µm — but the viewer must not use `result.x`.)
- **Importing `plotting2d.py`.** It is matplotlib (`tricontourf`), not usable in `st.plotly_chart`. Use it as label/colormap spec only.
- **`fill_value=0` on doping/potential griddata.** Would paint fake zeros; prefer NaN-transparent (A3).

## Don't Hand-Roll

| Problem                         | Don't Build                     | Use Instead                                                      | Why                                                          |
| ------------------------------- | ------------------------------- | ---------------------------------------------------------------- | ------------------------------------------------------------ |
| Irregular→grid interpolation    | Manual barycentric / Delaunay   | `scipy.interpolate.griddata((x,y), z, (Xi,Yi), method="linear")` | Locked (STATE.md:92); handles unstructured mesh; verified §3 |
| 2D heatmap rendering            | Custom canvas / manual colormap | `plotly.graph_objects.Heatmap`                                   | Renders NaN transparent, has colorbar, interactive zoom free |
| y-axis inversion                | Manual coordinate flip          | `fig.update_yaxes(autorange="reversed")`                         | One line; matches plotting2d.py `invert_yaxis()`             |
| Result caching                  | `@st.cache_data` / custom memo  | `st.session_state["field_result"]` (already exists)              | DeviceConfig unhashable (Phase 39 verified)                  |
| E-field magnitude / NaN cleanup | Recompute gradient in the page  | `mesh.node_values["ElectricField"]`                              | API already computed + NaN→0 (simulation.py:259-262)         |

**Key insight:** Everything physics/numeric (solve, E-field magnitude, NaN cleanup, node alignment) is done in the API layer. The viewer is a pure `MeshData → go.Figure` transform plus a Streamlit widget. Any recompute in the page is a bug.

## Common Pitfalls

### Pitfall 1: Leaving the 2D `st.stop()` guard in place (HIGHEST PRIORITY)

**What goes wrong:** `field_map.py` currently `st.stop()`s for any 2D config before `run_field`. If not removed, VIZ-01 is unreachable — the 2D heatmap can never render.
**How to avoid:** Remove the `if cfg.half_width_um is not None: st.warning(...); st.stop()` block. Route 2D through `run_field` (keep the `try/except RuntimeError`).
**Warning sign:** Selecting Dimensionality=2D shows the old "1D-only" warning instead of a heatmap.

### Pitfall 2: Feeding empty 2D arrays into the line-chart builder

**What goes wrong:** For 2D, `result.x`/`result.y`/`metadata["potential"]` are empty (or `metadata["potential"]` is a full node array but not a depth profile). `build_field_figures(result)` produces empty/meaningless line charts, and `metadata`-key access may `KeyError` if the 2D metadata differs.
**How to avoid:** Branch on `result.mesh.y_coords is None`. Call `build_field_figures` and `to_csv_bytes` only on the 1D branch. On the 2D branch, render only the geometry heatmap.
**Warning sign:** Empty/flat line charts appear above the 2D heatmap, or a `KeyError`/`ValueError` on a 2D run.

### Pitfall 3: cm vs µm axis mislabel

**What goes wrong:** `mesh.x_coords`/`y_coords` are **cm**; plotting them directly labels the axis 10⁴× too small. (Phase 39's line charts avoided this by using `result.x`, already µm — but the viewer uses `mesh`.)
**How to avoid:** In `geometry_viewer.py`, always `* 1e4` the mesh coordinates before building axes.
**Warning sign:** Depth axis reads 0–0.001 instead of 0–10 µm.

### Pitfall 4: log10 on signed / zero NetDoping

**What goes wrong:** `NetDoping` is signed (p+ vs n-) and can be 0 at the junction → `log10` gives NaN/-inf, and a raw signed doping heatmap is unreadable across orders of magnitude.
**How to avoid:** For the doping quantity, plot `np.log10(np.abs(z))` with a small floor for z==0 (mirror plotting2d.py `plot_doping_2d`). Only doping needs log; E-field and potential are linear.
**Warning sign:** Doping heatmap is uniform/blank or throws on log of zero.

### Pitfall 5: Testing against a real 2D devsim solve

**What goes wrong:** 2D `ramp_bias` convergence is unverified and devsim resource-exhaustion is a documented CI blocker (STATE.md). A test that runs `run_field` on a 2D config may hang, crash the pool, or hit non-convergence.
**How to avoid:** Test `build_geometry_figure` against a **hand-built synthetic `MeshData`** (irregular `x_coords`/`y_coords`, fabricated `node_values`) — pure, no devsim, no Streamlit. Same pattern as `tests/test_app_csv_export.py`.
**Warning sign:** A viewer test imports/calls `etna.run_field` or `devsim`.

## Runtime State Inventory

Not applicable — Phase 40 is a greenfield presentation component (one new pure module + one page edit), no rename/refactor/migration of stored data, services, OS state, secrets, or build artifacts.

## Environment Availability

| Dependency | Required By                        | Available | Version     | Fallback                                             |
| ---------- | ---------------------------------- | --------- | ----------- | ---------------------------------------------------- |
| plotly     | heatmap/bar (VIZ-01/02)            | ✓         | **6.8.0**   | — (Phase 39 installed it; **note: not 5.x**)         |
| scipy      | griddata interpolation (VIZ-01)    | ✓         | 1.18.0      | —                                                    |
| numpy      | arrays                             | ✓         | (installed) | —                                                    |
| streamlit  | selectbox + page                   | ✓         | 1.55.0      | —                                                    |
| devsim     | `run_field` solve (2D convergence) | ✓ (loads) | —           | **2D convergence UNVERIFIED — open risk, see below** |

**Missing dependencies with no fallback:** none.

**Open risk (not a missing dependency):** 2D `run_field` convergence was **not** exercised this session (devsim resource-exhaustion + documented `ramp_bias` non-convergence at deep bias — STATE.md blockers). The page already wraps `run_field` in `try/except RuntimeError` and the viewer should render only when `result.mesh is not None`. VIZ-01 is code-complete-able and unit-testable against synthetic MeshData regardless, but end-to-end 2D browser verification may surface a solver non-convergence that is **upstream** (same class as the Phase 39 `ramp_bias` blocker), not a Phase 40 UI defect. Flag this to the planner: acceptance for VIZ-01 should be the builder + wiring + synthetic-mesh test, with browser 2D verification best-effort (a non-converging 2D solve is out-of-scope physics work).

## Validation Architecture

`workflow.nyquist_validation` is **absent** in config.json → treat as **enabled**.

### Test Framework

| Property          | Value                                                                   |
| ----------------- | ----------------------------------------------------------------------- |
| Framework         | pytest ≥7.0 (dev extra) + Streamlit `AppTest` (streamlit.testing.v1)    |
| Config file       | `pyproject.toml` `[tool.pytest.ini_options]`, `pythonpath=["."]`        |
| Quick run         | `uv run pytest tests/test_app_geometry_viewer.py -x`                    |
| Full suite (page) | `uv run pytest tests/test_app_field_page.py tests/test_app_pages.py -x` |

> Use `uv run pytest` (STATE.md:106 — bare `pytest` may miss venv packages). Per-file isolation is mandatory for any devsim-touching test (STATE.md:102,124) — but the viewer tests touch **no** devsim.

### Phase Requirements → Test Map

| Req      | Behavior                                                                      | Test Type                                                         | Command                                              | Exists?            |
| -------- | ----------------------------------------------------------------------------- | ----------------------------------------------------------------- | ---------------------------------------------------- | ------------------ |
| VIZ-01   | 2D MeshData → `go.Heatmap` (trace type, z shape, axis titles)                 | Pure unit on `build_geometry_figure` (synthetic 2D MeshData)      | `uv run pytest tests/test_app_geometry_viewer.py -x` | ❌ Wave 0          |
| VIZ-02   | 1D MeshData (`y_coords=None`) → `go.Bar`, same builder                        | Pure unit (synthetic 1D MeshData)                                 | (same file)                                          | ❌ Wave 0          |
| VIZ-03   | selectbox present; changing it re-renders from cached result without re-solve | AppTest (`at.selectbox`) + monkeypatch `etna.run_field`       | `uv run pytest tests/test_app_field_page.py -x`      | ⚠️ extend existing |
| 2D route | 2D config no longer hits `st.stop()`; heatmap path reached                    | AppTest, monkeypatched `run_field` returns synthetic 2D SimResult | (in field page test)                                 | ⚠️ extend existing |

**AppTest accessor facts (verified this session, streamlit 1.55):**

- `at.selectbox` **exists** (verified — dropdown interaction IS testable: set `at.selectbox[0].select("Net doping").run()` and assert `at.exception == []`).
- There is still **NO** `plotly_chart` and **NO** `download_button` accessor. So AppTest asserts on `at.exception == []`, `at.selectbox` label/options/value, `at.warning`/`at.info`, and `at.session_state["field_result"]` — **not** on the chart itself.
- Chart-content assertions (heatmap vs bar, axis titles, z shape) go in the **pure unit test** on `build_geometry_figure`, which needs neither Streamlit nor devsim.

**Mocking:** Continue referencing the facade as `etna.run_field` (module attribute) so tests `monkeypatch.setattr(etna, "run_field", fake)` returning a synthetic `SimResult` with a populated `mesh` (STATE.md:103; proven in `tests/test_app_field_page.py`). For a 2D fake, set `x=np.array([])`, `y=np.array([])`, and `mesh.y_coords` populated.

### Wave 0 Gaps

- [ ] `app/components/geometry_viewer.py` — `build_geometry_figure(mesh, quantity)` + `QUANTITIES` (PURE)
- [ ] `tests/test_app_geometry_viewer.py` — pure unit tests: synthetic 1D MeshData → Bar, synthetic 2D MeshData → Heatmap, per-quantity scaling, cm→µm, doping log10
- [ ] Extend `tests/test_app_field_page.py` — 2D route (no `st.stop`), selectbox present + change re-renders from cache (mock `run_field` → synthetic 2D SimResult)
- [ ] Remove obsolete 2D-guard assertion from the existing field-page test if it asserts the "1D-only" warning for a 2D config

_(No framework install needed — pytest + AppTest already present.)_

## Security Domain

`security_enforcement` **absent** in config.json → enabled. Phase 40 has **no new input surface**: inputs are numeric sidebar widgets (bounded) plus a fixed-option `st.selectbox` (values hardcoded in `QUANTITIES`). No file upload, no free-text, no CSV write introduced by this phase (2D CSV export is out of scope). Data rendered is numeric arrays from a trusted internal solve.

| ASVS                | Applies | Control                                                                                    |
| ------------------- | ------- | ------------------------------------------------------------------------------------------ |
| V5 Input Validation | partial | Dropdown options are a fixed hardcoded set; numeric sidebar widgets bounded by `min_value` |
| V6 Cryptography     | no      | —                                                                                          |

**Threat patterns:** none material — no tainted data reaches a sink (no eval, no file path from user, no template). Risk: LOW.

## Assumptions Log

| #   | Claim                                                                                             | Section | Risk if Wrong                                                                   |
| --- | ------------------------------------------------------------------------------------------------- | ------- | ------------------------------------------------------------------------------- |
| A1  | 1D: geometry-viewer bar **supplements** (not replaces) the existing E-field/potential line charts | §2      | Low — UX layout; both satisfy VIZ-02, planner may prefer replace                |
| A2  | 2D branch **skips** CSV download (out of scope; `to_csv_bytes` has no 2D branch)                  | §2      | Low — no VIZ req asks for 2D CSV; if wanted, needs a new to_csv_bytes 2D branch |
| A3  | griddata leaves NaN outside convex hull (transparent gaps), not `fill_value=0`                    | §3      | Low — 0-fill would paint fake zero-field; NaN is safer/clearer                  |
| A4  | Per-quantity scaling: E-field linear/viridis, potential linear/RdBu_r, doping log10/plasma        | §4      | Low — mirrors plotting2d.py; doping-log10 is load-bearing for readability       |
| A5  | griddata grid resolution n_x≈200, n_y≈100 (lateral dense, depth thinner)                          | §3      | Low — cosmetic sharpness; matches API's own 100×200 order-of-magnitude          |
| A6  | Dropdown default = "Electric field" (first option)                                                | §5      | Low — cosmetic default                                                          |
| A7  | Coordinate mapping: 2D mesh x=lateral, y=depth (per STATE.md:79); heatmap x=lateral               | §6      | Low — locked convention; verify once in browser if 2D solve converges           |

## Open Questions (RESOLVED)

1. **Does 2D `run_field` converge for a realistic 2D `DeviceConfig`?**
   - What we know: 1D `ramp_bias` fails at deep bias (STATE.md); 2D path is untested this session; devsim CI resource-exhaustion blocks a quick check here.
   - What's unclear: whether a 2D solve completes at the default/shallow bias to actually populate a non-empty `mesh` in the browser.
   - Recommendation: build + unit-test the viewer against synthetic MeshData (fully deterministic); treat browser 2D verification as best-effort; a non-convergence is upstream physics, not a Phase 40 defect (same posture Phase 39 took for `ramp_bias`).

2. **Should regions/contacts be overlaid on the heatmap?**
   - What we know: `MeshData.regions`/`contacts` carry bounding boxes / positions; no VIZ requirement asks to draw them.
   - Recommendation: out of scope (Deferred). Add as a trivial bonus (`fig.add_shape` rectangle / `add_hline` at contact depth) only if planner wants it; keep default off.

## Sources

### Primary (HIGH — read verbatim this session)

- `etna/api/results.py:20-46` — `MeshData`, `SimResult`
- `etna/api/simulation.py:144-317` — `run_field` full body (1D + 2D mesh population, empty-array 2D return, E-field magnitude + NaN cleanup)
- `etna/core/plotting2d.py:106-269` — matplotlib label/colormap/coordinate spec (`plot_potential_2d`, `plot_efield_2d`, `plot_doping_2d`; `invert_yaxis`, log10 doping)
- `app/workflows/field_map.py` (full) — existing page, 2D `st.stop()` guard to remove, `try/except`, render flow
- `app/components/results.py` (full) — `build_field_figures`, `to_csv_bytes` (no 2D branch)
- `.planning/STATE.md` — locked decisions (griddata:92, no-devsim:91, x=lateral:79, no-physics:49, unhashable-config/caching, `run_field` referenced as module attr:103, `uv run pytest`:106)
- `.planning/ROADMAP.md:491-503`, `.planning/REQUIREMENTS.md:36-40` — Phase 40 success criteria + VIZ-01/02/03
- `.planning/phases/39-*/39-RESEARCH.md`, `39-PATTERNS.md`, `39-04-SUMMARY.md` — inherited Streamlit/Plotly/AppTest patterns

### Empirically verified this session

- `griddata((x,y),z,(Xi,Yi),"linear")` → `go.Heatmap(x=xi,y=yi,z=Zi)` renders; NaN outside hull = transparent; `update_yaxes(autorange="reversed")` works
- `plotly==6.8.0` (NOT 5.x), `scipy==1.18.0` installed in venv
- streamlit 1.55 AppTest: `at.selectbox` **exists**; still no `plotly_chart` / `download_button` accessor

## Metadata

**Confidence breakdown:**

- MeshData contract + `run_field` 2D behavior: HIGH — read verbatim from source
- Rendering pattern (griddata→Heatmap, bar, yaxis reverse): HIGH — verified locally + locked decision
- Integration seams (remove 2D guard, skip line charts for 2D, no-rerun mechanic): HIGH — grounded in read source
- Per-quantity scaling / 1D bar-vs-line / NaN fill / grid resolution: MEDIUM (`[ASSUMED]`) — mirror plotting2d.py, no CONTEXT.md; concrete defaults recommended
- 2D solver convergence end-to-end: LOW — unverified, documented upstream risk

**Research date:** 2026-07-13
**Valid until:** ~30 days (stable local codebase; risks are upstream Streamlit accessor changes or 2D solver work)
