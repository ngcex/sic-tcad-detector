# Phase 40: Geometry Viewer - Pattern Map

**Mapped:** 2026-07-13
**Files analyzed:** 4 (2 new, 2 modified)
**Analogs found:** 4 / 4 (structure); heatmap/griddata body has no codebase analog (see No Analog Found)

## File Classification

| New/Modified File                         | Role                            | Data Flow                        | Closest Analog                                                          | Match Quality                                                             |
| ----------------------------------------- | ------------------------------- | -------------------------------- | ----------------------------------------------------------------------- | ------------------------------------------------------------------------- |
| `app/components/geometry_viewer.py` (NEW) | component (pure figure builder) | transform (MeshData → go.Figure) | `app/components/results.py`                                             | role-match (structure/purity/signature exact; heatmap body has NO analog) |
| `app/workflows/field_map.py` (MODIFIED)   | workflow (Streamlit page)       | request-response                 | itself (edit in place) + `app/components/device_sidebar.py` (selectbox) | exact (self-edit); selectbox role-match                                   |
| `tests/test_app_geometry_viewer.py` (NEW) | test (pure unit)                | transform                        | `tests/test_app_csv_export.py`                                          | exact (pure-fn-on-synthetic-fixture)                                      |
| `tests/test_app_field_page.py` (MODIFIED) | test (AppTest)                  | request-response                 | itself (extend)                                                         | exact (self-edit)                                                         |

**Match-quality caveat (per advisor):** `results.py` is a strong analog for the _structure_ of `geometry_viewer.py` (pure module, no `st.*`, `build_*_figure(...) -> go.Figure`, unit-testable). The `griddata → go.Heatmap` **body** has **no codebase analog** — no Plotly heatmap or scipy.griddata code exists anywhere in the repo. Do not treat the heatmap internals as copy-from-existing; source them from RESEARCH §3 (verified sketch) + `plotting2d.py` (label/colormap spec). See **No Analog Found**.

## Pattern Assignments

### `app/components/geometry_viewer.py` (component, transform) — NEW

**Analog:** `app/components/results.py` (structure/purity/signature only)

**Module header + imports pattern** — copy the purity discipline and `from __future__` (results.py lines 1-22):

```python
"""Geometry viewer: pure MeshData → Plotly figure builder.

PURE — no `st.*` calls. Turns a petringa MeshData (irregular devsim node
coords + node_values) into a go.Figure: 2D griddata heatmap or 1D bar.
Unit-testable without Streamlit or devsim, exactly like results.to_csv_bytes.
"""
from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from scipy.interpolate import griddata   # NEW dep-use; no repo precedent (see No Analog Found)

from petringa import MeshData            # import the dataclass for the type hint
```

Note: `results.py` imports `petringa` and `from petringa import SimResult`. Mirror that with `from petringa import MeshData`. Confirm `MeshData` is re-exported from top-level `petringa` (results.py imports `SimResult` this way; MeshData lives in `petringa/api/results.py` alongside it).

**Pure builder signature pattern** — mirror `build_field_figures(result: SimResult) -> tuple[go.Figure, go.Figure]` (results.py:63) and `build_cv_figure(result: SimResult) -> go.Figure` (results.py:25):

```python
QUANTITIES = {                          # friendly label -> node_values key
    "Electric field": "ElectricField",
    "Net doping": "NetDoping",
    "Electrostatic potential": "Potential",
}

def build_geometry_figure(mesh: MeshData, quantity: str) -> go.Figure:
    key = QUANTITIES[quantity]
    z = mesh.node_values[key]
    if mesh.y_coords is None:
        return _build_bar(mesh, z, key)      # 1D
    return _build_heatmap(mesh, z, key)      # 2D
```

**go.Figure construction + update_layout pattern** — copy verbatim from results.py:70-75 (single-trace `go.Figure(data=go.Scatter(...))` then `.update_layout(title=, xaxis_title=, yaxis_title=)`). Replace `go.Scatter` with `go.Heatmap` / `go.Bar`.

**cm→µm conversion (CRITICAL, no analog to copy — RESEARCH Pitfall 3):** `results.build_field_figures` explicitly does NOT multiply (`result.x` is already µm — see its docstring lines 64-69). The viewer is the OPPOSITE: `mesh.x_coords`/`y_coords` are **cm** and MUST be `* 1e4`. This is the one place the analog's behavior is inverted; call it out in the plan.

**Heatmap body** — see **No Analog Found** (RESEARCH §3 verified sketch).

---

### `app/workflows/field_map.py` (workflow, request-response) — MODIFIED

**Analog:** itself (edit in place). The full current file was read; changes are surgical.

**Change 1 — REMOVE the 2D stop-guard** (current lines 32-38, the highest-priority edit per RESEARCH Pitfall 1):

```python
# DELETE this entire block:
    if cfg.half_width_um is not None:
        st.warning(
            "These workflows are 1D-only. 2D field visualization arrives in "
            "Phase 40 (geometry viewer). Set Dimensionality to 1D in the "
            "sidebar."
        )
        st.stop()
```

The module docstring (lines 7-13) also references "Phase 40" and "2D guard here is CRITICAL" — update it to describe the new 2D-through routing.

**Change 2 — KEEP the try/except RuntimeError** exactly as-is (current lines 40-48). 2D `ramp_bias` may not converge; this is the convergence-failure guard. Do NOT convert to a pre-check (the docstring warns the guard-vs-try/except distinction). Precedent for the error-render shape is already in this file.

**Change 3 — branch the render on `result.mesh.y_coords is None`** (current lines 50-76 render unconditionally):

- 1D branch: keep existing `build_field_figures(result)` line charts + `to_csv_bytes` download + net-doping expander (lines 52-76 all stay valid — they read `result.x`/`result.y`/`metadata`, populated for 1D). Then ADD the geometry viewer (bar) below (A1: supplement).
- 2D branch: SKIP `build_field_figures` and `to_csv_bytes` (empty arrays for 2D; `to_csv_bytes` raises ValueError on 2D — it has no 2D branch, see results.py:130-141). Render ONLY the geometry viewer. Guard with `if result.mesh is not None:`.

**Change 4 — add the quantity selectbox** feeding the builder. **Analog:** `app/components/device_sidebar.py:69-71`:

```python
profile = st.sidebar.selectbox(
    "Doping profile", ["graded", "uniform"], index=0, key="cfg_doping_profile"
)
```

Mirror this exact widget shape (options list, `index=0`, persistent `key=`) in `field_map.py` (main body, not sidebar) for the quantity dropdown:

```python
from app.components.geometry_viewer import build_geometry_figure, QUANTITIES
quantity = st.selectbox("Quantity", list(QUANTITIES.keys()), index=0, key="geo_quantity")
st.plotly_chart(build_geometry_figure(result.mesh, quantity))
```

The persistent `key=` (as `cfg_doping_profile` does) is what makes the selection survive reruns — no manual session_state bookkeeping. VIZ-03 "no re-run" is satisfied because the Run button returns `False` on a selectbox rerun and `result` is read from the already-cached `session_state["field_result"]` (established Phase 39). Do NOT add `@st.cache_data` (DeviceConfig unhashable — Phase 39).

---

### `tests/test_app_geometry_viewer.py` (test, pure unit) — NEW

**Analog:** `tests/test_app_csv_export.py` (exact — pure fn on hand-built fixture, no Streamlit, no devsim)

**Test structure pattern** — copy the shape of `test_app_csv_export.py`: import the pure fn directly, build a synthetic dataclass fixture, assert on the return. Replace `SimResult`+`to_csv_bytes` with `MeshData`+`build_geometry_figure`:

```python
from __future__ import annotations
import numpy as np
import pytest
from petringa import MeshData
from app.components.geometry_viewer import build_geometry_figure, QUANTITIES

def _mesh_1d():
    return MeshData(
        x_coords=np.array([0.0, 1e-4, 2e-4]),   # cm (depth for 1D)
        y_coords=None,
        node_values={
            "ElectricField": np.array([1e5, 8e4, 5e4]),
            "Potential": np.array([0.0, -20.0, -45.0]),
            "NetDoping": np.array([3e15, 2e15, -9e13]),   # signed
        },
        regions=[], contacts=[],
    )

def _mesh_2d():
    # irregular (non-grid) node scatter to exercise griddata
    return MeshData(
        x_coords=np.array([0.0, 5e-4, 1e-3, 2e-4, 8e-4]),
        y_coords=np.array([0.0, 1e-4, 2e-4, 5e-4, 3e-4]),
        node_values={...same three keys...},
        regions=[], contacts=[],
    )
```

**Assertion targets (from RESEARCH Validation — chart content is testable here, NOT in AppTest):**

- 1D → `fig.data[0].type == "bar"`; 2D → `fig.data[0].type == "heatmap"`.
- 2D heatmap `fig.data[0].z` has 2D shape; axis titles "Lateral position (µm)" / "Depth (µm)".
- cm→µm: assert an axis coordinate is ~1e4× the raw cm value (guards Pitfall 3).
- doping log10: `build_geometry_figure(mesh, "Net doping")` does not raise on signed/zero z and z-values are log-scaled (guards Pitfall 4).
- yaxis reversed for 2D (`fig.layout.yaxis.autorange == "reversed"`).

Use `pytest.approx` for float compares (test_app_csv_export.py:56-59 precedent) and `pytest.raises` for error paths (test_app_csv_export.py:155-156).

---

### `tests/test_app_field_page.py` (test, AppTest) — MODIFIED

**Analog:** itself (extend). Full file read.

**Reuse the mock seam** — `_fake_run_field` + `monkeypatch.setattr(petringa, "run_field", fake)` (lines 27-49). For a 2D fake, return a SimResult with `x=np.array([])`, `y=np.array([])`, and a **populated `mesh`** with `y_coords` non-None:

```python
def _fake_run_field_2d(cfg, **kwargs):
    return SimResult(
        config=cfg, sim_type="field",
        x=np.array([]), y=np.array([]),   # 2D returns empty x/y
        metadata={},
        mesh=MeshData(x_coords=..., y_coords=..., node_values={...}, regions=[], contacts=[]),
    )
```

**DELETE the obsolete guard test** — `test_2d_config_warns_and_skips` (lines 66-76) asserts `len(at.warning) >= 1` and `"field_result" not in at.session_state` for a 2D config. Phase 40 removes exactly that behavior, so this test now asserts the wrong thing. Replace it with a 2D-route test: 2D config → Run → `at.exception == []`, `"field_result"` IS cached, no warning about "1D-only".

**Add selectbox tests (VIZ-03)** — AppTest 1.55 `at.selectbox` EXISTS (RESEARCH verified). Pattern:

```python
at.selectbox[0].select("Net doping").run()
assert at.exception == []
# devsim did NOT re-run: session_state["field_result"] is the same cached object,
# and the fake run_field was called only on the button click, not on the selectbox rerun.
```

Keep the existing passing tests (`test_run_caches_field_result` line 48, `test_empty_state_guard` line 78, `test_solver_convergence_failure_shows_error_not_crash` line 90) — they remain valid for the 1D path. Assertions stay limited to `at.exception`, `at.session_state`, `at.button`, `at.selectbox`, `at.warning`, `at.info`, `at.error` — there is still NO `plotly_chart`/`download_button` accessor.

---

## Shared Patterns

### Pure-figure-builder discipline (no `st.*`)

**Source:** `app/components/results.py` lines 1-22 (module docstring states the rule) + every `build_*` fn
**Apply to:** `app/components/geometry_viewer.py`

```python
from __future__ import annotations
import numpy as np
import plotly.graph_objects as go
# ... build_x(...) -> go.Figure ; go.Figure(data=go.<Trace>(...)) ; .update_layout(title=, xaxis_title=, yaxis_title=)
```

The `st.selectbox` / `st.plotly_chart` calls live ONLY in `field_map.py`; the builder stays pure so it is unit-tested like `to_csv_bytes`.

### run_field mock seam (module attribute)

**Source:** `tests/test_app_field_page.py` lines 27-49
**Apply to:** all AppTest tests touching the field page

```python
monkeypatch.setattr(petringa, "run_field", _fake)   # patch the MODULE attribute
# _fake returns a hand-built SimResult; the real devsim solve never runs
```

Never call real `petringa.run_field`/devsim in a viewer or page test (RESEARCH Pitfall 5; 2D convergence unverified).

### Streamlit widget with persistent key

**Source:** `app/components/device_sidebar.py` lines 66-71 (`st.sidebar.radio` / `st.sidebar.selectbox` with options list, `index=0`, `key=`)
**Apply to:** the VIZ-03 quantity dropdown in `field_map.py`

```python
st.selectbox("Quantity", list(QUANTITIES.keys()), index=0, key="geo_quantity")
```

The `key=` makes the selection persist across reruns automatically — this is the mechanism (with the cached `field_result`) behind VIZ-03 "without re-running".

### Label / colormap / coordinate spec (matplotlib → Plotly translation)

**Source:** `petringa/core/plotting2d.py` (authoritative spec — do NOT import it; it is matplotlib `tricontourf`)
**Apply to:** `geometry_viewer.py` axis titles, colorbar labels, colorscales

| Quantity                                                                                                                                                                                                         | Title                         | Colorbar label                 | Colorscale | Scaling                                | plotting2d.py line |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------- | ------------------------------ | ---------- | -------------------------------------- | ------------------ |
| Potential                                                                                                                                                                                                        | "2D Potential Map"            | "Potential (V)"                | RdBu_r     | linear                                 | 138, 142           |
| ElectricField                                                                                                                                                                                                    | "2D Electric Field Magnitude" | "\|Electric Field\| (V/cm)"    | viridis    | linear                                 | 227, 231           |
| NetDoping                                                                                                                                                                                                        | "2D Doping Profile"           | "log10(\|NetDoping\|) (cm^-3)" | plasma     | `log10(abs(z))`, floor `abs<1.0 → 1.0` | 259-270            |
| y-axis inverted (surface at top): matplotlib `ax.invert_yaxis()` (line 143/232) → Plotly `fig.update_yaxes(autorange="reversed")`. x-axis "Lateral position (µm)", y-axis "Depth (µm)" (lines 140-141, 229-230). |

## No Analog Found

Files/patterns with no close match in the codebase (planner uses RESEARCH.md §3 verified sketch + plotting2d.py spec instead):

| Pattern                                                  | Role      | Data Flow | Reason                                                                                                                                                                                                                              |
| -------------------------------------------------------- | --------- | --------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `griddata → go.Heatmap` body inside `geometry_viewer.py` | component | transform | No Plotly heatmap and no `scipy.interpolate.griddata` usage exists anywhere in the repo. `results.py` supplies structure/purity/signature only, NOT the heatmap internals. Source: RESEARCH §3 (empirically verified this session). |
| `go.Bar` 1D depth-profile                                | component | transform | No `go.Bar` usage in repo (results.py uses only `go.Scatter`). Trivial; RESEARCH §3 sketch.                                                                                                                                         |
| `scipy` `griddata` interpolation                         | —         | transform | First scipy-in-app-layer use; API layer uses `LinearNDInterpolator` (plotting2d.py:200, simulation.py) but that is not in `app/`. griddata is the STATE.md-locked choice.                                                           |

**Heatmap body to copy (RESEARCH §3, verified scipy 1.18.0 / plotly 6.8.0):**

```python
x_um = mesh.x_coords * 1e4                     # lateral (µm)
y_um = mesh.y_coords * 1e4                     # depth (µm)
xi = np.linspace(x_um.min(), x_um.max(), 200)  # A5: n_x≈200
yi = np.linspace(y_um.min(), y_um.max(), 100)  # n_y≈100 (epi thin)
Xi, Yi = np.meshgrid(xi, yi)
Zi = griddata((x_um, y_um), z, (Xi, Yi), method="linear")   # A3: NaN outside hull → transparent (no fill_value)
fig = go.Figure(data=go.Heatmap(x=xi, y=yi, z=Zi, colorscale=..., colorbar_title=...))
fig.update_layout(xaxis_title="Lateral position (µm)", yaxis_title="Depth (µm)", title=...)
fig.update_yaxes(autorange="reversed")
```

## Metadata

**Analog search scope:** `app/components/`, `app/workflows/`, `tests/`, `petringa/core/plotting2d.py`, `petringa/api/results.py`
**Files scanned:** 6 read verbatim (results.py, field_map.py, test_app_field_page.py, test_app_csv_export.py, api/results.py, plotting2d.py:106-275) + grep for `st.selectbox`/dimensionality across `app/`
**Key finding:** selectbox has a real analog (`device_sidebar.py:69`), overturning the assumption it was analog-less. griddata/heatmap body confirmed analog-less.
**Pattern extraction date:** 2026-07-13
