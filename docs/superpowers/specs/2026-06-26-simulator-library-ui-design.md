# Design Spec: petringa — From Notebooks to Simulator Library + UI

**Date:** 2026-06-26  
**Status:** Approved for planning  
**Scope:** 3-phase transformation of notebook-based TCAD toolkit into an installable Python library with a Streamlit UI

---

## 1. Problem Statement

The project currently has ~15k lines of physics code spread across `src/` (20 modules), 20+ Jupyter notebooks, and 25 pytest modules. The structure has grown organically: notebooks import internal classes directly, there is no stable public API, no installable package, and no way for group members to run simulations without reading the Python source. The goal is to:

1. Make the simulator usable by physics group members who are not the original developers (group Petringa + students)
2. Provide a UI experience comparable to commercial TCAD tools (Sentaurus, Silvaco) — parameterized input, geometry visualization, interactive results
3. Keep the option open for external collaborators / foundries (FBK, ST) to access a hosted version later

---

## 2. Architecture Overview

```
petringa/                       ← installable Python package
├── petringa/                   ← renamed from src/ (public API lives here)
│   ├── __init__.py             ← exports stable public API
│   ├── api/                    ← thin facades over internal modules
│   │   ├── device.py           ← DeviceConfig + build_device()
│   │   ├── simulation.py       ← run_cv(), run_cce(), run_field()
│   │   ├── damage.py           ← run_radiation_damage()
│   │   ├── microdosimetry.py   ← run_microdosimetry()
│   │   └── results.py          ← SimResult dataclass
│   ├── core/                   ← existing modules (renamed, not rewritten)
│   │   ├── sic_material.py
│   │   ├── device.py
│   │   ├── device2d.py
│   │   ├── poisson.py
│   │   ├── drift_diffusion.py
│   │   ├── cv_analysis.py
│   │   ├── charge_collection.py
│   │   ├── charge_collection_2d.py
│   │   ├── generation_profiles.py
│   │   ├── single_particle.py
│   │   ├── mc_coupling.py
│   │   ├── microdosimetry.py
│   │   ├── radiation_damage.py
│   │   ├── dark_current.py
│   │   ├── flash_recombination.py
│   │   ├── temperature_sweep.py
│   │   ├── transient.py
│   │   ├── incomplete_ionization.py
│   │   ├── alternative_structures.py
│   │   ├── optimization.py
│   │   ├── validation.py
│   │   ├── analytical.py
│   │   ├── plotting.py
│   │   ├── plotting2d.py
│   │   └── devsim_reset.py
│   └── _version.py
├── app/                        ← Streamlit UI (Phase 2+)
│   ├── main.py                 ← entry point: streamlit run app/main.py
│   ├── pages/
│   │   ├── 01_device_config.py
│   │   ├── 02_cv_analysis.py
│   │   ├── 03_charge_collection.py
│   │   ├── 04_field_map.py
│   │   ├── 05_radiation_damage.py
│   │   ├── 06_dark_current.py
│   │   ├── 07_microdosimetry.py
│   │   └── 08_batch_sweep.py
│   └── components/
│       ├── geometry_viewer.py  ← 2D plotly mesh/field visualizer
│       ├── result_panel.py     ← standard result display widget
│       └── param_sidebar.py    ← reusable parameter input sidebar
├── notebooks/                  ← kept as-is; update imports to petringa.*
├── tests/                      ← kept as-is; update imports
├── pyproject.toml              ← replaces requirements.txt
└── docs/
```

---

## 3. Public API Design

The public API lives in `petringa/api/`. Internal modules in `petringa/core/` are not part of the public contract and may change. Users (including the Streamlit UI) only import from `petringa.*`.

### 3.1 DeviceConfig

```python
from petringa import DeviceConfig

cfg = DeviceConfig(
    # Geometry
    epi_thickness_um=10.0,          # float, default 10 um
    substrate_thickness_um=1.0,     # float, default 1 um
    half_width_um=None,             # float or None → 1D device if None, 2D if set

    # Doping
    N_A=1e19,                       # float cm^-3, p+ substrate
    doping_profile="graded",        # "uniform" | "graded"
    N_D=None,                       # float cm^-3 (uniform) or None
    N_D_junction=2.93e15,           # float cm^-3 (graded)
    N_D_bulk=8.82e13,               # float cm^-3 (graded)
    L_transition_um=0.987,          # float um (graded)

    # Operating conditions
    T=300.0,                        # float K
    area_cm2=1e-4,                  # float cm^2, detector area
)
```

### 3.2 SimResult

```python
@dataclass
class SimResult:
    config: DeviceConfig
    sim_type: str                   # "cv" | "cce" | "field" | "damage" | ...
    x: np.ndarray                   # primary axis (bias, depth, fluence, ...)
    y: np.ndarray                   # primary output
    metadata: dict                  # sim-type-specific extras
    mesh: MeshData | None           # populated after build, used by geometry viewer
```

### 3.3 MeshData (geometry viewer contract)

```python
@dataclass
class MeshData:
    x_coords: np.ndarray            # node x coordinates (cm)
    y_coords: np.ndarray | None     # node y coordinates (cm), None for 1D
    node_values: dict[str, np.ndarray]  # "NetDoping", "ElectricField", etc.
    regions: list[dict]             # [{name, x_min, x_max, y_min, y_max}]
    contacts: list[dict]            # [{name, position}]
```

`MeshData` is populated **after** `devsim` builds the device (post-build extraction via `devsim.get_node_model_values()`). The geometry viewer reads `MeshData` and renders via Plotly — it does not call `devsim` itself.

### 3.4 Simulation functions

```python
from petringa import run_cv, run_cce, run_field, run_radiation_damage, run_microdosimetry

# C-V sweep
result = run_cv(cfg, v_start=0, v_stop=-200, n_points=40)
# → SimResult with x=bias_V, y=capacitance_F_cm2, metadata includes W, 1/C^2

# CCE vs bias
result = run_cce(cfg, v_start=-10, v_stop=-200, n_points=30)
# → SimResult with x=bias_V, y=CCE (0–1)

# Field / potential profile at fixed bias
result = run_field(cfg, bias_V=-100)
# → SimResult with x=depth_um, y=field_V_cm, metadata includes potential, doping
# result.mesh is populated → geometry viewer can render

# Radiation damage sweep
result = run_radiation_damage(cfg, fluences=[1e12, 1e13, 1e14], proton_energy_MeV=5.6)
# → SimResult with x=fluence, y=CCE_post_damage

# Microdosimetry
result = run_microdosimetry(cfg, mc_csv_path, sv_thickness_um=10, sv_width_um=150)
# → SimResult with x=lineal_energy_keV_um, y=yd_y (dose spectrum)
```

### 3.5 Batch sweep

```python
from petringa import ParametricSweep

sweep = ParametricSweep(
    base_config=cfg,
    param="epi_thickness_um",
    values=[5, 10, 20, 50],
    sim_fn=run_cce,
    sim_kwargs={"v_start": -10, "v_stop": -200, "n_points": 20},
)
results = sweep.run()    # list[SimResult]
```

---

## 4. Phase Breakdown

### Phase 1 — Installable Library + Vertical Slice (C-V)

**Goal:** `pip install -e .` works; `from petringa import DeviceConfig, run_cv` works; all existing tests pass unchanged.

**Scope:**

- Add `pyproject.toml` (replaces `requirements.txt`); declare `petringa` as package
- Rename `src/` → `petringa/core/` (mechanical move, zero logic change)
- Create `petringa/api/` with `DeviceConfig`, `SimResult`, `MeshData` dataclasses
- Implement `petringa/api/simulation.py::run_cv()` as thin facade over `core/cv_analysis.py`
- Implement `petringa/api/simulation.py::run_field()` (populates `MeshData` via devsim node extraction)
- Update all `from src.X import Y` → `from petringa.core.X import Y` across `src/`, `tests/`, `notebooks/`
- `petringa/__init__.py` exports: `DeviceConfig`, `SimResult`, `run_cv`, `run_field`

**Acceptance gate:** `pytest -q` green (all 25 modules); `v3_frozen.json` baseline unchanged.

**Explicitly out of scope:** UI, any physics change, kappa/SRIM data issues.

**Vertical slice validation:** After Phase 1, write a one-file script `examples/cv_example.py` that does:

```python
from petringa import DeviceConfig, run_cv
cfg = DeviceConfig()
result = run_cv(cfg, v_start=0, v_stop=-200, n_points=20)
print(result.y)   # capacitance array
```

If this runs cleanly, the API design is validated before touching remaining modules.

---

### Phase 2 — Streamlit MVP + Geometry Viewer

**Goal:** `streamlit run app/main.py` launches a working UI covering the 3 core workflows (C-V, CCE, field map). Geometry viewer shows 2D mesh with field/doping overlay.

**Scope:**

- Expand library API: add `run_cce()` façade (over `core/charge_collection.py`)
- `app/main.py` — multi-page Streamlit shell
- `app/pages/01_device_config.py` — sidebar form for all `DeviceConfig` fields; stores config in `st.session_state`
- `app/pages/02_cv_analysis.py` — calls `run_cv()`, renders C-V and 1/C² vs V plots (matplotlib/plotly)
- `app/pages/03_charge_collection.py` — calls `run_cce()`, renders CCE vs bias
- `app/pages/04_field_map.py` — calls `run_field()`, renders electric field and potential depth profiles; triggers geometry viewer for 2D devices
- `app/components/geometry_viewer.py` — takes `MeshData`, renders 2D heatmap (Plotly `go.Heatmap` for field/doping on triangular mesh via interpolation to regular grid)
- `app/components/result_panel.py` — download button (CSV export of `SimResult.x`, `SimResult.y`)

**Geometry viewer detail:** For 1D devices, renders a depth-profile bar. For 2D devices (when `half_width_um` is set in `DeviceConfig`), renders a 2D colormap. Node coordinates and field values extracted from `MeshData` (populated by `run_field()`). Interpolated to a regular grid using `scipy.interpolate.griddata` for `go.Heatmap`.

**UI layout per page:**

```
┌─ Sidebar ──────────┐  ┌─ Main panel ─────────────────────┐
│ DeviceConfig form   │  │ [Run simulation] button           │
│ (all params)        │  │                                   │
│                     │  │ Plot area (plotly interactive)    │
│                     │  │                                   │
│                     │  │ [Download CSV]                    │
└────────────────────┘  └──────────────────────────────────┘
```

**Acceptance gate:** All 3 pages produce correct outputs matching notebook results for the default `DeviceConfig`. Geometry viewer renders without error for both 1D and 2D configs.

**Explicitly out of scope:** Radiation damage, microdosimetry, batch sweep, dark current, temperature sweep, FLASH.

---

### Phase 3 — Feature Complete

**Goal:** All simulation workflows from the 20 notebooks are accessible via the UI and the library API. Batch parametric sweeps launchable from the UI.

**Scope — library API additions:**

- `run_radiation_damage()` — facade over `core/radiation_damage.py` + `core/cv_analysis.py`
- `run_dark_current()` — facade over `core/dark_current.py`
- `run_temperature_sweep()` — facade over `core/temperature_sweep.py`
- `run_flash_recombination()` — facade over `core/flash_recombination.py`
- `run_transient()` — facade over `core/transient.py`
- `run_microdosimetry()` — facade over `core/microdosimetry.py` + `core/mc_coupling.py`
- `ParametricSweep` class (see §3.5)

**Scope — UI additions:**

- `app/pages/05_radiation_damage.py` — fluence sweep, CCE vs fluence, defect concentrations
- `app/pages/06_dark_current.py` — dark current vs temperature, trap contributions
- `app/pages/07_microdosimetry.py` — MC CSV upload, lineal energy spectrum (y·d(y) vs log y), y_F/y_D readout
- `app/pages/08_batch_sweep.py` — select parameter + range + sim type, run sweep, overlaid plots, bulk CSV download
- `app/pages/09_alternative_structures.py` — mesa, 3D-electrode, ΔE-E telescope configs

**Known data-blocked items (explicitly excluded from Phase 3 scope):**

- `kappa(E)` NIEL hardness factors — placeholders in `core/radiation_damage.py`; any absolute Phi_crit values are unvalidated. UI will display a warning banner on the radiation damage page.
- Graded doping calibration from SRIM data — Phase 3 uses calibrated defaults from Phase 26; new SRIM data is a separate future phase.

**Acceptance gate:** All 20 notebooks can be replaced by equivalent UI workflows. `ParametricSweep` tested with ≥2 parameters. Microdosimetry page produces y·d(y) spectrum matching notebook 18 for the default MC CSV.

---

## 5. Deployment

**Local (all phases):**

```bash
pip install -e .               # or: uv pip install -e .
streamlit run app/main.py      # Phase 2+
```

**Shared lab server (Phase 2+):**

```bash
streamlit run app/main.py --server.port 8501 --server.address 0.0.0.0
```

No Docker required. Streamlit's built-in server is sufficient for a group of <20 users.

**External / Streamlit Community Cloud (future Phase 4):**
devsim is a C extension not available on Streamlit Community Cloud's sandbox. External deployment requires a VPS with devsim pre-installed. Out of scope for Phase 3.

---

## 6. Dependency Management

Replace `requirements.txt` with `pyproject.toml`:

```toml
[project]
name = "petringa"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = [
    "devsim>=2.10.0",
    "numpy>=1.24",
    "scipy>=1.11",
    "matplotlib>=3.7",
    "plotly>=5.0",
    "streamlit>=1.35",
    "pandas>=2.0",
]

[project.optional-dependencies]
dev = ["pytest>=7.0", "jupyter>=1.0"]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"
```

Use `uv` for all dependency operations (project convention).

---

## 7. Testing Strategy

**Phase 1:**

- All existing 25 pytest modules pass without modification (import paths updated `src.X` → `petringa.core.X`)
- `v3_frozen.json` baseline regression unchanged
- One new test: `tests/test_api_cv.py` — calls `run_cv(DeviceConfig())`, checks output shape and that C values are in physically reasonable range

**Phase 2:**

- UI tested manually against notebook 02 (C-V) and notebook 03 (CCE) golden outputs
- `tests/test_api_field.py` — calls `run_field()`, checks `MeshData` is populated with correct node count
- No Streamlit unit tests (UI layer is thin, physics tested via API tests)

**Phase 3:**

- `tests/test_api_*.py` for each new facade function
- `tests/test_parametric_sweep.py` — sweep over 2 values, checks result list length
- Microdosimetry: `tests/test_api_microdosimetry.py` against `data/synthetic_mc_events.csv`

**Slow test marking:** Existing `@pytest.mark.slow` convention maintained. DD-heavy tests continue to be run one module at a time.

---

## 8. Key Constraints and Non-Goals

**Constraints:**

- Physics logic must not change during refactor (Phases 1–2). All physics changes are separate GSD phases.
- `devsim` CGS unit convention (cm, cm⁻³, V, s, F/cm) preserved throughout the public API.
- `MeshData` is extracted post-build from devsim; the geometry viewer never calls devsim directly.
- kappa(E) data-blocked issue documented with UI warning banner; no fabricated values.

**Non-goals:**

- 3D geometry viewer (devsim devices are 1D/2D; 2D Plotly heatmap is sufficient)
- Real-time live simulation (devsim runs are seconds to minutes; UI triggers on button click)
- Web deployment with devsim (devsim is a C extension; needs local/VPS install)
- Any new physics (new defect models, new materials, etc.)
- Mask layout / GDS output (foundry-facing work is in `deliverables/`, not this scope)

---

## 9. GSD Phase Mapping

| GSD Phase | Name              | Deliverable                                                     |
| --------- | ----------------- | --------------------------------------------------------------- |
| Phase A   | Library packaging | `petringa` installable, C-V API slice, all tests green          |
| Phase B   | Streamlit MVP     | Working UI: device config + C-V + CCE + field + geometry viewer |
| Phase C   | Feature complete  | All workflows in UI + batch sweep + microdosimetry              |

Each phase should have its own GSD PLAN.md. Phase A is a prerequisite for Phase B; Phase B is a prerequisite for Phase C.

**Suggested GSD phase goal statements:**

- **Phase A:** "Users can `pip install -e .` the petringa package and call `run_cv(DeviceConfig())` from a Python script. All existing tests pass unchanged."
- **Phase B:** "Users can `streamlit run app/main.py`, configure a device via a form, and interactively view C-V curves, CCE vs bias, electric field profiles, and a 2D geometry heatmap."
- **Phase C:** "All 20 notebook workflows are reproducible via the Streamlit UI. Parametric sweeps are launchable from the batch page. The microdosimetry page accepts a Geant4 CSV and outputs y·d(y) spectra."
