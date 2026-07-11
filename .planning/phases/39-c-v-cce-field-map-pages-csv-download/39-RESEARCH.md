# Phase 39: C-V, CCE, Field Map Pages + CSV Download — Research

**Researched:** 2026-07-11
**Domain:** Streamlit UI (v1.55.0) wiring petringa API facades to interactive Plotly charts + CSV export
**Confidence:** HIGH (all contracts read verbatim from source; all technical claims empirically verified in this session)

## Summary

Phase 39 wires three existing, fully-implemented petringa API facades (`run_cv`, `run_cce`, `run_field`) into three placeholder Streamlit pages that currently only echo the device config as JSON. The API contracts are precise and stable — every `run_*` returns a uniform `SimResult` dataclass (`x`, `y`, `metadata`, `mesh`), and I quote the exact fields below. The heavy lifting (devsim solves, unit conversions, physics) is done; Phase 39 is pure presentation + serialization glue.

Three findings materially shape the plan and are not in the success criteria: **(1)** the sidebar can produce a 2D config, but all three workflows are 1D-only — `run_cv`/`run_cce` raise `NotImplementedError` for 2D and `run_field` silently returns empty `x`/`y` — so every page needs a 1D-only pre-check guard. **(2)** `plotly` is declared in `pyproject.toml` (`plotly>=5.0`) but is **NOT installed** in the venv (verified: `ModuleNotFoundError`), and no code anywhere imports it yet — it must be installed via `uv` (per project memory) as Wave 0, or every AppTest breaks at import. **(3)** `DeviceConfig` is unhashable (verified `TypeError`), so `st.cache_data` keyed on config is impossible — result caching MUST be manual `st.session_state`, and this is load-bearing for UI-06 (the download button rerun would otherwise wipe the result).

**Primary recommendation:** Build a small shared `app/components/results.py` helper module (Plotly figure builders + a `to_csv_bytes()` serializer with commented metadata header) plus a uniform per-page pattern: 1D-guard → "Run simulation" button → cache `SimResult` in `st.session_state[f"{sim}_result"]` → render Plotly charts from the cached result → `st.download_button` reading pre-built CSV bytes from the cached result.

## Architectural Responsibility Map

| Capability             | Primary Tier                | Secondary Tier                                | Rationale                                                                          |
| ---------------------- | --------------------------- | --------------------------------------------- | ---------------------------------------------------------------------------------- |
| Physics solve (devsim) | petringa API (`run_*`)      | —                                             | Already implemented; pages call, never touch devsim                                |
| Result caching         | Streamlit session_state     | —                                             | Reruns are per-interaction; unhashable config forbids `st.cache_data`              |
| Chart construction     | app UI (Plotly builders)    | —                                             | No reusable Plotly code exists; matplotlib `plotting.py` not reusable in Streamlit |
| CSV serialization      | app UI helper               | pandas                                        | No API-level CSV export exists; pattern is `df.to_csv`                             |
| Provenance metadata    | app UI helper               | `dataclasses.asdict` + `petringa.__version__` | DeviceConfig carried on every SimResult; version available                         |
| 2D dispatch guard      | app UI (per-page pre-check) | —                                             | Workflows are 1D-only; sidebar can emit 2D config                                  |

## User Constraints

**No CONTEXT.md exists** — discuss-phase was deliberately skipped ("facciamo la plan direttamente"). The planner must make default decisions. This research provides concrete `[ASSUMED]` recommendations (CSV format, metadata, extra plots) that the planner should adopt as defaults; all such items are listed in the Assumptions Log for later user confirmation.

## Phase Requirements

| ID    | Description                                                                 | Research Support                                                                                       |
| ----- | --------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| UI-03 | C-V page: run `run_cv()`, show Plotly C-V curve + 1/C² vs V (Mott-Schottky) | `run_cv` contract (simulation.py:34-131); `metadata["one_over_C_squared"]` (line 108) already computed |
| UI-04 | CCE page: run `run_cce()`, show Plotly CCE vs bias                          | `run_cce` contract (simulation.py:320-411); `y` in [0,1] (line 372)                                    |
| UI-05 | Field map page: run `run_field()`, show Plotly E-field + potential vs depth | `run_field` contract (simulation.py:133-317); both E-field and Potential returned in one call (see §6) |
| UI-06 | Download results as CSV from any result page                                | No existing CSV export in app; `st.download_button` + session_state caching pattern (see §5)           |

---

## 1. API Contract (verbatim from source)

All three facades are re-exported at package top level (`petringa/__init__.py:8-18`), so pages import `from petringa import run_cv, run_cce, run_field, DeviceConfig, SimResult`. `petringa/api/__init__.py` is empty (just a docstring) — the real exports live in `petringa/__init__.py`.

### `SimResult` / `MeshData` envelope (`petringa/api/results.py:20-46`)

```python
@dataclass
class MeshData:
    x_coords: np.ndarray  # node x coordinates (cm)
    y_coords: "np.ndarray | None"  # node y coordinates (cm), None for 1D
    node_values: dict[str, np.ndarray]  # "NetDoping", "ElectricField", etc.
    regions: list[dict]  # [{name, x_min, x_max, y_min, y_max}]
    contacts: list[dict]  # [{name, position}]

@dataclass
class SimResult:
    config: "DeviceConfig"
    sim_type: str  # "cv" | "cce" | "field" | "damage" | ...
    x: np.ndarray  # primary axis (bias, depth, fluence, ...)
    y: np.ndarray  # primary output
    metadata: dict = field(default_factory=dict)  # sim-type-specific extras
    mesh: "MeshData | None" = None
```

### `run_cv` (`petringa/api/simulation.py:34-39`)

```python
def run_cv(
    config: DeviceConfig,
    v_start: float = 0.0,
    v_stop: float = -200.0,
    n_points: int = 40,
) -> SimResult:
```

Returns (docstring lines 64-69, code lines 106-119): `sim_type="cv"`, `x=cv_result["voltages"]` (bias V), `y=cv_result["capacitance"]` (F), `metadata` keys: `"depletion_widths"` (cm), `"one_over_C_squared"` (= `1.0 / capacitance**2`), `"area_cm2"`. `mesh=None`.

**2D guard (lines 82-88):** raises `NotImplementedError` if `config.half_width_um is not None`.

### `run_cce` (`petringa/api/simulation.py:320-325`)

```python
def run_cce(
    config: DeviceConfig,
    v_start: float = -10.0,
    v_stop: float = -200.0,
    n_points: int = 30,
) -> SimResult:
```

Returns (docstring lines 371-375, code lines 401-411): `sim_type="cce"`, `x=result["voltages"]` (bias V), `y=result["cce_values"]` (dimensionless, in [0,1]), `metadata` keys: `"I_collected"` (A/cm²), `"I_generated"` (A/cm²). `mesh=None`.

**2D guard (lines 377-383):** raises `NotImplementedError` if `config.half_width_um is not None`.

### `run_field` (`petringa/api/simulation.py:133`)

```python
def run_field(config: DeviceConfig, bias_V: float = -100.0) -> SimResult:
```

Returns (docstring lines 160-172, code lines 299-306). For **1D** (`config.half_width_um is None`):

- `sim_type="field"`
- `x = x_coords * 1e4` → **depth in µm** (line 296) — NOTE: `x` is µm, but `mesh.x_coords` is raw **cm**.
- `y = field_nodes` → node-aligned ElectricField magnitude (V/cm) (line 297)
- `metadata`: `"bias_V"` (line 279), `"potential"` (V, node array, line 280), `"net_doping"` (cm⁻³, node array, line 281)
- `mesh`: populated `MeshData` with `node_values` `{"NetDoping", "Potential", "ElectricField"}` (lines 266-276)

For **2D** (`config.half_width_um is not None`): `run_field` does **NOT raise** — it returns `x=np.array([])`, `y=np.array([])` (lines 293-294) and routes all data to `mesh`. This is the CR-01-bug-avoidance behavior (lines 285-292).

**DeviceConfig** (`petringa/api/device.py:25-43`) is a plain `@dataclass` with 11 fields, all defaulted. Key field for the 2D guard: `half_width_um: Optional[float] = None` (line 31) — `None`=1D, float=2D.

---

## 2. Existing Streamlit Patterns to Follow

**`render()` contract:** Every page module exposes a top-level `def render() -> None` with all `st.*` calls inside it (no module-level side effects), enabling both `st.Page(render, ...)` registration in `app/main.py:19-25` and headless `AppTest.from_function(render)` (home.py:1-6 documents this explicitly).

**Empty-state guard pattern** (identical in cv.py/cce.py/field_map.py/home.py):

```python
cfg = st.session_state.get("device_config")
if cfg is None:
    st.info("Configure a device in the sidebar to begin.")
    st.stop()
```

The exact string `"Configure a device in the sidebar to begin."` is asserted by `test_empty_state_guard` (test_app_pages.py:90) — preserve it verbatim.

**Config JSON echo idiom:** `st.json({k: getattr(cfg, k) for k in cfg.__dataclass_fields__})` — the established way to display config; reusable for a "device parameters used" expander on result pages.

**session_state key convention:** The only non-widget key is `st.session_state["device_config"]` (device_sidebar.py:150), written by `render_device_sidebar()` which runs on **every** page before `pg.run()` (main.py:36). Sidebar widget keys are all `cfg_*` prefixed (e.g. `cfg_epi_thickness_um`). Recommend new result keys namespaced distinctly: `cv_result`, `cce_result`, `field_result` (avoid `cfg_*`).

**Docstring style:** Module docstring explaining the page's role + phase provenance; the current placeholder caption `"Running this simulation is implemented in Phase 39."` should be removed.

**Test patterns to extend (`tests/test_app_pages.py`):**

- `test_first_edit_survives_rerun` (line 18) — boots real `app/main.py`, sets a sidebar widget, reruns, asserts persistence. Accessor: `at.sidebar.number_input[0]` is "Epi thickness (µm)".
- `test_no_magic_pages_directory` (line 53) — structural guard: `app/pages` must NOT exist, `app/workflows` must. **Do not create an `app/pages` dir.**
- `test_empty_state_guard` (line 78) — `AppTest.from_function` on a wrapper that imports `render` **inside its body** (required: `from_function` bodies must be self-contained), asserts no exception + info text.
- `test_nav_sidebar_smoke` (line 93) — boots `app/main.py`, asserts `at.exception == []` and `device_config` in session_state.

---

## 3. Plotly Usage Precedent

**There is NO existing Plotly code anywhere** in `petringa/` or `notebooks/` (verified: `grep -rln "plotly|go.Figure|px\." petringa/ notebooks/` returns nothing). The UI must build all charts from scratch.

`petringa/core/plotting.py` and `plotting2d.py` exist but are **matplotlib-only** (`import matplotlib.pyplot as plt`, return `matplotlib.axes.Axes`) — **not reusable** in a Streamlit `st.plotly_chart`. However, they are an authoritative reference for **data keys, axis labels, and unit conversions** to mirror in the Plotly builders:

- `plot_cv_curve` (plotting.py:339) consumes `{"voltages", "capacitance", "one_over_C_squared", "depletion_widths"}` — same keys `run_cv` produces. Labels: C-V → "Capacitance (F/cm²)"; Mott-Schottky → "1/C² (cm⁴/F²)", title "Mott-Schottky Plot (1/C² vs V)", x "Voltage (V)".
- `plot_cce_vs_bias` (plotting.py:472) consumes `{"voltages", "cce_values"}`, plots `|V|` on x (`np.abs`), y-limit `[0, 1.1]`, ref line at CCE=1.0, labels "|Reverse Bias| (V)" / "Charge Collection Efficiency".
- `plot_electric_field` (plotting.py:36) converts x cm→µm (`*1e4`) — but `run_field.x` is **already** in µm, so the Plotly field builder must NOT re-multiply.

**Recommendation:** Write Plotly `go.Figure` builders in a new `app/components/results.py`, using these labels/units as the spec. Do not import `plotting.py`.

---

## 4. CSV Export Precedent

**No CSV export exists in `app/`** (confirmed). The project-wide convention is plain `df.to_csv(path, index=False)`:

- `tests/test_mc_coupling.py:71,168,194`: `df.to_csv(path, index=False)`
- `scripts/create_notebook_17.py:285`: `synthetic_df.to_csv(csv_path, index=False)`

**No metadata-header (commented `#` lines) CSV precedent exists** anywhere (verified: grep for `comment=`, `StringIO`, header-comment patterns returns nothing). So the traceability-header format is a new decision → tag `[ASSUMED]`.

**No `to_dict()` on DeviceConfig** — but `dataclasses.asdict(config)` works cleanly (verified: returns all 11 fields as a plain dict). The app already uses `{k: getattr(cfg, k) for k in cfg.__dataclass_fields__}`. `petringa.__version__ = "5.0.0"` (`petringa/_version.py`) is available for the metadata header.

`petringa/api/sweep.py` (`ParametricSweep`) uses `dataclasses.replace` for config cloning but has **no** `to_csv`/`to_dataframe` method — no serialization convention to inherit there.

### Recommended CSV format (all `[ASSUMED]`)

Build CSV bytes in-memory (`df.to_csv(index=False).encode()` or a `StringIO`) — no temp file. Column recommendations per page:

- **C-V** (`cv_YYYYMMDD-HHMMSS.csv`): `bias_V, capacitance_F, one_over_C2_cm4_per_F2, depletion_width_cm`
  (from `result.x`, `result.y`, `metadata["one_over_C_squared"]`, `metadata["depletion_widths"]`)
- **CCE** (`cce_...csv`): `bias_V, CCE, I_collected_A_per_cm2, I_generated_A_per_cm2`
  (`I_generated` may be scalar — broadcast or place in metadata header)
- **Field 1D** (`field_...csv`): `depth_um, ElectricField_V_per_cm, Potential_V, NetDoping_cm-3`
  (`result.x` is µm; `metadata["potential"]`, `metadata["net_doping"]` are node arrays; all node-aligned & same length)

### Recommended traceability metadata header (`[ASSUMED]`)

Prepend commented lines before the column table — this matches the PSTAR/SRIM commented-header convention the project already _consumes_ elsewhere, and is the audit-flagged gap (a bare number CSV isn't traceable to its run):

```
# petringa SiC TCAD Simulator — <sim_type> result
# software_version: 5.0.0
# generated: 2026-07-11T14:32:05Z
# bias_range_V: <v_start> to <v_stop>, n_points: <n>
# device: epi_thickness_um=10.0, doping_profile=graded, N_D_junction=2.93e15, ... (all 11 DeviceConfig fields via asdict)
# units: <per-column units>
```

**Tradeoff to document:** re-reading needs `pandas.read_csv(..., comment='#')`; Excel shows `#` lines as literal rows. Alternative (simpler, less rich): a second "parameters" CSV or an in-UI expander showing `asdict(config)`. Recommend the commented-header approach for single-file traceability.

---

## 5. `st.download_button` Mechanics + Caching (load-bearing)

**Verified facts:**

- `DeviceConfig` is **unhashable** (`TypeError: unhashable type: 'DeviceConfig'`). Therefore `st.cache_data` keyed on a config argument is **impossible** — do not attempt `@st.cache_data` on `run_cv`/etc.
- devsim solves are expensive (each `run_*` calls `reset_devsim_fully()` + builds/ramps a device). Recomputing on every rerun is unacceptable.

**Required pattern (not an optimization):** Streamlit reruns the whole script on _every_ widget interaction. On a Download-button click, the script reruns → the "Run simulation" button returns `False` → without a cached result, the plot and download vanish. So:

```python
if st.button("Run simulation"):
    st.session_state["cv_result"] = run_cv(cfg)   # or petringa.run_cv (see §testing)

result = st.session_state.get("cv_result")
if result is not None:
    st.plotly_chart(build_cv_figure(result))
    st.plotly_chart(build_mott_schottky_figure(result))
    st.download_button(
        "Download CSV",
        data=to_csv_bytes(result),      # bytes prepared before rerun — no server temp file
        file_name="cv_result.csv",
        mime="text/csv",
    )
```

`st.download_button` requires the bytes to exist at render time (no server-side temp file). Because `to_csv_bytes` runs on the cached `SimResult`, clicking Download does NOT recompute the simulation.

---

## 6. Field Map — Both E-field AND Potential from ONE call

**Confirmed: `run_field` returns both in a single call, no second API call needed.** For 1D (simulation.py:266-306):

- ElectricField: `result.y` (node array, V/cm) AND `mesh.node_values["ElectricField"]`
- Potential: `metadata["potential"]` (node array, V) AND `mesh.node_values["Potential"]`
- Depth axis: `result.x` (µm, node-aligned to both arrays)

So success criterion 3 (E-field AND potential vs depth) is satisfied by one `run_field(cfg)` call. Plot both as two Plotly subplots/figures sharing the depth x-axis. `metadata["net_doping"]` is a bonus third quantity available for a doping-profile plot if desired.

---

## Extra Plots Worth Surfacing (per-page, `[ASSUMED]`)

- **C-V:** required = C-V curve + Mott-Schottky. Bonus available at zero extra cost: depletion width vs V (`metadata["depletion_widths"]`).
- **CCE:** required = CCE vs bias. `I_collected`/`I_generated` are in metadata — surface as a small caption/metric, not a required plot.
- **Field map:** required = E-field + potential vs depth. Bonus: net doping vs depth (`metadata["net_doping"]`) — cheap and physically informative. Full 2D field overlay is explicitly **Phase 40** (geometry viewer), not this phase.

Keep the required set minimal; offer bonuses as opt-in (e.g. an expander) so pages stay clean.

---

## Common Pitfalls

### Pitfall 1: 2D config crashes the workflow (HIGHEST PRIORITY)

**What goes wrong:** Sidebar `Dimensionality="2D"` sets `half_width_um=50.0` (device_sidebar.py:86-89). Then: C-V and CCE raise uncaught `NotImplementedError`; field map silently plots nothing (`x=y=[]`).
**How to avoid:** Add a pre-check at the top of every Run handler:

```python
if cfg.half_width_um is not None:
    st.warning("These workflows are 1D-only. 2D field visualization arrives in Phase 40 (geometry viewer). Set Dimensionality to 1D in the sidebar.")
    st.stop()
```

Must be a **pre-check, not try/except** — `run_field` doesn't raise, so an exception handler wouldn't catch its empty-array case.

### Pitfall 2: plotly not installed breaks ALL tests

**What goes wrong:** Pages will `import plotly.graph_objects`. Plotly is in `pyproject.toml` but absent from the venv. After Phase 39 edits, `AppTest.from_file("app/main.py")` fails at import, breaking the entire existing suite.
**How to avoid:** Wave 0 — install via **uv** (project uses uv, not pip): `uv sync` or `uv add plotly`. Verify `python -c "import plotly"` succeeds before writing page code.

### Pitfall 3: unit mismatch — `result.x` (µm) vs `mesh.x_coords` (cm)

**What goes wrong:** `run_field.x = x_coords * 1e4` (µm), but `mesh.x_coords` is raw cm. Mixing them mislabels the depth axis by 10⁴.
**How to avoid:** For 1D depth plots, use `result.x` (already µm) and `result.y` / `metadata["potential"]`. Reserve `mesh` for Phase 40.

### Pitfall 4: download rerun wipes results

**What goes wrong:** Without session_state caching, clicking Download reruns the script, the Run button returns False, and the plot+download disappear (and a naive re-run would re-solve devsim). See §5.
**How to avoid:** Cache the `SimResult` in `st.session_state` keyed per sim type.

---

## Don't Hand-Roll

| Problem              | Don't Build                          | Use Instead                                 | Why                                          |
| -------------------- | ------------------------------------ | ------------------------------------------- | -------------------------------------------- |
| CSV string           | Manual string joins                  | `pandas.DataFrame(...).to_csv(index=False)` | Project convention; handles quoting/escaping |
| Config serialization | Manual field listing                 | `dataclasses.asdict(config)`                | Verified working; auto-tracks all 11 fields  |
| Result caching       | Custom memoization / `st.cache_data` | Manual `st.session_state`                   | Config unhashable → `st.cache_data` raises   |
| Physics/units        | Any recompute in the page            | The `SimResult` from `run_*`                | All physics + unit conversion already done   |

---

## Validation Architecture

`workflow.nyquist_validation` is **absent** in config.json → treat as **enabled**.

### Test Framework

| Property    | Value                                                               |
| ----------- | ------------------------------------------------------------------- |
| Framework   | pytest ≥7.0 (dev extra), Streamlit `AppTest` (streamlit.testing.v1) |
| Config file | `pyproject.toml` `[tool.pytest.ini_options]`, `pythonpath = ["."]`  |
| Quick run   | `pytest tests/test_app_pages.py -x`                                 |
| Full suite  | `pytest`                                                            |

### Phase Requirements → Test Map

| Req      | Behavior                           | Test Type                          | Command                                  | Exists?   |
| -------- | ---------------------------------- | ---------------------------------- | ---------------------------------------- | --------- |
| UI-03    | C-V run → charts, no exception     | AppTest + monkeypatch              | `pytest tests/test_app_cv_page.py -x`    | ❌ Wave 0 |
| UI-04    | CCE run → chart, no exception      | AppTest + monkeypatch              | `pytest tests/test_app_cce_page.py -x`   | ❌ Wave 0 |
| UI-05    | Field run → E+potential charts     | AppTest + monkeypatch              | `pytest tests/test_app_field_page.py -x` | ❌ Wave 0 |
| UI-06    | CSV bytes correct (columns/header) | Pure unit test on `to_csv_bytes()` | `pytest tests/test_app_csv_export.py -x` | ❌ Wave 0 |
| 2D guard | 2D config → warning + no crash     | AppTest                            | (in page test files)                     | ❌ Wave 0 |

**Critical testing decision (verified):** AppTest 1.55 has **no `plotly_chart` and no `download_button` accessor** (verified: accessor list contains `button`, `caption`, `title`, `info`, `json`, `dataframe`, `exception`, etc. — no chart/download). `at.get("plotly_chart")` also fails (element_type must be an AppTest attribute). So AppTest tests assert on: `at.exception == []`, `at.button` label, `at.caption`/`at.markdown` marker text, and `session_state` result keys — NOT on the chart or download widget directly. Put the CSV-content assertions in a **pure unit test** on the serializer (`to_csv_bytes(SimResult) -> bytes`), which is fully testable without Streamlit or devsim.

**Mocking (avoid real devsim solves):** Reference the facades as `petringa.run_cv(...)` (i.e. `import petringa; petringa.run_cv`) so tests can `monkeypatch.setattr(petringa, "run_cv", fake_run_cv)`. `from petringa import run_cv` binds at page-import time and would require patching the page-module attribute instead — more fragile. Wave 0 must include a spike verifying AppTest + monkeypatch actually intercepts the call (this dictates page import structure).

### Wave 0 Gaps

- [ ] Install `plotly` via uv (blocking — all AppTests fail without it)
- [ ] Verify AppTest + monkeypatch interception of `petringa.run_*` (dictates page structure)
- [ ] `app/components/results.py` — shared Plotly builders + `to_csv_bytes()`
- [ ] `tests/test_app_csv_export.py` — pure serializer test (columns, metadata header, no-devsim)
- [ ] Per-page AppTest files (cv/cce/field) with mocked `run_*`

---

## Security Domain

`security_enforcement` absent → enabled. Phase 39 has **no file-input surface** (microdosimetry's CSV input is a different phase). Inputs are numeric sidebar widgets, bounded by `min_value` (device_sidebar.py). Only honest item: **CSV injection** — a leading `=`/`+`/`-`/`@` in a CSV cell can execute in Excel. But nothing user-controlled reaches the CSV as a free-text string: config values are numeric, `sim_type` is a fixed literal, column headers are hardcoded. Risk: LOW. Mitigation: if any string ever enters a cell, prefix-escape it; otherwise no action needed.

| ASVS                | Applies | Control                                              |
| ------------------- | ------- | ---------------------------------------------------- |
| V5 Input Validation | partial | Numeric widgets bounded by `min_value`; no free-text |
| V6 Cryptography     | no      | —                                                    |

---

## Environment Availability

| Dependency | Required By              | Available | Version                         | Fallback                                        |
| ---------- | ------------------------ | --------- | ------------------------------- | ----------------------------------------------- |
| devsim     | `run_*` solves           | ✓         | (loads UMFPACK 5.1)             | —                                               |
| streamlit  | all UI                   | ✓         | 1.55.0                          | —                                               |
| pandas     | CSV export               | ✓         | (installed)                     | —                                               |
| numpy      | arrays                   | ✓         | (installed)                     | —                                               |
| **plotly** | all charts (UI-03/04/05) | **✗**     | declared `>=5.0`, NOT installed | **none — MUST install via uv (Wave 0 blocker)** |
| matplotlib | (not used by Phase 39)   | ✓         | —                               | —                                               |

**Missing with no fallback:** `plotly` — hard blocker. Install via uv per project memory (`uv add plotly` / `uv sync`), not pip/venv.

---

## Assumptions Log

| #   | Claim                                                                | Section     | Risk if Wrong                                                                                                 |
| --- | -------------------------------------------------------------------- | ----------- | ------------------------------------------------------------------------------------------------------------- |
| A1  | CSV columns per page as listed in §4                                 | CSV Export  | Low — cosmetic; user may want different names/units                                                           |
| A2  | Commented `#` metadata header format                                 | §4          | Medium — affects re-readability (needs `comment='#'`); user may prefer separate params file or in-UI expander |
| A3  | Timestamp + `petringa.__version__` + full `asdict(config)` in header | §4          | Low — this is the audit-flagged traceability requirement; format is the only open part                        |
| A4  | Bonus plots (depletion width, net doping) offered opt-in             | Extra Plots | Low — additive, removable                                                                                     |
| A5  | 2D config → `st.warning` + `st.stop` (vs. disabling Run)             | Pitfall 1   | Low — UX choice; behavior (no crash) is required either way                                                   |
| A6  | Reference facades as `petringa.run_*` for mockability                | Validation  | Low — but changes page import style; verify in Wave 0 spike                                                   |

---

## Sources

### Primary (HIGH — read verbatim this session)

- `petringa/api/simulation.py` (877 lines) — `run_cv` (34-131), `run_field` (133-317), `run_cce` (320-411)
- `petringa/api/results.py` (20-46) — `SimResult`, `MeshData`
- `petringa/api/device.py` (25-115) — `DeviceConfig`, `build_device`
- `petringa/__init__.py` (8-36) — public exports
- `app/main.py`, `app/components/device_sidebar.py`, `app/workflows/{cv,cce,field_map,home}.py`, `tests/test_app_pages.py`
- `petringa/core/plotting.py` (matplotlib reference: 339-530)
- `pyproject.toml`, `.planning/config.json`, `.planning/ROADMAP.md` (450-480)

### Empirically verified this session

- plotly NOT installed (`ModuleNotFoundError`); pyproject declares `plotly>=5.0`
- `import devsim; import petringa` chain OK (petringa 5.0.0, UMFPACK loads)
- streamlit 1.55.0; AppTest accessor list confirms NO `plotly_chart`/`download_button` accessors
- `DeviceConfig` unhashable (`TypeError`); `dataclasses.asdict(DeviceConfig())` works
- `petringa.__version__ = "5.0.0"`

## Metadata

**Confidence breakdown:**

- API contract: HIGH — quoted verbatim from source
- Environment (plotly missing, unhashable config, AppTest accessors): HIGH — empirically verified
- CSV format / metadata / extra plots: MEDIUM (`[ASSUMED]`) — no CONTEXT.md, no in-repo precedent; concrete defaults recommended for planner
- Pitfalls (2D guard, caching): HIGH — grounded in verified contracts + behavior

**Research date:** 2026-07-11
**Valid until:** ~30 days (stable local codebase; only risk is upstream Streamlit accessor changes)
