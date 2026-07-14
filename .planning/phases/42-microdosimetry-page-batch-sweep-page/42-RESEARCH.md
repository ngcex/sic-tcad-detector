# Phase 42: Microdosimetry Page + Batch Sweep Page - Research

**Researched:** 2026-07-14
**Domain:** Streamlit page wiring (file upload + parametric sweep) over the existing `petringa` pure-Python/devsim facades
**Confidence:** HIGH

## Summary

Phase 42 fleshes out two already-registered placeholder pages — `app/workflows/microdosimetry.py` (FEAT-03) and `app/workflows/batch_sweep.py` (FEAT-04) — into the same Run → cache → render → download shape proven five times already (C-V, CCE, Field Map, Radiation Damage, Dark Current). No new packages, no new physics, no new API. Every underlying facade (`run_microdosimetry`, `ParametricSweep`, `run_cce`, etc.) already exists and is tested. This is pure UI-wiring against a frozen library contract. `[VERIFIED: codebase]`

The **microdosimetry page** is the simplest new surface: `run_microdosimetry` is a _pure data pipeline_ (no devsim, no solver, no convergence risk) that takes an MC-events CSV path and returns a `SimResult` with a 300-bin `y·d(y)` spectrum plus `y_F`/`y_D` scalar readouts. The only genuinely new mechanism is bridging Streamlit's `st.file_uploader` (which yields an in-memory `UploadedFile`) to the facade's `mc_csv_path: str` (a filesystem path) — solved with a server-generated `tempfile.NamedTemporaryFile` written from the uploaded bytes and removed in `finally`. There is **no path-traversal surface** (the user supplies bytes, not a path) and **no code-execution surface** (`load_mc_events_csv` is `pd.read_csv` only — no eval/pickle), consistent with threat T-37-02-V5. `[VERIFIED: codebase]`

The **batch sweep page** is a direct generalization of the already-shipped Dark Current page, which is _itself_ a `ParametricSweep(param="T", sim_fn=run_dark_current, ...)` page. The generalization: let the user (1) pick a swept `DeviceConfig` field from a **curated numeric-only selectbox** (not "all fields" — several fields break the 1D facades), (2) enter a numeric value list, (3) pick a simulation facade that returns an overlayable x/y curve, then render one Plotly trace per swept value and offer a **bulk CSV** with a run-identifier column. `run_cce` is the recommended default sim type because its default `v_stop=-40` stays inside the DD convergence envelope AND it _truncates gracefully_ (returns a shorter valid curve) instead of crashing — the safe demo path for success-criterion 4 (≥3 values render). `[VERIFIED: codebase + STATE.md]`

**Primary recommendation:** Copy the Dark Current page structure verbatim for batch_sweep (it is the canonical ParametricSweep-page precedent); add a `tempfile`-bridged `st.file_uploader` to the microdosimetry page; add two pure builders (`build_microdosimetry_figure`, `build_sweep_overlay_figure`) and one new bulk-CSV serializer to `app/components/results.py`; default the batch sweep demo to `run_cce` over `epi_thickness_um=[10,15,20]`.

<phase_requirements>

## Phase Requirements

| ID                    | Description                                                                                                                                                         | Research Support                                                                                                                                                                                                                                                                                                 |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FEAT-03               | User can upload a Geant4 MC CSV file on the microdosimetry page and see a y·d(y) vs log(y) lineal energy spectrum with y_F and y_D displayed                        | `run_microdosimetry` facade (pure, no devsim) confirmed working end-to-end against `data/synthetic_mc_events.csv`; y_F=17.23, y_D=53.22 keV/µm, 300 bins. File-upload bridge via `tempfile`. New `build_microdosimetry_figure` builder (`xaxis_type="log"`). y_F/y_D shown as `st.metric`/`st.caption` readouts. |
| FEAT-04               | User can configure and run a parametric sweep (select parameter, define value list, choose sim type), view overlaid results, and download all results as a bulk CSV | `ParametricSweep` already exists and is proven live on the Dark Current page. New `build_sweep_overlay_figure(list[SimResult], param, values)` builder + new `sweep_results_to_csv_bytes(...)` bulk serializer with a run-identifier column. Curated sweepable-field selectbox + facade selectbox.               |
| </phase_requirements> |

## Architectural Responsibility Map

| Capability                           | Primary Tier                                         | Secondary Tier | Rationale                                                                                                          |
| ------------------------------------ | ---------------------------------------------------- | -------------- | ------------------------------------------------------------------------------------------------------------------ |
| MC CSV upload + validation           | Frontend (Streamlit page)                            | —              | `st.file_uploader` is a UI widget; validation (schema/columns) happens page-side before calling the pure facade    |
| Microdosimetric spectrum computation | Library (`petringa.run_microdosimetry`)              | —              | Pure data pipeline already implemented; page never re-implements spectrum math                                     |
| File bytes → path bridge             | Frontend (Streamlit page)                            | —              | Temp-file lifecycle is a page concern; the facade contract is path-based and must not change                       |
| Parametric sweep orchestration       | Library (`petringa.ParametricSweep`)                 | —              | `.run()` clones config via `dataclasses.replace` and calls the facade N times; page must NOT hand-roll the loop    |
| Per-run TCAD solve                   | Library (facade `sim_fn`, e.g. `run_cce`)            | —              | Each facade builds + tears down its own devsim device; page is devsim-agnostic                                     |
| Overlay figure + bulk CSV            | Frontend (`app/components/results.py` pure builders) | —              | Pure `go.Figure`/`bytes` functions, no `st.*`, consumed by the page — matches the established `results.py` pattern |
| Result caching across reruns         | Frontend (`st.session_state`)                        | —              | Same cache-in-session_state pattern as all prior result pages                                                      |

## Standard Stack

### Core

| Library   | Version         | Purpose                                                                                           | Why Standard                                                                                                                  |
| --------- | --------------- | ------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| streamlit | 1.58.0          | Page UI, `st.file_uploader`, widgets, `st.session_state`, `st.plotly_chart`, `st.download_button` | Project-locked UI framework (STATE.md v5.0 decision); already the app runtime `[VERIFIED: uv run python -c import streamlit]` |
| plotly    | (installed)     | `go.Figure` builders (spectrum + overlay)                                                         | Established plot layer across all Phase 39–41 pages `[VERIFIED: codebase]`                                                    |
| pandas    | ≥2.0            | `to_csv` serialization; `load_mc_events_csv` uses `read_csv` internally                           | Already a runtime dep (STATE.md); used by existing `to_csv_bytes` `[VERIFIED: pyproject/STATE.md]`                            |
| numpy     | (installed)     | array assembly for sweep aggregation                                                              | Ubiquitous in facades `[VERIFIED: codebase]`                                                                                  |
| petringa  | 0.x (this repo) | `run_microdosimetry`, `ParametricSweep`, `run_cce`, `SimResult`, `DeviceConfig`                   | The library under UI — all facades already implemented and tested `[VERIFIED: codebase]`                                      |

**Standard library only (no install):** `tempfile` (upload bridge), `os`/`pathlib` (temp cleanup), `dataclasses.fields`/`asdict` (field introspection for the sweepable-param selectbox).

### Supporting

| Library                        | Version | Purpose                               | When to Use                                                                                                  |
| ------------------------------ | ------- | ------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `streamlit.testing.v1.AppTest` | 1.58.0  | AppTest harness for page verification | Every page test; **this phase's `file_uploader` IS drivable** (see Pitfall/testability finding) `[VERIFIED]` |

### Alternatives Considered

| Instead of                              | Could Use                              | Tradeoff                                                                                                                                                                                                                                                                |
| --------------------------------------- | -------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `tempfile.NamedTemporaryFile` + path    | `io.BytesIO` straight to `pd.read_csv` | **Rejected** — `run_microdosimetry`'s signature is `mc_csv_path: str` (a filesystem path). Passing a BytesIO would require changing the frozen library contract (no physics/API changes allowed in v5.0). Honor the contract; bridge at the UI. `[VERIFIED: signature]` |
| Curated field selectbox                 | Free-text field name                   | **Rejected** — an arbitrary string into `ParametricSweep.param` either `TypeError`s (unknown field) or, worse, selects a field that breaks the 1D facades (`half_width_um`, `doping_profile`). A selectbox of vetted numeric fields is safer and clearer.               |
| Text-input value list + `float()` parse | `st.data_editor` dynamic add-row       | Either works; text+`float()` is simpler, AppTest-drivable, and matches the numeric-widget precedent. `st.data_editor` is heavier and its AppTest support is less certain. Recommend text-input parse.                                                                   |

**Installation:**

```bash
# NONE — every dependency already present. No install step in any Phase 42 task.
```

**Version verification:** `streamlit 1.58.0` confirmed via `uv run python -c "import streamlit; print(streamlit.__version__)"`. All other libs already declared in `pyproject.toml` (STATE.md). `[VERIFIED]`

## Package Legitimacy Audit

> No external packages are installed by this phase. Every dependency (streamlit, plotly, pandas, numpy, petringa) is already present and vetted in prior phases; `tempfile`/`os`/`dataclasses` are Python standard library.

| Package | Registry | Age | Downloads | Source Repo | slopcheck | Disposition                |
| ------- | -------- | --- | --------- | ----------- | --------- | -------------------------- |
| (none)  | —        | —   | —         | —           | n/a       | No install step this phase |

**Packages removed due to slopcheck [SLOP] verdict:** none
**Packages flagged as suspicious [SUS]:** none

_No install step exists in this phase. Planner must NOT invent one._

## Architecture Patterns

### System Architecture Diagram

**Microdosimetry page (FEAT-03):**

```
[st.file_uploader]  --UploadedFile(bytes)-->  [page: validate columns]
        |                                              |
   (Run button)                                        v
        |                              [tempfile.NamedTemporaryFile(.csv, delete=False)]
        |                                    write uploaded.getvalue()
        v                                              |
   try/finally  ------------------------------------>  |  .name (server-generated path)
        |                                              v
        |                          petringa.run_microdosimetry(cfg, path, sv_t, sv_w)
        |                                              |  (PURE data pipeline, no devsim)
        |                                              v
        |                                SimResult(sim_type="microdosimetry",
        |                                    x=bin_centers, y=y·d(y),
        |                                    metadata={y_F, y_D, f_y, l_bar_um})
   finally: os.remove(path)                            |
        |                                              v
        v                          st.session_state["microdosimetry_result"] = result
   [build_microdosimetry_figure] --> st.plotly_chart (xaxis_type="log")
   [y_F / y_D readouts] ----------> st.metric / st.caption
```

**Batch sweep page (FEAT-04):**

```
[selectbox: swept field] [text_input: value list] [selectbox: sim type]
        |                        |                        |
        v                        v (parse: float() per token, reject non-numeric)
    param="epi_thickness_um"  values=[10,15,20]     sim_fn=petringa.run_cce
        |                        |                        |
        +------------------------+------------------------+
                                 |  (Run button, try/except RuntimeError)
                                 v
        petringa.ParametricSweep(base_config=cfg, param, values, sim_fn,
                                 sim_kwargs).run()   # real orchestration, NOT hand-rolled
                                 |
                                 v
                          list[SimResult]  (one per swept value; each builds+tears
                                 |          down its own devsim device sequentially)
              +------------------+------------------+
              v                                     v
  [build_sweep_overlay_figure(results,   [sweep_results_to_csv_bytes(results,
      param, values)] -> 1 trace/value       param, values)] -> ONE csv, run-id column
              |                                     |
              v                                     v
        st.plotly_chart                    st.download_button("Download all results as CSV")
```

### Recommended Project Structure

```
app/
├── workflows/
│   ├── microdosimetry.py   # flesh out placeholder: uploader + tempfile bridge + run + render
│   └── batch_sweep.py      # flesh out placeholder: selectbox×2 + value list + ParametricSweep + overlay
├── components/
│   └── results.py          # ADD: build_microdosimetry_figure, build_sweep_overlay_figure,
│                           #      sweep_results_to_csv_bytes  (all pure, no st.*)
tests/
├── test_app_microdosimetry_page.py   # NEW — uses AppTest file_uploader.set_value/upload
└── test_app_batch_sweep_page.py      # NEW — monkeypatch sim_fn, real ParametricSweep.run()
data/
└── synthetic_mc_events.csv            # EXISTING fixture (2000 events, 11116 rows) — reuse for tests
```

### Pattern 1: File-upload → path bridge (microdosimetry page)

**What:** Bridge Streamlit's in-memory `UploadedFile` to the facade's `mc_csv_path: str` without changing the frozen API.
**When to use:** The only file-input surface in this phase.
**Example:**

```python
# Source: derived from run_microdosimetry signature (petringa/api/simulation.py:843)
#         + Streamlit st.file_uploader contract (v1.58). Pattern verified safe: no
#         path traversal (server-generated temp path), no eval (pd.read_csv only).
import os
import tempfile
import streamlit as st
import petringa

uploaded = st.file_uploader("Upload MC events CSV", type=["csv"], key="micro_csv")
sv_thickness = st.number_input("Sensitive-volume thickness (µm)", value=10.0, key="micro_sv_t")
sv_width = st.number_input("Sensitive-volume width (µm)", value=150.0, key="micro_sv_w")

if st.button("Run simulation"):
    if uploaded is None:
        st.warning("Upload an MC events CSV to run.")
        st.stop()
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
            tmp.write(uploaded.getvalue())   # UploadedFile bytes
            tmp_path = tmp.name              # server-generated path — NOT user-controlled
        result = petringa.run_microdosimetry(
            cfg, mc_csv_path=tmp_path,
            sv_thickness_um=sv_thickness, sv_width_um=sv_width,
        )
        st.session_state["microdosimetry_result"] = result
    except (ValueError, KeyError) as e:      # malformed CSV -> load_mc_events_csv raises
        st.error(f"Could not parse the uploaded CSV: {e}")
    finally:
        if tmp_path is not None and os.path.exists(tmp_path):
            os.remove(tmp_path)
```

`petringa.run_microdosimetry` MUST be referenced as a **module attribute** (`import petringa; petringa.run_microdosimetry(...)`), never `from petringa import run_microdosimetry`, so tests can `monkeypatch.setattr(petringa, "run_microdosimetry", fake)` — the mockability seam locked in 39-01 and reused by every page. `[VERIFIED: STATE.md decision + all 5 prior pages]`

### Pattern 2: ParametricSweep page (batch sweep) — copy Dark Current verbatim

**What:** Drive `petringa.ParametricSweep(...).run()` and render one trace per swept value.
**When to use:** The batch sweep page IS the general case of the dark current page.
**Example:**

```python
# Source: app/workflows/dark_current.py (41-03, the canonical ParametricSweep page)
import petringa
from app.components.results import build_sweep_overlay_figure, sweep_results_to_csv_bytes

SWEEPABLE_FIELDS = [   # CURATED — numeric, 1D-facade-safe (see "Curated sweep fields" below)
    "epi_thickness_um", "substrate_thickness_um", "N_A",
    "N_D_junction", "N_D_bulk", "L_transition_um", "T", "area_cm2",
]
SIM_FACADES = {        # label -> module-attribute facade; each returns an overlayable x/y curve
    "CCE vs bias (run_cce)": "run_cce",
    "C-V (run_cv)": "run_cv",
    "CCE vs temperature (run_temperature_sweep)": "run_temperature_sweep",
    # run_field EXCLUDED: single-target ramp, raises with no partial fallback (see Pitfall 3)
}

param = st.selectbox("Sweep parameter", SWEEPABLE_FIELDS, key="sweep_param")
values_raw = st.text_input("Values (comma-separated)", value="10, 15, 20", key="sweep_values")
sim_label = st.selectbox("Simulation type", list(SIM_FACADES), key="sweep_sim")

if st.button("Run simulation"):
    try:
        values = [float(v.strip()) for v in values_raw.split(",") if v.strip()]
    except ValueError:
        st.error("Values must be a comma-separated list of numbers.")
        st.stop()
    if len(values) < 1:
        st.warning("Enter at least one value.")
        st.stop()
    sim_fn = getattr(petringa, SIM_FACADES[sim_label])
    try:
        results = petringa.ParametricSweep(   # REAL .run() — never hand-roll the loop
            base_config=cfg, param=param, values=values, sim_fn=sim_fn,
        ).run()
        st.session_state["sweep_results"] = results
        st.session_state["sweep_param"] = param
        st.session_state["sweep_values"] = values
    except RuntimeError as e:
        st.error(f"Simulation failed to converge: {e}\n\nTry a shallower value range.")
```

### Anti-Patterns to Avoid

- **Hand-rolling the sweep loop** instead of calling `ParametricSweep.run()`: defeats success-criterion 4's "uses ParametricSweep under the hood" and the security-hardened `dataclasses.replace` cloning. The Dark Current test explicitly asserts the _real_ `.run()` runs (only the facade is monkeypatched).
- **Passing `st.file_uploader`'s BytesIO straight into a path API**: doesn't work (`run_microdosimetry` opens a path); tempting but wrong.
- **`eval`/`exec` for the value list** or **arbitrary attribute injection** into config: banned by threat T-37-03-V5. Use `float()` per token and `ParametricSweep`'s `dataclasses.replace` (which `TypeError`s on unknown fields — do not bypass).
- **Leaving the temp file undeleted**: always `os.remove` in `finally`.
- **Adding `sv_thickness_um`/`sv_width_um` to the sidebar DeviceConfig**: they are `run_microdosimetry` params, NOT `DeviceConfig` fields — put their widgets on the microdosimetry page only.

## Don't Hand-Roll

| Problem                               | Don't Build                                              | Use Instead                                                                                                     | Why                                                                                                            |
| ------------------------------------- | -------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| Sweeping a config field over N values | A manual `for value in values:` loop that mutates config | `petringa.ParametricSweep(...).run()`                                                                           | Success-criterion 4 requires it; `dataclasses.replace` cloning is the security-hardened path (no setattr/eval) |
| MC CSV parsing/validation             | Custom CSV reader                                        | `petringa.core.mc_coupling.load_mc_events_csv` (called inside `run_microdosimetry`)                             | Already handles the column map + unit conversion; page only needs to catch its exceptions                      |
| Lineal-energy spectrum math           | Reimplement ICRU-36 binning                              | `run_microdosimetry`                                                                                            | Pure pipeline already implemented, tested, and frozen                                                          |
| Config cloning                        | `copy.deepcopy` + setattr                                | `dataclasses.replace` (inside ParametricSweep)                                                                  | Rejects unknown fields with `TypeError`; no attribute injection                                                |
| CSV metadata header / timestamp       | New header format                                        | Mirror existing `to_csv_bytes` header convention (`#` comment lines, ISO-8601, software_version, device fields) | Consistency with 5 existing CSV exports                                                                        |

**Key insight:** Every hard part of this phase is already solved in the library. The phase is 90% widget wiring + 3 pure `results.py` functions. The only _new_ logic is (a) the tempfile bridge and (b) the bulk-CSV run-identifier column — both trivial and testable.

## Common Pitfalls

### Pitfall 1: Assuming AppTest cannot drive file upload

**What goes wrong:** Prior pages (39/40/41) documented that AppTest 1.55 has no `plotly_chart`/`download_button` accessor, leading to an assumption that file upload is also untestable and the happy path can't be exercised.
**Why it happens:** Cargo-culting the "AppTest can't reach the chart" caveat to a different widget.
**How to avoid:** AppTest **1.58** exposes `at.file_uploader[...]` with `.set_value((name, bytes, mime))` and `.upload(...)` — **file upload IS drivable in tests.** `[VERIFIED: streamlit.testing.v1.element_tree.FileUploader has set_value/upload/value/clear]`. Use `data/synthetic_mc_events.csv`'s bytes as the injected upload to test the full Run→cache path. (Chart/download-button _rendering_ still isn't AppTest-accessible; assert on `at.session_state`, `at.exception`, `at.error`, `at.warning` as prior pages do.)
**Warning signs:** A plan that says "happy path can only be verified in the browser" for the microdosimetry page — that's stale; the upload path is unit-testable.

### Pitfall 2: `to_csv_bytes` has no `microdosimetry` branch; bulk CSV is NOT a `to_csv_bytes` branch

**What goes wrong:** Trying to route the batch-sweep bulk export through the existing `to_csv_bytes(result)` (which takes ONE `SimResult` and dispatches on `sim_type`).
**Why it happens:** `to_csv_bytes` is the obvious CSV helper, but it is single-result and has branches only for `cv|cce|field|damage|dark_current` — no `microdosimetry`, and no multi-result mode.
**How to avoid:** Add a **separate** `sweep_results_to_csv_bytes(results: list[SimResult], param: str, values: list) -> bytes` that emits ONE CSV with a leading run-identifier column (the swept `param` value) so all N curves live in one file (success-criterion 3). Keep the existing `to_csv_bytes` untouched. For the microdosimetry single-result download (not a stated criterion — discretion): either add a `microdosimetry` branch to `to_csv_bytes` OR skip a single-result download entirely. Flag for UI-SPEC.
**Warning signs:** `to_csv_bytes(result)` called with a `list` argument, or a `KeyError`/`ValueError: unknown sim_type` at export time.

### Pitfall 3: Choosing a demo sweep that hits the DD convergence wall

**What goes wrong:** The success-criterion-4 demo (≥3 values render) crashes because the default facade/config combination hits the known `ramp_bias` non-convergence at deep bias (STATE.md Phase 39: `run_cce` fails ~V≈60.5, `run_field` fails ~V=66 for the default config).
**Why it happens:** Deep-depletion punch-through of the default 10µm epi.
**How to avoid:**

- **Default sim type = `run_cce`** (default `v_stop=-40` stays inside the envelope AND truncates gracefully — returns a shorter valid curve with `truncated=True` rather than raising).
- **Default sweep = `param="epi_thickness_um", values=[10, 15, 20]`** — sweeping epi _thicker_ moves punch-through to deeper bias (safer). Thinner epi punches through at _shallower_ bias; but even then `run_cce` truncates rather than crashing, so the overlay still renders.
- **Exclude `run_field` from the sim-type selectbox** — it is a single-target ramp with NO partial-result fallback: a beyond-envelope value raises `RuntimeError` outright (no truncated curve to plot). `run_cce`/`run_cv` sweeps degrade; `run_field` does not.
- Keep the `try/except RuntimeError` → `st.error` guard (mirrors all 3 Phase 39 pages) so any residual non-convergence shows a friendly message, not a traceback.
  **Warning signs:** `run_field` in the sim-type list; a default value range that pushes bias past ~-50V for the default config. **[ASSUMED — recommend a planner Wave-0 live-devsim spike to confirm `[10,15,20]`+`run_cce` renders 3 clean curves, exactly as Phase 41 spiked its defaults.]**

### Pitfall 4: Offering un-sweepable / 1D-breaking config fields

**What goes wrong:** A user picks `half_width_um` (switches 1D→2D → `run_cce`/`run_cv` raise `NotImplementedError`), `doping_profile` (a `str`, not numeric), or `N_D` (`None` unless `doping_profile="uniform"`) and the sweep crashes.
**Why it happens:** `dataclasses.fields(DeviceConfig)` naively yields all 12 fields including non-numeric and dimension-switching ones.
**How to avoid:** Present a **curated selectbox** of the vetted numeric, 1D-safe fields only: `epi_thickness_um`, `substrate_thickness_um`, `N_A`, `N_D_junction`, `N_D_bulk`, `L_transition_um`, `T`, `area_cm2`. Exclude `half_width_um` (1D→2D), `doping_profile` (str), `N_D` (None-by-default / uniform-only).
**Warning signs:** `NotImplementedError` from `run_cce`/`run_cv`, or `TypeError` on non-numeric replace values.

### Pitfall 5: Empty-state guard on `device_config` for the microdosimetry page

**What goes wrong:** The microdosimetry spectrum does **not** depend on `device_config` (the config is carried on `SimResult` for provenance only; the spectrum depends solely on MC events + SV geometry). The existing placeholder still guards on `cfg is None` and `st.stop()`s — keeping that guard forces the user to configure a device that has no effect on the result.
**Why it happens:** Copy-paste of the standard empty-state guard from pages where config DOES matter.
**How to avoid:** This is a **UX decision for the UI-SPEC**, not a silent implementation choice. Keeping the guard is _consistent_ with every other page (config is still recorded as provenance and shown in the CSV header) but adds friction. Surface both options in the UI-SPEC; do not decide here. **[ASSUMED — flag for UI-SPEC/discuss.]**

## Runtime State Inventory

> Not applicable — this is a greenfield UI-wiring phase (two placeholder pages fleshed out), not a rename/refactor/migration. No stored data, live-service config, OS-registered state, secrets, or build artifacts are modified. The temp file created during upload is created _and deleted within the same request_ (in `finally`), so it leaves no persistent runtime state.

## Code Examples

### Concrete valid MC CSV row + fixture facts (for test fixture / in-UI format hint)

```
event_id,x,y,z,edep
0,0.0035860123023206346,0.0,0.0003046023557741905,353.3627594257208
```

`[VERIFIED: head -6 data/synthetic_mc_events.csv]`

- **Default column map** (`load_mc_events_csv`): `event_id`, `x`(→x_cm), `y`(→y_cm), `z`(→z_cm), `edep`(→edep_keV). Default `pos_unit="cm"`, `energy_unit="keV"`. `[VERIFIED: petringa/core/mc_coupling.py:107-158]`
- **Fixture stats:** 11116 rows, 2000 unique `event_id`s, edep range ~0.0002–2255 keV. `[VERIFIED: pandas describe]`
- **Facade output for this fixture** (`sv_thickness_um=10, sv_width_um=150`): `sim_type="microdosimetry"`, `x`/`y` length **300**, `y_F=17.23 keV/µm`, `y_D=53.22 keV/µm`, `l_bar_um=20.0`, x-range 0.0102 – 9772 keV/µm. `[VERIFIED: live uv run]`
- **In-UI format hint copy** (suggested): "CSV columns: `event_id, x, y, z, edep` — positions in cm, energy deposit in keV. One row per MC step; steps are summed per `event_id`."

### `build_microdosimetry_figure` (new pure builder)

```python
# Source: pattern derived from existing build_cce_figure (app/components/results.py)
def build_microdosimetry_figure(result: SimResult) -> go.Figure:
    """y·d(y) vs lineal energy y (LOG x-axis per ICRU-36 convention)."""
    fig = go.Figure(data=go.Scatter(x=result.x, y=result.y, mode="lines"))
    fig.update_layout(
        title="Microdosimetric Spectrum",
        xaxis_title="Lineal energy y (keV/µm)",
        xaxis_type="log",              # criterion says "vs log(y)" — x MUST be log
        yaxis_title="y · d(y)",
    )
    return fig
```

**Note:** x-range spans 0.01–9772 keV/µm → `xaxis_type="log"` is mandatory, not optional.

### `build_sweep_overlay_figure` (new pure builder)

```python
def build_sweep_overlay_figure(results: list[SimResult], param: str, values: list) -> go.Figure:
    """One trace per swept value; legend keyed by the parameter value.
    Reuse the Phase 41 data-viz qualitative palette for distinguishable traces."""
    PALETTE = ["#1F6FEB", "#D32F2F", "#2E7D32", "#1A1A1A", "#9AA0A6"]  # cycle if >5
    fig = go.Figure()
    for i, (val, res) in enumerate(zip(values, results)):
        # run_cce/run_cv use |bias| convention on the CCE page; keep the axis
        # consistent with the chosen facade (see per-facade x-axis note below).
        fig.add_trace(go.Scatter(
            x=res.x, y=res.y, mode="lines+markers",
            name=f"{param}={val}", line=dict(color=PALETTE[i % len(PALETTE)]),
        ))
    fig.update_layout(title=f"Parametric Sweep: {param}",
                      xaxis_title="(facade x-axis)", yaxis_title="(facade y-axis)")
    return fig
```

### `sweep_results_to_csv_bytes` (new bulk CSV — run-identifier column)

```python
def sweep_results_to_csv_bytes(results: list[SimResult], param: str, values: list) -> bytes:
    """ALL sweep runs in ONE CSV. Leading `<param>` column identifies each run.
    Mirror the existing to_csv_bytes header convention (# comment lines)."""
    frames = []
    for val, res in zip(values, results):
        df = pd.DataFrame({param: val, "x": res.x, "y": res.y})
        frames.append(df)
    combined = pd.concat(frames, ignore_index=True)
    header = [
        f"# petringa SiC TCAD Simulator — parametric sweep ({param})",
        f"# software_version: {petringa.__version__}",
        f"# generated: {datetime.now(timezone.utc).isoformat()}",
        f"# swept_values: {values}",
    ]
    return ("\n".join(header) + "\n" + combined.to_csv(index=False)).encode("utf-8")
```

## State of the Art

| Old Approach                              | Current Approach                                               | When Changed               | Impact                                                                                                                              |
| ----------------------------------------- | -------------------------------------------------------------- | -------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| AppTest can't reach result widgets (1.55) | AppTest **1.58** `file_uploader` supports `set_value`/`upload` | streamlit 1.58 (this repo) | The microdosimetry happy path is now unit-testable end-to-end via injected fixture bytes — no browser-only caveat needed for upload |

**Deprecated/outdated:** none relevant.

## Assumptions Log

| #   | Claim                                                                                                                                           | Section                    | Risk if Wrong                                                                                                       |
| --- | ----------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| A1  | `run_cce` + `epi_thickness_um=[10,15,20]` renders 3 clean/gracefully-truncated overlay curves for the default config                            | Pitfall 3                  | Demo (criterion 4) hits a convergence wall; a Wave-0 live spike (as Phase 41 did) confirms cheaply                  |
| A2  | Text-input + `float()` parse for the value list is AppTest-drivable and preferable to `st.data_editor`                                          | Standard Stack / Pattern 2 | If AppTest can't drive it, fall back to `st.number_input`-per-row; low risk (text_input is a proven AppTest widget) |
| A3  | Keeping the `device_config` empty-state guard on the microdosimetry page is the right UX (vs. removing it since spectrum is config-independent) | Pitfall 5                  | Extra user friction, or inconsistency with other pages — resolve in UI-SPEC, not code                               |
| A4  | The bulk-CSV run-identifier column should be the swept `param` value (e.g. a column literally named `epi_thickness_um`)                         | Pitfall 2 / code example   | Only a schema-shape choice; any clear identifier column satisfies criterion 3                                       |
| A5  | Overlaying `run_cce`/`run_cv` curves across configs shares a coherent x-axis (bias) so traces are comparable                                    | Architecture / builder     | True within one facade; DO NOT mix facades in one sweep (the page picks one sim type per run, so this holds)        |

**If this table were empty:** it is not — A1 and A3 in particular should be confirmed (A1 via a live spike, A3 via the UI-SPEC) before locking the plan.

## Open Questions

1. **Should the microdosimetry page keep the `device_config` empty-state guard even though the spectrum is config-independent?**
   - What we know: config is provenance-only on the `SimResult`; the spectrum depends only on MC events + SV geometry.
   - What's unclear: whether forcing device config adds unwanted friction vs. cross-page consistency.
   - Recommendation: surface both in the UI-SPEC; default to keeping the guard for consistency unless the UI-SPEC decides otherwise. (A3)

2. **Single-result CSV download on the microdosimetry page — include it or not?**
   - What we know: not a stated success criterion; `to_csv_bytes` currently lacks a `microdosimetry` branch.
   - Recommendation: discretion — if included, add a `microdosimetry` branch to `to_csv_bytes` (columns: `y_keV_per_um`, `y_times_d_y`, plus `y_F`/`y_D`/`l_bar_um` in the `#` header); else omit. Flag for UI-SPEC.

3. **Which sim facades to expose in the batch-sweep sim-type selectbox beyond `run_cce`?**
   - What we know: `run_cce`, `run_cv`, `run_temperature_sweep`, `run_radiation_damage`, `run_flash_recombination`, `run_transient` all return an overlayable x/y curve. `run_field` returns a single-bias depth profile with no partial-failure fallback (exclude). `run_microdosimetry` needs an mc_csv_path (doesn't fit the config-sweep model — exclude).
   - Recommendation: start with `run_cce` (default), `run_cv`, `run_temperature_sweep`; keep the list curated and explicitly EXCLUDE `run_field` and `run_microdosimetry`. Confirm in UI-SPEC.

## Environment Availability

| Dependency                     | Required By                        | Available                                                                                        | Version     | Fallback |
| ------------------------------ | ---------------------------------- | ------------------------------------------------------------------------------------------------ | ----------- | -------- |
| streamlit                      | both pages                         | ✓                                                                                                | 1.58.0      | —        |
| plotly                         | figure builders                    | ✓                                                                                                | installed   | —        |
| pandas                         | CSV load/serialize                 | ✓                                                                                                | ≥2.0        | —        |
| numpy                          | array assembly                     | ✓                                                                                                | installed   | —        |
| devsim                         | batch-sweep facades (run_cce etc.) | ✓ (with a MISSING-DLL warning for `libopenblas.dylib`, but `liblapack` loads and solves succeed) | installed   | —        |
| `data/synthetic_mc_events.csv` | microdosimetry test fixture        | ✓                                                                                                | 2000 events | —        |

**Missing dependencies with no fallback:** none.
**Missing dependencies with fallback:** none — the microdosimetry page (pure pipeline) has zero devsim dependency; the batch-sweep page uses the already-working facades.

## Validation Architecture

### Test Framework

| Property           | Value                                                                                                   |
| ------------------ | ------------------------------------------------------------------------------------------------------- |
| Framework          | pytest (with `streamlit.testing.v1.AppTest`)                                                            |
| Config file        | `pyproject.toml` `[project.optional-dependencies].dev` (pytest); no separate pytest.ini                 |
| Quick run command  | `uv run pytest tests/test_app_microdosimetry_page.py tests/test_app_batch_sweep_page.py -q`             |
| Full suite command | Per-file/per-class isolation (bare `pytest -q` is unsatisfiable — devsim resource exhaustion, STATE.md) |

### Phase Requirements → Test Map

| Req ID  | Behavior                                                                                  | Test Type                                                                                   | Automated Command                                        | File Exists?      |
| ------- | ----------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- | -------------------------------------------------------- | ----------------- |
| FEAT-03 | Upload CSV → Run → spectrum cached; y_F/y_D in session                                    | AppTest (inject fixture bytes via `file_uploader.upload`, monkeypatch `run_microdosimetry`) | `uv run pytest tests/test_app_microdosimetry_page.py -q` | ❌ Wave 0         |
| FEAT-03 | Empty state / no-file guard                                                               | AppTest                                                                                     | same file                                                | ❌ Wave 0         |
| FEAT-03 | Malformed CSV → `st.error`, no crash                                                      | AppTest (upload junk bytes)                                                                 | same file                                                | ❌ Wave 0         |
| FEAT-03 | `run_microdosimetry` pure-pipeline contract (x/y len, y_D≥y_F)                            | unit (existing)                                                                             | `uv run pytest tests/test_api_microdosimetry.py -q`      | ✅ exists         |
| FEAT-04 | Select param + values + sim type → Run → `ParametricSweep.run()` caches `list[SimResult]` | AppTest (monkeypatch `run_cce`, real `ParametricSweep`)                                     | `uv run pytest tests/test_app_batch_sweep_page.py -q`    | ❌ Wave 0         |
| FEAT-04 | ≥3 values produce ≥3 cached results                                                       | AppTest                                                                                     | same file                                                | ❌ Wave 0         |
| FEAT-04 | Bad value list → `st.error`, no crash                                                     | AppTest                                                                                     | same file                                                | ❌ Wave 0         |
| FEAT-04 | Bulk CSV serializer emits one CSV with run-id column for N results                        | unit (pure, no Streamlit)                                                                   | `uv run pytest tests/test_app_csv_export.py -q` (extend) | ⚠ extend existing |
| FEAT-04 | RuntimeError from facade → `st.error`, no crash                                           | AppTest                                                                                     | same file                                                | ❌ Wave 0         |

### Sampling Rate

- **Per task commit:** `uv run pytest tests/test_app_microdosimetry_page.py -q` (or the batch-sweep file for that task).
- **Per wave merge:** both new page test files + `tests/test_app_csv_export.py` + `tests/test_api_microdosimetry.py`.
- **Phase gate:** the four files above green (per-file isolation) before `/gsd:verify-work`.

### Wave 0 Gaps

- [ ] `tests/test_app_microdosimetry_page.py` — FEAT-03 (upload via `file_uploader.upload`, empty guard, malformed CSV, happy-path cache). Reuse `data/synthetic_mc_events.csv` bytes.
- [ ] `tests/test_app_batch_sweep_page.py` — FEAT-04 (monkeypatch `run_cce`, assert real `ParametricSweep.run()` executes, ≥3 results, bad-list guard, RuntimeError guard).
- [ ] Extend `tests/test_app_csv_export.py` (or a new pure test) — `sweep_results_to_csv_bytes` shape: one CSV, N runs, run-id column present.
- [ ] (Optional, planner Wave-0 spike) live-devsim confirm `run_cce`+`epi_thickness_um=[10,15,20]` renders 3 curves (A1).

## Security Domain

`security_enforcement` is absent from `.planning/config.json` — treat as enabled per the absent-key convention.

### Applicable ASVS Categories

| ASVS Category         | Applies | Standard Control                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| --------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| V2 Authentication     | no      | Local/shared-lab Streamlit app, no auth layer (Phase 38, out of scope)                                                                                                                                                                                                                                                                                                                                                                                                   |
| V3 Session Management | no      | `st.session_state` is Streamlit's in-process mechanism, not a security boundary                                                                                                                                                                                                                                                                                                                                                                                          |
| V4 Access Control     | no      | Single-user local tool, no roles                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| V5 Input Validation   | **yes** | **Two new input surfaces this phase:** (1) MC CSV upload — restrict `st.file_uploader(type=["csv"])`, catch `load_mc_events_csv`'s `pd.read_csv` exceptions; NO eval/pickle (it is `read_csv` only, matching T-37-02-V5). (2) Batch-sweep value list — parse with `float()` per token, reject non-numeric with `st.error`; NEVER `eval`/`exec`. Field selection is a curated selectbox (not free text) so `ParametricSweep.param` is always a real `DeviceConfig` field. |
| V6 Cryptography       | no      | No crypto operations                                                                                                                                                                                                                                                                                                                                                                                                                                                     |

### Known Threat Patterns for {Streamlit + pandas + devsim}

| Pattern                                                   | STRIDE                | Standard Mitigation                                                                                                                                                                                         |
| --------------------------------------------------------- | --------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Path traversal via uploaded file path                     | Tampering             | **None possible** — user uploads _bytes_, the temp path is server-generated by `tempfile.NamedTemporaryFile`; the user never supplies a path. `os.remove` in `finally`.                                     |
| Code execution from CSV content                           | Tampering / Elevation | `load_mc_events_csv` uses `pd.read_csv` only — no eval, no pickle, no binary deserialization (T-37-02-V5). CSV bytes are pure data.                                                                         |
| Arbitrary attribute injection into config via sweep param | Tampering             | `ParametricSweep` clones via `dataclasses.replace`, which `TypeError`s on unknown fields — no setattr/eval (T-37-03-V5). Curated selectbox further constrains `param` to vetted fields. Do NOT weaken this. |
| Malformed value list crashing the sweep                   | DoS (local)           | `float()` per token inside `try/except ValueError` → `st.error`; require ≥1 parsed value before running.                                                                                                    |
| Degenerate SV geometry / `n`-like inputs into devsim      | DoS (local)           | Streamlit `number_input` type coercion + the existing `try/except RuntimeError` on the facade call (mirrors all Phase 39–41 pages).                                                                         |
| Temp-file leak (disk exhaustion over many uploads)        | DoS (local)           | `os.remove` in `finally` on every run guarantees single-request lifetime; no accumulation.                                                                                                                  |

## Sources

### Primary (HIGH confidence)

- `petringa/api/simulation.py` (run_microdosimetry:843, run_cce:337, run_field:144, run_cv:34, run_dark_current:665) — facade signatures, sweep suitability, convergence-envelope docstrings `[VERIFIED: Read]`
- `petringa/api/sweep.py` — `ParametricSweep` contract, `dataclasses.replace` security note `[VERIFIED: Read]`
- `petringa/api/results.py` — `SimResult`/`MeshData` dataclasses `[VERIFIED: Read]`
- `petringa/api/device.py` — `DeviceConfig` full field list + types (sweepable-field curation) `[VERIFIED: Read]`
- `petringa/core/mc_coupling.py` (load_mc_events_csv:107) — CSV column map, units, `pd.read_csv` only `[VERIFIED: Read]`
- `petringa/core/microdosimetry.py` (lineal_energy_spectrum:117) — spectrum output keys `[VERIFIED: Read]`
- `app/components/results.py` — existing pure-builder + `to_csv_bytes` conventions (sim_type branches, header format) `[VERIFIED: Read]`
- `app/workflows/dark_current.py` (41-03) — canonical ParametricSweep-page precedent `[VERIFIED: Read]`
- `app/workflows/microdosimetry.py` / `batch_sweep.py` — current placeholders `[VERIFIED: Read]`
- `app/main.py` — both pages already registered in `st.navigation` `[VERIFIED: Read]`
- `data/synthetic_mc_events.csv` — fixture (2000 events, 11116 rows, columns event_id/x/y/z/edep) `[VERIFIED: head + pandas]`
- Live run of `run_microdosimetry` on the fixture — y_F=17.23, y_D=53.22, 300 bins `[VERIFIED: uv run]`
- `streamlit.testing.v1.element_tree.FileUploader` — `set_value`/`upload`/`value`/`clear` present in 1.58 `[VERIFIED: uv run inspect]`
- `.planning/STATE.md` — devsim resource exhaustion, ramp_bias convergence wall, page-import mockability decision `[VERIFIED: Read]`
- `.planning/phases/41-.../41-RESEARCH.md` Security Domain, `41-UI-SPEC.md`, `41-02-SUMMARY.md` — precedent for banner/guard/palette/test structure `[VERIFIED: Read]`

### Secondary (MEDIUM confidence)

- Streamlit `st.file_uploader` general contract (`type=`, `UploadedFile.getvalue()`) — inferred from installed version behavior + widget accessor `[CITED: streamlit 1.58 installed]`

### Tertiary (LOW confidence)

- none — every load-bearing claim is grounded in a codebase read or a live run.

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH — every dependency present and version-verified; zero new installs.
- Architecture: HIGH — batch sweep is a direct generalization of the shipped Dark Current page; microdosimetry facade run live.
- Pitfalls: HIGH — each grounded in a codebase read or STATE.md finding; the two [ASSUMED] items (A1 demo-convergence, A3 empty-state UX) are explicitly flagged for a Wave-0 spike / UI-SPEC.

**Research date:** 2026-07-14
**Valid until:** 2026-08-13 (stable — internal codebase contract, no fast-moving external deps)
