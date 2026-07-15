# Phase 38: Streamlit Shell + Device Config Page - Research

**Researched:** 2026-07-10
**Domain:** Streamlit multi-page app architecture, session_state persistence, dataclass-driven form UI
**Confidence:** HIGH

## Summary

Phase 38 establishes the Streamlit shell for the etna SiC TCAD simulator: a multi-page app whose sidebar collects **all** `DeviceConfig` fields into a single config object that persists across page navigation. It runs **no** simulations — every page is a placeholder/empty-state page. Because Phases 39-43 build on whatever navigation and session_state pattern this phase sets, the foundational contract must be right.

Three facts drive every recommendation below, all verified this session against official Streamlit docs (v1.58.0, the installed version) and the actual source code:

1. **`st.navigation` + `st.Page` is the current recommended multi-page API** and it is available in the installed Streamlit 1.58.0. Critically, with `st.navigation` the sidebar defined in the **entry script** runs on every page automatically — which directly satisfies UI-02 ("form on any page") and UI-07 ("persist across nav"). The legacy `pages/` directory forces per-page sidebar duplication. `[CITED: docs.streamlit.io/develop/api-reference/navigation/st.navigation]`
2. **Keyed widget state is garbage-collected when the widget is not rendered on a run** (confirmed verbatim in the docs). Therefore the assembled `DeviceConfig` must be stored under a **single non-widget key** (recommend `st.session_state["device_config"]`), not read directly from widget keys. `[CITED: docs.streamlit.io/develop/concepts/architecture/widget-behavior]`
3. **`DeviceConfig` has two mode switches that gate other fields** — `half_width_um` (None=1D / float=2D) and `doping_profile` ("graded" uses N_D_junction/N_D_bulk/L_transition_um; "uniform" uses N_D). `st.form` batches inputs and defers reruns until submit, so a mode selectbox **inside** `st.form` cannot toggle conditional field visibility before submit (confirmed in docs). The two mode selectors must be **reactive** (outside any `st.form`). `[CITED: docs.streamlit.io/develop/api-reference/execution-flow/st.form]`

**Primary recommendation:** Build `app/main.py` as an `st.navigation` entry script that renders the full device-config sidebar (mode selectors reactive/outside a form; conditional fields rendered per-mode), assembles a `DeviceConfig`, and stores it under `st.session_state["device_config"]`. Register one placeholder page per downstream workflow. Each placeholder guards on `"device_config" not in st.session_state` with an `st.info` prompt. Do **not** add `st.cache_resource`/`st.cache_data` in this phase — it is not warranted (see Don't Hand-Roll / State of the Art).

## Architectural Responsibility Map

| Capability                         | Primary Tier                               | Secondary Tier                      | Rationale                                                                                                             |
| ---------------------------------- | ------------------------------------------ | ----------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| Multi-page navigation              | Streamlit runtime (entry script)           | —                                   | `st.navigation` in `app/main.py` owns page registration and routing                                                   |
| Device config form controls        | Streamlit UI (sidebar, entry script)       | —                                   | Entry-script sidebar renders on every page, so config is editable "on any page" (UI-02) `[CITED: st.navigation docs]` |
| Config persistence across nav      | `st.session_state` (single non-widget key) | —                                   | session_state persists across pages within a session `[CITED: session-state docs]`                                    |
| DeviceConfig assembly / validation | App layer (`app/` helper)                  | `etna.DeviceConfig` (dataclass) | The dataclass is the contract; the app marshals form values into it                                                   |
| Empty-state guard on pages         | Streamlit UI (each page)                   | —                                   | Pages check for `device_config` and prompt if absent (success criterion 4)                                            |
| Simulation execution               | **NOT this phase**                         | `etna` facades                  | Phase 38 runs no simulations; facades are used starting Phase 39                                                      |

## Standard Stack

### Core

| Library   | Version                               | Purpose                                           | Why Standard                                                              |
| --------- | ------------------------------------- | ------------------------------------------------- | ------------------------------------------------------------------------- |
| streamlit | 1.58.0 (installed; declared `>=1.30`) | Multi-page app shell, sidebar form, session_state | Project-mandated UI framework (pyproject.toml) `[VERIFIED: pip list]`     |
| etna  | local editable install                | `DeviceConfig` dataclass consumed by the form     | The config contract this phase exposes `[VERIFIED: etna/__init__.py]` |

### Supporting

| Library | Version                  | Purpose                 | When to Use                                                                             |
| ------- | ------------------------ | ----------------------- | --------------------------------------------------------------------------------------- |
| plotly  | 6.8.0 (declared `>=5.0`) | Interactive plots       | NOT needed in Phase 38 (no results rendered); used from Phase 39 `[VERIFIED: pip list]` |
| pandas  | 3.0.3 (declared `>=2.0`) | CSV export / dataframes | NOT needed in Phase 38; used from Phase 39 `[VERIFIED: pip list]`                       |

### Alternatives Considered

| Instead of                       | Could Use                                            | Tradeoff                                                                                                                                                                                                                                                                                                                                                                                           |
| -------------------------------- | ---------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `st.navigation` + `st.Page`      | Legacy `pages/` directory convention                 | `pages/` auto-discovers files but forces the sidebar/config form to be duplicated (or imported) in every page file, and gives less control over ordering/grouping. `st.navigation` runs the entry-script sidebar on every page, which is exactly what UI-02/UI-07 need. `[CITED: st.navigation docs — "As soon as any session executes st.navigation, your app will ignore the pages/ directory"]` |
| Single non-widget key for config | Read `DeviceConfig` fields directly from widget keys | Widget-keyed state is GC'd when the widget isn't rendered (e.g. after navigating to a page that doesn't re-render that exact widget, or when a conditional field is hidden). A single non-widget key survives. `[CITED: widget-behavior docs]`                                                                                                                                                     |

**Installation:** No new packages. All dependencies already declared in `pyproject.toml` and installed. No `npm`/`pip install` step required for this phase.

**Version verification:**

- `streamlit 1.58.0` — verified via `pip list` this session. `st.navigation`/`st.Page` present in this version's docs. `[VERIFIED: pip list + docs.streamlit.io v1.58.0]`
- `plotly 6.8.0`, `pandas 3.0.3` — verified installed. `[VERIFIED: pip list]`

## Package Legitimacy Audit

> This phase installs **no** external packages. All required libraries (streamlit, plotly, pandas, etna) are already declared in `pyproject.toml` and installed in the environment. No slopcheck / registry verification needed — nothing is being added.

| Package                         | Registry | Age | Downloads | Source Repo | slopcheck | Disposition |
| ------------------------------- | -------- | --- | --------- | ----------- | --------- | ----------- |
| (none — no installs this phase) | —        | —   | —         | —           | —         | —           |

**Packages removed due to slopcheck [SLOP] verdict:** none
**Packages flagged as suspicious [SUS]:** none

## Phase Requirements

| ID    | Description                                                                                                    | Research Support                                                                                                                                                                                                                  |
| ----- | -------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| UI-01 | Launch with `streamlit run app/main.py`; multi-page navigation lists all simulation workflow pages             | `st.navigation` + `st.Page` API (System Architecture + Pattern 1). Placeholder pages registered for all downstream workflows (C-V, CCE, Field Map, Radiation Damage, Dark Current, Microdosimetry, Batch Sweep, Geometry Viewer). |
| UI-02 | Configure all `DeviceConfig` parameters via sidebar form controls on any page, persisted in `st.session_state` | Entry-script sidebar (renders on every page) + full field mapping table (all 11 DeviceConfig fields) + Pattern 2 (reactive mode selectors) + Pattern 3 (single non-widget key).                                                   |
| UI-07 | Navigate between all pages without losing device configuration                                                 | session_state persistence across pages `[CITED: session-state docs]` + single non-widget storage key (Pitfall 1).                                                                                                                 |

## Architecture Patterns

### System Architecture Diagram

```
                    streamlit run app/main.py
                              │
                              ▼
              ┌──────────────────────────────────┐
              │      app/main.py (ENTRY SCRIPT)   │
              │                                   │
              │  1. render_device_sidebar()  ─────┼──► SIDEBAR (runs every page)
              │       ├─ mode selectors (REACTIVE, no st.form)
              │       │    ├─ dimensionality: 1D / 2D  ──► gates half_width_um
              │       │    └─ doping_profile: graded/uniform ──► gates N_D vs N_D_junction/bulk/L_transition
              │       ├─ geometry / doping / operating inputs
              │       └─ assemble DeviceConfig(...)
              │                    │
              │                    ▼
              │     st.session_state["device_config"] = DeviceConfig(...)   ◄── single non-widget key
              │                                   │                             (persists across nav)
              │  2. pg = st.navigation({...pages...})                          │
              │  3. pg.run()  ────────────────────┼──► selected PAGE ──────────┘
              └──────────────────────────────────┘         │
                                                            ▼
                                          ┌─────────────────────────────────┐
                                          │  placeholder page (e.g. C-V)    │
                                          │  if "device_config" not in      │
                                          │     st.session_state:           │
                                          │       st.info("Configure a      │
                                          │         device in the sidebar") │
                                          │  else:                          │
                                          │       show config summary +     │
                                          │       "coming in Phase 39" note │
                                          └─────────────────────────────────┘
```

Data flow: user edits sidebar → Streamlit reruns entry script → sidebar re-renders and re-assembles `DeviceConfig` → stored under `device_config` key → `pg.run()` executes the selected page → page reads `device_config` from session_state. Navigation reruns the entry script (sidebar always renders), so config survives.

### Component Responsibilities

| File                                           | Responsibility                                                                                          |
| ---------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| `app/main.py`                                  | Entry script: `st.navigation` registration, `pg.run()`, calls the sidebar renderer                      |
| `app/components/device_sidebar.py` (suggested) | Renders all DeviceConfig controls, assembles `DeviceConfig`, writes `st.session_state["device_config"]` |
| `app/pages/*.py` (or page functions)           | One placeholder per workflow; empty-state guard + config summary                                        |
| `etna.DeviceConfig`                        | The dataclass contract the sidebar populates (do NOT redefine fields in the app)                        |

### Recommended Project Structure

```
app/
├── main.py                       # st.navigation entry script (streamlit run app/main.py)
├── components/
│   └── device_sidebar.py         # render_device_sidebar() -> writes st.session_state["device_config"]
└── pages/                        # placeholder page modules (or plain functions passed to st.Page)
    ├── cv.py                     # C-V (Phase 39)
    ├── cce.py                    # CCE (Phase 39)
    ├── field_map.py              # Field map / geometry (Phase 39/40)
    ├── radiation_damage.py       # (Phase 41)
    ├── dark_current.py           # (Phase 41)
    ├── microdosimetry.py         # (Phase 42)
    └── batch_sweep.py            # (Phase 42)
```

Note: with `st.navigation` these page modules do **not** need to live in a magic `pages/` directory — you can put them anywhere and reference them via `st.Page("app/pages/cv.py", ...)` or pass a callable. Keep the registration list in `main.py` **extensible** so Phases 39-42 append pages.

### DeviceConfig → Form Control Mapping (ALL fields — read from source)

Source of truth: `etna/api/device.py` (verified this session). All 11 fields:

| Field                    | Type / Default         | Gated by                                        | Suggested control                                     |
| ------------------------ | ---------------------- | ----------------------------------------------- | ----------------------------------------------------- |
| `epi_thickness_um`       | float = 10.0           | always                                          | `st.number_input` (µm)                                |
| `substrate_thickness_um` | float = 1.0            | always                                          | `st.number_input` (µm)                                |
| `half_width_um`          | Optional[float] = None | **dimensionality selector** (None=1D, float=2D) | shown only when "2D" selected; `st.number_input` (µm) |
| `N_A`                    | float = 1e19           | always                                          | doping input (see Pitfall 3 — wide range)             |
| `doping_profile`         | str = "graded"         | **mode selector** ("graded"/"uniform")          | `st.selectbox`/`st.radio` (REACTIVE, outside form)    |
| `N_D`                    | Optional[float] = None | `doping_profile == "uniform"`                   | shown only for "uniform"; doping input                |
| `N_D_junction`           | float = 2.93e15        | `doping_profile == "graded"`                    | shown only for "graded"; doping input                 |
| `N_D_bulk`               | float = 8.82e13        | `doping_profile == "graded"`                    | shown only for "graded"; doping input                 |
| `L_transition_um`        | float = 0.987          | `doping_profile == "graded"`                    | shown only for "graded"; `st.number_input` (µm)       |
| `T`                      | float = 300.0          | always                                          | `st.number_input` (K)                                 |
| `area_cm2`               | float = 1e-4           | always                                          | doping-style input (small value; see Pitfall 3)       |

**Two conditional dimensions to model in the form:**

1. **Dimensionality** (1D vs 2D) — a UI-only concept mapping to `half_width_um is None` (1D) vs a float (2D). There is no explicit "dimensionality" field on `DeviceConfig`; the app derives it. When "1D" is selected, set `half_width_um=None`; when "2D", collect the µm value.
2. **Doping profile** ("graded" vs "uniform") — an actual field. Graded uses `N_D_junction`, `N_D_bulk`, `L_transition_um` and leaves `N_D=None`. Uniform uses `N_D` and leaves the graded triplet at defaults (or hides them).

### Pattern 1: `st.navigation` multi-page shell

**What:** Register pages and run the selected one from a single entry script.
**When to use:** Every multi-page Streamlit app in current best practice; required here to satisfy UI-01/UI-02/UI-07 cleanly.
**Example:**

```python
# Source: docs.streamlit.io/develop/api-reference/navigation/st.navigation (v1.58.0)
import streamlit as st
from app.components.device_sidebar import render_device_sidebar

render_device_sidebar()  # runs on EVERY page; writes st.session_state["device_config"]

pages = [
    st.Page("app/pages/cv.py", title="C-V", icon="📉"),
    st.Page("app/pages/cce.py", title="CCE"),
    st.Page("app/pages/field_map.py", title="Field Map"),
    # ... append Phase 41/42 pages here later ...
]
pg = st.navigation(pages)   # sidebar nav by default
pg.run()
```

Note: `st.set_page_config()` (if used) must be called before other Streamlit commands in the entry script.

### Pattern 2: Reactive mode selectors (NOT inside `st.form`)

**What:** The dimensionality and doping-profile selectors trigger a rerun immediately so conditional fields appear/disappear.
**When to use:** Any time a control gates the visibility of other controls.
**Why:** Inside `st.form`, widget changes do not rerun until submit, so conditional fields cannot toggle before submit. `[CITED: st.form docs — "widgets inside a form cannot trigger conditional rendering before submit"]`
**Example:**

```python
# Mode selectors are plain widgets (reactive) — NOT wrapped in st.form
dim = st.sidebar.radio("Dimensionality", ["1D", "2D"])
profile = st.sidebar.selectbox("Doping profile", ["graded", "uniform"])

half_width_um = st.sidebar.number_input("Half-width (µm)", value=50.0) if dim == "2D" else None
if profile == "graded":
    N_D_junction = st.sidebar.number_input("N_D junction (cm⁻³)", value=2.93e15, format="%.3e")
    N_D_bulk     = st.sidebar.number_input("N_D bulk (cm⁻³)",     value=8.82e13, format="%.3e")
    L_transition_um = st.sidebar.number_input("Transition length (µm)", value=0.987)
    N_D = None
else:  # uniform
    N_D = st.sidebar.number_input("N_D (cm⁻³)", value=1e15, format="%.3e")
    N_D_junction, N_D_bulk, L_transition_um = 2.93e15, 8.82e13, 0.987  # defaults
```

### Pattern 3: Single non-widget key for the assembled config

**What:** Assemble a `DeviceConfig` from the widget values and store it under one dedicated key.
**Why:** This key is the **downstream contract** for Phases 39-43. Widget-keyed state is unreliable across navigation (GC when not rendered); a single non-widget key persists across pages within the session. `[CITED: widget-behavior docs + session-state docs]`
**Example:**

```python
from etna import DeviceConfig

st.session_state["device_config"] = DeviceConfig(
    epi_thickness_um=epi_thickness_um,
    substrate_thickness_um=substrate_thickness_um,
    half_width_um=half_width_um,
    N_A=N_A,
    doping_profile=profile,
    N_D=N_D,
    N_D_junction=N_D_junction,
    N_D_bulk=N_D_bulk,
    L_transition_um=L_transition_um,
    T=T,
    area_cm2=area_cm2,
)
```

### Pattern 4: Empty-state guard on every page (success criterion 4)

```python
import streamlit as st
cfg = st.session_state.get("device_config")
if cfg is None:
    st.info("Configure a device in the sidebar to begin.")
    st.stop()   # or simply return / skip the rest
# else: render config summary; "Simulation coming in Phase 39" placeholder
```

Because the entry-script sidebar always runs before `pg.run()`, `device_config` will normally be set on the very first render — but the guard is required by criterion 4 and defends against edge cases (e.g. a page loaded before the sidebar writes, or future refactors).

### Anti-Patterns to Avoid

- **Mode selector inside `st.form`:** conditional fields won't toggle until submit. Keep selectors reactive. `[CITED: st.form docs]`
- **Reading config from individual widget keys on result pages:** those keys are GC'd when not rendered on that run. Read the single `device_config` object instead. `[CITED: widget-behavior docs]`
- **Adding `st.cache_resource`/`st.cache_data` in Phase 38:** no expensive object creation happens here; module imports are already cached by `sys.modules`. Premature and unnecessary (see State of the Art).
- **Using the legacy `pages/` magic directory AND `st.navigation`:** once `st.navigation` runs, the `pages/` directory is ignored — mixing them causes confusion. Pick `st.navigation`. `[CITED: st.navigation docs]`
- **Redefining DeviceConfig fields in the app:** import `DeviceConfig` from `etna`; the dataclass is the single source of truth. If a field is added to the dataclass later, the form should be updated in lockstep (see Open Questions).

## Don't Hand-Roll

| Problem                              | Don't Build                                          | Use Instead                               | Why                                                                                                                 |
| ------------------------------------ | ---------------------------------------------------- | ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| Multi-page routing                   | Custom `st.radio` page switcher + `if/elif` dispatch | `st.navigation` + `st.Page`               | Native, gives real URLs, sidebar nav, and runs the entry-script sidebar on every page `[CITED: st.navigation docs]` |
| Cross-page state                     | Query params / cookies / files                       | `st.session_state["device_config"]`       | session_state persists across pages within a session by design `[CITED: session-state docs]`                        |
| Config schema                        | Hand-written dict of fields in the app               | `etna.DeviceConfig` dataclass         | Already the API contract; avoids drift between UI and library                                                       |
| Config cloning for sweeps (Phase 42) | `setattr`/`eval`                                     | `dataclasses.replace` / `ParametricSweep` | Already implemented; injection-safe `[VERIFIED: etna/api/sweep.py]`                                             |

**Key insight:** Streamlit's native primitives (`st.navigation`, `st.session_state`) already solve routing and state persistence. The only real design work in Phase 38 is faithfully mapping the 11 `DeviceConfig` fields (with two conditional dimensions) to reactive controls and choosing the single storage key — the rest is framework-provided.

## Common Pitfalls

### Pitfall 1: Config lost on navigation (widget-key GC)

**What goes wrong:** Result pages read `N_D`, `epi_thickness_um`, etc. from widget keys; after navigating away and back, some values are gone or reset to defaults.
**Why it happens:** "If a widget command for a specific widget instance isn't called during a script run, then none of its parts are retained, including its value in `st.session_state`." `[CITED: widget-behavior docs]`
**How to avoid:** Store the assembled `DeviceConfig` under a single non-widget key (`device_config`). Pages read that object, never the raw widget keys.
**Warning signs:** Config appears correct on the sidebar page but stale/default on a result page after navigation.

### Pitfall 2: Conditional fields don't appear (form batching)

**What goes wrong:** User switches doping profile to "uniform" but the `N_D` field never shows (or the graded fields never hide) until they click submit.
**Why it happens:** `st.form` batches all inputs and defers rerun to submit; widgets inside a form can't drive conditional rendering before submit. `[CITED: st.form docs]`
**How to avoid:** Keep the dimensionality and doping-profile selectors as plain reactive widgets (outside any `st.form`). Since Phase 38 runs no simulation, there is no batching benefit to using `st.form` at all — omitting it is the simplest correct choice.
**Warning signs:** Fields only update after a submit click; "uniform"/"graded" toggle feels one step behind.

### Pitfall 3: `st.number_input` across a 1e13–1e19 doping range

**What goes wrong:** Doping fields (`N_A`=1e19, `N_D_junction`=2.93e15, `N_D_bulk`=8.82e13, `N_D` uniform ~1e15, `area_cm2`=1e-4) span many orders of magnitude; default `st.number_input` step/format is awkward and can round or clamp.
**Why it happens:** `number_input` defaults assume integer-ish steps; scientific-notation values need explicit formatting.
**How to avoid:** Use `st.number_input(..., format="%.3e")` (scientific notation) for doping/area fields, or accept a text input parsed to float. Pick sensible `min_value`/`step`. This is a UI ergonomics decision, not a physics constraint. `[ASSUMED]` — control choice is at Claude's discretion; verify UX with the user if a specific input style is desired.
**Warning signs:** Users can't enter 8.82e13 cleanly; values snap to round numbers.

### Pitfall 4: devsim banner noise on startup

**What goes wrong:** `import etna` eagerly imports devsim, which prints a BLAS/LAPACK/UMFPACK banner to stderr and takes ~1.4s import time (measured this session). Streamlit will show this once at server start.
**Why it happens:** `etna/__init__.py` imports the facades, which import devsim at module load. Verified: `import etna` → `'devsim' in sys.modules` is `True`. `[VERIFIED: python -c timing this session]`
**How to avoid:** This is paid **once per server process** (Python caches modules in `sys.modules`), not per rerun — so it is not a performance problem for Phase 38. Do NOT wrap the import in `st.cache_resource` (that's for objects, not modules). If the banner is cosmetically undesirable, importing `etna` (or just `DeviceConfig`) lazily inside the sidebar function is optional but unnecessary. Since Phase 38 only needs `DeviceConfig`, `from etna import DeviceConfig` is sufficient and still triggers the eager devsim import — acceptable.
**Warning signs:** Slow first page load (~1-2s), devsim text in the server console — both benign.

## Code Examples

### Full sidebar assembly (reactive selectors + single-key storage)

```python
# app/components/device_sidebar.py
# Source: composed from docs.streamlit.io st.navigation / session-state / widget-behavior
import streamlit as st
from etna import DeviceConfig

def render_device_sidebar() -> None:
    st.sidebar.header("Device configuration")

    # --- reactive mode selectors (NOT inside st.form) ---
    dim = st.sidebar.radio("Dimensionality", ["1D", "2D"], key="ui_dim")
    profile = st.sidebar.selectbox("Doping profile", ["graded", "uniform"], key="ui_profile")

    # --- geometry ---
    epi = st.sidebar.number_input("Epi thickness (µm)", value=10.0, min_value=0.1)
    sub = st.sidebar.number_input("Substrate thickness (µm)", value=1.0, min_value=0.0)
    half_width = st.sidebar.number_input("Half-width (µm)", value=50.0, min_value=0.1) if dim == "2D" else None

    # --- doping ---
    N_A = st.sidebar.number_input("N_A substrate (cm⁻³)", value=1e19, format="%.3e")
    if profile == "graded":
        N_D_junction = st.sidebar.number_input("N_D junction (cm⁻³)", value=2.93e15, format="%.3e")
        N_D_bulk     = st.sidebar.number_input("N_D bulk (cm⁻³)",     value=8.82e13, format="%.3e")
        L_trans      = st.sidebar.number_input("Transition length (µm)", value=0.987, min_value=0.0)
        N_D = None
    else:
        N_D = st.sidebar.number_input("N_D uniform (cm⁻³)", value=1e15, format="%.3e")
        N_D_junction, N_D_bulk, L_trans = 2.93e15, 8.82e13, 0.987

    # --- operating conditions ---
    T = st.sidebar.number_input("Temperature (K)", value=300.0, min_value=1.0)
    area = st.sidebar.number_input("Area (cm²)", value=1e-4, format="%.3e")

    # --- assemble + store under the single downstream contract key ---
    st.session_state["device_config"] = DeviceConfig(
        epi_thickness_um=epi, substrate_thickness_um=sub, half_width_um=half_width,
        N_A=N_A, doping_profile=profile, N_D=N_D,
        N_D_junction=N_D_junction, N_D_bulk=N_D_bulk, L_transition_um=L_trans,
        T=T, area_cm2=area,
    )
```

### Placeholder page with empty-state guard

```python
# app/pages/cv.py
import streamlit as st

st.title("C-V Simulation")
cfg = st.session_state.get("device_config")
if cfg is None:
    st.info("Configure a device in the sidebar to begin.")
    st.stop()

st.write("Current device configuration:")
st.json({k: getattr(cfg, k) for k in cfg.__dataclass_fields__})
st.caption("Running C-V simulations is implemented in Phase 39.")
```

## State of the Art

| Old Approach                                 | Current Approach                               | When Changed                                                   | Impact                                                                                                         |
| -------------------------------------------- | ---------------------------------------------- | -------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `pages/` magic directory for multi-page apps | `st.navigation` + `st.Page` in an entry script | GA in Streamlit 1.36 (mid-2024); available in installed 1.58.0 | Entry-script sidebar runs on every page → clean cross-page config (UI-02/UI-07). `[CITED: st.navigation docs]` |
| Manual `st.experimental_*` state hacks       | Stable `st.session_state`                      | long stable                                                    | Store `DeviceConfig` under a non-widget key                                                                    |

**Deprecated/outdated:**

- Manually building page switchers with `st.radio` + `if/elif`: superseded by `st.navigation`.
- `st.cache` (unqualified): replaced by `st.cache_data` / `st.cache_resource` — but **neither is needed in Phase 38** (no expensive objects; module imports are `sys.modules`-cached). Introduce `st.cache_resource`/`st.cache_data` in Phase 39+ when facades actually build devsim devices and produce results.

## Assumptions Log

| #   | Claim                                                                                                                                                     | Section                       | Risk if Wrong                                                                                                                         |
| --- | --------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| A1  | `st.number_input(..., format="%.3e")` (or text-parse) is the right control for wide-range doping/area fields                                              | Pitfall 3 / Code Examples     | Low — UX ergonomics only; any numeric-entry control satisfies UI-02. Confirm with user if a specific input style is desired.          |
| A2  | Suggested `app/` directory layout (`components/`, `pages/`)                                                                                               | Recommended Project Structure | Low — organizational only; `st.navigation` doesn't mandate a directory. Planner may choose flat layout.                               |
| A3  | Placeholder page set = the eight downstream workflows (C-V, CCE, Field Map, Radiation Damage, Dark Current, Microdosimetry, Batch Sweep, Geometry Viewer) | System Architecture           | Low — UI-01 says "all simulation workflow pages"; derived from ROADMAP Phases 39-42. Confirm exact page list/titles in discuss-phase. |

## Open Questions

1. **Exact page list and titles/icons for the navigation.**
   - What we know: Phases 39-42 add C-V, CCE, Field Map, Geometry Viewer, Radiation Damage, Dark Current, Microdosimetry, Batch Sweep. UI-01 requires "all simulation workflow pages" listed.
   - What's unclear: Whether all 8 appear as placeholders now, or only the Phase 39 trio (C-V/CCE/Field) plus a stub for the rest. Whether pages are grouped (`st.navigation` supports section grouping).
   - Recommendation: Register **all** downstream workflow pages as placeholders now (satisfies UI-01 literally and makes the nav extensible), each with the empty-state guard. Keep the registration list in `main.py` easy to append to.

2. **Should "dimensionality" (1D/2D) be a labeled selector or inferred from a half-width toggle?**
   - What we know: `DeviceConfig` has no explicit dimensionality field; it's `half_width_um is None` vs float.
   - What's unclear: UX preference — explicit radio vs a "2D device?" checkbox.
   - Recommendation: Explicit `st.radio(["1D","2D"])` (clearest), mapping "1D"→`half_width_um=None`. Discretionary; confirm in discuss-phase if desired.

3. **DeviceConfig field drift.** If a future phase adds a field to `DeviceConfig`, the sidebar must be updated in lockstep. Not a blocker for Phase 38, but worth a note/test.
   - Recommendation: Consider a lightweight test asserting the sidebar assembles a `DeviceConfig` with all `__dataclass_fields__` populated (guards against silent field drift).

## Environment Availability

| Dependency              | Required By                         | Available | Version        | Fallback |
| ----------------------- | ----------------------------------- | --------- | -------------- | -------- |
| streamlit               | UI-01/02/07                         | ✓         | 1.58.0         | —        |
| etna (DeviceConfig) | UI-02                               | ✓         | local editable | —        |
| plotly                  | (Phase 39+) not this phase          | ✓         | 6.8.0          | —        |
| pandas                  | (Phase 39+) not this phase          | ✓         | 3.0.3          | —        |
| devsim                  | transitively imported by `etna` | ✓         | 2.10.0         | —        |

**Missing dependencies with no fallback:** none.
**Missing dependencies with fallback:** none.

## Validation Architecture

> `workflow.research: true` and `nyquist_validation` not disabled → included. Note: Phase 38 UI behavior is largely runtime/visual; automated tests cover the config-assembly and empty-state logic, while browser launch is manual/smoke.

### Test Framework

| Property           | Value                                                                                                            |
| ------------------ | ---------------------------------------------------------------------------------------------------------------- |
| Framework          | pytest (project standard; used by all 25 existing test modules) `[VERIFIED: pyproject.toml dev extras + tests/]` |
| Config file        | pyproject.toml (project uses `pytest -q`)                                                                        |
| Quick run command  | `pytest tests/test_app_device_sidebar.py -x` (new — Wave 0)                                                      |
| Full suite command | `pytest -q`                                                                                                      |

### Phase Requirements → Test Map

| Req ID | Behavior                                                                                                 | Test Type    | Automated Command                                                             | File Exists? |
| ------ | -------------------------------------------------------------------------------------------------------- | ------------ | ----------------------------------------------------------------------------- | ------------ |
| UI-02  | Sidebar assembles a `DeviceConfig` with ALL 11 fields populated from form values                         | unit         | `pytest tests/test_app_device_sidebar.py::test_assemble_config_all_fields -x` | ❌ Wave 0    |
| UI-02  | Graded profile → N_D=None + graded triplet set; uniform → N_D set                                        | unit         | `pytest tests/test_app_device_sidebar.py::test_doping_mode_mapping -x`        | ❌ Wave 0    |
| UI-02  | 1D → half_width_um=None; 2D → float                                                                      | unit         | `pytest tests/test_app_device_sidebar.py::test_dimensionality_mapping -x`     | ❌ Wave 0    |
| UI-07  | Config stored under single `device_config` key survives (logic-level: object round-trips through a dict) | unit         | `pytest tests/test_app_session.py::test_config_persistence_key -x`            | ❌ Wave 0    |
| UI-01  | Empty-state guard returns a prompt when `device_config` absent                                           | unit         | `pytest tests/test_app_pages.py::test_empty_state_guard -x`                   | ❌ Wave 0    |
| UI-01  | `streamlit run app/main.py` launches without error, nav lists pages                                      | manual/smoke | manual browser check (or `streamlit run` headless smoke)                      | manual       |

Testing note: Streamlit UIs are testable headlessly via `streamlit.testing.v1.AppTest` (the official app-testing harness). The planner should prefer `AppTest` for exercising sidebar → session_state → page-guard flows without a browser; keep pure config-assembly logic in a plain function (e.g. `assemble_config(values) -> DeviceConfig`) so it is unit-testable independent of Streamlit. `[CITED: docs.streamlit.io app testing]`

### Sampling Rate

- **Per task commit:** `pytest tests/test_app_*.py -x`
- **Per wave merge:** `pytest -q`
- **Phase gate:** Full suite green + manual `streamlit run app/main.py` smoke before `/gsd:verify-work`

### Wave 0 Gaps

- [ ] `tests/test_app_device_sidebar.py` — covers UI-02 (config assembly, doping/dimensionality mapping); requires refactoring config assembly into a pure `assemble_config()` function
- [ ] `tests/test_app_pages.py` — covers UI-01 empty-state guard
- [ ] `tests/test_app_session.py` — covers UI-07 persistence-key contract
- [ ] Consider `streamlit.testing.v1.AppTest`-based smoke test for nav + sidebar-on-every-page

## Security Domain

> `security_enforcement` not set to false in config → included. Phase 38 has a minimal attack surface: no file uploads, no simulation execution, no network, no secrets, no auth. It is a local single-user Streamlit tool assembling a dataclass from numeric form inputs.

### Applicable ASVS Categories

| ASVS Category         | Applies     | Standard Control                                                                                             |
| --------------------- | ----------- | ------------------------------------------------------------------------------------------------------------ |
| V2 Authentication     | no          | Local single-user tool; no auth in scope                                                                     |
| V3 Session Management | no          | Streamlit per-tab session; no credentials                                                                    |
| V4 Access Control     | no          | No multi-user access model                                                                                   |
| V5 Input Validation   | yes (light) | Numeric inputs constrained via `st.number_input` min/max; `DeviceConfig` field types. No eval/exec on input. |
| V6 Cryptography       | no          | No secrets or crypto in this phase                                                                           |

### Known Threat Patterns for Streamlit config UI

| Pattern                                                             | STRIDE          | Standard Mitigation                                                                                                                                                                                                                  |
| ------------------------------------------------------------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Arbitrary attribute injection into config (relevant Phase 42 sweep) | Tampering       | Already mitigated in library: `ParametricSweep` uses `dataclasses.replace`, never `setattr`/`eval` `[VERIFIED: etna/api/sweep.py]`. Phase 38 assembles `DeviceConfig` with explicit named kwargs — no dynamic attribute setting. |
| Malformed numeric input (NaN/inf/negative thickness)                | Tampering / DoS | `st.number_input` `min_value` guards; optional validation before storing config. No simulation runs in Phase 38 so blast radius is nil.                                                                                              |
| (No file upload in this phase)                                      | —               | File-input threat surface (`run_microdosimetry` CSV) is Phase 42, not here.                                                                                                                                                          |

## Sources

### Primary (HIGH confidence)

- `docs.streamlit.io/develop/api-reference/navigation/st.navigation` — st.navigation/st.Page, entry-script sidebar behavior, pages/ dir superseded
- `docs.streamlit.io/develop/concepts/architecture/session-state` — persistence across reruns and pages, no persistence on server crash
- `docs.streamlit.io/develop/concepts/architecture/widget-behavior` — keyed widget state GC when not rendered; placeholder-key pattern
- `docs.streamlit.io/develop/api-reference/execution-flow/st.form` — form batching, no conditional rendering before submit, callback constraints
- Source code (verified this session): `etna/api/device.py` (DeviceConfig 11 fields), `etna/api/results.py` (SimResult/MeshData), `etna/api/sweep.py` (ParametricSweep), `etna/__init__.py` (public API), `etna/api/simulation.py` (facade signatures)
- `pip list` this session: streamlit 1.58.0, plotly 6.8.0, pandas 3.0.3, devsim 2.10.0
- `python -c "import etna"` timing: eager devsim import confirmed (~1.4s, `'devsim' in sys.modules` True)

### Secondary (MEDIUM confidence)

- `docs.streamlit.io` app testing (`streamlit.testing.v1.AppTest`) — referenced for headless UI tests; API details to be confirmed by planner if AppTest tests are written

### Tertiary (LOW confidence)

- Exact GA version of st.navigation (stated 1.36 from general knowledge) — not load-bearing; availability confirmed empirically in installed 1.58.0 `[ASSUMED for the specific version number]`

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH — all versions verified via `pip list`; no new packages
- Architecture (st.navigation, session_state, form): HIGH — verified against official docs this session
- DeviceConfig field mapping: HIGH — read directly from source `etna/api/device.py`
- Pitfalls: HIGH — widget-GC and form-batching confirmed verbatim in docs; import cost measured empirically
- Test architecture: MEDIUM — AppTest specifics to be confirmed by planner

**Research date:** 2026-07-10
**Valid until:** 2026-08-09 (Streamlit is moderately fast-moving; st.navigation is stable but verify if Streamlit is upgraded past 1.58)
