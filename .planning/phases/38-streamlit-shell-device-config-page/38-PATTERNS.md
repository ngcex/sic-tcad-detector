# Phase 38: Streamlit Shell + Device Config Page - Pattern Map

**Mapped:** 2026-07-10
**Files analyzed:** 14 (11 app/config files + 3 test files)
**Analogs found:** 3 / 14 (all three test files map to `tests/test_api_sweep.py`; the 11 app/config files are greenfield — no Streamlit code exists anywhere in the repo)

> **Headline finding:** This is the FIRST UI phase. There is **no `app/` directory, no `.streamlit/` directory, and zero `import streamlit` anywhere in the codebase** (verified via `grep -rl "import streamlit"` → empty). Therefore the Streamlit UI files (`main.py`, `device_sidebar.py`, page modules, `config.toml`) have **no in-repo analog** — the planner must copy patterns from RESEARCH.md's Code Examples, which are already concrete and line-level. Do **not** force-fit a Python-library analog onto a Streamlit UI file.
>
> The real, strong in-repo analogs are: (1) the three **test files** → copy `tests/test_api_sweep.py` (pure-Python, no devsim, `DeviceConfig()` assembly, `pytest.raises(TypeError)`); and (2) the shared **`DeviceConfig` contract** + the **`build_device()` "map every field explicitly" convention** → the source-of-truth every UI file marshals into.

## File Classification

| New/Modified File                  | Role                 | Data Flow                           | Closest Analog                                            | Match Quality                      |
| ---------------------------------- | -------------------- | ----------------------------------- | --------------------------------------------------------- | ---------------------------------- |
| `.streamlit/config.toml`           | config               | —                                   | none (greenfield)                                         | no analog                          |
| `app/main.py`                      | entry-script / route | request-response (rerun)            | none (greenfield)                                         | no analog — use RESEARCH Pattern 1 |
| `app/components/device_sidebar.py` | component (form)     | transform (form values → dataclass) | `petringa/api/device.py` `build_device` (convention only) | partial — convention analog        |
| `app/pages/home.py`                | page                 | request-response                    | none (greenfield)                                         | no analog — use RESEARCH Pattern 4 |
| `app/pages/cv.py`                  | page                 | request-response                    | none (greenfield)                                         | no analog                          |
| `app/pages/cce.py`                 | page                 | request-response                    | none (greenfield)                                         | no analog                          |
| `app/pages/field_map.py`           | page                 | request-response                    | none (greenfield)                                         | no analog                          |
| `app/pages/radiation_damage.py`    | page                 | request-response                    | none (greenfield)                                         | no analog                          |
| `app/pages/dark_current.py`        | page                 | request-response                    | none (greenfield)                                         | no analog                          |
| `app/pages/microdosimetry.py`      | page                 | request-response                    | none (greenfield)                                         | no analog                          |
| `app/pages/batch_sweep.py`         | page                 | request-response                    | none (greenfield)                                         | no analog                          |
| `tests/test_app_device_sidebar.py` | test                 | —                                   | `tests/test_api_sweep.py`                                 | exact                              |
| `tests/test_app_session.py`        | test                 | —                                   | `tests/test_api_sweep.py`                                 | exact                              |
| `tests/test_app_pages.py`          | test                 | —                                   | `tests/test_api_sweep.py`                                 | exact                              |

**Note on page modules:** all 8 page modules (Home + 7 workflow pages) are structurally identical for Phase 38 — a title, the empty-state guard, and a config-summary/"coming in Phase N" caption. They share ONE pattern (see Shared Patterns → Empty-State Guard + Config Summary). The planner should treat them as instances of a single template, not 8 distinct patterns.

## Pattern Assignments

### `tests/test_app_device_sidebar.py`, `tests/test_app_session.py`, `tests/test_app_pages.py` (test)

**Analog:** `tests/test_api_sweep.py` — the project's canonical **fast, pure-Python, no-devsim** test module. This is the exact template the three Wave-0 app tests must follow: construct a `DeviceConfig()`, exercise pure logic, assert field values, use a **fake object** instead of the real solver, and assert `TypeError` on the injection guard.

**Imports pattern** (`tests/test_api_sweep.py` lines 15-19):

```python
import numpy as np
import pytest

from petringa import DeviceConfig, ParametricSweep
from petringa.api.results import SimResult
```

For the app tests, import the pure logic under test instead — e.g. `from app.components.device_sidebar import assemble_config` (RESEARCH recommends refactoring config assembly into a pure `assemble_config(values) -> DeviceConfig` so it is unit-testable without Streamlit; Validation Architecture note).

**Fake-object pattern — avoid the real solver** (`tests/test_api_sweep.py` lines 22-33):

```python
def fake_sim(cfg, **kw):
    """Canned sim_fn: echoes the swept field back through the SimResult.

    Never builds a devsim device — returns a SimResult whose x carries the
    injected `epi_thickness_um` so config injection can be asserted.
    """
    return SimResult(
        config=cfg,
        sim_type="fake",
        x=np.array([cfg.epi_thickness_um]),
        y=np.array([1.0]),
    )
```

App-test analog: exercise the pure `assemble_config()` (or the empty-state guard as a pure function) — never launch a live Streamlit server in a unit test. If a headless UI test is wanted, RESEARCH points at `streamlit.testing.v1.AppTest`, but keep the config-assembly logic in a plain function so the bulk of tests stay framework-free like this analog.

**Assembly + length assertion pattern** (`tests/test_api_sweep.py` lines 36-47) — maps to `test_assemble_config_all_fields`:

```python
def test_returns_list_of_correct_length():
    sweep = ParametricSweep(
        base_config=DeviceConfig(),
        param="epi_thickness_um",
        values=[5, 10, 20],
        sim_fn=fake_sim,
    )
    results = sweep.run()

    assert isinstance(results, list)
    assert len(results) == 3
    assert all(isinstance(r, SimResult) for r in results)
```

For `test_assemble_config_all_fields`: call `assemble_config(...)`, assert the result `isinstance(cfg, DeviceConfig)`, and assert **every** field in `cfg.__dataclass_fields__` is populated (guards against silent DeviceConfig field drift — RESEARCH Open Question 3).

**Field-value assertion pattern** (`tests/test_api_sweep.py` lines 50-63) — maps to `test_doping_mode_mapping` / `test_dimensionality_mapping`:

```python
def test_config_injection():
    values = [5, 10, 20]
    sweep = ParametricSweep(...)
    results = sweep.run()

    assert [r.config.epi_thickness_um for r in results] == values
```

For `test_doping_mode_mapping`: `assemble_config(profile="graded", ...)` → assert `cfg.N_D is None` and the graded triplet is set; `profile="uniform"` → assert `cfg.N_D` set. For `test_dimensionality_mapping`: `"1D"` → `cfg.half_width_um is None`; `"2D"` → `cfg.half_width_um == <float>`.

**Guard / raises pattern** (`tests/test_api_sweep.py` lines 66-73):

```python
def test_unknown_param_raises():
    with pytest.raises(TypeError):
        ParametricSweep(
            base_config=DeviceConfig(),
            param="not_a_field",
            values=[1],
            sim_fn=fake_sim,
        ).run()
```

Reuse the `pytest.raises(...)` structure for any validation guard the sidebar adds.

**CRITICAL marker constraint:** These app tests are **logic-only and must NOT carry `@pytest.mark.slow`.** The `slow` marker is reserved exclusively for devsim integration tests (see `tests/test_api_device.py` line 18 `@pytest.mark.slow` on `TestBuildDevice`, and `pytest.ini` line 3 which defines `slow` as ">10s devsim integration tests"). The app tests build no device — they belong in the fast default suite alongside `test_api_sweep.py` (which is unmarked). This is also what makes the RESEARCH command `pytest tests/test_app_*.py -x` fast.

---

### `app/components/device_sidebar.py` (component, transform)

**Analog:** None for the Streamlit rendering (greenfield). **Convention analog** for the field-marshalling half: `petringa/api/device.py` → `build_device()` (lines 46-115).

**Why this is the convention analog:** RESEARCH recommends refactoring config assembly into a pure `assemble_config(values) -> DeviceConfig`. That function should mirror `build_device`'s core convention: **map every DeviceConfig field explicitly with named kwargs; never rely on defaults; never use `setattr`/`eval`.** `build_device` states this verbatim ("Every DeviceConfig field is mapped explicitly — core constructor defaults are never relied upon", lines 54-56).

**Dispatch-on-mode convention** (`petringa/api/device.py` lines 75-96) — the sidebar's dimensionality logic mirrors this exact 1D/2D dispatch on `half_width_um is None`:

```python
if config.half_width_um is None:
    # 1D branch ...
    return device_info
# 2D branch: half_width_um is a float ...
```

Sidebar analog: `"1D"` selected → set `half_width_um=None`; `"2D"` → collect the float. Same None-vs-float switch the library already keys on.

**Explicit named-kwarg assembly** (`petringa/api/device.py` lines 79-90 show the named-kwarg style; the sidebar's final assembly should look like RESEARCH Pattern 3 / Code Example):

```python
# From RESEARCH.md Code Examples (no in-repo Streamlit analog exists):
st.session_state["device_config"] = DeviceConfig(
    epi_thickness_um=epi, substrate_thickness_um=sub, half_width_um=half_width,
    N_A=N_A, doping_profile=profile, N_D=N_D,
    N_D_junction=N_D_junction, N_D_bulk=N_D_bulk, L_transition_um=L_trans,
    T=T, area_cm2=area,
)
```

**Streamlit rendering pattern (NO in-repo analog):** copy directly from RESEARCH.md "Full sidebar assembly" Code Example (lines 305-347) and the UI-SPEC "Sidebar Device-Config Contract" (Groups 1-4). Key non-negotiables from both docs:

- Mode selectors (Dimensionality radio, Doping-profile selectbox) render **first** and **outside any `st.form`** (reactive — RESEARCH Pattern 2 / Pitfall 2).
- Doping + area fields use `format="%.3e"` (RESEARCH Pitfall 3; UI-SPEC Group 3/4).
- Every label carries its unit in parentheses (UI-SPEC "Formatting rule").
- Groups separated by `st.sidebar.header` + `st.sidebar.divider()` (UI-SPEC Field grouping).
- Store under the single non-widget key `st.session_state["device_config"]` (RESEARCH Pattern 3).

---

### `app/main.py` (entry-script / route)

**Analog:** None (greenfield). Copy from RESEARCH.md **Pattern 1** (`st.navigation` multi-page shell, lines 175-190) and UI-SPEC "Navigation & Page Contract".

Non-negotiables from the docs:

- `st.set_page_config(...)` called **first**, before any other Streamlit call, with `page_title="petringa — SiC TCAD Simulator"`, `layout="wide"`, `initial_sidebar_state="expanded"` (UI-SPEC).
- `render_device_sidebar()` called **before** `pg.run()` so the sidebar renders on every page (satisfies UI-02/UI-07).
- Register Home + all 7 workflow pages via `st.Page(...)`; keep the list easy to append to (UI-SPEC page table; RESEARCH Open Question 1). Page titles/icons/order per UI-SPEC "Page list" table.
- Do **not** use the legacy `pages/` magic directory (ignored once `st.navigation` runs — RESEARCH anti-pattern).

---

### `app/pages/*.py` (page × 8, request-response)

**Analog:** None (greenfield). All 8 pages share ONE template — copy from RESEARCH.md "Placeholder page with empty-state guard" (lines 351-363) + UI-SPEC "Empty-state guard contract" + "Main-panel config summary". See Shared Patterns below (the guard + summary is cross-cutting across all 8).

## Shared Patterns

### Empty-State Guard (applies to ALL 8 page modules — success criterion 4)

**Source:** RESEARCH.md Pattern 4 (lines 243-250) + UI-SPEC "Empty-state guard contract" (lines 144-149). No in-repo analog (greenfield). This is the SINGLE most-copied snippet in the phase — every page begins with it, verbatim:

```python
cfg = st.session_state.get("device_config")
if cfg is None:
    st.info("Configure a device in the sidebar to begin.")
    st.stop()
```

Copy string exactly as specified in UI-SPEC Copywriting Contract ("Configure a device in the sidebar to begin.").

### Config-Summary Readout (applies to Home + all placeholder pages)

**Source:** UI-SPEC "Main-panel config summary" (line 197) + RESEARCH page example (line 362). Renders the persisted config as a monospace JSON readout — visually confirms UI-07 persistence across nav:

```python
st.json({k: getattr(cfg, k) for k in cfg.__dataclass_fields__})
```

Iterating `cfg.__dataclass_fields__` (rather than hardcoding field names) is the same drift-resistant idiom the test's all-fields assertion uses — keep them consistent.

### DeviceConfig Contract (the source of truth — referenced by sidebar, pages, and all tests)

**Source:** `petringa/api/device.py` lines 25-43. **NEVER redefine these fields in the app** (RESEARCH anti-pattern; UI-SPEC line 157). Import via `from petringa import DeviceConfig`. The full 11-field contract with defaults:

```python
@dataclass
class DeviceConfig:
    epi_thickness_um: float = 10.0          # always
    substrate_thickness_um: float = 1.0     # always
    half_width_um: Optional[float] = None   # gated: 1D=None / 2D=float
    N_A: float = 1e19                       # always
    doping_profile: str = "graded"          # mode selector: "graded"/"uniform"
    N_D: Optional[float] = None             # gated: uniform only
    N_D_junction: float = 2.93e15           # gated: graded only
    N_D_bulk: float = 8.82e13               # gated: graded only
    L_transition_um: float = 0.987          # gated: graded only
    T: float = 300.0                        # always
    area_cm2: float = 1e-4                  # always
```

These defaults are exactly the sidebar widget default values specified in UI-SPEC Groups 2-4. Use them as the widget `value=` defaults so an untouched sidebar produces `DeviceConfig()`.

### Test Convention (applies to all 3 test files)

**Source:** `tests/test_api_sweep.py` (whole file) + `pytest.ini`. Fast, no-devsim, pure-Python; construct real `DeviceConfig()`; fake objects instead of the solver; `pytest.raises` for guards; **no `@pytest.mark.slow`** (reserved for devsim integration — `tests/test_api_device.py` line 18). Run via `pytest tests/test_app_*.py -x`.

## No Analog Found

The Streamlit UI + config files are genuinely greenfield — no Streamlit, no `app/`, no `.streamlit/` exists in the repo. The planner should copy patterns from **RESEARCH.md Code Examples** (concrete, line-level) and **UI-SPEC** (the visual/interaction contract), NOT from any Python-library file:

| File                                   | Role         | Data Flow        | Reason / Where to copy from instead                                                                                                                                                                     |
| -------------------------------------- | ------------ | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `.streamlit/config.toml`               | config       | —                | No theming config exists yet. Copy `[theme]` keys from UI-SPEC "Color" table (backgroundColor `#FFFFFF`, secondaryBackgroundColor `#F0F2F6`, primaryColor `#1F6FEB`, textColor `#1A1A1A`).              |
| `app/main.py`                          | entry-script | request-response | No Streamlit entry script exists. Copy RESEARCH Pattern 1 + UI-SPEC nav/page-config contract.                                                                                                           |
| `app/components/device_sidebar.py`     | component    | transform        | Rendering half is greenfield (copy RESEARCH "Full sidebar assembly" + UI-SPEC Sidebar Contract). Only the field-marshalling convention has an in-repo analog (`build_device`, see Pattern Assignments). |
| `app/pages/home.py` + 7 workflow pages | page         | request-response | No Streamlit pages exist. Copy the single shared template (RESEARCH placeholder-page example + Shared Patterns guard/summary above).                                                                    |

## Metadata

**Analog search scope:** whole repo (`grep -rl "import streamlit"` → empty; `grep -rl "st.navigation|st.Page|session_state"` → empty), `tests/` (34 test modules scanned; `test_api_sweep.py` and `test_api_device.py` read in full), `petringa/api/` (`device.py`, `results.py`, `sweep.py`, `__init__.py` read in full).
**Files scanned in depth:** 6 (`petringa/api/device.py`, `results.py`, `sweep.py`, `petringa/__init__.py`, `tests/test_api_sweep.py`, `tests/test_api_device.py`) + `pyproject.toml`, `pytest.ini`.
**Pattern extraction date:** 2026-07-10
