# Phase 39: C-V, CCE, Field Map Pages + CSV Download - Pattern Map

**Mapped:** 2026-07-11
**Files analyzed:** 8 (3 modified pages, 1 new component, 4 new test files)
**Analogs found:** 7 / 8 (the mocked-`run_*` AppTest pattern has no in-repo analog — see No Analog Found)

> Scope note: The `run_cv` / `run_cce` / `run_field` / `SimResult` API contract is
> owned verbatim by `39-RESEARCH.md` (§1). This document does NOT restate it — it
> only maps each new/modified file to the closest existing code to copy structure
> from, with line-anchored excerpts.

## File Classification

| New/Modified File                     | Role           | Data Flow            | Closest Analog                                                                                                            | Match Quality                         |
| ------------------------------------- | -------------- | -------------------- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------- |
| `app/workflows/cv.py` (modify)        | component/page | request-response     | `app/workflows/home.py`                                                                                                   | exact (same render() page shape)      |
| `app/workflows/cce.py` (modify)       | component/page | request-response     | `app/workflows/home.py`                                                                                                   | exact                                 |
| `app/workflows/field_map.py` (modify) | component/page | request-response     | `app/workflows/home.py`                                                                                                   | exact                                 |
| `app/components/results.py` (new)     | utility        | transform / file-I/O | `app/components/device_sidebar.py` (pure-seam structure); `petringa/core/plotting.py` (labels/units spec, NOT importable) | role-match                            |
| `tests/test_app_cv_page.py` (new)     | test           | request-response     | `tests/test_app_pages.py::test_empty_state_guard` (from_function shape)                                                   | role-match (mock layer has no analog) |
| `tests/test_app_cce_page.py` (new)    | test           | request-response     | `tests/test_app_pages.py`                                                                                                 | role-match                            |
| `tests/test_app_field_page.py` (new)  | test           | request-response     | `tests/test_app_pages.py`                                                                                                 | role-match                            |
| `tests/test_app_csv_export.py` (new)  | test           | transform            | `tests/test_app_device_sidebar.py` (pure-function unit test, no Streamlit)                                                | exact                                 |

---

## Pattern Assignments

### `app/workflows/cv.py` / `cce.py` / `field_map.py` (component/page, request-response)

All three pages share one skeleton. **Analog:** `app/workflows/home.py` (full file, 8-25) — the canonical `render()` page: module docstring documenting the "all `st.*` inside render(), no module-level side effects → importable for `st.Page` + `AppTest.from_function`" contract, then title, then empty-state guard.

**Module docstring + no-side-effects contract** — `app/workflows/home.py:1-6`:

```python
"""Home landing page: orientation line + current device config summary.

All Streamlit calls live inside render() — no module-level side effects —
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests.
"""
```

Copy this docstring style per page; **remove** the current placeholder caption
`"Running this simulation is implemented in Phase 39."` (present in cv.py:17, cce.py:17, field_map.py:19).

**Empty-state guard (preserve verbatim string)** — `app/workflows/home.py:20-23`:

```python
cfg = st.session_state.get("device_config")
if cfg is None:
    st.info("Configure a device in the sidebar to begin.")
    st.stop()
```

The string `"Configure a device in the sidebar to begin."` is asserted by
`tests/test_app_pages.py:90` — do not alter it.

**Config-echo idiom (reuse for a "device parameters used" expander)** — `app/workflows/home.py:25`:

```python
st.json({k: getattr(cfg, k) for k in cfg.__dataclass_fields__})
```

**2D pre-check guard (NEW pattern — not in analog; spec from RESEARCH Pitfall 1).**
No existing page has this; the closest structural precedent is the guard-then-`st.stop()`
shape above. Insert immediately after the empty-state guard, inside the Run handler
per RESEARCH §Pitfall-1. Must be a pre-check, NOT try/except (`run_field` returns
empty arrays instead of raising):

```python
if cfg.half_width_um is not None:
    st.warning("These workflows are 1D-only. 2D field visualization arrives "
               "in Phase 40 (geometry viewer). Set Dimensionality to 1D in the sidebar.")
    st.stop()
```

**Run → cache → render → download pattern (NEW; spec from RESEARCH §5, load-bearing).**
No existing page runs a simulation; there is no analog. session_state key convention
is inherited from `app/components/device_sidebar.py:150` (`st.session_state["device_config"] = ...`);
namespace new keys distinctly as `cv_result` / `cce_result` / `field_result` (avoid the
`cfg_*` widget prefix). Reference the facade as `petringa.run_cv` (module attribute) so
tests can monkeypatch it — see RESEARCH §Validation A6.

```python
import petringa
if st.button("Run simulation"):
    st.session_state["cv_result"] = petringa.run_cv(cfg)

result = st.session_state.get("cv_result")
if result is not None:
    st.plotly_chart(build_cv_figure(result))
    st.plotly_chart(build_mott_schottky_figure(result))
    st.download_button("Download CSV", data=to_csv_bytes(result),
                       file_name="cv_result.csv", mime="text/csv")
```

`st.cache_data` is forbidden (DeviceConfig unhashable — RESEARCH §5); caching MUST be
manual `st.session_state`.

---

### `app/components/results.py` (utility, transform / file-I/O)

Shared Plotly `go.Figure` builders + `to_csv_bytes(result: SimResult) -> bytes`.

**Structural analog:** `app/components/device_sidebar.py` — the established "pure seam +
thin Streamlit-facing wrapper" module. `assemble_config` (device_sidebar.py:24-50) is a
pure, Streamlit-free function imported directly by unit tests. Mirror that: `to_csv_bytes`
and the figure builders must be **pure** (take a `SimResult`, no `st.*` calls) so
`tests/test_app_csv_export.py` can import and test them without a Streamlit runtime.

Import convention (top of module) — mirror `device_sidebar.py:17-21`:

```python
from __future__ import annotations
import streamlit as st          # only if any st.* helper lives here; builders should stay pure
from petringa import SimResult  # type hints; also DeviceConfig via result.config
```

For serialization use `dataclasses.asdict(result.config)` (RESEARCH §4, verified) and
`petringa.__version__` for the metadata header; use `pandas.DataFrame(...).to_csv(index=False)`
(project convention, RESEARCH §Don't-Hand-Roll) — no manual string joins, no temp file
(`df.to_csv(index=False).encode()`).

**Plotly builders — labels/units SPEC (do NOT import `plotting.py`; it is matplotlib-only).**
`petringa/core/plotting.py` is the authoritative source for axis labels, titles, unit
conversions, and reference lines. Mirror these into `go.Figure`:

C-V + Mott-Schottky — `petringa/core/plotting.py:367-395`:

```python
# C_vs_V:   x=voltages,             ylabel "Capacitance (F/cm²)",  title "C-V Characteristic"
# 1/C2_vs_V: y=one_over_C_squared,  ylabel "1/C² (cm⁴/F²)",        title "Mott-Schottky Plot (1/C² vs V)"
# W_vs_V (bonus): depletion_widths * 1e4 (cm->µm), ylabel "Depletion Width (µm)"
# shared xlabel "Voltage (V)"
```

CCE — `petringa/core/plotting.py:491-517`:

```python
V = np.abs(np.asarray(cce_data["voltages"]))   # plot |V| on x
cce = np.asarray(cce_data["cce_values"])
# axhline y=1.0 reference; ylim [0, 1.1]
# xlabel "|Reverse Bias| (V)", ylabel "Charge Collection Efficiency", title "CCE vs Reverse Bias"
```

Field / E-field vs depth — `petringa/core/plotting.py:59-68`:

```python
x_um = np.asarray(x_cm) * 1e4  # cm -> um   <-- DO NOT replicate: run_field.x is ALREADY µm (RESEARCH §Pitfall 3, line 154)
# xlabel "Depth (µm)", ylabel "Electric Field (V/cm)"
```

For the Field page use `result.x` (µm, node-aligned) with `result.y` (E-field, V/cm) and
`metadata["potential"]` (V); bonus `metadata["net_doping"]` (cm⁻³). Both E-field AND
potential come from ONE `run_field(cfg)` call (RESEARCH §6). Ignore `mesh` (reserved for
Phase 40; `mesh.x_coords` is cm, not µm — Pitfall 3).

---

### `tests/test_app_csv_export.py` (test, transform)

**Analog:** `tests/test_app_device_sidebar.py:1-40` — the pure-function unit-test precedent.
It imports a Streamlit-free seam directly and asserts on its output, with no `AppTest` and
no server. Copy this exact shape for `to_csv_bytes`: construct a `SimResult` fixture in
memory, call `to_csv_bytes(result)`, assert on the decoded bytes (column headers, the
commented `#` metadata header lines, values). This is where CSV-content assertions live —
AppTest 1.55 has no `download_button` accessor (RESEARCH §Validation), so content cannot be
asserted through the page test.

Import style to mirror — `tests/test_app_device_sidebar.py:9-10`:

```python
from petringa import DeviceConfig
from app.components.device_sidebar import assemble_config   # -> from app.components.results import to_csv_bytes
```

### `tests/test_app_cv_page.py` / `test_app_cce_page.py` / `test_app_field_page.py` (test, request-response)

**Partial analog:** `tests/test_app_pages.py::test_empty_state_guard` (79-90) — the
`AppTest.from_function` wrapper whose body imports `render` INSIDE itself (required:
`from_function` bodies must be self-contained), runs it, then asserts `at.exception == []`
and inspects `at.info`:

```python
def test_empty_state_guard():
    def _run_cv_page():
        from app.workflows.cv import render
        render()
    at = AppTest.from_function(_run_cv_page)
    at.run()
    assert at.exception == [], f"page crashed on empty session_state: {at.exception}"
    info_texts = [el.value for el in at.info]
    assert "Configure a device in the sidebar to begin." in info_texts
```

Available assertion surface (RESEARCH §Validation, verified): `at.exception == []`,
`at.button` (label), `at.caption` / `at.markdown` marker text, `at.warning` (2D guard),
and `at.session_state` result keys. There is **NO** `plotly_chart` or `download_button`
accessor — assert the result landed in `session_state`, not on the chart/download widgets.

To seed a config, set `at.session_state["device_config"] = DeviceConfig(...)` before
`at.run()` (2D: pass `half_width_um=50.0` to exercise the guard).

---

## Shared Patterns

### session_state result caching

**Source convention:** `app/components/device_sidebar.py:150` (`st.session_state["device_config"] = assemble_config(values)`)
**Apply to:** all three pages
Cache each `SimResult` under a distinct, non-`cfg_`-prefixed key: `cv_result`, `cce_result`,
`field_result`. Required (not an optimization): a Download-button click reruns the script and
the Run button returns `False`; without the cache the chart+download vanish and a naive rerun
re-solves devsim. `st.cache_data` is impossible (DeviceConfig unhashable).

### render() / no-side-effects page contract

**Source:** `app/workflows/home.py:1-6, 13`
**Apply to:** all three pages
Every `st.*` call inside `render()`; module docstring documents the `st.Page` + AppTest dual use.

### Empty-state guard

**Source:** `app/workflows/home.py:20-23`
**Apply to:** all three pages
Verbatim string `"Configure a device in the sidebar to begin."` (asserted at `tests/test_app_pages.py:90`).

### Config serialization

**Source:** `app/workflows/home.py:25` (`{k: getattr(cfg, k) for k in cfg.__dataclass_fields__}`) and `dataclasses.asdict(config)` (RESEARCH §4)
**Apply to:** results.py CSV metadata header + an optional per-page "device parameters used" expander

### Facade reference for mockability

**Source:** RESEARCH §Validation A6 (no in-repo analog)
**Apply to:** all three pages
Call `petringa.run_cv(...)` (i.e. `import petringa; petringa.run_cv`), NOT `from petringa import run_cv`, so tests can `monkeypatch.setattr(petringa, "run_cv", fake)`.

### CSV serialization mechanics

**Source:** project convention `df.to_csv(index=False)` (RESEARCH §4, e.g. tests/test_mc_coupling.py:71)
**Apply to:** results.py `to_csv_bytes`
Build bytes in-memory (`df.to_csv(index=False).encode()`); no temp file. Prepend commented `#` metadata header (`[ASSUMED]` format per RESEARCH §4).

---

## No Analog Found

Files/patterns with no close match in the codebase (planner should use RESEARCH.md instead of copying existing code):

| Item                                                                                  | Role           | Data Flow        | Reason                                                                                                                                                                                                                            |
| ------------------------------------------------------------------------------------- | -------------- | ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Plotly `go.Figure` builders in `results.py`                                           | utility        | transform        | No Plotly code exists anywhere (RESEARCH §3, grep-verified). Build from scratch; use `plotting.py` labels/units as spec only — it is matplotlib, not importable into `st.plotly_chart`.                                           |
| Run → cache → download flow in the three pages                                        | component/page | request-response | No existing page runs a simulation or caches a result; the session_state write convention is the only inherited piece. Spec = RESEARCH §5.                                                                                        |
| 2D pre-check guard                                                                    | component/page | request-response | No page guards on `half_width_um`. Spec = RESEARCH §Pitfall 1.                                                                                                                                                                    |
| monkeypatch + `AppTest.from_function` intercepting `petringa.run_*` in the page tests | test           | request-response | Nothing in the repo monkeypatches `petringa.run_cv`. `test_empty_state_guard` gives the `from_function` wrapper shape only. Spec = RESEARCH §Validation (A6 + Wave-0 spike verifying interception before writing page structure). |
| Commented-`#` CSV metadata header                                                     | utility        | file-I/O         | No metadata-header CSV precedent in-repo (RESEARCH §4, grep-verified). Format is `[ASSUMED]` — RESEARCH §4 gives the recommended layout.                                                                                          |

---

## Metadata

**Analog search scope:** `app/workflows/`, `app/components/`, `tests/`, `petringa/core/plotting.py`, `petringa/__init__.py`
**Files scanned:** 9 (home.py, cv.py, cce.py, field_map.py, main.py, device_sidebar.py, test_app_pages.py, test_app_device_sidebar.py, plotting.py)
**Pattern extraction date:** 2026-07-11
