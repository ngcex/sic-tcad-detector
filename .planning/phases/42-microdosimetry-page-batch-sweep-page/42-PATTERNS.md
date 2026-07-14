# Phase 42: Microdosimetry Page + Batch Sweep Page - Pattern Map

**Mapped:** 2026-07-14
**Files analyzed:** 6 (2 pages modified, 1 components file modified, 3 test files new/extended)
**Analogs found:** 6 / 6 (4 exact/role-match in-codebase; 2 mechanisms are net-new, sourced to RESEARCH.md)

> **UI-SPEC supersedes RESEARCH.md.** Where RESEARCH.md's code sketches conflict with `42-UI-SPEC.md`, this map carries the **UI-SPEC** version. Three specific overrides are flagged inline: (1) y_F/y_D readout uses `st.caption`, not `st.metric`; (2) `build_sweep_overlay_figure` axis titles come from a per-facade label map, not the literal `"(facade x-axis)"` placeholder; (3) microdosimetry single-result Download CSV is **included** (Q2 resolved), so `to_csv_bytes` gains a new `microdosimetry` branch.

## File Classification

| New/Modified File                       | Change               | Role                  | Data Flow                   | Closest Analog                                                              | Match Quality                          |
| --------------------------------------- | -------------------- | --------------------- | --------------------------- | --------------------------------------------------------------------------- | -------------------------------------- |
| `app/workflows/microdosimetry.py`       | modify placeholder   | page (Streamlit)      | file-I/O + request-response | `app/workflows/radiation_damage.py` (skeleton) + tempfile bridge (new)      | role-match (skeleton) + new mechanism  |
| `app/workflows/batch_sweep.py`          | modify placeholder   | page (Streamlit)      | batch (parametric sweep)    | `app/workflows/dark_current.py`                                             | exact (canonical ParametricSweep page) |
| `app/components/results.py`             | modify (4 additions) | builders + serializer | transform (pure)            | in-file: `build_damage_figure`, `build_dark_current_figure`, `to_csv_bytes` | exact (same module)                    |
| `tests/test_app_microdosimetry_page.py` | new                  | test                  | AppTest                     | `tests/test_app_radiation_damage_page.py`                                   | role-match                             |
| `tests/test_app_batch_sweep_page.py`    | new                  | test                  | AppTest                     | `tests/test_app_dark_current_page.py`                                       | exact (sweep-page test)                |
| `tests/test_app_csv_export.py`          | extend               | test                  | pure unit                   | its own per-sim_type tests                                                  | exact (same file)                      |

**results.py = 4 changes (not 3):** `build_microdosimetry_figure`, `build_sweep_overlay_figure`, `sweep_results_to_csv_bytes`, AND a new `microdosimetry` branch in `to_csv_bytes` (Q2).

## Pattern Assignments

### `app/workflows/batch_sweep.py` (page, batch)

**Analog:** `app/workflows/dark_current.py` — copy verbatim; it IS the canonical `ParametricSweep`-driven page. The batch sweep page generalizes it from a hard-coded `param="T"` sweep to a user-selected curated field + facade.

**Module-docstring + imports pattern** (`dark_current.py:39-45`):

```python
from __future__ import annotations

import numpy as np
import streamlit as st

import petringa
from app.components.results import build_sweep_overlay_figure, sweep_results_to_csv_bytes
```

`petringa.ParametricSweep` and the facades MUST be referenced as **module attributes** (`petringa.run_cce`, `getattr(petringa, ...)`), never `from petringa import ...`, so tests can `monkeypatch.setattr(petringa, "run_cce", fake)` while real `ParametricSweep.run()` executes — the seam proven in `tests/test_app_run_mockability.py` (39-01).

**Empty-state + 1D guard pattern** (`dark_current.py:48-60`) — keep BOTH on this page (UI-SPEC line 121 puts the 1D guard on batch sweep):

```python
def render() -> None:
    st.title("Batch Sweep")

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    if cfg.half_width_um is not None:
        st.warning(
            "These workflows are 1D-only. Set Dimensionality to 1D in the sidebar."
        )
        st.stop()
```

**Curated selectbox constants** (NEW — module-level, from UI-SPEC line 144/146 + RESEARCH.md Pattern 2):

```python
SWEEPABLE_FIELDS = [   # CURATED — numeric, 1D-facade-safe (excludes half_width_um, doping_profile, N_D)
    "epi_thickness_um", "substrate_thickness_um", "N_A",
    "N_D_junction", "N_D_bulk", "L_transition_um", "T", "area_cm2",
]
SIM_FACADES = {        # label -> module-attribute facade name (Q1: run_field/run_microdosimetry EXCLUDED)
    "CCE vs bias (run_cce)": "run_cce",
    "C-V (run_cv)": "run_cv",
    "CCE vs temperature (run_temperature_sweep)": "run_temperature_sweep",
}
```

**Input widgets** (UI-SPEC line 138-142 — note defaults):

```python
param = st.selectbox("Sweep parameter", SWEEPABLE_FIELDS, key="sweep_param")  # default epi_thickness_um
values_raw = st.text_input("Values (comma-separated)", value="10, 15, 20", key="sweep_values")
sim_label = st.selectbox("Simulation type", list(SIM_FACADES), key="sweep_sim")  # default "CCE vs bias (run_cce)"
```

**Value-list parse + Run + partial-failure aggregation** — mirror `dark_current.py:97-149`'s `try/except RuntimeError` + `T_ok`/`I_total` skip-empty loop, adapted per-swept-value (UI-SPEC Interaction Contract line 160-168):

```python
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
        results = petringa.ParametricSweep(          # REAL .run() — never hand-roll the loop
            base_config=cfg, param=param, values=values, sim_fn=sim_fn,
        ).run()
        # per-swept-value skip-empty aggregation (mirror dark_current.py:121-128):
        ok_values, ok_results = [], []
        for val, res in zip(values, results):
            if len(res.x) < 1:
                continue
            ok_values.append(val); ok_results.append(res)
        st.session_state["sweep_results"] = ok_results
        st.session_state["sweep_param"] = param
        st.session_state["sweep_values"] = ok_values
        st.session_state["sweep_sim_label"] = sim_label   # needed for per-facade axis titles at render
        st.session_state["sweep_n_ok"] = len(ok_results)
        st.session_state["sweep_n_requested"] = len(values)
    except RuntimeError as e:
        st.error(f"Simulation failed to converge: {e}\n\nTry a shallower value range.")
```

**Partial-failure banner + render + bulk download** — mirror `dark_current.py:151-167` (`n_ok < n_requested` warning above chart), extended with the `n_ok == 0` → `st.error` case (UI-SPEC line 167) and the bulk-CSV button label (UI-SPEC line 123, must read exactly `"Download all results as CSV"`):

```python
results = st.session_state.get("sweep_results")
if results is not None:
    n_ok = st.session_state.get("sweep_n_ok", 0)
    n_requested = st.session_state.get("sweep_n_requested", 0)
    if n_ok == 0:
        st.error("All swept values failed to converge. Try a shallower value range.")
    else:
        if n_ok < n_requested:
            st.warning(
                f"{n_ok} of {n_requested} values completed successfully; the rest "
                "failed to converge or returned no data and are omitted from the plot below."
            )
        param = st.session_state["sweep_param"]
        values = st.session_state["sweep_values"]
        sim_label = st.session_state["sweep_sim_label"]
        st.plotly_chart(build_sweep_overlay_figure(results, param, values, sim_label))
        st.download_button(
            "Download all results as CSV",
            data=sweep_results_to_csv_bytes(results, param, values),
            file_name="batch_sweep_result.csv",
            mime="text/csv",
        )
```

> **Convergence-safe default (UI-SPEC line 156 / RESEARCH Pitfall 3):** default `sim_label="CCE vs bias (run_cce)"` + `param="epi_thickness_um"` + `values="10, 15, 20"`. `run_cce` (`v_stop=-40`) stays inside the DD convergence envelope AND truncates gracefully. Do NOT default to `run_field`. A planner Wave-0 live-devsim spike should confirm this exact combo renders 3 clean/truncated curves.

---

### `app/workflows/microdosimetry.py` (page, file-I/O + request-response)

**Analog:** `app/workflows/radiation_damage.py` for the overall Run→cache→render→download **skeleton** and the empty-state guard (kept per Q3, UI-SPEC line 21). The `tempfile` upload bridge itself is a **new mechanism** — see No Analog Found.

**COPY from `radiation_damage.py`, but explicitly DROP three things** (do NOT cargo-cult the whole file):

- **DROP the 1D-only `half_width_um` guard** (`radiation_damage.py:54-58`) — `run_microdosimetry` is a pure pipeline that touches no devsim; the spectrum is config-independent (UI-SPEC line 21). No dimensionality guard on this page.
- **DROP the NaN-tolerance handling** (`radiation_damage.py:114-119`) — `run_microdosimetry` has no partial-convergence failure mode (Plot Contract, UI-SPEC line 183). No `np.isnan` check.
- **DROP the kappa data-blocked banner** (`radiation_damage.py:39-47`) — that is radiation-damage-specific.

**KEEP from `radiation_damage.py`:** the module docstring "no module-level side effects" note, the module-attribute import seam, the empty-state guard, the `if st.button("Run simulation"):` → `try/except` → cache-in-session_state → `if result is not None:` render+download shape.

**Imports** (module-attribute seam — RESEARCH.md line 199):

```python
from __future__ import annotations

import os
import tempfile

import streamlit as st

import petringa
from app.components.results import build_microdosimetry_figure, to_csv_bytes
```

**Empty-state guard (kept, Q3)** — same as `radiation_damage.py:49-52`:

```python
cfg = st.session_state.get("device_config")
if cfg is None:
    st.info("Configure a device in the sidebar to begin.")
    st.stop()
```

**Uploader + SV-geometry inputs + format hint** (UI-SPEC line 115/132-134):

```python
uploaded = st.file_uploader("Upload MC events CSV", type=["csv"], key="micro_csv")
st.caption(
    "CSV columns: event_id, x, y, z, edep — positions in cm, energy deposit in keV. "
    "One row per MC step; steps are summed per event_id."
)
sv_thickness = st.number_input("Sensitive-volume thickness (µm)", value=10.0, key="micro_sv_t")
sv_width = st.number_input("Sensitive-volume width (µm)", value=150.0, key="micro_sv_w")
```

**Run + tempfile bridge + malformed-CSV guard** (NEW bridge — RESEARCH.md Pattern 1, lines 178-197; `run_microdosimetry` signature verified `petringa/api/simulation.py:843`, `mc_csv_path: str`):

```python
if st.button("Run simulation"):
    if uploaded is None:
        st.warning("Upload an MC events CSV to run.")
        st.stop()
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
            tmp.write(uploaded.getvalue())   # UploadedFile bytes; server-generated path (no traversal)
            tmp_path = tmp.name
        st.session_state["microdosimetry_result"] = petringa.run_microdosimetry(
            cfg, mc_csv_path=tmp_path,
            sv_thickness_um=sv_thickness, sv_width_um=sv_width,
        )
    except (ValueError, KeyError) as e:       # malformed CSV -> load_mc_events_csv raises
        st.error(f"Could not parse the uploaded CSV: {e}")
    finally:
        if tmp_path is not None and os.path.exists(tmp_path):
            os.remove(tmp_path)               # single-request temp lifetime; no disk accumulation
```

**Render + y_F/y_D readout + single-result download** — the readout uses **`st.caption`, NOT `st.metric`** (UI-SPEC Typography Discretion line 64 explicitly rejects `st.metric`; RESEARCH.md's `st.metric` sketch is superseded). Download button included (Q2):

```python
result = st.session_state.get("microdosimetry_result")
if result is not None:
    st.plotly_chart(build_microdosimetry_figure(result))
    st.caption(
        f"y_F = {result.metadata['y_F']:.2f} keV/µm    "
        f"y_D = {result.metadata['y_D']:.2f} keV/µm"
    )
    st.download_button(
        "Download CSV",
        data=to_csv_bytes(result),
        file_name="microdosimetry_result.csv",
        mime="text/csv",
    )
```

---

### `app/components/results.py` (builders + serializer, pure transform)

All four additions are PURE (no `st.*`), matching the module's existing contract (`results.py:1-11`). Imports already present: `petringa`, `numpy as np`, `pandas as pd`, `plotly.graph_objects as go`, `datetime`/`timezone`, `SimResult`.

**(1) `build_microdosimetry_figure` — analog `build_damage_figure` (`results.py:63-86`)** (single trace, log-x, `#1F6FEB`):

```python
def build_microdosimetry_figure(result: SimResult) -> go.Figure:
    """y·d(y) vs lineal energy y (LOG x-axis per ICRU-36; x spans 0.01–9772 keV/µm)."""
    fig = go.Figure(
        data=go.Scatter(x=result.x, y=result.y, mode="lines", line=dict(color="#1F6FEB"))
    )
    fig.update_layout(
        title="Microdosimetric Spectrum",
        xaxis_title="Lineal energy y (keV/µm)",
        xaxis_type="log",              # MANDATORY (UI-SPEC Plot Contract line 178)
        yaxis_title="y · d(y)",
    )
    return fig
```

**(2) `build_sweep_overlay_figure` — analog `build_dark_current_figure` (`results.py:89-132`)** (multi-trace palette loop). **SIGNATURE GAP:** RESEARCH.md's sketch is `(results, param, values)` (3 args) with a literal `"(facade x-axis)"` placeholder. The UI-SPEC per-facade axis map (line 193-197) requires knowing which facade ran, so the planner MUST add a **4th `sim_label` (or facade key) argument** and select axis titles from a map. Palette is cycled by trace order `i % len(PALETTE)` (value-keyed qualitative — NOT Phase 41's fixed-quantity mapping; UI-SPEC line 99):

```python
_SWEEP_PALETTE = ["#1F6FEB", "#D32F2F", "#2E7D32", "#1A1A1A", "#9AA0A6"]  # cycle if >5

_SWEEP_AXIS_TITLES = {   # keyed by SIM_FACADES label (UI-SPEC line 193-197)
    "CCE vs bias (run_cce)": ("Bias V (V)", "Charge Collection Efficiency"),
    "C-V (run_cv)": ("Bias V (V)", "Capacitance (F)"),
    "CCE vs temperature (run_temperature_sweep)": ("Temperature (K)", "Value"),
}

def build_sweep_overlay_figure(results, param, values, sim_label) -> go.Figure:
    """One trace per swept value; legend f"{param}={val}"; per-facade axis titles."""
    fig = go.Figure()
    for i, (val, res) in enumerate(zip(values, results)):
        fig.add_trace(go.Scatter(
            x=res.x, y=res.y, mode="lines+markers",
            name=f"{param}={val}",
            line=dict(color=_SWEEP_PALETTE[i % len(_SWEEP_PALETTE)]),
        ))
    x_title, y_title = _SWEEP_AXIS_TITLES.get(sim_label, ("(facade x-axis)", "Value"))
    fig.update_layout(title=f"Parametric Sweep: {param}",
                      xaxis_title=x_title, yaxis_title=y_title)
    return fig
```

**(3) `sweep_results_to_csv_bytes` — analog `to_csv_bytes` header convention (`results.py:234-248`)** — a SEPARATE bulk serializer (NOT a `to_csv_bytes` branch — RESEARCH Pitfall 2). One CSV, leading run-identifier column named after `param`:

```python
def sweep_results_to_csv_bytes(results, param, values) -> bytes:
    """ALL sweep runs in ONE CSV; leading `<param>` column identifies each run."""
    frames = [pd.DataFrame({param: val, "x": res.x, "y": res.y})
              for val, res in zip(values, results)]
    combined = pd.concat(frames, ignore_index=True)
    header_lines = [
        f"# petringa SiC TCAD Simulator — parametric sweep ({param})",
        f"# software_version: {petringa.__version__}",
        f"# generated: {datetime.now(timezone.utc).isoformat()}",
        f"# swept_values: {list(values)}",
    ]
    return ("\n".join(header_lines) + "\n" + combined.to_csv(index=False)).encode("utf-8")
```

**(4) NEW `microdosimetry` branch in `to_csv_bytes` — analog the existing per-sim_type branches (`results.py:169-232`)** (Q2 requires it; RESEARCH left it open). Insert a new `elif result.sim_type == "microdosimetry":` block before the final `else` raise, following the exact `df` + `extra_header_lines` shape of the sibling branches. Suggested columns `y_keV_per_um` (`result.x`), `y_times_d_y` (`result.y`); header lines carry `y_F`, `y_D`, `l_bar_um` from `result.metadata`. The shared `# device:` header + `datetime` footer (`results.py:234-248`) apply unchanged.

---

### `tests/test_app_batch_sweep_page.py` (test, new)

**Analog:** `tests/test_app_dark_current_page.py` — near-exact structural template.

**Imports + module-attribute monkeypatch seam** (`test_app_dark_current_page.py:21-27`):

```python
from __future__ import annotations
import numpy as np
from streamlit.testing.v1 import AppTest
import petringa
from petringa import DeviceConfig, SimResult
```

**Fake facade** — mirror `_fake_run_dark_current` (`test_app_dark_current_page.py:30-43`), but shaped for the swept facade (e.g. `run_cce` returns an x/y curve). Monkeypatch `petringa.run_cce` (the facade), NEVER `ParametricSweep` — the real `.run()` must execute (asserts genuine sweep wiring, `test_app_dark_current_page.py:8-13`).

**Page wrapper + tests to mirror:**

- `_run_batch_sweep_page()` importing `app.workflows.batch_sweep.render` (`test_app_dark_current_page.py:46-49`)
- empty-state guard (`test_app_dark_current_page.py:52-61`)
- 2D-config 1D-only warning (`test_app_dark_current_page.py:64-72`)
- Run uses real `ParametricSweep` + caches `list[SimResult]`, ≥3 values → ≥3 results (`test_app_dark_current_page.py:75-99`) — assert `len(at.session_state["sweep_results"]) >= 3`
- partial-value-failure warning, not crash (`test_app_dark_current_page.py:102-140`) — one fake returns empty `x` → `n_ok < n_requested` banner
- bad value-list → `st.error`, no crash (NEW vs dark_current — set `sweep_values` to non-numeric via `at.text_input`, assert `at.error`)
- facade `RuntimeError` → `st.error`, no cache (`test_app_dark_current_page.py:142-159`)

Assertions limited to `at.exception`, `at.session_state`, `at.button`, `at.warning`, `at.info`, `at.error` (no `plotly_chart`/`download_button` accessor — `test_app_dark_current_page.py:15-19`).

---

### `tests/test_app_microdosimetry_page.py` (test, new)

**Analog:** `tests/test_app_radiation_damage_page.py` for structure (fake facade, `_run_..._page()` wrapper, monkeypatch seam, assertion set). The **`file_uploader.upload((name, bytes, mime))` driving is a new mechanism** — see No Analog Found.

Mirror `test_app_radiation_damage_page.py`:

- module-attribute monkeypatch of `petringa.run_microdosimetry` (`test_app_radiation_damage_page.py:1-12`)
- `_fake_run_microdosimetry(cfg, mc_csv_path, sv_thickness_um=10.0, sv_width_um=150.0)` returning a `SimResult(sim_type="microdosimetry", x=..., y=..., metadata={"y_F":..., "y_D":..., "l_bar_um":...})` (shape from `test_app_radiation_damage_page.py:23-35`)
- empty-state guard test (`test_app_radiation_damage_page.py:72-81`)
- no-file-on-Run → `st.warning("Upload an MC events CSV to run.")` (NEW — assert `at.warning`)
- happy-path: inject `data/synthetic_mc_events.csv` bytes via `at.file_uploader[0].upload(...)`, click Run, assert `at.session_state["microdosimetry_result"].sim_type == "microdosimetry"` (RESEARCH Pitfall 1 — AppTest **1.58** file_uploader IS drivable)
- malformed CSV → upload junk bytes → `st.error`, no crash (mirror the RuntimeError-not-crash shape `test_app_radiation_damage_page.py:120-139`, but catching the parse `st.error`)
- **NO 1D-guard test and NO NaN-tolerance test** — those mechanisms are dropped from this page (see page section above).

---

### `tests/test_app_csv_export.py` (test, extend)

**Analog:** its own existing per-sim_type tests (`test_app_csv_export.py:24-198`) — pure `SimResult` fixture → `to_csv_bytes` → assert columns/header. Add:

- **`test_microdosimetry_csv_columns_and_header`** — mirror `test_damage_csv_columns_and_header` (`test_app_csv_export.py:145-172`): hand-build a `sim_type="microdosimetry"` `SimResult`, assert the new branch's columns + the `y_F`/`y_D`/`l_bar_um` header lines.
- **`test_sweep_results_to_csv_bytes_shape`** — pure test of the NEW bulk serializer (no Streamlit): build N `SimResult`s, assert ONE CSV, N runs concatenated, run-identifier column named after `param` present (RESEARCH Test Map line 458, 471). Import `sweep_results_to_csv_bytes` alongside the existing `to_csv_bytes`/`build_dark_current_figure` imports (`test_app_csv_export.py:17`).

## Shared Patterns

### Module-attribute mockability seam

**Source:** `app/workflows/dark_current.py:32-36` + `tests/test_app_run_mockability.py` (39-01)
**Apply to:** both new pages + both new page tests
`import petringa; petringa.run_X(...)` (never `from petringa import run_X`) so tests `monkeypatch.setattr(petringa, "run_X", fake)`. For batch sweep, monkeypatch the **facade** (`run_cce`), never `ParametricSweep` — the real `.run()` must execute.

### Run → cache → render → download shape

**Source:** `app/workflows/radiation_damage.py:96-126` and `app/workflows/dark_current.py:97-167`
**Apply to:** both pages
`if st.button("Run simulation"):` → `try/except RuntimeError` → `st.session_state[key] = result` → later `if result is not None:` → `st.plotly_chart(...)` + `st.download_button(...)`.

### Empty-state guard

**Source:** `app/workflows/radiation_damage.py:49-52` (identical across all pages)
**Apply to:** both pages (microdosimetry keeps it per Q3, even though config is provenance-only)

```python
cfg = st.session_state.get("device_config")
if cfg is None:
    st.info("Configure a device in the sidebar to begin.")
    st.stop()
```

### Partial-failure skip-empty aggregation + banner

**Source:** `app/workflows/dark_current.py:121-160` (`T_ok`/`I_total` loop + `n_ok < n_requested` warning)
**Apply to:** batch sweep page (per-swept-value, extended with the `n_ok == 0` → `st.error` case, UI-SPEC line 167)

### CSV metadata-header convention

**Source:** `app/components/results.py:234-248` (`#`-comment lines: title, `software_version`, ISO-8601 `generated`, `device:` fields)
**Apply to:** the new `microdosimetry` branch of `to_csv_bytes` AND `sweep_results_to_csv_bytes` (the latter omits `device:`, adds `swept_values:`)

### Pure-builder / pure-serializer contract

**Source:** `app/components/results.py:1-11` (no `st.*` anywhere in this module)
**Apply to:** all 4 results.py additions — keeps them testable without a Streamlit runtime or devsim build (`tests/test_app_csv_export.py`)

### AppTest assertion surface

**Source:** `tests/test_app_dark_current_page.py:15-19`
**Apply to:** both new page tests — assert on `at.exception`, `at.session_state`, `at.button`, `at.warning`, `at.info`, `at.error` (no `plotly_chart`/`download_button` accessor). Exception: `file_uploader` IS drivable in 1.58 (below).

## No Analog Found

Two genuinely new mechanisms — no codebase precedent. Planner should use RESEARCH.md rather than force-fit an analog:

| Mechanism                                                        | Where                                   | Reason no analog                                                                                                                                                                                         | Source                                                             |
| ---------------------------------------------------------------- | --------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| `tempfile.NamedTemporaryFile` upload → `mc_csv_path: str` bridge | `app/workflows/microdosimetry.py`       | First file-upload surface in the app; no prior page reads an uploaded file. Server-generated temp path (no traversal), `os.remove` in `finally`.                                                         | RESEARCH.md Pattern 1 (lines 159-197); threat T-37-02-V5           |
| `at.file_uploader[0].upload((name, bytes, mime))` test-driving   | `tests/test_app_microdosimetry_page.py` | Prior page tests (39/40/41) documented AppTest 1.55 could not reach result widgets; file-upload driving is a **1.58** capability with no in-repo precedent. Inject `data/synthetic_mc_events.csv` bytes. | RESEARCH.md Pitfall 1 (lines 270-275); State of the Art (line 391) |

## Metadata

**Analog search scope:** `app/workflows/`, `app/components/`, `tests/`, `petringa/api/` (signature verification only)
**Files scanned:** dark_current.py, radiation_damage.py, microdosimetry.py, batch_sweep.py (workflows); results.py (components); test_app_csv_export.py, test_app_dark_current_page.py, test_app_radiation_damage_page.py (tests); sweep.py, simulation.py (signatures)
**Signatures verified:** `ParametricSweep(base_config, param, values, sim_fn, sim_kwargs).run() -> list[SimResult]` (`petringa/api/sweep.py:31-70`); `run_microdosimetry(config, mc_csv_path: str, sv_thickness_um=10.0, sv_width_um=150.0) -> SimResult` (`petringa/api/simulation.py:843-848`)
**Pattern extraction date:** 2026-07-14
