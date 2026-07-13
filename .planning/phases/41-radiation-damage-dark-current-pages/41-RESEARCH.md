# Phase 41: Radiation Damage + Dark Current Pages - Research

**Researched:** 2026-07-13
**Domain:** Streamlit UI wiring over existing `petringa.run_radiation_damage()` / `petringa.run_dark_current()` facades — no new physics/core code.
**Confidence:** HIGH (facade behavior empirically verified this session by actually calling both facades against a live devsim); MEDIUM on the "vs temperature" wording question (flagged as an open question requiring user confirmation).

## Summary

Phase 41 is pure UI wiring: `petringa/api/simulation.py` already has fully-implemented,
tested `run_radiation_damage()` (line 447) and `run_dark_current()` (line 665) facades.
`app/main.py` already imports and registers `render_radiation_damage` / `render_dark_current`
from `app/workflows/radiation_damage.py` and `app/workflows/dark_current.py` — both files
currently exist ONLY as placeholders ("Running this simulation is implemented in Phase 41.").
This phase replaces those two placeholder bodies following the exact
Run→cache→render→download skeleton established in Phase 39 (`cv.py`/`cce.py`/`field_map.py`)
and the shared pure-figure-builder module `app/components/results.py`.

The single largest risk in this phase is **not** convergence (though that matters too) but
**silent partial-failure return shapes that differ from Phase 39's precedent**. Empirically
verified this session (see Pitfalls below): `run_dark_current` at its own defaults
(`v_stop=-100.0`) truncates the returned arrays with **no `"truncated"` metadata key at all**
(unlike `run_cv`/`run_cce`, which explicitly report `truncated`/`requested_v_stop`).
`run_radiation_damage` at its own defaults returns a **NaN embedded mid-array** with the
full-length `x` intact and **no exception raised** — a third, distinct failure mode from
Phase 39's `RuntimeError`. Both pages must be built to tolerate these shapes, not just
mirror Phase 39's `try/except RuntimeError` pattern verbatim.

The success-criterion wording "dark current vs temperature" does not match `run_dark_current`'s
actual signature (which sweeps bias, `x=voltages`) — this is flagged as Open Question #1 and
needs discuss-phase / user confirmation before the planner locks page architecture.

**Primary recommendation:** Build both pages as thin `render()` modules mirroring
`app/workflows/cce.py` exactly (guard → widget inputs → Run button → `try/except RuntimeError`
→ session_state cache → pure figure builder from `results.py` → CSV download), with two new
pure builders (`build_damage_figure`, `build_dark_current_figure`) and two new CSV branches
added to `to_csv_bytes`. Choose widget defaults that diverge from the facades' raw defaults
where the raw defaults are known to hit partial non-convergence (see Pitfall 2), and add a
`st.warning` for the kappa-data-blocked banner unconditionally at the top of the radiation
damage page.

## Architectural Responsibility Map

| Capability                                                                                               | Primary Tier                                                      | Secondary Tier | Rationale                                                                                               |
| -------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- | -------------- | ------------------------------------------------------------------------------------------------------- |
| Radiation-damage physics (κ(E), fluence sweep, CCE degradation)                                          | API/Backend (`petringa/core/`, `petringa/api/simulation.py`)      | —              | Already implemented; Phase 41 must not touch it                                                         |
| Dark-current physics (SRH/TAT/SRV decomposition, bias sweep)                                             | API/Backend (`petringa/core/dark_current.py`, `run_dark_current`) | —              | Already implemented; Phase 41 must not touch it                                                         |
| Widget inputs (fluence range, proton energy, dark-current v_start/v_stop/n_points, optional N_t/S_n/S_p) | Frontend Server (Streamlit `render()`)                            | —              | Page-local `st.*` widgets, not sidebar (sidebar is only for `DeviceConfig`, per Phase 38/39 convention) |
| Kappa-data-blocked warning banner                                                                        | Frontend Server (Streamlit `render()`)                            | —              | Static/unconditional `st.warning`, no data dependency                                                   |
| Result caching across reruns                                                                             | Frontend Server (`st.session_state`)                              | —              | `DeviceConfig` is unhashable so `st.cache_data` is impossible (established Phase 39)                    |
| Figure construction (CCE vs fluence, dark current decomposition)                                         | Frontend Server (`app/components/results.py`, pure)               | —              | Pure Plotly builders, no `st.*`, unit-testable without Streamlit runtime (established Phase 39/40)      |
| CSV export                                                                                               | Frontend Server (`app/components/results.py::to_csv_bytes`)       | —              | Extends the existing dispatcher; must not regress the `test_unknown_sim_type_raises_value_error` test   |

## Phase Requirements

| ID      | Description                                                                                                                                    | Research Support                                                                                                                                                                                                                                |
| ------- | ---------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FEAT-01 | User can run a radiation damage sweep specifying fluence range and proton energy, see CCE vs fluence curves; kappa-data-blocked warning banner | `run_radiation_damage()` signature/defaults documented below; kappa warning text sourced verbatim from the facade's own docstring/module comments; NaN-in-array failure mode empirically verified                                               |
| FEAT-02 | User can run a dark current simulation and see dark current vs temperature/bias with J_SRH/J_TAT/J_SRV decomposed and overlaid                 | `run_dark_current()` signature/defaults documented below; decomposition trace spec (log-abs, zero-guard) sourced from the existing `plot_dark_current_decomposition` matplotlib reference; "vs temperature" wording flagged as Open Question #1 |

## User Constraints

No CONTEXT.md exists for this phase (checked `.planning/phases/41-radiation-damage-dark-current-pages/` — empty prior to this research). No locked decisions, discretion areas, or deferred ideas to copy verbatim. If `/gsd:discuss-phase` runs before planning, Open Question #1 below should be the primary discussion topic.

## Project Constraints (from CLAUDE.md)

No `./CLAUDE.md` exists in the repository root. No project-specific directives to enforce beyond what STATE.md/ROADMAP.md/REQUIREMENTS.md already establish (no new physics in v5.0; UI-wiring-only for this phase).

No `.claude/skills/` or `.agents/skills/` directory with `SKILL.md` files was found (only `.claude/*.workflow.js` audit scripts, unrelated to this phase).

## Standard Stack

### Core

| Library                           | Version                            | Purpose                                                      | Why Standard                                                            |
| --------------------------------- | ---------------------------------- | ------------------------------------------------------------ | ----------------------------------------------------------------------- |
| `streamlit`                       | already pinned in `pyproject.toml` | Page rendering, widgets, session_state                       | Established in Phases 38-40; no new dependency                          |
| `plotly` (`plotly.graph_objects`) | already pinned                     | `go.Figure` construction in `results.py` builders            | Established in Phase 39; `app/components/results.py` already imports it |
| `pandas`                          | already pinned                     | `to_csv_bytes` DataFrame construction                        | Established in Phase 39                                                 |
| `numpy`                           | already pinned                     | Array construction for widget-derived fluence/voltage ranges | Already used throughout `petringa/api/simulation.py`                    |

No new packages are needed for this phase — it consumes facades and libraries already declared in `pyproject.toml`. **Package Legitimacy Audit is not applicable** (no new installs).

### Alternatives Considered

| Instead of                                                   | Could Use                                                      | Tradeoff                                                                                                                                                                                                                                                                               |
| ------------------------------------------------------------ | -------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `st.number_input` pairs for fluence min/max + `np.geomspace` | `st.slider` with log formatting                                | `st.slider` has no native log-scale mode in Streamlit 1.5x; two `number_input`s + explicit `np.geomspace(min, max, n)` construction (mirrors `run_radiation_damage`'s own default `np.geomspace(1e13, 1e16, 6)`) is simpler and matches the facade's own array-construction convention |
| `try/except RuntimeError` alone (Phase 39 pattern)           | `try/except RuntimeError` + explicit NaN/short-array tolerance | Verified this session that RuntimeError alone is insufficient for these two facades (see Pitfalls 2/3) — must additionally handle graceful partial-failure return shapes, not just hard exceptions                                                                                     |

## Package Legitimacy Audit

Not applicable — this phase installs no new packages. All imports (`streamlit`, `plotly`, `pandas`, `numpy`, `petringa`) are already present in `pyproject.toml` and used by prior phases.

## Architecture Patterns

### System Architecture Diagram

```
User (browser)
   |
   v
app/main.py (st.navigation)  --------------------------------+
   |  renders sidebar (device_sidebar.py) -> session_state["device_config"]
   |
   +--> app/workflows/radiation_damage.py :: render()
   |        |
   |        [1] st.warning(kappa banner)  <-- UNCONDITIONAL, every rerun
   |        [2] widget inputs: fluence_min, fluence_max, n_points, proton_energy_MeV, V_bias
   |        [3] st.button("Run simulation")
   |             -> petringa.run_radiation_damage(cfg, fluences=np.geomspace(...), V_bias=..., proton_energy_MeV=...)
   |             -> try/except RuntimeError -> st.error(...)
   |             -> session_state["damage_result"] = SimResult(sim_type="damage")
   |        [4] result = session_state.get("damage_result")
   |             -> build_damage_figure(result)  [app/components/results.py, PURE]
   |             -> st.plotly_chart(fig)
   |             -> st.download_button(to_csv_bytes(result))
   |
   +--> app/workflows/dark_current.py :: render()
            |
            [1] widget inputs: v_start, v_stop, n_points, (optional expander: N_t, S_n, S_p)
            [2] st.button("Run simulation")
                 -> petringa.run_dark_current(cfg, v_start=..., v_stop=..., n_points=..., N_t=..., S_n=..., S_p=...)
                 -> try/except RuntimeError -> st.error(...)
                 -> session_state["dark_current_result"] = SimResult(sim_type="dark_current")
            [3] result = session_state.get("dark_current_result")
                 -> build_dark_current_figure(result)  [app/components/results.py, PURE]
                      overlays I_total, I_SRH, I_TAT, I_SRV (abs, log-y, zero-guarded)
                 -> st.plotly_chart(fig)
                 -> st.download_button(to_csv_bytes(result))

app/components/results.py (PURE, no st.* calls)
   build_damage_figure(result) -> go.Figure          [NEW]
   build_dark_current_figure(result) -> go.Figure    [NEW]
   to_csv_bytes(result) -> bytes                     [EXTENDED: + "damage", "dark_current" branches]

petringa/api/simulation.py (UNCHANGED, already implemented)
   run_radiation_damage(config, fluences=None, V_bias=-40.0, proton_energy_MeV=5.6) -> SimResult
   run_dark_current(config, v_start=0.0, v_stop=-100.0, n_points=20, N_t=None, S_n=None, S_p=None) -> SimResult
```

### Recommended Project Structure

No new files/directories beyond what already exists:

```
app/
├── workflows/
│   ├── radiation_damage.py   # REPLACE placeholder body (FEAT-01)
│   └── dark_current.py       # REPLACE placeholder body (FEAT-02)
├── components/
│   └── results.py            # EXTEND: build_damage_figure, build_dark_current_figure, +2 to_csv_bytes branches
tests/
├── test_app_radiation_damage_page.py   # NEW (mirrors test_app_field_page.py / test_app_cce_page.py)
├── test_app_dark_current_page.py       # NEW
└── test_app_csv_export.py              # EXTEND: damage + dark_current CSV tests; FIX test_unknown_sim_type_raises_value_error (see Pitfall 4)
```

### Pattern 1: Run → cache → render → download (Phase 39 skeleton, verbatim)

**What:** Every result page follows the identical 5-step shape: empty-state guard → optional pre-check → Run button calling the facade referenced as a **module attribute** (`petringa.run_x`, never `from petringa import run_x`) → `session_state` cache → pure-builder render → CSV download.

**When to use:** Both new pages, no exceptions.

**Example (from `app/workflows/cce.py`, copy this shape exactly):**

```python
# Source: app/workflows/cce.py (Phase 39, verbatim structure)
from __future__ import annotations

import streamlit as st

import petringa
from app.components.results import build_cce_figure, to_csv_bytes


def render() -> None:
    st.title("Charge Collection (CCE)")

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    if cfg.half_width_um is not None:
        st.warning("These workflows are 1D-only. ...")
        st.stop()

    if st.button("Run simulation"):
        try:
            st.session_state["cce_result"] = petringa.run_cce(cfg)
        except RuntimeError as e:
            st.error(f"Simulation failed to converge: {e}\n\n...")

    result = st.session_state.get("cce_result")
    if result is not None:
        st.plotly_chart(build_cce_figure(result))
        st.download_button("Download CSV", data=to_csv_bytes(result),
                            file_name="cce_result.csv", mime="text/csv")
```

Both `run_radiation_damage` and `run_dark_current` are **1D-only** implicitly: `cce_vs_fluence`
(called by `run_radiation_damage`) uses `create_dd_device` (1D), and `create_dark_current_device`
(called by `run_dark_current`) also uses `create_dd_device`. Neither raises `NotImplementedError`
for 2D configs the way `run_cv`/`run_cce` do — **empirically unverified this session whether they
silently mishandle 2D `half_width_um`** (not tested; out of scope to test since `DeviceConfig`
sidebar always allows 2D). **Recommendation: add the same `cfg.half_width_um is not None` 1D-only
pre-check guard used in `cv.py`/`cce.py`** as a defensive measure, even though it is not proven
necessary — the cost of adding it is one `if` block, the cost of NOT adding it is an unverified
2D failure mode reaching users. Flag this as a pre-check to add, following the `cv.py`/`cce.py`
precedent rather than `field_map.py`'s (which deliberately removed its 2D guard because `run_field`
explicitly supports 2D — that is NOT the case here).

### Pattern 2: Widget inputs for array-shaped parameters (NEW — no direct analog)

**What:** `run_radiation_damage`'s `fluences` parameter and `run_dark_current`'s bias sweep are
both **arrays** constructed from simpler scalar inputs (`np.geomspace`/`np.linspace`), not
raw arrays a user would type in directly.

**When to use:** Both new pages need a "range → array" widget group.

**Radiation damage page — fluence range (mirrors `run_radiation_damage`'s own default
construction, `np.geomspace(1e13, 1e16, 6)`, confirmed in `petringa/api/simulation.py:496`):**

```python
col1, col2, col3 = st.columns(3)
with col1:
    fluence_min = st.number_input("Min fluence (p/cm²)", value=1e13, format="%.3e")
with col2:
    fluence_max = st.number_input("Max fluence (p/cm²)", value=1e16, format="%.3e")
with col3:
    n_points = st.number_input("Number of points", value=6, min_value=2, step=1)

proton_energy_MeV = st.selectbox(
    "Proton energy (MeV)", [30, 62, 70, 150], index=1,
)  # matches NIEL_HARDNESS_PROTON_SIC table keys exactly (petringa/core/radiation_damage.py:61-66)

V_bias = st.number_input("Reverse bias V_bias (V)", value=-40.0)

fluences = np.geomspace(fluence_min, fluence_max, int(n_points))
```

**Widget default rationale:** `proton_energy_MeV` as a `selectbox` (not free `number_input`) is
recommended because `get_hardness_factor` (petringa/core/radiation_damage.py:407) only has real
tabulated entries at {30, 62, 70, 150} MeV — any other value silently linearly interpolates
between table entries via `np.interp`, which is fine numerically but misleadingly precise given
the kappa values are placeholders. A `selectbox` of the 4 real table keys communicates this
implicitly; a free `number_input` (like `run_radiation_damage`'s own `proton_energy_MeV=5.6`
default, which is BELOW the table's lowest key of 30 and therefore clamps to the 30 MeV value
via `np.interp`'s edge behavior) invites false precision. **This is a UX judgment call, not a
hard requirement** — `number_input` with the 5.6 MeV facade default is also acceptable if the
planner prefers matching the facade default exactly; flag either choice to discuss-phase.

**Dark current page — bias sweep (mirrors `run_dark_current`'s own parameters exactly, they are
already scalars):**

```python
col1, col2, col3 = st.columns(3)
with col1:
    v_start = st.number_input("V start (V)", value=0.0)
with col2:
    v_stop = st.number_input("V stop (V)", value=-50.0)  # NOTE: widget default -50.0, NOT
                                                           # the facade's raw default -100.0 —
                                                           # see Pitfall 2 for why
with col3:
    n_points = st.number_input("Number of points", value=20, min_value=2, step=1)

with st.expander("Advanced (trap/surface parameters)"):
    N_t = st.number_input("N_t override (cm⁻³/s, blank = default)", value=None, ...)
    S_n = st.number_input("S_n override (cm/s, blank = default)", value=None, ...)
    S_p = st.number_input("S_p override (cm/s, blank = default)", value=None, ...)
```

### Pattern 3: Persistent warning banner (NEW pattern for this repo — no prior "persistent" banner exists)

**What:** A banner that must appear on every rerun of the radiation damage page, not gated
behind `if result is not None`.

**When to use:** Radiation damage page only (FEAT-01 criterion #2).

**Example:**

```python
def render() -> None:
    st.title("Radiation Damage")

    st.warning(
        "**Data-blocked placeholder:** kappa (NIEL hardness factor) values used here "
        "are unvalidated placeholders (see petringa/core/radiation_damage.py "
        "NIEL_HARDNESS_PROTON_SIC table). The energy TREND is physically motivated but "
        "no ABSOLUTE Phi_crit or defect-concentration number is citable until real "
        "SR-NIEL SiC proton NIEL data replaces these placeholders. Treat CCE-vs-fluence "
        "curves below as a relative sensitivity shape only."
    )

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()
    # ... rest of page
```

Because `render()` re-executes on every Streamlit rerun (button click, widget change,
navigation), placing `st.warning(...)` unconditionally at the top of the function — BEFORE the
empty-state guard, so it is visible even with no device configured — makes it "persistent"
with zero `session_state` bookkeeping required. This is simpler than the advisor's suggestion
to place it merely "outside the `if result is not None` block"; placing it before even the
`cfg is None` guard ensures FEAT-01 criterion #2 ("displays a persistent warning banner") holds
under every possible page state, including first load before any device is configured.

**Text sourcing:** Draw the banner text verbatim from two authoritative in-repo sources —
`run_radiation_damage`'s own docstring Warning section (petringa/api/simulation.py:486-493,
"RESEARCH Pitfall 4: the NIEL kappa hardness factors ... are DATA-BLOCKED placeholders") and
the `NIEL_HARDNESS_PROTON_SIC` module-level comment block (petringa/core/radiation_damage.py:46-59,
"AUDIT C-5 ... DATA-BLOCKED, NOT YET FIXED"). Both are `[CITED: petringa/api/simulation.py,
petringa/core/radiation_damage.py]` — this is in-repo authoritative documentation, not an
external source, but it establishes the exact scientific caveat that must reach the UI verbatim
per FEAT-01's explicit requirement ("kappa values are data-blocked placeholders ... absolute
Phi_crit numbers are unvalidated").

### Anti-Patterns to Avoid

- **Treating `try/except RuntimeError` as sufficient error handling for these two facades:**
  Phase 39's pattern catches hard solver-convergence exceptions. Verified this session that
  BOTH new facades have a DIFFERENT, more common failure mode that does NOT raise: silent
  truncation (`run_dark_current`) and silent NaN-in-array (`run_radiation_damage`). Both must
  additionally be handled at the render layer (Plotly renders NaN as a gap automatically; short
  arrays render fine but should not be silently unexplained to the user).
- **Reusing the facades' raw numeric defaults as widget defaults without checking convergence:**
  `v_stop=-100.0` (dark current) and the implicit chain through `cce_vs_fluence` are both known,
  empirically, to hit partial failure at their own stated defaults for the plain default
  `DeviceConfig()`. Pick widget defaults that converge cleanly on first Run (see Pitfall 2).
- **Building a truly free-text array input for `fluences`:** the facade already provides a
  sensible `np.geomspace` default; re-deriving a log-spaced array from 2 endpoints + point-count
  widgets is simpler and safer than parsing a free-text list.
- **Adding `@st.cache_data` to the Run button handler:** `DeviceConfig` is unhashable (established
  Phase 39) — this will raise. Caching is manual via `session_state`, as in all other pages.

## Don't Hand-Roll

| Problem                                      | Don't Build                                                                | Use Instead                                                                                                                                            | Why                                                                                                                                                                                                                                                                                                                                                                                |
| -------------------------------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| CCE-vs-fluence physics                       | A new fluence-sweep loop in the page                                       | `petringa.run_radiation_damage()`                                                                                                                      | Already implemented, tested (`test_api_facades.py`), and the ONLY sanctioned way to call this physics per phase scope ("no new physics/core code")                                                                                                                                                                                                                                 |
| Dark-current decomposition physics           | A new SRH/TAT/SRV extraction loop                                          | `petringa.run_dark_current()`                                                                                                                          | Same — already implemented, returns `metadata["I_SRH"]`/`I_TAT`/`I_SRV"]` pre-decomposed                                                                                                                                                                                                                                                                                           |
| CSV serialization                            | A new ad-hoc CSV writer per page                                           | Extend `app/components/results.py::to_csv_bytes`                                                                                                       | Established Phase 39 convention: one dispatcher function, `sim_type`-keyed branches, commented `#` metadata header with software version + device config                                                                                                                                                                                                                           |
| Figure construction                          | Inline `go.Figure(...)` calls directly in the page `render()`              | New pure builders in `app/components/results.py`                                                                                                       | Established Phase 39/40 convention: figure-building is pure and unit-testable without a Streamlit runtime; page code only calls `st.plotly_chart(build_x_figure(result))`                                                                                                                                                                                                          |
| Log-scale multi-trace decomposition plotting | A naive `go.Scatter(..., y=result.metadata["I_TAT"])` on a log-axis figure | `np.abs()` + zero-value filtering per trace, mirroring `petringa.core.dark_current.plot_dark_current_decomposition` (matplotlib reference, line ~1205) | `I_TAT` is negative (generation, not recombination) and `I_SRV` is exactly 0.0 at default `S_n`/`S_p` in the empirically-verified run this session — a raw `semilogy`/log-y Plotly trace either drops negative points silently or raises on all-zero log(0); the existing core matplotlib plotter already solved this exact problem with `abs()` + `if np.any(I[mask] > 0)` guards |

**Key insight:** Both facades are complete, tested, physically-reasoned black boxes from this
phase's perspective. The entire engineering surface is: (1) turn simple widget scalars into the
facade's expected array/kwarg shape, (2) render whatever comes back — including NaN and
truncation — without crashing, (3) serialize it. Nothing here should touch `petringa/core/` or
`petringa/api/simulation.py`.

## Common Pitfalls

### Pitfall 1: "Dark current vs temperature" wording does not match `run_dark_current`'s x-axis

**What goes wrong:** FEAT-02 and Success Criterion #3 both say "dark current vs temperature",
but `run_dark_current(config, v_start, v_stop, n_points, ...)` sweeps **bias** — its `SimResult.x`
is `voltages` (V), not temperature (K). There is no `run_dark_current_vs_temperature` facade.
The closest existing temperature-sweep facade, `run_temperature_sweep`, sweeps **CCE**, not dark
current (`sim_type="temperature"`, `y` = CCE values, not current).
**Why it happens:** the ROADMAP/REQUIREMENTS wording was likely written from an early design
intent (`petringa.core.dark_current.dark_current_vs_fluence` — note NOT vs-temperature either —
exists in core, and there's a `nt_temperature_scale(T)` helper in `dark_current.py` used
internally, so temperature dependence of dark current IS modeled physically, just not exposed
as a swept axis in any existing facade).
**How to avoid:** This is Open Question #1 (below) — surface it explicitly to discuss-phase
before the planner locks the page's x-axis. Two readings are both defensible:

- **Reading A (recommended default if no user input arrives):** Treat "vs temperature" as
  loose/informal phrasing for "temperature-dependent dark current model" (which IS true —
  `nt_temperature_scale` makes `run_dark_current`'s N_t effective rate temperature-dependent
  via `config.T`), and build the page with `x=voltage` (matching the facade exactly), while
  exposing `config.T` as a normal sidebar-configured `DeviceConfig` field that the user can
  change to see the bias-sweep curve shift. This requires NO new facade code and satisfies
  the facade's actual return shape.
- **Reading B:** Literally sweep temperature as the x-axis by looping `run_dark_current` at
  a fixed bias across several `config.T` values in the page/app layer (NOT in `petringa/core`
  or `petringa/api`, keeping "no new physics" intact) — this is more expensive (N devsim
  solves instead of 1) and duplicates ParametricSweep's job outside of `ParametricSweep`
  (which Phase 42's batch-sweep page is explicitly supposed to own).
  **Warning signs:** If the planner writes a task like "sweep dark current across temperatures
  using run_dark_current's built-in T axis" — that facade has no such axis; verify against the
  actual signature in `petringa/api/simulation.py:665` before writing any task.

### Pitfall 2: Both facades silently produce partial/degenerate results at their own default parameters — try/except RuntimeError is not sufficient

**What goes wrong (empirically verified this session, live devsim run):**

- `run_dark_current(DeviceConfig())` with its own default `v_stop=-100.0` fails to converge at
  V≈-62.5V and the underlying `dark_current_sweep` returns **early with shorter arrays** — NOT
  a raised exception, and **no `"truncated"` key in metadata at all** (contrast with
  `run_cv`/`run_cce`, which DO set `metadata["truncated"]`/`metadata["requested_v_stop"]`).
  Confirmed: with `n_points=5` requested, only 3 points (`V = 0, -25, -50`) were returned.
- `run_radiation_damage(DeviceConfig())` with its own default `fluences=np.geomspace(1e13, 1e16, 6)`
  succeeds in returning a **full-length 6-element array**, but ONE interior element (fluence
  ≈3.98e13) is `np.nan` — because `cce_vs_fluence` (petringa/core/charge_collection.py, catch
  block around line 801-803) wraps each per-fluence-point solve in `try/except Exception`,
  logging a warning and setting `cce_values[i] = np.nan`, never propagating the exception to the
  caller.
  **Why it happens:** Both underlying core functions were designed to be robust/non-fatal for
  long parameter sweeps (reasonable for a batch/notebook context), silently degrading rather than
  aborting the whole sweep on one bad point. This is good behavior for a sweep, but it means the
  UI cannot assume "no exception raised" implies "full clean data".
  **How to avoid:**
  1. Keep the `try/except RuntimeError` wrapper around the `petringa.run_*` call (mirrors
     Phase 39 exactly) — it is still necessary for cases where the FIRST point fails outright
     (observed for `run_radiation_damage` in Phase 39's field/CCE precedent, and possible here too).
  2. Additionally render tolerantly: Plotly's `go.Scatter`/`go.Bar` silently skip/gap NaN y-values
     without crashing — no special code needed for the NaN case beyond not calling `.dropna()`
     accidentally.
  3. For `run_dark_current`'s silent truncation, there is no metadata flag to key off — the ONLY
     way to detect it at the page level is comparing `len(result.x)` against the requested
     `n_points` (the page already knows what it asked for, since it constructed the sweep). If
     `len(result.x) < n_points`, show an informational `st.info`/`st.warning` similar in spirit
     to Phase 39's "Sweep stopped at ... — solver reached full depletion" message, but phrased
     for dark current (no direct "requested_v_stop" value available in metadata, only the
     locally-known requested `v_stop` widget value).
  4. Choose WIDGET defaults (not facade defaults) that avoid tripping either failure mode on
     first Run: for dark current, `v_stop=-50.0` (not `-100.0`) converged fully in this
     session's live test up to at least the tested range; a full sweep to `-50.0` should be
     safer than `-100.0` (though not exhaustively re-verified point-by-point in this research
     session — recommend the plan's Wave 0/spike step re-confirm with `n_points` matching the
     final widget default). For radiation damage, the observed NaN occurred at
     fluence≈3.98e13 with `V_bias=-40.0` (the facade's own default) — a shallower `V_bias`
     (e.g. -20.0 or -30.0) is worth testing as a safer widget default, mirroring Phase 39's own
     finding that "-20.0 converges fine" for `run_field`/`run_cce` at the plain default config.
     **Warning signs:** A plan/task that says "the try/except RuntimeError pattern from Phase 39
     covers all error cases" — it does not, for these two facades specifically.

### Pitfall 3: Log-scale y-axis for dark-current decomposition needs `abs()` + zero-guarding

**What goes wrong:** `metadata["I_TAT"]` was observed negative in this session's live test
(effective generation, sign convention is net current not magnitude), and `metadata["I_SRV"]`
was exactly `0.0` at default `S_n`/`S_p` values. A naive Plotly `yaxis_type="log"` trace either
silently omits negative/zero points (Plotly's default behavior — no crash, but a
confusingly-vanishing trace) or produces a broken-looking chart if all points in a trace are
non-positive.
**Why it happens:** `I_TAT`'s physical sign (generation vs recombination) is not guaranteed
positive; `I_SRV` is genuinely zero unless non-default `S_n`/`S_p` values are supplied.
**How to avoid:** Mirror the existing core matplotlib reference plotter exactly —
`petringa.core.dark_current.plot_dark_current_decomposition` (line ~1205-1251) does:

```python
# Source: petringa/core/dark_current.py plot_dark_current_decomposition (existing reference)
for key, label, color, lw, ls in components:
    I = np.abs(np.asarray(sweep_result[key]))
    if np.any(I > 0):          # zero-guard: skip an all-zero trace entirely
        ax.semilogy(V, I, ...)
```

Translate this guard into the new Plotly builder: take `np.abs()` of each of
`I_total`/`I_SRH`/`I_TAT`/`I_SRV` before plotting, and only add a `go.Scatter` trace for a
component if `np.any(component > 0)` (skip `I_SRV` entirely from the legend if it's all-zero,
rather than showing an empty/broken log trace).
**Warning signs:** A chart that "looks empty" for I_TAT/I_SRV in a manual browser check — check
the raw metadata sign/zero-ness before assuming a code bug.

### Pitfall 4: Extending `to_csv_bytes` will invert an existing passing test

**What goes wrong:** `tests/test_app_csv_export.py::test_unknown_sim_type_raises_value_error`
currently constructs a `SimResult(sim_type="damage", ...)` specifically BECAUSE `"damage"` is
NOT YET a handled branch in `to_csv_bytes`, and asserts `pytest.raises(ValueError)`. The moment
Phase 41 adds a `"damage"` branch to `to_csv_bytes` (which it must, per FEAT-01 criterion #4
"CSV download button"), that existing test's assertion inverts — the call will no longer raise.
**Why it happens:** The test was written in Phase 39 using `"damage"` as a stand-in for "any
sim_type not yet implemented", not realizing it would become a real, implemented type one phase
later.
**How to avoid:** The plan MUST include editing `test_unknown_sim_type_raises_value_error` to
use a genuinely-unimplemented sim_type string (e.g. `"nonexistent_type"` or `"microdosimetry"`
— check Phase 42 timing; `"microdosimetry"` is also not yet a `to_csv_bytes` branch as of this
research, so it is a safe, still-future-proof placeholder, OR use an obviously-fake string like
`"not_a_real_sim_type"` to avoid any future collision). Flag this file edit explicitly in the
plan's file list — it is easy to omit since it looks like an unrelated pre-existing test.
**Warning signs:** `pytest -q tests/test_app_csv_export.py` going red after Phase 41's
`to_csv_bytes` extension, if this edit is missed.

### Pitfall 5: `run_radiation_damage`/`run_dark_current` 2D-config behavior is unverified

**What goes wrong:** Unlike `run_cv`/`run_cce` (which explicitly raise `NotImplementedError` for
`config.half_width_um is not None`), neither `run_radiation_damage` nor `run_dark_current` has
any 2D guard in their source (`petringa/api/simulation.py:447-767`, read in full this session —
no `if config.half_width_um is not None: raise` block exists in either function). Both delegate
to `create_dd_device`/`create_dark_current_device`, which are 1D constructors
(`petringa.core.drift_diffusion.create_dd_device`) — passing a 2D-configured `DeviceConfig`
would presumably either be silently ignored (if `half_width_um` isn't read by the 1D path at
all) or cause an unrelated downstream error. **Not tested empirically this session** (would
require constructing a 2D `DeviceConfig` and running a live devsim solve, out of scope for the
research budget).
**Why it happens:** These two facades were implemented in Phase 37 without an explicit 2D guard,
possibly because 2D radiation-damage/dark-current was never in scope and nobody anticipated a
2D `DeviceConfig` reaching them.
**How to avoid:** Add the same defensive 1D-only pre-check guard used in `cv.py`/`cce.py`
(`if cfg.half_width_um is not None: st.warning(...); st.stop()`) to BOTH new pages, even though
its necessity is unverified — the guard is cheap insurance against an unknown failure mode.
**Warning signs:** A 2D `DeviceConfig` passed to either facade producing a confusing, hard-to-
diagnose devsim error deep in `create_dd_device` rather than a clean page-level message.

## Code Examples

### Pure figure builder: CCE vs fluence (radiation damage)

```python
# New addition to app/components/results.py, mirrors build_cce_figure's shape
# (results.py:48-60) but with a log-x axis since fluence spans orders of magnitude
# (facade default np.geomspace(1e13, 1e16, 6) — log-spaced by construction).
def build_damage_figure(result: SimResult) -> go.Figure:
    """CCE vs proton fluence (log-x), NaN-tolerant (Pitfall 2)."""
    fig = go.Figure(
        data=go.Scatter(x=result.x, y=result.y, mode="lines+markers")
    )
    fig.add_hline(y=1.0, line_dash="dash", line_color="gray", opacity=0.5)
    fig.update_layout(
        title="CCE vs Proton Fluence",
        xaxis_title="Proton Fluence (p/cm²)",
        yaxis_title="Charge Collection Efficiency",
        xaxis_type="log",
        yaxis_range=[0, 1.1],
    )
    return fig
```

### Pure figure builder: dark current decomposition

```python
# New addition to app/components/results.py.
# Source pattern: petringa/core/dark_current.py::plot_dark_current_decomposition
# (matplotlib reference, ~line 1205) — abs() + zero-guard translated to Plotly.
def build_dark_current_figure(result: SimResult) -> go.Figure:
    """|I| vs bias, log-y, 4 overlaid traces (total + 3 components), zero-guarded."""
    fig = go.Figure()
    components = [
        ("I_total_placeholder", "Total", result.y),
        ("I_SRH", "SRH (bulk)", result.metadata["I_SRH"]),
        ("I_TAT", "TAT (effective)", result.metadata["I_TAT"]),
        ("I_SRV", "SRV (surface)", result.metadata["I_SRV"]),
    ]
    for _, label, values in components:
        I = np.abs(np.asarray(values))
        if np.any(I > 0):
            fig.add_trace(go.Scatter(x=result.x, y=I, mode="lines+markers", name=label))
    fig.update_layout(
        title="Dark Current Decomposition",
        xaxis_title="Voltage (V)",
        yaxis_title="|Dark Current| (A)",
        yaxis_type="log",
    )
    return fig
```

### CSV export: new `to_csv_bytes` branches

```python
# Extend app/components/results.py::to_csv_bytes with 2 new elif branches,
# following the exact existing dispatch pattern (results.py:97-141).
elif result.sim_type == "damage":
    df = pd.DataFrame({"fluence_p_per_cm2": result.x, "CCE": result.y})
    extra_header_lines = [
        f"# V_bias: {result.metadata['V_bias']}",
        f"# energy_MeV: {result.metadata['energy_MeV']}",
        "# WARNING: kappa (NIEL hardness factor) is a data-blocked placeholder; "
        "absolute Phi_crit numbers are unvalidated (see RESEARCH.md Pattern 3).",
    ]
elif result.sim_type == "dark_current":
    df = pd.DataFrame(
        {
            "bias_V": result.x,
            "I_total_A": result.y,
            "I_SRH_A": result.metadata["I_SRH"],
            "I_TAT_A": result.metadata["I_TAT"],
            "I_SRV_A": result.metadata["I_SRV"],
        }
    )
    extra_header_lines = [f"# area_cm2: {result.metadata['area_cm2']}"]
```

### AppTest mock fixtures (must include the two verified partial-failure shapes)

```python
# tests/test_app_radiation_damage_page.py — the NaN-in-array fixture is the REAL
# shape returned by run_radiation_damage at defaults (verified this session),
# not a hypothetical edge case. Testing only the happy path misses the actual bug class.
def _fake_run_radiation_damage_with_nan(cfg, **kwargs):
    return petringa.SimResult(
        config=cfg, sim_type="damage",
        x=np.array([1e13, 3.98e13, 1.58e14, 6.31e14, 2.51e15, 1e16]),
        y=np.array([0.996, np.nan, 0.985, 0.957, 0.860, 0.596]),  # NaN mid-array
        metadata={"V_bias": -40.0, "energy_MeV": 5.6},
    )

# tests/test_app_dark_current_page.py — the short-array (truncated, no metadata flag)
# fixture is the REAL shape returned by run_dark_current at its own default v_stop=-100.0.
def _fake_run_dark_current_truncated(cfg, **kwargs):
    return petringa.SimResult(
        config=cfg, sim_type="dark_current",
        x=np.array([0.0, -25.0, -50.0]),          # requested n_points=5, got 3
        y=np.array([4.9e-17, 2.4e-13, 3.5e-13]),
        metadata={
            "I_SRH": np.array([-1.8e-60, 3.4e-17, 3.4e-17]),
            "I_TAT": np.array([-3.8e-14, -2.5e-13, -3.5e-13]),  # negative — Pitfall 3
            "I_SRV": np.array([0.0, 0.0, 0.0]),                 # all-zero — Pitfall 3
            "area_cm2": 1e-4,
        },
    )
```

## State of the Art

| Old Approach                                                                                  | Current Approach                                                                                         | When Changed                                                                                               | Impact                                                                          |
| --------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| Jupyter notebooks calling `cce_vs_fluence`/`dark_current_sweep` directly (v2.0, Phases 14-15) | Streamlit UI calling `run_radiation_damage`/`run_dark_current` facades (v5.0, Phase 37 API, Phase 41 UI) | Phase 37 (2026-07-09) introduced the facade layer; Phase 41 is the first phase to expose it in the browser | Non-developers (Petringa group) can now run these sweeps without editing Python |

**Deprecated/outdated:** None specific to this phase — the underlying physics modules
(`petringa/core/radiation_damage.py`, `petringa/core/dark_current.py`) are unchanged from v2.0/v1.1
and explicitly out of scope to modify.

## Assumptions Log

| #   | Claim                                                                                                                                                                            | Section                   | Risk if Wrong                                                                                                                                                                                                                                                                                                 |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| A1  | "Dark current vs temperature" is loose phrasing and Reading A (bias x-axis + `config.T` as a normal sidebar field) satisfies FEAT-02's intent                                    | Pitfall 1, Open Questions | If the user actually wants a literal temperature x-axis (Reading B), the page architecture and CSV schema must change — an app-layer loop over `run_dark_current` at varying `config.T`, which is more compute-heavy and changes `SimResult.x` semantics                                                      |
| A2  | Widget default `v_stop=-50.0` for dark current avoids the partial-truncation failure mode observed at the facade's own `-100.0` default                                          | Pitfall 2, Pattern 2      | Not exhaustively re-verified in this session at `n_points` matching the final widget default — if `-50.0` also truncates, the plan's Wave 0 spike step should catch it and adjust the default further (e.g. `-40.0`), following Phase 39's precedent of iterating the default down until convergence is clean |
| A3  | Shallower `V_bias` (e.g. -20.0/-30.0) would avoid the NaN observed at radiation damage's default `V_bias=-40.0`                                                                  | Pitfall 2                 | Not tested empirically this session (only the -40.0 default was run); if untrue, the NaN-tolerant rendering (Plotly gap) is still an acceptable fallback per FEAT-01's criteria (curves are still displayed, one point missing is not a hard failure)                                                         |
| A4  | 2D `DeviceConfig` passed to `run_radiation_damage`/`run_dark_current` produces an unverified failure mode requiring a defensive pre-check guard                                  | Pitfall 5                 | If untrue (i.e. both facades handle 2D configs gracefully or 2D configs never reach these pages because the sidebar's 1D default is rarely changed), the added guard is harmless over-caution, not a bug                                                                                                      |
| A5  | `proton_energy_MeV` as a `selectbox` restricted to the 4 NIEL table keys {30, 62, 70, 150} is preferable UX to a free `number_input` matching the facade's raw default (5.6 MeV) | Pattern 2                 | Low risk either way — both are valid implementations; flagged as a discretion area, not a correctness issue                                                                                                                                                                                                   |

## Decision Addendum (post-research, pre-planning)

**Open Question #1 RESOLVED: Reading B (literal temperature x-axis) — CONFIRMED.**

Decision: the dark current page sweeps **temperature**, not bias, as its primary x-axis, per FEAT-02's
and Success Criterion #3's literal wording ("dark current vs temperature"), which appears identically
in both authoritative sources. Rationale: the plain reading, the doubled wording across ROADMAP and
REQUIREMENTS, the domain norm (a detector-physics "dark current vs temperature" plot is a specific,
standard characterization plot with T on the x-axis), and Phase 43's integration audit (which will
check FEAT-02 against its literal wording) all point to B over the effort-based "zero new code"
tiebreaker that originally favored Reading A.

**Implementation path — use `petringa.ParametricSweep`, NOT a hand-rolled loop:**

`petringa/api/sweep.py::ParametricSweep` already exists and is exactly fit for this:
`ParametricSweep(base_config=cfg, param="T", values=temperatures, sim_fn=petringa.run_dark_current,
sim_kwargs={"v_start": ..., "v_stop": fixed_bias, "n_points": 2}).run()` — clones `cfg` per
temperature via `dataclasses.replace`, calls `run_dark_current` per clone, returns
`list[SimResult]` of length `len(temperatures)`. This keeps "no new physics/core code" fully intact
(ParametricSweep and run_dark_current are both already implemented and tested) and reuses tested
sweep code instead of a page-local loop.

**New page architecture (supersedes Pattern 1/2's dark-current section and the Code Examples above):**

- Widget inputs: `T_min`, `T_max`, `n_temperatures` (construct `np.linspace(T_min, T_max, n)`,
  mirroring the fluence-range widget-group pattern from Pattern 2), plus a single **fixed operating
  bias** input (e.g. `V_bias` number_input, default a shallow/safe value per Pitfall 2 — NOT
  `run_dark_current`'s own `v_stop=-100.0` default) at which dark current is evaluated for each T.
  Optional advanced expander for `N_t`/`S_n`/`S_p` overrides (unchanged from original Pattern 2).
- Run button calls `ParametricSweep(base_config=cfg, param="T", values=temperatures,
sim_fn=petringa.run_dark_current, sim_kwargs={"v_start": V_bias, "v_stop": V_bias, "n_points": 1}
or similar single-point-per-T kwargs).run()` — each per-temperature `run_dark_current` call should
  request a minimal bias sweep (e.g. `n_points=1` or `2`, `v_start=v_stop=V_bias`) since only the
  fixed-bias operating point is needed per temperature; verify `run_dark_current` tolerates
  `v_start == v_stop` / `n_points=1` without raising — flag this as a Wave 0 spike check if unverified.
- Result shape: `ParametricSweep.run()` returns a **list of `SimResult`** (one per temperature),
  not a single `SimResult`. The dark current **page** (not the builder) aggregates this list into
  a single `SimResult` — extracting one point (`result.y[0]` or an aggregate) per temperature from
  each per-T `SimResult` to build `x=temperatures`, `y`/`metadata["I_SRH"/"I_TAT"/"I_SRV"]` arrays
  of length `n_temperatures` — BEFORE calling `build_dark_current_figure`. This keeps
  `build_dark_current_figure` a pure single-`SimResult`-in function, consistent with every other
  builder in `results.py`, and matches how `to_csv_bytes` is called elsewhere (one `SimResult` in,
  one CSV out; columns: `T_K`, `I_total_A`, `I_SRH_A`, `I_TAT_A`, `I_SRV_A`, one row per
  temperature — NOT one row per bias point). Page-level caching stores the aggregated
  `SimResult` (`session_state["dark_current_result"]`), not the raw `list[SimResult]`.
- Pitfalls 2 (silent truncation/NaN) and 3 (log-scale abs+zero-guard) still apply per-temperature
  per-call — each per-T `run_dark_current` call can independently truncate or NaN; the page/builder
  must tolerate a per-temperature result being partial or missing a component.
- Pitfall 5 (2D-config guard) is unaffected — still add the defensive 1D-only pre-check.
- Pitfall 4 (`to_csv_bytes` regression) is unaffected — still applies, but the `"dark_current"`
  CSV branch's column schema changes to the per-temperature shape described above.

This addendum supersedes the dark-current-specific portions of Pattern 1 (skeleton is still valid
structurally — guard → widgets → Run → cache → render → download — but the Run step now invokes
`ParametricSweep` instead of a single `run_dark_current` call), Pattern 2 (dark current widget group,
now T-range + fixed bias instead of bias-range), the dark-current Code non-workingExample (now needs a
list-of-SimResult-aware builder), and Open Question #1/Assumption A1 (RESOLVED, no longer open).
The radiation damage page (FEAT-01) is UNAFFECTED by this addendum — it remains a single
`run_radiation_damage` call as originally researched.

## Open Questions (RESOLVED)

All three questions below are resolved: Q1 by the Decision Addendum (temperature x-axis via
`ParametricSweep`, confirmed with the advisor), Q2 by 41-01-PLAN.md's Task 1 live-devsim spike
(confirms convergent widget defaults before Wave 2 builds around them), and Q3 by 41-UI-SPEC.md's
Discretion Resolved section (`proton_energy_MeV` as a `st.selectbox` restricted to the 4 real NIEL
table keys). Retained below for provenance/history only — no further action needed.

1. **Does "dark current vs temperature" (FEAT-02, Success Criterion #3) mean a literal
   temperature x-axis, or is it loose phrasing for "the existing bias-sweep facade, whose
   underlying generation-rate model IS temperature-dependent via config.T"?**
   - What we know: `run_dark_current(config, v_start, v_stop, n_points, ...)` sweeps **bias**
     (`SimResult.x = voltages`). There is no facade that sweeps temperature and returns dark
     current. `run_temperature_sweep` exists but sweeps **CCE**, not dark current. The core
     module `petringa/core/dark_current.py` does have `nt_temperature_scale(T)`, making the
     effective generation rate `N_t` temperature-dependent — so `config.T` (already a normal
     `DeviceConfig` field, already exposed in the sidebar) genuinely does affect
     `run_dark_current`'s output, just not as the swept axis.
   - What's unclear: whether the ROADMAP/REQUIREMENTS wording is a literal spec requirement or
     informal shorthand written before the exact facade signature was finalized in Phase 37.
   - Recommendation: default to Reading A (x=bias voltage, `config.T` remains a normal sidebar
     field the user can change to see the curve shift) since it requires zero new
     facade/app-layer sweep logic and matches the facade's actual, tested return shape. Surface
     this explicitly to `/gsd:discuss-phase` before planning locks in the page's x-axis and CSV
     schema — if the user confirms Reading B is required, the plan needs an app-layer
     `for T in temperatures: run_dark_current(replace(cfg, T=T), ...)` loop (using
     `dataclasses.replace`, matching `ParametricSweep`'s own internal pattern per STATE.md) at a
     fixed bias, which is a materially different (and slower — N devsim solves) page design.

2. **Does the widget default `v_stop=-50.0` for dark current (recommended in Pitfall 2/Pattern 2)
   actually converge cleanly at the FINAL chosen `n_points` widget default, or does it also
   partially truncate?**
   - What we know: `-100.0` truncates around -62.5V (empirically verified, `n_points=5`
     requested → 3 points returned). `-50.0` was NOT independently re-verified at the sweep
     granularity a real widget default would use.
   - What's unclear: the exact convergence boundary for `dark_current_sweep` at the plain
     default `DeviceConfig()`.
   - Recommendation: the plan should include a cheap Wave 0 spike (single live `run_dark_current`
     call at the chosen widget defaults, mirroring Phase 39-01's spike precedent) to confirm the
     default converges cleanly before writing the full page, OR accept that the truncation-
     detection UI (Pitfall 2, item 3) makes even a partially-truncating default acceptable per
     FEAT-02's criteria (the decomposition still renders for whatever prefix DID converge).

3. **Should the radiation damage page's `proton_energy_MeV` widget be a `selectbox` restricted
   to the NIEL table's 4 real keys, or a free `number_input` matching the facade's raw
   5.6 MeV default?**
   - What we know: `get_hardness_factor` linearly interpolates via `np.interp` for any energy
     not in `{30, 62, 70, 150}`; the facade's own default (5.6 MeV) is BELOW the lowest table key.
   - What's unclear: whether the Petringa group has a specific experimental beam energy they
     always use (which might not be one of the 4 table keys either).
   - Recommendation: low-stakes; either widget type is fully functional. Flagged in the
     Assumptions Log (A5) as a discretion area for the planner, not a blocker.

## Environment Availability

Skipped — this phase has no external dependencies beyond what Phases 35-40 already installed
and verified (`streamlit`, `plotly`, `pandas`, `numpy`, `devsim`, all already present and used).
No new CLI tools, services, or runtimes are introduced.

## Validation Architecture

`workflow.nyquist_validation` is absent from `.planning/config.json` (checked directly this
session — the file has `mode`, `parallelization`, `commit_docs`, `model_profile`, `workflow.research`,
`workflow.plan_check`, `workflow.verifier`, `workflow._auto_chain_active`, `granularity`, but no
`nyquist_validation` key). Per the absent-key convention, treat as enabled and include this section.

### Test Framework

| Property           | Value                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Framework          | pytest 7.0+ (per `pyproject.toml` `pytest>=7.0`), `streamlit.testing.v1.AppTest` for page tests                                                                                                                                                                                                                                                                                                                                            |
| Config file        | `pytest.ini` (repo root; declares only the `slow` marker)                                                                                                                                                                                                                                                                                                                                                                                  |
| Quick run command  | `uv run pytest tests/test_app_radiation_damage_page.py tests/test_app_dark_current_page.py tests/test_app_csv_export.py -q`                                                                                                                                                                                                                                                                                                                |
| Full suite command | Per-file/per-class isolation is the durable convention (STATE.md Phase 35 finding — bare single-process `pytest -q` is unsatisfiable due to devsim resource exhaustion, confirmed pre-existing not caused by any phase). Run each new/modified test file individually: `uv run pytest tests/test_app_radiation_damage_page.py -q && uv run pytest tests/test_app_dark_current_page.py -q && uv run pytest tests/test_app_csv_export.py -q` |

### Phase Requirements → Test Map

| Req ID                                              | Behavior                                                                                             | Test Type                                            | Automated Command                                                                                                                                                                                              | File Exists?                      |
| --------------------------------------------------- | ---------------------------------------------------------------------------------------------------- | ---------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------- |
| FEAT-01 (criterion 1: run + CCE-vs-fluence plot)    | Radiation damage page Run button caches a `SimResult` and renders `build_damage_figure`              | AppTest (mocked `petringa.run_radiation_damage`)     | `uv run pytest tests/test_app_radiation_damage_page.py::test_run_caches_damage_result -x`                                                                                                                      | ❌ Wave 0 (new file)              |
| FEAT-01 (criterion 2: persistent kappa banner)      | `st.warning` with kappa-data-blocked text appears on every page load, even with no device configured | AppTest                                              | `uv run pytest tests/test_app_radiation_damage_page.py::test_kappa_banner_persistent -x`                                                                                                                       | ❌ Wave 0 (new file)              |
| FEAT-01 (NaN-tolerant rendering — Pitfall 2)        | Page does not crash when `run_radiation_damage` returns a NaN mid-array                              | AppTest (mocked with the verified NaN fixture)       | `uv run pytest tests/test_app_radiation_damage_page.py::test_nan_in_result_does_not_crash -x`                                                                                                                  | ❌ Wave 0 (new file)              |
| FEAT-02 (criterion 3: decomposition overlay)        | Dark current page renders 4 overlaid traces (total + SRH + TAT + SRV), zero/negative-guarded         | AppTest + pure builder unit test                     | `uv run pytest tests/test_app_dark_current_page.py::test_run_caches_dark_current_result -x` and `uv run pytest tests/test_app_results_builders.py::test_build_dark_current_figure_guards_zero_and_negative -x` | ❌ Wave 0 (both new files)        |
| FEAT-02 (truncation-tolerant rendering — Pitfall 2) | Page does not crash and informs the user when `run_dark_current` returns fewer points than requested | AppTest (mocked with the verified truncated fixture) | `uv run pytest tests/test_app_dark_current_page.py::test_truncated_result_shows_info -x`                                                                                                                       | ❌ Wave 0 (new file)              |
| FEAT-01 + FEAT-02 (criterion 4: CSV download)       | `to_csv_bytes` handles `"damage"` and `"dark_current"` sim_types with correct columns                | Pure unit test                                       | `uv run pytest tests/test_app_csv_export.py::test_damage_csv_columns_and_header tests/test_app_csv_export.py::test_dark_current_csv_columns_and_header -x`                                                     | ❌ Wave 0 (extend existing file)  |
| (regression guard)                                  | `to_csv_bytes`'s unknown-type test still passes after adding the "damage" branch                     | Pure unit test                                       | `uv run pytest tests/test_app_csv_export.py::test_unknown_sim_type_raises_value_error -x`                                                                                                                      | ✅ exists, needs edit (Pitfall 4) |

### Sampling Rate

- **Per task commit:** run the specific new/modified test file in isolation (per-file
  convention established Phase 35/39/40)
- **Per wave merge:** run all of `test_app_radiation_damage_page.py`,
  `test_app_dark_current_page.py`, `test_app_csv_export.py` sequentially (not concatenated —
  isolation convention)
- **Phase gate:** all listed files green individually before `/gsd:verify-work`; a manual browser
  smoke-check of both pages (mirroring the Phase 39/40 "browser sign-off" precedent recorded in
  STATE.md) is recommended given this phase's two confirmed silent-partial-failure modes —
  automated AppTest with mocked facades cannot prove the REAL facade's live behavior matches the
  fixtures used in tests; a real `streamlit run app/main.py` click-through on both new pages
  should be the final verification step, exactly as Phase 39/40 did.

### Wave 0 Gaps

- [ ] `tests/test_app_radiation_damage_page.py` — new file, covers FEAT-01
- [ ] `tests/test_app_dark_current_page.py` — new file, covers FEAT-02
- [ ] `tests/test_app_csv_export.py` — extend with 2 new test functions + fix
      `test_unknown_sim_type_raises_value_error` (Pitfall 4)
- [ ] Optional: `tests/test_app_results_builders.py` — pure unit tests for
      `build_damage_figure`/`build_dark_current_figure` isolated from AppTest (mirrors
      `tests/test_app_geometry_viewer.py`'s pure-builder-testing precedent from Phase 40); can
      also be folded into the existing `test_app_csv_export.py` file or a new
      `test_app_results.py` if the planner prefers consolidating all `results.py` unit tests
- No framework install needed — `pytest`, `streamlit.testing.v1.AppTest` already present and
  used by every prior Phase 38-40 test file.

## Security Domain

`security_enforcement` is absent from `.planning/config.json` — treat as enabled per the
absent-key convention.

### Applicable ASVS Categories

| ASVS Category         | Applies | Standard Control                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| --------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| V2 Authentication     | no      | Local/shared-lab-server Streamlit app, no auth layer exists anywhere in this project (established Phase 38, out of scope)                                                                                                                                                                                                                                                                                                                                     |
| V3 Session Management | no      | `st.session_state` is Streamlit's own in-process mechanism, not a security boundary here                                                                                                                                                                                                                                                                                                                                                                      |
| V4 Access Control     | no      | Single-user local tool, no roles/permissions concept in this codebase                                                                                                                                                                                                                                                                                                                                                                                         |
| V5 Input Validation   | yes     | Numeric widget inputs (`fluence_min`/`max`, `n_points`, `V_bias`, `v_start`/`v_stop`, optional `N_t`/`S_n`/`S_p`) should use Streamlit's built-in `number_input` `min_value`/type coercion rather than parsing free text — Streamlit's own widgets already provide this; no manual `eval`/`exec` or string-to-array parsing is needed anywhere in this phase (unlike Phase 42's CSV-upload page, which DOES need file-content validation — out of scope here) |
| V6 Cryptography       | no      | No cryptographic operations in this phase                                                                                                                                                                                                                                                                                                                                                                                                                     |

### Known Threat Patterns for this stack

| Pattern                                                                                                                                                         | STRIDE                    | Standard Mitigation                                                                                                                                                                                                                      |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Malformed numeric widget input causing an unhandled exception deep in devsim (e.g. `n_points=0` causing an empty/degenerate `np.geomspace`/`np.linspace` array) | Denial of Service (local) | Use `min_value=2` on all `n_points`-style `number_input` widgets (Streamlit enforces this client-side); the existing `try/except RuntimeError` + this research's NaN/truncation tolerance already covers the solver-side failure surface |
| No new file-upload or eval surface introduced in this phase                                                                                                     | —                         | Not applicable — this phase has no CSV upload (that's Phase 42's microdosimetry page) or any other user-supplied-code-path surface                                                                                                       |

## Sources

### Primary (HIGH confidence — verified this session via direct code read and/or live execution)

- `petringa/api/simulation.py` (full file read) — `run_radiation_damage` (line 447),
  `run_dark_current` (line 665) exact signatures, defaults, docstrings, warnings
- `petringa/core/radiation_damage.py` (full file read) — `NIEL_HARDNESS_PROTON_SIC` table,
  `get_hardness_factor`, `compute_phi_crit`, kappa-data-blocked provenance comments
- `petringa/core/dark_current.py` (full file read) — `dark_current_sweep`,
  `extract_dark_current_components`, `plot_dark_current_decomposition` (matplotlib reference),
  `nt_temperature_scale`
- `petringa/core/charge_collection.py` — grep-verified the `except Exception ... cce_values[i] =
np.nan` pattern (lines 801-803) underlying `cce_vs_fluence`'s silent-NaN behavior
- Live execution (this session): `petringa.run_dark_current(DeviceConfig(), n_points=5)` —
  confirmed truncation to 3 points at V≈-62.5V, no `"truncated"` metadata key present
- Live execution (this session): `petringa.run_radiation_damage(DeviceConfig())` — confirmed
  full-length array with NaN at fluence≈3.98e13, no exception raised
- `app/workflows/cv.py`, `cce.py`, `field_map.py` (full files read) — Run→cache→render→download
  skeleton, 1D-only pre-check pattern, try/except RuntimeError pattern
- `app/components/results.py` (full file read) — existing pure builders, `to_csv_bytes` dispatch
  structure, existing `sim_type` branches (`"cv"`, `"cce"`, `"field"`)
- `app/main.py` (full file read) — confirmed both pages are already registered in navigation
- `app/workflows/radiation_damage.py`, `app/workflows/dark_current.py` (full files read) —
  confirmed both are currently pure placeholders
- `tests/test_app_run_mockability.py`, `tests/test_app_field_page.py`,
  `tests/test_app_csv_export.py` (full files read) — mock-seam pattern, AppTest assertion
  surface, existing `test_unknown_sim_type_raises_value_error` (Pitfall 4 source)
- `tests/test_api_facades.py` (full file read) — confirmed `run_radiation_damage`/
  `run_dark_current` have ONLY contract-level (non-executing) tests in the existing suite; no
  prior end-to-end test exercised either facade against a live devsim solve
- `.planning/STATE.md`, `.planning/REQUIREMENTS.md`, `.planning/ROADMAP.md` — decisions,
  blockers, Phase 39 RuntimeError precedent, FEAT-01/FEAT-02 exact wording
- `.planning/phases/39-.../39-PATTERNS.md`, `.planning/phases/40-geometry-viewer/40-PATTERNS.md`
  — pattern-map precedent for pure-builder/page structure
- `app/components/device_sidebar.py` (full file read) — confirmed `config.T` is already a
  normal sidebar-exposed `DeviceConfig` field (relevant to Open Question #1 Reading A)

### Secondary (MEDIUM confidence)

None — no external web sources were needed for this phase; all facts were derivable from direct
code inspection and live execution of the already-implemented, in-repo facades.

### Tertiary (LOW confidence)

None.

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH — no new dependencies, all already declared and used
- Architecture: HIGH — direct extension of an already-established, 3-times-repeated Phase 39/40
  page pattern
- Pitfalls: HIGH — the two most important pitfalls (silent truncation, silent NaN) were verified
  by actually executing both facades against a live devsim solve this session, not inferred from
  reading code alone
- Open Question #1 (temperature wording): MEDIUM — the factual claim (no temperature-sweep dark
  current facade exists) is HIGH confidence (verified by reading the full simulation.py file);
  the interpretation of what the requirement author intended is inherently uncertain and flagged
  for user confirmation

**Research date:** 2026-07-13
**Valid until:** No expiry driver — this research is tied to the current, frozen state of
`petringa/api/simulation.py` and `petringa/core/{radiation_damage,dark_current}.py`, which are
explicitly out of scope to modify in this phase or any v5.0 phase. Re-research only if those
facades change (e.g. a future v6.0 phase adds real SR-NIEL data or a temperature-sweep dark
current facade).
