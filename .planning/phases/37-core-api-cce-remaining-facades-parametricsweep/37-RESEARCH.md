# Phase 37: Core API — CCE + Remaining Facades + ParametricSweep - Research

**Researched:** 2026-07-08
**Domain:** Python library facade/adapter layer over an existing devsim TCAD simulation core (refactor-only, no physics changes)
**Confidence:** HIGH (all findings verified against the in-repo codebase and design spec)

## Summary

Phase 37 completes the `petringa` public API by adding seven simulation facades (`run_cce`, `run_radiation_damage`, `run_dark_current`, `run_temperature_sweep`, `run_flash_recombination`, `run_transient`, `run_microdosimetry`) plus a `ParametricSweep` utility class. All underlying physics already exists in `petringa/core/*` and was written in v3.0/v4.0 — this phase wraps it, it does NOT reimplement it (STATE.md: "No physics changes allowed in any v5.0 phase — refactor only").

**The single most important finding:** Phase 37 CANNOT clone the Phase 36 `build_device(config) → core_fn(device_info)` pattern for most facades. The Phase 37 core functions split into three structurally incompatible shapes. Five of them (`cce_vs_bias`, `sweep_cce_vs_temperature`, `cce_vs_dose_rate`, `cce_vs_fluence`, `transient_cce_vs_dose_rate`) build their _own_ devsim device internally from scalar kwargs and do **not** accept a pre-built `device_info`. Refactoring them to accept `device_info` would break the 20 validated notebooks that call `petringa.core.*` directly and would itself constitute a structural change. The correct Phase 37 pattern is a **`DeviceConfig` → core-kwargs adapter**, not the Phase 36 device-passing pattern.

**Primary recommendation:** Classify each of the 7 facades into one of three buckets (see Architecture Patterns) and wrap the core function at its _existing_ entry point. `run_cce` is the flagship (only facade with a real physics gate: CCE ∈ [0,1]) — give it a slow integration test like `run_cv` got. The other 6 have a weaker acceptance bar (imports + accepts `DeviceConfig` first arg) — give them light contract tests. `ParametricSweep` should be tested with a **fake `sim_fn`** so its unit tests never touch devsim.

## Architectural Responsibility Map

| Capability                | Primary Tier            | Secondary Tier                                                                    | Rationale                                                                |
| ------------------------- | ----------------------- | --------------------------------------------------------------------------------- | ------------------------------------------------------------------------ |
| `run_cce` (CCE vs bias)   | API / Backend (facade)  | Core physics (`charge_collection.cce_vs_bias`)                                    | Facade adapts `DeviceConfig`→kwargs; core owns the devsim solve          |
| `run_radiation_damage`    | API / Backend (facade)  | Core (`charge_collection.cce_vs_fluence`)                                         | Facade maps fluence array + energy; core builds fresh device per fluence |
| `run_dark_current`        | API / Backend (facade)  | Core (`dark_current.create_dark_current_device` + `dark_current_sweep`)           | Facade must set up TAT+SRV models before sweeping                        |
| `run_temperature_sweep`   | API / Backend (facade)  | Core (`temperature_sweep.sweep_cce_vs_temperature`)                               | Facade maps temp array; core builds a device per T                       |
| `run_flash_recombination` | API / Backend (facade)  | Core (`flash_recombination.cce_vs_dose_rate`)                                     | Facade maps dose-rate array; core builds device internally               |
| `run_transient`           | API / Backend (facade)  | Core (`transient.TransientSolver` / `transient_cce_vs_dose_rate`)                 | Facade orchestrates solver init + pulse; underspecified in design spec   |
| `run_microdosimetry`      | API / Backend (facade)  | Core (`mc_coupling.load_mc_events_csv` + `microdosimetry.lineal_energy_spectrum`) | Pure data pipeline — NO devsim device built                              |
| `ParametricSweep`         | API / Backend (utility) | —                                                                                 | Pure Python orchestration over any `sim_fn`; devsim-agnostic             |

## Standard Stack

No new external packages are introduced by this phase. All work is internal Python over the already-declared dependency set (`numpy`, `scipy`, `pandas`, `devsim`). [VERIFIED: petringa/pyproject.toml + codebase]

### Core (already installed, already declared)

| Library              | Version     | Purpose                                                                                             | Why Standard                                                  |
| -------------------- | ----------- | --------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| numpy                | (installed) | Arrays, `np.linspace`, `np.geomspace` for sweep grids                                               | Already the array type of every `SimResult.x`/`.y`            |
| pandas               | >=2.0       | Several core sweep fns return `pd.DataFrame` (temperature, radiation hardness, transient dose-rate) | Facade must convert DataFrame → `SimResult.x`/`.y`/`metadata` |
| scipy                | (installed) | Used inside `run_field` (LinearNDInterpolator); no new use expected here                            | Already a runtime dep                                         |
| devsim               | (installed) | The C-extension TCAD solver wrapped by core physics                                                 | Process-global state — drives the lifecycle pitfall below     |
| dataclasses (stdlib) | —           | `ParametricSweep` uses `dataclasses.replace(base_config, **{param: value})` to clone+inject         | Idiomatic config cloning                                      |

### Supporting (stdlib only)

| Library            | Purpose                                                                           | When to Use                                                  |
| ------------------ | --------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| `uuid` (stdlib)    | Unique device names (only relevant for buckets b/c that build via `build_device`) | Buckets (a) already generate their own uuid names internally |
| `logging` (stdlib) | Facade-level warning on cleanup failure (Phase 36 convention)                     | In any `finally` cleanup block the facade adds               |

### Alternatives Considered

| Instead of                                 | Could Use                                                                                    | Tradeoff                                                                                      |
| ------------------------------------------ | -------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| Wrapping `cce_vs_bias` directly            | Hand-rolling via `add_generation_to_dd` + `compute_cce_from_dd` on a `build_device()` device | REJECTED — reconstructs physics orchestration = risk of physics change; forbidden by STATE.md |
| Refactoring core fns to take `device_info` | Adapter mapping `DeviceConfig`→core kwargs                                                   | REJECTED — refactoring core signatures breaks 20 notebooks + is a structural change           |

**Installation:** None required — all dependencies already declared and installed. [VERIFIED: pyproject.toml]

## Package Legitimacy Audit

> Not applicable — this phase installs **zero** external packages. All code is internal Python over already-declared, already-installed dependencies (`numpy`, `scipy`, `pandas`, `devsim`, stdlib). slopcheck / registry verification is moot because no `pip install` / `npm install` step exists in this phase.

**Packages removed due to slopcheck [SLOP] verdict:** none (no packages installed)
**Packages flagged as suspicious [SUS]:** none (no packages installed)

## Architecture Patterns

### The Phase 36 template (what run_cv / run_field established)

Located in `petringa/api/simulation.py`. The established shape is:

```python
def run_cv(config: DeviceConfig, ...) -> SimResult:
    if config.half_width_um is not None:
        raise NotImplementedError(...)          # 2D guard, mirrors design
    reset_devsim_fully()                        # clean global devsim state at entry
    device_info = build_device(config)          # DeviceConfig -> solved DD device
    try:
        result = cv_sweep(device_info, ...)     # call core physics, unchanged
        return SimResult(config=config, sim_type="cv", x=..., y=..., metadata={...}, mesh=None)
    finally:
        try:
            devsim.delete_device(device=device_info["device_name"])
        except Exception:
            logger.warning(..., exc_info=True)
            reset_devsim_fully()
```

`SimResult` and `MeshData` are frozen dataclasses in `petringa/api/results.py`:

```python
@dataclass
class SimResult:
    config: DeviceConfig
    sim_type: str      # "cv" | "cce" | "field" | "damage" | ...
    x: np.ndarray      # primary axis (bias, depth, fluence, ...)
    y: np.ndarray      # primary output
    metadata: dict = field(default_factory=dict)
    mesh: MeshData | None = None
```

[VERIFIED: petringa/api/simulation.py, petringa/api/results.py]

### The three-bucket classification (backbone of this phase)

Each Phase 37 core function falls into exactly one bucket. **This determines whether `build_device()` is used at all.**

#### Bucket (a): Core builds its own device from scalar kwargs — `build_device()` NOT used

The facade is a **`DeviceConfig` → core-kwargs adapter**. It maps config fields to the core fn's parameters (many accept a `device_kwargs=` dict or explicit `epi_thickness_cm=`, `N_D_junction=`, etc.), calls the core fn, and repackages the returned dict/DataFrame as a `SimResult`. The facade does **not** call `build_device()` and does **not** call `reset_devsim_fully()`/`delete_device` around the core call **if** the core fn already cleans up its own device (verified true for `cce_vs_bias` — see Pitfall 1).

Facades in bucket (a):

- **`run_cce`** → `charge_collection.cce_vs_bias(V_range, epi_thickness_cm, alpha_range_cm, generation_rate, device_kwargs)`. Returns dict with `"voltages"`, `"cce_values"`, `"I_collected"`, `"I_generated"`. 1D-only (uses `create_dd_device`) → mirror `run_cv`'s `NotImplementedError` 2D guard. **`cce_vs_bias` deletes its own device in a `finally` block (verified line 54-56), so the facade must NOT add a second `delete_device`.** [VERIFIED: charge_collection.py:432-597]
- **`run_radiation_damage`** → `charge_collection.cce_vs_fluence(fluence_range, V_bias, epi_thickness_cm, N_D_junction, N_D_bulk, L_transition, energy_MeV, ...)`. Returns dict with `"fluences"`, `"cce_values"` (creates a fresh device per fluence). Design spec §3.4: `run_radiation_damage(cfg, fluences=[...], proton_energy_MeV=5.6)` → `x=fluence, y=CCE_post_damage`. [VERIFIED: charge_collection.py:597, spec:156]
- **`run_temperature_sweep`** → `temperature_sweep.sweep_cce_vs_temperature(temperatures, voltages, method, **device_kwargs)`. Returns a **`pd.DataFrame`** (long-format columns `T, V, CCE`) — facade must extract arrays. [VERIFIED: temperature_sweep.py:215]
- **`run_flash_recombination`** → `flash_recombination.cce_vs_dose_rate(dose_rates_Gy_s, V_bias, epi_thickness_cm, E_MeV, ...)`. Builds device internally. NOTE: module carries an explicit SCOPE/LIMITATIONS block — outputs are a numerical sensitivity bound, not validated FLASH physics. [VERIFIED: flash_recombination.py:263]

#### Bucket (b): Core takes `device_info` but needs extra model setup — `build_device()` usable but insufficient alone

The facade builds a device (via `create_dark_current_device` / `build_device`), runs the required extra setup, then calls the sweep/solver. Wrap with the full Phase 36 `reset → try → finally delete` lifecycle.

Facades in bucket (b):

- **`run_dark_current`** → recommended entry is `dark_current.create_dark_current_device(T, N_t, S_n, S_p, **kwargs)` which internally does `create_dd_device` + `setup_tat_model` + `setup_surface_recombination`, THEN `dark_current.dark_current_sweep(device_info, V_range, area, V_step)`. Returns dict with `voltages, I_total, I_SRH, I_TAT, I_SRV, J_*`. Do NOT use bare `build_device()` here — it omits the TAT/SRV model setup the sweep requires. [VERIFIED: dark_current.py:435, 538]
- **`run_transient`** → `transient.TransientSolver(device_info, contact, method)` needs `.initialize()` then `.simulate_pulse(G_spatial, ...)`. A ready-made wrapper `transient.transient_cce_vs_dose_rate(V_bias, dose_rates, ...)` builds its own device (bucket-a shaped) and may be the simplest satisfier of the weak criterion-2 bar. **`run_transient` has NO signature in design spec §3.4** — it is underspecified; keep minimal. [VERIFIED: transient.py:172-545, spec §3.4 omits it]

#### Bucket (c): No devsim device at all — pure data pipeline

- **`run_microdosimetry`** → `mc_coupling.load_mc_events_csv(filepath, column_map, pos_unit, energy_unit)` → energies → `microdosimetry.mean_chord_length(sv_thickness_um, sv_width_um, sv_depth_um)` → `microdosimetry.lineal_energy_spectrum(collected_energies_keV, l_bar_um, ...)` which returns a dict with `bin_centers`, `f_y`, `d_y`, `y_F`, `y_D`. Design spec §3.4: `run_microdosimetry(cfg, mc_csv_path, sv_thickness_um=10, sv_width_um=150)` → `x=lineal_energy, y=yd_y`. No devsim, no `reset_devsim_fully`. [VERIFIED: mc_coupling.py:107, microdosimetry.py:36,117, spec:160]

### ParametricSweep design (LIB-07)

Design spec §3.5 fully specifies the API — implement exactly this:

```python
sweep = ParametricSweep(
    base_config=cfg,
    param="epi_thickness_um",       # a DeviceConfig field name (str)
    values=[5, 10, 20, 50],
    sim_fn=run_cce,                 # any facade
    sim_kwargs={"v_start": -10, "v_stop": -200, "n_points": 20},
)
results = sweep.run()               # -> list[SimResult], len == len(values)
```

Recommended `.run()` idiom: for each `value`, build a cloned config via `dataclasses.replace(base_config, **{param: value})`, then call `sim_fn(cfg_i, **sim_kwargs)`, collect into a list. [VERIFIED: spec §3.5:167-177]

### System Architecture Diagram

```
User code / Streamlit page
        │
        │  DeviceConfig(...)  +  facade-specific args
        ▼
┌─────────────────────────────────────────────────────────┐
│  petringa/api/simulation.py  (facade layer — Phase 37)   │
│                                                          │
│  run_cce ─────────► adapter (bucket a) ──► cce_vs_bias   │
│  run_radiation_damage ► adapter (a) ────► cce_vs_fluence │
│  run_temperature_sweep ► adapter (a) ──► sweep_cce_vs_T  │
│  run_flash_recombination ► adapter (a) ► cce_vs_dose_rate│
│  run_dark_current ► build+setup (b) ──► dark_current_sweep│
│  run_transient ► build+solver (b) ────► TransientSolver  │
│  run_microdosimetry ► data pipe (c) ──► load_mc_events + │
│                                          lineal_energy    │
└───────────────┬──────────────────────────┬──────────────┘
                │ (a,b) devsim solve         │ (c) numpy/pandas only
                ▼                            ▼
        petringa/core/*.py            data/synthetic_mc_events.csv
        (UNCHANGED physics)           (MC input fixture)
                │
                ▼
        SimResult(config, sim_type, x, y, metadata, mesh=None)

ParametricSweep(base_config, param, values, sim_fn, sim_kwargs)
        │  for v in values: dataclasses.replace(base_config, {param:v})
        ▼
        list[SimResult]   ← calls any facade above N times
```

### Recommended Project Structure

Design spec §2 shows facades split across `simulation.py`, `damage.py`, `microdosimetry.py`. **However, Phase 36 placed both `run_cv` and `run_field` in a single `simulation.py`.** Recommendation: follow the _actual_ Phase 36 convention (one `simulation.py`) unless the planner prefers the spec's split; the acceptance criteria only require `from petringa import run_X` to work, which is import-location-agnostic.

```
petringa/api/
├── simulation.py    # all run_* facades (extends existing file) — or split per spec §2
├── sweep.py         # ParametricSweep class (new)  [or inline in simulation.py]
├── device.py        # DeviceConfig, build_device (existing, unchanged)
└── results.py       # SimResult, MeshData (existing, unchanged)
petringa/__init__.py # add 7 run_* + ParametricSweep to imports + __all__
tests/
├── test_api_cce.py             # slow integration (physics gate)
├── test_api_facades.py         # light contract tests for the other 5 devsim facades
├── test_api_microdosimetry.py  # data-pipeline test (uses data/synthetic_mc_events.csv)
└── test_api_sweep.py           # ParametricSweep unit tests with FAKE sim_fn (fast, no devsim)
```

### Anti-Patterns to Avoid

- **Refactoring core physics fns to accept `device_info`:** breaks 20 notebooks + is a structural change. Adapt config→kwargs instead.
- **Hand-rolling CCE from `add_generation_to_dd`+`compute_cce_from_dd`:** reconstructs physics = forbidden. Wrap `cce_vs_bias`.
- **Double device cleanup:** bucket-(a) core fns already delete their device. Adding a facade-level `delete_device` on the same name will fail/log-warn spuriously.
- **Mislabeling axes (the CR-01 bug from Phase 36):** don't return a non-physical array as `SimResult.x`. For CCE, `x=voltages`, `y=cce_values`. For unsupported cases, return empty arrays + document, mirroring the Phase 36 fix.

## Don't Hand-Roll

| Problem                    | Don't Build                                    | Use Instead                                                      | Why                                                 |
| -------------------------- | ---------------------------------------------- | ---------------------------------------------------------------- | --------------------------------------------------- |
| CCE vs bias solve          | Custom devsim generation + current integration | `charge_collection.cce_vs_bias`                                  | Already validated in v3.0; physics-change ban       |
| Fluence damage sweep       | Manual carrier removal + re-solve loop         | `charge_collection.cce_vs_fluence`                               | Fresh-device-per-fluence pattern already correct    |
| Dark current decomposition | Manual SRH/TAT/SRV extraction                  | `dark_current.create_dark_current_device` + `dark_current_sweep` | TAT+SRV model setup is intricate                    |
| Transient pulse solve      | Custom BDF1 time-stepping                      | `transient.TransientSolver`                                      | `charge_error=1e10` + adaptive dt already tuned     |
| y-spectrum computation     | Manual histogram + normalization               | `microdosimetry.lineal_energy_spectrum`                          | ICRU-36 binning + normalization validation built in |
| MC CSV parsing             | Custom `pd.read_csv` + unit conversion         | `mc_coupling.load_mc_events_csv`                                 | Handles column mapping + pos/energy unit conversion |
| Config cloning in sweep    | Manual dict copy + setattr                     | `dataclasses.replace(base_config, **{param: value})`             | Immutable, typo-safe, idiomatic                     |

**Key insight:** Every physics operation this phase needs already exists and is validated. The phase's entire value is the thin adapter/repackaging layer. Any line of new _physics_ is a bug (violates the v5.0 refactor-only constraint).

## Runtime State Inventory

> This is a refactor/wrapping phase (no rename, no data migration). Runtime-state categories are assessed for completeness:

| Category            | Items Found                                                                                                                                                                                          | Action Required         |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------- |
| Stored data         | None — facades read/compute, do not persist keyed state.                                                                                                                                             | none                    |
| Live service config | None — no external services touched.                                                                                                                                                                 | none                    |
| OS-registered state | **devsim process-global device registry** — devices created by core fns live in the shared process. Buckets (a)/(c) self-clean or never build; bucket (b) facades must `delete_device` in `finally`. | facade lifecycle (code) |
| Secrets/env vars    | None.                                                                                                                                                                                                | none                    |
| Build artifacts     | None — no package rename; no egg-info churn.                                                                                                                                                         | none                    |

**Nothing found requiring data migration.** The only runtime concern is devsim's process-global device registry, addressed by facade-level cleanup (Pitfall 1).

## Common Pitfalls

### Pitfall 1: devsim process-global device exhaustion across facades + ParametricSweep

**What goes wrong:** devsim state is process-global (STATE.md documented blocker: bare single-process `pytest -q` is unsatisfiable due to device/resource exhaustion). Bucket-(a) core fns build a device internally; if a facade _also_ builds/leaks one, or if `ParametricSweep.run()` calls N devsim facades in one process, resources exhaust and later solves fail.
**Why it happens:** Each devsim device consumes process-global solver resources that don't fully release without explicit `delete_device` / `reset_devsim_fully`.
**How to avoid:**

- Verify per-facade whether the core fn cleans up. `cce_vs_bias` DOES (`finally: delete_device`, line 54-56) → facade must NOT double-delete. [VERIFIED]
- For bucket-(b) facades (`run_dark_current`, `run_transient`) that build via `build_device`/`create_dark_current_device`, wrap with `reset_devsim_fully()` at entry + `delete_device` in `finally` (exact Phase 36 pattern).
- **Test `ParametricSweep` with a fake `sim_fn`** (a lambda returning a canned `SimResult`) so its unit tests never build a device — this satisfies criterion 3 AND sidesteps exhaustion.
- All slow tests use per-file isolation (`.venv/bin/pytest tests/test_api_cce.py -q`), per STATE.md convention. Never run a monolithic `pytest -q`.
  **Warning signs:** A facade test passes alone but fails when run after another slow test in the same process.

### Pitfall 2: CCE facade doping calibration mismatch (calibration divergence)

**What goes wrong:** `cce_vs_bias` hardcodes `N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1.0e-4` in its internal `create_dd_device` call. These **differ** from `DeviceConfig` defaults (`2.93e15, 8.82e13, 0.987e-4`). [VERIFIED: charge_collection.py:479-483 vs device.py DeviceConfig defaults]
**Why it happens:** The v3.0 CCE calibration used slightly different values than the v4.0 C-V-calibrated DeviceConfig defaults.
**The dilemma ("no physics changes" cuts both ways):**

- If `run_cce` passes `config.N_D_junction` etc. through `device_kwargs`, results diverge from the CCE notebooks (a physics-behavior change).
- If `run_cce` ignores config doping and uses `cce_vs_bias` defaults, the facade silently ignores the user's `DeviceConfig` (a contract surprise).
  **How to avoid:** **This is a decision the planner/user must lock — do NOT choose silently.** See Open Question #1. Recommended default: pass config values through so `run_cce(DeviceConfig())` is self-consistent, and document that CCE defaults differ from the notebook's hardcoded calibration. [ASSUMED — needs confirmation]
  **Warning signs:** `run_cce(DeviceConfig())` produces CCE values that don't match notebook 02/03 outputs.

### Pitfall 3: Core fns return pandas DataFrames, not dicts

**What goes wrong:** `sweep_cce_vs_temperature`, `radiation_hardness_sweep`, and `transient_cce_vs_dose_rate` return `pd.DataFrame`, not the `{"voltages":..., "cce_values":...}` dict shape that `cce_vs_bias`/`cce_vs_fluence` return. A facade written assuming a dict will crash. [VERIFIED: temperature_sweep.py:215, transient.py:545]
**How to avoid:** Per-facade, check the actual return type and extract the right columns into `SimResult.x`/`.y` (e.g., temperature sweep: `x=T`, `y=CCE`, with `V` fixed or in metadata).

### Pitfall 4: Radiation damage kappa values are data-blocked placeholders

**What goes wrong:** `radiation_damage.py` NIEL hardness factors are explicit placeholders ("obtain from SR-NIEL", lines 62-65). Absolute Φ_crit numbers are unvalidated (STATE.md blocker). [VERIFIED]
**How to avoid:** `run_radiation_damage` must NOT present outputs as validated. This is not a Phase 37 blocker (the facade just wraps existing behavior), but the eventual UI (FEAT-01) shows a warning banner. Keep the facade honest; don't fabricate values.

### Pitfall 5: run_microdosimetry needs a real MC input file

**What goes wrong:** `run_microdosimetry` and its test are dead without an actual MC CSV. `lineal_energy_spectrum` takes energies, not a config alone.
**How to avoid:** Use the existing fixture `data/synthetic_mc_events.csv` (verified present) via `mc_coupling.load_mc_events_csv`. The test can point at it directly. [VERIFIED: data/synthetic_mc_events.csv exists; test_mc_coupling.py + test_microdosimetry.py already exercise this path]

## Code Examples

### run_cce (bucket a — flagship)

```python
# Source: pattern derived from petringa/api/simulation.py::run_cv (Phase 36)
#         wrapping petringa/core/charge_collection.py::cce_vs_bias
def run_cce(config, v_start=-10.0, v_stop=-200.0, n_points=30):
    if config.half_width_um is not None:
        raise NotImplementedError("run_cce: 2D CCE out of scope; pass half_width_um=None")
    bias_array = np.linspace(v_start, v_stop, n_points)
    # cce_vs_bias builds AND deletes its own device (finally block) -> no facade cleanup here
    result = cce_vs_bias(
        V_range=bias_array,
        epi_thickness_cm=config.epi_thickness_um * 1e-4,
        device_kwargs={  # OPEN QUESTION #1: whether to forward config doping — see below
            "N_D_junction": config.N_D_junction,
            "N_D_bulk": config.N_D_bulk,
            "L_transition": config.L_transition_um * 1e-4,
        },
    )
    return SimResult(
        config=config, sim_type="cce",
        x=result["voltages"], y=result["cce_values"],
        metadata={"I_collected": result["I_collected"], "I_generated": result["I_generated"]},
        mesh=None,
    )
```

Note: `cce_vs_bias` return keys are `"voltages"`, `"cce_values"`, `"I_collected"`, `"I_generated"`. [VERIFIED: charge_collection.py:432-597]

### run_microdosimetry (bucket c — data pipeline, no devsim)

```python
# Source: petringa/core/mc_coupling.py + petringa/core/microdosimetry.py
def run_microdosimetry(config, mc_csv_path, sv_thickness_um=10.0, sv_width_um=None):
    events = load_mc_events_csv(mc_csv_path)                 # DataFrame of MC events
    l_bar = mean_chord_length(sv_thickness_um, sv_width_um=sv_width_um)
    # collected energies come from events (column per mc_coupling schema)
    spec = lineal_energy_spectrum(collected_energies_keV, l_bar)
    return SimResult(
        config=config, sim_type="microdosimetry",
        x=spec["bin_centers"], y=spec["d_y"] * spec["bin_centers"],  # y*d(y)
        metadata={"y_F": spec["y_F"], "y_D": spec["y_D"], "f_y": spec["f_y"]},
        mesh=None,
    )
```

[VERIFIED: mc_coupling.py:107, microdosimetry.py:36,117]

### ParametricSweep

```python
# Source: design spec §3.5
from dataclasses import replace

@dataclass
class ParametricSweep:
    base_config: "DeviceConfig"
    param: str
    values: list
    sim_fn: "Callable"
    sim_kwargs: dict = field(default_factory=dict)

    def run(self) -> list:
        results = []
        for value in self.values:
            cfg_i = replace(self.base_config, **{self.param: value})
            results.append(self.sim_fn(cfg_i, **self.sim_kwargs))
        return results
```

### ParametricSweep unit test with fake sim_fn (satisfies criterion 3, no devsim)

```python
def test_parametric_sweep_returns_list_of_correct_length():
    def fake_sim(cfg, **kw):
        return SimResult(config=cfg, sim_type="fake",
                         x=np.array([cfg.epi_thickness_um]), y=np.array([1.0]))
    sweep = ParametricSweep(base_config=DeviceConfig(), param="epi_thickness_um",
                            values=[5, 10, 20], sim_fn=fake_sim)
    results = sweep.run()
    assert len(results) == 3
    assert [r.x[0] for r in results] == [5, 10, 20]   # config injection verified
```

## State of the Art

| Old Approach                                                     | Current Approach                                              | When Changed       | Impact                                                 |
| ---------------------------------------------------------------- | ------------------------------------------------------------- | ------------------ | ------------------------------------------------------ |
| Notebooks import `petringa.core.*` and call physics fns directly | Public facades in `petringa.api.*` return uniform `SimResult` | v5.0 (Phase 36-37) | UI + external users never touch core                   |
| Phase 36 `build_device → core(device_info)` device-passing       | Phase 37 config→kwargs adapter (buckets a/c)                  | Phase 37           | Core fns keep building own devices; notebooks unbroken |

**Deprecated/outdated:** Nothing deprecated. Core physics is frozen (refactor-only milestone).

## Assumptions Log

| #   | Claim                                                                                                                                                                                       | Section                       | Risk if Wrong                                                                                                                                                                     |
| --- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| A1  | `run_cce` should forward `config` doping values into `cce_vs_bias.device_kwargs` (making `run_cce(DeviceConfig())` self-consistent) rather than use the core fn's hardcoded CCE calibration | Pitfall 2, Code Examples      | If wrong, `run_cce` either silently ignores user config OR diverges from CCE notebook results — this is the one real design decision in the phase; must be locked by planner/user |
| A2  | Single-file `simulation.py` (Phase 36 actual convention) is preferred over the spec's `damage.py`/`microdosimetry.py` split                                                                 | Recommended Project Structure | Low — import location is invisible to acceptance criteria (`from petringa import run_X`)                                                                                          |
| A3  | `transient_cce_vs_dose_rate` (self-building) is an acceptable minimal satisfier for `run_transient` given the criterion-2 weak bar and the spec omitting a `run_transient` signature        | Bucket (b), Open Q #3         | Low — criterion only requires import + `DeviceConfig` first arg; but the "natural" signature is undefined                                                                         |
| A4  | `data/synthetic_mc_events.csv` is the correct default MC fixture for `run_microdosimetry` tests                                                                                             | Pitfall 5, Environment        | Low — file verified present and already used by existing MC/microdosim tests                                                                                                      |

## Open Questions

1. **CCE doping calibration (HIGH priority — the phase's one real decision).**
   - What we know: `cce_vs_bias` hardcodes `N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1.0e-4`; `DeviceConfig` defaults are `2.93e15, 8.82e13, 0.987e-4`.
   - What's unclear: Should `run_cce` forward `config` doping (self-consistent facade, but CCE numbers shift vs notebooks) or ignore it (matches notebooks, but ignores user config)?
   - Recommendation: Forward config values (A1) and document the difference; but flag for discuss-phase / user confirmation before locking. "No physics changes" is ambiguous here — surface it, don't decide silently.

2. **Facade file layout: single `simulation.py` vs spec's per-domain split (`damage.py`, `microdosimetry.py`).**
   - What we know: Phase 36 put both facades in one `simulation.py`; spec §2 shows a split.
   - Recommendation: Follow Phase 36 (single file) for consistency unless planner prefers spec; acceptance criteria are location-agnostic.

3. **`run_transient` signature — undefined.**
   - What we know: Design spec §3.4 lists every facade signature EXCEPT `run_transient`. Core offers both `TransientSolver` (device_info + init + pulse) and the self-building `transient_cce_vs_dose_rate`.
   - Recommendation: Implement the minimal satisfier of criterion 2 (import + `DeviceConfig` first arg). Consider `run_transient(config, dose_rate_Gy_s=..., ...)` wrapping `transient_cce_vs_dose_rate`. Flag as underspecified.

4. **`run_temperature_sweep`/`run_flash_recombination`/`run_transient` SimResult axis semantics.**
   - What we know: These return DataFrames or multi-column data (T×V grid for temperature). A single `x`/`y` pair must be chosen.
   - Recommendation: Pick the primary swept axis as `x` (T for temperature sweep, dose_rate for flash/transient), CCE as `y`, put the rest (fixed V, other columns) in `metadata`. Confirm with planner.

## Environment Availability

| Dependency                     | Required By                                                                                                    | Available                               | Version                    | Fallback                          |
| ------------------------------ | -------------------------------------------------------------------------------------------------------------- | --------------------------------------- | -------------------------- | --------------------------------- |
| devsim                         | run_cce, run_radiation_damage, run_dark_current, run_temperature_sweep, run_flash_recombination, run_transient | ✓ (installed, used by all prior phases) | as installed               | — (no fallback; core requirement) |
| numpy / scipy / pandas         | all facades + repackaging                                                                                      | ✓                                       | declared in pyproject.toml | —                                 |
| `data/synthetic_mc_events.csv` | run_microdosimetry + its test                                                                                  | ✓                                       | —                          | none needed (present)             |
| pytest 8.x + `slow` marker     | all tests                                                                                                      | ✓                                       | 8.x (`.venv/bin/pytest`)   | —                                 |

**Missing dependencies with no fallback:** none.
**Missing dependencies with fallback:** none.

## Validation Architecture

### Test Framework

| Property           | Value                                                                                                                                     |
| ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------- |
| Framework          | pytest 8.x (`.venv/bin/pytest`)                                                                                                           |
| Config file        | `pytest.ini` (registers `slow` marker only)                                                                                               |
| Quick run command  | `.venv/bin/pytest tests/test_api_<name>.py -q` (per-file isolation)                                                                       |
| Full suite command | N/A by project convention — monolithic `pytest -q` is unsatisfiable (devsim process exhaustion, STATE.md/PKG-03). Use per-file isolation. |

### Phase Requirements → Test Map

| Req ID | Behavior                                                                     | Test Type                                          | Automated Command                                                    | File Exists? |
| ------ | ---------------------------------------------------------------------------- | -------------------------------------------------- | -------------------------------------------------------------------- | ------------ |
| LIB-04 | `run_cce(DeviceConfig())` → SimResult, CCE ∈ [0,1]                           | slow integration                                   | `.venv/bin/pytest tests/test_api_cce.py -x`                          | ❌ Wave 0    |
| LIB-04 | `run_cce` rejects 2D config                                                  | unit (fast)                                        | `.venv/bin/pytest tests/test_api_cce.py::test_run_cce_rejects_2d -x` | ❌ Wave 0    |
| LIB-06 | all 6 facades import + accept DeviceConfig first arg                         | contract (fast, no solve — signature/import check) | `.venv/bin/pytest tests/test_api_facades.py -x`                      | ❌ Wave 0    |
| LIB-06 | `run_microdosimetry` produces y·d(y) spectrum from CSV                       | data-pipeline (fast, no devsim)                    | `.venv/bin/pytest tests/test_api_microdosimetry.py -x`               | ❌ Wave 0    |
| LIB-07 | `ParametricSweep.run()` over ≥2 values → `list[SimResult]` of correct length | unit (fast, fake sim_fn)                           | `.venv/bin/pytest tests/test_api_sweep.py -x`                        | ❌ Wave 0    |

### Sampling Rate

- **Per task commit:** run the single affected test file in isolation (e.g., `.venv/bin/pytest tests/test_api_sweep.py -q`).
- **Per wave merge:** re-run that wave's test file(s) + a regression spot-check of the closest core test (`tests/test_charge_collection.py -q` for CCE work).
- **Phase gate:** all new `tests/test_api_*.py` green in per-file isolation before `/gsd:verify-work`. NEVER a monolithic run.

### Wave 0 Gaps

- [ ] `tests/test_api_cce.py` — covers LIB-04 (slow integration + fast 2D-guard). Mirror `tests/test_api_cv.py` structure.
- [ ] `tests/test_api_facades.py` — covers LIB-06 for the 5 devsim facades. Prefer **fast** import/signature contract tests (criterion 2 is "imports + accepts DeviceConfig first arg" — a `inspect.signature` + import assertion satisfies it) plus optionally one slow smoke each if the planner wants runtime proof.
- [ ] `tests/test_api_microdosimetry.py` — covers LIB-06 microdosimetry via `data/synthetic_mc_events.csv` (fast, no devsim).
- [ ] `tests/test_api_sweep.py` — covers LIB-07 with a **fake sim_fn** (fast, no devsim).

_Recommendation: bias toward FAST contract/data tests. Only `run_cce` needs a slow integration test (its physics gate). Over-investing in slow validated-output tests for the other 6 exceeds the acceptance bar and worsens the devsim-exhaustion risk._

## Security Domain

> `security_enforcement` is absent from `.planning/config.json` → treated as enabled. This phase is an internal library API over an existing scientific compute core — no auth, no network, no untrusted input surface beyond a local MC CSV path.

### Applicable ASVS Categories

| ASVS Category         | Applies | Standard Control                                                                                                                                                                                                                                                                                           |
| --------------------- | ------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| V2 Authentication     | no      | No auth surface (library, local execution)                                                                                                                                                                                                                                                                 |
| V3 Session Management | no      | No sessions                                                                                                                                                                                                                                                                                                |
| V4 Access Control     | no      | No access control surface                                                                                                                                                                                                                                                                                  |
| V5 Input Validation   | partial | `run_microdosimetry` takes a filesystem path (`mc_csv_path`) — parsing is delegated to `mc_coupling.load_mc_events_csv` (pandas `read_csv`). ParametricSweep `param` is a config field name → `dataclasses.replace` raises on unknown field (safe, no `eval`). Facade numeric args should not be `eval`'d. |
| V6 Cryptography       | no      | No crypto                                                                                                                                                                                                                                                                                                  |

### Known Threat Patterns for this stack

| Pattern                                        | STRIDE            | Standard Mitigation                                                                                                                                                                |
| ---------------------------------------------- | ----------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Malicious/oversized MC CSV via `mc_csv_path`   | Tampering / DoS   | Path is caller-supplied local file; `load_mc_events_csv` uses pandas (no code exec). Not a new surface — existing MC tests already parse CSVs. No `pickle`/`eval` on file content. |
| `ParametricSweep.param` as attribute injection | Tampering         | Use `dataclasses.replace(cfg, **{param: value})` — raises `TypeError` on unknown field; NEVER `setattr`+`getattr` with `eval`.                                                     |
| devsim resource exhaustion (self-DoS)          | Denial of Service | Per-file test isolation + facade cleanup (Pitfall 1). Not adversarial; operational.                                                                                                |

No new external attack surface is introduced. The phase adds no network, no deserialization of untrusted binary, no dynamic code execution.

## Sources

### Primary (HIGH confidence — in-repo verification)

- `petringa/api/simulation.py`, `petringa/api/device.py`, `petringa/api/results.py` — Phase 36 facade template, `SimResult`/`MeshData`/`DeviceConfig` contracts
- `petringa/core/charge_collection.py` (cce_vs_bias:432, cce_vs_fluence:597, self-cleanup:54-56) — CCE + fluence entry points and device-lifecycle behavior
- `petringa/core/dark_current.py` (create_dark_current_device:538, dark_current_sweep:435) — dark current entry + required model setup
- `petringa/core/temperature_sweep.py` (sweep_cce_vs_temperature:215), `flash_recombination.py` (cce_vs_dose_rate:263), `transient.py` (TransientSolver:172, transient_cce_vs_dose_rate:545), `radiation_damage.py` (kappa placeholders:62-65)
- `petringa/core/mc_coupling.py` (load_mc_events_csv:107), `microdosimetry.py` (mean_chord_length:36, lineal_energy_spectrum:117) — microdosimetry pipeline
- `docs/superpowers/specs/2026-06-26-simulator-library-ui-design.md` §3.4-§3.5 — canonical facade + ParametricSweep signatures
- `tests/test_api_cv.py`, `tests/test_api_field.py`, `tests/test_api_device.py` — Phase 36 test conventions
- `.planning/REQUIREMENTS.md` (LIB-04/06/07), `.planning/STATE.md` (refactor-only constraint, devsim exhaustion blocker), `.planning/phases/36-*/36-REVIEW-FIX.md` + `36-VALIDATION.md` — Phase 36 lessons

### Secondary (MEDIUM confidence)

- None required — all claims verified against primary in-repo sources.

### Tertiary (LOW confidence)

- None.

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH — no new packages; deps verified in pyproject.toml + codebase
- Architecture (3-bucket classification): HIGH — each core fn signature read directly from source; device-lifecycle behavior verified
- Pitfalls: HIGH — CCE calibration mismatch, DataFrame returns, and self-cleanup all confirmed by reading exact source lines
- The one genuine ambiguity (CCE doping forwarding) is surfaced as Open Question #1 / A1, not resolved silently

**Research date:** 2026-07-08
**Valid until:** ~2026-08-08 (30 days — codebase is stable, refactor-only milestone; core physics frozen)
