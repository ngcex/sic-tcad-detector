# Phase 37: Core API — CCE + Remaining Facades + ParametricSweep - Pattern Map

**Mapped:** 2026-07-08
**Files analyzed:** 12 (2 modified: `simulation.py`, `__init__.py`; 2 new source: `sweep.py`, and 4 new test files; plus classification of the 7 facade functions + ParametricSweep)
**Analogs found:** 10 / 12 (2 partial/no-analog gaps flagged: `run_microdosimetry` events→energies bridge; `run_transient` signature)

> No CONTEXT.md exists for this phase (verified: only `37-RESEARCH.md` present). The `run_cce` doping-forwarding decision (RESEARCH Open Question #1 / A1) is therefore **unresolved** — this map surfaces both concrete facts for the planner rather than picking one.

## File Classification

The **Data Flow** column encodes the three-bucket spine from RESEARCH (not flattened to "request-response"):

- **bucket-a (config→kwargs adapter):** core fn builds AND cleans its own device; facade does NOT call `build_device`/`reset_devsim_fully`/`delete_device`.
- **bucket-b (build+setup lifecycle):** facade builds device + extra model setup, wraps with `reset → try → finally delete` (full Phase 36 pattern).
- **bucket-c (pure data pipeline):** no devsim device at all.

| New/Modified File / Symbol                              | Role    | Data Flow (bucket)                      | Closest Analog                                                                     | Match Quality                                    |
| ------------------------------------------------------- | ------- | --------------------------------------- | ---------------------------------------------------------------------------------- | ------------------------------------------------ |
| `etna/api/simulation.py :: run_cce`                 | facade  | CRUD / bucket-a                         | `simulation.py::run_cv` (structure) + `run_field` (2D guard)                       | role-match (device lifecycle differs — bucket-a) |
| `etna/api/simulation.py :: run_radiation_damage`    | facade  | CRUD / bucket-a                         | `simulation.py::run_cv` (shell)                                                    | role-match                                       |
| `etna/api/simulation.py :: run_temperature_sweep`   | facade  | transform / bucket-a (DataFrame return) | `simulation.py::run_cv` (shell)                                                    | role-match                                       |
| `etna/api/simulation.py :: run_flash_recombination` | facade  | CRUD / bucket-a                         | `simulation.py::run_cv` (shell)                                                    | role-match                                       |
| `etna/api/simulation.py :: run_dark_current`        | facade  | CRUD / bucket-b                         | `simulation.py::run_cv` (full reset→try→finally lifecycle)                         | exact (lifecycle)                                |
| `etna/api/simulation.py :: run_transient`           | facade  | streaming / bucket-b (or a)             | `simulation.py::run_cv` (shell); signature UNDEFINED in spec                       | partial (no signature analog)                    |
| `etna/api/simulation.py :: run_microdosimetry`      | facade  | file-I/O / bucket-c                     | `mc_coupling`/`microdosimetry` core fns (no facade analog)                         | partial (events→energies bridge has NO analog)   |
| `etna/api/sweep.py :: ParametricSweep`              | utility | batch                                   | none (new orchestration class); spec §3.5 is canonical                             | no analog (fully specified by spec)              |
| `etna/__init__.py` (modify)                         | config  | —                                       | existing `__init__.py` (add 7 fns + ParametricSweep)                               | exact                                            |
| `tests/test_api_cce.py` (new)                           | test    | slow integration + fast guard           | `tests/test_api_cv.py`                                                             | exact                                            |
| `tests/test_api_facades.py` (new)                       | test    | fast contract (import/signature)        | `tests/test_api_device.py` (fast-assertion style)                                  | role-match                                       |
| `tests/test_api_microdosimetry.py` (new)                | test    | fast data-pipeline                      | `tests/test_microdosimetry.py` (fixtures) + `tests/test_mc_coupling.py` (CSV load) | role-match                                       |
| `tests/test_api_sweep.py` (new)                         | test    | fast unit (fake sim_fn)                 | `tests/test_api_device.py` (fast unit style)                                       | role-match                                       |

## Pattern Assignments

### `run_cce` (facade, bucket-a) — flagship

**Analogs:** `etna/api/simulation.py::run_cv` (SimResult packaging + 2D guard), wrapping `etna/core/charge_collection.py::cce_vs_bias`.

**2D guard pattern** — copy from `run_cv` (simulation.py:75-81):

```python
if config.half_width_um is not None:
    raise NotImplementedError(
        "run_cce: 2D CCE is out of Phase 37 scope. cce_vs_bias operates on the "
        "1D DD device (create_dd_device); pass half_width_um=None."
    )
```

**CRITICAL — do NOT clone run_cv's device lifecycle.** `cce_vs_bias` builds AND deletes its own device in a `finally` block (charge_collection.py:583-587):

```python
    finally:
        try:
            devsim.delete_device(device=device)
        except Exception:
            pass
```

Therefore the facade must NOT call `reset_devsim_fully()`, `build_device()`, or `delete_device()`. Adding a facade-level cleanup double-deletes the same device (RESEARCH Anti-Pattern "Double device cleanup").

**Core call signature** (charge_collection.py:432-438):

```python
def cce_vs_bias(V_range, epi_thickness_cm=10e-4, alpha_range_cm=15e-4,
                generation_rate=1e18, device_kwargs=None):
```

Returns dict (charge_collection.py:589-594): `{"voltages", "cce_values", "I_collected", "I_generated"}`.

**SimResult packaging** — mirror `run_cv` (simulation.py:105-112): `x=result["voltages"]`, `y=result["cce_values"]`, `metadata={"I_collected", "I_generated"}`, `mesh=None`, `sim_type="cce"`.

**⚠ UNRESOLVED DECISION (planner must lock) — doping forwarding via `device_kwargs`:**
The one real design choice in this phase (RESEARCH Open Question #1 / A1), unresolved because no CONTEXT.md exists. The two concrete facts:

| Source                                                         | N_D_junction | N_D_bulk  | L_transition                            |
| -------------------------------------------------------------- | ------------ | --------- | --------------------------------------- |
| `cce_vs_bias` hardcoded default (charge_collection.py:481-483) | `2.90e15`    | `8.50e13` | `1.0e-4` cm                             |
| `DeviceConfig` defaults (device.py:37-41)                      | `2.93e15`    | `8.82e13` | `0.987e-4` cm (`L_transition_um=0.987`) |

- **Option A (forward config):** pass `device_kwargs={"N_D_junction": config.N_D_junction, "N_D_bulk": config.N_D_bulk, "L_transition": config.L_transition_um*1e-4}` → `run_cce(DeviceConfig())` is self-consistent, but CCE numbers shift vs the v3.0 CCE notebooks.
- **Option B (ignore config):** omit `device_kwargs` → matches notebooks, but silently ignores the user's `DeviceConfig` doping.

RESEARCH recommends Option A + documentation, but flags it for explicit lock. **Do not choose silently in the plan.** Note: `epi_thickness_cm=config.epi_thickness_um * 1e-4` is unambiguous and forwarded either way.

---

### `run_radiation_damage` (facade, bucket-a)

**Analog:** `run_cv` shell (SimResult packaging only — no device lifecycle in facade). Wraps `etna/core/charge_collection.py::cce_vs_fluence` (charge_collection.py:597-608):

```python
def cce_vs_fluence(fluence_range, V_bias=-40.0, epi_thickness_cm=10e-4,
                   N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1.0e-4,
                   alpha_range_cm=15e-4, generation_rate=1e18, energy_MeV=62.0,
                   lifetime_model="linear", damage_params=None):
```

Returns dict (charge_collection.py:648-654): `{"fluences", "cce_values", "V_bias", "energy_MeV", "lifetime_model"}`. Builds a fresh device per fluence internally (self-cleaning, bucket-a).

**Design spec §3.4:** `run_radiation_damage(cfg, fluences=[...], proton_energy_MeV=5.6)` → `x=fluence, y=CCE_post_damage`. Map: `x=result["fluences"]`, `y=result["cce_values"]`, `sim_type="damage"`, rest to metadata.

**Honesty note (RESEARCH Pitfall 4):** NIEL kappa factors in `radiation_damage.py:62-65` are data-blocked placeholders. Do NOT present outputs as validated; keep metadata factual. Same doping-default mismatch as `run_cce` applies (this fn also hardcodes `2.90e15/8.50e13/1.0e-4` as defaults) — same planner decision.

---

### `run_temperature_sweep` (facade, bucket-a, DataFrame return)

**Analog:** `run_cv` shell. Wraps `etna/core/temperature_sweep.py::sweep_cce_vs_temperature` (temperature_sweep.py:215-219):

```python
def sweep_cce_vs_temperature(temperatures, voltages=None, method="hecht", **device_kwargs):
```

**⚠ Returns `pd.DataFrame`, NOT a dict** (temperature_sweep.py:237-238): long-format columns `T, V, CCE` (RESEARCH Pitfall 3). Facade must extract arrays:

```python
df = sweep_cce_vs_temperature(temperatures=temps, voltages=voltages, method=method)
# choose primary swept axis as x=T, y=CCE; fixed V / other columns -> metadata
# (RESEARCH Open Q #4: confirm axis choice with planner)
x = df["T"].to_numpy(); y = df["CCE"].to_numpy()
```

`sim_type="temperature"`. Note `voltages` defaults to `[-10,-20,-30]` internally → if a single V is desired for a clean 1D `x=T` curve, pass a single-element `voltages`.

---

### `run_flash_recombination` (facade, bucket-a)

**Analog:** `run_cv` shell. Wraps `etna/core/flash_recombination.py::cce_vs_dose_rate` (flash_recombination.py:263-272):

```python
def cce_vs_dose_rate(dose_rates_Gy_s, V_bias=-30.0, epi_thickness_cm=10e-4, E_MeV=62,
                     n_continuation_steps=5, N_D_junction=2.90e15, N_D_bulk=8.50e13,
                     L_transition=1.0e-4):
```

Returns dict (flash_recombination.py:305+): `{"dose_rates", "cce_values", "cce_no_auger_ref", "V_bias", "epi_thickness_cm", "E_MeV", ...}`. Builds device internally (bucket-a). Map `x=result["dose_rates"]`, `y=result["cce_values"]`, `sim_type="flash"`.

**Honesty note:** module carries an explicit SCOPE/LIMITATIONS block — outputs are a numerical sensitivity bound, not validated FLASH physics. Do not overstate.

---

### `run_dark_current` (facade, bucket-b) — full Phase 36 lifecycle

**Analog:** `run_cv` device lifecycle (simulation.py:86-123) — this is an **exact** structural match. Two-step core call:

1. Build device WITH TAT+SRV setup (do NOT use bare `build_device()` — it omits TAT/SRV). `etna/core/dark_current.py::create_dark_current_device` (dark_current.py:538-564):

```python
def create_dark_current_device(T=300, N_t=None, S_n=None, S_p=None, **kwargs):
    device_info = create_dd_device(T=T, **kwargs)
    setup_tat_model(device_info, N_t=N_t)
    setup_surface_recombination(device_info, S_n=S_n, S_p=S_p)
    return device_info
```

2. Sweep: `dark_current.py::dark_current_sweep(device_info, V_range, area=0.05, V_step=0.5)` (dark_current.py:435). Returns dict: `{voltages, I_total, I_SRH, I_TAT, I_SRV, J_total, J_SRH, J_TAT, J_SRV}`.

**Lifecycle wrapper — copy exactly from run_cv (simulation.py:87-123):**

```python
reset_devsim_fully()
device_info = create_dark_current_device(T=config.T, ...)   # replaces build_device()
try:
    result = dark_current_sweep(device_info, V_range=bias_array, area=config.area_cm2)
    return SimResult(config=config, sim_type="dark_current",
                     x=result["voltages"], y=result["I_total"],
                     metadata={"I_SRH":..., "I_TAT":..., "I_SRV":...}, mesh=None)
finally:
    try:
        devsim.delete_device(device=device_info["device_name"])
    except Exception:
        logger.warning(..., exc_info=True)
        reset_devsim_fully()
```

---

### `run_transient` (facade, bucket-b) — ⚠ SIGNATURE UNDEFINED

**Analog:** `run_cv` shell for SimResult packaging; NO signature analog. Design spec §3.4 omits `run_transient` entirely (RESEARCH Open Q #3). Acceptance bar is only criterion-2: imports + accepts `DeviceConfig` first arg.

**Do not invent a signature in the plan without flagging it undefined.** Two core entry points exist:

- Minimal satisfier (self-building, bucket-a shaped): `etna/core/transient.py::transient_cce_vs_dose_rate(V_bias=-30.0, dose_rates=None, t_rise=1e-6, t_duration=1e-3, t_fall=1e-6, dt_min=1e-8, dt_max=1e-4, epi_thickness_cm=10e-4)` (transient.py:545). **⚠ Returns `pd.DataFrame`** with columns `["dose_rate_Gy_s", "transient_cce"]` (transient.py:582-583) — Pitfall 3 applies.
- Full solver (bucket-b): `transient.TransientSolver(device_info, contact, method)` needs `.initialize()` then `.simulate_pulse(...)` (transient.py:172-431).

**Recommendation (RESEARCH A3):** wrap `transient_cce_vs_dose_rate` for the minimal bar; leave final signature to planner. Map `x=df["dose_rate_Gy_s"]`, `y=df["transient_cce"]`, `sim_type="transient"`.

---

### `run_microdosimetry` (facade, bucket-c) — pure data pipeline, NO devsim

**Analogs:** `etna/core/mc_coupling.py::load_mc_events_csv` + `etna/core/microdosimetry.py::{mean_chord_length, lineal_energy_spectrum}`. No `reset_devsim_fully`, no `build_device`, no cleanup.

**Pipeline:**

```python
events = load_mc_events_csv(mc_csv_path)          # mc_coupling.py:107 -> DataFrame
l_bar = mean_chord_length(sv_thickness_um, sv_width_um=sv_width_um)   # microdosimetry.py:36
spec = lineal_energy_spectrum(collected_energies_keV, l_bar)         # microdosimetry.py:117
```

`lineal_energy_spectrum` returns dict (microdosimetry.py:208+): `{bin_edges, bin_centers, bin_widths, f_y, d_y, y_F, y_D, n_events, y_values}`.

**Design spec §3.4:** `run_microdosimetry(cfg, mc_csv_path, sv_thickness_um=10, sv_width_um=150)` → `x=lineal_energy (=bin_centers)`, `y=yd_y`. RESEARCH code example uses `y = spec["d_y"] * spec["bin_centers"]` (the y·d(y) representation). `sim_type="microdosimetry"`.

**⚠ NO ANALOG for the events→energies bridge — see "No Analog Found" below.**

`load_mc_events_csv` returns a per-**step** DataFrame (columns `event_id, x_cm, y_cm, z_cm, edep_keV`); the synthetic fixture has multiple rows per `event_id` (verified: `event_id=0` repeats in `data/synthetic_mc_events.csv`). `lineal_energy_spectrum` wants **per-event collected energy**. The natural bridge is `collected_energies_keV = events.groupby("event_id")["edep_keV"].sum().to_numpy()`, but this exact aggregation has NO existing analog (existing microdosim tests feed synthetic energies directly). `mc_coupling.process_mc_ensemble` (mc_coupling.py:392) exists but requires a `cce_interp` — not a drop-in. Planner must define this aggregation step explicitly.

---

### `ParametricSweep` (utility, batch) — `etna/api/sweep.py` (new)

**No codebase analog** — fully specified by design spec §3.5. Implement exactly (RESEARCH Code Examples):

```python
from dataclasses import dataclass, field, replace

@dataclass
class ParametricSweep:
    base_config: "DeviceConfig"
    param: str
    values: list
    sim_fn: "Callable"
    sim_kwargs: dict = field(default_factory=dict)

    def run(self) -> list["SimResult"]:
        results = []
        for value in self.values:
            cfg_i = replace(self.base_config, **{self.param: value})
            results.append(self.sim_fn(cfg_i, **self.sim_kwargs))
        return results
```

**Security (RESEARCH V5):** use `dataclasses.replace(cfg, **{param: value})` — raises `TypeError` on unknown field. NEVER `setattr`/`eval`.

---

### `etna/__init__.py` (modify)

**Analog:** current `__init__.py:1-17` (exact pattern). Add the 7 `run_*` facades + `ParametricSweep` to imports and `__all__`:

```python
from etna.api.simulation import (
    run_cv, run_field, run_cce, run_radiation_damage, run_dark_current,
    run_temperature_sweep, run_flash_recombination, run_transient, run_microdosimetry,
)
from etna.api.sweep import ParametricSweep
```

Add each name to `__all__`. Acceptance criteria require `from etna import run_X` to work (import-location-agnostic).

## Shared Patterns

### Device lifecycle (bucket-b facades only)

**Source:** `etna/api/simulation.py::run_cv` lines 86-123 (reset → build → try → finally delete → fallback reset).
**Apply to:** `run_dark_current`, `run_transient` (if using `TransientSolver` path).
**Do NOT apply to:** `run_cce`, `run_radiation_damage`, `run_temperature_sweep`, `run_flash_recombination` (bucket-a, core self-cleans), `run_microdosimetry` (bucket-c, no device).

```python
reset_devsim_fully()
device_info = <build+setup>(config)
try:
    ...
finally:
    try:
        devsim.delete_device(device=device_info["device_name"])
    except Exception:
        logger.warning("...: delete_device(%r) failed, falling back to full reset",
                       device_info["device_name"], exc_info=True)
        reset_devsim_fully()
```

### 2D NotImplementedError guard

**Source:** `run_cv` (simulation.py:75-81).
**Apply to:** `run_cce` (1D-only, uses `create_dd_device`). Guard is a plain `if config.half_width_um is not None:` check at function top → makes the rejection test FAST (no devsim).

### Unit conversion (um → cm)

**Source:** `build_device` (device.py:71-73): `epi_thickness_cm = config.epi_thickness_um * 1e-4`; `L_transition = config.L_transition_um * 1e-4`.
**Apply to:** all facades forwarding config geometry/doping into core `*_cm` kwargs.

### SimResult packaging

**Source:** `run_cv` (simulation.py:105-112). All facades return `SimResult(config=config, sim_type=<str>, x=..., y=..., metadata={...}, mesh=None)`. `SimResult`/`MeshData` are frozen-ish dataclasses in `results.py:36-45`. Axis rule (RESEARCH Anti-Pattern CR-01): `x` = primary swept physical axis; never mislabel. For unsupported cases return empty arrays, mirroring `run_field`'s 2D handling (simulation.py:277-290).

### Test conventions

**Slow integration test:** `tests/test_api_cv.py` — `@pytest.mark.slow` class, `run_X(DeviceConfig(), ...)`, assert `isinstance(SimResult)`, sim_type, shape, physics gate. Apply to `test_api_cce.py` (CCE ∈ [0,1] gate).
**Fast guard test:** `tests/test_api_cv.py::test_run_cv_rejects_2d_config` — `pytest.raises(NotImplementedError)`, no devsim. Apply to `test_api_cce.py` 2D-guard.
**Fast contract test:** `tests/test_api_device.py` style (plain asserts) + `inspect.signature`/import checks for `test_api_facades.py` (LIB-06: imports + accepts DeviceConfig first arg).
**Data-pipeline test:** `tests/test_microdosimetry.py` fixtures + `tests/test_mc_coupling.py` CSV load, pointed at `data/synthetic_mc_events.csv`.
**Fake sim_fn test:** RESEARCH §"ParametricSweep unit test" — lambda returning a canned `SimResult`, asserts `len(results)==len(values)` + config injection. Never touches devsim.
**Isolation:** run each test file alone (`.venv/bin/pytest tests/test_api_<name>.py -q`) — never monolithic `pytest -q` (devsim process exhaustion, STATE.md/Pitfall 1).

## No Analog Found

| File / Step                                      | Role        | Data Flow            | Reason                                                                                                                                                                                                                                                                                                                                                                                                          |
| ------------------------------------------------ | ----------- | -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `run_microdosimetry` events→energies aggregation | facade step | file-I/O / bucket-c  | `lineal_energy_spectrum` needs per-event collected energy, but `load_mc_events_csv` returns a per-step DataFrame with repeated `event_id`. No existing code performs the `groupby("event_id")["edep_keV"].sum()` bridge; existing microdosim tests feed synthetic energies directly. `mc_coupling.process_mc_ensemble` requires a `cce_interp` (not a drop-in). Planner must define the aggregation explicitly. |
| `run_transient` signature                        | facade      | streaming / bucket-b | Design spec §3.4 omits `run_transient` entirely. No signature analog exists. Core offers both `TransientSolver` (device_info + init + pulse) and self-building `transient_cce_vs_dose_rate` (DataFrame return). Planner must choose + define signature; RESEARCH recommends wrapping `transient_cce_vs_dose_rate` for the minimal criterion-2 bar.                                                              |
| `ParametricSweep`                                | utility     | batch                | No codebase analog, but NOT a gap — design spec §3.5 fully specifies it. Listed here only to note the source is the spec, not an existing file.                                                                                                                                                                                                                                                                 |

## Metadata

**Analog search scope:** `etna/api/` (simulation.py, device.py, results.py), `etna/core/` (charge_collection.py, dark_current.py, temperature_sweep.py, flash_recombination.py, transient.py, mc_coupling.py, microdosimetry.py), `tests/` (test_api_cv.py, test_api_field.py, test_api_device.py, test_microdosimetry.py, test_mc_coupling.py), `data/synthetic_mc_events.csv`.
**Files scanned:** ~15 source/test files + 1 data fixture.
**Pattern extraction date:** 2026-07-08
**Key unresolved item for planner:** `run_cce` doping-forwarding decision (no CONTEXT.md exists to resolve it).
