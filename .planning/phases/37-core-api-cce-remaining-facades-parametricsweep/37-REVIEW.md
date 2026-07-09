---
phase: 37-core-api-cce-remaining-facades-parametricsweep
reviewed: 2026-07-09T00:00:00Z
depth: standard
files_reviewed: 7
files_reviewed_list:
  - petringa/__init__.py
  - petringa/api/simulation.py
  - petringa/api/sweep.py
  - tests/test_api_cce.py
  - tests/test_api_facades.py
  - tests/test_api_microdosimetry.py
  - tests/test_api_sweep.py
findings:
  critical: 0
  warning: 3
  info: 2
  total: 5
status: issues_found
---

# Phase 37: Code Review Report

**Reviewed:** 2026-07-09
**Depth:** standard
**Files Reviewed:** 7
**Status:** issues_found

## Summary

Phase 37 adds `run_cce()` plus six more `run_*` facades and `ParametricSweep` on top of the existing devsim-based core. I traced every facade's dict/DataFrame consumption against the actual return statements of the wrapped core functions (`cce_vs_bias`, `cce_vs_fluence`, `sweep_cce_vs_temperature`, `cce_vs_dose_rate`, `dark_current_sweep`/`create_dark_current_device`, `transient_cce_vs_dose_rate`, `load_mc_events_csv`/`mean_chord_length`/`lineal_energy_spectrum`) rather than trusting the plan/summary claims, since 5 of the 9 facades have signature-only tests and are never executed by the test suite. All key/column accesses match the real core return shapes — no KeyError-class bug found. The bucket lifecycle classification (bucket-a self-cleaning / bucket-b facade-owned / bucket-c no-devsim) is applied correctly and consistently: no facade double-deletes a device, and `run_dark_current` is the only one that owns a full reset→build→try→finally-delete lifecycle. `ParametricSweep` correctly uses `dataclasses.replace` (no `setattr`/`eval`), verified to raise `TypeError` on an unknown field.

Two real defects were found that the phase's own test suite cannot catch (both hidden behind the signature-only facade tests / an under-specific integration assertion): a dead sensitive-volume parameter in `run_microdosimetry`, and a materially inconsistent config-forwarding contract across the 7 non-run_cv/run_field facades that silently breaks `ParametricSweep` for several DeviceConfig fields and for `run_transient` entirely. Neither causes a crash, so both are classified WARNING rather than BLOCKER, but they directly undermine the phase's stated goal ("self-consistent... whole API," "uniform config-forwarding").

## Warnings

### WR-01: `sv_width_um` in `run_microdosimetry` is a dead parameter — always silently falls back to the slab approximation

**File:** `petringa/api/simulation.py:810-863` (specifically line 862)
**Issue:** `run_microdosimetry` calls:

```python
l_bar = mean_chord_length(sv_thickness_um, sv_width_um=sv_width_um)
```

but never passes `sv_depth_um`. The core function's branch selection is:

```python
if sv_width_um is not None and sv_depth_um is not None:
    ... # 4V/S 3D rectangular-parallelepiped formula
else:
    l_bar = 2.0 * t   # slab approximation
```

(`petringa/core/microdosimetry.py:59-75`). Since `sv_depth_um` is never supplied by the facade, the 3D branch is **never** reachable through this public API — `sv_width_um` has zero effect on the returned spectrum no matter what value the caller passes. A caller who sets `sv_width_um=50.0` vs `sv_width_um=5000.0` gets an identical `l_bar`, `y_F`, `y_D`, and spectrum. The facade's own docstring (`Returns` section, and the `sv_width_um` parameter doc) implies the width parameter is meaningful microdosimetric-geometry input, which is misleading given the actual code path.
`tests/test_api_microdosimetry.py:29-34` calls the facade with `sv_width_um=150` but only asserts finiteness and `y_D >= y_F`, so this masking is invisible to the test suite — the test would pass identically with `sv_width_um` deleted from the call entirely.
**Fix:** Either (a) add an `sv_depth_um` parameter to `run_microdosimetry` and forward it so the 3D branch is reachable, or (b) if the slab approximation is intentionally the only supported geometry for this facade, drop `sv_width_um` from the signature (or document explicitly in the docstring that it is currently unused/reserved) to avoid silently misleading callers:

```python
def run_microdosimetry(
    config: DeviceConfig,
    mc_csv_path: str,
    sv_thickness_um: float = 10.0,
    sv_width_um: float = 150.0,
    sv_depth_um: float | None = None,
) -> SimResult:
    ...
    l_bar = mean_chord_length(
        sv_thickness_um, sv_width_um=sv_width_um, sv_depth_um=sv_depth_um
    )
```

### WR-02: Config-forwarding is inconsistent across the 7 new facades, silently breaking `ParametricSweep` for several `DeviceConfig` fields

**File:** `petringa/api/simulation.py:414-807` (all 6 Plan-02 facades), `petringa/api/sweep.py:70-81`
**Issue:** The phase explicitly frames "uniform config-forwarding (D-01)" as an established convention (37-01-SUMMARY.md, 37-02-SUMMARY.md, 37-02-PLAN.md `must_haves`), and `ParametricSweep` is designed to sweep _any_ `DeviceConfig` field through _any_ facade `sim_fn`. In practice, forwarding is partial and facade-specific:

- Every facade forwards only `epi_thickness_um`, `N_D_junction`, `N_D_bulk`, `L_transition_um` (and `run_dark_current` additionally forwards `T`/`area_cm2`). None forward `config.N_A`, `config.N_D`, `config.doping_profile`, or `config.substrate_thickness_um`, even though `DeviceConfig` defines these fields and the underlying core constructors (`create_dd_device`/`create_sic_device`) accept them.
- `run_transient` (`petringa/api/simulation.py:737-807`) forwards **only** `epi_thickness_cm` — `transient_cce_vs_dose_rate` (`petringa/core/transient.py:545-645`) has no `N_D_junction`/`N_D_bulk`/`L_transition` parameters at all and hardcodes them internally (`N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1.0e-4` at `transient.py:597-599`). So `run_transient(DeviceConfig(N_D_junction=5e15))` silently ignores the doping override entirely, unlike every sibling facade.
- Concretely, this means `ParametricSweep(base_config=DeviceConfig(), param="N_D_junction", values=[...], sim_fn=run_transient)` produces **identical output for every value in the sweep** — no error, no warning, just a flat sweep curve that looks like a legitimate (if boring) physics result. This is exactly the kind of silent logic bug a reviewer must flag: it isn't caught by any test (facade contract tests are signature-only; `ParametricSweep` tests use a fake `sim_fn` that echoes any field back), and a downstream user has no signal that their sweep parameter was dropped.
- Similarly, sweeping `param="doping_profile"` (e.g., "uniform" vs "graded") via `ParametricSweep` over any of the 7 facades has zero effect on any of them, since none forward it — but nothing raises or warns; `dataclasses.replace` happily creates the config, and the facade quietly discards the field.
  **Fix:** At minimum, document per-facade which `DeviceConfig` fields are actually forwarded (the docstrings currently assert "config-forwarding (D-01)" uniformly without listing exceptions beyond `run_temperature_sweep`'s documented `T` exclusion). Preferably, either extend `run_transient`'s core function to accept doping kwargs and forward them like its siblings, or explicitly and loudly document run_transient's doping-forwarding gap next to the existing D-03 minimal-satisfier note (it currently reads as an oversight, not a stated limitation). Consider having `ParametricSweep` (or a debug helper) warn when `param` doesn't appear anywhere in the resulting `SimResult.metadata`/traceable forwarding path, though that is a larger design change.

### WR-03: `run_dark_current` never supplies a `device_name`, unlike every device-owning sibling facade

**File:** `petringa/api/simulation.py:692-702`
**Issue:** `run_cv`/`run_field` build devices via `build_device()`, which always generates a `uuid4`-based unique name (`petringa/api/device.py:78,96`). `run_dark_current` instead calls `create_dark_current_device(T=..., N_t=..., S_n=..., S_p=..., epi_thickness_cm=..., N_D_junction=..., N_D_bulk=..., L_transition=...)` with no `device_name=` kwarg, so it falls through to `create_dd_device` → `create_sic_device`'s hardcoded default `device_name="sic_diode"` (`petringa/core/device.py:43`). This is safe under the documented single-threaded, `reset_devsim_fully()`-at-entry / `finally`-delete usage pattern this phase relies on everywhere, but it is an inconsistency relative to the established pattern (`build_device`, `cce_vs_bias`, `cce_vs_fluence`, `transient_cce_vs_dose_rate` all mint unique names) and is a latent trap if `run_dark_current` is ever called concurrently with another facade in the same process, or if a future refactor removes the blanket `reset_devsim_fully()` call.
**Fix:** Pass an explicit unique `device_name` (e.g. `f"dark_current_{uuid.uuid4().hex[:8]}"`) to `create_dark_current_device`, matching the naming discipline used by every other device-owning code path in this phase and Phase 36.

## Info

### IN-01: `sv_width_um`'s docstring doesn't disclose the geometry-branch dependency

**File:** `petringa/api/simulation.py:846-847`
**Issue:** The parameter docstring for `sv_width_um` reads simply "Sensitive-volume width (um). Default 150.0." with no mention that it is only used together with a (currently unexposed) `sv_depth_um` to select the 3D chord-length formula. Related to WR-01 — flagged separately since it's a documentation-only gap that would remain even if WR-01 is fixed by exposing `sv_depth_um` (the docstring should then explain the two-parameter-or-slab-fallback contract).
**Fix:** Once WR-01 is resolved, document the interaction explicitly, e.g.: "Only used (together with `sv_depth_um`) for the 3D rectangular-parallelepiped chord length; if `sv_depth_um` is None, the slab approximation `l_bar = 2 * sv_thickness_um` is used and `sv_width_um` has no effect."

### IN-02: Five of nine public facades have zero behavioral test coverage (expected/accepted, noted per review brief)

**File:** `tests/test_api_facades.py`, `petringa/api/simulation.py:414-807`
**Issue:** `run_radiation_damage`, `run_dark_current`, `run_temperature_sweep`, `run_flash_recombination`, and `run_transient` are covered only by signature-introspection tests (`callable()` + first-param-name checks); no test ever executes their devsim code paths or validates their `SimResult` output against a physics gate. This is explicitly called out in 37-02-SUMMARY.md as a deliberate choice to avoid devsim process-exhaustion, and the review brief itself says to treat this as expected/accepted rather than novel. Noted here only because it is the direct cause of WR-01/WR-02 being invisible to CI: any facade-level key/column mismatch, wrong-axis error, or dropped-forwarding bug in these five facades would only surface via manual execution or a future integration test, not via `pytest`.
**Fix:** No action required for this phase per the accepted test-coverage tradeoff; flagged for downstream phases planning end-to-end facade coverage to prioritize these five, especially `run_transient` and `run_dark_current` given WR-02/WR-03 above.

---

_Reviewed: 2026-07-09_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
