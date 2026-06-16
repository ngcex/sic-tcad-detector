# Phase 26: Graded Doping 2D Calibration — Research

**Researched:** 2026-05-17
**Domain:** devsim 2D Poisson + drift-diffusion convergence at high reverse bias; graded-doping calibration; devsim global-state hygiene
**Confidence:** HIGH

## Summary

Phase 26 fixes the known v3.0 tech debt that `device2d.py` cannot converge above approximately −15 V reverse bias, blocking the full clinical bias range (−50 V) needed by Phases 29 (noise) and 30 (build-up). The phase touches **existing files only** — no new modules. Three families of work: (1) re-calibrate `_N_D_JUNCTION_DEFAULT / _N_D_BULK_DEFAULT / _L_TRANSITION_DEFAULT` for the 2D mesh by adding a `calibrate_graded_doping_2d` function mirroring `device.py:468`, (2) extend `reset_devsim()` (currently inlined at `optimization.py:153`) to clear all devsim global state that leaks across 2D devices including cylindrical-axis parameters from `alternative_structures.py`, and (3) freeze and protect the existing 22 v3.0 notebooks with a baseline-frozen regression sweep.

**The calibration target is the validated 1D C-V curve at the device center, NOT a separate 2D experimental dataset.** Petringa measured C-V on a real device — the 1D simulator was already calibrated to this data (W(0V)=1.7 µm, W(−10V)=9.5 µm, W(−30V)=9.73 µm). Success criterion #2 (R² ≥ 0.99 between 2D-center and 1D C-V) means the 2D solver must reproduce its own validated 1D twin — a self-consistency target, not external data.

**Primary recommendation:** Sequence the work as ROOT-CAUSE DIAGNOSIS → CALIBRATION → CLEANUP → REGRESSION. Do not jump straight to refitting `{N_D_junction, N_D_bulk, L_transition}` for 2D — first confirm whether the failure mode is profile-driven (H1), solver/mesh/BC differential between 1D and 2D (H2), or both (H3). The same `{2.9e15, 8.5e13, 1e-4}` profile converges to −30 V in 1D today; that empirical fact is diagnostic.

## Phase Requirements

| ID      | Description                                                                                                                                                                                                    | Research Support                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| ------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| CONS-01 | User can run 2D device simulations at reverse biases beyond −15 V without solver divergence, using a re-calibrated graded epi doping profile in `device2d.py` that matches C-V data across the full bias range | Existing 1D infrastructure (`device.py:calibrate_graded_doping` lines 468–637) provides a verified pattern to port to 2D using `create_2d_dd_device` (charge_collection_2d.py:116). PITFALLS.md P27 (graded doping at node), P03/P30 (global state leaks), P20 (device-name collision), P24 (regression baseline policy) define the cleanup required to make convergence stable. The `cv_sweep` function (`cv_analysis.py:121`) already implements adaptive solve-with-fallback at high bias and can be reused as the 2D verification harness. |

## Standard Stack

No new packages. v3.0 stack is sufficient for this phase.

### Core (unchanged)

| Library        | Version       | Purpose                                | Why Standard                                                                                                            |
| -------------- | ------------- | -------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| devsim         | ≥2.10.0       | TCAD finite-volume Poisson + DD solver | [VERIFIED: existing code, calibrated through v3.0] Already the project's TCAD engine; 2D mesh API in use since Phase 19 |
| scipy.optimize | (scipy ≥1.11) | Nelder-Mead minimiser for calibration  | [VERIFIED: device.py:617] Already used by 1D `calibrate_graded_doping`; port pattern unchanged                          |
| numpy          | ≥1.24         | Array math                             | [VERIFIED]                                                                                                              |
| pytest         | ≥7.0          | Regression test harness                | [VERIFIED: existing tests/test_device2d.py]                                                                             |

### Files touched (existing only)

| File                            | Role                                                                         | Change Type                                                                                                                                              |
| ------------------------------- | ---------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/device2d.py`               | 2D mesh + doping setup; holds `_N_D_*_DEFAULT` constants                     | EDIT — update calibrated defaults; ensure graded-doping helper survives reverse-bias sweep                                                               |
| `src/cv_analysis.py`            | C-V sweep harness                                                            | REUSE — already supports 1D; assess whether `cv_sweep` works on a 2D device unchanged or needs a thin wrapper to pick the center column for W extraction |
| `src/charge_collection_2d.py`   | `create_2d_dd_device(half_width_um, V_bias, ...)`                            | REUSE for calibration trial loop; no edits expected                                                                                                      |
| `src/optimization.py`           | inline `reset_devsim()` call at line 153                                     | EXTRACT — promote to a shared utility (e.g. `src/devsim_reset.py` or top-level in `device2d.py`) that enumerates devices and clears global parameters    |
| `src/alternative_structures.py` | uses `raxis_zero`, `raxis_variable`, `cylindrical_*` models — leaks globally | NO EDIT — but the extended reset must clean up after it                                                                                                  |
| `tests/test_device2d.py`        | existing 4 test classes                                                      | EXTEND — add reverse-bias convergence tests, 2D-vs-1D C-V regression, reset-state regression                                                             |

### Alternatives Considered

| Instead of                       | Could Use                                         | Why Not                                                                                                                                                                                             |
| -------------------------------- | ------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Re-fitting graded profile in 2D  | Use 1D-fitted values as-is                        | The 1D profile DOES converge in 1D up to −30 V at least; assuming the 2D failure is purely profile-driven without testing the mesh/BC differential first risks treating a symptom rather than cause |
| New `src/devsim_reset.py` module | Keep `reset_devsim()` inline in `optimization.py` | Phase 31 (anisotropic) explicitly says it must also extend reset; this code will be called from multiple modules — promoting it now removes a forced refactor in Phase 31                           |
| Numerical-jacobian calibration   | Nelder-Mead simplex                               | Nelder-Mead is what 1D uses (`device.py:617`) and is robust for the ~3-parameter low-dim fit; no need to introduce gradient methods                                                                 |

**Installation:** None — phase uses existing stack only.

## Architecture Patterns

### Current Module Layout (no changes)

```
src/
├── device2d.py              # 2D mesh creation, graded doping (THIS PHASE EDITS)
├── charge_collection_2d.py  # create_2d_dd_device wrapper, ramps to V_bias
├── cv_analysis.py           # cv_sweep — reused for 2D verification
├── optimization.py          # holds the inline reset_devsim() pattern (THIS PHASE EXTRACTS)
├── alternative_structures.py # source of cylindrical-axis global-state leak
├── device.py                # 1D — frozen, source of calibration pattern to port
└── poisson.py, drift_diffusion.py  # dimension-agnostic, no changes
```

### Pattern 1: Calibration Function (port 1D → 2D)

**What:** Mirror `device.py:calibrate_graded_doping` (lines 468–637) for 2D, using `create_2d_dd_device` for trial builds and extracting W at the center column.
**When to use:** This phase only — calibration is one-shot; results bake into `_N_D_*_DEFAULT` constants in `device2d.py`.
**Reference (the 1D pattern to port):**

```python
# Source: src/device.py:534–615 (1D calibrate_graded_doping objective fn)
def objective(params_vec):
    N_D_j, N_D_b, L_t = params_vec
    # bounds penalty ...
    try:
        device_info = create_dd_device(
            device_name=trial_name,
            doping_profile="graded",
            N_D_junction=N_D_j, N_D_bulk=N_D_b, L_transition=L_t,
        )
        # ramp cathode in 0.5V steps to each target, extract W
        # cost = sum of squared relative errors vs W_exp
    finally:
        devsim.delete_device(device=trial_name)
        devsim.delete_mesh(mesh=mesh_name)
    return cost
```

**2D adaptation requirements** [VERIFIED: cv_analysis.py:121, charge_collection_2d.py:116]:

- Use `create_2d_dd_device(half_width_um=..., V_bias=0.0, doping_profile="graded", ...)` to instantiate.
- Extract W at center column: filter nodes by `x < 1e-6`, sort by `y`, find depletion edge from carrier profile (same logic as `extract_depletion_width_numerical` in `poisson.py`).
- Use `cv_sweep` (`cv_analysis.py:121`) or replicate its solve-with-fallback pattern — it already handles the high-bias convergence path:

```python
# Source: src/cv_analysis.py:188–207 (proven fallback pattern)
try:
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-10, maximum_iterations=40)
except devsim.error:
    devsim.solve(type="dc", absolute_error=1e12, relative_error=1e-8, maximum_iterations=100)
```

### Pattern 2: Extended reset_devsim() (PITFALLS P03, P20, P30)

**What:** A reset utility that nukes all devsim global state, not just devices.
**When to use:** Between alt-structure runs and planar runs in the same Python session; in pytest fixtures; in calibration trial loops.

```python
# Targets per PITFALLS P03 + P30 — verbatim global names that leak
def reset_devsim_fully():
    # 1. Delete all devices (enumerate, not hardcode — P20)
    for dev in list(devsim.get_device_list()):
        try:
            devsim.delete_device(device=dev)
        except Exception:
            pass
    # 2. Clear cylindrical-axis globals (set by alternative_structures.py:579–584)
    for name in ("raxis_zero", "raxis_variable",
                 "node_volume_model", "edge_couple_model",
                 "element_edge_couple_model", "surface_area_model"):
        try:
            devsim.set_parameter(name=name, value="")
        except Exception:
            pass
    # 3. Save+restore solver settings (per optimization.py:96–104 existing pattern)
    devsim.reset_devsim()
```

**Justification:** `alternative_structures.py:573–584` explicitly sets `raxis_zero`, `raxis_variable`, and replaces `node_volume_model`/`edge_couple_model`/`surface_area_model` with `cylindrical_*` variants. devsim's `reset_devsim()` alone does NOT unset these — they persist as process globals into the next planar run. This is PITFALLS P03 (severity HIGH).

### Pattern 3: Baseline-Frozen Regression Sweep

**What:** A single pytest that loads frozen reference values (CCE, dark current, W(V), key y-spectrum statistics) for the 22 v3.0 notebooks and asserts the new code reproduces them within tolerance.
**When to use:** Once per phase (and once again in Phase 34 milestone audit).

```python
# Pattern (Wave 0 must create the baseline JSON file first)
@pytest.mark.regression
def test_v3_baseline_preserved():
    baseline = json.load(open("tests/baselines/v3_frozen.json"))
    for key, expected in baseline.items():
        # rebuild device, extract metric, compare
        actual = rebuild_and_extract(key)
        assert math.isclose(actual, expected, rel_tol=1e-3), key
```

The "20 notebooks" count in the phase description is approximate — the actual notebook directory contains 22 files (`01_*` through `20_*` plus `03_executed.ipynb` and others). [VERIFIED: `ls notebooks/*.ipynb | wc -l` → 22] The planner should count notebooks that have a numbered scientific deliverable and exclude purely re-executed copies.

### Anti-Patterns to Avoid

- **Hardcoding device names** [PITFALLS P20]: never use literal `"device2d"` — always pass `device_name=f"cal_2d_{uuid.uuid4().hex[:8]}"`. The existing `test_device2d.py` `_unique_name()` factory is a good template.
- **Evaluating doping in Python and injecting via `set_node_values`** [PITFALLS P27]: stay with the `devsim.node_model(equation="...")` pattern already in `set_graded_doping_2d` — devsim handles per-node sampling correctly. Verified in `device2d.py:126–143`.
- **Adding a `src/calibration.py` module**: tempting but unnecessary — the calibration function lives next to `set_graded_doping_2d` in `device2d.py` and is run once, then the result bakes into the default constants.
- **Refitting with a different objective than 1D**: keep the cost function `sum((W_sim - W_exp)/W_exp)²` and the bounds (`N_D_j ∈ [1e14, 1e16]`, `N_D_b ∈ [1e12, 1e15]`, `L_t ∈ [0.5e-4, 5e-4]`, `N_D_b < N_D_j`) identical to 1D. Asymmetric bounds would make the 2D result physically inconsistent with the 1D twin.

## Don't Hand-Roll

| Problem                                  | Don't Build                                        | Use Instead                                                                                                        | Why                                                                                                        |
| ---------------------------------------- | -------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------- |
| Bias ramping to high reverse voltage     | Custom step-doubling logic                         | `ramp_bias()` from `drift_diffusion.py:245` + fallback pattern from `cv_analysis.py:188`                           | Already battle-tested in v3.0; handles convergence failures gracefully                                     |
| Depletion-width extraction at the center | New numerical method                               | `extract_depletion_width_numerical(device_info)` from `poisson.py` filtered to `x < 1e-6`                          | The 1D version exists and is already validated against C-V data                                            |
| Optimisation loop                        | scipy.curve_fit / lmfit                            | `scipy.optimize.minimize(method="Nelder-Mead")`                                                                    | Matches 1D exactly (`device.py:617`); no new dependency; well-behaved on 3-parameter low-dim problem       |
| Device-name uniqueness                   | `uuid.uuid4().hex[:8]` rolled fresh in each module | Existing `_unique_name()` pattern in `tests/test_device2d.py` or `uuid` directly                                   | Already in use across the codebase (16 sites use `uuid` for device names)                                  |
| C-V sweep over reverse-bias array        | New 2D sweep                                       | `cv_sweep(device_info, V_range, area=...)` from `cv_analysis.py:121`                                               | Already iterates voltages with adaptive step + retry; agnostic to 1D/2D since it operates on `device_info` |
| Calibration loop structure               | New optimisation harness                           | Port `device.py:calibrate_graded_doping` (lines 468–637) verbatim, swap `create_dd_device` → `create_2d_dd_device` | The 1D function is the proven pattern; same Nelder-Mead, same bounds, same cost, same cleanup              |

**Key insight:** This phase is calibration + cleanup, not new architecture. Every helper needed already exists. The work is **wiring existing pieces in a new sequence** and tightening the global-state contract.

## Common Pitfalls

(See `.planning/research/PITFALLS.md` for full v4.0 catalog. The phase-specific items below cite that document; one-line summaries follow.)

### Pitfall P27 — Graded doping at node, not mesh (HIGH)

**What goes wrong:** Evaluating the doping function on mesh-line coordinates instead of devsim node coordinates produces spurious gradients at refinement transitions.
**Status in our code:** Already correct in `device2d.py:126–143` — `set_graded_doping_2d` uses `devsim.node_model(equation="...")` which is the right pattern.
**Action this phase:** Verify the calibrated parameters still produce a smooth donor profile across the y-axis (no kinks) — add a `test_graded_doping_smoothness` to `test_device2d.py`.

### Pitfall P03 + P30 — Cylindrical-axis globals leak (HIGH)

**What goes wrong:** `alternative_structures.py:573–584` calls `devsim.set_parameter(name="raxis_zero", value=0.0)` and replaces `node_volume_model`/`edge_couple_model`/`surface_area_model` with `cylindrical_*` variants. These are devsim **process globals**, not per-device. After alt-structure execution, any subsequent planar 2D device silently inherits cylindrical assembly weights and produces wrong Poisson results.
**How to detect:** Success criterion #3 of this phase — "regression test that runs alt-structures then planar in sequence" is exactly the canary for this bug.
**Fix this phase:** Extend `reset_devsim()` per Pattern 2 above. Required globals to clear: `raxis_zero`, `raxis_variable`, `node_volume_model`, `edge_couple_model`, `element_edge_couple_model`, `surface_area_model`.

### Pitfall P20 — Hardcoded device names collide (HIGH)

**What goes wrong:** Notebooks re-run in the same kernel hit "device already exists". The calibration loop, which builds 30–100 trial devices, is especially vulnerable.
**Fix this phase:** Always use `device_name=f"cal_2d_{uuid.uuid4().hex[:8]}"` in the new `calibrate_graded_doping_2d` (matches `device.py:531`). Ensure existing notebooks 15, 17, 19 still pass — they currently use `_unique_name()` or literal strings; do not change their interface.

### Pitfall P24 — Regression baseline policy (MEDIUM)

**What goes wrong:** Notebook 14 (validation) snapshots reference values. Changing the calibrated `_N_D_*_DEFAULT` constants in `device2d.py` could shift 2D CCE / dark-current / W(V) outputs measurably, breaking the validation notebook.
**Decision needed:** Per the v4.0 SUMMARY.md decision D4: **"Freeze v3.0 baselines; `anisotropic=False` default preserves all existing notebooks"** — meaning the calibrated 2D defaults must shift 2D outputs **as little as possible** at biases ≤ −10 V (where v3.0 already converged). The fit objective should weight low-bias agreement (preserve v3.0) AND high-bias convergence (the new requirement) — not just high-bias.
**Implementation:** The cost function should include W(0V), W(−10V), W(−30V) AND verify convergence (not W value) at W(−50V). This preserves v3.0 behavior in the already-validated range while extending the operating envelope.

### Pitfall (new) — Air buffer interaction with high-field SCR (MEDIUM)

**What goes wrong:** `device2d.py:216` defines `air_buffer = 1e-8` cm — a 0.1 nm thick "air" region whose only purpose is to provide nodes for devsim contact detection. At high reverse bias (−50 V) the depletion region extends to ~10 µm; the field at the cathode contact may interact with the air-buffer interface in ways that don't occur in 1D (which has no buffer regions).
**Diagnostic test (Wave 0):** Run identical 1D and 2D devices with the same graded profile to −30 V and compare the cathode-side carrier profile. If 2D has anomalies near `y ≈ total_depth` that 1D does not, the air buffer is suspect.
**Source confidence:** [ASSUMED] — flagged because no existing test isolates this; the planner should add a diagnostic step.

## Runtime State Inventory

> Phase 26 is a calibration + cleanup phase, not a rename/refactor. The categories below answer the question "what runtime state must change beyond source files?"

| Category            | Items Found                                                                                                                                                                                                                                        | Action Required                                                         |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------- |
| Stored data         | **None** — no databases, no persisted state. Calibration writes results into `_N_D_*_DEFAULT` constants in `device2d.py` (source code) only.                                                                                                       | None                                                                    |
| Live service config | **None** — no external services (n8n, Datadog, Cloudflare, etc.) in this project.                                                                                                                                                                  | None                                                                    |
| OS-registered state | **None** — no Task Scheduler entries, no launchd plists, no systemd units.                                                                                                                                                                         | None                                                                    |
| Secrets/env vars    | **None** — no .env, no SOPS, no credentials. Pure scientific code.                                                                                                                                                                                 | None                                                                    |
| Build artifacts     | `__pycache__/` directories per PITFALLS P29; project memory `feedback_uv_venv.md` flags uv-managed venv. After changing `_N_D_*_DEFAULT` constants, stale `.pyc` files in `src/__pycache__/` could shadow new sources for non-uv invocation paths. | `find . -name "__pycache__" -exec rm -rf {} +` before re-running pytest |

## Validation Architecture

> `workflow.nyquist_validation` is not set in `.planning/config.json`; treated as enabled per default policy.

### Test Framework

| Property           | Value                                                                                                                                                                                                             |
| ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Framework          | pytest ≥7.0 [VERIFIED: existing 22 test files]                                                                                                                                                                    |
| Config file        | None at repo root; pytest discovers via convention (`tests/test_*.py`) [VERIFIED]                                                                                                                                 |
| Quick run command  | `uv run pytest tests/test_device2d.py -x`                                                                                                                                                                         |
| Full suite command | `uv run pytest -x --tb=short`                                                                                                                                                                                     |
| Slow tests         | None marked `@pytest.mark.slow`; full DD calibration loop in this phase will be slow (~30–60 s per trial × ~50 trials) — recommend marking the calibration test `@pytest.mark.slow` and excluding from quick runs |

### Phase Requirements → Test Map

| Req ID       | Behavior                                                                 | Test Type   | Automated Command                                                   | File Exists?                    |
| ------------ | ------------------------------------------------------------------------ | ----------- | ------------------------------------------------------------------- | ------------------------------- |
| CONS-01 (#1) | 2D device converges at V = −15, −30, −50 V on both 100 µm and 300 µm SVs | integration | `pytest tests/test_device2d.py::TestReverseBiasConvergence -x`      | ❌ Wave 0                       |
| CONS-01 (#2) | 2D C-V at center matches 1D C-V (R² ≥ 0.99)                              | integration | `pytest tests/test_device2d.py::test_2d_vs_1d_cv_centerline -x`     | ❌ Wave 0                       |
| CONS-01 (#3) | reset_devsim clears alt-structure global state                           | regression  | `pytest tests/test_device2d.py::test_reset_after_alt_structures -x` | ❌ Wave 0                       |
| CONS-01 (#4) | 22 v3.0 notebook outputs unchanged (within tolerance)                    | regression  | `pytest tests/test_v3_baseline_regression.py -x`                    | ❌ Wave 0 (needs baseline JSON) |

### Sampling Rate

- **Per task commit:** `uv run pytest tests/test_device2d.py -x` (~15 s without calibration; ~5 min with `--run-slow`)
- **Per wave merge:** `uv run pytest -x --tb=short` (~10–20 min — DD simulations dominate)
- **Phase gate:** Full suite green + the 4 new regression tests above pass + manual visual inspection of 2D-vs-1D C-V overlay plot

### Wave 0 Gaps

- [ ] `tests/test_device2d.py` — add `TestReverseBiasConvergence` class (3 cases: −15, −30, −50 V × 2 SV sizes = 6 cases)
- [ ] `tests/test_device2d.py` — add `test_2d_vs_1d_cv_centerline` (compares `cv_sweep` output on 2D-center vs 1D for V ∈ [0, −50] in 5 V steps)
- [ ] `tests/test_device2d.py` — add `test_reset_after_alt_structures` (run a mesa or 3D-electrode structure, run reset, run planar 2D, assert correctness)
- [ ] `tests/test_v3_baseline_regression.py` — new file
- [ ] `tests/baselines/v3_frozen.json` — Wave 0 must extract and freeze v3.0 reference outputs (CCE at 50 V for 100/300 µm SVs; dark current at −10 V; W(V) at 5 voltages) BEFORE any calibration changes
- [ ] `pytest.ini` or `conftest.py` — register `@pytest.mark.slow` marker if not already (verify in Wave 0)
- [ ] Optional: a single integration smoke that imports `device2d`, `charge_collection_2d`, `alternative_structures`, and `cv_analysis` and runs a 30-second end-to-end happy path — catches import-cycle issues introduced by extracting `reset_devsim_fully`

## Code Examples

### Pattern: 2D calibration objective (port from 1D)

```python
# Source: adapt from src/device.py:534–615
# Key swap: create_dd_device → create_2d_dd_device; extract W at center column

import uuid
import numpy as np
import devsim
import scipy.optimize

from src.charge_collection_2d import create_2d_dd_device

def calibrate_graded_doping_2d(
    target_W_data=None,
    half_width_um=50.0,
    epi_thickness_cm=10e-4,
    N_A=1e19,
    T=300,
    x0=None,
    maxiter=80,
):
    """2D analog of src/device.py:calibrate_graded_doping.

    Target voltages include extended reverse-bias range (−50 V) to address CONS-01.
    Cost function blends low-bias agreement (preserves v3.0) with high-bias
    convergence (new requirement).
    """
    if target_W_data is None:
        # Petringa 1D-calibrated targets + extended reverse-bias points.
        # The −50 V target is "must converge"; for the W-value we accept the 1D
        # twin's W value as ground truth at center column (success criterion #2).
        target_W_data = {0.0: 1.7e-4, -10.0: 9.5e-4, -30.0: 9.73e-4, -50.0: None}

    if x0 is None:
        x0 = [2.9e15, 8.5e13, 1e-4]  # v3.0 1D-calibrated values as starting point

    run_id = uuid.uuid4().hex[:8]
    counter = [0]

    def objective(p):
        N_D_j, N_D_b, L_t = p
        if not (1e14 <= N_D_j <= 1e16 and 1e12 <= N_D_b <= 1e15
                and 0.5e-4 <= L_t <= 5e-4 and N_D_b < N_D_j):
            return 1e6

        name = f"cal2d_{run_id}_{counter[0]}"
        counter[0] += 1
        cost = 0.0
        try:
            # Build device at 0 V then ramp via cv_sweep (handles fallback tolerances)
            dev = create_2d_dd_device(
                device_name=name,
                half_width_um=half_width_um,
                V_bias=0.0,                          # ramp via cv_sweep below
                doping_profile="graded",
                N_D_junction=N_D_j, N_D_bulk=N_D_b, L_transition=L_t,
                epi_thickness_cm=epi_thickness_cm, N_A=N_A, T=T,
            )
            from src.cv_analysis import cv_sweep
            voltages = sorted(target_W_data.keys(), reverse=True)  # 0, -10, -30, -50
            result = cv_sweep(dev, V_range=voltages)
            # Extract W at center column only (success criterion #2)
            W_2d_center = extract_W_at_center(dev, voltages)  # helper from poisson.py pattern
            for v, W_exp in target_W_data.items():
                i = list(result["voltages"]).index(v)
                W_sim = W_2d_center[i]
                if W_exp is None:
                    # −50 V: only convergence requirement, no W target
                    if W_sim <= 0 or not np.isfinite(W_sim):
                        cost += 1e3
                else:
                    cost += ((W_sim - W_exp) / W_exp) ** 2
        except Exception as e:
            cost = 1e6
        finally:
            try:
                devsim.delete_device(device=name)
            except Exception:
                pass
        return cost

    result = scipy.optimize.minimize(
        objective, x0, method="Nelder-Mead",
        options={"maxiter": maxiter, "xatol": 1e-10, "fatol": 1e-6},
    )
    # write result.x into src/device2d.py _N_D_*_DEFAULT constants (manual step)
    return {"N_D_junction": result.x[0], "N_D_bulk": result.x[1],
            "L_transition": result.x[2], "final_cost": result.fun}
```

### Pattern: Extended reset utility

```python
# Add to src/device2d.py OR new src/devsim_reset.py (planner decides)
# Source: derived from PITFALLS P03/P20/P30 + existing optimization.py:96–158 pattern

import devsim

# Globals that leak from alternative_structures.py:573–584
_CYLINDRICAL_GLOBALS = (
    "raxis_zero", "raxis_variable",
    "node_volume_model", "edge_couple_model",
    "element_edge_couple_model", "surface_area_model",
)

def reset_devsim_fully(preserve_solver=True):
    """Clear all devsim global state including cylindrical-axis parameters.

    Replaces the inline `devsim.reset_devsim()` call in optimization.py:153.
    Must be used after any cylindrical 2D run (alternative_structures.py 3D electrode)
    before instantiating a planar 2D device, or the planar mesh assembly will be wrong.
    """
    saved_solver = None
    saved_callback = None
    if preserve_solver:
        try:
            saved_solver = devsim.get_parameter(name="direct_solver")
        except devsim.error:
            saved_solver = "custom"
        try:
            saved_callback = devsim.get_parameter(name="solver_callback")
        except devsim.error:
            saved_callback = None

    # Enumerate-don't-hardcode (P20)
    for dev in list(devsim.get_device_list()):
        try:
            devsim.delete_device(device=dev)
        except Exception:
            pass

    # Clear cylindrical-axis + assembly-model globals (P03, P30)
    for name in _CYLINDRICAL_GLOBALS:
        try:
            devsim.set_parameter(name=name, value="")
        except Exception:
            pass

    devsim.reset_devsim()

    if preserve_solver:
        if saved_solver is not None:
            devsim.set_parameter(name="direct_solver", value=saved_solver)
        if saved_callback is not None:
            devsim.set_parameter(name="solver_callback", value=saved_callback)
```

### Pattern: Regression-canary test for cylindrical state leak (success criterion #3)

```python
# tests/test_device2d.py (NEW test)
def test_reset_after_alt_structures():
    """alt-structures → reset → planar should match alt → planar (no reset = wrong)."""
    from src.alternative_structures import create_3d_electrode_structure  # cylindrical
    from src.charge_collection_2d import create_2d_dd_device
    from src.device2d import reset_devsim_fully

    # Reference: planar in a fresh session
    reset_devsim_fully()
    dev_ref = create_2d_dd_device(device_name="planar_ref", V_bias=10.0)
    I_ref = _extract_dark_current(dev_ref)
    devsim.delete_device(device="planar_ref")

    # Now: alt-structures FIRST, then reset, then planar
    reset_devsim_fully()
    dev_alt = create_3d_electrode_structure(device_name="alt_test")
    # ... run a brief equilibrium solve ...
    reset_devsim_fully()
    dev_planar = create_2d_dd_device(device_name="planar_after_reset", V_bias=10.0)
    I_after = _extract_dark_current(dev_planar)

    assert abs(I_after - I_ref) / abs(I_ref) < 1e-3, (
        "reset_devsim_fully failed to clear cylindrical globals"
    )
```

## State of the Art

| Old Approach (v3.0)                                                                | Current Approach (Phase 26)                                               | When Changed        | Impact                                                                                  |
| ---------------------------------------------------------------------------------- | ------------------------------------------------------------------------- | ------------------- | --------------------------------------------------------------------------------------- |
| 2D uses 1D-calibrated `{2.9e15, 8.5e13, 1e-4}` graded profile, fails > −15 V       | 2D-specific calibration (or proven-equivalent profile after H2 diagnosis) | This phase          | Unblocks Phases 29 (noise) and 30 (build-up) which need stable field at clinical biases |
| `reset_devsim()` is inline in `optimization.py` and clears devices only            | Extended utility clears cylindrical-axis globals + assembly model names   | This phase          | Eliminates a known silent-corruption path for alt-structure → planar notebook sequences |
| Notebook 14 validation snapshots v3.0 outputs but no programmatic regression sweep | Frozen JSON of v3.0 outputs + pytest regression                           | This phase (Wave 0) | Catches downstream physics drift in Phases 27–34                                        |

**Deprecated / outdated:**

- Inlined `try: devsim.reset_devsim() ... except: pass` blocks scattered through `optimization.py` and notebooks — replace with `reset_devsim_fully()` import once available.

## Assumptions Log

| #   | Claim                                                                                                                                                                                                                                              | Section                 | Risk if Wrong                                                                                                                                                                                                      |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| A1  | The 2D solver divergence above −15 V is caused EITHER by the 2D mesh interacting with the existing graded profile OR by the profile itself being marginal — **NOT** by a devsim version-specific solver bug                                        | Summary, Pitfalls       | If devsim itself is the issue, the calibration approach won't fix it. Mitigation: Wave 0 diagnostic step that runs a 1D device with identical settings ramped to −50 V; if 1D also fails, escalate                 |
| A2  | The Petringa C-V dataset already encoded in `device.py:519` (W(0V)=1.7 µm, W(−10V)=9.5 µm, W(−30V)=9.73 µm) is authoritative for the 2D calibration target. No separate 2D experimental C-V exists.                                                | Summary                 | If experimental 2D C-V data does exist, the target dataset shifts. Mitigation: explicitly ask the user in `/gsd-discuss-phase` whether external 2D C-V data is available                                           |
| A3  | The "20 v3.0 notebooks" baseline-freeze target in success criterion #4 should be interpreted as "all numbered scientific notebooks 01–20" — the actual notebooks/ directory contains 22 files including `03_executed.ipynb` and one or two re-runs | Validation Architecture | If the intent was literally 20 distinct deliverables, baseline selection differs. Mitigation: planner verifies notebook count and identifies which are "frozen scientific outputs"                                 |
| A4  | The −50 V target is achievable with a graded epi profile (not requiring a guard-ring or mesa termination) — this is consistent with PROJECT.md framing 100×100×10 µm SV as the planar target                                                       | Code Examples           | If breakdown at −50 V requires structural changes (guard ring), the phase scope is wrong. Mitigation: discuss-phase should confirm with PI that −50 V on a planar diode without guard ring is the realistic target |
| A5  | Air-buffer regions (`air_buffer = 1e-8 cm`, `device2d.py:216`) do not contribute to the −15 V divergence — the buffer is far from the SCR                                                                                                          | Pitfalls                | If the air buffer is implicated, the fix differs (mesh redesign, not profile refit). Mitigation: Wave 0 diagnostic comparing carrier profiles near `y = total_depth` in 1D vs 2D                                   |
| A6  | `cv_sweep` from `cv_analysis.py:121` works on a 2D `device_info` dict unchanged — the function operates on `device_info["device_name"]` and uses contact bias, both dimension-agnostic                                                             | Architecture Patterns   | If `cv_sweep` has 1D assumptions in `extract_depletion_width_numerical`, a thin wrapper is needed. Mitigation: planner verifies by inspection or smoke test in Wave 0                                              |
| A7  | Calibration result will fit within the existing parameter bounds (`N_D_j ∈ [1e14, 1e16]`, `L_t ∈ [0.5e-4, 5e-4]`)                                                                                                                                  | Code Examples           | If 2D optimum lies outside, bounds need widening — may indicate a non-physical optimum and require revisiting H2/H3                                                                                                |

## Open Questions (RESOLVED)

All five open questions below have been resolved during planning. The RESOLVED: line in
each item gives the answer that downstream plans (Plans 01–04) operate against.

1. **Hypothesis confirmation: H1 (profile) vs H2 (mesh/BC) vs H3 (both).**
   - What we know: 1D `{2.9e15, 8.5e13, 1e-4}` converges to at least −30 V; 2D same profile fails > −15 V. 1D and 2D share identical y-axis mesh spacing at the junction (both use `ps=1e-7` at junction, `ps=1e-6` 5 µm before, `ps=5e-6` in epi mid).
   - What's unclear: which difference (air buffer regions, lateral x-mesh interaction, symmetry BC at `x=0`) causes the divergence.
   - Recommendation: Wave 0 diagnostic — port the failing 2D bias point to 1D with matched parameters; if 1D succeeds, the 2D-specific differential (air buffer, x-mesh) is the culprit. Plan accordingly.
   - **RESOLVED:** Defer the H1/H2/H3 classification to runtime. Plan 26-01 Task 1 (`scripts/diagnose_1d_2d_parity.py`) empirically classifies the hypothesis by running 1D and 2D devices with the SAME graded profile to −50 V and writing the result to `26-DIAGNOSIS.md`. Plan 26-03 (`scripts/run_calibration_2d.py`) reads the YAML `hypothesis:` key and gates the calibration: H1 or H3 → proceed with Nelder-Mead; H2 → abort with exit code 2 and a message naming the suspect (air buffer, x-mesh, x-symmetry BC) so the user can re-scope before sinking 30–60 min into a doomed fit. The plan does not pre-commit to an answer — the diagnostic IS the answer.

2. **Cost function weighting between low-bias preservation and high-bias convergence.**
   - What we know: Phase success criterion #4 demands ZERO regression on 22 notebooks (low-bias must be preserved); criterion #1 demands −50 V convergence (new).
   - What's unclear: how to weight `W(0V), W(−10V), W(−30V)` residuals vs the new `W(−50V)` requirement.
   - Recommendation: Start with equal weights (1/N) for known-W targets plus a hard convergence penalty (1e3) for `V = −50 V`; tune in `/gsd-discuss-phase` based on initial fit quality.
   - **RESOLVED:** Adopt the recommended scheme verbatim. The Plan 26-03 `calibrate_graded_doping_2d` cost function uses equal weights `(1/N)` over the three known-W targets `{0 V, −10 V, −30 V}` from the Petringa dataset (cost contribution `sum(((W_sim - W_exp)/W_exp)**2)`) PLUS a hard `divergence_penalty=1e3` added whenever the cathode ramp fails to reach `V_target_for_convergence_only=-50.0 V`. This preserves v3.0 low-bias behavior (PITFALLS P24 protected by the regression sweep in Plan 26-04) while making non-convergence at −50 V dominate the cost surface so Nelder-Mead favors solutions that converge. The two knobs (`V_target_for_convergence_only`, `divergence_penalty`) are exposed as keyword arguments for retuning if Plan 26-03 Task 3's full-range R² test fails on the first calibration pass.

3. **Where does the extended `reset_devsim_fully()` live?**
   - What we know: Currently `reset_devsim` is inline at `optimization.py:153`; Phase 31 will extend it again for tensor-mobility globals.
   - What's unclear: New module (`src/devsim_reset.py`) vs add to `src/device2d.py` vs add to `src/__init__.py`.
   - Recommendation: `src/devsim_reset.py` — single-responsibility, easy to import from any module including alt-structures' own teardown, easy to extend in Phase 31. Defer to planner; this is a `Claude's discretion` style choice.
   - **RESOLVED:** New module `src/devsim_reset.py`. Plan 26-02 Task 1 creates the file with two exports: the public `reset_devsim_fully(preserve_solver=True)` function and the module-level constant tuple `_CYLINDRICAL_GLOBALS` (7 strings: `raxis_zero`, `raxis_variable`, `node_volume_model`, `edge_couple_model`, `element_edge_couple_model`, `element_node0_volume_model`, `element_node1_volume_model`). The same task refactors `src/optimization.py` to `from src.devsim_reset import reset_devsim_fully` and removes the inline `devsim.reset_devsim()` block at lines 96-158. Single-responsibility module makes the Phase 31 tensor-mobility extension a one-tuple append rather than a code-archaeology exercise.

4. **Notebook scope for the regression: 20 or 22?**
   - What we know: notebooks/ has 22 files. PROJECT.md and prior milestones reference "20 notebooks."
   - Recommendation: Planner counts notebooks numbered 01–20 with scientific content; treats `03_executed.ipynb` as a re-run not a deliverable. Wave 0 produces the explicit list in `tests/baselines/v3_frozen.json`.
   - **RESOLVED:** 20 notebooks, scientific deliverables only. Plan 26-01 Task 2 hardcodes the canonical 20-notebook list (`01_phase1_validation.ipynb` through `20_feasibility_report.ipynb`) into `scripts/freeze_v3_baselines.py` and writes it as the `notebook_list` field of `tests/baselines/v3_frozen.json`. Duplicates / re-runs (`03_executed.ipynb`, `05_parametric_studies.ipynb`) are explicitly EXCLUDED — these are not deliverables per the v3.0 milestone audit. Plan 26-04 Task 2 reads `notebook_list` and asserts `len == 20`; the regression sweep re-executes exactly those 20 files via nbclient.

5. **Should this phase also re-snap notebook 14 (validation) baselines?**
   - What we know: `notebook 14` (`14_validation.ipynb`) holds v3.0 validation snapshots; PITFALLS P24 says "baseline freeze vs rebaseline" is a milestone-level decision.
   - Recommendation: Freeze for this phase; defer any rebaseline to Phase 34 (milestone audit). Plan accordingly — no edits to notebook 14 in Phase 26.
   - **RESOLVED:** Freeze, do NOT re-snap. Phase 26 makes no edits to `notebooks/14_validation.ipynb` (or any other notebook in the canonical 20). Plan 26-04 Task 2's `scripts/regression_sweep_v3_notebooks.py` re-EXECUTES each notebook via nbclient but explicitly does NOT call `nbformat.write` — the on-disk notebooks are read-only during the sweep. The acceptance criterion `grep -c "nbformat.write" scripts/regression_sweep_v3_notebooks.py` returns 0 enforces this contract. Any decision to rebaseline notebook 14 is deferred to the Phase 34 milestone audit (per the v4.0 milestone roadmap).

## Project Constraints (from CLAUDE.md)

> `./CLAUDE.md` does **not exist** in the working directory. [VERIFIED: `ls /Users/ngcex/projects/physics/petringa/CLAUDE.md` → "No such file or directory"]

Project-wide constraints discovered from auto-memory and STATE.md instead:

- **Use uv, not pip/venv, for all Python dependency management.** [Source: `memory/feedback_uv_venv.md`] — Calibration scripts run with `uv run pytest ...`; if `physdata` ever needs adding (Phase 27, not this phase), use `uv add` not `uv pip install`.
- **`device.py` (1D) is frozen** to protect 20 validated notebooks. [Source: STATE.md "Decisions"] — Phase 26 must NOT touch `src/device.py`. Calibration ports the _pattern_ (Nelder-Mead loop, bounds, cost) but the new function lives in `device2d.py` or a sibling.
- **`charge_error=1e10` required for all BDF1 transient solves.** [Source: STATE.md] — Calibration uses `type="dc"` only, so this constraint does not bite here; flag for Phase 30/29 transient work.
- **x = lateral, y = depth coordinate convention for all 2D modules.** [Source: STATE.md] — Existing in `device2d.py`; preserve.
- **Doping-profile uniform N_D fails at reverse bias in 2D.** [Source: `memory/project_doping_profile.md`] — This is the exact gap Phase 26 closes; confirms problem statement.

## Sources

### Primary (HIGH confidence — internal code + project planning docs)

- `src/device.py:468–637` — 1D `calibrate_graded_doping` reference implementation to port
- `src/device2d.py:55–467` — existing 2D mesh + graded doping; site of the calibration result
- `src/charge_collection_2d.py:50–166` — `_robust_dc_solve()` + `create_2d_dd_device` wrapper to reuse
- `src/cv_analysis.py:121–227` — `cv_sweep` C-V harness, already has the high-bias fallback pattern
- `src/optimization.py:96–158` — current `reset_devsim()` site with solver save/restore pattern
- `src/alternative_structures.py:573–584` — exact location where cylindrical globals are set (the leak source)
- `tests/test_device2d.py` — existing 4 test classes; pattern to extend
- `.planning/research/PITFALLS.md` — full v4.0 pitfalls catalog (P03, P20, P24, P27, P30 directly relevant)
- `.planning/research/SUMMARY.md` — v4.0 synthesis, decision D4 (regression baseline policy)
- `.planning/research/STACK.md` — confirms no new packages required for this phase
- `.planning/research/FEATURES.md` — confirms Phase 26 is "CALIBRATE existing" (not new module)
- `.planning/v3.0-MILESTONE-AUDIT.md` — confirms 22 notebooks in v3.0 directory; no audit gap mentions −15 V ceiling explicitly (this is "known tech debt from STATE.md")
- `.planning/STATE.md` lines 90–95 — explicit recognition of "N_D uniform in 2D fails at reverse bias → Phase 26"

### Secondary (MEDIUM confidence — external documentation, used to corroborate stack assumptions)

- [devsim Models Documentation](https://devsim.net/models.html) — confirms `edge_model` + `unitx/unity` framework; relevant for understanding why graded-doping `node_model` evaluation is correct in current code [CITED]
- [devsim forum thread on global state](https://forum.devsim.org/) — corroborates PITFALLS P03 finding that `raxis_*` parameters are global [CITED via SUMMARY.md]

### Tertiary (LOW confidence — none required for this phase)

None — the entire phase scope is covered by primary sources.

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH — no new dependencies, all helpers exist and are battle-tested through v3.0
- Architecture: HIGH — files-to-edit verified by direct inspection; calibration pattern is a direct port
- Pitfalls: HIGH — drawn from project-specific `.planning/research/PITFALLS.md` (severity ratings already triaged)
- Validation: MEDIUM-HIGH — pytest infrastructure exists; the 4 new tests are well-specified but the v3 baseline JSON does not yet exist (Wave 0 work)
- Root cause hypothesis: MEDIUM — three plausible failure modes (H1/H2/H3) and the differential between 1D and 2D is not yet diagnosed; this is intentionally an open question for Wave 0 to resolve before committing to a calibration

**Research date:** 2026-05-17
**Valid until:** 2026-06-17 (30 days — stable scientific domain; no fast-moving APIs in scope; devsim 2.10.0 stable since Oct 2025)
