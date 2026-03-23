# Phase 7: Solver Robustness - Research

**Researched:** 2026-03-21
**Domain:** devsim equation re-registration / documentation alignment
**Confidence:** HIGH

## Summary

Phase 7 addresses two tech debt items identified in the v1.0 milestone audit. Both are well-scoped and low-risk.

**Bug 1 -- `time_node_model` dropped on equation re-registration:** Both `add_generation_to_dd` (in `src/charge_collection.py`, lines 291-308) and `add_auger_recombination` (in `src/flash_recombination.py`, lines 104-121) call `devsim.equation()` to re-register the ElectronContinuityEquation and HoleContinuityEquation after updating node models. Neither call passes `time_node_model`. The original registration in `src/drift_diffusion.py` (lines 143-168) sets `time_node_model="NCharge"` for electrons and `time_node_model="PCharge"` for holes. Since `devsim.equation()` performs a full replacement when re-registering (confirmed by devsim's "Replacing equation with equation of the same name" warning), omitting `time_node_model` clears the transient term. This is a latent bug -- it has no effect on the current DC-only solver workflow, but would break any future transient (`transient_bdf1`, `transient_bdf2`, `transient_tr`) solves.

**Bug 2 -- ROADMAP SC-3 wording:** The Phase 4 success criterion SC-3 (ROADMAP.md line 112) currently reads: "CCE vs dose-rate curve ... shows physically meaningful trend (flat CCE confirms Auger is negligible -- valid null result)". This was partially updated during Phase 4 execution to acknowledge the null result, but the wording still contains the ambiguous phrase "shows physically meaningful trend". Phase 7 should make the wording unambiguously state that flat CCE is the accepted scientific finding, not a failure to observe degradation.

**Primary recommendation:** Add `time_node_model` parameters to all four `devsim.equation()` re-registration calls, add a regression test verifying transient capability is preserved after generation/Auger setup, and update ROADMAP SC-3 wording.

## Standard Stack

### Core

| Library | Version     | Purpose                | Why Standard                                                            |
| ------- | ----------- | ---------------------- | ----------------------------------------------------------------------- |
| devsim  | (installed) | TCAD device simulation | Project's core solver; `devsim.equation()` API is the target of the fix |

No new libraries needed. This phase is purely internal code and documentation fixes.

## Architecture Patterns

### Pattern 1: devsim.equation() Re-Registration with Full Parameter Set

**What:** When re-registering an equation (to update `node_model` after changing generation/recombination expressions), ALL parameters from the original registration must be re-specified, including `time_node_model`.

**When to use:** Any time `devsim.equation()` is called on an already-registered equation name.

**Example:**

```python
# CORRECT: preserves time_node_model
devsim.equation(
    device=device,
    region=region,
    name="ElectronContinuityEquation",
    variable_name="Electrons",
    node_model="ElectronGeneration",
    edge_model="ElectronCurrent",
    time_node_model="NCharge",         # MUST include
    variable_update="positive",
)

# WRONG: strips time_node_model (current buggy code)
devsim.equation(
    device=device,
    region=region,
    name="ElectronContinuityEquation",
    variable_name="Electrons",
    node_model="ElectronGeneration",
    edge_model="ElectronCurrent",
    variable_update="positive",
    # time_node_model MISSING -- transient capability lost
)
```

**Source:** devsim API (`help(devsim.equation)`) confirms `time_node_model` is optional, meaning omission = cleared.

### Anti-Patterns to Avoid

- **Partial equation re-registration:** Never assume devsim preserves unspecified optional parameters on re-registration. It performs a full replacement. Always pass the complete set.

## Don't Hand-Roll

| Problem               | Don't Build                                       | Use Instead                                 | Why                                                                                                                                                                   |
| --------------------- | ------------------------------------------------- | ------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Transient model names | Storing model names in variables or a config dict | Hardcode `"NCharge"` and `"PCharge"` inline | These are stable constants defined once in `setup_sic_drift_diffusion`; indirection adds complexity for no benefit. The original registration is the source of truth. |

## Common Pitfalls

### Pitfall 1: Forgetting time_node_model on Re-Registration

**What goes wrong:** Transient solves silently produce wrong results or fail because the charge storage term (dQ/dt) is missing from the continuity equation.
**Why it happens:** devsim's `equation()` API treats all keyword arguments as optional. Omitting `time_node_model` doesn't raise an error -- it simply clears the time-dependent term.
**How to avoid:** Always pass `time_node_model` when re-registering continuity equations. A unit test should verify transient capability is preserved.
**Warning signs:** DC solves work perfectly fine; only transient solves fail or produce unphysical results.

### Pitfall 2: Testing Only DC After Equation Changes

**What goes wrong:** The bug is latent because the entire current test suite uses `type="dc"` solves. A DC solve doesn't exercise `time_node_model`.
**Why it happens:** The project doesn't currently use transient solves (continuation solver uses DC steps).
**How to avoid:** Add a minimal regression test that verifies `time_node_model` is still set after calling `add_generation_to_dd` and `add_auger_recombination`. This can be done by checking the equation attributes or attempting a trivial transient step.

## Code Examples

### Fix for add_generation_to_dd (charge_collection.py, lines 291-308)

Current code (buggy):

```python
devsim.equation(
    device=device, region=region,
    name="ElectronContinuityEquation",
    variable_name="Electrons",
    node_model="ElectronGeneration",
    edge_model="ElectronCurrent",
    variable_update="positive",
)
devsim.equation(
    device=device, region=region,
    name="HoleContinuityEquation",
    variable_name="Holes",
    node_model="HoleGeneration",
    edge_model="HoleCurrent",
    variable_update="positive",
)
```

Fixed code:

```python
devsim.equation(
    device=device, region=region,
    name="ElectronContinuityEquation",
    variable_name="Electrons",
    node_model="ElectronGeneration",
    edge_model="ElectronCurrent",
    time_node_model="NCharge",
    variable_update="positive",
)
devsim.equation(
    device=device, region=region,
    name="HoleContinuityEquation",
    variable_name="Holes",
    node_model="HoleGeneration",
    edge_model="HoleCurrent",
    time_node_model="PCharge",
    variable_update="positive",
)
```

### Fix for add_auger_recombination (flash_recombination.py, lines 104-121)

Identical pattern: add `time_node_model="NCharge"` and `time_node_model="PCharge"` to the two `devsim.equation()` calls.

### ROADMAP SC-3 Wording Fix

Current (line 112):

```
3. CCE vs dose-rate curve spanning 20 to 230 Gy/s at reference conditions (-30V, 10 um epi, 62 MeV protons) shows physically meaningful trend (flat CCE confirms Auger is negligible — valid null result)
```

Suggested replacement:

```
3. CCE vs dose-rate curve spanning 20 to 230 Gy/s at reference conditions (-30V, 10 um epi, 62 MeV protons) produces flat CCE (~1.0) confirming Auger recombination is negligible at therapeutic FLASH dose rates — an accepted null result consistent with delta_n << Auger threshold (~1e16 cm^-3)
```

## State of the Art

Not applicable -- this phase is internal tech debt, not a new feature or technology adoption.

## Open Questions

1. **Regression test approach for transient capability**
   - What we know: DC solves don't exercise `time_node_model`. A transient solve test would be the most thorough verification, but the project doesn't currently have transient solve infrastructure.
   - What's unclear: Whether a minimal transient solve test (single BDF1 step) is stable enough for CI, or whether checking equation metadata is sufficient.
   - Recommendation: Use a single `transient_bdf1` step with relaxed tolerances as the regression test. If devsim doesn't expose equation metadata for direct inspection, the transient solve test is the only way to verify. If the transient test is too fragile, fall back to a source-level assertion (grep the code for `time_node_model` in equation calls -- less robust but better than nothing).

## Affected Files

| File                                               | Change                                                | Lines   |
| -------------------------------------------------- | ----------------------------------------------------- | ------- |
| `src/charge_collection.py`                         | Add `time_node_model` to 2 `devsim.equation()` calls  | 291-308 |
| `src/flash_recombination.py`                       | Add `time_node_model` to 2 `devsim.equation()` calls  | 104-121 |
| `.planning/ROADMAP.md`                             | Update SC-3 wording                                   | 112     |
| `tests/test_charge_collection.py` or new test file | Regression test for transient capability preservation | New     |

## Sources

### Primary (HIGH confidence)

- `devsim.equation()` API documentation (`help(devsim.equation)`) -- confirms `time_node_model` is optional, omission = no transient term
- Source code inspection: `src/drift_diffusion.py` lines 143-168 (original equation registration with `time_node_model`)
- Source code inspection: `src/charge_collection.py` lines 291-308 (re-registration without `time_node_model`)
- Source code inspection: `src/flash_recombination.py` lines 104-121 (re-registration without `time_node_model`)
- `.planning/v1.0-MILESTONE-AUDIT.md` -- tech debt items documenting both issues
- `.planning/ROADMAP.md` line 112 -- current SC-3 wording
- STATE.md decision [04-02] -- "Null result (no CCE degradation) is valid scientific finding"

## Metadata

**Confidence breakdown:**

- Bug identification: HIGH -- directly confirmed by reading source code and devsim API docs
- Fix approach: HIGH -- straightforward parameter addition; exact code changes identified
- ROADMAP wording: HIGH -- clear intent documented in milestone audit and Phase 7 description
- Regression test approach: MEDIUM -- transient solve test feasibility not fully verified

**Research date:** 2026-03-21
**Valid until:** indefinite (internal code fix, not dependent on external ecosystem changes)
