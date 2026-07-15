---
phase: 36-core-api-deviceconfig-c-v-field-vertical-slice
reviewed: 2026-07-02T00:00:00Z
depth: standard
files_reviewed: 8
files_reviewed_list:
  - examples/cv_example.py
  - etna/__init__.py
  - etna/api/__init__.py
  - etna/api/device.py
  - etna/api/results.py
  - etna/api/simulation.py
  - tests/test_api_cv.py
  - tests/test_api_field.py
findings:
  critical: 1
  warning: 3
  info: 2
  total: 6
status: issues_found
---

# Phase 36: Code Review Report

**Reviewed:** 2026-07-02T00:00:00Z
**Depth:** standard
**Files Reviewed:** 8
**Status:** issues_found

## Summary

Reviewed the `DeviceConfig` -> `build_device` -> `run_cv`/`run_field` -> `SimResult` vertical slice. The 1D path (`run_cv`, `run_field` with `half_width_um=None`) is correctly wired: unit conversions (um -> cm, cm -> um) are consistent, the reverse-bias sign convention matches `core.cv_analysis.cv_sweep`'s documented contract, device cleanup is structured so a `build_device()` failure can't cause a `NameError` in `finally`, and `N_D=None` is handled safely by both 1D and 2D core constructors.

The 2D path in `run_field()` has one contract-breaking defect: it returns the device's _lateral_ coordinate as `SimResult.x` labeled and documented as "depth," which the design spec explicitly requires to be `x=depth_um` (`docs/superpowers/specs/2026-06-26-simulator-library-ui-design.md:152`). No test exercises the 2D `run_field` result's `x`/`y` fields against this contract — `tests/test_api_field.py::TestRunFieldIntegration2D` only checks `mesh`, never `result.x`/`result.y`.

Two further robustness issues affect any process embedding `etna` as a library: `run_cv`/`run_field` call `reset_devsim_fully()` unconditionally at entry, which deletes _every_ devsim device in the process (not just etna's own), and the cleanup `finally` blocks swallow `delete_device` failures with no logging, silently falling back to a full reset.

## Critical Issues

### CR-01: `run_field()` returns the wrong array as `SimResult.x`/`.y` for 2D devices, breaking the design-spec depth-profile contract

**File:** `etna/api/simulation.py:236-243`
**Issue:**
The design spec (`docs/superpowers/specs/2026-06-26-simulator-library-ui-design.md:151-152`) specifies the `run_field()` contract as:

```
result = run_field(cfg, bias_V=-100)
# → SimResult with x=depth_um, y=field_V_cm, metadata includes potential, doping
```

and further states that `04_field_map.py` "renders electric field and potential depth profiles" from this result, with 2D rendering handled separately via `MeshData` on a regular grid (spec line 229).

For 1D devices, `x_coords` genuinely is depth (the 1D mesh in `etna/core/device.py` runs along depth only), so `x = x_coords * 1e4` correctly yields `depth_um`. For 2D devices, however, `create_sic_2d_device` builds a mesh where `x` is the _lateral_ half-width coordinate and `y` is depth (`etna/core/device2d.py:169-170`: "Structure uses symmetry: x = 0 (center/symmetry plane) ... Depth: y = 0 (anode, top) to total_depth (cathode, bottom)"). `run_field()` unconditionally sets:

```python
return SimResult(
    ...
    x=x_coords * 1e4,      # lateral position for 2D, NOT depth
    y=field_nodes,         # per-node E magnitude, unsorted, multivalued vs x
    ...
)
```

For a 2D device this means `result.x` is lateral position mislabeled as depth, and `result.y` is a flat, unordered per-node array — for a given lateral `x` there are many nodes at different depths all sharing that `x`, so `(x, y)` is not a valid function/profile at all (any naive `plot(result.x, result.y)` produces multivalued, effectively random-looking scatter, not a depth profile). Any downstream consumer following the documented contract (e.g., a future `04_field_map.py` depth-profile renderer that isn't 2D-aware) will silently plot garbage instead of raising an error. `run_field`'s own docstring (`simulation.py:133`, `"x=depth (um)"`) does not caveat that this is false for the 2D branch — it directly asserts an incorrect contract.

**Fix:** Either (a) restrict `SimResult.x`/`.y` to a genuine depth cut for 2D devices (e.g. along the lateral symmetry line `x=0`, mirroring how `MeshData` is the authoritative source for full 2D rendering), or (b) explicitly document/raise that `SimResult.x`/`.y` are not meaningful for 2D `run_field()` results and downstream code must use `mesh` instead — mirroring the `NotImplementedError` pattern `run_cv()` already uses for its unsupported 2D case:

```python
if is_2d:
    # SimResult.x/.y are not a valid depth profile for 2D devices (x is the
    # lateral mesh coordinate, not depth); a real depth cut isn't defined
    # without picking a lateral slice. Consumers must use `mesh` instead.
    x_out = np.array([])
    y_out = np.array([])
else:
    x_out = x_coords * 1e4
    y_out = field_nodes
...
return SimResult(..., x=x_out, y=y_out, ..., mesh=mesh)
```

and update the docstring to state this explicitly. At minimum, add a regression test asserting `result.x`/`result.y` semantics (or their documented emptiness) for the 2D case in `tests/test_api_field.py::TestRunFieldIntegration2D`.

## Warnings

### WR-01: `reset_devsim_fully()` at the top of `run_cv`/`run_field` deletes all devsim devices in the process, not just etna's own

**File:** `etna/api/simulation.py:71`, `etna/api/simulation.py:138`
**Issue:** Both public facades unconditionally call `reset_devsim_fully()` before building their own device. Per `etna/core/devsim_reset.py:100-109`, this enumerates and deletes _every_ device currently registered with devsim in the process — including devices created by other, unrelated code sharing the same process (devsim state is process-global, not etna-scoped). For a one-shot script this is harmless, but `run_cv`/`run_field` are public library API entry points (`etna/__init__.py` exports them at top level); a caller that holds another live devsim device across two `etna.run_cv()` calls (e.g. in a notebook, a web app request handler, or a test suite that builds a device via `etna.core.*` directly and then calls `run_cv()`) will have that unrelated device silently deleted with no warning.
**Fix:** At minimum, document this destructive side effect prominently in both docstrings ("Calling this function deletes all devsim devices currently in the process, not just those created by etna"). Consider scoping cleanup to only the device this function itself creates (already handled correctly in the `finally` block) and dropping the blanket `reset_devsim_fully()` call at entry, reserving the full reset for cases where cylindrical-coordinate leakage is actually suspected (per `devsim_reset.py`'s own stated purpose).

### WR-02: Cleanup `finally` blocks swallow `delete_device` failures silently

**File:** `etna/api/simulation.py:97-102`, `etna/api/simulation.py:244-249`
**Issue:**

```python
finally:
    try:
        devsim.delete_device(device=device_info["device_name"])
    except Exception:
        reset_devsim_fully()
```

If `devsim.delete_device` raises, the exception is discarded and silently replaced with a full `reset_devsim_fully()` call — with no logging of what went wrong. This makes debugging device-cleanup failures (e.g., a device name collision, or an unexpected devsim internal state) effectively impossible from the caller's side, and additionally amplifies WR-01's blast radius: a failed single-device delete escalates straight to deleting every device in the process without any log line.
**Fix:**

```python
finally:
    try:
        devsim.delete_device(device=device_info["device_name"])
    except Exception:
        logger.warning(
            "run_cv: delete_device(%r) failed, falling back to full reset",
            device_info["device_name"], exc_info=True,
        )
        reset_devsim_fully()
```

(requires adding `logger = logging.getLogger(__name__)` to `etna/api/simulation.py`, consistent with the logging pattern already used in `etna/core/cv_analysis.py` and `etna/core/devsim_reset.py`).

### WR-03: No test coverage for `run_cv()`'s documented 2D `NotImplementedError` or for `build_device()` directly

**File:** `tests/test_api_cv.py`, `etna/api/simulation.py:59-65`
**Issue:** `run_cv()`'s docstring and implementation explicitly guard against 2D configs:

```python
if config.half_width_um is not None:
    raise NotImplementedError(...)
```

but no test in `tests/test_api_cv.py` (or elsewhere) exercises this branch. Likewise, `build_device()` (`etna/api/device.py`) — the dispatch point between the 1D and 2D construction paths — has no direct unit test; it is only exercised indirectly through `run_cv`/`run_field`. A regression that silently removed the 2D guard in `run_cv`, or that broke the `device_name` uniqueness/uuid logic in `build_device`, would not be caught by the current test suite.
**Fix:** Add a fast (non-`slow`) test:

```python
def test_run_cv_rejects_2d_config():
    with pytest.raises(NotImplementedError):
        run_cv(DeviceConfig(half_width_um=50.0))
```

and a direct `build_device()` unit test asserting `device_info["dd_initialized"] is True` and the expected `dimension` key for both 1D and 2D configs.

## Info

### IN-01: `run_cv`/`run_field` do not guard against `n_points <= 0` or non-finite bias arguments

**File:** `etna/api/simulation.py:23-27`, `etna/api/simulation.py:67`
**Issue:** `run_cv(config, v_start=0, v_stop=-200, n_points=0)` silently produces `np.linspace(0, -200, 0) == []`, so `cv_sweep` runs its loop zero times and `run_cv` returns a `SimResult` with empty `x`/`y`/`metadata["depletion_widths"]` arrays and no error or warning. This is a valid-but-degenerate input that a caller could pass accidentally (e.g. from a UI slider defaulting to 0); the resulting empty `SimResult` looks superficially like a success.
**Fix:** Add a lightweight guard, e.g. `if n_points < 1: raise ValueError("n_points must be >= 1")`, at the top of `run_cv()`.

### IN-02: `run_field`'s `contacts` list reports the same scalar `position` for 2D devices, ignoring lateral extent

**File:** `etna/api/simulation.py:213-216`
**Issue:** `contacts = [{"name": "anode", "position": 0.0}, {"name": "cathode", "position": float(device_info["total_length"])}]` reports only a depth-axis scalar position for both contacts, regardless of dimensionality. For 2D devices this is likely intentional (contacts span the full lateral extent at fixed depth), but nothing in `MeshData`'s `contacts: list[dict]` field documents this convention, so a geometry-viewer consumer has no way to know whether `position` means "depth only, full lateral span" versus a specific (x, y) point.
**Fix:** Document the convention inline (e.g. a comment: `# position is depth only; contacts span the full lateral extent at this depth for 2D devices`), or extend the contact dict with an explicit `x_min`/`x_max` for the 2D case to mirror the `regions` dict's own explicit min/max fields.

---

_Reviewed: 2026-07-02T00:00:00Z_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
