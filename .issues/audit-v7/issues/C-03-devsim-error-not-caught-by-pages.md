---
id: C-03
source: raw-C F3
title: devsim.error is not a RuntimeError; pages only catch RuntimeError, so sidebar-permitted inputs crash with a raw traceback
severity: medium
category: bug
status: fixed
location:
  - app/components/device_sidebar.py:80-85 (substrate min_value=0.0), 97-103, 109-114, 132-134 (N_A, N_D, N_D_junction, N_D_bulk, area: no min_value)
  - app/workflows/cv.py:43, cce.py:55, field_map.py:47, dark_current.py:144, batch_sweep.py:131, radiation_damage.py:105 (`except RuntimeError` only)
labels: [app, error-handling, devsim]
---

## Summary

`devsim.error.__mro__ == (devsim_py3.error, Exception, BaseException, object)`. Only
`solve_equilibrium` and `ramp_bias` translate it into RuntimeError (poisson.py:165-166).
Mesh and doping construction failures propagate as raw `devsim.error`. The sidebar allows
`substrate_thickness_um = 0.0` (`min_value=0.0`). Clicking Run on C-V / CCE / Field /
Dark Current with that value raises an uncaught `devsim_py3.error: DEVSIM FATAL ...`, and
Streamlit shows a red traceback.

Correction to the raw finding: the raw file says the process-global devsim state is "undefined"
after the FATAL. In my run, P4 (`run_cv(area_cm2=0)`) succeeded in the same process right after
P3's FATAL, so subsequent runs are not obviously broken. That sub-claim is softened here.

## Reproduction / Evidence

```
uv run python -c "import devsim; print(devsim.error.__mro__)"
-> (<class 'devsim_py3.error'>, <class 'Exception'>, <class 'BaseException'>, <class 'object'>)

uv run python <scratchpad>/repro_api.py   (P3)
-> (<class 'devsim_py3.error'>, <class 'Exception'>, ...) DEVSIM FATAL: There was a fatal
   exception in the program.  Shutting down.  Please inspect program output.
```

`grep -n "except" app/workflows/*.py`: every devsim-calling page uses `except RuntimeError` only.
microdosimetry catches `(ValueError, KeyError)` and is not devsim-backed.
`app/components/device_sidebar.py:80-85` has `min_value=0.0` for the substrate.

## Impact

Inputs that the UI explicitly permits produce a stack trace instead of an error message. The same
path is hit by the NaN/inf batch-sweep values (C-10) and by concurrent sessions (C-04).

## Proposed fix

```diff
# etna/api/__init__.py / errors.py
+class SimulationError(RuntimeError):
+    """Raised by etna facades for any solver/mesh failure (wraps devsim.error)."""

# etna/api/simulation.py: in each facade (or a decorator)
+    except devsim.error as e:
+        raise SimulationError(str(e)) from e

# app/components/device_sidebar.py
-        min_value=0.0,
+        min_value=0.01,
# and min_value=1e10 (or any >0) for N_A / N_D / N_D_junction / N_D_bulk, min_value=1e-8 for area_cm2
```

Pages keep `except RuntimeError` and catch `SimulationError` for free, so the app never needs to
import devsim. Also see C-12 for `DeviceConfig.__post_init__` validation.

## Test plan

- AppTest for each page with a fake facade raising `devsim.error("x")`, or `SimulationError` after the fix. Assert that `st.error` is shown and `at.exception` is empty.
- Unit test: `run_cv(DeviceConfig(substrate_thickness_um=0.0))` raises `SimulationError`, or `ValueError` once C-12 lands.

## Effort

S

## Status 2026-10-18 (final-c): fixed (cf64bef)

- Every devsim-calling page now catches `Exception` and routes it through `app.components.page_utils.show_error`, which treats the devsim `error` class (by name) as a solver failure (friendly st.error + details expander). DeviceConfig.__post_init__ rejects substrate_thickness_um <= 0 and the other non-positive inputs, and the sidebar shows the validation error instead of building.
- Evidence: `grep -c "except RuntimeError" app/workflows/*.py` -> 0 in every file; `uv run pytest -q tests/test_app_robustness.py -k "devsim_error or stale or invalid_sidebar or rejects_out_of_range"` -> 15 passed (test_devsim_error_shows_friendly_error_with_details).
