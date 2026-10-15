---
id: C-08
source: raw-C F8
title: Radiation Damage page accepts fluence_min <= 0 or min > max: zero crashes with an uncaught ValueError, negative silently yields NaN fluences
severity: low
category: bug
status: open
location:
  - app/workflows/radiation_damage.py:62-74 (fluence inputs, no min_value)
  - app/workflows/radiation_damage.py:96-110 (geomspace outside any ValueError handler)
labels: [app, input-validation]
---

## Summary

The `fluence_min` / `fluence_max` number_inputs have no `min_value`, and there is no validation
before `np.geomspace(fluence_min, fluence_max, n)`:

- min = 0 raises `ValueError: Geometric sequence cannot include zero`. It is not caught (`except RuntimeError` only), so the user sees a raw traceback.
- min < 0 does not raise, contrary to the raw finding's "negative values likewise". geomspace returns NaN interior points (`[-1e13, nan, nan, 1e16]`), which are passed silently into `run_radiation_damage`. This is worse than a crash.
- min > max gives a descending fluence axis with no message.

## Reproduction / Evidence

```
uv run python -c "import numpy as np; ..."
(0.0, 1e16)     ValueError Geometric sequence cannot include zero
(-1e13, 1e16)   [-1.e+13     nan     nan  1.e+16]
(1e16, 1e13)    [1.e+16 1.e+15 1.e+14 1.e+13]
```

Code read: radiation_damage.py:62-74 has no min_value. Lines 96-105 show `try: fluences =
np.geomspace(...); ... except RuntimeError`.

## Impact

Raw traceback for zero. NaN fluences flow into the NIEL damage model with undefined downstream
behavior for negative values.

## Proposed fix

```diff
         fluence_min = st.number_input(
             "Min fluence (p/cm²)",
             format="%.3e",
+            min_value=1e8,
             value=1e13,
 ...
         fluence_max = st.number_input(
             "Max fluence (p/cm²)",
             format="%.3e",
+            min_value=1e8,
             value=1e16,
 ...
     if st.button("Run simulation"):
+        if fluence_min >= fluence_max:
+            st.error("Min fluence must be smaller than max fluence.")
+            st.stop()
         try:
             ...
-        except RuntimeError as e:
+        except (RuntimeError, ValueError) as e:
```

## Test plan

AppTest cases: fluence_min=0 (rejected by the widget, or st.error shown, with no exception);
fluence_min=1e16, fluence_max=1e13 (st.error shown and the facade is not called).

## Effort

S
