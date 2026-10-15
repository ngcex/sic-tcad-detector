---
id: C-10
source: raw-C F10
title: Batch Sweep value parser accepts nan, inf and non-positive values, which are injected straight into DeviceConfig
severity: low
category: bug
status: open
location:
  - app/workflows/batch_sweep.py:93-100
labels: [app, input-validation]
---

## Summary

`values = [float(v.strip()) for v in values_raw.split(",") if v.strip()]` only guards against
`ValueError`. `float("nan")`, `float("inf")` and negative numbers parse successfully and are
passed to `ParametricSweep`. ParametricSweep validates only the field name (via `dataclasses.replace`).
`grep -rn isfinite etna/api app` finds no downstream finiteness check.

## Reproduction / Evidence

- Own run of repro_api.py, P8: `[float(v.strip()) for v in "10, nan, inf".split(",")]` -> `[10.0, nan, inf]`.
- Code read of batch_sweep.py:93-100. The only checks are ValueError and an empty list.
- Not probed: what devsim does with `epi_thickness_um=nan`. The raw file's "devsim.error (uncaught, C-03) or a hang" is inferred, not observed.

## Impact

Unvalidated non-finite or non-physical geometry, doping or T reaches the mesh builder. The
likely outcomes are an uncaught devsim error (C-03) or a pathological solve. Low severity, because
it requires deliberate input.

## Proposed fix

```diff
         try:
             values = [float(v.strip()) for v in values_raw.split(",") if v.strip()]
         except ValueError:
             st.error("Values must be a comma-separated list of numbers.")
             st.stop()
+        if not all(math.isfinite(v) and v > 0 for v in values):
+            st.error("Values must be finite, positive numbers.")
+            st.stop()
```

All current SWEEPABLE_FIELDS are strictly positive physical quantities. L_transition_um=0 may be
legitimate (abrupt junction); if so, allow `>= 0` for that field only. Long-term, C-12's
`DeviceConfig.__post_init__` provides the same guarantee for library users.

## Test plan

AppTest: values "10, nan" -> st.error shown and no facade call. "10, -5" -> st.error shown.
"10, 15" -> runs as before.

## Effort

S
