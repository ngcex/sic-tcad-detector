---
id: B-04
title: compare_cce_hecht_vs_dd raises ValueError (length mismatch) whenever the DD sweep truncates
severity: medium
category: bug
status: open
location: etna/core/charge_collection.py:1471-1505 (compare_cce_hecht_vs_dd); etna/core/charge_collection.py:612-623 (cce_vs_bias truncated return); etna/core/validation.py:62; etna/core/plotting.py:541-549 (plot_cce_comparison)
labels: [cce, hecht, truncation, robustness, raw-F4]
---

## Summary

`cce_vs_bias` returns only the bias points that converged (`V_range[reached_mask]`, `truncated=True`). `compare_cce_hecht_vs_dd` evaluates Hecht and partial-depletion Hecht on the **full** `V_range`. It then computes `np.abs(cce_dd - cce_hecht)` and calls `compute_agreement_metrics(cce_dd, cce_hecht)`, both of which need equal lengths. If any requested bias lies past the punch-through wall (about 60-66 V for the default device, lower for thinner or more lightly doped devices), the call raises instead of returning a partial result. `plot_cce_comparison` plots `cce_dd` against the full `voltages`, so it would also misalign if it were ever handed a truncated result.

## Reproduction / Evidence

- call chain (read, not re-run): charge_collection.py:1474-1475 sets `cce_dd = dd_result["cce_values"]` (reached points only). Line 1478 computes `hecht_cce(V_range)` (full length). Line 1498 does `np.abs(cce_dd - cce_hecht)`, which is a broadcast error unless len(cce_dd) is 1 or matches. validation.py:62 raises `ValueError("Array length mismatch: sim=..., exp=...")`.
- From the raw audit, not re-run: `compare_cce_hecht_vs_dd(np.array([-30.0, -90.0]))` raises `ValueError: Array length mismatch: sim=1, exp=2`.
- The only test (`TestCCEDDvsHecht`, slow) uses `[-30, -40, -60]`, which sits close to the wall.
- Notebook 03 calls this function with a user-chosen `V_range_comp`.

## Impact

A user who asks for a high-bias comparison gets an exception instead of a partial comparison flagged `truncated`. The existing test is fragile, because -60 V is near the punch-through wall.

## Proposed fix

```diff
     dd_result = cce_vs_bias(V_range, epi_thickness_cm=epi_thickness_cm, **kwargs)
     cce_dd = dd_result["cce_values"]
+    V_used = np.asarray(dd_result["voltages"], dtype=float)
-    cce_hecht = hecht_cce(V_range, d=epi_thickness_cm)
+    cce_hecht = hecht_cce(V_used, d=epi_thickness_cm)
 ...
-    cce_hecht_partial = hecht_cce_partial_depletion(V_range, epi_thickness_cm, W_func)
+    cce_hecht_partial = hecht_cce_partial_depletion(V_used, epi_thickness_cm, W_func)
 ...
     return {
-        "voltages": V_range,
+        "voltages": V_used,
+        "truncated": dd_result["truncated"],
+        "requested_v_stop": dd_result["requested_v_stop"],
```

Guard the metrics when `len(V_used) < 2`, where R^2 is undefined, by returning NaN metrics.

## Test plan

Monkeypatch `cce_vs_bias` to return a one-point truncated dict for a two-point request. Assert that the comparison returns arrays of length 1, `truncated=True`, and no exception. Run `plot_cce_comparison` on that result.

## Effort

S
