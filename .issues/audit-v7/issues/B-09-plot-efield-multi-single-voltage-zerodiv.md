---
id: B-09
title: plot_electric_field_multi raises ZeroDivisionError for a single-voltage sweep
severity: low
category: bug
status: fixed 2026-10-18 (wp-c2 8c8dec4: plot_electric_field_multi removed, no caller)
location: etna/core/plotting.py:100-104
labels: [plotting, edge-case, raw-F9]
---

## Summary

`n_curves = min(8, len(voltages))`, followed by `colors = [cmap(i / (n_curves - 1)) ...]`. With one voltage this divides by zero. The parametric plots in the same file already use the guarded pattern `i / max(n - 1, 1)` (lines 765, 828, 893).

## Reproduction / Evidence

- call chain (read): plotting.py:101 `n_curves = min(8, len(voltages))`, then :104 `colors = [cmap(i / (n_curves - 1)) for i in range(n_curves)]`. With `n_curves == 1`, `range(1)` yields i=0 and `0 / 0` raises ZeroDivisionError.
- ran: `grep -rn "cmap(.*max\|/ max(" etna/core/plotting*.py` finds the guarded form at lines 765, 828 and 893.
- ran: `grep -rln plot_electric_field_multi etna app tests scripts notebooks` finds only `etna/core/plotting.py` and `tests/test_plotting.py`, so there are no production callers (see B-12).
- From the raw audit, not re-run: `plot_electric_field_multi({"voltages":[0.0], "E_fields":[(x, E)]})` raises `ZeroDivisionError`.

## Impact

Low. A one-voltage call crashes. The function currently has no callers outside tests.

## Proposed fix

```diff
-    colors = [cmap(i / (n_curves - 1)) for i in range(n_curves)]
+    colors = [cmap(i / max(n_curves - 1, 1)) for i in range(n_curves)]
```

Alternatively, delete the function as part of B-12.

## Test plan

`plot_electric_field_multi({"voltages": [0.0], "E_fields": [(x, E)]})` returns an Axes with one line.

## Effort

S
