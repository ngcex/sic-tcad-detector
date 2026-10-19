---
id: C-07
source: raw-C F7
title: sim_type="dark_current" has two incompatible shapes; the facade's native bias-sweep result is mislabeled as T_K by to_csv_bytes / build_dark_current_figure, and truncation is not reported
severity: low
category: bug
status: fixed
location:
  - etna/api/simulation.py:742-754 (run_dark_current: x = bias V)
  - app/components/results.py:89-132 (figure x-axis "Temperature (K)"), :286-296 (CSV first column "T_K")
  - app/workflows/dark_current.py:130-140 (page re-packs a second "dark_current" SimResult with x = T)
  - etna/core/dark_current.py:499-505 (early `return results` on convergence failure)
labels: [api, app, labels, truncation]
---

## Summary

`run_dark_current` returns `sim_type="dark_current"` with x = bias voltages. The Dark Current
page aggregates several of those into another `sim_type="dark_current"` SimResult with
x = temperatures. The shared serializer and figure builder assume the page shape. A library user
running `to_csv_bytes(run_dark_current(cfg))` gets a first column named `T_K` that contains
bias voltages.

In addition, `dark_current_sweep` returns a partial array on a convergence failure
(dark_current.py:499-505, logger.warning and early return). run_dark_current then exposes no
`truncated`/`requested_v_stop` metadata, unlike run_cv/run_cce. The raw file cites lines 481-486
for this; the current location is 499-505.

Correction to the raw finding: the raw file says the page aggregate drops area provenance
("the dark-current CSV has no area provenance"). That is false. The aggregate carries
`config=cfg`, and `to_csv_bytes` writes `asdict(result.config)`, which includes `area_cm2`, into
the `# device:` header. Only `metadata["area_cm2"]` is dropped, and nothing reads it. That
sub-claim is removed.

## Reproduction / Evidence

Own run of repro_api.py, P5:

```
x: [-5.] sim_type: dark_current
T_K,I_total_A,I_SRH_A,I_TAT_A,I_SRV_A
-5.0,6.99...
```

Code read: results.py:286-296 hard-codes `"T_K": result.x`. results.py:127-128 titles the plot
"vs Temperature" / "Temperature (K)". dark_current.py (app):130-140 builds the aggregate with
`x=np.array(T_ok)`.

## Impact

Library users get a CSV and plot with a wrong axis label. Silent truncation of a bias sweep is
indistinguishable from a complete one. The app path itself is labeled correctly today.

## Proposed fix

```diff
# app/workflows/dark_current.py
-                sim_type="dark_current",
+                sim_type="dark_current_vs_T",

# app/components/results.py
-    elif result.sim_type == "dark_current":
-        df = pd.DataFrame({"T_K": result.x, ...})
+    elif result.sim_type == "dark_current_vs_T":
+        df = pd.DataFrame({"T_K": result.x, ...})
+    elif result.sim_type == "dark_current":
+        df = pd.DataFrame({"bias_V": result.x, "I_total_A": result.y, ...})
# build_dark_current_figure: dispatch the x-axis title on sim_type the same way

# etna/api/simulation.py run_dark_current metadata
+                "truncated": len(result["voltages"]) < n_points,
+                "requested_v_stop": v_stop,
```

## Test plan

- `to_csv_bytes(run_dark_current(DeviceConfig(), v_start=-5, v_stop=-5, n_points=1))`: the first header column is `bias_V`.
- Existing dark-current page AppTests still pass with the renamed aggregate sim_type, and the page CSV still starts with `T_K`.
- Unit test with a monkeypatched `dark_current_sweep` returning a short array: `metadata["truncated"] is True`.

## Effort

S

## Status 2026-10-18 (final-c): fixed (5ca8e98 + 1d57122)

- 5ca8e98: run_dark_current reports `truncated`, `requested_v_start/stop`, `reached_v_stop`.
- 1d57122: the Dark Current page aggregate is `sim_type="dark_current_vs_T"` (CSV `T_K,I_total_A`); run_dark_current keeps `"dark_current"` (x = bias) and to_csv_bytes exports it as `bias_V,I_total_A` (+ `past_reach_through` column for the legacy model, `# truncated:` and `# validity:` header lines); build_dark_current_figure titles the axis Temperature or Bias by sim_type.
- Evidence: `uv run pytest -q tests/test_app_csv_export.py tests/test_app_dark_current_page.py` -> passed (test_run_dark_current_shape_is_labelled_bias_not_temperature, test_run_dark_current_truncation_and_legacy_validity_are_exported; page CSV still starts with T_K).
