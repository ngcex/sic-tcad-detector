---
id: C-06
source: raw-C F6
title: C-V and Mott-Schottky figures and the 1/C^2 CSV column carry per-area units (F/cm^2, cm^4/F^2) while run_cv returns absolute Farads
severity: low
category: bug
status: fixed
location:
  - app/components/results.py:31 (`yaxis_title="Capacitance (F/cm²)"`)
  - app/components/results.py:43 (`yaxis_title="1/C² (cm⁴/F²)"`)
  - app/components/results.py:240 (CSV column `one_over_C2_cm4_per_F2`)
  - etna/api/simulation.py:111-117 (area forwarded; `1.0 / C**2`, so units are F^-2)
labels: [app, units, labels]
---

## Summary

`run_cv` passes `area=config.area_cm2` to `cv_sweep` -> `junction_capacitance(W, area)` =
eps*area/W. Its docstring (cv_analysis.py:60) says "F if area given, F/cm^2 if area=1.0", so the
result is in Farads. The comment at simulation.py:107-110 states this explicitly. The figure
axes still say F/cm² and cm⁴/F², and the CSV 1/C² column name encodes cm⁴/F².

Clarifications to the raw finding:

- The CSV capacitance column is already correctly named `capacitance_F`. Only the 1/C² column and the two figure axis titles are wrong.
- The raw fix says "keep `# area_cm2:` in the CSV header (already in metadata but not written)". In fact `area_cm2` is already written, as part of the `# device:` line built from `asdict(result.config)` (results.py:307-316). No header change is needed.

## Reproduction / Evidence

- P2b (own run): default area 1e-4 cm² gives C = 4.935e-13 F, about 0.5 pF, which is consistent with Farads. A per-area value would be about 4.9e-9 F/cm².
- `sed -n 235,250p app/components/results.py` shows columns `bias_V, capacitance_F, one_over_C2_cm4_per_F2, depletion_width_cm`.
- `grep -rn "cm4" tests` finds tests/test_app_csv_export.py:55 and :69, which pin the column name.

## Impact

A reader taking the 1/C² column or the plot axis at face value mis-scales by 1/area (1e4 at the
default area) when extracting N_D from a Mott-Schottky slope. That is exactly the doping
extraction workflow this tab exists for.

## Proposed fix

```diff
# app/components/results.py
-        yaxis_title="Capacitance (F/cm²)",
+        yaxis_title="Capacitance (F)",
 ...
-        yaxis_title="1/C² (cm⁴/F²)",
+        yaxis_title="1/C² (F⁻²)",
 ...
-                "one_over_C2_cm4_per_F2": result.metadata["one_over_C_squared"],
+                "one_over_C2_per_F2": result.metadata["one_over_C_squared"],

# tests/test_app_csv_export.py:55, :69
-        == "bias_V,capacitance_F,one_over_C2_cm4_per_F2,depletion_width_cm"
+        == "bias_V,capacitance_F,one_over_C2_per_F2,depletion_width_cm"
```

## Test plan

- Update test_app_csv_export assertions as above.
- Add an assertion on `build_cv_figure(r).layout.yaxis.title.text == "Capacitance (F)"` and on the Mott-Schottky title.

## Effort

S

## Status 2026-10-18 (final-c): fixed (00ce048 + f88e0bb)

- 00ce048: figure axes "Capacitance (F)" and "1/C² (F⁻²)". f88e0bb: CSV column renamed to `one_over_C2_per_F2`.
- Evidence: `uv run pytest -q tests/test_app_csv_export.py` -> 22 passed (test_cv_csv_columns_and_header with the new header, test_cv_figure_axes_are_absolute_farads).
