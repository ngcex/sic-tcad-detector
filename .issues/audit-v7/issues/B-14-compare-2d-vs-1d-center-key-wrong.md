---
id: B-14
title: compare_cce_2d_vs_1d returns "cce_2d_center" equal to the full-area CCE, not the center CCE
severity: low
category: docs
status: fixed 2026-10-18 (wp-c2 3d7a375: cce_2d_center = centre-column CCE; equal to the full-area value to ~1e-9 on the laterally uniform device, notebook 15 output unchanged at 4 decimals)
location: etna/core/charge_collection_2d.py:497-515 (compare_cce_2d_vs_1d return dict); notebooks/15_2d_electrostatics_cce.ipynb; scripts/create_notebook_15*.py
labels: [2d, cce, wrong-key, raw-F14]
---

## Summary

`"cce_2d_center"` and `"cce_2d_full_area"` are both set to `float(cce_2d_full)`, the CCE for uniform generation over the full area. The true center CCE, `lateral_scan["cce_values"][0]`, is already computed as the local `cce_center` (line ~498) but is never exposed under the "center" key. Notebook 15 prints this key as "2D center".

The raw audit filed this under "docs". It is really a wrong-value bug, but low severity, because the center value can still be read from `lateral_scan`.

## Reproduction / Evidence

- ran: `sed -n 470,515p etna/core/charge_collection_2d.py`:
  - `cce_center = cce_lat[0]`
  - `"cce_2d_center": float(cce_2d_full)`
  - `"cce_2d_full_area": float(cce_2d_full)`
- ran: `grep -rln cce_2d_center ...` finds charge_collection_2d.py, test_charge_collection_2d.py, create_notebook_15.py, create_notebook_15_v2.py and notebooks/15_2d_electrostatics_cce.ipynb.

## Impact

Notebook 15 shows the full-area value labelled as the center value. The "2D center vs 1D" comparison is therefore not the comparison it claims to be.

## Proposed fix

```diff
-            "cce_2d_center": float(cce_2d_full),
+            "cce_2d_center": float(cce_center),
             "cce_2d_full_area": float(cce_2d_full),
```

Re-execute notebook 15, or at least regenerate the affected cell.

## Test plan

Assert `result["cce_2d_center"] == pytest.approx(result["lateral_scan"]["cce_values"][0])` in `test_charge_collection_2d.py`. Check whether existing assertions on `cce_2d_center` assumed the full-area value.

## Effort

S
