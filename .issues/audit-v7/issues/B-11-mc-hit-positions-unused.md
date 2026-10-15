---
id: B-11
title: MC event hit positions are loaded and unit-converted but never used; per-event CCE ignores lateral position (edge effects absent from the spectrum)
severity: low
category: improvement
status: open
location: etna/core/mc_coupling.py:92-99 (convert_units), 314-384 (events_to_charge_profiles), 392-462 (process_mc_ensemble, LET warning at ~l.440); etna/core/charge_collection_2d.py:203+ (cce_lateral_scan), 326-394 (cce_heatmap_2d)
labels: [microdosimetry, mc-coupling, 2d, improvement, raw-F11]
---

## Summary

Both the "full path" (`events_to_charge_profiles`) and the "fast path" (`process_mc_ensemble`) reduce each event to `LET = sum(edep) / sv_thickness` and apply a CCE that does not depend on position. The loaders require `x_cm`, `y_cm` and `z_cm`, and `convert_units` converts them, but nothing reads them afterwards. The module docstring ("Per-event mesh charge profile mapping") and the Phase-20 lateral CCE map suggest that edge effects reach the microdosimetric spectrum. They do not. Also, `process_mc_ensemble` warns only for LET outside [0.01, 1e4] keV/um, while the default CCE(LET) table spans roughly [0.1, 1000] and `np.interp` clamps silently outside the table.

## Reproduction / Evidence

- ran: `grep -n "x_cm\|y_cm\|z_cm" etna/core/mc_coupling.py etna/core/microdosimetry.py`. Every hit is in a loader column mapping, `convert_units` (l.92) or a docstring. None is a consumer.
- ran: `grep -n "cce_heatmap_2d\|cce_lateral\|x_positions_cm" etna/core/mc_coupling.py` returns nothing.
- ran: `grep -n "0.01\|1e4\|warn" etna/core/mc_coupling.py` shows line 440, `if np.any(LET < 0.01) or np.any(LET > 10000)`.
- The CCE table range [0.1, 1000] comes from the raw audit and was not re-checked.

## Impact

Microdosimetric spectra from MC input leave out the lateral and edge-collection losses that the 2D work characterised. For small SVs, where the edge fraction is large, y_F and y_D are biased high. The current docs overstate what the coupling does.

## Proposed fix

```diff
-def process_mc_ensemble(events_df, cce_table, sv_thickness_um, ...):
+def process_mc_ensemble(events_df, cce_table, sv_thickness_um, ..., cce_lateral=None):
     ...
     cce = np.interp(LET, table_let, table_cce)
+    lo, hi = table_let.min(), table_let.max()
+    if np.any((LET < lo) | (LET > hi)):
+        logger.warning("LET outside CCE table [%g, %g]: clamped", lo, hi)
+    if cce_lateral is not None:
+        x_evt = per_event_edep_weighted(events_df, "x_cm")  # |x| centroid
+        f_edge = cce_lateral["cce_values"] / cce_lateral["cce_values"][0]
+        cce = cce * np.interp(np.abs(x_evt), cce_lateral["x_positions_cm"], f_edge)
```

Update the module docstring to describe current behaviour if this is not implemented.

## Test plan

Two synthetic events with the same edep, one at x=0 and one at x=edge, with a lateral map in which the edge CCE is half the centre CCE. Assert the collected energies differ by 2x. Add a caplog test for the out-of-table warning.

## Effort

M
