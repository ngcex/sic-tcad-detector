---
id: B-02
title: _compute_node_efield assumes a 1D edge chain; on 2D meshes it gives nodes the field of unrelated edges (get_dark_current_2d / notebook 20 noise floor)
severity: high
category: bug
status: fixed  # 2026-10-18 sync: 9d19d9c: TAT node E-field on 2D meshes uses |grad Potential| (dimension branch, dark_current.py get_dimension == 2)
location: etna/core/dark_current.py:233-265 (_compute_node_efield), 94-230 (setup_tat_model); etna/core/optimization.py:380-433 (get_dark_current_2d); scripts/create_notebook_20.py:428
labels: [physics, dark-current, 2d, tat, microdosimetry, raw-F2]
---

## Summary

For node `i`, `_compute_node_efield` averages `abs_E[i-1]` and `abs_E[i]` taken from the `ElectricField` edge model. That only makes sense for a 1D chain with `n_edges == n_nodes - 1` and edges stored in node order. On a 2D triangular mesh there are about 3x more edges than nodes, and edge index has no positional relation to node index. `E_field_node`, and therefore `f_depl = min(E_field_node / E_ref, 1)` in `setup_tat_model`, is meaningless. The loop also never reads edges past index `n_nodes - 1`. Neither function has a `dimension` guard. `optimization.get_dark_current_2d` calls `setup_tat_model` on a 2D device and feeds the resulting current into `estimate_noise_floor`, which notebook 20 uses for its noise floor and structure scoring.

## Reproduction / Evidence

- ran: `create_2d_dd_device(half_width_um=25, V_bias=20)` followed by `setup_tat_model(d2)` gives **n_nodes = 3924, n_edges = 11095**. `n_edges != n_nodes - 1`, so the 1D index formula cannot hold.
- call chain: dark_current.py:252-258 contains the loop `E_node[i] = 0.5*(abs_E[i-1] + abs_E[i])`. `grep -n dimension etna/core/dark_current.py` has no hits.
- call chain: `scripts/create_notebook_20.py:80` imports `get_dark_current_2d` and line 428 calls `I_dark = get_dark_current_2d(...)`. `notebooks/20_feasibility_report.ipynb` has 4 occurrences.
- From the raw audit, not re-run: the Pearson correlation between log10(E_field_node) and the node-averaged triangle |grad V| is 0.056. Mean f_depl in the epi is 0.379 from this function vs 0.776 from the true gradient, and max E is 4.76e4 vs 9.34e4 V/cm.
- `tests/test_optimization.py` does not exercise `get_dark_current_2d`.

## Impact

The 2D dark current, and with it the noise floor and microdosimetry structure scores in notebook 20, is computed from an arbitrary generation volume. The numbers are not physically meaningful, and the error is silent.

## Proposed fix

Use a primitive that works in any dimension, for example:

```diff
 def _compute_node_efield(device, region):
-    E_edge = np.array(devsim.get_edge_model_values(..., name="ElectricField"))
-    ... 1D index loop ...
-    devsim.node_model(device=device, region=region, name="E_field_node", equation="1e3")
-    devsim.set_node_values(...)
+    devsim.edge_model(device=device, region=region, name="absE_edge",
+                      equation="max(abs(ElectricField), 1e3)")
+    devsim.edge_average_model(device=device, region=region, edge_model="absE_edge",
+                              node_model="E_field_node", average_type="arithmetic_mean")
```

For an exact 2D |E|, use `element_from_edge_model` to get Ex and Ey, then `|E| = sqrt(Ex^2 + Ey^2)`. As a minimum stopgap, `raise NotImplementedError` when `device_info.get("dimension") == 2`. That requires passing device_info, or checking the edge count, e.g. `n_edges != n_nodes - 1`.

## Test plan

- 2D: on `create_2d_dd_device(25, 20)`, assert corr(E_field_node, node-averaged triangle-gradient |E|) > 0.95.
- 1D: assert `E_field_node` is unchanged within 1e-6 relative on the default 1D device, so the calibration is preserved.
- Smoke test that `get_dark_current_2d` returns a finite value of plausible magnitude.

## Effort

M
