---
id: B-01
title: sweep_iv_vs_temperature(use_tat=True) ramps with a frozen zero-bias E_field_node, so leakage comes out ~7.6x too low
severity: high
category: bug
status: fixed  # 2026-10-18 sync: a7eea2b: sweep_iv_vs_temperature refreshes TAT field/gamma per bias step (PERF-02 section 7)
location: etna/core/temperature_sweep.py:88-102; etna/core/dark_current.py:141,210,233-265; etna/core/drift_diffusion.py ramp_bias (~l.245-310)
labels: [physics, dark-current, temperature, tat, silent-wrong-result, raw-F1]
---

## Summary

`setup_tat_model` builds `U_TAT = ... - G0_TAT * Gamma_n * min(E_field_node / E_ref_TAT, 1.0)`. `E_field_node` is **not** a devsim expression. It is a static node array written once by `_compute_node_efield()` (`devsim.node_model(..., equation="1e3")` followed by `set_node_values`). All sweeps inside `dark_current.py` call `_compute_node_efield` and `_compute_gamma_factors` again at every bias step. `temperature_sweep.sweep_iv_vs_temperature` does not. It calls `create_dark_current_device()`, which snapshots the field at 0 V, and then `drift_diffusion.ramp_bias()`, which never refreshes that array. The depletion indicator `f_depl` therefore stays at its 0 V width while the real depletion region grows, and the generation current comes out several times too low.

## Reproduction / Evidence

- ran: `grep -n "E_field_node" etna/core/drift_diffusion.py`. No hits: ramp_bias only sets the bias parameter and solves.
- ran: scratch script, `uv run python`, using the graded device as in temperature_sweep (N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1e-4, T=300):
  - `create_dark_current_device(...)`, then `ramp_bias(d, 30.0, contact="cathode")`, then `np.array_equal(E_field_node_before, E_field_node_after)` gives **True**, with I = 3.698e-10 A/cm^2
  - `sweep_iv_vs_temperature([300.0], V_reverse=-30)["I_reverse"]` = **3.698e-10 A/cm^2**
  - `dark_current_sweep(create_dark_current_device(...), [-30.0], area=0.05)["J_total"]` = **2.825e-09 A/cm^2**
  - Ratio: **7.64x**.
- From the raw audit, not re-run: the integrated f_depl is 1.065e-4 cm stale vs 8.148e-4 cm refreshed, the same 7.6x.
- call chain: `tests/test_temperature_sweep.py::TestSweepIV` checks only the columns and I(310K) > I(290K). The n_i(T) scale factor keeps the result monotonic, so the bug never shows up in the suite.
- Consumers: notebook 06 (`sweep_iv_vs_temperature(T_clinical, V_reverse=-30)`) and tests. The Streamlit dark-current page uses `dark_current_sweep`, which is correct.

## Impact

The absolute temperature-dependent leakage I_reverse(T) from the "C8" path, including notebook 06's leakage coefficient dI/dT, is not on the same footing as the 300 K calibration it claims to preserve. It is roughly 7.6x low at -30 V. The error depends on bias, so dI/dT is also biased.

## Proposed fix

Minimal option: reuse the path that refreshes the field.

```diff
-            cathode_V = -V_reverse
-            ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5)
-            I_total = extract_contact_current(device_info, contact="cathode")
+            if use_tat:
+                res = dark_current_sweep(device_info, [V_reverse], area=1.0, V_step=0.5)
+                I_total = float(res["J_total"][-1])
+            else:
+                ramp_bias(device_info, -V_reverse, contact="cathode", V_step=0.5)
+                I_total = extract_contact_current(device_info, contact="cathode")
```

Structural option (recommended; it also fixes B-02 and the lag within each step): make `E_field_node` a real devsim model. For example, define an edge model `absE = abs(ElectricField)`, then `devsim.edge_average_model(device, region, edge_model="absE", node_model="E_field_node", average_type="arithmetic_mean")`. Either register the `Potential` derivatives or accept a lag within each step, and drop the static `set_node_values`.

## Test plan

- New slow test: `sweep_iv_vs_temperature([300], V_reverse=-30)["I_reverse"]` must match `dark_current_sweep(create_dark_current_device(same kwargs), [-30])["J_total"]` within 5%.
- Regression: after any ramp in the TAT path, assert that `E_field_node` differs from its 0 V snapshot.

## Effort

M
