---
id: B-05
title: cce_vs_bias_at_fluence loses every converged point on a mid-sweep failure; sweep_cce_vs_temperature("dd") hides truncation
severity: medium
category: bug
status: fixed 2026-10-18 (final-abd 317d2c4: cce_vs_bias_at_fluence keeps the converged prefix with truncated / requested_v_stop; sweep_cce_vs_temperature hecht and dd pair CCE with the returned voltages and add a truncated column. Evidence: pytest tests/test_truncation_b05.py -> 3 passed)
location: etna/core/charge_collection.py:819-1015 (cce_vs_bias_at_fluence, try/finally with no except at ~973-1010); etna/core/temperature_sweep.py:303-309; etna/api/simulation.py:528 (run_temperature_sweep method param)
labels: [cce, radiation-damage, temperature, truncation, raw-F5]
---

## Summary

1. `cce_vs_bias` catches the `RuntimeError` from `ramp_bias` and returns the converged subset with `truncated=True`. `cce_vs_bias_at_fluence` has only `try/finally`. The first bias it cannot reach raises, `finally` deletes the device, and every CCE value computed so far is lost. Damaged devices are where this is most likely, because carrier removal lowers the effective N_D and moves punch-through to lower |V|.
2. `sweep_cce_vs_temperature(method="dd")` loops `for j, V in enumerate(voltages): records.append({... result["cce_values"][j]})`. When `cce_vs_bias` truncates, `cce_values` is shorter than `voltages`. The rows for the leading converged indices are appended. Then `IndexError` is raised, caught by `except Exception`, and logged only as a warning. **Correction to the raw finding:** the raw text says the whole temperature point is dropped. For voltages in ramp order, the converged prefix is in fact kept and the remaining voltages are silently missing, with no `truncated` flag. For voltages not in ramp order, the index `j` into `voltages` does not line up with the reached-mask order, so CCE values could be assigned to the wrong voltage.

## Reproduction / Evidence

- call chain: `ramp_bias` raises `RuntimeError("ramp_bias: failed to converge ...")` (drift_diffusion.py, ramp_bias raise at ~l.310). In `cce_vs_bias_at_fluence`, `grep` over lines 960-1015 shows `try:` at the loop, `finally:`, and no `except`.
- call chain: `cce_vs_bias` returns `V_range[reached_mask]` and `cce_sorted[reached_mask]` (charge_collection.py:616-623). temperature_sweep.py:303-309 indexes by the requested `voltages`.
- ran: `grep -l cce_vs_bias_at_fluence notebooks/*.ipynb` lists `10_cce_vs_fluence.ipynb` and `14_validation.ipynb`.
- ran: `grep -n "def run_temperature_sweep" -A30 etna/api/simulation.py` shows `method: str = "hecht"` and forwarding to `sweep_cce_vs_temperature`, so `"dd"` is reachable through the API.
- The earlier "truncation reported, not silent" fix covered `cv_sweep` and `cce_vs_bias` only.

## Impact

Fluence studies in notebooks 10 and 14 fail outright, discarding valid points, when one bias is unreachable. Temperature DD sweeps silently return fewer rows than requested, with no flag.

## Proposed fix

```diff
 # cce_vs_bias_at_fluence
+    reached_mask = np.zeros(len(V_range), dtype=bool)
+    truncated = False
     try:
         for idx, V_target in zip(sorted_indices, sorted_V):
-            ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5)
+            try:
+                ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5)
+            except RuntimeError as exc:
+                logger.warning("cce_vs_bias_at_fluence truncated at %s V: %s", V_target, exc)
+                truncated = True
+                break
             ...
+            reached_mask[idx] = True
     finally: ...
-    return {"voltages": V_range, "cce_values": cce, ...}
+    return {"voltages": V_range[reached_mask], "cce_values": cce[reached_mask], ...,
+            "truncated": truncated, "requested_v_stop": float(V_range.min())}
```

```diff
 # sweep_cce_vs_temperature, dd branch
-                for j, V in enumerate(voltages):
-                    records.append({"T": T, "V": V, "CCE": result["cce_values"][j]})
+                for V, c in zip(result["voltages"], result["cce_values"]):
+                    records.append({"T": T, "V": V, "CCE": c,
+                                    "truncated": result.get("truncated", False)})
```

## Test plan

- Monkeypatch `ramp_bias` to raise `RuntimeError` on the second target. Assert that `cce_vs_bias_at_fluence` returns one point with `truncated=True`.
- Monkeypatch `cce_vs_bias` to return a truncated dict. Assert that the DataFrame from `sweep_cce_vs_temperature(method="dd")` has rows only for the reached voltages, `truncated=True`, and each V paired with the right CCE.

## Effort

S
