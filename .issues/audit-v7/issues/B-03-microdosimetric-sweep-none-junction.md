---
id: B-03
title: microdosimetric_sweep fails every configuration when N_D_junction / L_transition are not passed (create_sic_2d_device stores raw None)
severity: high
category: bug
status: fixed  # 2026-10-18 sync: 845a66b: create_sic_2d_device returns the resolved N_D_junction/N_D_bulk/L_transition, never None (device2d.py 'audit B-03')
location: etna/core/optimization.py:116-144 (Mj-3 gate), 180-208 (broad except plus reset); etna/core/device2d.py:405-410 (resolves defaults), 462-464 (returns the raw args)
labels: [microdosimetry, 2d, optimization, silent-failure, raw-F3]
---

## Summary

`microdosimetric_sweep` calls `create_2d_dd_device(half_width_um, V_bias, epi_thickness_cm, N_D_bulk=nd)` without passing `N_D_junction` or `L_transition`. `create_sic_2d_device` resolves both to `_N_D_JUNCTION_DEFAULT` and `_L_TRANSITION_DEFAULT` for the doping it builds (device2d.py:407-409). It then returns the **raw arguments**, so `device_info["N_D_junction"]` and `["L_transition"]` are `None`. The audit Mj-3 full-depletion gate calls `built_in_potential(N_A_ionized, None, n_i)`, which raises `TypeError`. The broad `except Exception` catches it, logs "Sweep i/N failed", writes a row of NaN with `is_valid=False`, and calls `reset_devsim_fully()`, all after a DD solve and lateral scan that had succeeded. With its default call signature, which is also the one notebook 20 uses, the ranking feature cannot work.

## Reproduction / Evidence

- ran (scratch script, `uv run python`):
  ```
  d = optimization.create_2d_dd_device(half_width_um=25, V_bias=0.0, epi_thickness_cm=10e-4, N_D_bulk=8.5e13)
  d["N_D_junction"], d["L_transition"]   -> None None
  built_in_potential(d["N_A_ionized"], d["N_D_junction"], d["n_i"])
     -> TypeError '<=' not supported between instances of 'NoneType' and 'int'
  ```
- call chain: optimization.py:133-144 passes `device_info["N_D_junction"]` and `["L_transition"]` to `built_in_potential` and `full_depletion_voltage_graded`. Lines 180-208 contain `except Exception` and `reset_devsim_fully`.
- call chain: `scripts/create_notebook_20.py:194-200` passes only half_widths, epi_thicknesses, N_D_bulks and V_biases.
- From the raw audit, not re-run: a full `microdosimetric_sweep(...)` with one config returns `center_cce=nan, V_fd=nan, is_valid=False` and the same TypeError in the log.
- `tests/test_optimization.py` does not exercise the sweep.

## Impact

If notebook 20 were re-executed, the sweep would produce 0 valid rows (see B-15). The ranking and scoring functions are dead in practice, and the failure is buried in a WARNING log.

## Proposed fix

Store the resolved values at the source:

```diff
 # etna/core/device2d.py (return dict of create_sic_2d_device)
-        "N_D_junction": N_D_junction,
-        "N_D_bulk": N_D_bulk,
-        "L_transition": L_transition,
+        "N_D_junction": _N_D_junction,
+        "N_D_bulk": _N_D_bulk,
+        "L_transition": _L_transition,
```

Optionally, narrow the `except Exception` in `microdosimetric_sweep` to solver errors (`RuntimeError`, `devsim.error`) so a `TypeError` coding bug propagates instead of becoming a NaN row.

## Test plan

- Unit (needs devsim): `create_sic_2d_device(...)` with N_D_junction and L_transition omitted returns floats equal to the module defaults.
- Slow: `microdosimetric_sweep(half_widths_um=(25,), epi_thicknesses_cm=(10e-4,), N_D_bulks=(8.5e13,), V_biases=(40,), n_lateral_points=2)` gives `np.isfinite(df["V_fd"]).all()` and no "failed" in the caplog.

## Effort

S
