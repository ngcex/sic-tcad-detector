---
id: C-13
source: raw-C F13
title: Dead or test-only surface in the API layer (MeshData.regions/contacts, build_device device_name, run_flash_recombination/run_transient not exportable)
severity: info
category: dead-code
status: fixed
location:
  - etna/api/results.py:32-33 and etna/api/simulation.py (run_field populates MeshData.regions / .contacts)
  - etna/api/device.py:46 (`build_device(..., device_name=None)`)
  - etna/api/simulation.py:592-662 (run_flash_recombination), 770-840 (run_transient)
  - app/components/results.py:226-305 (to_csv_bytes has no "temperature" / "flash" / "transient" branch)
labels: [api, dead-code, cleanup]
---

## Summary

- (a) `MeshData.regions` and `MeshData.contacts` are computed in run_field but read only by tests/test_api_field.py:47-53. The geometry viewer (app/components/geometry_viewer.py) does not reference them, although results.py describes them as the "geometry viewer contract".
- (b) `build_device`'s `device_name` kwarg is never passed by any caller in etna, app, tests or scripts.
- (c) `run_flash_recombination` and `run_transient` are exported from `etna` but referenced only in `etna/__init__.py`, `etna/api/simulation.py` and `tests/test_api_facades.py` (import/signature smoke). `to_csv_bytes` raises `ValueError: unknown sim_type` for their results and for `run_temperature_sweep`'s, so 3 of the 9 public facades cannot be exported with the shared serializer.

## Reproduction / Evidence

```
grep -rn "\.regions\|\.contacts" app etna tests   -> tests/test_api_field.py:47,48,53 only
grep -n "regions\|contacts" app/components/geometry_viewer.py  -> (none)
grep -rn "build_device(" etna app tests scripts     -> only positional single-arg calls (simulation.py:104, :204)
grep -rln "run_flash_recombination\|run_transient" etna app scripts tests notebooks
   -> etna/__init__.py, etna/api/simulation.py, tests/test_api_facades.py
sed -n 226,305p app/components/results.py  -> branches: cv, cce, field, damage, dark_current, microdosimetry; else raise
```

## Impact

Maintenance cost and a misleading "contract" docstring. Library users cannot serialize 3 facade
results with the shared helper. There is no correctness impact.

## Proposed fix

- Either draw region bounds and contacts in geometry_viewer (which makes the fields live), or drop them from MeshData and run_field and update test_api_field.
- Add to_csv_bytes branches:

```diff
+    elif result.sim_type in ("temperature", "flash", "transient"):
+        xname = {"temperature": "T_K", "flash": "dose_rate_Gy_s", "transient": "dose_rate_Gy_s"}[result.sim_type]
+        df = pd.DataFrame({xname: result.x, "cce": result.y})
+        extra_header_lines = [f"# {k}: {v}" for k, v in result.metadata.items() if np.isscalar(v)]
```

- Keep `device_name` (it is harmless, and useful for deterministic names in tests) but add one test that uses it, or remove it.

## Test plan

`to_csv_bytes` round-trip tests for the temperature, flash and transient sim_types, using
synthetic SimResults (fast, no devsim). If regions/contacts are dropped, update
test_api_field accordingly.

## Effort

S

## Status 2026-10-18 (final-c): fixed (f88e0bb + 9122689)

- (a) 9122689: MeshData.regions/.contacts are LIVE: app/components/geometry_viewer.py draws each contact as a dashed line at its depth (vertical in the 1D bar chart, horizontal on the reversed 2D depth axis) and outlines 2D regions; empty lists draw nothing. MeshData docstring rewritten (units, dict keys, devsim and engine2d producers). Live check: run_field 1D -> lines at 0 and 11.3 um; 2D (engine2d) -> rect + lines at 0 and 11.3 um.
- (b) build_device(device_name=...) kept; 9122689 adds tests/test_api_device.py::test_build_device_forwards_explicit_device_name (1D and 2D, no devsim).
- (c) f88e0bb: to_csv_bytes branches for "temperature" (T_K, CCE, solver from solver_per_T; voltage_V/method header), "flash" (dose_rate_Gy_s, CCE; V_bias/E_MeV + scope WARNING) and "transient" (dose_rate_Gy_s, CCE_clipped_0_1, CCE_raw from transient_cce_raw; current_contact header).
- Evidence: `uv run pytest -q tests/test_app_geometry_viewer.py tests/test_app_csv_export.py tests/test_api_device.py -k "not slow"` pass (3 viewer tests, 3 round-trip tests, 1 device_name test).
