---
id: C-12
source: raw-C F12
title: DeviceConfig has no validation; an unknown doping_profile silently becomes uniform, and uniform with N_D=None silently uses a hidden default
severity: low
category: bug
status: open
location:
  - etna/api/device.py:25-43 (DeviceConfig, no __post_init__)
  - etna/core/device.py:92-93 (`if N_D is None: N_D = DEFAULT_N_D`, DEFAULT_N_D = 1.07e15 at :33)
  - etna/core/device.py:228-247 (`if doping_profile == "graded": ... else: uniform`)
labels: [api, input-validation, public-contract]
---

## Summary

- `DeviceConfig(doping_profile="Graded")`, `"step"` or any typo builds a uniform device with no error.
- `DeviceConfig(doping_profile="uniform")` with N_D left at None builds with the hidden `DEFAULT_N_D = 1.07e15`.
- Geometry, doping, T and area can be zero, negative or non-finite.

DeviceConfig is the public contract and is serialized into every CSV header.

Line-number correction: the raw file cites core/device.py:51-52 and :187-206. The current
locations are :92-93 and :228-247.

Scope note: the app always supplies N_D when "uniform" is selected (device_sidebar.py:100-103)
and restricts the profile to two strings. The N_D=None and typo cases therefore affect library
and notebook users only. The zero/negative cases affect both (C-03, C-09).

## Reproduction / Evidence

```
uv run python -c "from etna.core.device import create_sic_device; \
  d = create_sic_device(device_name='x1', doping_profile='Graded'); print('built ok', d['doping_profile'])"
-> built ok with Graded: Graded          # no error; took the uniform else-branch
grep -n "__post_init__\|raise" etna/api/device.py  -> (no output)
```

## Impact

Typos produce a physically different device that is labeled with the user's intended name.
Missing N_D silently substitutes a value the user never chose.

## Proposed fix

```diff
 @dataclass
 class DeviceConfig:
     ...
+    def __post_init__(self):
+        if self.doping_profile not in ("graded", "uniform"):
+            raise ValueError(f"doping_profile must be 'graded' or 'uniform', got {self.doping_profile!r}")
+        if self.doping_profile == "uniform" and self.N_D is None:
+            raise ValueError("doping_profile='uniform' requires N_D")
+        for name in ("epi_thickness_um", "substrate_thickness_um", "N_A", "N_D_junction",
+                     "N_D_bulk", "T", "area_cm2"):
+            v = getattr(self, name)
+            if not (math.isfinite(v) and v > 0):
+                raise ValueError(f"{name} must be finite and > 0, got {v!r}")
+        if not (math.isfinite(self.L_transition_um) and self.L_transition_um >= 0):
+            raise ValueError("L_transition_um must be finite and >= 0")
+        if self.half_width_um is not None and not self.half_width_um > 0:
+            raise ValueError("half_width_um must be > 0 or None")
```

`dataclasses.replace` (used by ParametricSweep) re-runs `__post_init__`, so sweeps are covered
automatically.

## Test plan

A parametrized `pytest.raises(ValueError)` case for each rule. Also check that `DeviceConfig()`
and the sidebar's `assemble_config` outputs still construct, and that the existing
test_api_* / test_app_* suites pass.

## Effort

S
