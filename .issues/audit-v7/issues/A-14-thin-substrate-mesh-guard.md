---
id: A-14
title: 1D mesh builder has no guard for thin substrates (negative mesh line -> opaque devsim error) while the 2D builder does; DeviceConfig unvalidated
severity: low
category: improvement
status: open
location: etna/core/device.py:108 (pos=junction_pos - 5e-6 unconditional) vs etna/core/device2d.py:244-247 (guarded); etna/api/device.py:30 (substrate_thickness_um unvalidated); app/components/device_sidebar.py:80-85 (min_value=0.0)
labels: [audit-v7, slice-A, input-validation, mesh]
source: raw-A-core-solvers.md F14
---

## Summary

`create_sic_device` always adds a mesh line at `junction_pos - 5e-6`. When
`substrate_thickness_cm < 5e-6` (0.05 µm), that position is negative or out of order. The
2D builder guards the same line with `if junction_pos - 5e-6 > 1e-8`. `DeviceConfig`
validates no field, and the sidebar allows `min_value=0.0` for substrate thickness.

The raw audit marked this "plausible" because it was never executed. I executed it and
upgraded it to confirmed, so the status is open. The failure is loud, not silent, but the
error message is opaque.

## Reproduction / Evidence

Command I ran:

```
create_sic_device('thin', substrate_thickness_cm=2e-6)
-> devsim error: "Contacts must be at both ends of 1d device"
```

I read the code at device.py:108 (no guard), device2d.py:244-247 (guarded) and
device_sidebar.py:80-85 (`min_value=0.0`).

## Impact

A user who types a substrate thickness of 0.04 µm or less in the sidebar, or passes it
through `DeviceConfig`, gets an unintelligible devsim error instead of a validation
message. No silent wrong results.

## Proposed fix

```diff
 # etna/core/device.py
-    devsim.add_1d_mesh_line(mesh=mesh_name, pos=junction_pos - 5e-6, ps=1e-6)
+    if junction_pos - 5e-6 > 1e-8:
+        devsim.add_1d_mesh_line(mesh=mesh_name, pos=junction_pos - 5e-6, ps=1e-6)
```

```python
# etna/api/device.py
def __post_init__(self):
    for f in ("epi_thickness_um", "substrate_thickness_um", "area_cm2"):
        if getattr(self, f) <= 0:
            raise ValueError(f"{f} must be > 0")
    if self.half_width_um is not None and self.half_width_um <= 0:
        raise ValueError("half_width_um must be > 0")
```

Also raise the sidebar `min_value` to 0.1.

## Test plan

- `create_sic_device(substrate_thickness_cm=2e-6)` builds, and the minimum node position
  is 0.
- `DeviceConfig(epi_thickness_um=-1)` raises `ValueError`.
- `DeviceConfig(substrate_thickness_um=0)` raises `ValueError`.

## Effort

S
