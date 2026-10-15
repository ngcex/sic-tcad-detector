---
id: A-09
title: delete_device leaves the mesh alive; per-call devices in cv_at_fluence, transient_cce_vs_dose_rate and API finally-blocks leak a mesh each
severity: low
category: improvement
status: open
location: etna/core/cv_analysis.py:~332-336, ~413-417; etna/core/transient.py:~644-648; etna/api/simulation.py:~133-141, ~326-334, ~758-767; contrast etna/core/device.py:~605-613 (calibrate_graded_doping deletes both)
labels: [audit-v7, slice-A, resource-leak, devsim]
source: raw-A-core-solvers.md F9
---

## Summary

devsim stores meshes separately from devices. After `devsim.delete_device(device=X)`, the
mesh `X_mesh` remains, and a new device can still be created from it. The per-call
devices use unique uuid names, so there are no collisions. However, long fluence or
dose-rate loops accumulate meshes for the lifetime of the process. The API entry points
call `reset_devsim_fully()` on entry, which cleans up only on the next call.

## Reproduction / Evidence

Command I ran:

```
create_sic_device('mp'); devsim.delete_device(device='mp')
devsim.create_device(mesh='mp_mesh', device='mp2')   -> succeeded ("mesh survived delete_device")
```

From the raw audit, not re-run: the specific `finally` blocks listed under location call
only `delete_device`. `calibrate_graded_doping` deletes both the device and the mesh.

## Impact

1D meshes have about 330 nodes, so the memory per call is negligible. Growth is unbounded
only in long loops such as notebooks 11, 12 and 14 and `run_radiation_damage`. This is a
hygiene and consistency issue.

## Proposed fix

```python
finally:
    for fn, kw in ((devsim.delete_device, {"device": name}),
                   (devsim.delete_mesh, {"mesh": f"{name}_mesh"})):
        try:
            fn(**kw)
        except devsim.error:
            pass
```

Better still, add `device_info["mesh_name"]` in `create_sic_device` and
`create_sic_device_2d`, plus a helper `delete_device_and_mesh(device_info)` that every
`finally` uses.

## Test plan

After `cv_at_fluence(0, [0, -5])` (using a patched uuid, or by reading the name from a
log), `devsim.create_device(mesh="<name>_mesh", device="x")` must raise `devsim.error`.
Add the same check for one API function.

## Effort

S
