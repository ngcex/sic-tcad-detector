---
id: A-15
title: calibrate_graded_doping_2d destroys every devsim device in the process on each trial (undocumented at the call site)
severity: low
category: improvement
status: fixed 2026-10-18 (WP-A 93944c4: Cartesian globals restored once at entry, each trial deletes only its own device and mesh)
location: etna/core/device2d.py:552 (import), ~664-665 (finally: reset_devsim_fully(preserve_solver=True) inside objective), 509-512 (docstring)
labels: [audit-v7, slice-A, side-effect, devsim]
source: raw-A-core-solvers.md F15
related: A-08
---

## Summary

The 2D calibration objective calls `reset_devsim_fully(preserve_solver=True)` in a
`finally` block on every trial, roughly 80 times per run. `reset_devsim_fully` enumerates
and deletes every device in the process. Any caller-held device is therefore silently
deleted after the first trial, for example a 1D twin kept for comparison or an app
session device.

The docstring says the reset clears "all devsim global state", but it does not warn that
caller-owned devices are destroyed. `run_cv` documents the same hazard explicitly. The
docstring also cites the stale `src.devsim_reset` path (see A-08). By contrast, the 1D
calibrator deletes only its own trial device and mesh.

## Reproduction / Evidence

- `grep -n reset_devsim_fully etna/core/device2d.py` returns 509 (docstring), 552
  (import) and 665 (`finally` in the objective).
- I read device2d.py:655-668 and confirmed the call runs unconditionally on every trial.

From the raw audit, not re-run: `scripts/run_calibration_2d.py` is the only caller.

## Impact

This is a surprising side effect for library users. A full reset per trial is also
unnecessary: per A-08, `devsim.reset_devsim()` already clears the cylindrical globals,
and deleting the trial device and mesh is enough between planar trials.

## Proposed fix

```python
def calibrate_graded_doping_2d(...):
    """...
    Warning
    -------
    Calls reset_devsim_fully() once on entry: ALL existing devsim devices
    in the process are deleted. Build reference devices after calibration.
    """
    reset_devsim_fully(preserve_solver=True)          # once, at entry
    def objective(p):
        ...
        finally:
            for fn, kw in ((devsim.delete_device, {"device": trial_name}),
                           (devsim.delete_mesh, {"mesh": f"{trial_name}_mesh"})):
                try: fn(**kw)
                except devsim.error: pass
```

## Test plan

Patch the objective to run 2 trials (maxiter=1 with a small simplex, or call `objective`
directly). Then check:

- a sentinel device created after the entry reset survives both trials
- `devsim.get_device_list()` holds only the sentinel afterwards
- the calibration result is unchanged versus the current implementation for a fixed seed
  and starting point

## Effort

S
