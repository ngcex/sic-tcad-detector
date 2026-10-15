---
id: A-08
title: devsim_reset.py rationale is false for devsim 2.10 (reset_devsim() does clear globals) and docs cite stale src/ paths
severity: low
category: docs
status: open
location: etna/core/devsim_reset.py:1-25 (header), 35, 48, 63-66, step 3 (~111-117); etna/core/device2d.py:484, 509; tests/test_devsim_reset.py:1-5
labels: [audit-v7, slice-A, docs, devsim]
source: raw-A-core-solvers.md F8
related: A-15
---

## Summary

The module header says: "devsim's own reset_devsim() deletes devices but does NOT unset
the process-global parameters". On the pinned devsim 2.10.0, `reset_devsim()` behaves as
follows:

- It removes custom globals such as `raxis_zero`.
- It restores `node_volume_model` to "NodeVolume".

Step 3 of `reset_devsim_fully` writes the Cartesian defaults before step 4 calls
`reset_devsim()`, and step 4 resets those same values, so step 3 is redundant. The real
value of the wrapper is twofold:

- It enumerates and deletes devices.
- It preserves `direct_solver` and `solver_callback` across the reset (steps 1 and 5).

The docs also reference `src/...` paths from a package layout that no longer exists.

## Reproduction / Evidence

Command I ran:

```
devsim.set_parameter(name='raxis_zero', value=0.0)
devsim.set_parameter(name='node_volume_model', value='CylindricalNodeVolume')
devsim.reset_devsim()
get_parameter('raxis_zero')        -> error: Cannot find parameter "raxis_zero"
get_parameter('node_volume_model') -> 'NodeVolume'
```

Command I ran: `grep -n "src/" etna/core/devsim_reset.py etna/core/device2d.py`.

- devsim_reset.py: 5, 18, 25, 35, 48, 64, 65
- device2d.py: 484
- Also device2d.py:509, ":func:`src.devsim_reset.reset_devsim_fully`" (a dotted form the
  grep missed)

## Impact

The code is harmless, but the rationale will mislead maintainers. For example, it was
used to justify the per-trial full reset in A-15. The stale paths break navigation. No
numerical impact.

## Proposed fix

- Rewrite the header along these lines: "belt-and-braces wrapper: enumerates and deletes
  devices and preserves solver parameters that `reset_devsim()` drops. devsim ≥ 2.10
  already clears custom globals and restores the default assembly models. Step 3 is kept
  as defence in case a future devsim regresses."
- Alternatively, move step 3 after step 4.
- Replace `src/` with `etna/core/` and `src.devsim_reset` with
  `etna.core.devsim_reset` in both files.

## Test plan

Add `test_plain_reset_devsim_clears_cylindrical_globals`:

- Set `node_volume_model` to cylindrical and `raxis_zero`.
- Call plain `devsim.reset_devsim()`.
- Assert `node_volume_model == "NodeVolume"` and that `raxis_zero` is absent.

This documents the devsim behaviour the module relies on and flags any devsim upgrade
that changes it. Keep `TestResetStateLeak`.

## Effort

S
