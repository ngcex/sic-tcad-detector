---
id: A-03
title: Scalar mobilities are evaluated at the other region's doping (epi hole mobility at N_A=1e19; substrate electron mobility at N_D~1e15)
severity: medium
category: bug
status: partial (hole half fixed 2026-10-17, feat/cv-observable; electron mobility in the p+ and graded local doping still open)
location: etna/core/device.py:184-185, etna/core/device2d.py:362-363 (same pattern in etna/core/alternative_structures.py:103-104, outside slice); consumer etna/core/drift_diffusion.py (CreateElectronCurrent/CreateHoleCurrent reading mu_n/mu_p)
labels: [audit-v7, slice-A, physics, mobility, transient]
source: raw-A-core-solvers.md F3
---

## Summary

`mu_n = mobility_caughey_thomas_T(N_D, ...)` and `mu_p = mobility_caughey_thomas_T(N_A, ...)`
are each set as one region-wide devsim parameter. This produces two errors:

- **Holes in the n-epi:** holes are the collected minority carrier in the sensitive
  volume, but their mobility is evaluated at the substrate doping of 1e19.
- **Electrons in the p+ substrate:** their mobility is evaluated at the epi doping.

For graded profiles, `N_D` is the uniform fallback of 1.07e15, not the local doping.

## Reproduction / Evidence

Command I ran: `uv run python -c "from etna.core.sic_material import mobility_caughey_thomas_T as m; ..."`

- `m(1e19,300,'hole') = 75.7` cm²/Vs (used now). Compare `m(1e15,300,'hole') = 121.2`, which is 38 % higher.
- `m(1.07e15,300,'electron') = 913.4` (used now in the p+ substrate). Compare `m(1e19,300,'electron') = 115.3`, which is about 8x lower.

I read the code at device.py:183-185 and device2d.py:361-363 and confirmed both pass
N_D/N_A into one scalar per region.

From the raw audit, not re-run:

- No `node_model` or `edge_model` named mu_n/mu_p exists.
- `charge_collection.py:115-119` (Hecht path) uses `mu_p_max`=125 for the same epi.
- `tests/test_material.py::test_device_mobility_at_300k` pins the current values.

## Impact

The error affects:

- hole transit time and drift in transient pulses
- high-injection and ambipolar behaviour in the FLASH module
- the minority-diffusion dark-current component

It also leaves the DD path inconsistent with the Hecht CCE path. Steady-state CCE at full
depletion is largely insensitive, which is why the existing tests do not catch it.

## Proposed fix

**Full fix:** make mobility a node model driven by local doping, and edge-average it for
the current models.

```python
expr_n = f"{mu_min_n} + ({mu_max_n_T}-{mu_min_n})/(1+pow((Donors+Acceptors)/{N_ref_n},{alpha_n}))"
devsim.node_model(device=dev, region=reg, name="mu_n_node", equation=expr_n)
devsim.edge_average_model(device=dev, region=reg, node_model="mu_n_node", edge_model="mu_n")
# same for mu_p; CreateElectronCurrent/CreateHoleCurrent keep referencing "mu_n"/"mu_p"
```

Check that the simple_dd expressions resolve `mu_n` as an edge model, taking precedence
over the parameter.

**Minimal fix:** evaluate both mu_n and mu_p at the epi doping (N_D, or N_D_bulk for
graded profiles) and document that substrate values are approximate.

## Test plan

- On the default graded device, the effective hole mobility on epi nodes is in the range
  118-125 cm²/Vs.
- Hecht (`charge_collection`) and DD use the same hole mobility within 5 %.
- Update `test_device_mobility_at_300k` to the new semantics.
- Re-run the transient and CCE suites and record any baseline shifts.

## Effort

M
