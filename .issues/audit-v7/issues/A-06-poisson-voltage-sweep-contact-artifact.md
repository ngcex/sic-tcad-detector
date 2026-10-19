---
id: A-06
title: poisson.ramp_voltage / voltage_sweep return contact-artifact E-fields, use an inverted sign convention, and report analytic uniform-N_D widths for graded devices
severity: medium
category: bug
status: fixed 2026-10-18 (final-abd 97e6986: Poisson-only ramp_voltage / voltage_sweep removed, no caller outside test_voltage_sweep_short; bias sweeps are ramp_bias / cv_sweep. Evidence: grep -rn 'voltage_sweep|ramp_voltage' etna app tests scripts notebooks examples -> no hit; pytest tests/test_poisson.py -> 13 passed)
location: etna/core/poisson.py:169-249 (ramp_voltage), 282-332 (extract_depletion_width uses device_info["N_D"]), 541-607 (voltage_sweep); tests/test_poisson.py:213-236 (test_voltage_sweep_short); etna/core/plotting.py:82
labels: [audit-v7, slice-A, poisson, misleading-output]
source: raw-A-core-solvers.md F6
related: A-07
---

## Summary

The Poisson-only model ties carriers to the local potential and has no quasi-Fermi
levels. The quasi-neutral bulk therefore cannot follow the contact, and the applied bias
drops across the first mesh cell at the biased contact. `voltage_sweep` still returns
these fields as `E_fields`/`E_max` against bias. It pairs them with
`extract_depletion_width(V)`, which is the analytic one-sided formula using
`device_info["N_D"]`. For graded devices that value is the uniform fallback of 1.07e15.

There are two further problems:

- **Inverted polarity:** the documented convention ("for reverse bias on the cathode, use
  negative V_end") is the forward-bias polarity used everywhere else in the codebase.
- **Silent truncation:** `ramp_voltage` exits with `break` on solver failure and sets no
  truncation flag.

## Reproduction / Evidence

Command I ran: `create_sic_device('pv'); setup_poisson; solve_equilibrium; ramp_voltage(di,'cathode',0,-5,-0.5)`

```
W_num at 0 V        = 1.682 um
V=-5.0: |E|max      = 7.47e5 V/cm at x = 10.986 um   (cathode at 11.0 um, junction at 1.0 um)
W_num at -5 V       = 1.682 um   (unchanged)
create_sic_device('pg', doping_profile='graded'): N_D = 1.07e15, extract_depletion_width(-30 V) = 5.747 um
```

From the raw audit, not re-run:

- At +5 V the peak is 1.75e6 V/cm at the same contact node.
- The calibrated W(-30 V) is 9.73 um.

Grep shows `voltage_sweep` is used only by tests/test_poisson.py and documented in
plotting.py. `ramp_voltage` is used only by notebook 19 (see A-07). Neither is reachable
from the API or app.

## Impact

Anyone plotting `voltage_sweep` output sees a field spike at the back contact and widths
that do not match the device. `test_voltage_sweep_short` passes only because the contact
artifact grows with |V|. The impact is limited to tests, notebooks and plotting.

## Proposed fix

Choose one:

**(a) Deprecate:** have `ramp_voltage` and `voltage_sweep` emit a `DeprecationWarning`
pointing to `drift_diffusion.ramp_bias` plus `extract_electric_field`.

**(b) Reimplement on DD:**

```python
def voltage_sweep(device_info, voltages=None):
    # device_info must come from create_dd_device
    for V in voltages:                              # conventional, V<=0 = reverse
        ramp_bias(device_info, V_target=-V, contact="cathode")
        x, E = extract_electric_field(device_info)
        W = extract_depletion_width_numerical(device_info)
        ...
    return {..., "truncated": truncated}
```

Under either option:

- Fix the docstring sign.
- Return or log `truncated` when `ramp_voltage` breaks out early.

## Test plan

- Replace `test_voltage_sweep_short` with two checks: the |E| peak lies within 0.5 µm of
  the junction at every bias, and `W_numerical` increases strictly with reverse |V|.
- Graded device: the width reported by the sweep is within 5 % of the
  `calibrate_graded_doping` target at -30 V.

## Effort

M
