---
id: A-13
title: iv_sweep assumes the device starts at 0 V instead of reading the live contact bias
severity: low
category: improvement
status: fixed 2026-10-18 (final-abd 54c7135: iv_sweep starts from devsim.get_parameter(<contact>_bias), returns truncated, docstring states the anode sign. Evidence: pytest tests/test_iv_sweep_balanced.py::test_iv_sweep_starts_from_the_live_bias -> the only anode_bias set after ramp_bias(-10) is -9.5)
location: etna/core/drift_diffusion.py:~360-361 (current_V = 0.0 in iv_sweep); contrast ~276-281 (ramp_bias reads devsim.get_parameter)
labels: [audit-v7, slice-A, solver, convergence]
source: raw-A-core-solvers.md F13
---

## Summary

`iv_sweep` initialises `current_V = 0.0`. It does not read the contact's present bias,
which is what `ramp_bias` does. If the device was already biased, for example after
`ramp_bias` or a previous `iv_sweep`, the first Newton solve can jump by the full prior
bias. Avoiding that kind of jump is the whole reason the stepping exists.

There is also a secondary inconsistency. The default `contact="anode"` differs from the
cathode convention used by `cv_sweep`, `ramp_bias` and `run_field`.

## Reproduction / Evidence

I read the code:

- drift_diffusion.py:276-281: `ramp_bias` does `current_V = devsim.get_parameter(device=device, name=bias_name)` with a fallback to 0.0.
- In the `iv_sweep` body (around lines 321-363): the signature defaults
  `contact="anode"`, and the body sets `current_V = 0.0` before computing `bias_name`.

From the raw audit, not re-run: callers are notebook 02 and plotting.py. No caller
pre-biases the device, and no API caller exists. The issue is latent.

## Impact

Convergence failures or wrong branches appear only if `iv_sweep` is chained after another
bias operation. The recorded voltages themselves are correct.

## Proposed fix

```diff
-    current_V = 0.0
     bias_name = simple_physics.GetContactBiasName(contact)
+    try:
+        current_V = float(devsim.get_parameter(device=device, name=bias_name))
+    except devsim.error:
+        current_V = 0.0
```

Document the anode sign convention in the docstring: positive anode voltage means forward
bias.

## Test plan

- `ramp_bias(dev, 5.0, contact="anode")`, then
  `iv_sweep(dev, [4.0, 3.0], contact="anode")`.
- Spy on `devsim.set_parameter`. The first set value must be within `V_step` of 5.0, not
  of 0.0.

## Effort

S
