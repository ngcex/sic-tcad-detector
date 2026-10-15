---
id: A-07
title: Notebook 19 calls poisson.ramp_voltage with ramp_bias's kwargs (V_target=, contact=) -> guaranteed TypeError, masked by an earlier convergence failure
severity: low
category: bug
status: open
location: scripts/create_notebook_19.py:513-514; notebooks/19_alternative_structures.ipynb cell 12; etna/core/poisson.py:169-178 (signature)
labels: [audit-v7, slice-A, notebook, latent-crash]
source: raw-A-core-solvers.md F7
related: A-06
---

## Summary

The notebook calls:

```
ramp_voltage(dee_info, V_target=V_BIAS, contact="de_anode", V_step=5.0)
```

The signature is `ramp_voltage(device_info, contact_name="cathode", V_start=0.0, V_end=-60.0, V_step=-0.5, ...)`,
with no `**kwargs`. `V_target` and `contact` cannot bind, so this call is a guaranteed
`TypeError`. The three sibling structures in the same script correctly use
`ramp_bias(..., V_target=, contact=, V_step=)`.

## Reproduction / Evidence

- `grep -n "ramp_voltage\|ramp_bias" scripts/create_notebook_19.py` shows lines 83 and
  513 (imports from etna.core.poisson) and line 514 (the bad call). Lines 232, 293 and
  398 are the correct `ramp_bias` calls.
- I read `etna/core/poisson.py:169-178` and confirmed there is no `V_target`, no
  `contact` parameter and no `**kwargs`.
- JSON inspection of notebooks/19 cell 12 shows an output tail of
  "Delta-E/E structure failed: Convergence failure! Continuing with other structures...".
  The cell dies at the equilibrium solve before it reaches line 514.

## Impact

The Delta-E/E study is currently broken by convergence. Once that is fixed, it will fail
with a TypeError instead. A second problem: even with the right signature, `ramp_voltage`
is a Poisson-only ramp, which has the contact artifact described in A-06.

## Proposed fix

```diff
-    from etna.core.poisson import ramp_voltage
-    ramp_voltage(dee_info, V_target=V_BIAS, contact="de_anode", V_step=5.0)
+    ramp_bias(dee_info, V_target=V_BIAS, contact="de_anode", V_step=5.0)
```

This requires the Delta-E/E device to have DD equations set up (`setup_sic_drift_diffusion`),
as the sibling structures do. Regenerate the notebook afterwards.

## Test plan

- Add a `pytest.mark.slow` smoke test that runs the Delta-E/E builder from the
  notebook-19 generator once and asserts that it reaches the "Bias ramped" print. The
  test may be `xfail` until the convergence failure is fixed.
- Add a quick static check: an AST-scan test asserting every `ramp_voltage(...)` call in
  `scripts/` uses only valid keyword names.

## Effort

S
