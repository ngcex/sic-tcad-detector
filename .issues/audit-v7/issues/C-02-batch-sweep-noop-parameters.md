---
id: C-02
source: raw-C F2
title: Batch Sweep offers sweep parameters the selected facade ignores, producing N identical overlaid curves
severity: medium
category: bug
status: fixed
location:
  - app/workflows/batch_sweep.py:49-69 (SWEEPABLE_FIELDS x SIM_FACADES)
labels: [app, ui-contract, depends-on-C-01]
---

## Summary

`SWEEPABLE_FIELDS` = epi_thickness_um, substrate_thickness_um, N_A, N_D_junction, N_D_bulk,
L_transition_um, T, area_cm2. `SIM_FACADES` = run_cce, run_cv, run_temperature_sweep. Several
pairs are no-ops, and the page then plots overlapping lines and exports a CSV whose leading
parameter column implies a dependence that does not exist:

- run_cce x {T, N_A, substrate_thickness_um}: ignored because of C-01. `area_cm2` is also a no-op, but it is physically correct for a dimensionless CCE.
- run_temperature_sweep x T: T is the swept axis, so this is by design. x {N_A, substrate_thickness_um}: ignored because of C-01. x area_cm2: no-op.
- run_cv honours all fields through build_device, so it is not affected.

## Reproduction / Evidence

- P2 (own run of repro_api.py): `run_cce(DeviceConfig(T=250))` and `T=350` give bit-identical y. See C-01 for the output.
- Read `app/workflows/batch_sweep.py:49-69`. There is no per-facade allow-list. The module comment says "1D-facade-safe", meaning safe from crashing, not meaningful.
- `_SWEEP_AXIS_TITLES` in app/components/results.py is labels only.

## Impact

The user sees "T=250, T=300, T=350" as three perfectly overlapping CCE curves and may conclude
that CCE is temperature-independent. That is a misleading scientific conclusion with no warning.

## Proposed fix

Fixing C-01 removes most of the no-op pairs. Independently, encode which parameters each facade honours:

```diff
-SIM_FACADES = {"CCE vs bias (run_cce)": "run_cce", ...}
+SIM_FACADES = {
+    "CCE vs bias (run_cce)": ("run_cce", {"epi_thickness_um", "N_D_junction", "N_D_bulk",
+                                          "L_transition_um"}),   # widen after C-01
+    "C-V (run_cv)": ("run_cv", set(SWEEPABLE_FIELDS)),
+    "CCE vs temperature (run_temperature_sweep)": ("run_temperature_sweep",
+        {"epi_thickness_um", "N_D_junction", "N_D_bulk", "L_transition_um"}),
+}
 ...
-param = st.selectbox("Sweep parameter", SWEEPABLE_FIELDS, key="sweep_param")
-sim_label = st.selectbox("Simulation type", list(SIM_FACADES), key="sweep_sim")
+sim_label = st.selectbox("Simulation type", list(SIM_FACADES), key="sweep_sim")
+attr, honoured = SIM_FACADES[sim_label]
+param = st.selectbox("Sweep parameter", [f for f in SWEEPABLE_FIELDS if f in honoured], key="sweep_param")
```

Alternatively, keep the full list and `st.warning(...)` when `param not in honoured`.

## Test plan

AppTest: select facade=run_cce. Assert "T" is absent from the param options, or that a warning
element is present when it is selected. Update this when C-01 lands and widens the honoured set.

## Effort

S

## Status 2026-10-18 (final-c): fixed

- C-01 (930ffdf + 8eff930) made run_cce / run_temperature_sweep honour T, N_A and the substrate; T x run_temperature_sweep is disabled with a warning since cf64bef.
- 4d5fb13: the last no-op pair, area_cm2 x {run_cce, run_temperature_sweep} (CCE is dimensionless), now shows a warning and disables Run; area x run_cv still runs.
- Evidence: `uv run pytest -q tests/test_app_batch_sweep_page.py` -> 9 passed (incl. test_area_sweep_with_cce_facade_is_disabled_with_warning); `uv run pytest -q tests/test_api_config_forwarding.py` -> 19 passed.
