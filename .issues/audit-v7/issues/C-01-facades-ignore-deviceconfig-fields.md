---
id: C-01
source: raw-C F1
title: Most simulation facades silently ignore DeviceConfig fields (T, N_A, substrate, doping_profile/N_D); results and CSV provenance are mislabeled
severity: high
category: bug
status: fixed  # 2026-10-18 sync: 8fd9602 + a0bd37c (file status already fixed)
location:
  - etna/api/simulation.py:422-430 (run_cce)
  - etna/api/simulation.py:501-509 (run_radiation_damage)
  - etna/api/simulation.py:568-576 (run_temperature_sweep)
  - etna/api/simulation.py:642-650 (run_flash_recombination)
  - etna/api/simulation.py:726-735 (run_dark_current)
  - etna/api/simulation.py:816-826 (run_transient)
  - app/components/results.py:307-316 (CSV "# device:" header serializes the full config)
labels: [api, correctness, provenance, physics-honesty]
---

## Summary

`DeviceConfig` is documented as the canonical device contract, and `build_device()` maps
every field. Only `run_cv` and `run_field` go through `build_device`. The other facades
hand-pick kwargs and drop the rest:

| facade                  | forwards                               | silently ignores                                    |
| ----------------------- | -------------------------------------- | --------------------------------------------------- |
| run_cce                 | epi, N_D_junction/bulk, L_transition   | T, N_A, substrate_thickness_um, doping_profile/N_D  |
| run_radiation_damage    | epi, graded triplet                    | T, N_A, substrate, doping_profile/N_D               |
| run_temperature_sweep   | epi, graded triplet                    | N_A, substrate, doping_profile/N_D (T is by design) |
| run_flash_recombination | epi, graded triplet                    | T, N_A, substrate, doping_profile/N_D               |
| run_dark_current        | T, epi, graded triplet (area to sweep) | N_A, substrate, doping_profile/N_D                  |
| run_transient           | epi only                               | everything else, including doping                   |

Correction to the raw finding: `area_cm2` is listed as "ignored" by run_cce in the raw file.
CCE is dimensionless and area-independent in 1D, so it is not a real omission and is removed here.

The app sidebar lets the user choose `doping_profile="uniform"` + `N_D`, `T`, `N_A` and
substrate. CCE / Radiation Damage / Dark Current / Batch Sweep then simulate the default
graded device at default N_A/substrate, while the plot and the CSV header
(`# device: ... doping_profile=uniform, N_D=1e+15 ...`, written from `asdict(result.config)`)
claim otherwise.

## Reproduction / Evidence

Ran the auditor's probe script myself (output saved to
`/private/tmp/claude-501/-Users-ngcex-projects-physics-petringa/df2eeea8-54e5-4d29-83f3-c7a1539b1919/scratchpad/repro_C_out.txt`):

```
uv run python <scratchpad>/repro_api.py
== P1: run_cce ignores doping_profile/N_D (uniform vs graded) ==
graded  y: [0.963266   0.99611564]
uniform y: [0.963266   0.99611564]
identical: True
== P2: run_cce ignores config.T (250K vs 350K) ==
T=250 y: [0.963266 0.99611564]  T=350 y: [0.963266 0.99611564]  identical: True
== P2b: run_cv honours T (sanity, should differ) ==
cv T=250 y: [4.935e-13 1.469e-13]  T=350 y: [5.198e-13 1.517e-13]  identical: False
```

Code read (own verification):

- `etna/api/simulation.py:422-430`: `cce_vs_bias(..., device_kwargs={"N_D_junction", "N_D_bulk", "L_transition"})`. No T, N_A, substrate or profile.
- `:726-735`: `create_dark_current_device(T=..., N_t, S_n, S_p, epi, graded triplet)`. No N_A, substrate, profile or N_D.
- `:822-826`: `transient_cce_vs_dose_rate(V_bias, dose_rates, epi_thickness_cm)` only.
- Core signatures: `cce_vs_fluence` (charge_collection.py:626), `cce_vs_dose_rate` (flash_recombination.py:263) and `transient_cce_vs_dose_rate` (transient.py:545) have fixed kwargs with no T/N_A/substrate/profile. `sweep_cce_vs_temperature(**device_kwargs)` (temperature_sweep.py:215) and `create_dark_current_device(**kwargs)` (dark_current.py:532) already pass kwargs through, so those two facades can forward today.
- `etna/api/device.py:79-90`: `build_device` forwards every field, which is why run_cv/run_field are correct.
- The default `DeviceConfig` values (substrate 1 um, N_A 1e19, T 300) coincide with the core defaults. That is why the existing default-config tests pass.

## Impact

The tool's stated goal is foundry-ready design numbers. Any non-default sidebar choice of
profile, N_D, N_A, substrate or (for CCE/damage/flash) T produces a result that is labeled
with a configuration that was never simulated. That includes the exported CSV provenance.
This silently invalidates parameter studies.

## Proposed fix

```diff
# etna/api/device.py
+def device_kwargs(config: DeviceConfig) -> dict:
+    """Single source of truth for DeviceConfig -> create_dd_device kwargs (1D)."""
+    return dict(
+        epi_thickness_cm=config.epi_thickness_um * 1e-4,
+        substrate_thickness_cm=config.substrate_thickness_um * 1e-4,
+        N_A=config.N_A, N_D=config.N_D, T=config.T,
+        doping_profile=config.doping_profile,
+        N_D_junction=config.N_D_junction, N_D_bulk=config.N_D_bulk,
+        L_transition=config.L_transition_um * 1e-4,
+    )
# build_device 1D branch: create_dd_device(device_name=name, **device_kwargs(config))

# etna/api/simulation.py  run_cce
-    result = cce_vs_bias(V_range=bias_array, epi_thickness_cm=..., device_kwargs={...triplet...})
+    kw = device_kwargs(config); epi = kw.pop("epi_thickness_cm")
+    result = cce_vs_bias(V_range=bias_array, epi_thickness_cm=epi, device_kwargs=kw)

# run_temperature_sweep: pass **{k: v for k, v in device_kwargs(config).items() if k != "T"}
# run_dark_current: create_dark_current_device(N_t=..., S_n=..., S_p=..., **device_kwargs(config))

# core: add `device_kwargs=None` passthrough to cce_vs_fluence, cce_vs_dose_rate,
# transient_cce_vs_dose_rate (mirroring cce_vs_bias), then forward from the facades.
```

Interim option if the core change is deferred: in the facades that cannot forward yet, raise
`NotImplementedError` (or log a warning and set `metadata["ignored_config_fields"]`) when a
non-forwardable field differs from its core default, and have the app surface that warning.

## Test plan

- Fast: monkeypatch `create_dd_device` (or each core entry point) to capture kwargs. For each facade, assert that every DeviceConfig field arrives with the right unit conversion.
- Slow (`@pytest.mark.slow`): `run_cce(DeviceConfig(T=250)).y != run_cce(DeviceConfig(T=350)).y`, and `DeviceConfig(doping_profile="uniform", N_D=1e15)` differs from default for run_cce, run_dark_current and run_radiation_damage.
- Regression: `run_cce(DeviceConfig())` is unchanged vs current output (defaults coincide).

## Effort

M (helper plus 3 easy facades is S; the core passthrough for the 3 fixed-signature functions and their tests pushes it to M)
