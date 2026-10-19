---
id: B-06
title: Dark-current calibration claims (18 pA / 8 pA / N_t=1e12) contradict each other and the code (~141 pA at -30 V, 0.05 cm^2); area defaults differ
severity: medium
category: docs
status: fixed  # 2026-10-18 sync: 0af94c1: N_t recalibrated to 2.8e12, 18 pA is the calibration input, not a prediction (PERF-02 section 7)
location: etna/core/dark_current.py:11-20 (module docstring), 115-117 (setup_tat_model docstring), 184-188 (C8 comment), 353/429/562 (area=0.05), 700/925 (area=0.04); etna/core/sic_material.py:68; tests/test_dark_current.py:228-247
labels: [dark-current, calibration, documentation, physics-honesty, raw-F6]
---

## Summary

The calibration is described four different ways:

- Module docstring: "Default N_t = 1e12 cm^-3 s^-1 targets ~18 pA".
- `setup_tat_model` docstring: "Default 1e12 targets ~8 pA at -30V".
- `SiC4H_Parameters.N_t = 2.2e13`, annotated "(calibrated ~18 pA at -30V)".
- C8 comment: "the existing 18 pA calibration at 300 K is preserved".

The code actually gives about **141 pA** at -30 V over 0.05 cm^2, which is 7.8x the stated target. The regression test accepts anything from 0.1x to 10x of 18 pA, so nothing pins the claimed value. The `area` defaults also differ: 0.05 cm^2 in `extract_dark_current_components`, `dark_current_sweep` and `sensitivity_sweep`, but 0.04 cm^2 in `dark_current_vs_fluence` and `dark_current_post_anneal`. Ampere figures from different functions are therefore not comparable.

## Reproduction / Evidence

- ran (scratch script, `uv run python`): `dark_current_sweep(create_dark_current_device(T=300, doping_profile="graded", N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1e-4), [-30.0], area=0.05)` gives `J_total = 2.825e-09 A/cm^2` and **`I_total = 1.413e-10 A` (141 pA)**.
- ran: `sed -n 1,25p; 108,120p; 180,192p etna/core/dark_current.py` and `grep -n N_t etna/core/sic_material.py` show the four statements quoted above.
- ran: `grep -n "area=0.0" etna/core/dark_current.py` shows lines 353, 429 and 562 with 0.05, and lines 700 and 925 with 0.04.
- call chain: `test_dark_current_order_of_magnitude` asserts `1.8e-12 <= I <= 180e-12`. It uses the default `create_dark_current_device()` kwargs, which may differ slightly from the graded kwargs used above. The magnitude still sits well above 18 pA. The 18 pA figure would imply an area of about 0.0064 cm^2, and no function uses that.

## Impact

Anyone reading the docs believes the model is calibrated to 18 pA when it outputs about 140 pA. Noise-floor and SNR conclusions that cite "~18 pA" are off by about 8x. This is a physics-honesty issue of the same kind as the earlier v6 items.

## Proposed fix

- Choose one of two options:
  - (a) Keep N_t = 2.2e13 and document it honestly: "N_t = 2.2e13 gives I(-30 V, 0.05 cm^2) ~ 1.4e-10 A (J ~ 2.8e-9 A/cm^2); the 18 pA experimental target is not reproduced".
  - (b) Re-calibrate N_t to about 2.8e12 to hit 18 pA at the experimental area, once that area is confirmed.
- Remove "1e12" from both docstrings.
- Add `DEFAULT_AREA_CM2 = 0.05` to the module and use it in every function signature.

```diff
-    voltage-dependent scaling. Default N_t = 1e12 cm^-3 s^-1 targets ~18 pA.
+    voltage-dependent scaling. Default N_t = params.N_t (2.2e13 cm^-3 s^-1),
+    giving J(-30 V) ~ 2.8e-9 A/cm^2 (~140 pA over 0.05 cm^2).
```

## Test plan

Tighten `test_dark_current_order_of_magnitude` to the documented value within 30%. Add a test asserting that every public function in dark_current.py has the same `area` default (signature introspection).

## Effort

S for option (a), the docs-only fix. M for option (b), recalibration.
