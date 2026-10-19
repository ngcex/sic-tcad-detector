---
id: D-02
title: DESIGN-1 quotes model leakage "~18 pA at -30 V" with no area; the code gives 141 pA (0.05 cm2) or 0.28 pA (API default area)
severity: high
category: docs
status: fixed  # 2026-10-18 sync: 0af94c1: see B-06 (PERF-02 section 7)
location: deliverables/DESIGN-1-dosimetry-pn.md:103-106; etna/core/dark_current.py:13,20,68,117,187; etna/core/sic_material.py:68; notebooks/07_dark_current.ipynb (committed output); .planning/MILESTONES.md:35,61; .planning/PHYSICS_AUDIT_v5.md
labels: [docs, foundry-deliverable, dark-current, calibration]
source: raw-D F2
---

## Summary

The foundry-facing DESIGN-1 says "Our model leakage (~18 pA at −30 V, 300 K) is a single-point
calibration". The current model does not reproduce that number at any stated area. The current density
is about 2.8 nA/cm² on both paths, so 18 pA would need about 6.4e-3 cm², and neither code path uses that area.
The docstrings also contradict each other and the actual default:

- dark_current.py:20 says "Default N_t = 1e12 targets ~18 pA".
- dark_current.py:117 says "Default 1e12 targets ~8 pA".
- sic_material.py:68 sets `N_t = 2.2e13` with the comment "calibrated ~18 pA".

Correction to the raw report: README.md:25-26 and 156-157 do **not** quote 18 pA. They only say
"calibrated, single-point" / "single-point fit", which is accurate, so README is not affected.

## Reproduction/Evidence

```
# core defaults (area 0.05 cm^2)
dark_current_sweep(create_dark_current_device(), V_range=[0,-10,-20,-30])
-> V=-30: I_total=1.4125e-10 A, I_SRH=1.70e-14, I_TAT=-1.436e-10, I_SRV=0
# public API (DeviceConfig default area 1e-4 cm^2)
etna.run_dark_current(DeviceConfig(), v_start=0, v_stop=-30, n_points=4)
-> V=-30: I_total=2.824e-13 A
$ grep -o "Simulated:[^\"]*pA" notebooks/07_dark_current.ipynb
Simulated:    1.41e-10 A (141.3 pA
$ grep -o "Target:[^\"]*pA" notebooks/07_dark_current.ipynb
Target:       1.80e-11 A (18.0 pA
```

1.41e-10/0.05 = 2.82e-9 A/cm², and 2.82e-13/1e-4 = 2.82e-9 A/cm²: the density is the same on both paths.

## Impact

The foundry would receive a leakage figure that is 8x off the model at the notebook's area, with no
area stated so it cannot be checked. Notebook 07 already shows the mismatch next to the target. This
undermines the credibility of the deliverable.

## Proposed fix

Option A (docs, preferred given the "estimate only" framing):

```diff
- Our model leakage (~18 pA at −30 V, 300 K) is a single-point calibration
+ Our model leakage density is ~2.8 nA/cm² at −30 V, 300 K (≈141 pA for a 0.05 cm² diode,
+ ≈0.28 pA for 1e-4 cm²). The trap term was originally fitted to an 18 pA measurement; the
+ current defaults no longer reproduce it. Single-point calibration only.
```

Fix the docstrings (`dark_current.py:20,117` → "Default N_t = 2.2e13 (params.N_t)"; drop the pA claims
or state them as a density). Update `.planning/MILESTONES.md` and `PHYSICS_AUDIT_v5.md` the same way.
Option B: recalibrate `N_t` (about 2.2e13 × 18/141 ≈ 2.8e12, if the response is linear) so the −30 V
value at 0.05 cm² is 18 pA, then rerun notebook 07.

## Test plan

Add `test_dark_current_calibration_point`: at −30 V and area 0.05 cm², assert that `I_total` matches the
documented value within 20%. This pins the docs and code together.

## Effort

S (docs) / M (recalibration and notebook rerun)
