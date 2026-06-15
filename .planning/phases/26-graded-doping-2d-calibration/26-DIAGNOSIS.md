---
phase: 26
document: DIAGNOSIS
created: 2026-06-15T12:02:00.435038+00:00
hypothesis: H1
---

# Phase 26 — 1D vs 2D Reverse-Bias Convergence Diagnosis

Profile under test (v3.0 defaults): `N_D_junction=2.9e15`, `N_D_bulk=8.5e13`,
`L_transition=1e-4` cm. Each device was ramped from 0 V toward -50 V in 0.5 V
cathode steps using the cv_sweep adaptive-fallback solver pattern.

## Results

| Device          | V_max_reverse_bias (V) | W(-10 V) (cm) | W(-30 V) (cm) | Failure mode |
|-----------------|------------------------|---------------|---------------|--------------|
| 1D (device.py)  | -50.00             | 9.593e-04    | 9.984e-04    | reached -50 V without failure    |
| 2D 100 µm SV    | -50.00         | n/a           | n/a           | reached -50 V without failure|
| 2D 300 µm SV    | -50.00         | n/a           | n/a           | reached -50 V without failure|

Machine-readable evidence:

- `V_1d_max_converged`: -50.00 V
- `V_2d_max_converged_100um`: -50.00 V
- `V_2d_max_converged_300um`: -50.00 V

## Hypothesis: H1

Both the 1D and/or 2D devices fail to reach -30 V with the v3.0 `{2.9e15, 8.5e13, 1e-4}` profile, OR both reach -30 V comfortably. In either reading the *profile* is the dominant lever: if it fails, the profile is marginal and a Nelder-Mead refit (Plan 03) is justified; if both converge deep, the profile is already adequate and calibration is fine-tuning only. Either way the lever is the doping profile, not the mesh/BC.

## Consumed by

- .planning/phases/26-graded-doping-2d-calibration/26-03-PLAN.md (Task 1 reads `hypothesis:` and gates calibration)
