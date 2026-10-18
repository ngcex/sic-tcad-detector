# Engine validated boxes and devsim ramp steps (WP-B, 2026-10-18)

## 1D engine T box: 250-400 K -> 150-450 K

run_cv (C, W) and run_cce (CCE) on the batched engine vs the devsim path
(`sim._run_cv_devsim` / `sim._run_cce_devsim`), engine box T bound lifted for the run.
21-point C-V 0..-100 V (legacy stack: 0..-40 V), 7-point CCE -10..V_stop.

| stack                  | T (K)   | max abs(dC/C) | max abs(dW) (cm) | max abs(dCCE) |
| ---------------------- | ------- | ------------- | ---------------- | ------------- |
| default (measured C-V) | 150-500 | <= 2.9e-8     | <= 1.3e-18       | <= 7.8e-16    |
| nplus_calibrated       | 150-500 | <= 3.2e-8     | <= 1.6e-18       | <= 5.6e-16    |
| legacy (no n+)         | 150-500 | <= 3.2e-8     | <= 2.6e-18       | <= 1.3e-13    |
| microdosimeter         | 150-500 | <= 2.6e-8     | <= 4.3e-19       | <= 1.1e-15    |

(abs(dC/C) ~ 3e-8 is the devsim central-difference capacitance, same as inside the old box.)
At 100 K devsim itself fails on the legacy stack (cv_sweep convergence failure at -30 V),
so no reference there. Adopted box: 150-450 K (500 K stays outside: several tests use it
as the "devsim path" example). Pinned by tests/test_engine_T_box.py.
run_cv / run_cce at 200 K: devsim 2.8 / 3.3 s -> engine 0.02 / 0.12 s.

Bench: lab/engine_box/bench/t_box_parity1d.py 150,200,225,425,450,500

## 2D engine epi box: 5-20 um -> 5-30 um (ENGINE2D_BOX and ENGINE2D_DARK_BOX)

lab/engine_box/bench/epi30_parity2d.py: dark_state_2d vs create_2d_dd_device (Potential on
every node, the run_field path) and lateral_cce_scan_2d vs cce_lateral_scan (3 points),
device keywords from device_kwargs(DeviceConfig(...)) as the facades pass them.

| case                                              | max abs(dpsi) (V) | max abs(dCCE) | engine | devsim  |
| ------------------------------------------------- | ----------------- | ------------- | ------ | ------- |
| default, epi 22, hw 50, 50 V                      | 7.8e-14           | 3.3e-14       | 0.8 s  | 68.1 s  |
| default, epi 25, hw 25, 20 V                      | 1.4e-14           | 3.1e-14       | 0.4 s  | 26.5 s  |
| default, epi 30, hw 25, 60 V                      | 1.4e-13           | 1.1e-14       | 0.4 s  | 48.4 s  |
| default, epi 30, hw 25, 100 V                     | 1.4e-13           | 1.2e-14       | 0.5 s  | 53.1 s  |
| n+ without graded transition, epi 30, hw 25, 60 V | 1.1e-13           | 1.4e-14       | 0.5 s  | 41.3 s  |
| N_D_bulk 2e13, epi 30, hw 150, 80 V               | 3.1e-13           | 2.5e-14       | 7.9 s  | 461.0 s |
| legacy (no n+), epi 30, hw 25, 60 V               | 1.4e-14           | 3.8e-15       | 0.1 s  | 13.6 s  |

run_field(DeviceConfig(half_width_um=50, epi_thickness_um=22), -50): devsim 68-90 s ->
engine2d 0.8 s. validated_domain_reason now reports the 2D engine envelope for 2D configs
(was the fixed string "2D device"). Pinned by tests/test_engine2d_epi30.py.

## Dark-current ramp step (dark_current_sweep V_step 0.5 V -> 2 V)

lab/engine_box/bench/dark_vstep.py, default config, run_dark_current's grid (0..-100 V,
20 points). With the old single refresh + re-solve per bias, the result depends on the
ramp path (one-step lag of the frozen TAT field), max |dI/I| vs the 0.5 V reference:

| T (K) | 1 V     | 2 V     | 3 / 5 V | 6 V     |
| ----- | ------- | ------- | ------- | ------- |
| 250   | 0       | 0       | 2e-16   | 2e-16   |
| 300   | 9.1e-15 | 2.9e-14 | 4.2e-14 | 5.1e-14 |
| 350   | 1.6e-10 | 4.9e-10 | 7.3e-10 | 8.5e-10 |
| 400   | 2.6e-7  | 8.1e-7  | 1.2e-6  | 1.4e-6  |

and on a 51-point grid at 400 K the 2 V ramp is 4.3e-6 off the 0.5 V one. Iterating the
refresh to the fixed point (E_field_node = E(state), `_resolve_self_consistent`) removes
the path dependence: the fixed points of the 0.5 V and 2 V ramps agree to 4e-16, and every
step size above agrees with the 0.5 V fixed point to <= 3.1e-15. The old 0.5 V single-refresh
reference is itself off the fixed point by 3.2e-7 (400 K, 20 points; 1.1e-6 on 51 points),
<= 1.9e-10 at <= 350 K, 9e-15 at 300 K (the 18 pA N_t calibration is unaffected).

Adopted: fixed-point iteration + DARK_RAMP_STEP_V = 2 V, the fastest step measured
(3 / 5 V take more Newton iterations per step and are slower than 2 V).
run_dark_current(DeviceConfig()): 3.45 s -> 2.10 s (same machine, back to back).
The three other copies of the ramp (sensitivity_sweep, dark_current_vs_fluence,
dark_current_post_anneal) now share dark_current_at_bias.

Refactored callers vs the pre-change module (git show 2f86879:etna/core/dark_current.py),
lab/engine_box/bench/dark_refactor_parity.py, area 0.04:

- dark_current_vs_fluence, fluences 0 / 1e12 / 1e13 / 2e13 (the last fully compensated,
  invalid in both), -30 and -80 V: |dI/I| <= 4.7e-15, same validity flags.
- dark_current_post_anneal (1e13, 600 K, 1800 s): 4.4e-15; (2e13, ...) invalid in both.
- sensitivity_sweep N_t / epi_thickness_cm at -30 V: 300 K <= 4.1e-15, 400 K <= 1.1e-7.
  sensitivity_sweep could not be called twice in one process (fixed device names; devsim
  keeps the mesh after delete_device): names now carry a uuid.
