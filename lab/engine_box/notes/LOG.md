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

| case | max abs(dpsi) (V) | max abs(dCCE) | engine | devsim |
|---|---|---|---|---|
| default, epi 22, hw 50, 50 V | 7.8e-14 | 3.3e-14 | 0.8 s | 68.1 s |
| default, epi 25, hw 25, 20 V | 1.4e-14 | 3.1e-14 | 0.4 s | 26.5 s |
| default, epi 30, hw 25, 60 V | 1.4e-13 | 1.1e-14 | 0.4 s | 48.4 s |
| default, epi 30, hw 25, 100 V | 1.4e-13 | 1.2e-14 | 0.5 s | 53.1 s |
| n+ without graded transition, epi 30, hw 25, 60 V | 1.1e-13 | 1.4e-14 | 0.5 s | 41.3 s |
| N_D_bulk 2e13, epi 30, hw 150, 80 V | 3.1e-13 | 2.5e-14 | 7.9 s | 461.0 s |
| legacy (no n+), epi 30, hw 25, 60 V | 1.4e-14 | 3.8e-15 | 0.1 s | 13.6 s |

run_field(DeviceConfig(half_width_um=50, epi_thickness_um=22), -50): devsim 68-90 s ->
engine2d 0.8 s. validated_domain_reason now reports the 2D engine envelope for 2D configs
(was the fixed string "2D device"). Pinned by tests/test_engine2d_epi30.py.
