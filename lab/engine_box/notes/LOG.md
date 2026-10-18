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
