---
id: A-04
title: Transient time stepper overshoots phase boundaries; the fall edge (and post-pulse) can be skipped with generation held at 100 %
severity: medium
category: bug
status: fixed 2026-10-18 (WP-A 03c7c17: pulse_time_step clips to phase edges, step-end stamps, BDF1 estimator)
location: etna/core/transient.py:125-169 (adaptive_dt), ~289-336 (simulate_pulse loop)
labels: [audit-v7, slice-A, transient, flash, numerics]
source: raw-A-core-solvers.md F4
---

## Summary

On the plateau, `adaptive_dt` returns `dt_max` without clipping it to the next envelope
breakpoint. The loop clips only to `t_end`. When `t_rise + t_duration` is not a multiple
of the accumulated steps, the last plateau step runs past the plateau end. The whole step
uses `env = pulse_envelope(t_start) = 1.0`, so the fall ramp is swallowed, and sometimes
the post-pulse window too.

Secondary issues, both one-step lags:

- The current solved at `t+dt` is appended with timestamp `t`.
- Generation for the implicit BDF1 step is evaluated at the start of the step.

## Reproduction / Evidence

Command I ran: a pure-Python replica of the loop that imports the real `pulse_envelope`
and `adaptive_dt` and returns (applied ∫env·dt, analytic trapezoid, number of steps
inside the fall window).

```
(1e-6, 1e-3,   1e-6) -> (1.001000e-3, 1.001000e-3, 10)   # defaults: exact
(1e-6, 0.55e-3, 1e-6) -> (5.5645e-4,  5.51e-4,    0)     # +1.0 %, fall never resolved
(1e-6, 50e-6,  1e-6) -> (5.645e-5,   5.1e-5,     0)     # +10.7 %
```

I read the loop code (transient.py ~289-336). `env` is evaluated at `t` before the solve,
and `times.append(t)` follows the solve.

From the raw audit, not re-run:

- `compute_transient_cce` normalises by the analytic trapezoid.
- Tests use only commensurate durations.
- `t_duration` is not exposed in the API or app.

## Impact

Any non-default `t_duration` passed to `simulate_pulse`, `simulate_pulse_train` or
`transient_cce_vs_dose_rate` over-generates charge and never resolves the trailing edge.
The CCE result is biased by the ratio of applied to analytic integral (+1 % to +11 % in
the examples above). The shipped defaults happen to be exact.

## Proposed fix

```python
# etna/core/transient.py, in simulate_pulse loop after dt = adaptive_dt(...)
for b in (t_rise, t_rise + t_duration, t_rise + t_duration + t_fall, t_end):
    if t < b - 1e-15:
        dt = min(dt, b - t)
        break
...
t_next = t + dt
env = pulse_envelope(t_next, t_rise, t_duration, t_fall)   # BDF1: evaluate at step end
add_generation_to_dd(self.device_info, G_spatial * env)
devsim.solve(..., tdelta=dt, ...)
times.append(t_next)
t = t_next
```

## Test plan

- Replica test for `(1e-6, 0.55e-3, 1e-6)`: `sum(env*dt)` equals `t_rise/2 + t_duration + t_fall/2`
  within 1e-9 relative. Use the trapezoid rule on the step-end envelope.
- At least 10 steps fall inside `[t_rise+t_duration, t_rise+t_duration+t_fall)`.
- Add a unit test that `adaptive_dt` combined with the clipping never crosses a
  breakpoint.
- Default-pulse transient tests remain within their current tolerances.

## Effort

S
