---
id: B-16
title: simulate_single_particle masks a first-step transient failure with IndexError (times[-1] on an empty list)
severity: low
category: bug
status: closed 2026-10-18 (verified by WP-A: already resolved by c1d719c-era refactor; times starts with the injection step, a failed injection raises the devsim error, build_cce_let_table records NaN)
location: etna/core/single_particle.py:226-287 (transient loop and summary log); etna/core/single_particle.py:428-429 (build_cce_let_table generic except)
labels: [single-particle, transient, error-handling, raw-F16]
---

## Summary

`times = []` is initialised before the loop (l.226). If the first post-injection BDF1 step fails at both tolerance levels, the inner `except devsim.error` logs a warning and `break`s before anything is appended. Execution continues as follows:

1. `np.trapezoid` on the empty arrays returns `Q_collected = 0.0`.
2. `logger.info(f"... t_final={times[-1]:.3e} ...")` raises `IndexError`.
3. In `build_cce_let_table`, the generic `except Exception as exc` logs this as "Simulation failed for LET=...: index -1 is out of bounds", which hides the real devsim error.

If the log line were removed or guarded, `Q_collected = 0.0` would be recorded as CCE = 0, a plausible-looking but false value.

## Reproduction / Evidence

- call chain (read): single_particle.py:226-228 initialises `times = []`. Lines 247-261 have the nested `except devsim.error:` and then `break`, before `times.append(t)` at ~l.264. Line ~281 computes `Q_collected = float(np.trapezoid(...))`, which is 0.0 for empty input. Line ~285 uses `times[-1]` and raises IndexError. Lines 428-429 hold the `except Exception as exc: logger.warning(f"Simulation failed for LET=...: {exc}")` handler.
- The branch is reachable. The code already has relaxed-tolerance retries because 2D transients fail at strict tolerances (raw audit, lines 170-207).
- Not executed. Triggering it requires a forced devsim failure; see the test plan.

## Impact

Low. The diagnosis is misleading when a transient fails at the first step. The latent CCE = 0 recording would corrupt the CCE(LET) table if the log were changed.

## Proposed fix

```diff
     times = np.array(times)
     currents = np.array(currents)
+    if times.size == 0:
+        raise RuntimeError(
+            "simulate_single_particle: transient failed at the first post-injection "
+            f"step (t={dt_inject:.3e} s); see preceding devsim warning"
+        )
```

Alternatively, return `Q_collected=np.nan` and have `build_cce_let_table` record NaN rather than 0.

## Test plan

Monkeypatch `devsim.solve` to raise `devsim.error` for `type="transient_bdf1"` after the injection step. Assert that `simulate_single_particle` raises `RuntimeError` with "first post-injection step" in the message, or returns NaN. Assert that `build_cce_let_table` records NaN, not 0.

## Effort

S
