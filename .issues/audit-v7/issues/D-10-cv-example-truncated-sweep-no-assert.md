---
id: D-10
title: examples/cv_example.py sweeps to -200 V, silently gets 7/20 points, ignores truncation and asserts nothing; spec still shows -200 V
severity: low
category: improvement
status: fixed 2026-10-18 (final-abd 36083d1: cv_example uses run_cv(cfg) (0..-50 V), prints truncated / requested_v_stop, raises on truncation or non-monotonic C; tests/test_examples.py runs it in a subprocess with cwd and ETNA_FIGURE_DIR = tmp_path; spec examples -50 / -40 V; truncation warnings name the first bias not reached. Note: the '7 of 20 points at -200 V' premise is stale, the batched run_cv reaches -200 V. Evidence: pytest tests/test_examples.py -> 2 passed; python examples/cv_example.py -> truncated=False, C 2.0761e-10 -> 3.4910e-11 F)
location: examples/cv_example.py:20-36; docs/superpowers/specs/2026-06-26-simulator-library-ui-design.md:143,147,174,206; etna/core/cv_analysis.py:233 (warning wording)
labels: [examples, docs, punch-through, truncation]
source: raw-D F10
---

## Summary

The "LIB-05 vertical-slice validation script" calls `run_cv(cfg, v_start=0, v_stop=-200, n_points=20)`.
The default device stops converging at about −63 V, so it prints 7 points as though they were the full
sweep. It never reads `metadata["truncated"]`. Its docstring says "Capacitance must decrease
monotonically", but nothing checks that. The design spec still uses `v_stop=-200` in four examples,
while the facade default is now `-50`.

Side issue: the `cv_sweep` warning says "stopped sweep at V=-73.68V", which is the first point that was
_not_ reached. The last converged point is −63.16 V.

## Reproduction/Evidence

```
$ uv run python examples/cv_example.py
cv_sweep: failed at V_cathode=65.66V: Convergence failure!
cv_sweep: stopped sweep at V=-73.68V (requested down to -200.00V)
Bias (V): [0. -10.53 -21.05 -31.58 -42.11 -52.63 -63.16]
C at 0.0V = 5.0608e-13 F, C at -63.2V = 8.6003e-14 F (decreasing with reverse bias)
$ grep -n v_stop docs/superpowers/specs/2026-06-26-simulator-library-ui-design.md  -> 143,147,174,206 all -200
$ sed -n 34,38p etna/api/simulation.py  -> v_stop: float = -50.0
```

## Impact

The example teaches users to request sweeps the solver cannot deliver and hides the truncation. Its
"validation" is cosmetic, and the warning text misreports where the sweep ended.

## Proposed fix

```diff
-    result = run_cv(cfg, v_start=0, v_stop=-200, n_points=20)
+    result = run_cv(cfg)  # facade default: 0 .. -50 V
+    print(f"truncated={result.metadata.get('truncated')} "
+          f"requested_v_stop={result.metadata.get('requested_v_stop')}")
 ...
+    assert np.all(np.diff(result.y) <= 0), "C must decrease monotonically with reverse bias"
```

In `cv_analysis.py:233`: `f"cv_sweep: stopped sweep after V={solved_voltages[-1]:.2f}V (failed at {V_target:.2f}V; requested down to ...)"`.
Update the four spec examples to `v_stop=-50` (and `-40` for CCE, matching that facade's default).

## Test plan

Run `uv run python examples/cv_example.py`: no truncation warning and the assertion passes. Optionally
add `tests/test_examples.py`, marked slow, that runs the example.

## Effort

S
