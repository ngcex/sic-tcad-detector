---
id: D-13
title: Five tests have no assertions (incl. print-only test_full_pipeline_diagnostic and a slow transient smoke test)
severity: low
category: improvement
status: open
location: tests/test_analytical.py:328; tests/test_device2d.py:475; tests/test_devsim_reset.py:39; tests/test_flash_recombination.py:323; tests/test_poisson.py:55; tests/test_poisson.py:164 (W_60 unasserted)
labels: [tests, weak-assertions]
source: raw-D F13
---

## Summary

An AST scan for `test_*` functions with no `assert` statement, no `assert*` call and no `pytest.raises`/`warns`
finds exactly five:

- `test_analytical.py:328 test_full_pipeline_diagnostic`: prints a table and is always green. It is a diagnostic, not a test.
- `test_device2d.py:475 test_importable`: smoke test, acceptable.
- `test_devsim_reset.py:39 test_idempotent_on_empty_state`: checks only that nothing raises. It also leaves global state broken (see D-01).
- `test_flash_recombination.py:323 test_transient_solve_after_generation_and_auger`: a slow devsim transient run whose only check is that it did not raise.
- `test_poisson.py:55 test_equilibrium_converges`: smoke test, acceptable.

Separately, `W_60` at `test_poisson.py:164` is computed and never used (also listed in D-09).

## Reproduction/Evidence

```
$ uv run python -I scratchpad/noassert.py   (AST scan, counts Assert nodes, assert*/raises/warns calls, with pytest.raises)
tests/test_analytical.py 328 test_full_pipeline_diagnostic
tests/test_device2d.py 475 test_importable
tests/test_devsim_reset.py 39 test_idempotent_on_empty_state
tests/test_flash_recombination.py 323 test_transient_solve_after_generation_and_auger
tests/test_poisson.py 55 test_equilibrium_converges
```

## Impact

Coverage looks better than it is. The slow transient test costs CI time without checking any physics.

## Proposed fix

- `test_full_pipeline_diagnostic`: convert it to assertions on W at the listed biases (against analytical
  values with tolerance), or move it to `scripts/`.
- Transient smoke test: add `assert np.all(np.isfinite(result["I"]))` and a charge-conservation check (∫I dt ≈ q·N_gen within a tolerance).
- `test_poisson.py:164`: add the missing `W_60` assertion.
- `test_idempotent_on_empty_state`: add `assert devsim.get_device_list() == []` and restore the solver (D-01).

## Test plan

Re-run the AST scan: only the deliberate smoke tests (`test_importable`, `test_equilibrium_converges`)
remain, or they are explicitly allowlisted.

## Effort

S
