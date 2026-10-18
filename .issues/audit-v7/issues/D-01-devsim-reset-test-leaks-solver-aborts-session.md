---
id: D-01
title: test_devsim_reset.py leaks a non-existent direct_solver into global devsim state and SIGABRTs the pytest session
severity: high
category: bug
status: fixed  # 2026-10-18 sync: 306b657 (file status already fixed)
location: tests/test_devsim_reset.py:39-46 (test_idempotent_on_empty_state), tests/test_devsim_reset.py:75-82 (test_preserves_direct_solver); etna/core/devsim_reset.py:60-139
labels: [tests, devsim, global-state, crash]
source: raw-D F1
---

## Summary

`TestResetClearsState.test_preserves_direct_solver` sets `direct_solver="mkl_pardiso"`, calls
`reset_devsim_fully(preserve_solver=True)` (which restores it by design) and has no teardown. MKL Pardiso
is not available in this devsim build, so the next `devsim.solve` in the same interpreter hits a C++
assert and the process aborts (exit 134). The pytest session dies with no summary.
`test_idempotent_on_empty_state` also ends in `reset_devsim_fully(preserve_solver=False)`, which leaves
`direct_solver == "unknown"`. Any later solve would fail, but the next file test masks this.
Nothing in `etna/` ever puts the solver back.

## Reproduction/Evidence

```
$ uv run pytest tests/test_devsim_reset.py tests/test_drift_diffusion.py::test_dd_equilibrium_convergence -q -p no:cacheprovider
......Fatal Python error: Aborted
  File "etna/core/poisson.py", line 140 in solve_equilibrium
  File "etna/core/drift_diffusion.py", line 237 in create_dd_device
  File "tests/test_drift_diffusion.py", line 83 in test_dd_equilibrium_convergence
EXIT=134
$ uv run pytest tests/test_devsim_reset.py -q      # alone
6 passed in 0.99s
$ grep -rn "direct_solver\|solver_callback" etna/ | grep -v devsim_reset.py
(no output)
$ ls tests/conftest.py
No such file or directory
```

Code check: `devsim_reset.py` step 5 restores whatever snapshot it took, so a bogus value set by the test
carries over. I did not re-run the raw report's separate snippet for the "unknown" path. That claim
follows from the code: with `preserve_solver=False`, `reset_devsim()` runs and nothing restores the solver.

## Impact

`pytest tests/`, or any run that lists this file before another devsim module, crashes the whole run
partway through. CI gives no report, and later modules silently never run.

## Proposed fix

```python
# tests/test_devsim_reset.py
@pytest.fixture(autouse=True)
def _restore_solver():
    devsim = pytest.importorskip("devsim")
    saved = {k: devsim.get_parameter(name=k) for k in ("direct_solver", "solver_callback")}
    yield
    for k, v in saved.items():
        devsim.set_parameter(name=k, value=v)

-        devsim.set_parameter(name="direct_solver", value="mkl_pardiso")
+        original = devsim.get_parameter(name="direct_solver")
+        devsim.set_parameter(name="direct_solver", value="custom")  # exists everywhere
         reset_devsim_fully(preserve_solver=True)
-        assert devsim.get_parameter(name="direct_solver") == "mkl_pardiso"
+        assert devsim.get_parameter(name="direct_solver") == "custom"
```

Optional hardening in `etna/core/devsim_reset.py`: when `preserve_solver=False`, put back devsim's
import-time default (`custom` plus the UMFPACK callback) instead of leaving `"unknown"`.

## Test plan

- Add `test_solve_after_reset_suite`: after the reset tests, build a 1D DD device and assert equilibrium converges.
- Re-run the repro command above. Expect exit 0 and all tests passing.
- `uv run pytest tests/ -q` reaches the summary line.

## Effort

S
