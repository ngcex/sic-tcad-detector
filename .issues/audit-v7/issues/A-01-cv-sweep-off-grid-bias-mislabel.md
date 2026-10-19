---
id: A-01
title: cv_sweep ramps on a fixed 0.5 V lattice and records the requested bias, so off-grid sweeps are mislabelled and sub-0.5 V spacing never solves
severity: high
category: bug
status: fixed  # 2026-10-18 sync: 5964fac: cv_sweep solves off-lattice targets (tests/test_cv_sweep_lattice.py); run_cv on the batched engine is exact (PERF-02 item checked)
location: etna/core/cv_analysis.py:185-226 (cv_sweep ramp loop); callers etna/api/simulation.py:99,113 (run_cv), app/workflows/cv.py:42, scripts/validation/cv_vs_bruzzi2023.py:28
labels: [audit-v7, slice-A, cv, solver, silent-failure]
source: raw-A-core-solvers.md F1
---

## Summary

The inner ramp in `cv_sweep` starts at `V = current_V_cathode + 0.5` and steps by 0.5 V
while `V <= V_cathode_target`. On success it unconditionally sets
`current_V_cathode = V_cathode_target` and records W under `V_target`. So:

- The last bias actually applied is the largest point of the lattice `current + k*0.5`
  that does not exceed the target. It is not the target itself.
- If the next target is less than 0.5 V beyond the previous one, the loop body never
  runs: `set_parameter` is never called, but the bookkeeping still advances. The device
  stays at the last bias that was really applied, and the stale W is recorded under the
  new label.
- `truncated` stays `False` in every case.

`device2d._ramp_to` (`V = min(V + 0.5, target)`) and `drift_diffusion.ramp_bias` both
land exactly on the target. Only `cv_sweep` has this bug.

## Reproduction / Evidence

Commands I ran (devsim 2.10.0, `uv run python`):

1. A mixed grid on a real DD device:
   ```
   di = create_dd_device("cva"); r = cv_sweep(di, [0,-0.25,-0.75,-1.2,-1.5])
   devsim.get_parameter(device="cva", name=GetContactBiasName("cathode"))  -> 0.75
   W (um) = [1.682 1.682 1.891 1.891 1.891], truncated=False
   ```
   The "-1.5 V" point was solved at 0.75 V.
2. Sub-0.5 V spacing on a real DD device:
   ```
   di = create_dd_device("cvb"); r = cv_sweep(di, np.linspace(0,-5,40))
   applied cathode bias after sweep = 0.0, truncated=False, n=40, unique W values = 1 (1.682 um)
   ```
   The solver is never called. The C-V curve is flat and the result does not report it.
3. Default `run_cv` grid, `linspace(0,-50,40)`, 1.282 V step (pure-Python replica of the
   loop): the final applied bias is 49.718 V against the -50 V label. Every point lags its
   label by a constant 0.282 V. The replica is valid for grids with spacing of 0.5 V or
   more.
4. Call sites (grep):
   - `etna/api/simulation.py:99` builds `bias_array = np.linspace(v_start, v_stop, n_points)` with defaults 0 / -50 / 40.
   - `app/workflows/cv.py:42` calls `etna.run_cv(cfg)` with defaults.
   - `scripts/validation/cv_vs_bruzzi2023.py:28` calls `run_cv(cfg)` with defaults.
5. `scripts/validation/cv_vs_bruzzi2023.py:79` fits the 1/C² slope only, to recover N_D.
   A constant voltage offset does not change the slope, so the published N_D recovery is
   unaffected. No code in the repo extracts V_bi from the 1/C² intercept. The raw audit's
   "~10 % V_bi shift" is therefore a potential consequence, not an observed one.

From the raw audit, not re-run: existing tests (`tests/test_cv.py`, `tests/test_api_cv.py`,
the v3 baseline) use only grids whose spacing is a multiple of 0.5 V.

## Impact

- **Default public grid:** C-V points from `run_cv`, the Streamlit page and the validation
  script are labelled 0.282 V deeper than the bias actually applied.
- **`n_points` large enough for spacing below 0.5 V** (for example `run_cv(cfg, n_points=200)`):
  the C-V curve is flat and silent. W(0) is reported at every voltage with
  `truncated=False`.
- **Any downstream intercept analysis** (V_bi) would be biased.

## Proposed fix

```python
# etna/core/cv_analysis.py, cv_sweep inner ramp (mirror device2d._ramp_to)
V = current_V_cathode
converged = True
while abs(V - V_cathode_target) > 1e-10:
    if V_cathode_target > V:
        V = min(V + 0.5, V_cathode_target)
    else:
        V = max(V - 0.5, V_cathode_target)
    V = round(V, 10)
    devsim.set_parameter(device=device, name=bias_name, value=V)
    try:
        devsim.solve(...)            # unchanged
    except devsim.error:
        try:
            devsim.solve(...)        # unchanged relaxed retry
        except devsim.error as e:
            converged = False
            break
if converged:
    current_V_cathode = V            # bias actually applied (== target)
```

Also fix the docstring wording "should be <= 0" (see A-10).

## Test plan

- `cv_sweep(dev, [0,-0.25,-0.75,-1.2])`: afterwards the cathode bias parameter equals
  1.2, and W is strictly increasing.
- `cv_sweep(dev, np.linspace(0,-5,40))`: W is strictly increasing and has 40 unique
  values.
- `cv_sweep(dev, np.linspace(0,-50,40))`: the final cathode bias parameter equals 50.0
  within 1e-9.
- Re-run `tests/test_cv.py`, `tests/test_api_cv.py` and the v3 baseline regression. These
  use on-lattice grids, so their results should not change.

## Effort

S
