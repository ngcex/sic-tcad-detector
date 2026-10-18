---
id: A-12
title: validate_cv / validate_iv silently clamp or snap when the sweep does not reach the experimental voltages
severity: low
category: improvement
status: fixed 2026-10-18 (final-abd 7d7d026: validate_cv NaNs and excludes experimental voltages outside the simulated range, returns extrapolated / n_extrapolated, raises if none inside; validate_iv takes a target only within 0.25 V, else NaN, pass False, missing_targets. Evidence: pytest tests/test_validation.py -> 23 passed, 6 new)
location: etna/core/validation.py:~126-131, ~145-148 (argmin nearest-point picks in validate_iv), ~206 (np.interp clamp in validate_cv)
labels: [audit-v7, slice-A, validation, silent-failure]
source: raw-A-core-solvers.md F12
---

## Summary

- **`validate_cv`:** it calls `np.interp(V_exp, sort(V_sim), W_sim[...])`, and NumPy
  clamps values outside the simulated range. With a truncated C-V (a punch-through wall,
  or a sweep that stops short of the experimental voltages), the last available W is
  reported as `W_sim` at the unreached experimental voltage.
- **`validate_iv`:** it picks `argmin(|V - (-60)|)` and `argmin(|V ∓ 2|)` regardless of
  distance. It then labels whatever point is nearest as `dark_current_60V` and builds the
  rectification ratio from it.

Both functions return pass/fail booleans and R² with no extrapolation flag. `cv_sweep`
already returns `truncated` and `requested_v_stop`, but these functions ignore them.

## Reproduction / Evidence

- I read validation.py around 120-150: `idx_60 = np.argmin(np.abs(V - (-60.0)))` and
  `idx_p2/idx_m2` have no distance check.
- I read validation.py around 200-207: `W_sim_at_exp = np.interp(V_exp, np.sort(V_sim), W_sim[np.argsort(V_sim)])`.
- NumPy documents that `np.interp` returns the edge value outside `[xp[0], xp[-1]]`. For
  example, `np.interp(-30, [-20, 0], [W20, W0]) == W20`.

From the raw audit, not re-run: only notebook 02 calls these functions, and tests always
pass exact voltages.

## Impact

Validation metrics can look acceptable while they are actually computed from unreached
biases. This risk matters given the project's history with circular and overstated
validation claims.

## Proposed fix

```python
# validate_cv
lo, hi = V_sim.min(), V_sim.max()
out = (V_exp < lo - 1e-9) | (V_exp > hi + 1e-9)
result["extrapolated"] = bool(out.any() or cv_data.get("truncated", False))
result["n_extrapolated"] = int(out.sum())
W_sim_at_exp = np.where(out, np.nan, np.interp(...))   # or raise if strict=True

# validate_iv
def _pick(V, target, tol=0.25):
    i = int(np.argmin(np.abs(V - target)))
    return i if abs(V[i] - target) <= tol else None
# None -> value NaN, *_pass False, result["missing_targets"].append(target)
```

## Test plan

- `validate_cv({"voltages":[0,-10], "depletion_widths":[...]})` sets `extrapolated=True`,
  with NaN for the points beyond -10 V, or raises when `strict=True`.
- `validate_iv` with V that lacks -60 V (max reverse -20 V) returns NaN for
  `dark_current_60V` and lists -60 in `missing_targets`.
- The existing tests with exact voltages are unchanged.

## Effort

S
