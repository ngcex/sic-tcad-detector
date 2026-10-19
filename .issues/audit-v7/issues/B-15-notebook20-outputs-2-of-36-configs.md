---
id: B-15
title: Notebook 20 (feasibility report) committed outputs rest on 2 of 36 sweep configurations (stale May-2026 run with a solver_callback cascade)
severity: medium
category: docs
status: closed 2026-10-18 (5ad2ec8 re-executed notebook 20 after the B-02/B-03 fixes and retracted the claims; 8498e1e strips all notebook outputs. Evidence: notebook 20 has 0 outputs and 0 'Valid results: 2' / 'solver_callback' hits; its cells print the valid count of the sweep explicitly)
location: notebooks/20_feasibility_report.ipynb (sweep, ranking and scoring cells); scripts/create_notebook_20.py:190-202; etna/core/optimization.py:180-208; etna/core/devsim_reset.py:60-136
labels:
  [
    notebooks,
    microdosimetry,
    stale-outputs,
    feasibility,
    physics-honesty,
    raw-F15,
  ]
---

## Summary

The committed executed output of notebook 20 prints "Sweep complete: 36 configurations / Valid results: 2 / Failed: 34". The failure log contains 34 "failed (hw=" lines. According to the raw audit, failure 3 is a ramp convergence failure at 80 V (expected punch-through). The rest are `DEVSIM FATAL: Solver "custom" specified, but "solver_callback" not set`, a cascade caused by a post-failure reset that lost the solver configuration. The notebook file contains 148 `solver_callback` hits. The ranking, structure scoring and conclusions were therefore computed from 2 configurations.

**Correction to the raw finding (dates and history):**

- The raw audit says the notebook and `devsim_reset.py` were committed together on 2026-07-26. In fact the notebook was last touched in 986869e (2026-07-26) and `devsim_reset.py` in b337371 (2026-07-15, rename).
- Pickaxe (`git log -S'Valid results: 2'`) shows the 2/36 output **first appeared in f434c64 (2026-05-17)**.
- `git log --all -S'saved_callback'` shows that the solver-preserving reset arrived in that same commit f434c64.
- The executed outputs therefore predate the reset fix taking effect, or at least were not regenerated after it. They also predate the audit Mj-3 full-depletion gate. That explains why the logged failures are the callback cascade and not B-03's `TypeError`.
- With current code, re-execution would fail 36/36 because of B-03. The status here is `open` rather than the raw audit's "plausible", because every fact used is confirmed from the repo. Only the effect of re-execution is inferred.

## Reproduction / Evidence

- ran: `grep -c "failed (hw=" notebooks/20_feasibility_report.ipynb` gives **34**.
- ran: `grep -o "Valid results: [0-9]*\|Failed: [0-9]*\|Sweep complete: [0-9]* configurations"` gives "Sweep complete: 36 configurations", "Valid results: 2" and "Failed: 34".
- ran: `grep -o 'solver_callback...' | uniq -c` gives 148 hits.
- ran: `git log --format='%h %cd %s' -- notebooks/20_feasibility_report.ipynb` shows 986869e (Jul 26), 2cb7ae6 (Jul 15), ec13554, a4f5341, f434c64 (May 17), 2df010b.
- ran: `git log -S'Valid results: 2' -- notebooks/20_feasibility_report.ipynb` gives f434c64. `git log --all -S'saved_callback'` gives f434c64 and b5046a5.
- B-03 shows that current code makes every configuration fail.

## Impact

The headline feasibility conclusions for the microdosimetry design (Giada's second fab deliverable) are based on 2 valid points out of a 36-point grid. Anyone reading the committed notebook cannot tell this without scrolling through the logs.

## Proposed fix

1. Fix B-03 first, then B-02, since the noise floor uses `get_dark_current_2d`.
2. Re-execute notebook 20 end to end.
3. Add a guard cell after the sweep:

```python
n_valid = int(sweep_df["is_valid"].sum())
assert n_valid >= 0.5 * len(sweep_df), f"only {n_valid}/{len(sweep_df)} valid configs"
```

4. Add notebook 20 to `scripts/regression_sweep_v3_notebooks.py`, or add a slow test that runs a reduced grid of 2-4 configs.

## Test plan

After the fixes, re-executed notebook 20 shows at least N valid configs. The remaining failures are only expected punch-through cases, with no "solver_callback not set" and no TypeError in the output.

## Effort

M
