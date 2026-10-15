---
id: D-12
title: v3 baseline guard is degenerate (300 um CCE == 100 um CCE) and self-baselined from a dirty tree
severity: low
category: improvement
status: open
location: tests/baselines/v3_frozen.json:2-3,12,20; tests/test_v3_baseline_regression.py:68-105; scripts/freeze_v3_baselines.py:176-182
labels: [tests, regression-baseline, test-cost]
source: raw-D F12
---

## Summary

`cce_center_100um` and `cce_center_300um` are both exactly `0.96085509`. In a laterally uniform planar 2D
device with uniform generation, center CCE does not depend on width. The 300 µm test therefore solves a
larger mesh only to assert the same value: it adds cost and no coverage. The metadata has
`git_commit_sha: "WORKING_TREE"`, because `freeze_v3_baselines.py` writes that whenever
`git status --porcelain` is non-empty, and the values were rebaselined (`rebaselined_plan_26_04`). The
file is a forward guard on current behaviour, not a v3.0 reference as the docstring implies.

## Reproduction/Evidence

```
$ head -25 tests/baselines/v3_frozen.json
"cce_center_100um": 0.96085509,
"cce_center_300um": 0.96085509,
"git_commit_sha": "WORKING_TREE",
"rebaselined_plan_26_04": "cce_center_*um updated from 0.96446952 (v3.0 snapshot) to 0.96085509 ..."
$ sed -n 174-182 scripts/freeze_v3_baselines.py  -> WORKING_TREE when porcelain != ""
```

The width-independence argument is physical reasoning. I did not re-simulate.

## Impact

Wasted CI time on a slow 2D solve. The docstring overstates what the baseline guards.

## Proposed fix

- Remove `test_v3_cce_center_300um_preserved`, or replace it with an assertion that does depend on width
  (for example, the edge/center CCE ratio from `cce_lateral_scan`).
- Change the module docstring to "post-Plan-26-04 forward-guard baseline".
- Re-freeze from a clean tree so `git_commit_sha` is a real SHA, or make `freeze_v3_baselines.py` refuse
  to run on a dirty tree unless `--allow-dirty` is passed.

## Test plan

`uv run pytest tests/test_v3_baseline_regression.py -q` is green, with fewer slow 2D solves or a new
width-sensitive assertion. The JSON `git_commit_sha` matches `git rev-parse`.

## Effort

S
