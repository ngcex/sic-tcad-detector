---
id: D-11
title: Stale src/ module names, /physics/etna cd paths, "46 pytest modules" count and xfail-skeleton comments
severity: low
category: docs
status: partial  # 2026-10-18 sync: 7acbe10 (wp-c1): scripts/ usage lines and src/ refs fixed (remaining 'src/device2d.py' mentions are historical, past tense); open: README.md:123 (WP-D), tests/ comments (WP-C2), device2d.py refs
location: README.md:123; scripts/freeze_v3_baselines.py:23; scripts/diagnose_1d_2d_parity.py:27; src/-references in scripts/freeze_v3_baselines.py, scripts/run_calibration_2d.py, scripts/regression_sweep_v3_notebooks.py, tests/test_devsim_reset.py, tests/test_v3_baseline_regression.py, tests/test_plotting.py, tests/test_device2d.py, etna/core/devsim_reset.py, etna/core/device2d.py; tests/test_device2d.py:196-212; tests/test_v3_baseline_regression.py:7
labels: [docs, comment-rot, src-rename]
source: raw-D F11
---

## Summary

Comments and docs were not updated after the `src` → `etna.core` rename and the move of the project to
`petringa`:

- README.md:123 says "46 pytest modules", but there are 45.
- The usage lines in `freeze_v3_baselines.py:23` and `diagnose_1d_2d_parity.py:27` say
  `cd /Users/ngcex/projects/physics/etna`, which does not exist.
- Docstrings and comments in 9 files still name `src/...` modules.
- `test_device2d.py:196-212` and `test_v3_baseline_regression.py:7` describe `@pytest.mark.xfail(strict=True)`
  skeletons that no longer exist. The only `xfail` matches in the repo are these comments.

None of this affects runtime except D-05, where the stale import is executable code.

## Reproduction/Evidence

```
$ ls tests/test_*.py | wc -l            -> 45
$ sed -n 123p README.md                 -> "46 pytest modules — ..."
$ grep -n "physics/etna" scripts/*.py   -> diagnose_1d_2d_parity.py:27, freeze_v3_baselines.py:23 (+ create_notebook_15.py:479, see D-06)
$ grep -rln "src/\|src\." <listed files> -> all 9 files match
$ grep -rn xfail tests/                 -> only comment lines in test_device2d.py:199-212, test_v3_baseline_regression.py:7
```

## Impact

Misleading for maintainers. A wrong count in a document a foundry may read.

## Proposed fix

- README: drop the number ("pytest modules — ...") or change it to 45.
- Replace both `cd` lines with `cd <repo root> && uv run python scripts/...`.
- In comments and docstrings, replace `src/` with `etna/core/` and `src.` with `etna.core.` in the listed files (review each match; do not blanket-sed code).
- Delete the xfail-skeleton comment blocks, or rewrite them in the past tense.

## Test plan

`grep -rn "physics/etna\|from src\|src/" scripts tests etna README.md` returns only intentional matches.
The test suite is unchanged.

## Effort

S
