---
id: D-04
title: test_full_suite_collects_with_module_present always fails — subprocess cwd hardcoded to nonexistent /Users/ngcex/projects/physics/etna
severity: medium
category: bug
status: fixed  # 2026-10-18 sync: fd81f26: stale cwd in the collect meta-test fixed
location: tests/test_radiation_damage.py:423-440
labels: [tests, stale-path, always-red]
source: raw-D F4
---

## Summary

The meta-test runs `pytest tests/ --collect-only` with `cwd="/Users/ngcex/projects/physics/etna"`. The
repo lives in `.../physics/petringa`, so `subprocess.run` raises `FileNotFoundError` before pytest starts.
The guard against breaks in the import chain can never pass.
Severity is set to medium (raw: high). This is one permanently failing test with a one-line fix, and
the collection behaviour it guards is still exercised by every normal run.

## Reproduction/Evidence

```
$ ls -d /Users/ngcex/projects/physics/etna
ls: /Users/ngcex/projects/physics/etna: No such file or directory
$ uv run pytest "tests/test_radiation_damage.py::TestRegressionSafety::test_full_suite_collects_with_module_present" -q
/opt/miniconda3/lib/python3.13/subprocess.py:1972: FileNotFoundError
FAILED ... 1 failed in 0.88s
```

The path is a string literal with no marker, skip or environment override.

## Impact

The suite is always red, which teaches people to ignore failures. The intended guard does nothing.

## Proposed fix

```diff
+from pathlib import Path
 ...
-            cwd="/Users/ngcex/projects/physics/etna",
+            cwd=str(Path(__file__).resolve().parents[1]),
```

## Test plan

Run the node id above and expect it to pass. Running it from another directory (`cd /tmp && pytest <abs path>`) should also pass.

## Effort

S
