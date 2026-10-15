---
id: D-14
title: cwd-relative data paths in three tests silently skip (or fail) when pytest is not launched from repo root
severity: low
category: improvement
status: open
location: tests/test_v3_baseline_regression.py:23,27-28; tests/test_mc_coupling.py:426-428; tests/test_app_microdosimetry_page.py:91
labels: [tests, paths, silent-skip]
source: raw-D F14
---

## Summary

Three tests build data paths relative to the current working directory. Two of them turn a missing
committed fixture into `pytest.skip`, so the regression guards disappear quietly when the suite runs from
another directory or an IDE runner. `tests/test_api_microdosimetry.py` already does this correctly with
`Path(__file__).parent.parent / "data"`.

## Reproduction/Evidence

```
$ grep -n 'Path("' tests/test_v3_baseline_regression.py tests/test_mc_coupling.py tests/test_app_microdosimetry_page.py
tests/test_app_microdosimetry_page.py:91:    data = Path("data/synthetic_mc_events.csv").read_bytes()
tests/test_v3_baseline_regression.py:23:BASELINE_PATH = pathlib.Path("tests/baselines/v3_frozen.json")
tests/test_mc_coupling.py:426:        cce_path = Path("data/cce_let_table_100um.json")
test_v3_baseline_regression.py:27-28  -> if not BASELINE_PATH.exists(): pytest.skip(...)
test_mc_coupling.py:427-428           -> if not cce_path.exists(): pytest.skip(...)
```

## Impact

Guards can silently switch off. The third test fails with a confusing `FileNotFoundError` instead.

## Proposed fix

```diff
-BASELINE_PATH = pathlib.Path("tests/baselines/v3_frozen.json")
+BASELINE_PATH = pathlib.Path(__file__).resolve().parent / "baselines" / "v3_frozen.json"
 ...
-    if not BASELINE_PATH.exists():
-        pytest.skip(f"Baseline file not found at {BASELINE_PATH}")
+    assert BASELINE_PATH.exists(), f"committed baseline missing: {BASELINE_PATH}"
```

Apply the same change to `test_mc_coupling.py:426` (`parents[1] / "data" / ...`, assert instead of skip)
and `test_app_microdosimetry_page.py:91`.

## Test plan

`cd /tmp && uv run --project <repo> pytest <repo>/tests/test_v3_baseline_regression.py <repo>/tests/test_mc_coupling.py <repo>/tests/test_app_microdosimetry_page.py -q`
runs with no skips caused by missing files.

## Effort

S
