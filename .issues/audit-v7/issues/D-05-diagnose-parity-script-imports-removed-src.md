---
id: D-05
title: scripts/diagnose_1d_2d_parity.py imports removed `src` package unguarded and crashes in run_1d()
severity: medium
category: bug
status: open
location: scripts/diagnose_1d_2d_parity.py:137 (guarded), :196 (unguarded), :269 (call from run_1d), :363-371 (only devsim.error caught)
labels: [scripts, stale-import, src-rename]
source: raw-D F5
---

## Summary

`_W_at_voltages_1d()` runs `from src import cv_analysis as _cv` (line 196) with no try/except. `run_1d()`
calls it at line 269 after the slow 1D ramp. `main()` only catches `devsim.error`, so the script dies
with `ModuleNotFoundError`. Line 137 has the same import, but it is wrapped in `except Exception` with a
fallback that happens to be correct (`"cathode_bias"`).

Correction to the raw report: `.planning/phases/26-graded-doping-2d-calibration/26-DIAGNOSIS.md` **already
exists**, so `scripts/run_calibration_2d.py` is not blocked today. The diagnosis just cannot be
regenerated. For that reason severity is downgraded to medium (raw: high).

## Reproduction/Evidence

```
$ ls -d src                  -> No such file or directory
$ grep -n "_W_at_voltages_1d\|def run_1d" scripts/diagnose_1d_2d_parity.py
188:def _W_at_voltages_1d(voltages):
252:def run_1d():
269:    W = _W_at_voltages_1d([-10.0, -30.0])
$ sed -n 192,198p ...   -> `from src import cv_analysis as _cv` (no try)
$ sed -n 363,372p ...   -> `except devsim.error as e:` only
$ ls .planning/phases/26-graded-doping-2d-calibration/26-DIAGNOSIS.md  -> exists
```

I did not run the full script because the ramp is slow. The crash path follows from the code above.

## Impact

The documented diagnostic workflow (`run_calibration_2d.py` tells users to run this script first) is
broken. A rerun after any change to the doping defaults is impossible without editing the script.

## Proposed fix

```diff
-        try:
-            from src import cv_analysis as _cv  # for GetContactBiasName
-            bias_name = _cv.simple_physics.GetContactBiasName("cathode")
-        except Exception:
-            bias_name = "cathode_bias"
+        from etna.core import cv_analysis as _cv
+        bias_name = _cv.simple_physics.GetContactBiasName("cathode")
 ...
-    from src import cv_analysis as _cv
+    from etna.core import cv_analysis as _cv
```

Also drop the unused `numpy` import (line 48, see D-09) and fix the stale `cd .../etna` usage line (see D-11).

## Test plan

Add a cheap smoke test: `importlib` the script module and assert that `_W_at_voltages_1d` resolves its
import (or run `python -c "import ast; ..."` to check there are no `from src` imports anywhere in `scripts/`).
`grep -rn "from src" scripts/ etna/ tests/` should be empty.

## Effort

S
