---
id: D-16
title: Security sweep of scripts/, examples/, tests/, app/ — no findings (informational)
severity: info
category: security
status: rejected
location: scripts/freeze_v3_baselines.py:177-180; tests/test_radiation_damage.py:423; scripts/run_calibration_2d.py:187-221; app/workflows/microdosimetry.py
labels: [security, informational, no-action]
source: raw-D F16
---

## Summary

The raw report recorded a negative result: no `eval`/`exec`/`pickle`/`yaml.load`/`shell=True`/`os.system`
in this slice, `subprocess` called with list arguments, `ast.literal_eval` used for cache keys, and the only
upload surface is `st.file_uploader(type=["csv"])` → `pandas.read_csv`.
**Status rejected (as an issue), not because the claim is wrong.** My spot-check agrees with it, but a
negative result has no defect to track. It is kept here so the finding is not dropped silently.

## Reproduction/Evidence

```
$ grep -rnE "eval\(|exec\(|pickle|os\.system|yaml\.load|shell=True" scripts examples tests app | grep -v literal_eval
app/workflows/microdosimetry.py:19:no pickle), and st.file_uploader(type=["csv"]) restricts the accepted   (comment only)
$ sed -n 174-182 scripts/freeze_v3_baselines.py  -> subprocess.check_output([...list...])
```

`run_calibration_2d.py` rewrites `etna/core/device2d.py` in place with an anchored regex. This is
self-modifying by design and is a maintainability concern, not a security one. I did not re-check the
claims about secrets or tracked `__pycache__`. Dependency CVEs were out of scope.

## Impact

None.

## Proposed fix

None.

## Test plan

N/A.

## Effort

S (none required)
