---
id: C-04
source: raw-C F4
title: No serialization of devsim access; two concurrent facade calls (two Streamlit sessions or tabs) segfault the process or delete each other's devices
severity: high
category: bug
status: fixed  # 2026-10-18 sync: 5d20634 (file status already fixed)
location:
  - etna/api/simulation.py:102, ~203, 725 (`reset_devsim_fully()` at entry of run_cv / run_field / run_dark_current)
  - etna/core/devsim_reset.py:102-107 (deletes every device in `devsim.get_device_list()`)
  - all bucket-a facades (run_cce, run_radiation_damage, run_temperature_sweep, run_flash_recombination, run_transient), whose core functions build and delete their own devices
  - app/main.py (no lock; each page calls facades directly)
labels: [api, concurrency, crash, streamlit]
---

## Summary

Streamlit runs each browser session's script in its own thread inside one process. devsim
state is process-global and not thread-safe. `grep -rn "threading\|Lock(\|st.cache" app etna/api`
finds nothing. When two sessions click Run at overlapping times, session B's
`reset_devsim_fully()` deletes session A's device mid-solve. Worse, concurrent devsim calls
crash the interpreter with SIGSEGV, which kills the Streamlit server for every user.

The raw file had this as "plausible / not reproduced". I reproduced it, so the status is
upgraded to open and the severity raised to high.

## Reproduction / Evidence

Probe script `/private/tmp/claude-501/-Users-ngcex-projects-physics-petringa/df2eeea8-54e5-4d29-83f3-c7a1539b1919/scratchpad/c04_threads.py`
starts two threads, each calling `run_cv(DeviceConfig(epi_thickness_um=10|15), v_start=0, v_stop=-20, n_points=12)`.

```
for i in 1 2 3; do uv run python -I c04_threads.py > c04_run$i.txt 2>&1; echo "run$i exit=$?"; done
run1 exit=139
run2 exit=139
run3 exit=139            # 139 = SIGSEGV; no RESULT line printed
```

Earlier non-crashing runs of the same probe produced:

```
t0: RuntimeError 'Equilibrium solve failed to converge: DEVSIM FATAL: Unrecognized "direct_solver" parameter value "unknown"...'
t1: same
---
t0: devsim.error 'Device "device1d_15b2cf3b" does not exist...'     <- deleted by the other thread
t1: RuntimeError 'Equilibrium solve failed to converge: Convergence failure!'
```

Fix validation: the same probe with a module-level `threading.RLock()` around each `run_cv` call
(`c04_threads_locked.py`) gave 2/2 runs `exit=0` and
`RESULT {'t0': ('ok', 12, False), 't1': ('ok', 12, False)}`.

## Impact

- Any two overlapping solves crash the whole server, either from two tabs of one researcher or from LAN users (the server binds all interfaces, see C-16).
- Non-crashing interleavings yield wrong-device results or uncaught `devsim.error` (C-03).
- Library users with threads or Jupyter background tasks are equally exposed.

## Proposed fix

```diff
# etna/api/simulation.py
+import functools, threading
+_DEVSIM_LOCK = threading.RLock()   # devsim is process-global; one solver at a time
+
+def _serialized(fn):
+    @functools.wraps(fn)
+    def wrapper(*a, **kw):
+        with _DEVSIM_LOCK:
+            return fn(*a, **kw)
+    return wrapper
+
+@_serialized
 def run_cv(...):
+@_serialized
 def run_field(...):
 ... and EVERY devsim-touching facade: run_cce, run_radiation_damage,
     run_temperature_sweep, run_flash_recombination, run_dark_current, run_transient
```

The bucket-a facades need the lock too. Their core functions' devices are wiped by another
thread's `reset_devsim_fully()`. An RLock lets ParametricSweep or nested calls re-enter. Export
the lock (or a `devsim_session()` context manager) for direct core users, and document
"one solver per process". In the app, optionally wrap calls in `st.spinner("Waiting for solver...")`.

## Test plan

- `tests/test_api_concurrency.py` (slow): start two threads running `run_cv(DeviceConfig(), n_points=3)`. Both must return without exception, and the test process must survive. Run it in a subprocess so a regression shows up as a non-zero exit rather than killing pytest.
- Repeat with `run_cv` plus `run_cce` to cover the bucket-a path.

## Effort

S (decorator plus test). Throughput is serialized, which is acceptable for this workload.
