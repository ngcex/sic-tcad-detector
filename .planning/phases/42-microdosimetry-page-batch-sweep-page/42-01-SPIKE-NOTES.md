# Phase 42 Plan 01 — Live-devsim Spike Notes

Recorded 2026-07-14. The single check below was run against a live devsim build via
`uv run python -c "..."` (throwaway one-liner, NOT committed as a test — devsim CCE solves are
too slow/resource-heavy for the pytest suite per STATE.md). This finding is load-bearing for
Wave 2's batch sweep page (42-03) — it reads this file for its hardcoded default instead of
re-deriving or re-running the sweep.

Environment note: the benign `Skipping libblas.dylib` / `loading UMFPACK 5.1 as direct solver`
lines appear at import (the libopenblas/libblas DLL warning from RESEARCH Environment
Availability) — devsim still solves through it, as confirmed by this spike returning full curves.

## Check A1: `ParametricSweep(param="epi_thickness_um", values=[10,15,20])` + `run_cce`

**Confirmed working, no exception, on the first attempt.** All three swept values returned a
full (non-truncated) CCE curve.

Exact call that succeeded:

```python
petringa.ParametricSweep(
    base_config=petringa.DeviceConfig(),
    param="epi_thickness_um",
    values=[10, 15, 20],
    sim_fn=petringa.run_cce,
).run()
```

Result: `list` of length 3, no exception raised. Each element is a `SimResult` with
`sim_type="cce"`. Per-value observed shape:

| epi_thickness_um | x length | y length | truncated | x-range (V)    | renderable? |
| ---------------- | -------- | -------- | --------- | -------------- | ----------- |
| 10               | 30       | 30       | `False`   | [-40.0, -10.0] | ✅ full     |
| 15               | 30       | 30       | `False`   | [-40.0, -10.0] | ✅ full     |
| 20               | 30       | 30       | `False`   | [-40.0, -10.0] | ✅ full     |

**Renderable / truncated / empty count: 3 renderable, 0 truncated, 0 empty.**

None of the three results set `metadata["truncated"]` — the default `run_cce` `v_stop=-40` stays
inside the convergence envelope for all three epi thicknesses (10, 15, 20 µm). Sweeping epi
_thicker_ moves punch-through to deeper bias, so the thinnest value (10 µm) is the closest to the
wall and it still solved cleanly to -40 V with no truncation. No `run_field` was involved (it is
excluded from the batch-sweep sim-type selectbox per RESEARCH Pitfall 3 — it raises rather than
truncating).

## VERDICT

**The locked batch-sweep default is CONFIRMED — no revision needed.** 42-03 may hardcode:

- `sim_label = "CCE vs bias (run_cce)"`
- `param = "epi_thickness_um"`
- `values = "10, 15, 20"` (parsed to `[10, 15, 20]`)

This combination over the default `DeviceConfig()` yields 3 fully renderable overlay curves (30
points each, x-range -40 to -10 V, none truncated). `build_sweep_overlay_figure(results,
"epi_thickness_um", [10, 15, 20], "CCE vs bias (run_cce)")` will draw three complete
CCE-vs-bias traces with x-axis "Bias V (V)" and y-axis "Charge Collection Efficiency". Keep the
`try/except RuntimeError → st.error` guard from the Phase 39 pages as defence-in-depth, but no
truncation banner is expected to fire for this default.

## Summary for Wave 2

| Page                | Confirmed default(s)                                                                 | Source   |
| ------------------- | ------------------------------------------------------------------------------------ | -------- |
| Batch sweep (42-03) | `sim_label="CCE vs bias (run_cce)"`, `param="epi_thickness_um"`, `values=[10,15,20]` | Check A1 |
| Batch sweep (42-03) | Each result: `len(x)==len(y)==30`, `truncated=False`, x-range [-40, -10] V           | Check A1 |

No scratch Python script was committed — the spike was run via `uv run python -c "..."` inline,
not saved to a file in the repository.
