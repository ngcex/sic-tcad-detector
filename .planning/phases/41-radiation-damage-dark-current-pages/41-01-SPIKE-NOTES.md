# Phase 41 Plan 01 — Live-devsim Spike Notes

Recorded 2026-07-13. Both checks run against a live devsim build via `uv run python -c "..."`
(throwaway script, not committed). These findings are load-bearing for Wave 2 (41-02, 41-03) —
use the confirmed values below instead of re-deriving or re-testing.

## Check A: `ParametricSweep(param="T")` + `run_dark_current` at a minimal single-point bias

**Confirmed working, no exception, on the first attempt** — the UI-SPEC's assumed `n_points=1`
with `v_start == v_stop` did NOT need the `n_points=2` fallback described in the plan's Task 1
contingency.

Exact call that succeeded:

```python
petringa.ParametricSweep(
    base_config=petringa.DeviceConfig(),
    param="T",
    values=[250.0, 325.0, 400.0],
    sim_fn=petringa.run_dark_current,
    sim_kwargs={"v_start": -20.0, "v_stop": -20.0, "n_points": 1},
).run()
```

Result: `list` of length 3, no exception raised. Each element is a `SimResult` with
`sim_type="dark_current"`, `len(result.x) == 1` (i.e. `result.x == [-20.0]`), and
`result.metadata["I_SRH"]` / `["I_TAT"]` / `["I_SRV"]` each of length 1 (matching `result.x`'s
length). Observed `y` values (`I_total`, one per temperature): `2.076e-17` (250 K),
`3.203e-11` (325 K), `3.080e-06` (400 K) — physically sane monotonic increase with temperature.

**Confirmed safe kwargs for Wave 2's dark current page:** `v_start=-20.0, v_stop=-20.0,
n_points=1` (matches the UI-SPEC's `V_bias=-20.0` operating-bias default exactly — no amendment
needed). Because `n_points=1`, each per-temperature `SimResult` has arrays of length 1; Wave 2's
aggregation step must read `result.y[0]` (and `result.metadata["I_SRH"][0]` /
`["I_TAT"][0]` / `["I_SRV"][0]`) per temperature to build the aggregated length-`n_temperatures`
`SimResult`, NOT assume any other index or length.

## Check B: `run_radiation_damage(DeviceConfig(), V_bias=-20.0)`

**No exception raised.** Returned a full 6-element `result.y` array (matching the facade's
default `fluences=np.geomspace(1e13, 1e16, 6)`).

`result.x` (fluences, cm^-2): `[1.000e13, 3.981e13, 1.585e14, 6.310e14, 2.512e15, 1.000e16]`

`result.y` (CCE): `[0.9928, nan, 0.9706, 0.9208, 0.7564, 0.4191]`

**One NaN present**, at index 1, `fluence ≈ 3.98107e13`. This DISCONFIRMS RESEARCH.md's
Assumption A3 — the shallower `V_bias=-20.0` (vs. the facade's own default `V_bias=-40.0`) does
**not** avoid the NaN; it still occurs at the same fluence value observed at `-40.0` in the
research session. `metadata = {"V_bias": -20.0, "energy_MeV": 5.6}` — no `"truncated"` key, no
other anomaly.

**Decision for Wave 2: keep `V_bias=-20.0` as the confirmed default anyway.** The single NaN at
fluence≈3.98e13 is an accepted, expected outcome — `build_damage_figure` is explicitly designed
to be NaN-tolerant (per its `<action>` spec: do not call `.dropna()`; Plotly renders a NaN
y-value as a native line gap). Do not chase a bias value that "fixes" this NaN — none was tested
that eliminates it, and per RESEARCH.md Pitfall 2 item 2, tolerant rendering is an acceptable
fallback for a single missing point, not a bug to fix. The UI-SPEC's copywriting contract's
partial-failure info banner (`st.info(...)`, shown when `np.any(np.isnan(result.y))`) is the
correct, already-planned handling for this — 41-02 should implement that banner and move on.

## Summary for Wave 2

| Page                             | Confirmed default(s)                                                                                                   | Source                                                                                    |
| -------------------------------- | ---------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| Radiation damage (41-02)         | `V_bias=-20.0` (UI-SPEC default, unchanged)                                                                            | Check B — NaN persists at this value; accepted via NaN-tolerant rendering per plan design |
| Dark current (41-03)             | `V_bias=-20.0` fixed operating bias; per-temperature `sim_kwargs={"v_start": V_bias, "v_stop": V_bias, "n_points": 1}` | Check A — confirmed working with no fallback needed                                       |
| Dark current (41-03) aggregation | Each per-T `SimResult` has `len(x) == len(metadata["I_*"]) == 1`; use index `[0]`                                      | Check A                                                                                   |

No scratch Python script was committed — both checks were run via `uv run python -c "..."`
inline, not saved to a file in the repository.
