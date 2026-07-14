---
phase: 42-microdosimetry-page-batch-sweep-page
plan: 01
subsystem: ui
tags:
  [
    plotly,
    pandas,
    csv-export,
    streamlit,
    microdosimetry,
    parametric-sweep,
    devsim,
  ]

# Dependency graph
requires:
  - phase: 41-radiation-damage-dark-current-pages
    provides: app/components/results.py pure builder/serializer blueprint (build_damage_figure, build_dark_current_figure, to_csv_bytes damage/dark_current branches) + spike-notes precedent
provides:
  - build_microdosimetry_figure — pure single-trace log-x Plotly builder (y·d(y) vs lineal energy)
  - build_sweep_overlay_figure — pure 4-arg overlay builder with per-facade axis-title map (_SWEEP_AXIS_TITLES) and value-keyed palette (_SWEEP_PALETTE)
  - sweep_results_to_csv_bytes — standalone bulk serializer emitting N results as ONE CSV with a leading param run-identifier column
  - to_csv_bytes microdosimetry branch — y_keV_per_um/y_times_d_y columns + y_F/y_D/l_bar_um header
  - 42-01-SPIKE-NOTES.md — confirmed convergence-safe batch-sweep default (run_cce + epi_thickness_um=[10,15,20]) for the Wave 2 batch sweep page
affects: [42-02, 42-03, microdosimetry-page, batch-sweep-page]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "4-arg overlay builder: sim_label selects axis titles from a per-facade map while param names legend/run-id (axes and legend decouple by design)"
    - "Bulk multi-run CSV is a SEPARATE function (sweep_results_to_csv_bytes), NOT a to_csv_bytes branch — to_csv_bytes stays single-result / sim_type-dispatched"
    - "Spike-first: confirm live-devsim convergence before a downstream page hardcodes the default"

key-files:
  created:
    - .planning/phases/42-microdosimetry-page-batch-sweep-page/42-01-SPIKE-NOTES.md
  modified:
    - app/components/results.py
    - tests/test_app_csv_export.py

key-decisions:
  - "build_sweep_overlay_figure ships the 4-arg (results, param, values, sim_label) signature — the 3-arg RESEARCH.md sketch is superseded by UI-SPEC's per-facade axis map"
  - "Batch-sweep default CONFIRMED unchanged (run_cce + epi_thickness_um=[10,15,20]): live spike returned 3 full 30-point curves, none truncated, so no thicker-epi revision needed"
  - "microdosimetry to_csv_bytes header keys are y_F_keV_per_um/y_D_keV_per_um/l_bar_um but read from metadata keys y_F/y_D/l_bar_um (key name != metadata key)"

patterns-established:
  - "Pattern 1: per-facade axis-title lookup via module-level _SWEEP_AXIS_TITLES dict keyed by selectbox label, with a defensive-only fallback"
  - "Pattern 2: value-keyed qualitative palette cycled by trace order (i % len), distinct from Phase 41's fixed-quantity color mapping"

requirements-completed: [FEAT-03, FEAT-04]

# Metrics
duration: 8 min
completed: 2026-07-14
---

# Phase 42 Plan 01: Microdosimetry + Batch-Sweep results.py Foundation Summary

**Two new pure Plotly builders (log-x microdosimetry spectrum + 4-arg parametric overlay), a standalone bulk-sweep CSV serializer, and a microdosimetry to_csv_bytes branch in app/components/results.py — plus a live-devsim spike confirming the run_cce + epi_thickness_um=[10,15,20] batch-sweep default renders 3 full curves.**

## Performance

- **Duration:** 8 min
- **Started:** 2026-07-14T00:00:00Z (approx)
- **Completed:** 2026-07-14
- **Tasks:** 3
- **Files modified:** 3 (2 modified, 1 created)

## Accomplishments

- Added `build_microdosimetry_figure` (single log-x trace, color #1F6FEB) and `build_sweep_overlay_figure` (4-arg, per-facade axis map, value-keyed palette) — both pure, no Streamlit runtime
- Added `sweep_results_to_csv_bytes` bulk serializer (N results → one CSV, leading `<param>` run-id column) and a `microdosimetry` branch in `to_csv_bytes`, leaving the five existing branches and the shared header untouched
- Extended `tests/test_app_csv_export.py` with two pure tests (9 pass total)
- Ran the live-devsim A1 spike: `ParametricSweep(param="epi_thickness_um", values=[10,15,20], sim_fn=run_cce)` over the default `DeviceConfig()` returned 3 full non-truncated CCE curves — batch-sweep default CONFIRMED for 42-03

## Task Commits

Each task was committed atomically:

1. **Task 1: Two pure figure builders** - `3abea6d` (feat)
2. **Task 2: Bulk-CSV serializer + microdosimetry to_csv_bytes branch** - `254f909` (feat)
3. **Task 3: Extend CSV tests + batch-sweep spike** - `218523b` (test) + `177935e` (docs — SPIKE-NOTES.md)

_Task 3 is a tdd task, but its implementation (the serializer + branch) already landed in Task 2 by plan design — so its test commit is test-after and passed immediately (no RED phase); the spike output is committed as docs._

## Files Created/Modified

- `app/components/results.py` - Added `build_microdosimetry_figure`, `build_sweep_overlay_figure` (+ `_SWEEP_PALETTE`, `_SWEEP_AXIS_TITLES` constants), `sweep_results_to_csv_bytes`, and a `microdosimetry` branch in `to_csv_bytes`
- `tests/test_app_csv_export.py` - Added `test_microdosimetry_csv_columns_and_header` + `test_sweep_results_to_csv_bytes_shape`; imported `sweep_results_to_csv_bytes`
- `.planning/phases/42-microdosimetry-page-batch-sweep-page/42-01-SPIKE-NOTES.md` - Recorded the confirmed convergence-safe batch-sweep default for 42-03

## Decisions Made

- Shipped the 4-arg `build_sweep_overlay_figure(results, param, values, sim_label)` signature per UI-SPEC (superseding RESEARCH.md's 3-arg sketch), since axis titles require knowing which facade ran.
- Kept the batch-sweep default (`run_cce`, `epi_thickness_um=[10,15,20]`) unchanged — the live spike showed all 3 values solve fully to -40 V with no truncation, so no thicker-epi revision was needed.
- Microdosimetry CSV header keys (`y_F_keV_per_um`, `y_D_keV_per_um`, `l_bar_um`) intentionally differ from the source metadata keys (`y_F`, `y_D`, `l_bar_um`) — the header labels carry explicit units while the metadata keys stay as the library emits them.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None. A repository formatter (PostToolUse hook) reformatted `results.py` and `42-01-SPIKE-NOTES.md` after edits (collapsed a multi-line DataFrame call, cosmetic only) — no functional impact, all verifications green.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- `app/components/results.py` foundation is complete and tested; both Wave 2 page plans (42-02 microdosimetry, 42-03 batch sweep) can build against it without duplicating figure/serializer code.
- 42-03 has a confirmed, convergence-safe default in `42-01-SPIKE-NOTES.md` — no live re-testing needed before hardcoding it.

## Self-Check: PASSED

- Files exist: app/components/results.py, tests/test_app_csv_export.py, 42-01-SPIKE-NOTES.md — all FOUND
- Commits exist: 3abea6d, 254f909, 218523b, 177935e — all FOUND
- `uv run pytest tests/test_app_csv_export.py -q` → 9 passed
- Module purity: no `st.*` calls in results.py

---

_Phase: 42-microdosimetry-page-batch-sweep-page_
_Completed: 2026-07-14_
