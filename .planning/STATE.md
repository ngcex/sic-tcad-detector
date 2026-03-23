---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: Realistic Device Physics
status: unknown
last_updated: "2026-03-23T16:17:15.217Z"
progress:
  total_phases: 10
  completed_phases: 10
  total_plans: 23
  completed_plans: 23
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-22)

**Core value:** Predict how CCE in 4H-SiC detectors degrades under FLASH dose rates, providing the first TCAD-based explanation of plasma recombination effects in SiC dosimeters.
**Current focus:** Phase 11 — Dark Current Modeling

## Current Position

Phase: 11 of 12 (Dark Current Modeling)
Plan: 01 of 02
Status: Plan 11-01 complete
Last activity: 2026-03-23 — Completed 11-01 (Hurkx TAT dark current module)

Progress: [██████████] 100%

## Performance Metrics

**Velocity:**

- Total plans completed: 4 (v1.1)
- Average duration: 8.5 min
- Total execution time: 0.57 hours

**By Phase:**

| Phase | Plans | Total  | Avg/Plan |
| ----- | ----- | ------ | -------- |
| 10    | 3/3   | 13 min | 4.3 min  |
| 11    | 1/2   | 21 min | 21 min   |

**Recent Trend:**

- Last 5 plans: 10-01 (3 min), 10-02 (5 min), 10-03 (5 min), 11-01 (21 min)
- Trend: 11-01 took longer due to model reformulation (n_i^2 bottleneck discovery)

_Updated after each plan completion_

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [Roadmap]: 3 phases derived from 21 requirements — Temperature (9 reqs) -> Dark Current (6 reqs) -> Transient (6 reqs)
- [Roadmap]: Notebooks bundled with their physics phase rather than separated into a combined-analysis phase
- [10-01]: Calibrated E_g_0 to 3.2965625 eV so Varshni gives exactly 3.26 eV at 300K
- [10-01]: Used calibration factor n_i(T) = n_i_300 \* compute_ni(T)/compute_ni(300) to anchor n_i
- [10-02]: Used \_UNSET sentinel pattern in charge_collection.py for backward-compatible T-dependent defaults
- [10-02]: hecht_cce T-dependent defaults use low-doping limit mu_max\*(T/300)^gamma since Hecht assumes bulk mobility
- [10-03]: Used Hecht method as default CCE sweep (faster than DD, adequate for T-sensitivity analysis)
- [10-03]: Pandas DataFrames for sweep results (long format for CCE, wide format for IV)
- [11-01]: Used effective generation rate N_t=2.2e13 cm^-3/s (not physical trap density) because n_i^2 bottleneck prevents pA-level dark current in 1D SiC
- [11-01]: E-field-weighted depletion region selector (E/E_ref clamped to 1) for voltage-dependent generation
- [11-01]: Gamma=1 at SiC detector fields is correct physics; Hurkx enhancement requires MV/cm fields

### Pending Todos

None.

### Blockers/Concerns

- [Phase 10]: RESOLVED in 10-01 — n_i(300K) discrepancy handled via calibration factor: n_i(T) = n_i_300 \* compute_ni(T)/compute_ni(300).
- [Phase 11]: Dark current mechanism ambiguity — 18 pA may be perimeter leakage (inherently 2D, unmodellable in 1D). TAT with effective parameters is the fallback.
- [Phase 12]: Transient computational cost unverified — estimated ~30s/pulse but adaptive dt across 6-order timescale gap may be expensive.

## Session Continuity

Last session: 2026-03-23
Stopped at: Completed 11-01-PLAN.md
Resume file: None
