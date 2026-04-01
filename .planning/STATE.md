---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: SiC Microdosimeter Design Study
status: in-progress
last_updated: "2026-04-01T07:56:28.000Z"
progress:
  total_phases: 20
  completed_phases: 17
  total_plans: 39
  completed_plans: 38
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-27)

**Core value:** TCAD-based feasibility study for a novel 4H-SiC microdosimeter -- first open-source 2D simulation with microdosimetric spectra computation and design optimization guidance
**Current focus:** Phase 24 - Alternative Structures

## Current Position

Phase: 24 (6 of 8 in v3.0) — Alternative Structures [IN PROGRESS]
Plan: 1 of 2 in current phase [COMPLETE]
Status: Plan 24-01 complete -- 4 alternative structure mesh builders with Poisson validation
Last activity: 2026-04-01 — Plan 24-01 executed, alternative_structures.py with mesa/3D/delta-E/E/guard ring

Progress: [████████░░] ~82% (v3.0, 12 of ~15 plans)

## Performance Metrics

**Velocity:**

- Total plans completed: 52 (v1.0: 20, v1.1: 7, v2.0: 13, v3.0: 12)
- Average duration: ~14 min
- Total execution time: ~10 hours

**Recent Trend:**

- v2.0 plans: 13 plans in 3 days
- Trend: Stable

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [v3.0]: New device2d.py module; frozen device.py to protect 14 validated notebooks
- [v3.0]: devsim physics modules are dimension-agnostic (poisson, drift-diffusion, transient, CCE)
- [v3.0]: Two new dependencies: gmsh (>=4.15.1) for 2D mesh, uproot (>=5.6) for Geant4 ROOT files
- [v3.0]: Graded epi doping must be extended to 2D (lateral uniformity, correct junction position)
- [v3.0]: CCE(LET) lookup table pattern: 30-50 TCAD transients, then apply to 1000+ MC events
- [19-01]: x=lateral, y=depth coordinate convention for all 2D modules
- [19-01]: Air buffer regions use SiC material for devsim contact detection
- [19-01]: Graded doping default for 2D (same calibrated defaults as 1D)
- [19-02]: E-field 2D visualization via scipy LinearNDInterpolator + np.gradient on regular grid
- [19-02]: validate_2d_vs_1d uses center-slice potential gradient for E-field comparison
- [19-02]: poisson.py confirmed dimension-agnostic (works on 1D and 2D without modification)
- [20-01]: \_robust_dc_solve with fallback to relaxed tolerances for 2D DD convergence
- [20-01]: Delete 2D device before creating 1D device to avoid devsim global solver coupling
- [21-01]: charge_error=1e10 required for all devsim BDF1 transient solves (disables step rejection)
- [21-01]: Robust transient_dc init with fallback to relaxed tolerances for 2D devices
- [21-01]: Generation-pulse injection: G/dt_inject for one BDF1 step, then zero gen and collect
- [21-02]: 300um SV LET sweep reduced to 10 points (from 15) due to ~3x longer sim time on larger mesh
- [21-02]: CCE ~1.01 at 50V full depletion is numerically valid (1% overcollection from generation-pulse method)
- [22-01]: column_map keys are standard names, values are source names -- consistent across CSV and ROOT loaders
- [22-01]: uproot imported lazily inside ROOT functions so module works without uproot for CSV-only workflows
- [22-01]: process_mc_ensemble filters zero-energy events with logging rather than raising errors
- [22-02]: Demo CCE curve added for visual comparison when real CCE table is flat at full depletion
- [22-02]: 2000 synthetic events with bimodal distribution (proton + heavy-ion) for mixed-field demonstration
- [23-01]: PSTAR water stopping powers bundled as CSV (37 points, 0.1-1000 MeV) for kappa computation
- [23-01]: SiC stopping powers from Bethe-Bloch scaling of SRIM data, ~1.7x water values
- [23-01]: Constant kappa fallback (0.58) with warning when no energy-dependent table provided
- [23-02]: 22-cell notebook structure for full microdosimetric pipeline: MC events -> lineal energy -> spectra -> tissue equivalence
- [23-02]: 4 publication-quality figures: y*d(y), y*f(y), SiC vs tissue overlay, bar chart comparison
- [24-01]: devsim cannot have contacts at interface boundary -- delta-E/E uses 2 contacts with interface continuity
- [24-01]: Guard ring acceptor doping via explicit expression to avoid devsim cyclic model dependency
- [24-01]: Cylindrical coordinate lifecycle: create, activate, physics, delete, restore_cartesian_coords

### Pending Todos

None.

### Blockers/Concerns

- Geant4 TTree naming conventions from INFN-LNS group unknown -- need sample ROOT file before Phase 22
- devsim cylindrical coordinate API has limited community examples -- verify at Phase 19 start
- SiC-specific kappa tissue-equivalence factor computed from PSTAR/SRIM stopping power ratios in Phase 23-01

## Session Continuity

Last session: 2026-04-01
Stopped at: Completed 24-01-PLAN.md
Resume file: None
