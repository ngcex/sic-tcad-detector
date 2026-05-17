---
gsd_state_version: 1.0
milestone: v4.0
milestone_name: Scientific Validation & Extended Physics
status: planning
last_updated: "2026-05-17T00:00:00.000Z"
progress:
  total_phases: 0
  completed_phases: 0
  total_plans: 0
  completed_plans: 0
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-05-17)

**Core value:** TCAD-based feasibility study per microdosimetro 4H-SiC — prima simulazione 2D open-source con spettri microdosimetrici, ottimizzazione parametrica, e validazione scientifica paper-ready
**Current focus:** Defining requirements for v4.0

## Current Position

Phase: Not started (defining requirements)
Plan: —
Status: Defining requirements
Last activity: 2026-05-17 — Milestone v4.0 started

Progress: [░░░░░░░░░░] 0%

## Performance Metrics

**Velocity (storico):**

- Total plans completed: 55 (v1.0: 20, v1.1: 7, v2.0: 13, v3.0: 15)
- Average duration: ~14 min per plan
- Total execution time: ~10 hours

## Accumulated Context

### Decisions (v3.0 — portate avanti)

- New device2d.py module; frozen device.py to protect 20 validated notebooks
- devsim physics modules are dimension-agnostic (poisson, drift-diffusion, transient, CCE)
- Dependencies: gmsh (>=4.15.1) per 2D mesh, uproot (>=5.6) per Geant4 ROOT files
- CCE(LET) lookup table pattern: 30-50 TCAD transients → apply to 1000+ MC events
- x=lateral, y=depth coordinate convention per tutti i moduli 2D
- charge_error=1e10 required per tutti i BDF1 transient solves (disables step rejection)
- uproot importato lazily per backward-compatibility con workflow CSV-only
- Guard ring raccomandato come primo upgrade pratico per il gruppo Petringa

### Tech Debt da Risolvere in v4.0

- N_D uniforme in 2D fallisce a reverse bias — serve profilo graded epi in 2D (vedi device2d.py)
- ROOT/uproot integration mock-only — serve file sintetico Geant4-compatibile come fixture reale
- Kappa da scaling analitico Bethe-Bloch — serve tabella PSTAR+SRIM tabulata per confronto quantitativo
- CCE(LET) flat a 1.0 — valid physics (full depletion) ma necessita validazione umana prima pubblicazione
- t_collection anomalously fast (0.03 ns) — possibile artifact early termination, da verificare
- score_structures usa metriche hardcoded invece di output TCAD live

### Pending Todos

None.

### Blockers/Concerns

- Geant4 TTree naming conventions da INFN-LNS sconosciute — usiamo file sintetico in v4.0
- devsim 3D mesh API poco documentata — richede ricerca prima della fase 3D

## Session Continuity

Last session: 2026-05-17
Stopped at: Starting v4.0 milestone definition
Resume file: None
