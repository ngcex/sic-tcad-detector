# Research Summary — v4.0 Scientific Validation & Extended Physics

**Project:** 4H-SiC TCAD Simulator — Petringa Group (INFN-LNS)
**Milestone:** v4.0 — Scientific Validation & Extended Physics
**Synthesized:** 2026-05-17
**Sources:** STACK.md, FEATURES.md, ARCHITECTURE.md, PITFALLS.md

---

## Stack Additions

One new Python package is required. Everything else is custom code on the existing stack.

| Package / File | Version            | Purpose                                                        | Action                     |
| -------------- | ------------------ | -------------------------------------------------------------- | -------------------------- |
| `physdata`     | `>=0.2.0` (~50 KB) | Fetch NIST PSTAR proton stopping power tables programmatically | `uv add "physdata>=0.2.0"` |

Existing stack unchanged: devsim >=2.10, gmsh >=4.15.1, uproot >=5.6, numpy >=1.24, scipy >=1.11, matplotlib >=3.7.

**Vendored data files** (commit to repo):

| File                                   | Source                               | Notes                                                |
| -------------------------------------- | ------------------------------------ | ---------------------------------------------------- |
| `data/srim/sic_proton.txt`             | SRIM/TRIM one-time run               | Proton stopping power in SiC                         |
| `data/srim/water_proton.txt`           | SRIM/TRIM one-time run               | Proton stopping power in water                       |
| `tests/fixtures/synthetic_geant4.root` | `scripts/generate_geant4_fixture.py` | uproot-written flat ntuple for ROOT integration test |

**What NOT to add:** PyROOT (500 MB), pysrim (Windows .exe wrapper), pyvista (200 MB), pandas, FEniCS/FiPy.

**SiC stopping power note:** PSTAR does NOT include SiC as a named compound. Use Bragg additivity (PSTAR Si + PSTAR C, mass-weighted) as primary; vendor SRIM-tabulated SiC as higher-fidelity option. Also: w_SiC = 7.83 eV (Ricossa 2024) vs 8.4 eV currently in code — worth re-validating.

---

## New vs Extended Modules

| #   | Feature                      | Category                      | Module(s)                                                                                      | Key constraint                                                                                                         |
| --- | ---------------------------- | ----------------------------- | ---------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| F1  | Graded epi doping 2D         | **CALIBRATE existing**        | `src/device2d.py`                                                                              | `set_graded_doping_2d` already exists; re-tune constants against 2D C-V data at reverse bias                           |
| F2  | ROOT/Geant4 real integration | **EXTEND + FIXTURE**          | `src/mc_coupling.py`, `tests/fixtures/`, `scripts/generate_geant4_fixture.py`                  | Reader exists mock-only; add real fixture + end-to-end test; do NOT create separate `src/root_io.py`                   |
| F3  | PSTAR+SRIM kappa             | **DATA REPLACE**              | `data/` CSV files, `src/stopping_power.py` (new utility)                                       | Replace placeholder CSVs producing flat κ [0.575, 0.587]; `microdosimetry.py` already loads them                       |
| F4  | Complete noise analysis      | **NEW module**                | `src/noise.py`                                                                                 | `optimization.py` has shot noise only; pure-computation, no devsim dependency; easily testable                         |
| F5  | Build-up over-response 2D    | **NEW module**                | `src/build_up.py`                                                                              | Post-processing correction on 2D CCE output; depends on stable F1 near-surface field                                   |
| F6  | Azimuthal / angular response | **NEW module**                | `src/azimuthal.py` (2D polar-θ); F6B needs 3D for true azimuthal φ on square SV                | 2D polar-θ is tractable; true φ-sweep on 100×100/300×300 µm SV requires 3D                                             |
| F7  | Anisotropic mobility tensor  | **NEW module + opt-in flag**  | `src/mobility_tensor.py`, `src/sic_material.py`, `src/device2d.py` (kwarg `anisotropic=False`) | devsim has NO native tensor mobility; custom edge models via `unitx`/`unity`; default off to preserve all 20 notebooks |
| F8  | Full 3D simulation           | **NEW module (stretch goal)** | `src/device3d.py`, `src/charge_collection_3d.py`                                               | gmsh MSH v2.2 ONLY (devsim rejects v4.1); physically inconsistent without F7; confirm with PI                          |

**6 new source files:** `noise.py`, `build_up.py`, `azimuthal.py`, `mobility_tensor.py`, `device3d.py`, `charge_collection_3d.py`

---

## Recommended Phase Order

```
TIER 1 — Parallel unblockers (no devsim changes; can develop simultaneously)
=============================================================================
Phase 26 | F1 | Graded doping 2D calibration          → unblocks F4, F5, F7, F8
Phase 27 | F3 | PSTAR+SRIM real stopping power data    → independent; fixes flat-κ bug
Phase 28 | F2 | Geant4 ROOT fixture + integration test → independent; closes MCCP-02 gap

TIER 2 — Analysis modules (after Phase 26 ships)
=================================================
Phase 29 | F4 | src/noise.py — complete noise model (needs stable dark current from F1)
Phase 30 | F5 | src/build_up.py — surface dead-zone (needs stable near-surface field from F1)

TIER 3 — New devsim physics (sequential)
=========================================
Phase 31 | F7 | src/mobility_tensor.py — anisotropic mobility (prerequisite for F8)
Phase 32 | F6A | src/azimuthal.py 2D — polar-angle θ sweep on 2D mesh

TIER 4 — Stretch goals (F8 is PI decision)
===========================================
Phase 33 | F8  | src/device3d.py — full 3D simulation (confirm scope with PI first)
Phase 34 | F6B | Azimuthal 3D — true φ-sweep on square SV (depends on F8)
Phase 35 | Audit | v4.0 milestone audit, paper figures, integration regression sweep
```

---

## Top Risks (Watch Out For)

### Risk 1 — Kappa unit mismatch produces wrong tissue-equivalence (HIGH)

PSTAR outputs MeV·cm²/g; SRIM defaults to MeV/(mg/cm²). Mixing without conversion gives κ wrong by SiC density factor (~3×). Compound with existing flat-κ bug (coarse sampling + linear interpolation). **Fix:** `StoppingPowerTable` class with explicit source-unit declarations; log-log interpolation; ≥50 energy points per decade; unit test against published 100 MeV proton value in SiC.

### Risk 2 — devsim global-state contamination (HIGH)

Parameters `raxis_zero`, `node_volume_model`, `edge_couple_model` are process-global. Adding 3D introduces `element_edge_couple_model`. Hardcoded device names collide on notebook re-run. **Fix:** Extend `reset_devsim()` to clear all global parameters; `make_unique_device_name()` factory; regression test running alt-structures then planar notebooks in sequence.

### Risk 3 — 3D mesh format rejection + node-count explosion (HIGH)

gmsh ≥4 defaults to MSH 4.1; devsim rejects it. Uniform 0.5 µm refinement on 300×300×10 µm → >100 M tetrahedra. **Fix:** Always set `Mesh.MshFileVersion = 2.2`; graded mesh (0.2 µm near junction, 5 µm bulk); warn if >1 M nodes; assert .msh first-line version.

### Risk 4 — No native tensor mobility in devsim (HIGH)

Built-in DD uses scalar μ only. Custom edge models required. Labeling error (μ∥ vs μ⊥) silently inverts physics. Note: 4H-SiC electrons have μ∥c > μ⊥c (opposite to holes AND opposite to 6H-SiC). **Fix:** Naming `mu_n_c_axis` / `mu_n_basal_plane`; explicit `c_axis_direction` parameter; validate isotropic limit reproduces v3.0 within 0.1%.

### Risk 5 — Undefined Geant4 ROOT schema (HIGH)

Every group uses different branch names; energy often in keV not MeV. Real INFN-LNS file will differ from synthetic fixture. **Fix:** `RootSchemaMap` dataclass for configurable branch-name mapping; golden fixture (CSV + ROOT for same events); KS test asserting both readers produce equivalent y-spectra.

---

## Key Decisions Needed

| #   | Decision                                              | Recommendation                                                                                                       |
| --- | ----------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| D1  | Real ROOT file from INFN-LNS availability             | Build Phase 28 against synthetic fixture; do not block on external deliverable                                       |
| D2  | Is full 3D (Phase 33/F8) required for the v4.0 paper? | Treat as stretch goal; confirm with PI which result genuinely needs Cartesian 3D                                     |
| D3  | Hooge α policy for noise module                       | Explicit parameter; 3 named presets (`sic_best=2e-5`, `typical=1e-4`, `worst=1e-3`); notebook shows sensitivity band |
| D4  | Regression baseline after F1 + F7 physics changes     | Freeze v3.0 baselines; `anisotropic=False` default preserves all existing notebooks                                  |
| D5  | SiC stopping power: Bragg additivity vs SRIM table    | Implement both behind `source=` switch; default `pstar_bragg` for CI; SRIM for publication                           |

---

_Synthesized from: STACK.md, FEATURES.md, ARCHITECTURE.md, PITFALLS.md_
_Ready for: gsd-roadmapper (phases 26–35)_
