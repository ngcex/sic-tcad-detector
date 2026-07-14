# Phase 43: Integration Audit — All 20 Notebook Workflows - Research

**Researched:** 2026-07-14
**Domain:** Verification / integration audit (no new physics, no new packages)
**Confidence:** HIGH (all findings grep- or source-verified against the working tree)

## Summary

Phase 43 is a **documentation + verification audit**, not a build phase. ROADMAP explicitly
tags it "audit — no new code," and Success Criterion 3 states that any gaps are _logged as v6.0
tech debt, not silently omitted (nor fixed here)_. The single phase requirement is **FEAT-05**
("all 20 notebook workflows have an equivalent UI workflow"). This research catalogs the 20
notebooks, the 8 Streamlit workflow pages, produces a notebook→page coverage matrix, audits the
25 v5.0 requirements against actual source, and identifies the concrete gaps.

**Key tension (surface to planner):** FEAT-05 as written requires _all 20_ notebooks to have UI
equivalents. They demonstrably do not — several v3.0 microdosimetry workflows (single-particle
2D transient, alternative structures, optimization/feasibility, MC-coupling PHD, published-data
validation, multi-defect comparison) and two facades (`run_flash_recombination`, `run_transient`)
have **no page**. Therefore FEAT-05 will be marked **"partially satisfied — residual gaps logged
as v6.0 tech debt."** The phase deliverable is an audit document (recommended: a v5.0 milestone
audit note following the existing `v3.0-MILESTONE-AUDIT.md` template) plus a verification pass —
**not** feature code to close the gaps.

**Second tension (Success Criterion 4):** SC4 requires a final `pytest -q` to pass. STATE.md and
`35-02-SUMMARY.md` prove that **bare single-process `pytest -q` is unsatisfiable on this machine**
due to devsim resource exhaustion (reproduces on the pre-refactor commit too). The established
project convention (PKG-03) is **per-file/per-class test isolation**. See Open Questions Q1.

**Primary recommendation:** Plan this phase as (1) a per-notebook UI-reproduction verification
matrix (executor click-tests each workflow at a converging bias), (2) a source-verified 25-item
requirements re-check that also corrects three stale `[ ]` checkboxes (UI-01/UI-02/UI-07, actually
SATISFIED per Phase 38), (3) a per-file/per-class `pytest` isolation gate (not bare `pytest -q`),
and (4) a written v6.0 tech-debt log for every uncovered workflow. No source edits to close gaps.

## User Constraints

No `CONTEXT.md` exists for Phase 43 (`/gsd:discuss-phase` was not run — the phase directory is
empty). Constraints are therefore taken from ROADMAP.md, REQUIREMENTS.md, and STATE.md:

### Locked Decisions (from ROADMAP / STATE)

- **Audit phase — no new code.** ROADMAP line 112: "verify all 20 notebook workflows are
  reproducible via UI; confirm all 25 v5.0 requirements satisfied _(GROUP C, audit — no new code)_."
- **Gaps are logged, not fixed.** SC3: "any gaps are logged as v6.0 tech debt, not silent omissions."
- **No physics changes in any v5.0 phase.** STATE.md line 49: "No physics changes allowed in any
  v5.0 phase — refactor only." This forbids touching the `ramp_bias` convergence bug here.
- **`device.py` (1D) is frozen** to protect the 20 validated notebooks (STATE decisions).

### Claude's Discretion

- Exact format/location of the audit deliverable (a `43-*-SUMMARY.md`, a
  `.planning/milestones/v5.0-MILESTONE-AUDIT.md`, or a dedicated tech-debt file). Recommendation:
  follow the existing `v3.0-MILESTONE-AUDIT.md` template (see "v6.0 Tech Debt Logging" below).
- The converging bias value used when click-testing CCE/field workflows (`bias_V=-20.0` is
  documented as converging; the default `-100.0` does not — see Pitfall 2).

### Deferred Ideas (OUT OF SCOPE for Phase 43)

- Fixing the `ramp_bias` deep-bias non-convergence (physics/numerics — separate follow-up).
- Building any missing page (alt-structures, optimization, single-particle transient, etc.).
- Fixing `kappa(E)` data-blocked placeholders (v4.0 Phase 27 work).

## Phase Requirements

| ID      | Description                                                                                                                             | Research Support                                                                                                                                                                                                                                                                                                                                                                                                      |
| ------- | --------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FEAT-05 | All 20 notebook workflows (`notebooks/01_*` through `notebooks/20_*`) have an equivalent UI workflow accessible via the Streamlit pages | The Coverage Matrix below maps every notebook to page(s) with FULL/PARTIAL/NONE. FEAT-05 is **partially satisfiable**: core electrical + damage + microdosimetry-upload + batch-sweep workflows are covered; several v3.0 2D/optimization workflows and two facades (`run_flash_recombination`, `run_transient`) have no page. Planner must decide the FEAT-05 disposition (partial-satisfied + v6.0 log) explicitly. |

## Architectural Responsibility Map

| Capability                | Primary Tier                      | Secondary Tier                 | Rationale                                                                        |
| ------------------------- | --------------------------------- | ------------------------------ | -------------------------------------------------------------------------------- |
| Notebook workflow catalog | Documentation / audit             | —                              | Read-only inspection of `notebooks/*.ipynb`                                      |
| UI page inventory         | Streamlit app (`app/workflows/`)  | —                              | Each page is a `render()` callable registered in `app/main.py`                   |
| Simulation execution      | Public API (`petringa.run_*`)     | devsim core (`petringa/core/`) | Pages call facades as module attributes; facades wrap core physics               |
| Requirement verification  | Audit doc + source read + AppTest | Live browser click-through     | SC2 requires checking against the _running_ app, not just checkboxes             |
| Test gate                 | pytest (per-file isolation)       | —                              | devsim single-process exhaustion forbids a bare `pytest -q`                      |
| Gap logging               | `.planning/` audit doc            | —                              | v6.0 tech-debt destination (no TECH-DEBT.md exists; use milestone-audit pattern) |

## Standard Stack

**N/A — audit phase, no new packages.** The existing stack is Streamlit + Plotly + devsim +
numpy/scipy/pandas, all already installed via `pyproject.toml` (hatchling). No `npm`/`pip`/`cargo`
installs occur in this phase.

## Package Legitimacy Audit

**N/A — this phase installs no external packages.** (One _optional, unverified_ avenue —
`pytest-forked` — is discussed in Open Questions Q1 as a possible way to make a single `pytest`
command satisfy SC4 literally; it is **not** a recommendation and would require slopcheck +
registry verification before adoption.)

## The 20 Notebooks (Catalog)

**Notebook count reconciliation [VERIFIED: `ls notebooks/`]:** There are 21 `.ipynb` files. One is
a duplicate execution artifact (`03_executed.ipynb`) — not a distinct workflow. There are **two
`05_` files** (`05_dark_current_vs_fluence.ipynb` and `05_parametric_studies.ipynb`). Excluding the
`03_executed` duplicate leaves exactly **20 canonical workflow notebooks** matching the `01–20`
prefix range. Filenames do **not** align 1:1 with phase numbers (e.g. `06` = "Phase 10", `15`+ =
v3.0 phases 19–25). `[ASSUMED]` that FEAT-05's "20 notebooks" = these 20 (the two `05_` files both
count; `03_executed` does not) — planner should confirm this framing.

| #   | Notebook                         | Workflow                                                       | Key `petringa.core` imports                                                 | Nearest UI facade                                             |
| --- | -------------------------------- | -------------------------------------------------------------- | --------------------------------------------------------------------------- | ------------------------------------------------------------- |
| 01  | `01_phase1_validation`           | 1D material params, analytical electrostatics, Poisson E-field | `analytical`, `poisson`, `device`, `incomplete_ionization`, `sic_material`  | (field map / no direct page for analytical)                   |
| 02  | `02_electrical_characterization` | C-V analysis, graded doping, IV/DD                             | `cv_analysis`, `drift_diffusion`, `validation`                              | `run_cv` (C-V page)                                           |
| 03  | `03_charge_collection`           | CCE vs bias, generation profiles, Hecht                        | `charge_collection`, `generation_profiles`                                  | `run_cce` (CCE page)                                          |
| 04  | `04_flash_recombination`         | FLASH Auger recombination, CCE vs dose rate                    | `flash_recombination`, `generation_profiles`                                | `run_flash_recombination` — **NO page**                       |
| 05a | `05_dark_current_vs_fluence`     | Dark current vs proton fluence (delta-J)                       | `dark_current`, `radiation_damage`                                          | `run_radiation_damage` / `run_dark_current` (partial)         |
| 05b | `05_parametric_studies`          | Multi-dim CCE vs dose-rate sweep (epi×doping×bias)             | `flash_recombination` (`parametric_cce_sweep`), `plotting`                  | Batch Sweep (partial — FLASH sweep not a facade)              |
| 06  | `06_temperature_dependence`      | Material props + device response vs T (280–350K)               | `sic_material`, `temperature_sweep`                                         | `run_temperature_sweep` (via Batch Sweep)                     |
| 07  | `07_dark_current`                | Dark current decomposition (SRH/TAT/SRV) vs T                  | `dark_current`, `sic_material`                                              | `run_dark_current` (Dark Current page)                        |
| 08  | `08_transient_flash`             | Transient I(t) waveform, single FLASH pulse                    | `transient`, `drift_diffusion`, `flash_recombination`                       | `run_transient` — **NO page**                                 |
| 09  | `09_radiation_damage`            | Defect physics (Z1/2/EH6/7/EH4), lifetime degradation          | `radiation_damage`                                                          | `run_radiation_damage` (partial — physics view)               |
| 10  | `10_cce_vs_fluence`              | CCE degradation vs fluence (62 MeV proton)                     | `charge_collection` (`cce_vs_fluence`), `radiation_damage`                  | `run_radiation_damage` (Radiation Damage page)                |
| 11  | `11_dark_current_cv_evolution`   | Dark current + C-V evolution vs fluence                        | `cv_analysis` (`cv_at_fluence`), `dark_current`, `radiation_damage`         | partial (`run_dark_current`; no CV-vs-fluence UI)             |
| 12  | `12_multi_defect_comparison`     | 3-defect vs single-defect model comparison                     | `charge_collection`, `cv_analysis`, `dark_current`, `radiation_damage`      | **NO page** (comparison view)                                 |
| 13  | `13_parametric_optimization`     | Radiation-hardness sweep (epi×doping×bias)                     | `charge_collection`, `radiation_damage` (`radiation_hardness_sweep`)        | Batch Sweep (partial)                                         |
| 14  | `14_validation`                  | Validation vs published irradiation data                       | `charge_collection`, `cv_analysis`, `dark_current`, `validation`            | **NO page** (validation-vs-literature)                        |
| 15  | `15_2d_electrostatics_cce`       | 2D electrostatics + CCE, 2D vs 1D validation                   | `device2d`, `charge_collection_2d`, `plotting2d`, `poisson`                 | `run_field` 2D heatmap (Field Map + Geometry Viewer, partial) |
| 16  | `16_single_particle_cce`         | Single-particle 2D transient, CCE(LET) table                   | `single_particle`, `charge_collection_2d`                                   | **NO page** (single-particle transient)                       |
| 17  | `17_mc_coupling`                 | MC coupling: energy deposition → PHD                           | `mc_coupling`, `single_particle` (`load_cce_let_table`)                     | partial (microdosimetry page does spectrum, not PHD)          |
| 18  | `18_microdosimetric_spectra`     | y·d(y) spectrum, f(y)/d(y), y_F/y_D, kappa                     | `mc_coupling`, `microdosimetry`, `single_particle`                          | `run_microdosimetry` (Microdosimetry page)                    |
| 19  | `19_alternative_structures`      | Mesa / 3D-electrode / delta-E/E structure comparison           | `alternative_structures`, `charge_collection_2d`, `microdosimetry`          | **NO page** (alt structures)                                  |
| 20  | `20_feasibility_report`          | Parametric optimization + noise + fab recommendations          | `optimization`, `charge_collection_2d`, `microdosimetry`, `single_particle` | **NO page** (optimization/feasibility)                        |

## The 8 UI Workflow Pages (Inventory)

[VERIFIED: source read of `app/main.py` + each `app/workflows/*.py`]

| Page              | File                  | Facade(s) called                                                  | Notes                                                                                                                                      |
| ----------------- | --------------------- | ----------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| Home              | `home.py`             | none                                                              | Landing / config summary                                                                                                                   |
| C-V Analysis      | `cv.py`               | `petringa.run_cv`                                                 | Plotly C-V + Mott-Schottky, CSV download. Works end-to-end (browser-confirmed)                                                             |
| Charge Collection | `cce.py`              | `petringa.run_cce`                                                | CCE vs bias + CSV. **Blocked at default config by `ramp_bias` non-convergence** (graceful `st.error`)                                      |
| Field Map         | `field_map.py`        | `petringa.run_field`                                              | 1D depth profiles + 2D geometry heatmap (Geometry Viewer). 2D at shallow bias browser-confirmed; 1D default deep bias hits convergence bug |
| Radiation Damage  | `radiation_damage.py` | `petringa.run_radiation_damage`                                   | CCE vs fluence + persistent kappa-placeholder warning banner + CSV                                                                         |
| Dark Current      | `dark_current.py`     | `ParametricSweep(param="T", sim_fn=run_dark_current)`             | J_SRH/J_TAT/J_SRV decomposition vs T + CSV                                                                                                 |
| Microdosimetry    | `microdosimetry.py`   | `petringa.run_microdosimetry`                                     | CSV upload → tempfile bridge → y·d(y) spectrum + y_F/y_D + CSV                                                                             |
| Batch Sweep       | `batch_sweep.py`      | `ParametricSweep` over `run_cce`/`run_cv`/`run_temperature_sweep` | Curated `SWEEPABLE_FIELDS` (8 numeric fields) × 3 `SIM_FACADES`; overlay + bulk CSV                                                        |

**Public API vs pages [VERIFIED: `petringa/__init__.py` + grep of `app/workflows/`]:** The public
API exports 9 `run_*` facades + `ParametricSweep`. Of these, **`run_flash_recombination` and
`run_transient` are surfaced by NO page** (grep of `app/workflows/` returns zero hits). These are
notebooks 04 (FLASH), 05b (FLASH parametric), and 08 (transient).

## Notebook → Page Coverage Matrix

**Rubric:**

- **FULL** — the notebook's primary observable is reproducible via a UI page for equivalent inputs.
- **PARTIAL** — a related page exists but the specific workflow/plot/comparison is not reproducible,
  or is blocked by the convergence bug at default config.
- **NONE** — no page exposes the underlying facade/workflow at all.

**Confidence:** HIGH for NONE rows (grep-verified: no facade/page). MEDIUM for FULL/PARTIAL
(title-vs-facade inference; executor must click-test to confirm the observable renders).

| #   | Notebook                    | Coverage | Page                            | What's missing (for PARTIAL/NONE)                                                                                                                             |
| --- | --------------------------- | -------- | ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 01  | phase1_validation           | PARTIAL  | Field Map                       | Analytical electrostatics (built-in potential, depletion width, Hecht) not surfaced; only numerical field map                                                 |
| 02  | electrical_characterization | FULL     | C-V                             | C-V + Mott-Schottky + CSV; browser-confirmed working                                                                                                          |
| 03  | charge_collection           | PARTIAL  | CCE                             | CCE vs bias exists but **blocked at default deep bias** by `ramp_bias`; reproducible only at shallow bias. Generation-profile / Hecht comparison plots absent |
| 04  | flash_recombination         | NONE     | —                               | `run_flash_recombination` has no page (FLASH Auger CCE vs dose rate)                                                                                          |
| 05a | dark_current_vs_fluence     | PARTIAL  | Radiation Damage / Dark Current | No single page does dark-current-_vs-fluence_; damage page does CCE-vs-fluence, dark page does J-vs-T                                                         |
| 05b | parametric_studies          | PARTIAL  | Batch Sweep                     | FLASH multi-dim (`parametric_cce_sweep`) is not one of the 3 `SIM_FACADES`; only CCE/CV/T sweeps                                                              |
| 06  | temperature_dependence      | PARTIAL  | Batch Sweep                     | `run_temperature_sweep` available via Batch Sweep; material-props-vs-T plots (E_g, n_i, mobility) not surfaced                                                |
| 07  | dark_current                | FULL     | Dark Current                    | J_SRH/J_TAT/J_SRV decomposition vs T + CSV — matches notebook 07                                                                                              |
| 08  | transient_flash             | NONE     | —                               | `run_transient` has no page (I(t) waveform, single-pulse dynamics)                                                                                            |
| 09  | radiation_damage            | PARTIAL  | Radiation Damage                | Page runs the damage _sweep_; the defect-physics inspection plots (lifetime degradation curves) are notebook-only                                             |
| 10  | cce_vs_fluence              | FULL     | Radiation Damage                | CCE vs fluence for specified fluence range + energy + CSV                                                                                                     |
| 11  | dark_current_cv_evolution   | PARTIAL  | Dark Current                    | Dark-current side partially covered; **C-V-vs-fluence evolution has no UI**                                                                                   |
| 12  | multi_defect_comparison     | NONE     | —                               | No page compares 3-defect vs single-defect models                                                                                                             |
| 13  | parametric_optimization     | PARTIAL  | Batch Sweep                     | `radiation_hardness_sweep` isn't a `SIM_FACADE`; single-param sweeps only, no optimization scoring                                                            |
| 14  | validation                  | NONE     | —                               | No page validates predictions against published/literature data                                                                                               |
| 15  | 2d_electrostatics_cce       | PARTIAL  | Field Map + Geometry Viewer     | 2D field/doping heatmap works (browser-confirmed at shallow bias); 2D **CCE** + 2D-vs-1D validation not surfaced                                              |
| 16  | single_particle_cce         | NONE     | —                               | No page for single-particle 2D transient / CCE(LET) table generation                                                                                          |
| 17  | mc_coupling                 | PARTIAL  | Microdosimetry                  | Microdosimetry page consumes an MC CSV → spectrum, but the PHD (pulse-height-distribution) intermediate view and CCE(LET) coupling aren't separately shown    |
| 18  | microdosimetric_spectra     | FULL     | Microdosimetry                  | CSV upload → y·d(y) spectrum + y_F/y_D + CSV download — the design-spec acceptance gate for this page                                                         |
| 19  | alternative_structures      | NONE     | —                               | No page for mesa / 3D-electrode / delta-E/E structure comparison                                                                                              |
| 20  | feasibility_report          | NONE     | —                               | No page for parametric optimization + noise + fab recommendations                                                                                             |

**Coverage tally (MEDIUM confidence, executor to confirm):** FULL ≈ 4 (02, 07, 10, 18);
PARTIAL ≈ 9; NONE ≈ 7 (04, 08, 12, 14, 16, 19, 20). This is the concrete evidence that FEAT-05
is **partially satisfied** — the audit's job is to record each PARTIAL/NONE cell as a v6.0
tech-debt item, not to close it.

## 25 v5.0 Requirements — Source-Verified Status

SC2 requires checking each requirement **against the running app/source**, not the checkbox. Three
checkboxes are **stale** (marked `[ ]` in REQUIREMENTS.md but actually SATISFIED per Phase 38's
browser-confirmed verification). The audit should correct these.

| Req         | REQUIREMENTS.md box | Actual (source-verified)       | Evidence                                                                                                     |
| ----------- | ------------------- | ------------------------------ | ------------------------------------------------------------------------------------------------------------ |
| PKG-01      | `[x]`               | SATISFIED                      | `pyproject.toml` hatchling; editable install works                                                           |
| PKG-02      | `[x]`               | SATISFIED                      | runtime + `[dev]` deps declared                                                                              |
| PKG-03      | `[x]`               | SATISFIED (per-file isolation) | see Open Questions Q1                                                                                        |
| LIB-01..07  | `[x]`×7             | SATISFIED                      | `petringa/__init__.py` exports all 9 facades + `DeviceConfig`/`SimResult`/`MeshData`/`ParametricSweep`       |
| **UI-01**   | **`[ ]`**           | **SATISFIED (stale box)**      | `38-VERIFICATION.md`: `st.navigation` w/ 8 pages boots cleanly, browser-confirmed                            |
| **UI-02**   | **`[ ]`**           | **SATISFIED (stale box)**      | `38-VERIFICATION.md`: all 11 `DeviceConfig` fields as sidebar controls, persistence fixed                    |
| UI-03       | `[x]`               | SATISFIED                      | C-V page browser-confirmed (chart + CSV)                                                                     |
| UI-04       | `[~]`               | PARTIAL (upstream blocker)     | `run_cce` non-convergence at default bias; graceful `st.error`; works at shallow bias                        |
| UI-05       | `[~]`               | PARTIAL (upstream blocker)     | `run_field` non-convergence at default `bias_V=-100`; graceful `st.error`; 2D shallow-bias browser-confirmed |
| UI-06       | `[~]`               | PARTIAL                        | CSV download browser-confirmed on C-V; CCE/field blocked by UI-04/05                                         |
| **UI-07**   | **`[ ]`**           | **SATISFIED (stale box)**      | `38-VERIFICATION.md`: cross-page persistence, Playwright-confirmed (commit `f1bc8b7`)                        |
| VIZ-01      | `[ ]`               | SATISFIED (stale box)          | `40-VERIFICATION.md`: real 2D solve → heatmap, browser-confirmed 2026-07-13                                  |
| VIZ-02      | `[ ]`               | SATISFIED (stale box)          | `40-VERIFICATION.md`: 1D → depth-profile bar, browser-confirmed                                              |
| VIZ-03      | `[ ]`               | SATISFIED (stale box)          | `40-VERIFICATION.md`: quantity dropdown, no re-solve (call-counter test + browser)                           |
| FEAT-01..04 | `[x]`×4             | SATISFIED                      | Radiation Damage, Dark Current, Microdosimetry, Batch Sweep pages all built + verified                       |
| FEAT-05     | `[ ]`               | PARTIAL — this phase           | Coverage Matrix above; ~7 NONE + ~9 PARTIAL cells → v6.0 log                                                 |

**Note:** VIZ-01/02/03 are also `[ ]` in REQUIREMENTS.md's checkbox list but the traceability table
(line 207–209) says "Pending" while `40-VERIFICATION.md` reports all three SATISFIED and
browser-confirmed. These are stale as well — the audit should reconcile REQUIREMENTS.md against the
Phase 40 verification. So **6 requirement boxes are stale** (UI-01, UI-02, UI-07, VIZ-01, VIZ-02,
VIZ-03), all actually SATISFIED. The only genuinely-open items are UI-04/05/06 (upstream solver)
and FEAT-05 (this phase).

## Known Unresolved Issues from Phases 35–42

| Issue                                 | Source                                  | Impact on Phase 43                                                                                                                                                                                                                                                                                                                  |
| ------------------------------------- | --------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ramp_bias` deep-bias non-convergence | `39-VERIFICATION.md`, STATE.md line 130 | `run_cce` fails ~V≈60.5V, `run_field` fails ~V=66V for default `DeviceConfig()` (default `bias_V=-100`). Shallow bias (`-20V`) converges. Blocks UI-04/05/06 and makes notebooks 03/15/16 workflows only reproducible at shallow bias. **Out of scope to fix** (physics/numerics, v5.0 forbids physics changes) → log as v6.0 debt. |
| Distinct from Phase 26 fix            | `39-VERIFICATION.md` line 118           | Phase 26 fixed _2D uniform-doping_ divergence via graded epi. This is a _1D graded-doping deep-bias_ failure — a separate, still-open issue.                                                                                                                                                                                        |
| `kappa(E)` data-blocked               | STATE.md line 128, v4.0 audit C-5       | Radiation Damage page shows a warning banner; absolute Φ_crit unvalidated. Not a v5.0 gap (v4.0 Phase 27 work). Note in audit for completeness.                                                                                                                                                                                     |
| Bare `pytest -q` unsatisfiable        | `35-02-SUMMARY.md`, STATE.md line 129   | devsim single-process resource exhaustion (reproduces pre-refactor). PKG-03 gate = per-file/per-class isolation. Directly affects SC4 — see Open Questions Q1.                                                                                                                                                                      |

## Runtime State Inventory

**N/A — Phase 43 is a read-only audit.** It writes no stored data, registers no OS state, changes
no secrets/env vars, and produces no build artifacts. The only deliverable is Markdown audit
documentation under `.planning/`. Verified: the phase requirement (FEAT-05) and all success
criteria are verification/documentation actions, not migrations. **Nothing to migrate in any of
the 5 categories — verified by reading the phase goal and success criteria.**

## v6.0 Tech-Debt Logging (SC3 destination)

[VERIFIED: `find .planning -iname "*tech*debt*"` + `.planning/v3.0-MILESTONE-AUDIT.md`]

**There is no `TECH-DEBT.md` or `backlog.md`** in `.planning/`. The established project convention
for logging accumulated tech debt at a milestone boundary is a **milestone-audit document** with a
`tech_debt:` YAML frontmatter section, per `.planning/v3.0-MILESTONE-AUDIT.md` (which has
`status: tech_debt`, a "## Tech Debt Summary" section, and per-phase item lists totaling "10 items
across 5 phases"). `.planning/milestones/` exists as a directory for archived roadmaps.

**Recommendation for the planner:** the SC3 gap log should live in a
`.planning/milestones/v5.0-MILESTONE-AUDIT.md` (or `43-*-SUMMARY.md`) that mirrors the v3.0 audit
structure: a Requirements-Coverage cross-reference, the notebook→page coverage matrix, and a
"Tech Debt Summary" listing each NONE/PARTIAL cell as a discrete v6.0 backlog item (e.g. "v6.0:
add single-particle transient page (notebook 16)", "v6.0: fix `ramp_bias` deep-bias convergence
blocking UI-04/05"). This satisfies SC3's "logged as tech debt, not silent omission."

## Common Pitfalls

### Pitfall 1: Trusting the REQUIREMENTS.md checkboxes over source

**What goes wrong:** UI-01/02/07 and VIZ-01/02/03 are marked `[ ]`/"Pending" in REQUIREMENTS.md but
are actually SATISFIED (browser-confirmed in Phases 38 and 40).
**How to avoid:** SC2 mandates checking the _running app/source_. The audit must read
`38-VERIFICATION.md` and `40-VERIFICATION.md` and reconcile REQUIREMENTS.md, not copy its boxes.
**Warning sign:** an audit that reports "6 requirements pending" when 6 boxes are merely stale.

### Pitfall 2: Click-testing CCE/field workflows at the default bias

**What goes wrong:** `run_cce` and `run_field` raise devsim `RuntimeError` ("ramp_bias: failed to
converge") for the plain default `DeviceConfig()` (default `bias_V=-100`). An executor click-testing
notebook 03/15/16 at defaults will see the graceful `st.error` and wrongly log "page broken."
**How to avoid:** click-test these workflows at a **converging shallow bias** (`bias_V=-20.0` is
documented as converging in `40-VERIFICATION.md`). The page correctly renders at shallow bias.
**Warning sign:** conflating the upstream solver bug (UI-04/05, already known) with a UI defect.

### Pitfall 3: Trying to fix gaps in an audit phase

**What goes wrong:** the coverage matrix reveals 7 NONE cells; the temptation is to build the missing
pages or fix `ramp_bias`. ROADMAP ("no new code"), SC3 ("logged, not fixed"), and STATE ("no physics
changes in v5.0") all forbid this.
**How to avoid:** the deliverable is documentation. Every gap → a v6.0 tech-debt line item.

### Pitfall 4: Running a bare `pytest -q` and reporting failure

**What goes wrong:** SC4 literally says `pytest -q`; running it single-process crashes from devsim
resource exhaustion — but this is a _known machine limitation_, not a regression.
**How to avoid:** use the PKG-03 convention — per-file/per-class isolation. See Open Questions Q1.

## Validation Architecture

`workflow.nyquist_validation` is **absent** from `config.json` → treated as **enabled**.

### Test Framework

| Property           | Value                                                                        |
| ------------------ | ---------------------------------------------------------------------------- |
| Framework          | pytest (declared in `[dev]` extras; run via `uv run pytest`)                 |
| Config file        | none dedicated; `pyproject.toml` + `uv sync --extra dev` (STATE.md line 110) |
| Quick run command  | `uv run pytest tests/test_api_*.py -q` (fast API unit tests, no heavy DD)    |
| Full suite command | **per-file/per-class isolation** (NOT bare `pytest -q`) — see Q1             |

### Phase Requirements → Test Map

| Req       | Behavior                         | Test type              | Automated command                                                                                                                                                                                 | Exists?                          |
| --------- | -------------------------------- | ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------- |
| FEAT-05   | 20 notebooks reproducible via UI | manual (click-through) | human browser verification per workflow                                                                                                                                                           | ❌ manual-only by nature         |
| SC4       | all `tests/test_api_*.py` pass   | unit                   | `uv run pytest tests/test_api_cv.py tests/test_api_cce.py tests/test_api_device.py tests/test_api_facades.py tests/test_api_field.py tests/test_api_microdosimetry.py tests/test_api_sweep.py -q` | ✅ all 7 exist                   |
| App pages | pages render / mockable          | AppTest                | `uv run pytest tests/test_app_*.py -q`                                                                                                                                                            | ✅ 13 `test_app_*` modules exist |

### Sampling Rate

- **Per task:** `uv run pytest tests/test_api_*.py -q` (fast, no devsim DD).
- **Phase gate:** per-file/per-class isolation across all test modules (the DD-heavy `test_cv`,
  `test_charge_collection*`, `test_device2d`, `test_single_particle`, `test_transient`, etc. must be
  run in isolation, not in one process). Plus AppTest suite for the app pages.

### Wave 0 Gaps

- None — all 7 `tests/test_api_*.py` and 13 `tests/test_app_*.py` modules already exist. This phase
  adds **no new tests**; it runs the existing suite as an acceptance gate.

## Security Domain

**N/A — audit phase, no new code, no new attack surface.** The one existing input surface (the
Microdosimetry CSV upload, `st.file_uploader(type=["csv"])` → server-side tempfile) was introduced
in Phase 42, is restricted to CSV, and is out of scope for a read-only audit. No V2–V6 ASVS
controls change. (`security_enforcement` not set in config; noting N/A is the honest disposition.)

## Environment Availability

**N/A for external installs** — but the audit does depend on the existing local toolchain being
runnable: `uv` + the editable `petringa` install + devsim (C extension, local-only). All are
already present (every prior v5.0 phase ran against them). No fallback needed; if `uv run pytest`
silently falls back to a non-project `pytest`, sync with `uv sync --extra dev` first (STATE.md
line 110).

## Assumptions Log

| #   | Claim                                                                                        | Section                | Risk if Wrong                                                                                                 |
| --- | -------------------------------------------------------------------------------------------- | ---------------------- | ------------------------------------------------------------------------------------------------------------- |
| A1  | FEAT-05's "20 notebooks" = the 20 canonical files (both `05_*`; excluding `03_executed`)     | Notebook Catalog       | If the intended set differs, the matrix rows shift; low risk (any reasonable reading gives ~20)               |
| A2  | FULL/PARTIAL classifications (non-NONE rows) hold when click-tested in the browser           | Coverage Matrix        | MEDIUM — executor must confirm each observable renders; NONE rows are grep-hard, FULL/PARTIAL are inference   |
| A3  | The de-facto meaning of SC4's `pytest -q` is per-file/per-class isolation (PKG-03 precedent) | Open Questions Q1      | If the user insists on a literal single `pytest -q`, it is unsatisfiable without `pytest-forked` (unverified) |
| A4  | v6.0 gaps should be logged in a milestone-audit doc (v3.0 template), no TECH-DEBT.md exists  | v6.0 Tech-Debt Logging | Low — planner may choose a different doc; the _content_ is what SC3 requires                                  |

## Open Questions

1. **How is SC4's `pytest -q` satisfied given devsim resource exhaustion?**
   - **What we know:** Bare single-process `pytest -q` is _proven unsatisfiable_ on this machine
     (reproduces on pre-refactor commit `fe3b43c`) — `35-02-SUMMARY.md`, STATE.md line 129. PKG-03
     was formally redefined to accept **per-file/per-class isolation** (all 25 modules pass
     individually).
   - **What's unclear:** SC4 is worded literally ("a final `pytest -q` run … passes all tests").
     The audit cannot unilaterally reinterpret a success criterion.
   - **Recommendation:** Plan the gate as per-file/per-class isolation (the PKG-03 convention) and
     **flag this reinterpretation for user/planner confirmation** in the audit. One _unverified_
     avenue to satisfy SC4 more literally is `pytest-forked` (`--forked` subprocess-isolates each
     test, so a single command could pass) — **not installed, not verified, would need slopcheck +
     registry check before use.** Do not test or adopt it in research; it is a planner option only.

2. **Does the audit correct REQUIREMENTS.md's stale checkboxes, or only report them?**
   - **What we know:** 6 boxes (UI-01/02/07, VIZ-01/02/03) are stale-but-satisfied.
   - **Recommendation:** The audit should (a) document actual status against source, and (b) update
     REQUIREMENTS.md's traceability + checkboxes to match verified reality — but this is a doc edit,
     not "new code," so it is in scope. Confirm with planner whether editing REQUIREMENTS.md is
     desired within this phase.

3. **What is the FEAT-05 disposition — "partial + v6.0 log," or is FEAT-05 held open?**
   - **What we know:** ~7 notebooks have NONE coverage; SC3 says log gaps, don't fix.
   - **Recommendation:** Mark FEAT-05 "partially satisfied — residual gaps logged as v6.0 tech
     debt," consistent with SC3. This is the crux decision the planner must make explicit.

## Sources

### Primary (HIGH confidence)

- `.planning/REQUIREMENTS.md` — full v5.0 section + traceability tables (read verbatim)
- `.planning/ROADMAP.md` — Phase 35–43 descriptions + success criteria (read verbatim)
- `.planning/STATE.md` — accumulated decisions, blockers (ramp_bias, pytest exhaustion)
- `app/main.py`, `app/workflows/*.py`, `app/components/*.py`, `petringa/__init__.py` — source read
- `notebooks/01_*` … `notebooks/20_*` — programmatic extraction of titles + `petringa.core` imports
- `.planning/phases/38-*/38-VERIFICATION.md` — UI-01/02/07 actual SATISFIED status
- `.planning/phases/39-*/39-VERIFICATION.md` — ramp_bias convergence blocker, UI-04/05/06
- `.planning/phases/40-*/40-VERIFICATION.md` — VIZ-01/02/03 SATISFIED (browser-confirmed)
- `.planning/phases/35-*/35-02-SUMMARY.md` — pytest single-process unsatisfiability + PKG-03 gate
- `.planning/v3.0-MILESTONE-AUDIT.md` — tech-debt logging template (SC3 destination)

### Secondary (MEDIUM confidence)

- `docs/superpowers/specs/2026-06-26-simulator-library-ui-design.md` — Phase C acceptance gate
  ("all 20 notebooks replaceable by equivalent UI workflows")

### Tertiary (LOW confidence)

- None — no WebSearch used (this is a codebase-internal audit).

## Metadata

**Confidence breakdown:**

- Notebook catalog: HIGH — programmatic extraction from all 20 files.
- Page inventory: HIGH — direct source read of all 8 workflow modules + main.py.
- Coverage matrix (NONE cells): HIGH — grep-verified missing facades/pages.
- Coverage matrix (FULL/PARTIAL cells): MEDIUM — inference; executor click-test confirms.
- Requirement status: HIGH — cross-checked against phase VERIFICATION docs, not just checkboxes.
- pytest gate: HIGH — established constraint in STATE + 35-02-SUMMARY.

**Research date:** 2026-07-14
**Valid until:** 2026-08-14 (stable — internal codebase audit; only shifts if pages are added
or the ramp_bias bug is fixed in a future phase)
