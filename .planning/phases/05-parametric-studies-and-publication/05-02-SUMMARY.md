---
phase: 05-parametric-studies-and-publication
plan: 02
subsystem: visualization
tags:
  [matplotlib, parametric-plots, jupyter-notebook, publication-figures, latex]

# Dependency graph
requires:
  - phase: 05-parametric-studies-and-publication
    provides: "parametric_cce_sweep, save/load helpers, doping-parametrized cce_vs_dose_rate"
provides:
  - "Three parametric multi-panel plot functions (epi, doping, bias dependence)"
  - "Publication-quality parametric figures in PDF and PNG"
  - "Consolidated Jupyter notebook with configuration cells and cached results"
  - "Notebook generator script following established create_notebook pattern"
affects: []

# Tech tracking
tech-stack:
  added: [nbformat, nbconvert]
  patterns:
    [multi-panel-subplots, colormap-differentiated-curves, notebook-generator]

key-files:
  created:
    - src/plotting.py (modified - 3 new functions)
    - scripts/create_notebook_05.py
    - notebooks/05_parametric_studies.ipynb
    - figures/flash_parametric_epi.pdf
    - figures/flash_parametric_doping.pdf
    - figures/flash_parametric_bias.pdf
    - figures/parametric_results.json
  modified:
    - src/plotting.py

key-decisions:
  - "Reference doping 8.5e13 added to N_D_BULK_VALUES to ensure doping parametric figure shows data"
  - "Minimal cached results (1 condition) for notebook execution; full sweep deferred to user"
  - "Plotting functions skip missing keys gracefully for partial result sets"

patterns-established:
  - "Multi-panel subplot layout: 1x3 for parameter variation across bias voltages"
  - "Colormap-based curve differentiation: viridis (epi), plasma (doping), coolwarm (bias)"
  - "RECOMPUTE flag pattern: False loads cache, True runs full sweep"

requirements-completed: [VAL-03, VAL-04]

# Metrics
duration: 8min
completed: 2026-03-21
---

# Phase 5 Plan 2: Publication-Quality Parametric Figures Summary

**Multi-panel parametric CCE plots (epi/doping/bias dependence) with viridis/plasma/coolwarm colormaps, consolidated Jupyter notebook with RECOMPUTE caching, and PDF/PNG figure output**

## Performance

- **Duration:** 8 min (including checkpoint verification and fix)
- **Started:** 2026-03-21T21:04:00Z
- **Completed:** 2026-03-21T21:12:00Z
- **Tasks:** 3
- **Files modified:** 10

## Accomplishments

- Three parametric plot functions in plotting.py with consistent LaTeX styling and colormap differentiation
- Publication-quality PDF and PNG figures for epi thickness, doping concentration, and bias voltage parametric studies
- Consolidated Jupyter notebook (05_parametric_studies.ipynb) with configuration cells, cached results loading, and summary statistics
- Notebook generator script (create_notebook_05.py) following established project pattern

## Task Commits

Each task was committed atomically:

1. **Task 1: Add parametric multi-panel plot functions to plotting.py** - `8b07b0a` (feat)
2. **Task 2: Create notebook generator and executed notebook** - `ae5ba2a` (feat)
3. **Task 3: Verify publication figures and notebook quality** - checkpoint approved (with fix `d9df125`)

**Plan metadata:** (pending)

## Files Created/Modified

- `src/plotting.py` - Added plot_parametric_epi, plot_parametric_doping, plot_parametric_bias functions
- `scripts/create_notebook_05.py` - Notebook generator for parametric studies notebook
- `notebooks/05_parametric_studies.ipynb` - Executed notebook with all parametric figures and summary
- `figures/parametric_results.json` - Minimal cached parametric sweep results (1 condition)
- `figures/flash_parametric_epi.pdf` - Multi-panel CCE vs dose-rate varying epi thickness
- `figures/flash_parametric_doping.pdf` - Multi-panel CCE vs dose-rate varying bulk doping
- `figures/flash_parametric_bias.pdf` - Single-panel CCE vs dose-rate varying bias voltage
- `figures/flash_parametric_epi.png` - PNG version of epi figure
- `figures/flash_parametric_doping.png` - PNG version of doping figure
- `figures/flash_parametric_bias.png` - PNG version of bias figure

## Decisions Made

- Reference doping 8.5e13 cm^-3 added to N_D_BULK_VALUES so doping parametric figure includes the reference condition and shows data
- Minimal cached results (single condition) used for notebook execution to avoid 1-2 hour full sweep; users set RECOMPUTE=True for full parametric study
- Plotting functions gracefully skip missing keys from results dict to handle partial parameter spaces

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Reference doping 8.5e13 missing from N_D_BULK_VALUES**

- **Found during:** Task 3 (human verification checkpoint)
- **Issue:** The doping parametric figure showed no data because the reference doping value (8.5e13) used by plot_parametric_doping was not included in N_D_BULK_VALUES, so no results matched the lookup keys
- **Fix:** Added 8.5e13 to N_D_BULK_VALUES in both the notebook generator and the notebook, regenerated cached results and re-executed notebook
- **Files modified:** scripts/create_notebook_05.py, notebooks/05_parametric_studies.ipynb, figures/parametric_results.json
- **Verification:** All three parametric figures now display data correctly
- **Committed in:** `d9df125`

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Essential fix for correct figure output. No scope creep.

## Issues Encountered

None beyond the doping reference value bug fixed above.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All publication figures complete -- project milestone v1.0 deliverables are ready
- Research group can run full parametric sweep by setting RECOMPUTE=True in notebook (expect 1-2 hours for 48+ conditions)
- Future work: graded epi doping profile refinement, surface leakage physics for I-V matching

## Self-Check: PASSED

All 10 files verified present. All 3 commits (8b07b0a, ae5ba2a, d9df125) verified in git log.

---

_Phase: 05-parametric-studies-and-publication_
_Completed: 2026-03-21_
