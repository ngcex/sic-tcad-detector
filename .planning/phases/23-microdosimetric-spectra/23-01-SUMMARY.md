---
phase: 23-microdosimetric-spectra
plan: 01
subsystem: simulation
tags:
  [
    microdosimetry,
    lineal-energy,
    icru-36,
    tissue-equivalence,
    kappa,
    stopping-power,
    f-y,
    d-y,
  ]

requires:
  - phase: 22-monte-carlo-coupling
    provides: "process_mc_ensemble output (event_collected_keV, event_energies_keV) for y-spectra computation"
provides:
  - "Lineal energy spectrum computation: f(y), d(y), y_F, y_D per ICRU Report 36"
  - "Mean chord length calculation for rectangular parallelepiped and slab geometries"
  - "Energy-dependent tissue-equivalence correction via kappa = S_water / S_SiC"
  - "PSTAR water and SRIM SiC stopping power lookup tables"
  - "Publication-quality y*d(y) and y*f(y) spectrum plotting functions"
affects: [23-02, phase-24, notebooks, microdosimetric-spectra]

tech-stack:
  added: []
  patterns:
    [
      "300 log-spaced bins (50/decade, 6 decades) per ICRU 36 convention",
      "f(y) normalization: counts / (N_total * dy) so integral f(y)*dy = 1",
      "d(y) = y * f(y) / y_F with Jensen inequality validation y_D >= y_F",
      "Stopping power CSV lookup tables in data/ directory for kappa interpolation",
    ]

key-files:
  created:
    - "src/microdosimetry.py"
    - "tests/test_microdosimetry.py"
    - "data/stopping_power_water.csv"
    - "data/stopping_power_sic.csv"
  modified: []

key-decisions:
  - "PSTAR water stopping powers bundled as CSV with 37 energy points spanning 0.1-1000 MeV"
  - "SiC stopping powers from Bethe-Bloch scaling of SRIM data, ~1.7x water values"
  - "Constant kappa fallback (0.58) with warning when no energy-dependent table provided"

patterns-established:
  - "Stopping power CSV format: energy_MeV, stopping_power_MeV_cm2_per_g columns"
  - "Kappa table dict format: energy_MeV and kappa arrays for np.interp lookup"
  - "Spectrum result dict: bin_edges, bin_centers, bin_widths, f_y, d_y, y_F, y_D, n_events, y_values"

requirements-completed: [MDOS-01, MDOS-02, MDOS-03, MDOS-04, MDOS-05]

duration: 3min
completed: 2026-03-31
---

# Plan 23-01: Microdosimetry Module Summary

**ICRU Report 36 lineal energy spectra with f(y)/d(y) normalization, energy-dependent kappa tissue-equivalence correction from PSTAR/SRIM stopping powers, and publication-quality y\*d(y) plotting -- all 25 tests pass**

## Performance

- **Duration:** ~3 min
- **Started:** 2026-03-31T23:00:19Z
- **Completed:** 2026-03-31T23:03:08Z
- **Tasks:** 2
- **Files created:** 4

## Accomplishments

- 8 public functions: mean_chord_length, make_y_bins, lineal_energy_spectrum, compute_microdosimetric_means, compute_kappa_table, tissue_equivalence_correction, plot_yd_spectrum, plot_yf_spectrum
- f(y) normalization validated within 5% of 1.0 for 2000 lognormal events; y_D >= y_F (Jensen's inequality) confirmed for all test distributions
- Energy-dependent kappa tissue-equivalence correction with per-event interpolation; kappa values in [0.3, 1.0] range as expected
- All 25 tests pass in 0.56 seconds

## Task Commits

1. **Task 1: Create microdosimetry.py module** - `5e3f338` (feat)
2. **Task 2: Create test_microdosimetry.py** - `8b8a99b` (test)

## Files Created/Modified

- `src/microdosimetry.py` - Full microdosimetric computation pipeline: chord length, y-binning, f(y)/d(y), kappa correction, plotting (480 lines)
- `tests/test_microdosimetry.py` - 25 tests covering all functions: geometry, binning, normalization, Jensen inequality, kappa table, tissue correction, plotting (299 lines)
- `data/stopping_power_water.csv` - PSTAR proton stopping powers in water, 37 points from 0.1-1000 MeV
- `data/stopping_power_sic.csv` - Proton stopping powers in SiC (Bethe-Bloch scaling from SRIM), 37 points matching water grid

## Decisions Made

- Bundled PSTAR water stopping powers as CSV (37 energy points, 0.1-1000 MeV) rather than computing at runtime -- lookup tables are standard practice for stopping power data
- SiC stopping powers scaled from SRIM data at approximately 1.7x water values due to higher Z, consistent with published compound material scaling
- Constant kappa fallback (0.58) provided with explicit logging warning about the approximation, as recommended in research document

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- microdosimetry.py ready for notebook demonstration (Plan 23-02)
- All 8 functions importable and tested
- Stopping power tables loaded and producing kappa in physically reasonable range
- Plot functions verified with Agg backend for non-interactive use

---

_Phase: 23-microdosimetric-spectra_
_Completed: 2026-03-31_
