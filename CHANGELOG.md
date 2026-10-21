# Changelog

All notable changes to ETNA are documented in this file.

The format is based on [Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html) for the
version in `etna/_version.py`. Each release has a section here. The release workflow
(`.github/workflows/release.yml`) copies the section matching the tag into the GitHub
release notes. See `docs/RELEASING.md`.

ETNA is distributed from GitHub only (no PyPI release). Issue numbers refer to
<https://github.com/ngcex/sic-tcad-detector/issues>.

## [Unreleased]

## [5.0.0] - 2026-10-21

First versioned release. `5.0.0` has been the package version since the v5.0 library
milestone (2026-07-15), but no tag or GitHub release existed before this one. This
section therefore covers everything up to the first tag. The internal milestones v1.0 to
v6.0, listed under "Development history" at the end of `CHANGELOG.md`, were never tagged
or released.

### Added

- **Packaging and distribution (#3 P0).**
  - The wheel now ships the runtime data (NIST PSTAR/ASTAR/ESTAR tables, SR-NIEL CSVs,
    validation points, foundry layout defaults, measurement templates) as `etna/_data/`.
  - All data lookups go through `etna._paths.data_path()`. The lookup order is
    `$ETNA_DATA_DIR`, then the packaged copy, then the repository `data/`.
  - Raw provenance files (scraped HTML, fetch scripts, synthetic MC events) ship in the
    sdist only.
  - New `wheel-smoke` CI job installs the built wheel outside the checkout and runs
    `scripts/wheel_smoke.py`.
  - `tests/test_packaging.py` guards the force-include table.
- **Project metadata (#3 P1, P3, P4).**
  - Release and citation: `CITATION.cff`, this changelog, `docs/RELEASING.md`, and a
    tag-driven release workflow (sdist + wheel, wheel smoke test, GitHub release, no
    PyPI).
  - CI and dependencies: an Ubuntu + macOS-arm64 smoke matrix, a Dependabot
    configuration, and a dependency licence inventory (`scripts/license_inventory.py`,
    `THIRD_PARTY_LICENSES.md`).
  - Community: contribution, conduct, security and support files, plus issue and PR
    templates.
  - Documentation: an MkDocs site configuration (`mkdocs.yml`) and a manual-only GitHub
    Pages workflow.
- **General 2D engine (#4).** `fv2d` is a numba drift-diffusion engine for laterally
  non-uniform and axisymmetric 2D devices. Contact-cell graded meshes resolve the cathode
  contact cell of the contact-on-epi alternative structures.
- **Radiation-damage kinetics (#5).**
  - Opt-in time-dependent (kinetic) deep-trap occupancy and frequency-aware C-V.
  - Two-level (negative-U) traps and direct donor removal.
  - Calibration of the removal rate on published C-V, plus literature presets (negative-U
    Z1/Z2, EH6/7 emission cross-section).
  - Mechanism notes in `lab/radiation_traps/`.
  - Opt-in charged radiation traps in Poisson beyond the compensation fluence.
- **Physical-model library (#6).** Opt-in 4H-SiC models: high-field (velocity-saturation)
  and anisotropic mobility, impact ionization with avalanche breakdown, and a c-axis
  low-field mobility for vertical transport.
- **Monte Carlo coupling (#7).**
  - A condensed-history SiC deposition sampler with straggling and delta rays.
  - Allpix Squared exchange and event-by-event pulse heights.
  - Depth-resolved hit export.
- **Neutron NIEL and measurement intake (#8).**
  - SiC neutron NIEL from NJOY HEATR (ENDF/B-VIII.0), with opt-in NIEL conventions.
  - `etna.measurements`: editable C-V, SIMS/SRP and I(T) inputs, with unit and
    plausibility checks, model comparison and C-V inversion.
  - A Measurements app page, synthetic round-trip examples and
    `docs/DATA_REQUEST_CHECKLIST.md`.
- **Single-ion track physics (#9).**
  - An axisymmetric single-ion track engine (numba, coupled Newton, TR-BDF2) and
    pair-balance CCE.
  - Opt-in electrothermal coupling with an independent heat-sink ledger.
  - Excess-density unknowns, transverse slices, a 2.5D tilted-track model, tau
    calibration and track facades.
- **Stopping power and kappa.**
  - Real NIST PSTAR/ASTAR tables (Si, C, water, ICRU muscle) replace the placeholder
    path.
  - kappa(E) is now the SiC -> water LET ratio (C-1).
  - The NIST ASTAR Bragg curve is the default alpha profile (Mj-5), with PSTAR energy
    dependence for protons (B-10).
  - MC events carry the primary kinetic energy (C-2).
- **Proton NIEL.** SR-NIEL proton hardness factors for Si and SiC (1-1000 MeV) replace the
  C-5 placeholders (B-07). Out-of-table energies raise.
- **Opt-in incomplete ionization** (A-02), carrier-dependent, in devsim and both numba
  engines.
- **Visualisation and export.**
  - Figure layer `figure_for(result, kind)` with PNG/SVG/PDF export and a Simulate ->
    Figures app page.
  - Engine2d carrier maps and a matplotlib-frames animation fallback.
  - A user-editable layout parameter file read by the GDSII, spec-sheet and cross-section
    exports.
  - Mesh export with VTU.
- **Lateral CCE weighting** of MC events (B-11) and a finite-SV Cauchy mean chord length
  in `run_microdosimetry` (D-15).

### Changed

- **Default damage model (#5, #8).** The default NIEL convention is now `sic-damage`,
  with the calibrated kinetic trap model. The previous behaviour stays available as an
  explicit option.
- **Default mobility (#6).** The default low-field mobility is the electron-only c-axis
  set. Holes keep the legacy isotropic set; holes x 0.90 is the opt-in
  `c_axis_holes_0p9`. CCE(LET) tables were rebuilt on the merged code.
- **Dark current (B-08).** The Hurkx 1992 trap-assisted-tunnelling field enhancement
  replaces the inert Kt branch, and N_t is refitted to the 18 pA leakage bound.
- **Transient (Mj-4).** The terminal current is now conduction + displacement, and
  `t_collection` is measured on it.
- **FLASH.** Auger recombination is genuinely active in `cce_vs_dose_rate` and
  `transient_cce_vs_dose_rate`.
- **`DeviceConfig`** has a strict doping-profile contract (C-12), and facades re-validate
  their input (C-03).
- **Performance.**
  - Engine-seeded devsim start (`seeded_bias`) replaces the 0 V bias ramp.
  - Compiled Gummel loops, a fused carrier solve and serial/prange dispatch.
  - The 1D field runs on the batched engine.
  - Lazy imports.
  - A cached UMFPACK setup in devsim's solver callback, which is bit-identical.
- **Package description.** It no longer names an institutional group. The sole copyright
  holder is Cesare Scalia (#3 H1).

### Fixed

- **Notebooks.** Fluence grids are expressed as fractions of each device's own
  compensation fluence. In notebook 14, the unsourced reference curve and wrong citations
  are removed.
- **Streamlit app.** It is aligned with the current library, with per-page model and
  validity notes. Every `PLACEHOLDER` keeps `REFERENCE-NOT-FAB`.

## Development history (before 5.0.0)

These were internal milestones in `.planning/MILESTONES.md`. They were never tagged or
released, and are listed only for provenance. Dates are those of the milestone commits.

- **v6.0 Audit-v7 hardening and performance engines (2026-10-19, PR #1).**
  - A batched numba float64 1D engine (C-V, CCE, dark current, FLASH) and a
    laterally-uniform 2D engine, with devsim as the fallback outside the validated box.
  - A new default stack (p+ 0.3 um / graded n- epi / n+) with a fitted n-/n+ transition,
    and the `measured_cv_tabulated_config`, `nplus_calibrated_config`, `legacy_config`
    and `microdosimeter_config` presets.
  - C = dQ/dV, `etna.run_all`, the Run All and Export pages, and `docs/USAGE.md`.
  - Public git history rewritten to remove copyrighted PDFs (map in
    `docs/history-rewrite-commit-map.txt`).
- **v5.0 Simulator library and Streamlit UI (2026-07-15).**
  - The installable `etna` package (`etna/core`, `etna/api`), renamed from `petringa`.
  - The Streamlit multi-page app and geometry viewer.
  - Notebook workflows reproducible from the UI.
- **v4.0 Scientific validation and extended physics (2026-05-17 to 2026-06-16, partially
  executed).**
  - Graded-doping 2D calibration (phase 26).
  - Stopping-power machinery.
  - Lazy `uproot` reading of Geant4 ROOT files.
- **v3.0 SiC microdosimeter design study (2026-03-27 to 2026-05-17).**
- **v2.0 Radiation damage modelling (2026-03-24 to 2026-03-27).** Burin 2024 defect
  constants, CCE vs fluence, damage dark current, carrier removal and C-V evolution, and
  annealing.
- **v1.1 Realistic device physics (2026-03-22 to 2026-03-24).** Temperature dependence,
  dark current with Hurkx TAT and surface recombination, and transient FLASH pulses.
- **v1.0 SiC TCAD simulator MVP (2026-03-20 to 2026-03-22).** 4H-SiC material module,
  devsim Poisson and drift-diffusion, C-V, CCE vs bias, and the first FLASH dose-rate
  study.

[Unreleased]: https://github.com/ngcex/sic-tcad-detector/compare/v5.0.0...HEAD
[5.0.0]: https://github.com/ngcex/sic-tcad-detector/releases/tag/v5.0.0
