# Project Research Summary

**Project:** Petringa SiC TCAD Simulator v1.1 -- Realistic Device Physics
**Domain:** 1D 4H-SiC semiconductor device simulation for radiation detector dosimetry
**Researched:** 2026-03-23
**Confidence:** MEDIUM-HIGH

## Executive Summary

The v1.1 milestone extends the validated 1D 4H-SiC TCAD simulator with three capabilities: temperature-dependent material parameters (303-313K clinical range), realistic dark current modeling (matching the experimental 18 pA), and transient FLASH pulse dynamics (10-200 ms pulses at 20-230 Gy/s). The existing v1.0 stack (devsim 2.10.0, numpy, scipy, matplotlib, lmfit) is entirely sufficient -- no new packages are needed. All new work is physics modeling code: extending `sic_material.py` with T-dependent functions, creating new modules for surface physics and transient simulation, and wiring these into the existing drift-diffusion pipeline. The architecture is additive: three new modules (`temperature.py`, `surface_physics.py`, `transient.py`) plug into the existing `device.py` -> `drift_diffusion.py` -> `charge_collection.py` flow without restructuring it.

The recommended approach is a four-phase build: (1) temperature-dependent parameters as foundation, (2) dark current physics (TAT + surface recombination), (3) transient FLASH dynamics, (4) combined analysis and publication figures. This ordering is driven by hard dependencies -- surface recombination and TAT both require T-dependent n_i(T), and transient simulation benefits from having all steady-state physics stabilized first. The v1.0 steady-state FLASH result serves as the critical validation anchor for transient work (transient CCE at long times must converge to the DC steady-state solution).

Three risks dominate. First, the dark current modeling: the experimental 18 pA cannot be explained by any n_i-proportional mechanism because SiC's n_i (~5e-9 cm^-3) makes all such currents negligibly small (~1e-25 A). The 1D model must rely on field-enhanced generation (trap-assisted tunneling) with effective fitting parameters, and the fit may not be physically unique. Second, transient solver stability across a 6-order timescale gap (ns carrier dynamics vs ms pulse durations) requires user-implemented adaptive time-stepping since devsim has no built-in adaptive dt. Third, threading temperature through the codebase risks silently breaking the validated 300K baseline -- the `compute_ni(300)` function returns ~6.5e-9 vs the calibrated constant 5e-9, requiring explicit reconciliation before use.

## Key Findings

### Recommended Stack

No new packages. The entire v1.1 feature set runs on the existing devsim 2.10.0 + scipy + numpy + matplotlib stack. devsim natively supports transient simulation (BDF1/BDF2/TR-BDF2 integrators), custom node models for TAT, and custom contact equations for surface recombination. The existing `time_node_model` registrations in `setup_sic_drift_diffusion()` already make the DD equations transient-capable -- calling `solve(type="transient_bdf1")` instead of `solve(type="dc")` is sufficient.

**Core technologies (unchanged from v1.0):**

- **devsim 2.10.0:** Drift-diffusion solver with transient BDF1/BDF2 already wired via `time_node_model` -- no equation changes needed
- **scipy/numpy:** Parameter fitting (curve_fit/minimize), ODE validation (solve_ivp for simplified transient cross-checks)
- **lmfit:** Constrained parameter fitting for dark current calibration (tau_SRH, S, m_t against 18 pA target)

### Expected Features

**Must have (table stakes):**

- T-dependent bandgap E_g(T), n_i(T), mobility mu(T), NC/NV(T) -- Varshni + Caughey-Thomas with 4H-SiC exponents from Ayalew thesis
- Dark current decomposition: bulk SRH + surface + TAT contributions identified and plotted separately
- Reverse I-V matching 18 pA experimental value at -60V with physics-based models
- Temperature coefficient of sensitivity matching experimental -0.079 +/- 0.005 %/C
- Transient current pulse I(t) with intra-pulse carrier build-up dynamics
- Validation against analytical limits (transient CCE converges to v1.0 steady-state CCE at low dose rate)

**Should have (differentiators -- publishable novelty):**

- Hurkx trap-assisted tunneling as field-enhanced SRH -- dominant dark current mechanism in SiC
- Surface recombination velocity at contacts -- SRV-based boundary condition replacing Ohmic pinning
- Inter-pulse memory effects in multi-pulse FLASH trains -- no published SiC TCAD study exists
- Plasma build-up and sweep-out time constants extracted from transient simulations
- T-dependent CCE reproducing the -0.079%/C coefficient from combined T + transient physics
- CCE per pulse in multi-pulse trains showing pulse-resolved dosimetric response

**Defer (v2+):**

- Coupled thermal simulation (clinical T range 30-40C too narrow; SiC thermal conductivity 3-5 W/cm-K prevents self-heating)
- Multi-level trap model (single midgap trap captures 90%+ of physics; multi-trap adds 15+ unconstrained parameters)
- Band-to-band tunneling (max field ~0.1 MV/cm, far below the ~3 MV/cm threshold for SiC BTBT)
- 2D/3D effects, noise simulation, radiation damage evolution, readout electronics modeling

### Architecture Approach

The v1.1 architecture adds three new modules to the existing linear pipeline without restructuring it. Temperature becomes a parameter threaded through all material computations via a `TemperatureDependentParams` frozen dataclass stored in `device_info["tdep"]`. Surface physics attaches at the contact boundary layer using devsim's `contact_equation()` API, following the proven pattern from the existing Auger implementation in `flash_recombination.py`. Transient solving wraps devsim's built-in BDF integrators with a Python-side adaptive time-stepping loop.

**Major components:**

1. **`temperature.py`** (NEW) -- Single source of truth for all T-dependent parameters; returns frozen dataclass with E_g, n_i, NC, NV, mu_n, mu_p, tau_n, tau_p, V_t at any T
2. **`surface_physics.py`** (NEW) -- Surface SRH recombination at contacts via `contact_equation()` + field-enhanced generation (TAT) in depletion region via modified USRH node model
3. **`transient.py`** (NEW) -- Pulse ON/OFF time-stepping loop with adaptive dt, CCE extraction per pulse, multi-pulse train orchestration
4. **`device.py`** (MODIFY) -- Wire T-dependent params from `temperature.py` instead of hardcoded 300K constants; store `tdep` in `device_info`
5. **`drift_diffusion.py`** (MINOR) -- Parameterize n1/p1 with T-dependent n_i; no equation changes

**Key architectural patterns:**

- Parameter threading via `device_info["tdep"]` dict prevents T-inconsistency between modules
- Contact equation modification follows the proven Auger pattern: CreateNodeModel -> CreateNodeModelDerivative -> re-register equation
- All new physics is additive and toggleable: disabled by default, v1.0 behavior preserved when off
- Device recreated at each T for parametric sweeps (simplest approach, avoids mutable state)

### Critical Pitfalls

1. **Scattered 300K constants (Critical)** -- n_i, mu, tau hardcoded at 300K in 6+ locations across `device.py`, `charge_collection.py`, `sic_material.py`. Missing any one produces silent physics errors. Prevention: centralize all T-dependent computation in `temperature.py`, full codebase audit of every `params.n_i_300` / `params.mu_n_max` usage, verify bit-identical results at T=300K.

2. **Surface recombination cannot explain 18 pA (Critical)** -- At n_i=5e-9, surface SRH gives ~1e-25 A -- 16 orders of magnitude below the 18 pA target. No physically reasonable SRV bridges this gap. Prevention: use field-enhanced generation (TAT) as the primary dark current mechanism. Do the quantitative math BEFORE implementing any dark current model.

3. **TAT complexity explosion (Critical)** -- Full multi-trap Hurkx model has 15-20 fitting parameters for a single I-V curve -- massively underconstrained. Prevention: start with simplest 2-parameter field-dependent generation (prefactor G0 + field scale E0). Add complexity only if I-V shape demands it. Ensure all Jacobian derivatives are analytically correct.

4. **Breaking validated 300K baseline (Critical)** -- `compute_ni(300)` returns ~6.5e-9 vs calibrated 5e-9 (30% discrepancy from different DOS masses). Switching silently changes all equilibrium results. Prevention: extract golden reference values BEFORE code changes, build `test_300K_regression()`, reconcile DOS masses to reproduce 5e-9 or re-validate.

5. **Transient timestep across 6-order timescale gap (Moderate)** -- ns carrier dynamics vs ms pulse durations with no built-in adaptive stepping. Prevention: implement adaptive dt in Python driver (start at 0.1 ns, increase 1.5x when Newton converges in <10 iterations, halve when >30). Use BDF1 at pulse on/off discontinuities, BDF2 during quasi-steady.

## Implications for Roadmap

Based on combined research findings, suggested phase structure:

### Phase 1: Temperature Foundation

**Rationale:** Every other v1.1 feature depends on T-dependent parameters being available. This is a pure computation module with no solver changes -- low risk, immediately testable against published literature tables. Must come first because both dark current (Phase 2) and transient (Phase 3) need n_i(T), mu(T), tau(T).
**Delivers:** `temperature.py` module with `compute_material_params(T, N_D)` returning all T-dependent quantities; modified `device.py` wiring T into the pipeline; modified `sic_material.py` with T-exponent constants; T-dependent I-V and C-V demonstrating expected shifts; `test_300K_regression()` ensuring no baseline breakage.
**Addresses:** T-dependent E_g, n_i, mu, NC/NV (all table stakes); foundation for temperature coefficient.
**Avoids:** Pitfall 1 (scattered 300K constants) via centralized computation; Pitfall 2 (wrong T-exponents) by using verified 4H-SiC Ayalew values (gamma_n=-2.40, gamma_p=-2.15); Pitfall 15 (n_i mismatch) by reconciling `compute_ni` with calibrated 5e-9; Pitfall 3 (regression) by extracting golden values first.
**Estimated scope:** ~300 LOC new + ~150 LOC modifications. Low-medium complexity.

### Phase 2: Dark Current Physics

**Rationale:** Depends on T-dependent n_i from Phase 1. The 40-order-of-magnitude gap between bulk SRH prediction (6.7e-49 A) and experiment (18 pA) is the central physics problem. Must be solved before transient work to have a complete, credible device model for publication.
**Delivers:** `surface_physics.py` with field-enhanced generation (TAT, Klaassen approximation) and surface recombination contact equations; dark current calibration matching 18 pA at -60V; decomposition plot showing I_bulk + I_surface + I_TAT vs voltage; I_dark(T) validation.
**Addresses:** Reverse I-V matching experiment (table stakes); dark current decomposition (table stakes); Hurkx TAT (differentiator); SRV at contacts (differentiator).
**Avoids:** Pitfall 4 (surface recombination too weak) by leading with TAT as primary mechanism; Pitfall 5 (TAT complexity) by starting with 2-parameter model; Pitfall 8 (contact vs region model confusion) by implementing surface physics at contact level only and verifying spatial localization.
**Estimated scope:** ~300 LOC new. Medium-high complexity -- custom devsim contact_equation patterns with analytical Jacobian derivatives.

### Phase 3: Transient FLASH Dynamics

**Rationale:** Depends on working DD pipeline from Phases 1-2. Most complex new feature -- benefits from having all steady-state physics stabilized. The v1.0 steady-state FLASH result provides the validation anchor (transient CCE at t >> tau_SRH must match DC CCE). This is the highest-novelty phase.
**Delivers:** `transient.py` with single-pulse I(t) waveforms, multi-pulse train simulation (N=10+), inter-pulse memory quantification, adaptive time-stepping; transient CCE extraction per pulse; plasma build-up and sweep-out time constant extraction; validation against v1.0 steady-state at low dose rate.
**Addresses:** Transient I(t) waveform (table stakes); intra-pulse dynamics (table stakes); inter-pulse memory (differentiator); plasma build-up time constants (differentiator); carrier sweep-out dynamics (differentiator).
**Avoids:** Pitfall 6 (timestep selection) via adaptive dt with convergence monitoring; Pitfall 13 (initial condition discontinuity) via mandatory DC init + transient_dc; Pitfall 14 (charge conservation drift) via per-step conservation accounting.
**Estimated scope:** ~400 LOC new. High complexity -- stiff multi-scale time-stepping over 6 orders of magnitude.

### Phase 4: Combined Analysis and Publication Figures

**Rationale:** All physics modules must be complete and individually validated before parametric studies combining T + transient + dark current. This phase is orchestration and visualization, not new physics.
**Delivers:** T-dependent CCE vs bias at T = 300, 303, 306, 310, 313 K; temperature coefficient extraction and validation against -0.079%/C; combined parametric sweeps (CCE vs dose_rate, T, bias, pulse structure); publication-quality transient waveform figures; comparison table of transient vs steady-state CCE across full parameter space.
**Addresses:** T-dependent CCE (differentiator); temperature coefficient validation (table stakes); CCE per pulse in multi-pulse train (differentiator); pulse structure optimization (differentiator).
**Avoids:** No new pitfalls -- uses validated pieces from earlier phases.
**Estimated scope:** 2-3 Jupyter notebooks. Medium complexity (orchestration and visualization).

### Phase Ordering Rationale

- **Hard dependency chain:** Phase 1 (T-params) feeds Phase 2 (dark current needs n_i(T)) and Phase 3 (transient needs T-dependent material properties). Phase 4 integrates all three.
- **Validation anchoring:** Each phase has an unambiguous pass/fail criterion: Phase 1 matches Ayalew/Ioffe tables and preserves v1.0 baseline; Phase 2 matches 18 pA at -60V; Phase 3 converges to v1.0 steady-state at low dose rate; Phase 4 reproduces -0.079%/C.
- **Risk front-loading:** The hardest physics question (can a 1D model match 18 pA with physical parameters?) lands in Phase 2, early enough to pivot strategy if TAT alone is insufficient.
- **Regression safety:** Phase 1 mandates golden value extraction and regression test creation before any physics changes touch the codebase.
- **Scientific narrative:** The phase ordering mirrors the paper narrative -- establish T-dependent model, validate against dark current, predict novel transient behavior, combine for complete dosimetric characterization.

### Research Flags

Phases likely needing deeper research during planning:

- **Phase 2 (Dark Current):** The quantitative analysis proves surface SRH is 16 orders of magnitude too weak. TAT implementation in devsim requires custom Jacobian derivatives with no built-in support, no examples, and sparse forum discussion. The Klaassen approximation for Gamma(F) needs numerical validation. Key decision: fit to 18 pA with effective parameters, or attempt full physical decomposition (which may require accepting that perimeter leakage is unmodellable in 1D)?
- **Phase 3 (Transient):** Adaptive time-stepping strategy needs empirical validation -- devsim's Python-level time loop has no precedent for FLASH-duration (200 ms) pulses. Need to determine whether TR-BDF2 composite method is necessary or if BDF1 alone suffices for the generation on/off discontinuity. Computational cost estimate (30s/pulse) is unverified.

Phases with standard patterns (skip research-phase):

- **Phase 1 (Temperature):** Varshni, Caughey-Thomas T-scaling, and DOS T^1.5 are textbook models with well-documented 4H-SiC parameters from the Ayalew thesis. The implementation is parameter wiring into an existing pipeline.
- **Phase 4 (Combined Analysis):** Pure orchestration of validated components using existing parametric sweep patterns from v1.0. No new physics or solver patterns needed.

## Confidence Assessment

| Area         | Confidence  | Notes                                                                                                                                                                                                                                 |
| ------------ | ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Stack        | HIGH        | No new packages needed; all devsim transient capabilities verified from installed examples (tran_diode.py, transient_rc.py) and documentation                                                                                         |
| Features     | HIGH        | Feature set derived from experimental validation targets in Petringa papers and Lopez Paz 2024; clear quantitative pass/fail criteria for every feature                                                                               |
| Architecture | HIGH        | Additive design following existing proven patterns (Auger implementation as template); devsim transient API confirmed from source code                                                                                                |
| Pitfalls     | MEDIUM-HIGH | Codebase-specific pitfalls verified against actual source (line-level references); physics pitfalls confirmed by quantitative estimates (surface SRH math); TAT fitting risk is real but bounded by the 2-parameter starting approach |

**Overall confidence:** MEDIUM-HIGH

The stack, features, and architecture are well-characterized with high-quality sources. The main uncertainties are: (1) whether a 1D model can physically reproduce 18 pA dark current with meaningful (not purely curve-fitting) parameters, and (2) whether devsim's transient solver handles the FLASH timescale gap (ns to ms) without prohibitive computational cost or convergence failure at high injection.

### Gaps to Address

- **n_i(300K) discrepancy:** `compute_ni(300)` returns ~6.5e-9 vs calibrated 5e-9. Must reconcile (adjust DOS effective masses to reproduce 5e-9, or accept new value and fully re-validate all v1.0 results). This is a Phase 1 blocker -- resolve before any T-dependent work.
- **Dark current mechanism ambiguity:** The 18 pA could be perimeter leakage (inherently 2D, unmodellable in 1D), TAT through Z1/2 centers (1D-capable), or measurement artifact (cable/probe leakage). If TAT with physical parameters cannot match 18 pA, strategy must shift to effective parameters with documented limitations.
- **SRH lifetime T-dependence model:** Literature is contradictory for 4H-SiC -- the Z1/2 center shows sample-dependent tau(T) behavior. For the 303-313K clinical range, the effect is <5% regardless of model choice. Recommend: keep tau constant for v1.1, document the assumption.
- **Donor incomplete ionization:** N donors at E_C - 0.05/0.09 eV are ~85% ionized at 300K but not modeled (full ionization assumed). Effect is <3% in clinical range. Add `ionized_donor_concentration()` for completeness but do not iterate self-consistently.
- **Transient computational cost:** Estimated ~30s per pulse, ~5 min for 10-pulse train at one dose rate. Acceptable but unverified. If convergence requires smaller dt than estimated, cost could increase 10x, making parametric sweeps over dose rate and T expensive.

## Sources

### Primary (HIGH confidence)

- devsim 2.10.0 installed examples (`tran_diode.py`, `transient_rc.py`, `simple_physics.py`) -- transient solver patterns, contact equation API, time_node_model usage
- devsim documentation (devsim.net/solver.html, devsim.net/CommandReference.html, devsim.net/models.html) -- solve types, tdelta/charge_error parameters, custom equation API
- TU Wien Ayalew thesis, Table 3.5 -- 4H-SiC Caughey-Thomas T-dependent parameters: gamma_n=-2.40, gamma_p=-2.15, beta=-0.5
- Petringa group papers (SiC_Photons_MedicalPhysics, Microdosimetry.pdf, Flash.pdf) -- experimental validation targets (18 pA dark current, CCE, C-V)
- Lopez Paz et al., Med Phys 2024 -- 11 Gy/pulse linearity, -0.079%/C temperature coefficient, <2% sensitivity drift after 100 kGy, 4 MGy/s capability
- Hurkx et al., IEEE TED 39(2), 1992 -- original field-enhanced SRH (TAT) model formulation
- TU Wien Schleich thesis -- TAT current implementation for SiC TCAD

### Secondary (MEDIUM confidence)

- Kimoto et al., J. Appl. Phys. 127, 195702 (2020) -- SRV for SiO2-passivated 4H-SiC Si-face: 150-5000 cm/s
- arXiv:2503.09016 (2025) -- 4H-SiC TAT tunneling effective mass m_t = 0.25\*m_0
- Burin et al., CERN RD50 (2024) -- comprehensive 4H-SiC TCAD parameter survey
- IEEE Access 2024 (10538275) -- carrier lifetime vs T in 4H-SiC, power-law dependence
- Frontiers in Sensors 2025 (1622153) -- SiC detector FLASH challenges review
- Hatakeyama et al., Materials Science Forum (2013) -- alternative 4H-SiC mobility T-exponents
- Rakheja et al., Semicond. Sci. Technol. (2023) -- SRV vs carrier concentration
- devsim Forum -- surface recombination implementation discussion (sparse)

### Tertiary (LOW confidence)

- SRH lifetime T-scaling model (tau ~ T^1.5 phonon-assisted) -- generic semiconductor model, not SiC-specific validated
- Auger coefficient T-dependence for 4H-SiC -- poorly characterized in literature, assumed negligible for v1.1

---

_Research completed: 2026-03-23_
_Ready for roadmap: yes_
