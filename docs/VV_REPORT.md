# ETNA verification and validation report

Status: 2026-10-21, etna 5.0.0 (main after the issue #9 merge). Issue #3, P2.1 / P2.2.

**Read this first.** ETNA is a research TCAD code for 4H-SiC detectors. It is not
clinical software. It makes no medical-device, certification, dosimetry-protocol or
regulatory claim, and no standard certifies a TCAD code for dosimetry. This report
collects the credibility evidence that exists today, in the vocabulary of two
frameworks:

- **NASA-STD-7009** (Standard for Models and Simulations, Rev. B; listing opened, see
  issue #3) for the structure: the credibility factors of section 2.
- **ASME V&V 20** for what "validation" means. A validation comparison is a measured
  model error, E = S − D (simulation minus data), stated with the uncertainty of the
  comparison. It holds at the validation point only, and it is **not a pass/fail
  verdict**. Source: Eça, Dowding, Roache, ASME VVS2020, doi:10.1115/VVS2020-8826
  (abstract opened, as cited in issue #3).

The status classes used below and in the registry are:

| class             | meaning                                                                                               |
| ----------------- | ----------------------------------------------------------------------------------------------------- |
| **VERIFIED**      | code vs an analytic limit (correctness of the numerics)                                               |
| **CROSS-CHECKED** | code vs code: an independent solver or reference implementation                                       |
| **CALIBRATED**    | agreement with the data the model, or its input, was fitted to. This is fit quality, never validation |
| **VALIDATED**     | comparison with a measurement that was not fitted; E is reported, not judged                          |
| **EXPECTED-FAIL** | a documented disagreement with a measurement                                                          |

No CALIBRATED entry is relabelled VALIDATED anywhere in this report.
`tests/test_benchmarks.py` freezes those labels and checks that the evidence table in
section 5 is generated from the registry.

Every number in sections 4 and 5 comes from one of two places. Some come from a
benchmark-registry entry recomputed on 2026-10-21 (`uv run python -m benchmarks`). The
rest come from a lab LOG or issue comment named next to the number; those were read
again for this report but not re-run, and are marked as such.

---

## 1. Intended use

**May be used for:**

- Design studies of the two Petringa-group 4H-SiC devices (the 4 mm² p⁺/n⁻/n⁺ dosimeter
  and the microdosimeter), including:
  - relative comparisons between geometries and dopings;
  - electrostatics: C-V, depletion and field;
  - drift-diffusion CCE under low injection;
  - order-of-magnitude dark current.
- Ranking alternatives and planning measurements (what to measure, and which measurement
  would discriminate between models: [USAGE §8](USAGE.md#8-validity-and-limits)).
- Code-to-code reference values inside the cross-checked parameter boxes
  ([USAGE §7](USAGE.md#7-solver-provenance-and-fast-paths)).
- Teaching and method development (track engine, damage kinetics, avalanche library), as
  opt-in models with stated assumptions.

**Must not be used for:**

- Clinical decisions, treatment planning, patient dosimetry or beam calibration. The
  same applies to any IEC 60731 / IAEA TRS-398 dosimeter classification: no output of
  ETNA has been validated against a measurement of the dosimeter it models.
- Absolute predictions of these quantities, which are calibrations or extrapolations:
  - leakage current, a one-point fit to an upper bound;
  - radiation-damaged CCE, a model prediction with an unidentified removal mechanism;
  - breakdown voltage, a parallel-plane upper bound with no edge termination;
  - the pulse-height defect, conditional on an unknown track-core radius;
  - FLASH CCE, a numerical sensitivity bound.
- Foundry sign-off of mask dimensions. The reference layout is `REFERENCE-NOT-FAB`
  while PLACEHOLDER dimensions remain (README "Foundry layout parameters").

The per-output trust table is in [README "What to trust"](../README.md#what-to-trust)
and is not repeated here.

---

## 2. Credibility factors (NASA-STD-7009)

The eight-factor list below follows NASA-STD-7009A's credibility assessment. Issue #3
found the list only in a search snippet of a NASA NTRS paper, so its exact wording is
**unverified**. The factors are verification, validation, input pedigree, results
uncertainty, results robustness, use history, M&S management and people
qualifications. The assessment is qualitative. ETNA does not compute a numeric
credibility score (issue #3, open question 5).

| factor                    | evidence                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | assessment                                                                                                                                                                                       |
| ------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Verification**          | 12 analytic-limit and 12 code-to-code registry entries (sections 4 and 5). devsim parity gates across the unit tests. The mesh and time-step studies of section 4                                                                                                                                                                                                                                                                                                                                              | Strong for the numerics inside the boxes. Gaps: no method-of-manufactured-solutions (MMS) test for the 1D Poisson / drift-diffusion solvers; the 1D production mesh is first order (section 4.2) |
| **Validation**            | One electrostatic fit (CALIBRATED), one leakage bound (CALIBRATED), damage removal fitted to other devices (CALIBRATED). Independent comparisons exist only for the opt-in breakdown library (other diodes) and qualitative damage CCE on other devices. The one device-level comparison that was not fitted, microdosimeter CCE(V), fails (EXPECTED-FAIL)                                                                                                                                                     | **Weak.** No output is validated on the modelled devices                                                                                                                                         |
| **Input pedigree**        | NIST PSTAR/ASTAR tables with SHA-256 of the raw responses ([data/nist_star/README.md](../data/nist_star/README.md)). SR-NIEL and NJOY SiC NIEL ([data/niel/README.md](../data/niel/README.md)). Literature parameter sets labelled primary / secondary ([USAGE §18](USAGE.md#provenance-what-was-opened)). Doping from one digitised C-V figure. Layout dimensions classed DOCUMENTED / MEASURED-FROM-FIGURE / PLACEHOLDER ([deliverables/DIMENSIONS-provenance.md](../deliverables/DIMENSIONS-provenance.md)) | Mixed. Material data are traceable. Device inputs (doping, lifetimes, area, mask) are inferred or placeholders, so a SIMS profile and mask data are the main missing items (section 7)           |
| **Results uncertainty**   | Numerical error estimates exist for C, W and CCE (section 4) and for V_BD (`mesh_error_V`). Parameter spreads are documented for impact ionization (about ±20 % between sets), the NIEL convention, the doping inversion (V_fd 2.7 vs 3.9-7.6 V if the bulk is at the documented doping) and τ (17-42 ns band)                                                                                                                                                                                                 | Partial. No propagated input-uncertainty analysis (UQ) of any output                                                                                                                             |
| **Results robustness**    | Default flips are pinned before and after: the mobility, damage-model and NIEL-convention flips of 2026-10-21 ([lab/README.md](../lab/README.md)). Sensitivity studies in the lab LOGs                                                                                                                                                                                                                                                                                                                         | Partial. Sensitivities are documented per study, not as a systematic analysis                                                                                                                    |
| **Use history**           | One research group; git history from 2026-03-20 (982 commits). No external users                                                                                                                                                                                                                                                                                                                                                                                                                               | Low                                                                                                                                                                                              |
| **M&S management**        | Git history, CI (lock check, ruff, notebook hygiene, the default suite), a nightly exhaustive run with pip-audit, the wheel-smoke job. Provenance recorded per calibrated constant in `lab/*/notes/LOG.md`. This registry                                                                                                                                                                                                                                                                                      | Good for a research code. No tagged release, CHANGELOG or DOI yet (issue #3 P1)                                                                                                                  |
| **People qualifications** | Developed with AI coding agents under one maintainer's review. Physics audits in `.planning/` and `.issues/`                                                                                                                                                                                                                                                                                                                                                                                                   | Not independently reviewed. A JOSS or methods-paper review (issue #3 P5) would address this                                                                                                      |

---

## 3. What is verified, calibrated, validated: summary by observable

| observable (default device unless stated)           | numerics                                                                               | physics                                                                                                                                                  | registry                                                                                                         |
| --------------------------------------------------- | -------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| C-V, W, field                                       | VERIFIED vs abrupt-junction and Mott-Schottky limits; CROSS-CHECKED vs devsim (1D, 2D) | **CALIBRATED**: doping fitted to Photons Fig. 4b. Plateau rms 0.013 µm; knee **EXPECTED-FAIL** (rms 0.24 µm against σ 0.03)                              | AN-W-UNIFORM, AN-MS-ND, AN-VFD-KNEE, CC-CV-DEVSIM, CC-E2D-PSI, CA-PHOTONS-PLATEAU, CA-W-LABEL-*, EF-PHOTONS-KNEE |
| electrostatics solver scale on another device       | —                                                                                      | **CALIBRATED** (circular: the paper's own N_D fit of the same C-V); +1.0 % with V_bi matched                                                             | CA-BRUZZI-W40-*                                                                                                  |
| CCE vs bias (steady, low injection)                 | CROSS-CHECKED vs devsim (1D, 2D, fv2d)                                                 | not compared on the dosimeter. Microdosimeter tilted alpha: **EXPECTED-FAIL** (shape rms 0.035 against σ 0.0005)                                         | CC-CCE-DEVSIM, CC-FV2D-_, EF-TILTED-ALPHA-_                                                                      |
| Hecht CCE (`run_temperature_sweep(method="hecht")`) | closed form VERIFIED                                                                   | a formula. Drift-diffusion with SRH does not reduce to it away from CCE → 1 (section 4.3)                                                                | AN-HECHT                                                                                                         |
| dark current                                        | VERIFIED: SRH generation limit; CROSS-CHECKED vs devsim                                | **CALIBRATED**: N_t fitted to the < 18 pA bound at −60 V. Bias and T shape not validated                                                                 | AN-SRH-GEN, CC-DARK-DEVSIM, CA-DARK-60                                                                           |
| damage: C-V removal                                 | —                                                                                      | **CALIBRATED** to the Radmanovac Table 2 slopes (other devices, 252.7 MeV)                                                                               | CA-RADMANOVAC-ONSEMI                                                                                             |
| damage: CCE vs fluence                              | engine vs devsim (tests/test_fluence_engine.py)                                        | **VALIDATED (bound only)**: Gsponer-like device consistent with "> 50 % at 1e15 neq". **EXPECTED-FAIL**: De Napoli-like device, loss about 10x too small | VA-GSPONER-1E15, EF-DENAPOLI-LOSS-*                                                                              |
| breakdown (opt-in library)                          | VERIFIED: uniform-field limit; McIntyre closed forms; DD gain vs McIntyre              | hatakeyama2004 **CALIBRATED** (fitted on those diodes); niwa2015 **VALIDATED** on Hatakeyama's diodes (−6.4 / +5.2 %)                                    | AN-BD-UNIFORM, AN-MCINTYRE-_, CC-AVAL-MCINTYRE, CA-HATAKEYAMA-1E17, VA-NIWA-_                                    |
| stopping power, κ(E)                                | CROSS-CHECKED: CSDA integral vs NIST, κ pipeline                                       | data (NIST) + **unvalidated** Bragg additivity for SiC                                                                                                   | CC-CSDA-*, CC-KAPPA-62                                                                                           |
| energy-loss straggling (MC sampler)                 | VERIFIED: Bohr variance (+3.9 %, model form), CSDA mean, Vavilov cumulant              | not compared with a measured spectrum                                                                                                                    | AN-BOHR-VAR, AN-MEAN-DEP, AN-VAVILOV                                                                             |
| single-ion track, PHD (opt-in)                      | ledgers, energy theorem, dark state vs batched (section 4.4)                           | Fig. 4b tail **CALIBRATED** (τ_scale 0.039); all 7 points **EXPECTED-FAIL**                                                                              | CA-TILTED-TAU-SCALE, EF-TILTED-TRACK-ALL7                                                                        |

Mapping to the characterisation quantities of the IEC 60731 draft (prEN IEC 60731:2026,
EVS listing opened in issue #3). The draft wording is not a requirement on ETNA:

| IEC 60731 draft quantity             | ETNA output                                         | status                                                                 |
| ------------------------------------ | --------------------------------------------------- | ---------------------------------------------------------------------- |
| leakage current                      | `run_dark_current`                                  | CALIBRATED (one bound)                                                 |
| stability under bias and temperature | `run_temperature_sweep`, `run_cv`/`run_cce` vs bias | model only; T shape not validated                                      |
| radiation-quality dependence         | CCE(LET) tables, κ(E), track engine                 | κ: data + approximation; CCE(LET): model; PHD: conditional on `radial` |
| zero-bias operation                  | `run_cce(v_start=0)`                                | model; 0 V tilted-alpha data not reproduced (EXPECTED-FAIL)            |

---

## 4. Numerical verification

### 4.1 Code-to-code (devsim parity)

The engines solve devsim's discretisation (box method, Scharfetter-Gummel, SRH) on the
same mesh, so parity tests show that the same equations are solved. They say nothing
about the device. Measured in this run (registry, `uv run python -m benchmarks`):

| comparison                                          | measured                                                      | gate   |
| --------------------------------------------------- | ------------------------------------------------------------- | ------ |
| 1D C-V, batched vs devsim, 0..−20 V                 | max dC/C 4.0e-8 (devsim's ±1 mV central difference); dW 4e-16 | 1e-6   |
| 1D C-V to −100 V (slow)                             | 9.0e-9                                                        | 1e-6   |
| 1D CCE, −10..−40 V                                  | 3.0e-15                                                       | 1e-5   |
| 1D dark current, 0..−100 V                          | 2.0e-14 relative                                              | 1e-9   |
| 2D potential, engine2d vs devsim, hw 25 µm, −20 V   | 3.6e-15 V (dE/E_max 8.9e-15)                                  | 1e-6 V |
| fv2d vs stored devsim, mesa, 50 V                   | 3.3e-16                                                       | 1e-10  |
| fv2d vs stored devsim, guard ring, 10 V (slow)      | 8.9e-16                                                       | 1e-10  |
| fv2d vs engine2d, default n+ diode, 20 V            | 3.4e-15                                                       | 1e-9   |
| DD avalanche gain vs McIntyre M, 0.99 V_BD (M ≈ 34) | 4.3e-6 in (1 − 1/M)                                           | 1e-4   |

Other parity results, recorded but not re-run here:

- 3D electrode (r, z): fv2d vs devsim 2.4e-12 (issue #4).
- FLASH CCE: ≤ 7.5e-12 (USAGE §7).
- Physics-model engine vs devsim: 3e-12 at 300 V (USAGE §18).
- Track engine dark state vs the batched engine: ψ 1e-14 V; steady CCE ≤ 3e-11 (issue #9).

### 4.2 Solution verification: mesh and tolerance studies

| solver / quantity                             | study                                                          | result                                                                                                                                                                                                                                       | source                                                                                              |
| --------------------------------------------- | -------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| 1D production mesh (684 nodes), default stack | uniform bisection, 4 levels, Richardson                        | **first order**, at the junction and depletion-edge cells. Level-0 error: CCE +1.7e-4 (−1 V) to < 1e-6 (≥ −5 V); C −2.7e-3 (−1 V) to 5e-5 (−50 V); W 3e-3 to 2e-6. Below the C-V digitisation σ (0.02 µm); the C-V fit was made on this mesh | [lab/mesh_contact/notes/LOG.md](../lab/mesh_contact/notes/LOG.md) §1                                |
| 1D contact cells, default stack               | bisection 1-7x                                                 | CCE ≤ 2.3e-10, C 0: not affected (no generation on contacts; the n⁻/n⁺ junction reflects holes)                                                                                                                                              | same, §1                                                                                            |
| legacy 1D model (cathode on the epi)          | same                                                           | CCE about −1e-3 at full depletion. Kept for reproducibility                                                                                                                                                                                  | same; issue #4                                                                                      |
| 2D contact-on-epi cathode (ALT_STACK)         | uniform bisection to 17.5k lines from two starting meshes      | first order in the contact cell. Converged mesa centre CCE 0.9630651 (10 V), 0.9990906 (50 V), ±1e-7. The legacy mesh was +6.17e-3 / +7.86e-4 off; the graded default mesh is +5.4e-5 / −2.9e-6                                              | lab/mesh_contact §2; [lab/geometry2d/notes/LOG.md](../lab/geometry2d/notes/LOG.md) §§5, 9; issue #4 |
| fv2d adaptive refinement                      | AMR with a goal-oriented estimate                              | the first-pass AMR "converged" values were off by −1.2e-4 / −2.6e-4. The geometric-extrapolation estimate underestimated the error 15-20x (caveat in `refine.py`). AMR needs the contact-cell rule                                           | issue #4                                                                                            |
| fv2d box scheme                               | coaxial Laplace, uniform r                                     | observed order 1.99 (registry AN-FV2D-ORDER); exact on a geometric mesh                                                                                                                                                                      | tests/test_fv2d.py                                                                                  |
| default stack 2D (engine2d / devsim)          | 0.1 nm cathode cell                                            | CCE moves 4e-10: not affected                                                                                                                                                                                                                | lab/mesh_contact §1                                                                                 |
| breakdown V_BD                                | full re-solve h, h/2, h/4                                      | about first order (ratio 1.5-2.3). Default dosimeter mesh error ≤ 0.11 V; the 2h quadrature estimate underestimates it                                                                                                                       | [lab/physics_models/notes/LOG.md](../lab/physics_models/notes/LOG.md); USAGE §18                    |
| c-axis mobility change                        | h, h/2, h/4                                                    | the change is converged to 0.3 % (−2 V) and 0.06 % (−10 V)                                                                                                                                                                                   | [lab/mobility_caxis/notes/LOG.md](../lab/mobility_caxis/notes/LOG.md) §10                           |
| small-signal C                                | engine linearisation vs devsim central difference and AC solve | 2.5e-8 (0 V) to 1.5e-12 (−30 V); Richardson 5e-11                                                                                                                                                                                            | [lab/microdosimeter/notes/LOG.md](../lab/microdosimeter/notes/LOG.md)                               |
| track engine                                  | time step: TR-BDF2 local-truncation-error control              | dCCE 1.9e-7 between rtol 1e-4 and 1e-6. Poynting residual 1.4e-4 → 1.8e-6 as rtol goes from 1e-3 to 1e-6 (time-discretisation error)                                                                                                         | issue #9; [USAGE §19](USAGE.md#19-single-ion-track-engine-opt-in-etnacoretrack)                     |
| track engine                                  | mesh FAST → fine                                               | PHD changes 5-8 %: use DEFAULT or finer for PHD numbers. Sub-nm cores are not mesh-converged and are outside drift-diffusion validity                                                                                                        | issue #9; [lab/track_plasma/notes/LOG.md](../lab/track_plasma/notes/LOG.md)                         |
| track engine                                  | float64 floor                                                  | excess unknowns: pair-ledger floor 9.6e-7 → 9e-14 (LET 1, 0 V); CCE vs the absolute form ≤ 6e-9 (LET 1000)                                                                                                                                   | USAGE §19.1                                                                                         |
| heat solver                                   | analytic line source                                           | within 3 %                                                                                                                                                                                                                                   | issue #9                                                                                            |

### 4.3 Analytic limits (registry, recomputed today)

| id                  | check                                                                         | E               | gate  |
| ------------------- | ----------------------------------------------------------------------------- | --------------- | ----- |
| AN-W-UNIFORM        | ε·A/C at −40 V vs the abrupt junction with the 2kT/q and p⁺-share corrections | +4.7e-6         | 1e-3  |
| AN-MS-ND            | N_D from the 1/C² slope                                                       | −1.6e-3         | 5e-3  |
| AN-VFD-KNEE         | Mott-Schottky knee vs q N_D d²/2ε − V_bi + 2kT/q (1e14, 10 µm)                | +1.2e-2         | 3e-2  |
| AN-MCINTYRE-AB / -N | M = 1/(1 − αW), M_n = e^(αW)                                                  | ≤ 1.6e-15       | 1e-12 |
| AN-BD-UNIFORM       | V_BD = bW/ln(aW), nearly uniform field                                        | +2.5e-3         | 3e-2  |
| AN-HECHT            | Hecht closed form vs quadrature                                               | +4.4e-9         | 1e-6  |
| AN-SRH-GEN          | q n_i W_eff/(τ_n + τ_p), pure midgap SRH                                      | +6.9e-3         | 3e-2  |
| AN-BOHR-VAR         | sampler variance vs ξ T_max (1 − β²/2), κ ≈ 22                                | +3.9e-2         | 0.2   |
| AN-MEAN-DEP         | sampler mean vs CSDA integral                                                 | −5.0e-6         | 2e-3  |
| AN-VAVILOV          | second cumulant of the Vavilov density                                        | −4.2e-9         | 1e-5  |
| AN-FV2D-ORDER       | observed order, coaxial Laplace                                               | 1.99 (E −0.012) | 0.2   |

Findings worth knowing:

- **SRH drift-diffusion does not reduce to Hecht** away from CCE → 1. On a nearly
  uniform-field diode with τ = 2 ps, drift-diffusion gives 1.75-2.6x the Hecht CCE,
  because SRH recombination needs both carriers in the same place while Hecht traps each
  carrier independently. `run_temperature_sweep(method="hecht")` is the closed form,
  verified as a formula only. The existing DD-vs-Hecht test compares them only near
  CCE = 1 (abs 0.10).
- **The naive q n_i W/2τ** overestimates the generation current by 11 % (−40 V) to 42 % (−5 V).
  Generation stops where n or p exceeds n_i near the depletion edges; the edge-corrected
  width agrees to 1-2.5 %. With the real n_i (5e-9 cm⁻³) the SRH current is about 1e-49 A,
  which is why the leakage model is a calibrated trap-assisted term.
- **The sampler's energy-loss variance** is 3.9 % above Bohr in the Gaussian regime. That
  is about 6x the sampling error at this seed, so it is a model-form property of the
  Geant4 Urban/Universal port. At κ < 1 the ratio is 0.76-0.85 (heavy tail).
- **MMS is missing.** No method-of-manufactured-solutions or order-of-accuracy test
  exists for the 1D Poisson / drift-diffusion solvers. Their order is known only from the
  production-mesh bisection study (first order, section 4.2). This is an open item
  (section 6).

### 4.4 Conservation and internal consistency

- fv2d three-contact current conservation: 2e-16 (issue #4).
- Track pair ledger: 1e-9 to 1e-10 at LET ≥ 100; thermal ledger ≤ 1e-7 (issue #9,
  USAGE §19).
- Ionization-integral identity 1 − I_n = e^(−φ(W)) (1 − I_p): holds to round-off
  (tests/test_physics_models.py).

---

## 5. Validation evidence (registry)

The table below is generated from `benchmarks/registry.toml` by
`uv run python -m benchmarks --update-report`. Do not edit it by hand; the test fails
if it drifts.

Columns:

- **reference D (u_D)**: D is `computed` when the evaluator computes it (an analytic
  formula or a devsim run).
- **E recorded**: the comparison error pinned on 2026-10-21 for the non-verification
  classes. For VERIFIED / CROSS-CHECKED the measured E is in `benchmarks/results/`.
- **gate**: for CALIBRATED / VALIDATED / EXPECTED-FAIL this is only a regression guard
  on E. It is not an acceptance criterion.

<!-- BEGIN benchmark-registry (generated by `python -m benchmarks --update-report`) -->

| id | status | quantity | reference D (u_D) | E recorded | gate | run |
|---|---|---|---|---|---|---|
| AN-W-UNIFORM | **VERIFIED** | small-signal depletion width eps*A/C at -40 V, uniform 1e15, 30 um [um] | computed | see results | \|E\| <= 1e-03 (rel) | fast |
| AN-MS-ND | **VERIFIED** | N_D from the Mott-Schottky 1/C^2 slope (-5..-40 V), uniform 1e15 [cm^-3] | computed | see results | \|E\| <= 5e-03 (rel) | fast |
| AN-VFD-KNEE | **VERIFIED** | full-depletion voltage of a uniform 1e14, 10 um epi (no n+) [V] | computed | see results | \|E\| <= 3e-02 (rel) | fast |
| AN-MCINTYRE-AB | **VERIFIED** | multiplication M(x) for alpha = beta, alpha W = 0.5, constant field | computed | see results | \|E\| <= 1e-12 (rel) | fast |
| AN-MCINTYRE-N | **VERIFIED** | electron multiplication M_n, only electrons ionize, alpha W = 0.5 | computed | see results | \|E\| <= 1e-12 (rel) | fast |
| AN-BD-UNIFORM | **VERIFIED** | avalanche breakdown voltage, p+/1e12/n+ (nearly constant field), alpha = beta = a exp(-b/E) [V] | computed | see results | \|E\| <= 3e-02 (rel) | fast |
| AN-HECHT | **VERIFIED** | two-carrier depth-averaged Hecht CCE, tau 2 ps, 50 V, 10 um | computed | see results | \|E\| <= 1e-06 (rel) | fast |
| AN-SRH-GEN | **VERIFIED** | midgap SRH generation current density at -40 V, p+/1e15/10 um, n_i inflated to 1e8 [A/cm^2] | computed | see results | \|E\| <= 3e-02 (rel) | fast |
| AN-BOHR-VAR | **VERIFIED** | energy-loss variance, 5.486 MeV alpha, 2 um SiC (kappa ~ 22), 50k events, seed 2 [keV^2] | computed | see results | \|E\| <= 2e-01 (rel) | fast |
| AN-MEAN-DEP | **VERIFIED** | mean energy deposit, same sample as AN-BOHR-VAR [keV] | computed | see results | \|E\| <= 2e-03 (rel) | fast |
| AN-VAVILOV | **VERIFIED** | second cumulant of the Vavilov density, kappa 0.5, beta^2 0.3 | computed | see results | \|E\| <= 1e-05 (rel) | fast |
| AN-FV2D-ORDER | **VERIFIED** | observed order of accuracy, fv2d (r, z) coaxial Laplace on uniform r meshes 41 -> 81 | computed | see results | \|E\| <= 2e-01 (abs) | fast |
| CC-CV-DEVSIM | **CROSS-CHECKED** | max \|C_engine/C_devsim - 1\|, default stack, 0..-20 V (2.5 V lattice) | computed | see results | \|E\| <= 1e-06 (rel) | fast |
| CC-CCE-DEVSIM | **CROSS-CHECKED** | max \|CCE_engine - CCE_devsim\|, default stack, -10..-40 V, alpha | 0 | see results | \|E\| <= 1e-05 (abs) | fast |
| CC-DARK-DEVSIM | **CROSS-CHECKED** | max \|I_engine/I_devsim - 1\|, default stack dark current 0..-100 V | computed | see results | \|E\| <= 1e-09 (rel) | fast |
| CC-E2D-PSI | **CROSS-CHECKED** | max \|psi_engine2d - psi_devsim\|, 2D default stack, half width 25 um, -20 V [V] | 0 | see results | \|E\| <= 1e-06 (abs) | fast |
| CC-FV2D-MESA | **CROSS-CHECKED** | max \|dCCE\| fv2d vs stored devsim, mesa (ALT_STACK, legacy depth mesh), 50 V, 9-point scan | 0 | see results | \|E\| <= 1e-10 (abs) | fast |
| CC-FV2D-E2D | **CROSS-CHECKED** | max \|dCCE\| fv2d vs engine2d, default n+ diode, hw 25 um, 20 V | 0 | see results | \|E\| <= 1e-09 (abs) | fast |
| CC-AVAL-MCINTYRE | **CROSS-CHECKED** | (1 - 1/gain_DD) / (1 - 1/M_McIntyre), default dosimeter at 0.99 V_BD (M ~ 34) | computed | see results | \|E\| <= 1e-04 (rel) | fast |
| CC-CSDA-P62 | **CROSS-CHECKED** | CSDA range of 62 MeV protons in Si [g/cm^2] | computed | see results | \|E\| <= 2e-03 (rel) | fast |
| CC-CSDA-A55 | **CROSS-CHECKED** | CSDA range of 5.5 MeV alphas in Si [g/cm^2] | computed | see results | \|E\| <= 2e-03 (rel) | fast |
| CC-KAPPA-62 | **CROSS-CHECKED** | SiC -> water linear LET ratio kappa at 62 MeV (protons) | computed | see results | \|E\| <= 1e-06 (rel) | fast |
| CC-CV-DEVSIM-DEEP | **CROSS-CHECKED** | max \|C_engine/C_devsim - 1\|, default stack, 0..-100 V (5 V lattice, 21 points) | computed | see results | \|E\| <= 1e-06 (rel) | slow |
| CC-FV2D-GUARD-RING | **CROSS-CHECKED** | max \|dCCE\| fv2d vs stored devsim, guard ring (ALT_STACK, legacy mesh), 10 V, 11-point anode scan | 0 | see results | \|E\| <= 1e-10 (abs) | slow |
| CA-PHOTONS-PLATEAU | **CALIBRATED** | rms residual of eps*A/C vs the digitised Photons Fig. 4b plateau (-5..-59 V, 208 points) [um] | 0 (+-0.02) | 0.01326 | regression +-2e-03 (abs) | fast |
| CA-W-LABEL-0 | **CALIBRATED** | eps*A/C at 0 V vs the printed Fig. 4b label [um] | 1.68 (+-0.005) | -0.01504 | regression +-5e-04 (rel) | fast |
| CA-W-LABEL-10 | **CALIBRATED** | eps*A/C at -10 V vs the printed Fig. 4b label [um] | 9.52 (+-0.005) | -0.000155 | regression +-5e-04 (rel) | fast |
| CA-W-LABEL-30 | **CALIBRATED** | eps*A/C at -30 V vs the printed Fig. 4b label [um] | 9.73 (+-0.005) | 0.002085 | regression +-5e-04 (rel) | fast |
| CA-DARK-60 | **CALIBRATED** | dark current at -60 V, 300 K, 4 mm2, default stack [A] | 1.8e-11 | -6.3e-07 | regression +-1e-04 (rel) | fast |
| CA-RADMANOVAC-ONSEMI | **CALIBRATED** | C-V removal slope g (per p/cm2) at fixed width 10-40 um, onsemi 50 um 4e13, 252.7 MeV protons | 5.6 (+-0.81) | -0.1327 | regression +-5e-03 (rel) | fast |
| CA-HATAKEYAMA-1E17 | **CALIBRATED** | avalanche V_BD, p+/10 um/n+ at N_D 1e17, hatakeyama2004 set [V] | 307.5 | 0.01767 | regression +-1e-03 (rel) | fast |
| CA-BRUZZI-W40-MATCHED | **CALIBRATED** | eps*A/C at the paper's -40 V band bending (V_bi matched), Bruzzi diode configuration [um] | 4.2 (+-0.05) | 0.01001 | regression +-1e-03 (rel) | fast |
| CA-BRUZZI-W40-RAW | **CALIBRATED** | eps*A/C at -40 V (no V_bi matching), Bruzzi diode configuration [um] | 4.2 (+-0.05) | 0.03103 | regression +-1e-03 (rel) | fast |
| CA-TILTED-TAU-SCALE | **CALIBRATED** | plasma-lifetime scale tau_scale fitted to the Micro Fig. 4b -20..-50 V tail (2.5D tilted track) | computed | S = 0.039 | regression +-0e+00 (abs) | record |
| VA-NIWA-3E16 | **VALIDATED** | avalanche V_BD, p+/10 um/n+ at N_D 3e16, niwa2015 set (opt-in) [V] | 805.6 | -0.06444 | regression +-1e-03 (rel) | fast |
| VA-NIWA-1E17 | **VALIDATED** | avalanche V_BD at N_D 1e17, niwa2015 set (opt-in) [V] | 307.5 | 0.05236 | regression +-1e-03 (rel) | fast |
| VA-GSPONER-1E15 | **VALIDATED** | alpha CCE at 1e15 neq, -300 V, Gsponer-like diode (50 um, uniform 1.5e14), default damage model | 0.5 | 0.4545 | regression +-1e-02 (abs) | fast |
| EF-PHOTONS-KNEE | **EXPECTED-FAIL** | rms residual of eps*A/C vs the 13 digitised knee markers (-3.25..-0.25 V) [um] | 0 (+-0.03) | 0.2365 | regression +-1e-02 (abs) | fast |
| EF-TILTED-ALPHA-RMS | **EXPECTED-FAIL** | rms shape residual, microdosimeter tilted-alpha CCE(V) normalised at -50 V vs Micro Fig. 4b (7 points) | 0 (+-0.0005) | 0.03538 | regression +-2e-03 (abs) | fast |
| EF-TILTED-ALPHA-0V | **EXPECTED-FAIL** | normalised tilted-alpha CCE at 0 V | 0.9083 (+-0.0005) | 0.08449 | regression +-3e-03 (abs) | fast |
| EF-TILTED-TRACK-ALL7 | **EXPECTED-FAIL** | chi2 of the best single tau_scale over all 7 Micro Fig. 4b points (2.5D tilted track) | 0 | 371 | regression +-0e+00 (abs) | record |
| EF-DENAPOLI-LOSS-729 | **EXPECTED-FAIL** | relative alpha-CCE loss at 7.29e14 n/cm2 (14 MeV neutrons), 20 um uniform 1e14, -300 V, default damage model | 0.225 | -0.9064 | regression +-2e-02 (rel) | fast |
| EF-DENAPOLI-LOSS-131 | **EXPECTED-FAIL** | relative alpha-CCE loss at 1.31e14 n/cm2 (14 MeV neutrons), same device and model | 0.07 | -0.9375 | regression +-2e-02 (rel) | fast |

<!-- END benchmark-registry -->

The validation uncertainty u_val of V&V 20 is **not available** for the three VALIDATED
rows:

- niwa2015's D is Hatakeyama's power-law fit to measured V_BD, with no stated scatter;
- Gsponer's D is a bound.

E is therefore a model-error estimate without an interval.

How to read the non-verification rows:

- **Photons Fig. 4b** (CA-PHOTONS-PLATEAU, CA-W-LABEL-*, EF-PHOTONS-KNEE). The default
  doping is an inversion of this one curve, so plateau and label agreement is fit quality.
  The knee misses by 0.24 µm rms (13 markers, max 3.9σ): the exponential epi cannot make
  the 1.3e13 dip at 4-6 µm ([lab/buffer_profile/notes/LOG.md](../lab/buffer_profile/notes/LOG.md) §3).
- **Bruzzi & Verroi 2023** (CA-BRUZZI-*). The N_D fed to ETNA is the paper's own fit of
  the C-V from which W comes, so the +1.0 % (V_bi matched) and +3.1 % (raw) are a
  solver-scale consistency on a different device (a Schottky diode), not a validation.
  [.planning/VALIDATION_LITERATURE.md](../.planning/VALIDATION_LITERATURE.md) has the
  full argument.
- **Breakdown** (CA-HATAKEYAMA-1E17, VA-NIWA-*). Hatakeyama's coefficients were fitted on
  the same diodes, so the +1.8 % at 1e17 is consistency. niwa2015 was not fitted to them:
  −6.4 % at 3e16 and +5.2 % at 1e17 are model errors on those diodes. Two caveats:
  - Its coefficients are secondary-source values (the Burin review).
  - The default dosimeter's V_BD (1826-2226 V across the sets) is an extrapolation in
    doping and an upper bound with no edge termination
    ([USAGE §18](USAGE.md#validation-and-calibration-labels)).
- **Radiation damage** (CA-RADMANOVAC-ONSEMI, VA-GSPONER-1E15, EF-DENAPOLI-LOSS-*).
  - The removal slope is fitted to Radmanovac Table 2 (other devices). The −13 % on the
    onsemi slope is inside its ±0.81 data spread (±14 %).
  - Gsponer: the default model gives CCE 0.955 at 1e15 neq and −300 V. This is
    consistent with the abstract's bound "CCE better than 50 % up to 1e15 n_eq/cm²"
    (arXiv:2310.02047, abstract opened 2026-10-21; the abstract states no bias). A bound
    is a weak comparison.
  - **De Napoli** (Front. Phys. 10 (2022) 898833, opened 2026-10-21): 14 MeV neutrons,
    20 µm / 1e14 Schottky diodes, ²³⁹Pu alphas at 300 V. The measured CCE losses are
    about 7.0 % and 22.5 % at 1.31e14 and 7.29e14 n/cm². The default model (NJOY SiC
    NIEL, `particle="neutron"`) gives 0.44 % and 2.1 %, about 10x less loss.
  - Before the 2026-10-21 default flip, the lab bench (si-neq convention, steady traps,
    fluence read directly as neq) gave about 10 % at 7.3e14. That is the basis of
    USAGE's "same decade". **That statement no longer holds for the default model**
    (section 6).
- **Microdosimeter tilted alpha** (EF-TILTED-ALPHA-*, CA-TILTED-TAU-SCALE,
  EF-TILTED-TRACK-ALL7). The steady drift-diffusion shape misses the data by rms 0.035
  against σ 0.0005: model 0.993 vs data 0.908 at 0 V.
  - The 2.5D track model fits the −20..−50 V tail only with a fitted plasma lifetime of
    about 23 ns (about 100x below good-epi values). That fit is CALIBRATED.
  - It misses all 7 points (χ² 371 for 5 dof).
  - The track-core radius is not constrained by the data
    ([USAGE §19.3](USAGE.md#193-fig-4b-with-the-tilted-model-a-calibration-not-a-reproduction)).

---

## 6. Known limitations and open items

Physics limits are listed in [USAGE §8 "Not modelled"](USAGE.md#8-validity-and-limits)
and are not repeated here. Open items for credibility:

1. **No validation on the modelled devices.** The only device-level comparison that was
   not fitted fails (microdosimeter CCE(V)). A CCE-vs-bias measurement on the dosimeter,
   SIMS, and C-V to −100..−200 V would change this (section 7).
2. **Damage CCE vs De Napoli** (new, from this registry): under the defaults of
   2026-10-21 the modelled 14 MeV-neutron loss is about 10x below the measurement, where
   the pre-flip model was within 2x. Possible causes, none tested here:
   - the neutron neq normalisation;
   - the calibrated slow-removal component;
   - trapping not represented by the steady CCE;
   - Schottky vs p-n geometry.
     This belongs to issues #5 / #8. The "same decade" sentence in USAGE (radiation-traps
     section) describes the steady opt-in model with si-neq fluences; it now says so.
3. **Gsponer bound quoted as "above 80 %".** `lab/radiation_traps/notes/LOG.md` §3 (and
   USAGE until this report) quote "above 80 %". The abstract states "better than 50 %" up
   to 1e15 n_eq/cm², and `MECHANISMS.md` agrees. USAGE now quotes the abstract. The 80 %
   figure stays unverified unless a source in the paper body is cited.
4. **MMS / order-of-accuracy test for the 1D Poisson and DD solvers** is missing (issue #3
   P2.1 (i)). The 1D production mesh is first order (section 4.2).
5. **No propagated uncertainty (UQ)** on any output. The input spreads are documented
   (section 2) but not pushed through.
6. **Bragg additivity for SiC stopping power** is unvalidated (no SiC measurement in the
   repository). κ(E) inherits it.
7. **Not re-run here** (recorded entries): the tilted-track calibration
   (`lab/track_plasma/bench/tilted_bench.py`, minutes) and the long devsim 2D baselines.
8. **Platform coverage.** The goldens were measured on macOS arm64. CI runs ubuntu only.
   The verification gates are the existing test tolerances, not the measured values.
   The regression guards on E (1e-3 relative and looser) leave room for BLAS
   differences, but this has not been measured on Linux yet (issue #3 P3.1).
9. **Benchmark CI job.** The fast subset runs inside the default suite
   (`tests/test_benchmarks.py`, `-m benchmark`); there is no separate CI job or
   cross-platform matrix yet.

---

## 7. Data that would raise credibility

The full list, with templates, formats and loader checks, is
[docs/DATA_REQUEST_CHECKLIST.md](DATA_REQUEST_CHECKLIST.md). In V&V terms, the items that
would turn a CALIBRATED or EXPECTED-FAIL row into a validation are:

| measurement                                         | registry rows it would affect                                                               |
| --------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| SIMS / SRP through the epi and the n⁻/n⁺ transition | CA-PHOTONS-*, EF-PHOTONS-KNEE: the doping becomes an input, so the C-V becomes a validation |
| C-V on a second die size, and to −100..−200 V       | the fringe-capacitance and transition-form assumptions                                      |
| leakage vs bias and temperature                     | CA-DARK-60: shape and activation energy                                                     |
| CCE vs bias on the dosimeter (alpha, proton)        | the first device-level transport validation                                                 |
| C-V and CCE vs fluence on these devices             | CA-RADMANOVAC-ONSEMI, EF-DENAPOLI-*: removal and trapping on the actual epi                 |
| mask data                                           | the layout's `REFERENCE-NOT-FAB` state                                                      |

---

## 8. Reproducing this report

```bash
uv run python -m benchmarks                   # fast subset: benchmarks/results/fast.{json,md}
uv run python -m benchmarks --slow            # + slow entries:  benchmarks/results/all.{json,md}
uv run python -m benchmarks --update-report   # regenerate the table of section 5
uv run pytest -m benchmark -q                 # the fast subset as tests
```

Run time of the fast subset on the 10-core M4 (2026-10-21): 7-15 s with a warm numba
cache. With an empty `NUMBA_CACHE_DIR` it took 162 s at load average ~50, nearly all JIT
compilation of the engines. The `--slow` entries add about 12 s.

Each JSON result records the environment: Python, platform, and the etna, numpy, scipy,
numba and devsim versions, plus the `uv.lock` hash.
