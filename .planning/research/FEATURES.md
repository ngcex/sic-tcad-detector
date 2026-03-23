# Feature Landscape

**Domain:** v1.1 realistic device physics for 4H-SiC TCAD radiation detector simulator
**Researched:** 2026-03-23
**Scope:** NEW features only -- temperature-dependent simulation, realistic dark current, transient FLASH dynamics. Does not re-document v1.0 features (I-V, C-V, CCE, steady-state FLASH, parametric sweeps).

## Table Stakes

Features that a reviewer expects when a paper claims "realistic device simulation" or "temperature-dependent modeling." Missing any of these undermines credibility of the v1.1 claims.

| Feature                                     | Why Expected                                                                                                                                                                                                                    | Complexity  | Notes                                                                                                                                                                                                    |
| ------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **T-dependent bandgap E_g(T)**              | Varshni model is universal in TCAD. Any paper mentioning temperature must use it. Already partially coded in `compute_ni()` but not wired into device pipeline.                                                                 | Low         | Varshni params already in `SiC4H_Parameters`: E_g_0=3.265, alpha=6.5e-4, beta=1300. Just needs wiring.                                                                                                   |
| **T-dependent n_i(T)**                      | Intrinsic carrier concentration drives all junction physics. n_i(300K)~5e-9 but n_i(310K) is ~3x larger due to exponential Eg/2kT sensitivity. Must be self-consistent with E_g(T).                                             | Low         | `compute_ni(T)` already exists but returns a float not wired into devsim. Need to set as device parameter at simulation start.                                                                           |
| **T-dependent mobility mu(T)**              | Caughey-Thomas with T-exponents is the standard TCAD model. Phonon scattering increases with T, reducing mobility. Required for I-V and CCE accuracy at T != 300K.                                                              | Low-Medium  | mu_n_max(T) = 950*(T/300)^(-2.40), mu_p_max(T) = 125*(T/300)^(-2.15). TU Wien Ayalew Table 3.5. Need to parameterize existing Caughey-Thomas with T.                                                     |
| **T-dependent NC(T), NV(T)**                | Density of states scales as T^(3/2). Feeds n_i(T) computation. Standard in every TCAD tool.                                                                                                                                     | Low         | NC(T) = NC_300*(T/300)^1.5, NV(T) = NV_300*(T/300)^1.5. Trivial.                                                                                                                                         |
| **Dark current decomposition**              | Any paper claiming to match experimental I_dark must decompose contributions: bulk SRH generation, surface recombination, trap-assisted tunneling. Reviewers expect this breakdown to validate the dominant mechanism.          | Medium      | Current I_dark = 6.71e-49 A (bulk SRH only) vs experimental 18 pA. The ~40 orders of magnitude gap proves surface/TAT physics is essential, not optional.                                                |
| **Reverse I-V matching experiment**         | The Petringa photons paper reports I_dark < 18 pA up to -60V. Reproducing this with physics-based models (not just fitting) validates the entire device model.                                                                  | Medium-High | Requires surface recombination + Hurkx TAT working together. The 18 pA target provides a clear pass/fail criterion.                                                                                      |
| **Temperature coefficient of dark current** | Experimental value: (-0.079 +/- 0.005)%/C for SiC dosimeter sensitivity. A T-dependent model must reproduce the correct sign and approximate magnitude of I_dark(T) variation.                                                  | Medium      | Dominated by n_i(T) exponential increase but offset by mobility decrease. Net effect should match experimental trends from Lopez Paz 2024.                                                               |
| **Transient current pulse I(t)**            | Any paper claiming transient simulation must show the current waveform vs time. Reviewers expect to see carrier drift, diffusion tail, and collection time. The TCT (Transient Current Technique) is the experimental analogue. | Medium      | devsim transient solver (BDF1/BDF2) already wired. Need to implement time-step loop, extract current at each step, plot I(t). Time resolution ~1 ns for SiC (transit time ~50 ps across 10 um at v_sat). |
| **Intra-pulse carrier dynamics**            | During a FLASH pulse (10-200 ms), carriers build up toward steady state. Must show time-resolved n(x,t), p(x,t) during pulse. Without this, "transient" simulation is not transient.                                            | Medium-High | Requires stable time-stepping over 6+ orders of magnitude (ns to ms). Adaptive dt essential. devsim charge_error parameter controls this.                                                                |
| **Validation against analytical limits**    | Low dose-rate transient must recover steady-state CCE. Zero-field limit must show full recombination. High-field limit must show full collection. These limiting cases validate the transient solver independently.             | Low-Medium  | Plot CCE_transient(dose_rate -> 0) vs CCE_steady_state. Should agree to <1%. This is a necessary sanity check.                                                                                           |

## Differentiators

Features that add scientific novelty beyond standard TCAD. These are what makes v1.1 publishable rather than merely incremental.

| Feature                                              | Value Proposition                                                                                                                                                                                                                                                                      | Complexity  | Notes                                                                                                                                                                                                                                             |
| ---------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **T-dependent SRH lifetimes tau(T)**                 | SRH lifetime in 4H-SiC is governed by the Z1/2 center (carbon vacancy at EC-0.65 eV). tau(T) is relatively flat below ~500K then thermally activated above. For the clinical range (303-313K), the effect is small but measurable. Including this is beyond standard 300K-only TCAD.   | Medium      | Model: tau(T) = tau_0 * [1 + C*exp(-E_a/kT)] where E_a ~ 0.1-0.3 eV relates to Z1/2 emission. tau_0 needs calibration against dark current. IEEE Access 2024 (10538275) reports power-law T dependence.                                           |
| **Hurkx trap-assisted tunneling model**              | TAT via midgap traps is the dominant reverse leakage mechanism in SiC at room temperature (not band-to-band tunneling -- Eg too large). Implementing Hurkx as a field-enhanced SRH correction is the standard TCAD approach. Differentiates from ideal-SRH-only models.                | Medium-High | Klaassen approximation for Gamma(F). Key params: E_t = E_g/2 (midgap), m_t = 0.25*m_0 (from arXiv:2503.09016). F_crit = sqrt(24*m_t*E_t^3)/(q*hbar). Implemented as node model modifying existing USRH.                                           |
| **Surface recombination velocity (SRV) at contacts** | SiO2-passivated 4H-SiC surfaces have SRV = 150-5000 cm/s depending on preparation. This is likely the dominant contributor to the 18 pA dark current, not bulk SRH. Implementing SRV as a contact boundary condition is non-trivial in devsim.                                         | Medium-High | Replace Ohmic carrier pinning at anode with SRH surface recombination flux: R_s = S*(n*p - n_i^2)/(n + p + 2\*n_i). S is the primary fitting parameter. Kimoto et al. JAP 127 (2020) provides reference values.                                   |
| **Inter-pulse memory effects**                       | After a FLASH pulse ends, residual carriers decay via SRH/Auger over ~100 ns to ~10 us. If the next pulse arrives before full decay, the starting carrier density is elevated, reducing CCE for subsequent pulses. No TCAD study has quantified this for SiC.                          | High        | Requires multi-pulse simulation: pulse ON -> pulse OFF -> (inter-pulse gap) -> pulse ON. Sweep inter-pulse gap from 1 us to 100 ms. Plot CCE(pulse_number) to show memory effect. Key FLASH physics insight.                                      |
| **Plasma build-up time constant**                    | During a pulse, carriers accumulate from zero toward a steady-state density set by generation-recombination balance. The characteristic time to reach steady state depends on tau_SRH, dose rate, and electric field. Extracting this timescale is novel.                              | Medium      | From transient n(t) during pulse, fit n(t) = n_ss \* (1 - exp(-t/tau_buildup)). Plot tau_buildup vs dose rate. Compare to tau_SRH. At low dose rates, tau_buildup ~ tau_SRH. At high dose rates (Auger-limited), tau_buildup shortens.            |
| **Carrier sweep-out dynamics**                       | After pulse OFF, the E-field sweeps collected carriers out while residual carriers recombine. The competition between drift sweep-out and recombination determines the "afterglow" current. Time-resolved I(t) after pulse-off is experimentally measurable via fast electronics.      | Medium      | Plot I(t) after pulse-off. Extract sweep-out time constant. Compare to transit time d/v_sat ~ 50 ps and tau_SRH. The ratio determines whether the detector is "transit-time limited" or "recombination limited".                                  |
| **CCE per pulse in multi-pulse train**               | FLASH delivers dose in multiple pulses. Real CCE is per-pulse, not time-averaged. If inter-pulse memory is significant, CCE degrades pulse-by-pulse until reaching a new steady state. This pulse-resolved CCE is what the experimental SiC dosimeter actually measures.               | High        | Simulate N=10+ pulse train. Extract CCE for each pulse. Plot CCE(n) for pulse number n. This directly maps to the pulse-resolved measurement configuration used by Lopez Paz et al. 2024 (linear up to 11 Gy/pulse, <3% deviation up to 4 MGy/s). |
| **T-dependent CCE**                                  | Temperature affects mobility (reduces drift velocity), lifetimes (changes recombination), and n_i (changes equilibrium). Net effect on CCE at clinical temperatures (30-40C) is small but experimentally relevant for dosimetric accuracy. No published SiC TCAD study combines T+CCE. | Medium      | Run CCE vs bias at T = 300, 303, 306, 310, 313 K. The (-0.079%/C) temperature coefficient from experiment should emerge from the simulation naturally if physics is correct.                                                                      |

## Anti-Features

Features to explicitly NOT build in v1.1. These are tempting scope expansions that would delay delivery without adding to the core scientific contribution.

| Anti-Feature                                    | Why Avoid                                                                                                                                                                                                                               | What to Do Instead                                                                                                                |
| ----------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| **Coupled thermal simulation (self-heating)**   | Clinical T range is 30-40C; self-heating from radiation is negligible (< 0.01 K). SiC thermal conductivity (3-5 W/cm-K) makes thermal runaway impossible. Coupling heat equation to DD adds solver complexity for zero physics insight. | Parameterize T as a fixed input. Run separate simulations at discrete T values.                                                   |
| **Multi-level trap model**                      | Some SiC papers model Z1/2, EH6/7, and additional deep levels separately. For the 18 pA dark current match, a single effective midgap trap captures 90%+ of the physics. Multi-trap adds fitting parameters without unique insight.     | Single midgap trap (E_t = E_g/2) for Hurkx TAT. Effective tau_SRH calibrated to experiment. Note multi-trap as future refinement. |
| **Band-to-band tunneling (BTBT)**               | BTBT requires E-fields > 3 MV/cm to be significant in 4H-SiC (Eg = 3.26 eV). Max field in this device at -60V is ~0.1 MV/cm. BTBT is irrelevant.                                                                                        | Do not implement. Note in paper that BTBT is negligible for this field range.                                                     |
| **Non-local tunneling model**                   | Local Hurkx is sufficient for uniform junctions. Non-local models (WKB path integrals) are needed for heterostructures or very thin barriers. This is a simple p-n junction -- local model is appropriate.                              | Use local Hurkx (node model). Note non-local as unnecessary for this geometry.                                                    |
| **Noise simulation (1/f, shot, thermal)**       | Noise is small-signal; FLASH is large-signal. Dark current noise floor is irrelevant to CCE degradation at MGy/s dose rates.                                                                                                            | Report dark current value. Do not simulate noise spectrum.                                                                        |
| **2D/3D field effects during pulse**            | Edge effects and non-uniform carrier generation add 2D complexity. The 1D model captures the dominant physics (field along depth direction). Edge effects contribute ~2% (build-up over-response).                                      | Stay 1D. Note edge effects as v2 item.                                                                                            |
| **Radiation damage evolution**                  | Cumulative damage (sensitivity drift after 100 kGy) is a separate physics problem from instantaneous pulse dynamics. CERN RD50 groups are working on this.                                                                              | Assume pristine device. Cite Lopez Paz 2024: sensitivity reduction < 2% after 100 kGy.                                            |
| **Readout electronics simulation (RC circuit)** | Experimental SiC systems use RC shaping circuits. Simulating the electronics chain is outside device physics scope.                                                                                                                     | Report raw detector current I(t). Leave electronics integration to experimentalists.                                              |

## Feature Dependencies

```
Existing v1.0 Foundation
  |
  +---> Poisson + DD solver (working)
  |     |
  |     +---> SRH recombination (working)
  |     |
  |     +---> Auger recombination (working)
  |     |
  |     +---> Steady-state FLASH (working)
  |
  +---> sic_material.py params (300K only)
  |
  v
v1.1 Feature Groups (three semi-independent branches)

[Branch 1: Temperature Dependence]
  sic_material.py: T-dependent functions
    |
    +---> E_g(T) via Varshni  ----+
    |                              |
    +---> NC(T), NV(T) ~ T^1.5 --+--> n_i(T) = sqrt(NC*NV)*exp(-Eg/2kT)
    |                              |
    +---> mu_n(T), mu_p(T)         |
    |     Caughey-Thomas with      |
    |     T-exponents              |
    |                              |
    +---> tau_SRH(T)               |
          Z1/2 center model        |
                                   |
    device.py: wire T into --------+
    create_sic_device()
      |
      v
    T-dependent I-V, C-V, CCE
      |
      v
    Temperature coefficient validation (-0.079%/C)

[Branch 2: Realistic Dark Current]
  Hurkx TAT model (new module)
    |
    +---> Field enhancement Gamma(F)
    |     Klaassen approximation
    |     REQUIRES: ElectricField from Poisson (existing)
    |
    +---> Modified USRH with Gamma
          REQUIRES: existing SRH model in drift_diffusion.py
    |
    v
  Surface recombination (new module)
    |
    +---> SRV contact equation
    |     REQUIRES: custom contact_equation replacing Ohmic pinning
    |     REQUIRES: n_i(T) from Branch 1
    |
    v
  Dark current orchestrator
    |
    +---> I_dark = I_bulk_SRH + I_surface + I_TAT
    |     REQUIRES: T-dependent n_i from Branch 1
    |     FIT TARGET: 18 pA at -30V, 300K
    |
    v
  Validation: I_dark(V) and I_dark(T)

[Branch 3: Transient FLASH Dynamics]
  Transient solver loop (new module)
    |
    +---> Pulse ON: generation rate active
    |     REQUIRES: existing gen profiles from generation_profiles.py
    |     REQUIRES: devsim transient_bdf1/bdf2 (already wired)
    |
    +---> Pulse OFF: generation zeroed, carrier decay
    |
    +---> I(t) extraction at each time step
    |
    v
  Single-pulse analysis
    |
    +---> Intra-pulse carrier build-up n(t)
    +---> Post-pulse sweep-out I(t)
    +---> CCE = integral(I(t)) / Q_generated
    |
    v
  Multi-pulse train
    |
    +---> Inter-pulse memory: CCE(pulse_number)
    |     REQUIRES: single-pulse working
    |
    +---> Steady-state CCE vs dose rate (validates against v1.0)
    |     REQUIRES: low-dose-rate limit matches existing steady-state result
    |
    v
  Parametric studies
    |
    +---> CCE_transient vs dose rate
    +---> CCE_transient vs T (connects Branch 1 + Branch 3)
    +---> Pulse structure optimization (gap, duration, DPP)
```

**Key dependency insight:** Branches 1 and 2 share a critical coupling: surface recombination and TAT both depend on n_i(T). Branch 3 is mostly independent but needs Branch 1 for T-dependent transient studies. The v1.0 steady-state FLASH result serves as the validation anchor for Branch 3 (transient CCE at long times must converge to steady-state CCE).

**Cross-branch validation:** The CCE temperature coefficient emerges from Branch 1 + Branch 3 combined: run transient FLASH at multiple temperatures and extract sensitivity(T). This is the highest-value combined result.

## MVP Recommendation

**Prioritize (Phase 1 -- Temperature Foundation + Dark Current Physics):**

1. T-dependent material parameters: E_g(T), n_i(T), mu(T), NC/NV(T) in `sic_material.py`
2. Wire T into `create_sic_device()` and `setup_sic_drift_diffusion()`
3. T-dependent I-V validation (I_dark should increase with T)
4. Hurkx TAT node model (field-enhanced SRH)
5. Surface recombination contact equation with SRV
6. Calibrate {tau_SRH, S, m_t} to match 18 pA dark current
7. Dark current decomposition plot: I_bulk + I_surface + I_TAT vs V

**Prioritize (Phase 2 -- Transient FLASH Dynamics):**

1. Single-pulse transient: I(t) during and after pulse
2. CCE from time-integrated current
3. Validate: transient CCE at low dose rate matches v1.0 steady-state
4. Multi-pulse train: CCE(pulse_number) for N=10 pulses
5. Inter-pulse gap sweep: CCE vs gap duration
6. Carrier build-up/decay time constants

**Prioritize (Phase 3 -- Combined Analysis + Publication Figures):**

1. T-dependent CCE: CCE vs bias at T = 300-313K
2. Temperature coefficient extraction and validation (-0.079%/C)
3. Combined parametric: CCE vs {dose_rate, T, bias, pulse_structure}
4. Publication-quality transient waveform figures
5. Comparison table: transient CCE vs steady-state CCE across parameter space

**Defer to v2:**

- Build-up over-response (2D effects, lower novelty)
- Azimuthal response (needs 2D/3D mesh)
- Multi-trap dark current model (diminishing returns over single midgap)
- Coupled thermal simulation (T range too narrow to matter)

## Complexity Budget

| Phase                      | Features                                                        | Estimated Complexity                                                                                   | Confidence |
| -------------------------- | --------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------ | ---------- |
| Temperature + Dark Current | T-dependent params, Hurkx TAT, SRV, calibration                 | Medium-High -- new devsim equation patterns (contact_equation, modified USRH), fitting loop            | MEDIUM     |
| Transient FLASH            | Single-pulse I(t), multi-pulse train, inter-pulse memory        | High -- time-stepping stability over 6 orders of magnitude, adaptive dt, convergence at high injection | MEDIUM-LOW |
| Combined Analysis          | T-dependent transients, parametric studies, publication figures | Medium -- uses existing pieces, mainly orchestration and visualization                                 | HIGH       |

**Critical risks:**

1. **Dark current fitting (Phase 1):** Three free parameters (tau_SRH, S, m_t) with one target (18 pA). Risk of non-unique fit. Mitigate by constraining each to literature ranges and decomposing I-V shape (SRV dominates at low bias, TAT at high bias).

2. **Transient convergence (Phase 2):** devsim transient solver at high carrier injection (FLASH dose rates) may require very small initial time steps and careful ramp-up. The v1.0 steady-state FLASH study already showed convergence challenges at high generation rates. Transient adds time-stepping instability risk.

3. **Inter-pulse timescale gap (Phase 2):** Pulse duration (~10-200 ms) vs carrier dynamics (~ns-us) creates a 6-order timescale gap. Must use adaptive time-stepping to avoid either (a) millions of small steps or (b) missing the physics. The TR-BDF2 composite method in devsim handles stiff problems but needs careful gamma parameter tuning.

## Validation Targets

| Feature                      | Validation Method                                          | Pass Criterion                                     | Source                   |
| ---------------------------- | ---------------------------------------------------------- | -------------------------------------------------- | ------------------------ |
| T-dependent n_i              | Compare to Ioffe NSM Archive values                        | Within 10% at 300-400K                             | Ioffe NSM                |
| T-dependent mobility         | Compare to TU Wien Ayalew Table 3.5                        | Within 5% at 300K, correct T-exponent sign         | Ayalew thesis            |
| Dark current magnitude       | Match Petringa experimental value                          | I_dark = 18 pA (+/- factor 2) at -30V, 300K        | Photons paper p.6        |
| Dark current vs voltage      | I_dark(V) shape matches experiment                         | Monotonic increase, plateau above full depletion   | Photons paper Fig. 6     |
| T coefficient of sensitivity | Compare to Lopez Paz 2024                                  | (-0.079 +/- 0.005)%/C                              | Lopez Paz, Med Phys 2024 |
| Transient CCE low-dose limit | Converge to v1.0 steady-state CCE                          | Agreement within 1%                                | Internal consistency     |
| Transient CCE high-dose      | CCE flat at ~1.0 across 20-230 Gy/s (Auger negligible)     | Consistent with v1.0 null result                   | v1.0 finding             |
| SiC dose-rate linearity      | Linear up to 11 Gy/pulse, <3% deviation up to 4 MGy/s      | Consistent with Lopez Paz 2024 experimental        | Lopez Paz, Med Phys 2024 |
| Pulse-resolved operation     | Detector tracks pulse structure including irregular pulses | Qualitative agreement with experimental capability | Lopez Paz 2024           |

## Sources

- Petringa group papers (SiC_Photons_MedicalPhysics, Microdosimetry.pdf, Flash.pdf) -- experimental validation targets (HIGH confidence)
- [TU Wien Ayalew thesis, Table 3.5](https://www.iue.tuwien.ac.at/phd/ayalew/node65.html) -- Caughey-Thomas T-dependent mobility parameters for 4H-SiC (HIGH confidence)
- [Silicon Carbide Sensors in Radiotherapy Dosimetry: Review (Frontiers 2025)](https://www.frontiersin.org/journals/sensors/articles/10.3389/fsens.2025.1622153/full) -- SiC detector FLASH challenges, thermal stability (HIGH confidence)
- [State-of-the-art SiC diode dosimeters for UHDR (Lopez Paz 2024)](https://pubmed.ncbi.nlm.nih.gov/38530300/) -- 11 Gy/pulse linearity, 4 MGy/s, -0.079%/C temperature coefficient, <2% sensitivity drift after 100 kGy (HIGH confidence)
- [Carrier Lifetime vs Temperature in 4H-SiC (IEEE Access 2024)](https://ieeexplore.ieee.org/document/10538275) -- power-law tau(T) dependence, Arrhenius damage coefficient (MEDIUM confidence)
- [Hurkx TAT model original paper (IEEE TED 1992)](https://www.semanticscholar.org/paper/A-new-recombination-model-for-device-simulation-Hurkx-Klaassen/4e0ad76a1a7d0e1b4db5f1e48bc05a6f16614337) -- field-enhanced SRH formulation (HIGH confidence)
- [4H-SiC TAT with m_t = 0.25\*m_0 (arXiv:2503.09016)](https://arxiv.org/pdf/2503.09016) -- tunneling effective mass for SiC (MEDIUM confidence)
- [Kimoto et al., SRV for 4H-SiC (JAP 2020)](https://pubs.aip.org/aip/jap/article/127/19/195702/153502/) -- surface recombination velocity 150-5000 cm/s (MEDIUM confidence)
- [Multi-Level TAT for SiC-JBS Reverse Leakage (Scientific.Net)](https://www.scientific.net/MSF.924.601) -- Z1/2 and EH6/7 trap levels, field-dependent TAT (MEDIUM confidence)
- [Modified TCAD approach for 4H-SiC JBS diodes (ScienceDirect 2025)](https://www.sciencedirect.com/science/article/abs/pii/S0026271425002896) -- electron trapping at EC-0.19, EC-0.65, EC-1.65 eV levels (MEDIUM confidence)
- [Carrier lifetime modulation in SiC PiN (Discover Nano 2023)](https://link.springer.com/article/10.1186/s11671-023-03905-6) -- Z1/2 center lifetime control, pulsed operation physics (MEDIUM confidence)
- [First Characterization SiC Detectors UHDR Electron Beams (MDPI 2023)](https://www.mdpi.com/2076-3417/13/5/2986) -- SiC FLASH experimental context (MEDIUM confidence)
- [RASER: Simulation of radiation damage in SiC detectors (arXiv 2025)](https://arxiv.org/html/2504.20463) -- transient current simulation methodology, SRH trapping time (MEDIUM confidence)
- [DEVSIM: A TCAD Semiconductor Device Simulator](https://www.researchgate.net/publication/358631134_DEVSIM_A_TCAD_Semiconductor_Device_Simulator) -- transient solver capabilities (HIGH confidence)
- [TCAD Parameters for 4H-SiC: A Review (Burin, CERN RD50)](https://jburin.web.cern.ch/pdfs/4HSiC_review.pdf) -- comprehensive parameter survey (MEDIUM confidence, PDF unreadable)
