# Project Research Summary

**Project:** SiC TCAD Simulator v2.0 — Radiation Damage Modeling
**Domain:** Physics-based TCAD simulation / semiconductor detector science (4H-SiC)
**Researched:** 2026-03-24
**Confidence:** MEDIUM-HIGH

## Executive Summary

This project extends a validated 1D devsim-based 4H-SiC drift-diffusion simulator (v1.1) with fluence-dependent radiation damage physics. The target application is predicting performance degradation of the Petringa group's thin-epi PIN detector under clinical proton beams at INFN-LNS (62 MeV). The state-of-the-art approach is the Burin et al. (2024, CERN RD50) three-defect TCAD model: Z1/2, EH6/7, and EH4 deep levels whose concentrations scale linearly with fluence (N_i = g_i \* Phi). This is not empirical curve fitting — the physical defects are DLTS-measured with known energy levels and capture cross-sections. Critically, the entire v2.0 feature set is achievable with the existing stack (devsim + scipy + numpy + matplotlib); no new packages are needed.

The recommended architecture treats fluence exactly as temperature is treated today: a first-class parameter passed to `create_sic_device()` that modifies material properties before the devsim solver runs. This design means the drift-diffusion solver, CCE computation, and dark current modules require zero changes — they consume devsim parameters regardless of origin. The key new module is `radiation_damage.py`, a pure-Python physics layer with no devsim dependency, which computes degraded lifetime (1/tau = 1/tau_0 + K_tau*Phi), effective doping (N_eff = N_D - c_r*Phi), and fluence-dependent trap densities. A `fluence_sweep.py` module orchestrates per-fluence device creation following the established `temperature_sweep.py` pattern.

The dominant risk is regression: the radiation damage module touches SRH recombination, doping, and trap parameters that the v1.1 validated baseline depends on. A non-additive implementation will silently break C-V (R^2=0.998), CCE, and the 18 pA dark current calibration. A second major risk is specific to the Petringa device's low bulk doping (N_D_bulk = 8.5e13 cm^-3): with a carrier removal rate of ~5 cm^-1, full doping compensation occurs at Phi_crit ~ 1.6e13 p/cm^2 — a surprisingly low fluence that falls within the operational range and can cause Newton solver divergence if not handled explicitly.

## Key Findings

### Recommended Stack

No new Python packages are required. All v2.0 physics is pure modeling code on top of devsim 2.10.0 + scipy + numpy + matplotlib. The only "dependency" is a curated set of physical constants from peer-reviewed literature, hardcoded as Python dataclasses. scipy.integrate.solve_ivp (BDF method) handles the stiff annealing kinetics ODEs and is already in the stack. NIEL scaling values for protons in SiC should be pre-computed from the SR-NIEL web calculator and stored as a hardcoded lookup table — no web API integration is needed.

**Core technologies:**

- `devsim 2.10.0`: drift-diffusion + Poisson solver with custom node models — already validated, supports arbitrary trap levels via CreateNodeModel without solver changes
- `scipy.integrate.solve_ivp` (BDF): annealing kinetics ODE system — already in stack, handles stiff first-order equations
- `numpy`: defect parameter tables, NIEL lookup, fluence sweep arithmetic — already in stack
- Literature constants (Burin et al. 2024): three-defect model parameters (g_int, sigma_e, sigma_h per defect) — the core "data dependency", not a software package

**Key parameter values (HIGH confidence, from Burin et al. 2024):**

- Z1/2: E_C - 0.67 eV, acceptor, g = 5.0 cm^-1, sigma_e = 2.0e-14 cm^2, sigma_h = 3.5e-14 cm^2
- EH6/7: E_C - 1.60 eV, donor, g = 1.6 cm^-1, sigma_e = 9.0e-12 cm^2, sigma_h = 3.8e-14 cm^2
- EH4: E_C - 1.03 eV, acceptor, g = 2.4 cm^-1, sigma_e = 5.0e-13 cm^2, sigma_h = 5.0e-14 cm^2
- Carrier removal rate: c_r = 4.2–6.4 cm^-1 for ~250 MeV protons (NIEL-scale upward for 62 MeV)

### Expected Features

**Must have (table stakes for "radiation damage modeling" claim):**

- Defect introduction model (N_defect = g \* Phi) — foundation for all downstream effects
- Carrier lifetime degradation tau(Phi) — primary CCE loss mechanism; every radiation damage paper plots this
- Carrier removal / effective doping N_eff(Phi) — primary C-V observable; expected by every reviewer
- CCE vs fluence prediction — the headline scientific result
- C-V shift with fluence — direct experimental validation target (C-V flattening indicates full compensation)
- Dark current vs fluence — second validation target (note: SiC current does NOT increase monotonically like Si; may decrease at high fluence due to carrier removal)

**Should have (differentiators for publication):**

- Multi-defect trap model (Z1/2 + EH6/7 + EH4) — most published SiC papers use a single effective trap; three-defect model is what commercial TCAD does and enables accurate simultaneous I-V and C-V matching
- NIEL scaling for proton energy dependence — generalizes tool beyond one calibration point; essential for translating neutron-calibrated constants to 62 MeV Petringa conditions
- Forward I-V degradation — independent damage model validation that falls out from multi-defect model without extra code
- Parametric radiation hardness optimization — sweep epi thickness, doping, bias to find operating conditions maximizing CCE at target fluence; the design guidance the group needs for next-generation detectors

**Defer to v3+:**

- Annealing kinetics — complex, limited quantitative SiC data, low practical relevance (Z1/2 is thermally stable to >1000C; room-temperature annealing is negligible)
- Monte Carlo defect cascade simulation — handled by Geant4 group; NIEL scaling from tabulated values is sufficient
- Polarization / space charge buildup dynamics — non-equilibrium trapping, poorly understood, out of scope
- Surface damage / TID effects — bulk displacement damage dominates for this detector geometry

### Architecture Approach

The core architectural decision is to treat fluence as a first-class parameter parallel to temperature. Fluence flows through `radiation_damage.py` (new, pure physics, no devsim dependency) into `sic_material.py` (minimal backward-compatible additions) into `create_sic_device()` (gains Phi and damage_params parameters), where computed scalar irradiated values are set as devsim parameters. All downstream modules — `drift_diffusion.py`, `dark_current.py`, `charge_collection.py` — remain completely unchanged. A fresh devsim device is created per fluence point following the established `temperature_sweep.py` pattern, avoiding in-place parameter mutation issues.

**Major components:**

1. `radiation_damage.py` (NEW) — pure physics: defect concentrations, degraded lifetime, effective doping, NIEL lookup table, Arrhenius annealing ODE; no devsim imports
2. `sic_material.py` (MODIFIED) — adds optional `Phi` parameter to `srh_lifetime()`; backward-compatible at Phi=0 by construction
3. `device.py` (MODIFIED) — `create_sic_device()` gains `Phi` and `damage_params`; wires irradiated scalar values to devsim; stores radiation state in device_info dict
4. `fluence_sweep.py` (NEW) — orchestrates CCE/dark current/C-V vs fluence sweeps, mirrors `temperature_sweep.py` using UUID device names and try/finally cleanup
5. `annealing.py` (NEW, deferred) — ODE-based annealing kinetics using scipy.integrate.solve_ivp with BDF method

### Critical Pitfalls

1. **Fluence=0 regression failure** — the damage module touches the same SRH/doping/trap physics the v1.1 validated baseline depends on. Any non-additive modification silently breaks C-V (R^2=0.998), CCE=100% at V>-40V, and the 18 pA dark current calibration. Prevention: design the damage API as purely additive — at Phi=0 the function does nothing; lifetime degradation updates existing taun/taup rather than replacing the lifetime model; carrier removal creates a new Donors_eff node model rather than overwriting the graded profile. Build the regression test BEFORE writing any damage code.

2. **Carrier removal reaches full compensation at Phi_crit ~ 1.6e13 p/cm^2** — for the Petringa device (N_D_bulk = 8.5e13 cm^-3, c_r ~ 5 cm^-1), this is within the operational fluence range. At this boundary the Newton solver can diverge, the Hecht equation becomes invalid, and the electric field profile changes qualitatively. Prevention: compute Phi_crit as the first task of the carrier removal phase; use logarithmic fluence steps that straddle this boundary; monitor Newton iteration count.

3. **Wrong damage constants — orders-of-magnitude scatter in literature** — g_Z1/2 values range from 0.44 to 5.0 cm^-1 across studies; carrier removal rates range from 5 to 260 cm^-1 depending on particle type, energy, and measurement technique. Using neutron-calibrated introduction rates for 62 MeV proton predictions introduces systematic error. Prevention: all constants must be explicitly named with provenance in a `DamageParameters` dataclass (e.g., `g_Z12: float = 5.0  # Burin 2024, neutron`); use directly measured c_r from C-V (not derived from g values); implement sensitivity analysis varying each constant by 2x.

4. **Dark current double-counting** — the v1.1 Hurkx TAT model uses an effective trap density N_t = 2.2e13 cm^-3 calibrated to 18 pA. This already implicitly includes generation from pre-existing Z1/2 centers. Adding radiation-induced Z1/2 as a separate SRH term on top will double-count at Phi=0. Prevention: use the additive delta model: J_dark(Phi) = J_dark(0) + delta_J_dark(Phi), where delta_J represents only the additional generation from radiation-introduced defects beyond the calibrated baseline.

5. **Linear vs logarithmic lifetime degradation — disputed for SiC** — Luo et al. (2025) found a logarithmic dependence (1/tau = a\*ln(Phi) + b) rather than the standard linear model used in Si. The logarithmic model predicts significantly less CCE degradation at high fluence. Neither model is directly validated for 62 MeV proton conditions. Prevention: implement both models behind a `lifetime_model: "linear" | "logarithmic"` flag; show both predictions in CCE vs fluence plots.

## Implications for Roadmap

The strict physics dependency chain (defect introduction feeds lifetime and doping, which feed CCE and C-V) and the regression risk together dictate the phase structure. The architecture's clean separation — `radiation_damage.py` has no devsim dependency — enables early unit testing before any solver integration.

### Phase 1: Damage Physics Foundation

**Rationale:** Pure-Python physics with no devsim dependency; fully testable as standalone unit tests. Must exist before any device simulation can use irradiated parameters. Writing the regression test first enforces the additive API constraint and catches integration issues before they reach the solver.
**Delivers:** `radiation_damage.py` with `RadiationDamageParams` dataclass (provenance-tagged constants from Burin 2024), `defect_concentration()`, `degraded_lifetime()` (linear + logarithmic models behind flag), `effective_doping()` with floor at zero, NIEL lookup table stub. Regression test suite confirming Phi=0 returns bit-identical v1.1 results for C-V, CCE, dark current.
**Addresses:** Defect introduction model, carrier lifetime degradation (both table stakes)
**Avoids:** Pitfall 1 (regression — build test first), Pitfall 2 (wrong constants — provenance-tagged dataclass), Pitfall 7 (g vs c confusion — separate named variables), Pitfall 3 (SiC != Si — use Burin not Hamburg model)

### Phase 2: Device Integration

**Rationale:** Connects radiation physics to the devsim solver. Requires Phase 1. This is the only phase where new code touches the validated devsim device creation path — treat it with the same care as the regression phase.
**Delivers:** `create_sic_device()` extended with `Phi` and `damage_params` parameters (backward-compatible defaults). Irradiated tau_n, tau_p, N_eff, N_t wired to devsim parameters. device_info dict extended with radiation state. Phi_crit computed and logged for the Petringa device geometry. Carrier removal implemented as a node model (position-dependent) applying N_eff(x) = N_D(x) - c_r\*Phi to the graded profile.
**Uses:** radiation_damage.py from Phase 1, graded doping node model (Approach A: uniform carrier removal across epi)
**Implements:** Parameter passthrough architecture (fluence-as-temperature pattern), fresh-device-per-sweep-point (UUID device names + try/finally cleanup)
**Avoids:** Pitfall 4 (operating point shift — must handle Phi_crit), Pitfall 9 (position-dependence — use node model not scalar), Pitfall 1 (additive API — Phi=0 must change nothing)

### Phase 3: CCE vs Fluence

**Rationale:** The primary scientific deliverable. Requires Phases 1-2. This is the headline result that justifies the v2.0 milestone and enables comparison with the Hecht analytical model.
**Delivers:** `fluence_sweep.py` with `sweep_cce_vs_fluence()`. Notebook: CCE vs Phi at multiple bias voltages. Notebook: CCE vs bias at multiple fluence levels. Hecht equation cross-check with fluence-dependent parameters (DD solver is primary; Hecht is cross-check only). Logarithmic vs linear lifetime model comparison.
**Avoids:** Pitfall 13 (Hecht breakdown at high damage — rely on DD not Hecht beyond Phi_crit), Pitfall 10 (wrong comparison data — document device/energy mismatch when comparing to literature)

### Phase 4: Dark Current vs Fluence

**Rationale:** Second key observable. Requires Phases 1-2. Validates the damage model independently from CCE. Note the counterintuitive SiC behavior (dark current does not increase monotonically with fluence as in Si) must be reproduced.
**Delivers:** `sweep_dark_current_vs_fluence()`. Additive delta_J model preserving the v1.1 18 pA calibration at Phi=0. Notebook: I_dark vs Phi and I-V curves at multiple fluences.
**Avoids:** Pitfall 5 (double-counting — additive delta_J_dark, not merged trap model)

### Phase 5: Carrier Removal / C-V Shift + NIEL Scaling

**Rationale:** C-V is the most directly comparable experimental observable for validation against literature data. NIEL scaling is a prerequisite before claiming 62 MeV predictions since all primary literature data is for different particle energies.
**Delivers:** `sweep_cv_vs_fluence()`. N_eff(x, Phi) node model. NIEL lookup table populated from SR-NIEL web calculator (manual task: proton values at 30, 62, 70, 150 MeV for SiC target). Notebook: C-V curves shifting with fluence, N_D_eff vs Phi. Hardness factor k(62 MeV) documented with uncertainty.
**Avoids:** Pitfall 10 (wrong comparison data — note Luo et al. data is closest match by doping but has 10x different epi thickness)

### Phase 6: Multi-Defect Model + Parametric Optimization

**Rationale:** Replaces the single-effective-defect approximation with the full Burin three-defect TCAD model. Required for publication-quality simultaneous I-V and C-V matching. Parametric optimization synthesizes all damage features into concrete design guidance. Deliberately deferred until the single-defect pipeline is fully validated end-to-end.
**Delivers:** Full multi-level SRH trap model in devsim (Z1/2 + EH6/7 + EH4 as separate CreateNodeModel expressions). Forward I-V degradation vs fluence. Parametric sweep over (epi thickness, N_D, V_bias) at target fluence to find CCE-maximizing geometry. Publication-quality combined figures.
**Avoids:** Pitfall 6 (Newton convergence — start with single defect first, add levels incrementally with iteration monitoring), Pitfall 3 (EH6/7 donor assignment — validate C-V shift direction against literature)

### Phase 7: Annealing Kinetics (Optional / Deferred)

**Rationale:** Deferred because: (1) Z1/2 does not anneal at room temperature (activation energy ~4 eV, stable to >1000C), making this irrelevant for in-service operation; (2) quantitative SiC annealing kinetics data is sparse (LOW-MEDIUM confidence); (3) all other damage features are more scientifically urgent. Build only if the Petringa group explicitly requests thermal recovery predictions.
**Delivers:** `annealing.py` with `annealing_trajectory()` (Arrhenius first-order ODE) and `multi_step_annealing()`. Notebook: defect recovery at elevated temperatures.
**Avoids:** Pitfall 11 (Si activation energies — use SiC-specific E_a; Z1/2 does not anneal at room temperature)

### Phase Ordering Rationale

- Phases 1-2 are strictly required before any fluence simulation: defect physics must exist before it can be wired to the devsim solver.
- Phase 3 before 4-5 because CCE vs fluence is the primary result; dark current and C-V are secondary validation targets.
- Phase 6 is deliberately last: multi-defect convergence issues (Pitfall 6) are much easier to debug after the single-effective-defect pipeline is validated across the full fluence range. The Burin model's anomalously large EH6/7 sigma_e (9e-12 cm^2) will stress the Newton solver.
- The regression constraint drives Phase 1 to include the test infrastructure before any simulation code — not as a post-hoc addition.
- Phase 7 is fully decoupled from Phases 1-6 and can be dropped entirely without affecting the core deliverables.

### Research Flags

Phases likely needing deeper research during planning:

- **Phase 5 (NIEL Scaling):** NIEL values for 62 MeV protons in SiC must be obtained from the SR-NIEL web calculator (www.sr-niel.org) before implementation begins. This is a manual 10-minute task, not a software integration, but must not be forgotten. The scaling factor from 252 MeV (MedAustron data) to 62 MeV (Petringa) significantly changes c_r and Phi_crit.
- **Phase 6 (Multi-Defect Model):** EH6/7 charge character (donor vs acceptor) is actively debated in literature. The Burin model assigns it as a donor with anomalously large sigma_e (9e-12 cm^2 — 450x larger than Z1/2). If the donor assignment is wrong, the sign of its doping contribution reverses. Validate against the C-V shift direction in the first available Petringa irradiation data.
- **Phase 3-5 (Carrier Removal Rate at 62 MeV):** No direct measurement of c_r at 62 MeV exists in the literature. Published values are for 252 MeV (4.2-6.4 cm^-1) and 55 MeV (~75 cm^-1). The 62 MeV value must be NIEL-interpolated and carries significant uncertainty. Phi_crit for the Petringa device is highly sensitive to this value.

Phases with standard patterns (research not needed):

- **Phase 1 (Damage Physics Foundation):** Defect introduction (N = g*Phi), degraded lifetime (1/tau = 1/tau_0 + K*Phi), and effective doping (N_eff = N_D - c\*Phi) are textbook TCAD physics. Implement directly from Burin et al. parameter tables.
- **Phase 2 (Device Integration):** Follows the exact pattern of the existing temperature integration. The create_sic_device() + devsim.set_parameter() pipeline is already established.
- **Phase 7 (Annealing):** Arrhenius first-order ODE with scipy.integrate.solve_ivp is a standard numerical problem. Uncertainty is in the SiC-specific activation energies, not the solver approach.

## Confidence Assessment

| Area         | Confidence | Notes                                                                                                                                                                                                                                                                                                                                                              |
| ------------ | ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Stack        | HIGH       | No new packages needed; all v2.0 physics achievable with existing devsim + scipy + numpy. Zero dependency risk. devsim's CreateNodeModel supports arbitrary trap levels without solver changes.                                                                                                                                                                    |
| Features     | HIGH       | Table stakes features are well-defined by domain conventions; deferral decisions are clear and well-justified. Feature dependency chain is fully mapped.                                                                                                                                                                                                           |
| Architecture | HIGH       | Fluence-as-temperature-analogue pattern maps cleanly onto the existing codebase. Anti-patterns are explicitly identified. Build order is dependency-driven with no ambiguity. MEMORY note about graded doping profile is addressed by Approach A (node model).                                                                                                     |
| Pitfalls     | MEDIUM     | Integration pitfalls (regression, double-counting) are HIGH confidence from codebase analysis. Physics pitfalls (g vs c confusion, linear vs log lifetime) are MEDIUM from literature. Numerical edge cases near Phi_crit are inferred from commercial TCAD reports (GTS TCAD convergence issues explicitly reported by Burin et al.), not yet verified in devsim. |

**Overall confidence:** MEDIUM-HIGH

### Gaps to Address

- **NIEL values for 62 MeV protons in SiC:** Must run SR-NIEL calculator before Phase 5 begins. Populate the `NIEL_PROTON_SIC` dict in `radiation_damage.py` (currently stubbed with `None` placeholders). Estimated effort: 10 minutes. This is a blocker for any proton-specific predictions.
- **Carrier removal rate at 62 MeV:** No direct measurement exists. Must NIEL-scale from 252 MeV MedAustron data (c_r = 4.2–6.4 cm^-1). At 62 MeV, NIEL is substantially higher, so c_r will be higher, making Phi_crit for the Petringa device lower than 1.6e13 p/cm^2. Compute the range and document it as an uncertainty band on all Phi_crit-dependent predictions.
- **Linear vs logarithmic lifetime model for 62 MeV protons:** Empirically unresolved. Implement both behind a flag; present both curves in notebooks. State this uncertainty explicitly in any publication draft.
- **EH6/7 donor vs acceptor assignment:** Literature is split; the Burin model uses donor with anomalously large sigma_e. This must be validated against Petringa group experimental data as soon as any irradiated C-V or I-V measurements are available.
- **Graded doping profile under damage:** Confirmed from MEMORY: uniform N_D fails at reverse bias, graded profile is already implemented. Radiation damage carrier removal must be applied via the node model (Approach A: N_eff(x) = N_D(x) - c_r\*Phi), not as a scalar parameter override. This is explicitly addressed in the architecture.

## Sources

### Primary (HIGH confidence)

- Burin et al., arXiv:2407.16710 (2024) — three-defect model (Z1/2 + EH6/7 + EH4): introduction rates, capture cross-sections, TCAD methodology, validation against neutron-irradiated pad diodes
- Burin et al., arXiv:2407.11776v3 (2024) — extended five-defect model, GTS Sentaurus TCAD, forward/reverse I-V matching across full fluence range
- arXiv:2510.11304 (2025, MedAustron) — carrier removal rate 4.2–6.4 cm^-1 for 252.7 MeV clinical proton beams; in-situ C-V measurement
- Chen et al., CPB 28(1):010701 (2019) — Z1/2 capture cross-sections validated by DLTS, CCE model
- SR-NIEL web calculator (www.sr-niel.org) — NIEL(E) tabulations for protons in SiC

### Secondary (MEDIUM confidence)

- arXiv:2503.09016 (Luo et al., 2025) — logarithmic lifetime degradation model, EH3 introduction rate 1.48 cm^-1, counterintuitive dark current decrease; 80 MeV protons, different energy from Petringa
- IEEE Access 10538275 (2024) — temperature and fluence dependence of carrier lifetime in 4H-SiC (K_tau Arrhenius behavior)
- Hazdra, pss-a 2100218 (2021) — non-linear lifetime degradation in 4H-SiC bipolar devices, Z1/2 as dominant recombination center
- ScienceDirect S136980012300464X (2023) — predictive analytical model for carrier removal using NIEL; power diodes, not detectors
- Tudisco et al., Front. Phys. 10:898833 (2022) — SiC detector radiation hardness review
- Hiyoshi & Kimoto (2009); Hornos et al. (2011) — Z1/2 thermal stability (>1150 C annealing onset)
- Hazdra et al. (2019), IET — comprehensive displacement damage + TID review for 4H-SiC power devices

### Tertiary (LOW-MEDIUM confidence)

- Storasta & Bergman (2004) — EH1/EH3 low-temperature annealing at 300–500 C; limited quantitative data for proton-specific conditions
- Annealing kinetics parameters for room-temperature operation — synthesized from multiple sources; no single validated study for 62 MeV proton-irradiated thin-epi detectors

---

_Research completed: 2026-03-24_
_Ready for roadmap: yes_
