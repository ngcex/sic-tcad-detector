# Physics Audit v5 — 4H-SiC TCAD Radiation-Detector Simulator (Petringa / INFN-LNS)

**Method:** Full independent re-audit across 7 physics domains (material/electrostatics, drift-diffusion + charge collection, dark current + temperature, transient + FLASH, radiation damage + annealing, microdosimetry + MC coupling, alternative structures + optimization + 2D). Every raw finding was double-verified by two adversarial reviewers (literature-check + code-reality-check). This document is the PI-level synthesis: I resolve the two verifiers' disagreements with my own physics judgment, reclassify each surviving finding against the v4 fix history (REGRESSION / STILL-OPEN / NEWLY-FOUND / ALREADY-FIXED), and correct severities where a verifier over- or under-stated impact.

**Date:** 2026-06-14

**One-line vs v4:** v5 confirms all six v4 fixes (C1, C4, C7/C8, M4, M5, M6) held with **zero regressions**, re-confirms the v4 still-open criticals (FLASH plasma physics, EH6/7 sigma_n, NIEL kappa), and surfaces several **newly-found** defects — most importantly a paper-blocking inverted/fabricated tissue-equivalence kappa in the microdosimetry path that v4 never examined.

---

## 1. Executive Verdict

The simulator's **core device physics is sound and the v4 fixes all held** — electrostatics, calibrated dark current (18 pA at -30 V, 300 K), and CCE-vs-bias are defensible as a calibrated TCAD model, and the CCE>1 / Hecht double-count / n_i / E_pair regressions from v4 are genuinely closed. However, the project is **not paper-ready for its headline claim**: the "first TCAD explanation of plasma recombination in 4H-SiC under FLASH" remains **unsupported by the code** (no plasma screening, ambipolar transport, conductivity modulation, or track funnelling; the only high-injection channel is an Auger term ~9-14 orders below SRH at FLASH densities), exactly as v4 documented. Two **newly-found critical** defects in the microdosimetry domain (which v4 did not audit) corrupt the tissue-equivalent output: the tissue-equivalence factor kappa is **both inverted and fabricated** (~0.58 from placeholder CSVs vs the physical ~1.13–1.24), and kappa is indexed by deposited energy instead of particle kinetic energy. Two **major** structural defects in the alternative-structures module (3D-electrode graded-profile override; 3D-electrode anode on n-) make those device variants unphysical, and the two known still-open radiation-damage parameters (EH6/7 sigma_n, NIEL kappa) keep all absolute fluence/Phi_crit numbers unreliable.

**Overall confidence for paper-readiness: LOW.** The device-physics chapters (electrostatics, dark current, CCE-vs-bias) can be published as a _calibrated_ model with honesty caveats. The microdosimetry/tissue-equivalence and FLASH chapters cannot be published as predictions until the kappa data and the missing plasma physics are addressed. With the critical fixes below (mostly data + relabeling, modest code), confidence rises to MODERATE for a scoped, honestly-framed paper.

---

## 2. Verification of the v4 Fixes

| Fix                                                                                                                                                                  | What v4 claimed                                                      | v5 independent check                                                                                                                                                                                                                                                                                                                                                                                            | Verdict                                                                                                       |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| **C1** — Hecht two-carrier double-count → depth-averaged single-pair form, intrinsically bounded ≤1                                                                  | `charge_collection.py`                                               | Re-derived Hecht; confirmed single-pair form. No double-count surfaced in any surviving finding. Related Hecht findings (`hecht-tau-role-mismatch`, `Ldiff-inconsistent`) concern lifetime _labeling_/consistency in the analytical benchmark, not the double-count.                                                                                                                                            | **HELD — correctly fixed**                                                                                    |
| **C4** — CCE>1 root-caused: single-particle ~0.6% trapezoidal-quadrature; transient normalization bug fixed; blind clips → physical-ceiling [0,1] + warning if >1.05 | `transient.py`, `single_particle.py`                                 | Re-derived: trapezoidal envelope normalization (`generated_charge_trapezoidal_pulse`) matches the collection window; unit tests sound. `cce-still-clipped` confirms the residual overshoot is quadrature (convex peak → trapezoid over-estimates), **refuting** the auditor's displacement-current hypothesis on sign grounds (displacement integrates to ~0; truncation only under-counts).                    | **HELD — fix correct; v4 root-cause attribution (quadrature) is right**                                       |
| **C7/C8** — temperature-dependent dark current via `nt_temperature_scale(T)=n_i(T)/n_i(300)`, E_a=E_g/2, factor=1.0 at 300 K preserves 18 pA                         | `dark_current.py`                                                    | Confirmed the scaling is implemented and physically defensible (Sah-Noyce-Shockley J*gen ∝ n_i ⇒ E_a≈E_g/2). The T-coefficient is now a \_model*, not invalid. Residual issue: the validating test re-derives the input function (`tshape-imposed-not-predicted`) and the dominant-term activation energy is an _assumption_ about a near-midgap center, not Z1/2 physics (`z12-srh-negligible-vs-ni-scaling`). | **HELD with residual** — fix works; honesty labeling of the T-curve as calibrated (not predicted) is required |
| **M4** — `compute_ni` removed spurious valley factor (NC 5.08e19 → 1.69e19)                                                                                          | `sic_material.py`                                                    | Independently reproduced NC=1.695e19 from m_e=0.77 m0. Correct.                                                                                                                                                                                                                                                                                                                                                 | **HELD — correctly fixed**                                                                                    |
| **M5** — `compute_ni` E_g(0) unified to 3.2965625 eV (consistent with `bandgap()`)                                                                                   | `sic_material.py:144, 26`                                            | Verified both `bandgap()` and `compute_ni` use the identical calibrated triple; E_g(300)=3.26 to machine precision.                                                                                                                                                                                                                                                                                             | **HELD — correctly fixed**                                                                                    |
| **M6** — E_pair 8.4 → 7.8 eV across material/generation/single-particle                                                                                              | `sic_material.py:94`, `generation_profiles.py`, `single_particle.py` | Confirmed 7.8 eV propagates through the core charge-generation path. **BUT** two files were missed: `mc_coupling.py:314` and `optimization.py:41` still hardcode 8.4 eV (`epair-8p4-mc`, `opt-epair-wrong`).                                                                                                                                                                                                    | **HELD in core; INCOMPLETE in two secondary paths**                                                           |

**Bottom line:** no v4 fix regressed. The only fix that was not fully propagated is M6 (two stale 8.4 eV constants in MC-coupling and optimization).

---

## 3. Critical Issues (current) — must fix before publication

### C-1 [NEWLY-FOUND] Tissue-equivalence kappa is inverted AND fabricated (~0.58 vs physical ~1.13–1.24)

- **What is wrong:** `compute_kappa_table` (`microdosimetry.py:302`) correctly computes kappa = S*water/S_SiC, but the input CSVs (`data/stopping_power_water.csv`, `data/stopping_power_sic.csv`) are placeholder/fabricated. CSV water proton stopping at 1 MeV is 51.55 MeV·cm²/g vs NIST PSTAR ~260.8 (off ~5×, non-constant ratio 3–5× across the curve), and SiC is listed \_higher* than water. Result: kappa ≈ 0.575–0.587, flat to <1%. Physics requires kappa > 1 (water Z/A=0.555 vs SiC 0.499 → 1.11×; water lower I ~78 eV vs SiC ~136 eV further raises the ratio), with real PSTAR giving ~1.24 (1 MeV) → ~1.13 (100 MeV), a ~10% monotonic decrease. The code's own docstring warning "kappa varies ~20-30%" is contradicted by its own data (<1%).
- **Severity reconciliation:** both verifiers confirmed **critical**; I concur. This is a sign-inverted, ~2× wrong, mis-shaped correction on the detector's headline clinical output.
- **Downstream affected:** every tissue-equivalent lineal-energy spectrum, y_F, y_D, quality-factor/RBE estimate after `tissue_equivalence_correction`; Notebooks 18/19 microdosimetry outputs; the `kappa_constant=0.58` default path.
- **Fix:** replace both CSVs with real NIST PSTAR liquid-water proton mass stopping powers and SiC (Bragg additivity from PSTAR Si and C, mass fractions Si 0.700 / C 0.300, or ASTAR/SRIM). Set `kappa_constant` default to ~1.2. Correct the docstring to "~1.13–1.32, increasing at lower energy." **Effort: M** (data sourcing + regen).

### C-2 [NEWLY-FOUND] kappa(E) indexed by deposited energy instead of particle kinetic energy

- **What is wrong:** `tissue_equivalence_correction` (`microdosimetry.py:347-355`) interpolates kappa(E) using per-event _deposited_ energy in the thin SV, but stopping-power tables are indexed by _kinetic_ energy. A penetrating proton deposits ~10⁴× less than its kinetic energy in a 10 µm SV, so the lookup samples the wrong (low-energy) end for every penetrating particle.
- **Severity reconciliation:** auditor called it critical; both verifiers downgraded to **minor** _today_ because the shipped kappa table is flat (<2% over the whole range, and deposited energies fall below the 0.1 MeV table floor → clamped). I rate it **critical-conditional**: it is currently masked only by the _other_ bug (the fabricated flat table). Once C-1 is fixed and kappa(E) has its true ~10% slope, this wrong-variable lookup becomes a real error. It must be fixed _together with_ C-1, hence listed as critical.
- **Downstream affected:** energy-dependent tissue-equivalent spectra (same notebooks as C-1).
- **Fix:** pass the per-event _primary kinetic energy_ from MC truth (available upstream of the `groupby('event_id')` in `mc_coupling.py`), not `edep_keV.sum()`. If kinetic energy is unrecoverable for a given event, fall back to the energy-averaged kappa and document it. Fix the misleading docstring at `microdosimetry.py:331-332`. **Effort: M.**

### C-3 [STILL-OPEN, v4 C2/C3] FLASH plasma physics not implemented — headline claim unsupported

- **What is wrong:** the only high-injection loss channel is Auger; no e-h plasma field screening, ambipolar transport, conductivity modulation, or track funnelling. At 20–230 Gy/s the steady excess density (5.9e8–3.5e11 cm⁻³) is 240–8000× below doping — firmly low injection, nowhere near the plasma regime. Auger/SRH ≈ 1e-14–1e-9; any dose-rate CCE trend is numerical, not physical.
- **Severity reconciliation:** literature-check kept **critical** (the unsupported physics IS the project's stated core value); code-reality downgraded to major because the module honestly self-discloses. I keep it **critical** for _publication purposes_: a paper cannot make the headline claim from this code, regardless of how honest the docstring is.
- **Downstream affected:** `cce_vs_dose_rate`, `parametric_cce_sweep`, Notebooks covering FLASH; the project Core Value statement in v1.0/v1.1 requirements.
- **Fix (choose one):** (a) implement field-dependent screening / ambipolar mobility / conductivity modulation and validate against measured roll-off (**Effort: L**), or (b) remove the "first TCAD explanation of plasma recombination" claim from all project-facing materials and relabel the dose-rate outputs as exploratory sensitivity bounds (**Effort: S**).

### C-4 [STILL-OPEN, v4 C5] EH6/7 electron capture cross-section dominates electron K_tau

- **What is wrong:** `sigma_n_EH67 = 9e-12 cm²` (`radiation_damage.py:87`) implies a 17 nm capture radius (~55 lattice constants) and contributes ~92% of electron K_tau, so every post-irradiation tau_n and electron-side CCE-vs-fluence is off by ~an order of magnitude. The test asserts the wrong value, locking it into CI.
- **Severity reconciliation:** here the two verifiers genuinely split, and the split is important. Literature-check **retrieved Burin et al. arXiv:2407.16710 Table I directly** and found it prints EH6,7 sigma*e = 9e-12 — i.e. the code is a \_faithful transcription*, NOT a transcription error, and the auditor's proposed "9e-15 per Burin" is wrong. Code-reality (no paper access) kept critical on pure physical-implausibility grounds. **My resolution:** the value is physically implausible for a point defect (~92% K*tau dominance, 17 nm radius), so the \_symptom* is real and critical for any electron-lifetime result — but the _root cause_ is an uncertain value inherited from the source paper (which itself flags its cross-sections as deviating widely in the literature), not a code slip. Net severity: **critical** for downstream electron-lifetime predictions; **but** the fix is "replace with a primary-measured value," not "fix a typo."
- **Downstream affected:** all tau_n-vs-fluence, electron-limited CCE-vs-fluence, EH6/7 annealing recovery of tau_n. Hole side unaffected.
- **Fix:** replace with a primary DLTS-measured EH6/7 electron cross-section (~1–2e-14 cm², capture radius ~0.6 nm, consistent with the same defect's hole sigma_p=3.8e-14 and neighboring Z1/2/EH4). Keep the AUDIT C5 caveat; re-attribute the suspect value to the source paper rather than calling it a transcription slip; update the test. **Effort: S** (one constant + test), but requires a literature decision.

### C-5 [STILL-OPEN, v4 C6] NIEL hardness factors are unvalidated placeholders, magnitude uncertain

- **What is wrong:** `NIEL_HARDNESS_PROTON_SIC` = {30:0.50, 62:0.35, 70:0.33, 150:0.22} are self-labeled placeholders; kappa multiplies fluence everywhere (`fluence_neq = fluence_proton * kappa`), so all absolute proton-damage and Phi_crit numbers ride on these.
- **Severity reconciliation:** the two verifiers split hard. Code-reality (anchoring to _silicon_ RD48 values) kept it major, claiming ~3–5× too small. Literature-check **refuted** the 3–5× claim: the code normalizes to the _SiC_ 1-MeV-neutron NIEL (correct denominator for a SiC simulator), whereas the auditor compared against _silicon_-normalized factors; the SiC neutron denominator is ~2.5–3× larger, so SiC-normalized proton kappa of ~0.3–0.5 is _plausible_, not 3–5× low. **My resolution:** the specific "3–5× too small / should be ~1–2" claim is **not established** — it conflated Si and SiC normalization. What survives is the _original v4 truth_: these are unvalidated placeholders that must be replaced with SR-NIEL SiC data, and may be off by tens of percent (not 3–5×). I keep this **STILL-OPEN** and rate it **major→critical-for-absolute-numbers**: until replaced, no _absolute_ Phi*crit or defect-concentration number is citable, though the energy \_trend* is correct.
- **Downstream affected:** CCE-vs-proton-fluence, N_Z12/N_EH67/N_EH4, Phi_crit headline, 62-MeV default predictions.
- **Fix:** replace with SR-NIEL SiC proton NIEL at 30/62/70/150 MeV, expressed consistently as kappa(E)=NIEL_proton(E)/NIEL_neutron_SiC(1 MeV). Document the normalization convention explicitly to avoid the Si/SiC confusion. **Effort: M.**

---

## 4. Major Issues

### Mj-1 [NEWLY-FOUND] 3D-electrode device discards the graded epi profile (`alternative_structures.py:770-775`)

For default `doping_profile='graded'`, `_apply_doping` installs the graded Donors field but returns the scalar `DEFAULT_N_D=1.07e15`; line 774 then overwrites Donors with a uniform 1.07e15 step, silently discarding the graded profile (intended bulk 8.5e13, ~12.6× lower). Both verifiers confirmed **major**. Bulk N_D error of 12.6× changes depletion width (~3.5×), field shape, and full-depletion voltage for this variant. **Fix:** preserve the graded field (rename to `Donors_epi` before overwrite) and set `Donors = max(Donors_column, Donors_epi)`, not the scalar. **Effort: S.** Isolated to `create_3d_electrode_device` (telescope path is correct).

### Mj-2 [NEWLY-FOUND] 3D-electrode anode sits on n- epi, not p+ (`alternative_structures.py:717-728`)

No radial p+ is added at the outer wall; acceptors are depth-only `step(junction_pos - y)`, so the "anode" contacts p+ over only the top ~9% (1 µm of 11 µm) and is a metal-on-n- ohmic contact over ~91%. The intended radial p⁺/n⁻ junction does not exist; the real junction is the inherited horizontal plane. Both verifiers confirmed **major**. The codebase already knows how to dope a wall (guard-ring variant at lines 1289-1305). **Fix:** add a full-depth x-dependent p+ acceptor overlay at the outer radius. **Effort: S.** Isolated to the 3D-electrode variant.

### Mj-3 [NEWLY-FOUND] CCE-uniformity sweep ranks non-fully-depleted configs (`optimization.py:44-169`)

`microdosimetric_sweep` ranks by edge/center CCE ratio alone, with no full-depletion check. For graded doping the true V*fd (integrating N_D(y)) reaches ~120–191 V for several grid points, exceeding the swept biases, so partially-depleted devices with both-low CCE can show a deceptively uniform ratio and rank artificially high. Both verifiers confirmed **major**; code-reality correctly notes the \_production* notebook grid (epi=10 µm, V up to 80 V) limits the blast radius to ~2 of 9 points (N_D=5e14 at 20/50 V), but the function defaults reach the deeply non-depleted regime. **Fix:** compute V_fd by integrating the graded N_D(y), flag/exclude V_bias < V_fd, and add an absolute-CCE floor before ranking. **Effort: M.**

### Mj-4 [STILL-OPEN] Transient I(t) omits displacement/Ramo current (`drift_diffusion.py:187-212`)

`extract_contact_current` queries only the Electron/Hole continuity equations; the displacement term lives on the unqueried PotentialEquation contact charge model. Integrated charge Q (hence CCE) is fine (displacement integrates to ~0 over an equilibrium→equilibrium transient), but the reported I_peak and t_collection are conduction-arrival current, not the physical Ramo-induced signal. Both verifiers confirmed **major**. **Fix:** add `I_disp = get_contact_current(..., equation="PotentialEquation")` so terminal current = conduction + displacement; no separate weighting-field solve needed. **Effort: M.** Pulse-shape/timing observables must not be compared to measured waveforms until fixed.

### Mj-5 [STILL-OPEN, v4 M7/M8 adjacent] Generation profiles omit the Bragg peak / use wrong depth law

Three related defects: (a) `alpha_generation_profile` is a flat erfc box with no Bragg rise (`generation_profiles.py:107-130`); (b) the 2D lateral-scan alpha depth profile _decays_ exponentially — backwards for LET and inconsistent with the 1D model (`charge_collection_2d.py:262`); (c) proton profile is always flat, computed SiC range discarded. For the standard 10 µm fully-depleted epi these are **minor** (peak lies beyond the epi; normalization cancels in CCE), but they become **major** for `cce_vs_epi_thickness` sweeps ≥15 µm where the full Bragg peak falls inside the device. I rate the cluster **major** because at least one shipped study (epi sweep to ~20 µm) is biased, and the two alpha models are mutually inconsistent. **Fix:** replace boxes with a Bragg-Kleeman LET profile `~(R-x)^(1/p-1)`, p≈1.75, truncated at the epi edge; unify the 1D and 2D alpha depth models; amend docstrings (they currently frame the box only as anti-ringing smoothing). **Effort: M.**

---

## 5. Minor Issues and Observations

**Dark-current model honesty (all ALREADY-corrected-in-behavior, documentation/labeling only):**

- `gamma-identically-unity` (**major→labeling**): Hurkx Γ≡1 for all fields below SiC breakdown (Kt=180–510 at operating fields; Kt<4 needs ~4.5 MV/cm). TAT field enhancement is inoperative; voltage dependence comes entirely from the `min(E/E_ref,1)` depletion selector and width growth. Γ≈1 is _physically correct_ for a 0.65 eV trap at these fields — the defect is the docstring advertising active TAT. **Fix: relabel.**
- `gamma-wrong-placement-srh` (**major→minor**): Hurkx Γ misplaced (multiplies n1/p1 instead of dividing lifetimes, wrong sign), but the affected SRH term is ~31 orders below the dominant generation term and Γ≡1 anyway. Dormant. **Fix: rewrite for hygiene.**
- `gamma-subunity-discontinuous` (**minor→observation**): Schenk Γ<1 and discontinuous at Kt=4, but unreachable below breakdown. Note: the auditor's proposed "correct" asymptotic also falls below 1 at low field; the right fix is `Γ=max(1, asymptote)`. **Fix: hygiene.**
- `z12-srh-negligible-vs-ni-scaling` (**major→minor**): the E*a=E_g/2 scaling is a defensible near-midgap-center \_assumption*, not Z1/2 physics (a true Z1/2-limited term would carry E_a=2.61 eV; the Z1/2 SRH branch is ~1e-17 of the total). **Fix: docstring honesty.**
- `tshape-imposed-not-predicted` (**major, but methodological**): the T-coefficient is correct physics but the validating test re-derives its own input. **Fix: relabel test as a unit test of the scaling function; add a real device-solved Arrhenius slope test.**
- `fdepl-heuristic-field-weight` (**minor**): `min(E/E_ref,1)` is an engineering depletion selector (the code's own comments say so — the "mislabeled as field-dependent generation" framing is a strawman). Acceptable if disclosed.
- `srv-decoupled-and-inert` (**minor**): the SRV model is never attached to a contact equation and evaluates to ~0 at the ohmic contact (np pinned to n_i²) at all biases; I_SRV does not feed I_total. **Fix: remove the misleading I_SRV column or implement SRV at a free surface.**
- `docstring-Nt-value-inconsistency` (**observation**): docstrings say N_t=1e12 / ~8 pA while the code uses 2.2e13 → ~18 pA. **Fix: update lines 18, 115.**

**Material/electrostatics (all observation-level, no result impact):**

- `srh-n1-p1-midgap`: base SRH n1=p1=n*i mislabeled "midgap" for a 0.65 eV trap; over-predicts bulk SRH generation ~16 orders \_relatively*, but absolute leakage is ~3e-26 A (far below 18 pA) and the leakage-relevant path already uses the correct Z1/2 n1/p1. **Observation** (the auditor's "major" overstates absolute impact). **Fix: correct comment, optionally set base n1/p1 to Z1/2 values.**
- `ni-internal-inconsistency`: stored n_i=5e-9 vs first-principles 8.5e-9 (V_bi shifts <0.03 V); a literature anchor, not a fudge. **Observation.**
- `gibbs-missing-pNv`: simplified Gibbs ionization drops p/N_V; mooted by the empirical high-doping branch at device doping. **Observation.**
- `varshni-params-nonstandard`: **REFUTED on premise** — (6.5e-4, 1300) IS the standard Levinshtein/Ioffe NSM value; the auditor's "Choyke 8.2e-4/1.8e3" reference is itself wrong. No issue.
- `numerical-W-50pct-threshold`: depletion edge at n/N_D=0.5 is a defensible engineering definition (~one Debye length, ~8% of W); under bias the code uses the analytical W anyway. **Observation.** (Side note: `poisson.py:317` docstring wrongly says "1% threshold.")

**Charge-collection benchmark consistency (minor):**

- `Ldiff-inconsistent`: hardcoded L*p=7 µm vs self-consistent sqrt(D_p·tau_p)≈14 µm; affects only the partial-depletion \_analytical* benchmark, not the DD solver. **Fix: default L_diff_p = sqrt(mu_p·(kT/q)·tau_p).**
- `hecht-tau-role-mismatch`: **non-issue** — the code correctly pairs each carrier's drift length with its own capture lifetime; the "in p-type" suffix is just a measurement-context label. The conservative direction (1 ns underestimates lambda_e) makes CCE→1 even safer. Docstring-clarity at most.

**FLASH module (minor, self-disclosed):**

- `auger-negligible-at-flash`: Auger ~9–14 orders below SRH (auditor's "17 orders" / n³ framing is wrong — in the doped low-injection regime Auger scales ~N_D²·δp, ratio ~2e-9). Conclusion (negligible) holds. **Observation.**
- `radiative-B-not-wired`: B=1.5e-12 defined but never used; would exceed Auger (crossover n≈2.1e18) but both << SRH. **Minor.** Self-disclosed.
- `flash-tests-trivially-satisfied`: `assert cce_with_auger <= cce_no_auger + 1e-10` passes for zero effect; no test asserts a measurable dose-rate roll-off; `[0.5,1.5]` CCE window exceeds physical ceiling. **Major as a test-validity/false-confidence issue.** **Fix: add a real dose-rate-dependence test or remove the FLASH-validation framing.**
- `cce-still-clipped`: CCE is clip-bounded not intrinsically bounded; sub-5% bias clipped silently. The displacement-current-as-overshoot-source hypothesis is **refuted** (wrong sign). **Observation:** persist `cce_raw` for auditability.
- `proton-range-unused-flat`, `proton-range-scaling-unused`: dead `_R_sic_mm` (noqa F841), flat profile physically valid for thin detectors; tests validate a dead value. **Observation:** remove dead code or implement Bragg-Kleeman.

**Radiation damage (minor/observation):**

- `phicrit-test-docstring-direction`: test docstring mislabels 30 MeV as "harder" and quotes kappa~1.0 vs table 0.50. **Minor docstring fix.**
- `eta-provenance-neq-vs-proton`: **REFUTED** — Burin et al. is a _neutron_-equivalent study; the per-neq labeling and the kappa NIEL-scaling are correct. No double-conversion.

**Microdosimetry / MC (minor):**

- `kappa-range-test-tautology` (**major as a test issue**): `test_kappa_range` asserts 0.3<kappa<1.0, codifying the inverted value; a correct table (>1) would _fail_ this test. **Fix with C-1:** assert ~1.0<kappa<1.4.
- `kappa-test-not-validated`: energy-dependent kappa test asserts nothing about magnitude or trend; lets flat/inverted data pass. **Minor; fix with C-1.**
- `epair-8p4-mc`, `opt-epair-wrong`: stale 8.4 eV (M6 not propagated); ~8% charge under-estimate in MC full-path, ~7.7% high noise floor in optimization. **Minor.** **Fix: set both to 7.8 eV (import from `SiCParameters`).**

**Alternative structures (observation):**

- `mesa-subtrench-pplus`: **REFUTED as critical** — the sub-trench p+ region is electrically inert (no interface to the pillar, solve runs only on the pillar region), so no parasitic junction or cathode short exists. **Observation:** misleading comment/dead doping; set sub-trench to n-/insulating.
- `deltaee-internal-junction`: the n-/p+ adjacency at the de/E interface is benign for the intended two-independent-detector telescope (no transport couples the regions; `CreateSiliconOxideInterface` enforces _potential_ continuity, not blocking). **Observation:** fix the docstring at `alternative_structures.py:862` which falsely claims "current continuity."

---

## 6. Calibration vs Prediction — Updated Honesty Ledger

| Quantity                                                | Status                                                                 | Honest label for the paper?                                                                                          | Change vs v4                                                         |
| ------------------------------------------------------- | ---------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------- |
| Electrostatics / C–V (V_bi, W, field)                   | **Calibrated** (N_D tuned to W(0V)≈1.7 µm; n_i anchored to literature) | Prediction OK _qualitatively_; absolute W is a calibrated/defined quantity (50% edge convention)                     | unchanged                                                            |
| Dark current I(–30 V, 300 K) ≈ 18 pA                    | **Calibrated** single point (N_t=2.2e13 fit)                           | NOT a prediction — single-point calibration; voltage shape is heuristic (depletion selector + width growth, not TAT) | unchanged; TAT framing must be dropped                               |
| Dark-current T-coefficient (Arrhenius E_a≈E_g/2)        | **Calibrated/assumed** (imposed via n_i(T) scaling)                    | NOT a prediction — assumed near-midgap E_a; report as "assumed midgap-generation E_a," not "Z1/2 activation energy"  | **v4 had this 'invalid' after C7/C8; now a valid but assumed model** |
| CCE vs bias (high-bias, 10 µm epi)                      | **Prediction** (DD solver; CCE→1 robust)                               | YES — defensible                                                                                                     | confirmed (C1/C4 held)                                               |
| CCE pulse shape: I_peak, t_collection                   | **Not physical** (conduction-only, no Ramo/displacement)               | NO — do not compare to measured waveforms (Mj-4)                                                                     | unchanged                                                            |
| FLASH CCE vs dose-rate                                  | **Not modeled** (Auger negligible; no plasma physics)                  | NO — exploratory bound only; headline plasma claim unsupported (C-3)                                                 | unchanged (v4 C2/C3 still open)                                      |
| CCE vs fluence / Phi_crit (absolute)                    | **Unreliable** (EH6/7 sigma_n + NIEL kappa placeholders)               | NO for absolute numbers; energy _trend_ OK                                                                           | unchanged (v4 C5/C6 still open)                                      |
| CCE vs fluence — hole side, trend                       | **Qualitative prediction**                                             | Trend OK; absolute depends on NIEL                                                                                   | unchanged                                                            |
| Tissue-equivalence kappa                                | **Wrong** (inverted + fabricated CSVs; wrong lookup variable)          | NO — paper-blocking until C-1/C-2 fixed                                                                              | **NEW — not in v4 ledger**                                           |
| Microdosimetric y_F / y_D (raw, pre-tissue-equivalence) | **Prediction** (from deposited-energy spectra)                         | YES for SiC-medium spectra; NO once multiplied by the broken kappa                                                   | **NEW**                                                              |

---

## 7. Severity-Ranked Action List

1. **[C-1] Replace fabricated stopping-power CSVs** with real NIST PSTAR water + Bragg-additivity/ASTAR SiC; fix `kappa_constant` default and docstring; update `test_kappa_range` bounds to ~1.0–1.4. **Effort: M.** _(paper-blocking)_
2. **[C-2] Use kinetic energy (not deposited) for kappa(E) lookup**; pass primary energy from MC truth; fix docstring. **Effort: M.** _(do with C-1)_
3. **[C-3] Resolve the FLASH headline claim:** either implement plasma/ambipolar/screening physics and validate (**L**), or remove the "first TCAD explanation of plasma recombination" claim and relabel dose-rate outputs as sensitivity bounds (**S**). _(paper-blocking for the headline)_
4. **[C-4] Replace EH6/7 sigma_n** (9e-12 → ~1e-14 primary-measured); re-attribute as inherited-suspect, update test, keep caveat. **Effort: S.**
5. **[C-5] Replace NIEL kappa table** with SR-NIEL SiC proton values at the four beam energies; document the SiC-neutron normalization explicitly. **Effort: M.**
6. **[Mj-1] Fix 3D-electrode graded-profile override** (`max(Donors_column, Donors_epi)`). **Effort: S.**
7. **[Mj-2] Add radial p+ wall to 3D-electrode anode.** **Effort: S.**
8. **[Mj-4] Add displacement current** to `extract_contact_current`; restores physical I_peak/t_collection. **Effort: M.**
9. **[Mj-3] Add full-depletion gate + absolute-CCE floor** to `microdosimetric_sweep`. **Effort: M.**
10. **[Mj-5] Replace flat/box generation profiles** with Bragg-Kleeman LET; unify 1D/2D alpha models; fix docstrings. **Effort: M.**
11. **[M6 cleanup] Propagate E_pair=7.8 eV** to `mc_coupling.py:314` and `optimization.py:41` (import from `SiCParameters`). **Effort: S.**
12. **[Honesty pass] Relabel dark-current TAT/Z1/2/T-coefficient docstrings** as calibrated/assumed; fix the circular `test_activation_energy_is_half_bandgap`; fix N_t docstrings; remove/relabel inert I_SRV column; fix `cce_lateral_scan` depth law; persist `cce_raw`. **Effort: S–M.**
13. **[Hygiene] Fix dead code & comments:** mesa sub-trench p+, delta-E/E "current continuity" docstring, proton dead range, base SRH "midgap" comment, Γ sub-unity clamp, n_i comment, 50%-threshold docstring. **Effort: S.**

---

---

## 8. Fix Log — code-only critical fixes (2026-06-14)

Applied after this audit, each gated by an advisor (PI-level) confirmation pass
BEFORE implementation. Data-dependent criticals (C-1 kappa CSVs, C-5 NIEL table)
were deliberately NOT changed — fixing them requires real PSTAR/ASTAR/SR-NIEL
data we do not have, and fabricating it is the exact defect this audit found.

| ID                     | Action                                                                                                                                                                                                                                                                                                                         | Status                                                              |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------- |
| **M6 residual**        | E_pair 8.4 → 7.8 eV in `mc_coupling.py` and `optimization.py`, now imported from `SiC4H_Parameters` (single source of truth).                                                                                                                                                                                                  | **Fixed** — 32 tests pass; deliberate ~8% Q/noise shift.            |
| **C-3** (FLASH claim)  | Removed "first TCAD explanation of plasma recombination" headline from `PROJECT.md`, `v1.0/v1.1-REQUIREMENTS.md`; relabeled as exploratory sensitivity study. Amended `cce_vs_dose_rate` docstring + notebook 04 summary overclaim ("shows no plasma effects" → "null result of the implemented model"). v4 caveat block kept. | **Fixed (relabel)** — 9 flash tests pass.                           |
| **C-4** (EH6/7 σ_n)    | 9e-12 → 1.4e-14 cm² (order-of-magnitude, primary-DLTS-consistent). Comment re-attributed: faithful transcription of Burin Table I (NOT a typo), overridden as physically implausible (~17 nm radius). Test changed from `== 9e-12` point-lock to plausibility band `1e-15 < σ < 1e-13`.                                        | **Fixed** — 78 tests pass; electron τ_n shifts ~1 order (intended). |
| **C-2** (kappa lookup) | Advisor REJECTED a numeric fix: kinetic energy is NOT in the MC DataFrame schema, and fixing the lookup while C-1's flat table persists would manufacture false correctness. Applied **docstring-honesty only** at `tissue_equivalence_correction`; real fix gated behind C-1 + MC-ingestion extension.                        | **Documented; deferred with C-1**                                   |
| **C-1** (kappa CSVs)   | DATA-BLOCKED. Added a prominent `.. danger::` block to `compute_kappa_table` stating the CSVs are fabricated/inverted (~0.58 vs physical ~1.13–1.24) and the exact external-data fix.                                                                                                                                          | **Documented; needs PSTAR/ASTAR data**                              |
| **C-5** (NIEL table)   | DATA-BLOCKED. Expanded the placeholder comment: recorded that v5 REFUTED the "3–5× too small" claim (Si-vs-SiC normalization confusion); kept open as unvalidated.                                                                                                                                                             | **Documented; needs SR-NIEL SiC data**                              |

**Not yet addressed** (honesty/hygiene batch): transient displacement/Ramo current (Mj-4), Bragg-peak generation profiles (Mj-5), dark-current TAT/Z1/2 docstring relabeling. These remain open for a follow-up pass.

---

## 9. Fix Log — major structure/optimization fixes (2026-06-15)

On the critical path to the foundry deliverable (two detector designs). These
three corrupt the geometry/optimization the design hand-off depends on. Each
advisor-confirmed BEFORE implementation, TDD, then tested.

| ID                                      | Action                                                                                                                                                                                                                                                                                                                                                                               | Status                                                                                                                                                                                                            |
| --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Mj-3** (CCE-uniformity sweep)         | Added `full_depletion_voltage_graded()` to `analytical.py` (integrates the graded N_D(y) moment — the uniform-N_D_bulk estimate is anti-conservative). `microdosimetric_sweep` now computes V_fd per config, flags `is_fully_depleted`/`passes_cce_floor` (CCE_FLOOR=0.95)/`is_valid`, and ranks valid configs above invalid (retained, not dropped) via pure `_rank_sweep_results`. | **Fixed** — 7 new tests; V_fd numerically verified (scipy.quad). Advisor's "120-191 V standard config" range CORRECTED: standard 10µm config = 10.3 V; the ~188 V case is the 20µm/5e14 corner (gate catches it). |
| **Mj-1** (3D-electrode graded override) | `create_3d_electrode_device` rebuilt `Donors_epi` from an EXPLICIT graded expression (not `equation="Donors"`, which caused a real devsim cyclic-dependency crash — caught by the Poisson-solve test), then `Donors = max(Donors_column, Donors_epi)`. Removed dead `existing_donors_expr` call.                                                                                     | **Fixed** — bulk Donors now ~N_D_bulk not 1.07e15; regression test asserts it.                                                                                                                                    |
| **Mj-2** (3D-electrode anode on n−)     | Added a full-depth radial p+ shell (`Acceptors_wall`, 1µm, `max()` with substrate term) at the outer wall + a mesh line at its inner edge so it's resolved. NetDoping re-derived after.                                                                                                                                                                                              | **Fixed** — regression test asserts NetDoping<0 full-depth at the wall; Poisson solves.                                                                                                                           |

**Deferred (advisor: off critical path for the deliverable):** Schottky/JBS device model (microdosimetry = geometry-only per PI decision), Mj-4 displacement current, Mj-5 Bragg profiles, kappa/Phase 27.

**Still the gating item for both designs:** Phase 26 graded-doping 2D calibration (planned+approved, not executed).

---

_End of Physics Audit v5._
