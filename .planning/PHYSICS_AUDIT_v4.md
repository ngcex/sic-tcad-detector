# Deep Physics Audit — SiC TCAD Simulator (Petringa / INFN-LNS)

**Method:** 7 domain auditors → 83 adversarial verifications (2 independent skeptical reviewers per finding) → advisor synthesis.
**Result:** 43 raw findings → **41 survived** double verification, 1 refuted, 1 unverified (rate-limited).
**Date:** 2026-06-13. 94 agents, ~3.6M tokens.

---

## 1. Executive verdict

The simulator is a **well-engineered numerical framework with sound electrostatics**, but the **headline scientific claim — "first TCAD explanation of plasma recombination in SiC under FLASH" — is not supported by the current code.** Three independent confirmations establish this: (a) there is **no plasma/ambipolar/funneling physics** anywhere; the only high-injection term is Auger, which is **~14 orders of magnitude too small** to affect CCE at FLASH dose rates; (b) radiative recombination is described but **never added to the equations**; (c) CCE artifacts (CCE>1) are **clipped rather than corrected** in both transient and single-particle paths, masking the absence of real physics. In its current state the FLASH dose-rate dependence of CCE is an **artifact, not a prediction.**

Separately, several **calibration knobs are presented as physics**: the dominant dark-current term `G_eff = N_t` is a single-point fudge factor (and its Hurkx Γ enhancement is identically 1 and applied with the wrong sign), the NIEL hardness factors are placeholders likely wrong in magnitude _and inverted_, and the tissue-equivalence κ is computed from non-physical stopping-power CSVs (matches the known `project_kappa_flat` issue). There are also **two concrete computational bugs**: the two-carrier Hecht equation double-counts each e-h pair (CCE→~2.0 before clipping), and `compute_ni()` triple-counts conduction-band valleys.

**Overall confidence for paper-readiness: LOW.** The v1.0 electrostatics/C-V foundation is solid and publishable; the FLASH, radiation-damage, and microdosimetry chapters need real physics before any of their results can be called predictions.

---

## 2. Critical issues (must fix before any publication)

### C1 — Two-carrier Hecht double-counts each e-h pair

`src/charge_collection.py:83-131`. CCE sums electron + hole contributions as if they were independent charges, giving CCE ≈ 2.0 clipped to 1.0. **Fix:** the two-carrier Hecht for a pair created at depth x is `CCE(x) = (λ_e/d)[1−exp(−(d−x)/λ_e)] + (λ_h/d)[1−exp(−x/λ_h)]` where the two terms share the _single_ pair — already normalized to ≤1, no clip needed. **Affects:** every analytical-Hecht benchmark, the "CCE validated against Hecht" claim (notebooks 03, 16), all DD-vs-Hecht comparisons.

### C2 — No FLASH plasma physics; Auger 14 orders too small

`src/flash_recombination.py` (whole module), `:245 cce_vs_dose_rate`. The "plasma recombination" is just an Auger term that is negligible at the relevant carrier densities, and there is no ambipolar/screening/funnel model. **Fix:** either (a) implement the actual high-injection physics (ambipolar diffusion, conductivity modulation, field screening by the e-h plasma) and show it dominates at FLASH rates, or (b) **retract the "plasma recombination explanation" claim** and reframe as a sensitivity study. **Affects:** the entire core-value proposition, notebooks 04 and 08.

### C3 — Radiative recombination described but not implemented

`src/sic_material.py:77 (B=1.5e-12)`, `flash_recombination.py:74,91,98`. The B coefficient exists and is documented as part of the high-injection model but is never wired into the recombination equation. **Fix:** add `R_rad = B(np − n_i²)` to the continuity recombination, or remove the claim. **Affects:** any statement that FLASH modeling includes radiative recombination.

### C4 — CCE>1 artifacts clipped, not corrected

`single_particle.py:382-387` (clip [0,1]), `transient.py:356-357` (clip [0,2]). Clipping hides a charge-conservation/normalization error. **Fix:** find the root cause (almost certainly the same pair-double-counting family as C1, or a generation-integral normalization), correct it, and remove the clip so CCE is intrinsically bounded. **Affects:** all single-particle and transient CCE numbers.

### C5 — EH6/7 electron capture cross-section ~1000× too large

`src/radiation_damage.py:81 (sigma_n_EH67 = 9e-12 cm²)`. Physical capture cross-sections are ~1e-14–1e-16 cm²; 9e-12 cm² is unphysical and makes EH6/7 dominate K_τ. **Fix:** correct to the Burin-cited value (verify against the paper; likely 9e-15 or 9e-16). **Affects:** carrier-lifetime degradation, CCE-vs-fluence (notebook 10), Φ_crit.

### C6 — NIEL hardness factors placeholder, wrong magnitude & inverted

`src/radiation_damage.py:49-54`, applied at `:638-639, :429, :482`. Comments admit placeholders; verification finds they are likely both wrong in magnitude and inverted relative to real proton / 1-MeV-neutron-equivalent NIEL. **Fix:** load real NIEL from SR-NIEL for the proton energies used; the module already accepts a table, so this is data not code. **Affects:** every cross-energy fluence comparison, scale_to_proton_energy, Φ_crit, radiation-hardness optimization (notebook 13).

### C7 — Temperature I-V sweep uses midgap-SRH-only device → leakage 13 orders too low

`src/temperature_sweep.py:79`, `drift_diffusion.py:114-117`. The T-sweep builds a device with midgap SRH (n1=p1=n_i) and **no TAT/N_t term**, giving ~1e-25 A vs the real ~18 pA. The temperature coefficient extracted from this is meaningless. **Fix:** use the calibrated dark-current model (the `N_t`/TAT path) in the T-sweep, with a physically T-dependent N_t. **Affects:** the clinical temperature-coefficient result (notebook 06).

### C8 — Dominant dark-current term N_t has zero temperature dependence

`src/dark_current.py:143`, `sic_material.py:68`. `G0_TAT = N_t` is a constant, so the model has the wrong activation energy for Z1/2-limited SiC leakage. **Fix:** give the generation rate the correct Arrhenius T-dependence (`∝ exp(−E_a/kT)` with E_a tied to the Z1/2 level / `n1`), then recalibrate the single prefactor. **Affects:** any temperature-dependent dark-current claim; couples with C7.

---

## 3. Major issues

- **M1 `Nt-fudge`** — `G_eff=N_t` is a single-point empirical fit mislabeled as trap density / TAT physics (`dark_current.py:14-18,140-162`). Defensible as a calibration, **not** as predictive TAT. Relabel honestly.
- **M2 `gamma-inert`** — Hurkx Γ ≡ 1 at all detector fields; the module presents inert TAT as the mechanism (`dark_current.py:218-248`).
- **M3 `gamma-wrong-placement`** — Γ multiplies n1/p1 _inside the SRH denominator_ (wrong sign of effect) instead of enhancing the emission rate (`dark_current.py:157-162`).
- **M4 `nc-valley-doublecount`** — `compute_ni()` triple-counts conduction valleys: NC=5.1e19 vs correct 1.69e19, contradicting the stored `NC_300` (`sic_material.py:131-141`).
- **M5 `varshni-eg0-inconsistent`** — two E_g(0) values: `compute_ni` uses 3.265, `bandgap()` uses 3.2965625 (`sic_material.py:136 vs :26,205`).
- **M6 `pair-energy-too-high`** — E_pair = 8.4 eV vs measured ~7.7–7.8 eV for 4H-SiC; scales all deposited-charge → CCE/microdosimetry by ~8% (`sic_material.py:78`, `single_particle.py:43`, `generation_profiles.py`).
- **M7 `kappa-flat-nonphysical-data`** — stopping-power CSVs non-physical → flat κ≈0.58 mislabeled as energy-dependent (known `project_kappa_flat`). Needs real PSTAR+SRIM data.
- **M8 `kappa-wrong-energy-variable`** — κ looked up with _deposited_ energy, not particle _kinetic_ energy (`microdosimetry.py:347-355`).
- **M9 `eh67-` family / `carrier-removal-equals-z12`** — removal rate set = Z1/2 introduction rate 1:1 with no justification (`radiation_damage.py:93 vs :72`).
- **M10 `incomplete-ionization-fudge`** — 10–30% ionization target at N_A=1e19 is ~10–30× above the defensible value (`incomplete_ionization.py:66-109`).
- **M11 `3d-electrode-junction-mismatch`** — anode placed on n-type; junction geometry inconsistent with contacts (`alternative_structures.py:741-781`).
- **M12 `mesa-subtrench-ptype-epi`** — mesa sub-trench doped p+ through the full epi depth (`alternative_structures.py:519-540`).
- **M13 `sweep-undepleted-configs`** — parametric sweep includes non-fully-depleted configs, corrupting CCE-uniformity ranking (`optimization.py:47,107-130`).
- **M14 `no-plasma-ambipolar-physics`** — (see C2) flagged major in transient domain.
- **M15 `radiative-claimed-not-implemented`** — (see C3).
- **M16 `cce-clip-masks-artifacts`** — (see C4) 1 confirm / 1 refute → real but lower confidence.
- **M17 `auger-test-cannot-detect-effect`** — Auger validation tests are trivially satisfied (`tests/test_flash_recombination.py:109-170`).

## 4. Minor & observations (18)

n_i internal inconsistency (5e-9 vs 8.5e-9 first-principles); v_th uses DOS mass (~33% low); alpha range 15 µm vs ~18–20 µm; noise floor shot-only mislabeled "intrinsic limit"; SRV is post-hoc add-on not a BC; effective-LET = edep/thickness conflation; f(y) under-normalization biases y_F low; Gibbs ionization missing N_V factor; score_structures unit-mixing rank-instability; binned moments from geometric bin centers; E_g(300)=3.26 forced vs Varshni 3.23; "harder proton" terminology inverted in a test; SRH lifetime T-sign worth checking; proton range scaling neglects stopping-power; δE/E thickness bookkeeping can go negative; compute_ni Eg(0) mismatch distorts n_i(T) slope.

## 5. Calibration vs prediction — honesty ledger

| Quantity                                            | Status                                      | Can the paper call it a prediction? |
| --------------------------------------------------- | ------------------------------------------- | ----------------------------------- |
| Electrostatics, V_bi, depletion, **C-V (R²=0.998)** | First-principles + validated                | **Yes**                             |
| Dark current ~18 pA                                 | Single-point fit (`N_t`), Γ inert, no T-dep | **No** — calibration only           |
| Temperature coefficient of dark current             | Built on wrong (midgap) device              | **No** — invalid                    |
| CCE vs bias (steady-state)                          | Sound DD, but Hecht benchmark bugged        | Qualified yes after C1              |
| **FLASH CCE vs dose-rate**                          | No plasma physics; Auger negligible         | **No** — artifact                   |
| CCE vs fluence / Φ_crit                             | EH6/7 σ wrong, NIEL placeholder             | **No** until C5+C6                  |
| Tissue-equivalence κ                                | Non-physical CSV, flat                      | **No** until real data              |
| Microdosimetric y_F/y_D                             | ICRU-correct math, but E_pair & κ off       | Qualified after M6/M8               |

## 6. Severity-ranked action list

1. **C1** Fix two-carrier Hecht double-count — **S** (single function + test)
2. **C4** Root-cause CCE>1 and remove clips — **M**
3. **M4/M5** Fix `compute_ni` valley count + unify E_g(0) — **S**
4. **M6** E_pair 8.4 → 7.8 eV — **S**
5. **C5** EH6/7 σ_n correction — **S** (one constant, verify vs Burin)
6. **C6** Load real NIEL table — **S–M** (data)
7. **C2/C3/M14/M15/M17** FLASH: implement real plasma+radiative physics _or_ retract claim & reframe — **L**
8. **C7/C8** Wire calibrated dark-current model into T-sweep + add Arrhenius N_t(T) — **M**
9. **M1/M2/M3** Relabel N_t honestly; fix Γ placement & range — **M**
10. **M7/M8** Real κ data + fix energy variable — **M** (data + S code)
11. **M11/M12/M13** Alternative-structure doping/contact geometry + depletion-filtered sweep — **M**
12. Minor batch (n_i, v_th, alpha range, noise label, etc.) — **S each**

**Refuted (do not act):** `hecht-electron-lifetime` (both reviewers: the lifetime usage is correct for this geometry).

---

## 7. Fix log (2026-06-13)

Fixes applied this session, each double-audited by two independent agents
(one physics-correctness, one over-engineering/hallucination):

| ID | Fix | Status |
|---|---|---|
| C1 | Hecht two-carrier double-count → depth-averaged single-pair form (`charge_collection.py`). Re-derived by symbolic integration; intrinsically bounded ≤1. | **Fixed + verified** |
| M4 | `compute_ni` removed spurious ×M_c valley factor (NC 5.08e19 → 1.69e19, matches stored/Ioffe). | **Fixed + verified** |
| M5 | `compute_ni` E_g(0) 3.265 → 3.2965625 eV, consistent with `bandgap()` (E_g(300)=3.260). | **Fixed + verified** |
| M6 | E_pair 8.4 → 7.8 eV (measured 4H-SiC W-value ≈7.83 eV) across material/generation/single-particle. | **Fixed + verified** |
| C2/C3 | FLASH module + PROJECT.md: honest caveat that plasma/ambipolar/radiative physics is NOT implemented and the "plasma recombination explanation" claim is unsupported. | **Documented (no physics invented)** |
| C5 | EH6/7 σ_n = 9e-12 cm² flagged as suspect transcription error (correct value needs Burin Table I). | **Flagged, not changed** |
| C8 | `dark_current` N_t T-independence flagged (wrong activation energy; needs Arrhenius + recalibration). | **Flagged, not changed** |

**Deliberately NOT fixed** (need external data or substantial re-derivation — would risk fabricating physics):
- **C5** EH6/7 σ value — needs Burin et al. 2024 Table I.
- **C6** NIEL hardness factors — needs SR-NIEL calculator data.
- **C7/C8** dark-current T-dependence — needs Arrhenius model + single-point recalibration.
- **C2/C3** real FLASH plasma/ambipolar physics — L-effort new modeling.

**C4 correction (from fix-audit):** the over-engineering reviewer hypothesized C4's
CCE>1 clips in `single_particle.py:387` and `transient.py:357` shared C1's root cause
and might become removable. On inspection they do **not** — those clips mask a *separate*
BDF1 generation-pulse displacement-current overshoot (CCE from Q_collected/Q_generated,
not from the Hecht formula). C4 remains a genuine open critical needing its own fix
(root-cause the transient charge-conservation artifact rather than clipping).

**C7 manifested as a flaky test (discovered during fix-verification):**
`test_temperature_sweep::test_iv_current_increases_with_temperature` asserted that
reverse leakage rises with T. The midgap-SRH-only sweep device produces I_reverse
~1e-12..1e-14 A (solver residual, not physics), which is non-monotonic in T both
BEFORE (I=[9.8e-13, 2.7e-12, 1.9e-12]) and AFTER (I=[2.3e-12, 2.7e-12, 3.6e-14])
the n_i fix — it passed pre-fix only by luck. n_i(T) itself is correctly monotonic
(8.8e-10 -> 8.5e-9 -> 7.1e-8). Marked `xfail` (strict=False) referencing C7; it will
pass once the calibrated TAT/N_t(T) model is wired into the sweep. This is evidence
FOR the audit's C7 finding, not a regression from the fixes.

**Test updates** (consequences of the correct E_pair 8.4->7.8 fix, not bugs):
`test_generation_profiles::test_dose_rate_conversion_value`,
`test_single_particle::test_generation_integral_matches_expected` (+2 tautological
LET tests) updated to the corrected constant.
