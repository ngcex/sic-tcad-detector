# Feature Landscape — v4.0 Extended Physics

**Domain:** 4H-SiC TCAD microdosimeter simulator — new physics capabilities for v4.0
**Researched:** 2026-05-17
**Scope:** NEW features only. v1.0/v1.1/v2.0/v3.0 capabilities (1D electrostatics, FLASH, radiation damage, 2D mesh, single-particle transients, MC coupling, y-spectra, alternative structures, parametric optimization with shot-noise floor) are NOT re-researched here.

Each section is labeled **table stakes**, **differentiator**, or **anti-feature** for the v4.0 milestone. Complexity is rated relative to the existing devsim codebase (~20k LOC, 27 modules), not absolute.

---

## 1. Complete Noise Analysis (Beyond Shot Noise)

**Classification:** Table stakes for any paper claiming microdosimetric "energy resolution" or "minimum detectable lineal energy". v3.0 has only shot-noise floor estimation — reviewers will reject this as incomplete.

### Physics

A SiC PIN microdosimeter coupled to a charge-sensitive preamplifier has four physically distinct noise sources that contribute incoherently (sum of variances) to the Equivalent Noise Charge (ENC):

1. **Shot noise (parallel current noise)** — already in v3.0. Spectral density `S_I,shot = 2·q·I_leak`. Dominant when I_leak is high (post-irradiation, high temperature).
2. **Thermal / Johnson-Nyquist noise (series voltage noise)** — Johnson noise on bulk resistance and FET input. `S_V,thermal = 4·k_B·T·R_eq`. Dominant when detector capacitance C_d is large or shaping time is short.
3. **1/f noise (flicker)** — mobility-fluctuation origin in SiC Schottky/PIN, characterized by the Hooge parameter α_H. `S_V,1/f = α_H · V² / (N·f)` where N is the total carrier number in the active region. In 4H-SiC Schottky diodes the noise spectral density follows `S_I(I) ~ I^β` with β ≈ 1.6–1.7 at 300 K (Levinshtein et al.). Independent of shaping time in ENC.
4. **Generation-Recombination (G-R) noise from deep traps** — Z1/2 (E_c – 0.65 eV) and EH6/7 (E_c – 1.55 eV) act as G-R centers. Each trap contributes a Lorentzian spectrum `S(f) = 4·N_t·τ / (1 + (2πfτ)²)` with τ set by the trap emission rate. After irradiation, trap concentrations rise (Burin params already in `radiation_damage.py`) and G-R noise can dominate.

### ENC Formulation (CR-RC shaping)

For a Gaussian/CR-RC shaper of peaking time τ_s, the standard ENC² (Spieler, LBL ICFA notes) decomposes as:

```
ENC²_total = a · (S_V) · C_d² / τ_s     (series / voltage)
           + b · (S_I) · τ_s             (parallel / current)
           + c · A_f · C_d²              (1/f, ~independent of τ_s)
```

where a, b, c are shaping-dependent form factors of order unity, A_f is the 1/f coefficient, and C_d is the **detector capacitance** at operating bias. This is a direct dependency on v3.0's C-V module — without C_d(V) from `cv_analysis.py`, the voltage-noise term cannot be computed.

### Figures of Merit

| Quantity                             | Symbol           | Typical SiC microdosimeter target                      | Computed from              |
| ------------------------------------ | ---------------- | ------------------------------------------------------ | -------------------------- |
| Equivalent Noise Charge              | ENC (e⁻ rms)     | 60–200 e⁻ at room temp (cf. silicon SOI µdos: ~150 e⁻) | Sum of variance terms      |
| Noise Equivalent Energy              | NEE = ENC · ε_eh | < 1 keV (ε_eh = 8.4 eV in SiC → 60 e⁻ ≈ 0.5 keV)       | ENC × pair-creation energy |
| Energy resolution                    | ΔE/E or FWHM     | ~12% at 660 keV imparted (silicon ref)                 | 2.355 · NEE / E_signal     |
| Minimum detectable y                 | y_min (keV/µm)   | 0.1–0.4 keV/µm (silicon ref); lower for SiC            | NEE / mean chord length    |
| Lower-level discrimination threshold | LLD              | Set above 3·σ_noise; affects y_F, y_D                  | NEE × safety factor        |

The choice of LLD is non-trivial: Bianchi 2021 (`The effect of different lower detection thresholds in microdosimetric spectra`) shows y_F can shift by tens of percent depending on cutoff. This must be a knob in v4.0, not a hard-coded value.

### Complexity: MEDIUM

- Formulas are textbook (Spieler/Radeka). Implementation is algebraic, not PDE.
- The hard part is parameter sourcing: α_H for 4H-SiC PIN (not Schottky) is poorly tabulated; literature values for α span 10⁻⁹ to 10⁻³. Recommend treating α_H as a calibration knob with a literature-bounded prior.
- G-R noise: trap densities already come from `radiation_damage.py`; need to add emission-time-constant τ_t(T) per defect (Arrhenius from Burin activation energies — already in v2.0).

### Dependencies on v3.0

- **`cv_analysis.py`** → C_d(V) at operating bias (REQUIRED — voltage-noise term cannot be computed without it)
- **`dark_current.py`** → I_leak(V, T, Φ) for shot-noise term (already used by v3.0 shot-noise floor)
- **`radiation_damage.py`** → trap densities N_t(Φ) and emission constants for G-R noise
- **`optimization.py`** → existing scoring function gets new term: replace `shot_noise_only` with `total_ENC`

### Anti-features

- Don't try to model amplifier 1/f noise from first principles — make it an input parameter (datasheet-driven).
- Don't model microphonic / pickup noise. These are EMC, not device physics.

---

## 2. Build-up / Over-response in 2D (Near-Surface CCE Deficit)

**Classification:** Table stakes. Without this, the simulated y-spectra will systematically over-predict low-y events compared to experiment, and reviewers will ask "where is your dead-layer correction?"

### Terminology Clarification (CRITICAL)

"Build-up over-response" has two distinct meanings in the literature, and the project's Italian description ("campo elettrico near-surface, zona di non raccolta, correzione") points to definition (b):

- **(a) Photon dosimetry sense:** Electronic disequilibrium near the surface in MV photon beams; build-up cap thickness matters. **NOT what is meant here.**
- **(b) Charged-particle / microdosimetry sense:** The first ~µm below the p+ contact has reduced CCE because (i) the layer is undepleted at low bias, (ii) heavy doping increases SRH recombination, (iii) surface recombination velocity S_eff drains carriers to the contact, (iv) electric field is weak so drift is slow. Protons crossing this region deposit energy that is partially collected → distorts the low-y tail of the spectrum.

This document assumes definition (b) throughout. State this explicitly in the user-facing notebook.

### Physics

A 1D profile of CCE(z) starting from the p+/n-epi metallurgical junction typically shows:

- **z ∈ [0, w_dead]:** CCE < 1, monotonically rising. w_dead is set by depletion depth into p+ side (~10–100 nm in heavily doped p+), plus a recombination-limited zone whose extent depends on S_eff and τ_p in the p+.
- **z ∈ [w_dead, w_depl]:** CCE ≈ 1 (full drift in field).
- **z > w_depl (undepleted epi at low V):** CCE < 1, diffusive collection with characteristic length L_n = √(D_n · τ_n).

The 2D extension adds lateral structure:

- Near the SV perimeter (the corner where the top metal ends), the electric field has a 2D fringe pattern. Below the edge of the contact, the depletion extends laterally and CCE is reduced not only by depth but by lateral position too.
- Result: a 2D CCE map CCE(x, z) that is NOT separable as CCE_x(x) · CCE_z(z).

For 4H-SiC specifically:

- Pair-creation energy ε_eh = 8.4 eV → low-y tail is at ~0.1–1 keV/µm where ENC matters too. Build-up correction and noise threshold are coupled, not independent.
- SRV at the SiC surface is poorly characterized; v1.1 has an empirical S_eff already calibrated to the 18 pA dark current target.

### "Build-up correction": Experimental vs Simulated

- **Experimental:** Take a CCE(depth) profile from ion-microbeam scanning (e.g., proton microprobe IBIC) → divide raw pulse-height distribution by CCE(z_track). Requires assumption about where ion entered.
- **Simulated:** Use the 2D CCE map as a _weighting function_ during MC coupling. Each Geant4 step (x_i, z_i, ΔE_i) contributes ΔQ_i = ΔE_i · CCE_2D(x_i, z_i) / ε_eh to the event charge. Then y = Q_event · ε_eh / (q · ℓ_chord).

The v3.0 `mc_coupling.py` already does the event-by-event integration but with CCE=1 (or a coarse depth-only profile). The v4.0 upgrade is the 2D weighting function.

### Figures of Merit

| Quantity                                     | Typical SiC value                           | Notes                                                   |
| -------------------------------------------- | ------------------------------------------- | ------------------------------------------------------- |
| Dead-layer thickness w_dead                  | 0.1–1 µm                                    | Function of p+ doping, bias                             |
| Surface recombination velocity S_eff         | 10³–10⁶ cm/s                                | Calibration parameter                                   |
| Edge CCE deficit (lateral)                   | 10–30% reduction within 1–5 µm of perimeter | From silicon SOI µdos literature; SiC should be similar |
| Over-response at low y                       | "Bump" in y-spectrum below 1 keV/µm         | Compare with/without correction                         |
| Active fraction (active vol / geometric vol) | > 0.95 desired                              | Determines effective chord length                       |

### Complexity: MEDIUM-HIGH

- The physics is in `device2d.py` already — need mesh refinement near top contact (current mesh is uniform-ish in z).
- 2D CCE map extraction requires running the 2D adjoint or beam-induced charge calc at many (x_inj, z_inj) points → potentially slow. Consider sparse sampling + interpolation.
- The coupling back to `mc_coupling.py` is a one-line API change but requires re-running all v3.0 y-spectra.
- Anisotropic mobility (feature 4) couples here: vertical drift in dead layer uses µ*⊥ (holes) and µ*∥ (electrons) — get the tensor right or the dead-layer thickness will be off.

### Dependencies on v3.0

- **`device2d.py`** → mesh refinement near z=0; existing SRV boundary condition stays
- **`charge_collection_2d.py`** → extend to return full CCE_2D(x, z) map, not just integrated CCE
- **`mc_coupling.py`** → ingest 2D CCE map as weighting function
- **`single_particle.py`** → already does drift-diffusion of point charges; reuse for CCE map sampling

### Anti-features

- Don't try to model the SiO2/SiC interface trap density from first principles. Treat S_eff as one calibration knob.
- Don't model the metal contact as anything but ohmic. Schottky-effect on the p+ contact is second-order for a properly fabricated device.

---

## 3. Azimuthal / Angular Response (CCE vs Angle of Incidence)

**Classification:** Differentiator for proton/ion microdosimetry papers. Most TCAD design papers stop at normal incidence; including angular response makes the work clinically relevant.

### Terminology Clarification (verify with PI before coding)

The Italian phrase `risposta azimutale` and `sweep angolare in 2D con approssimazione 3D` is ambiguous:

- **Polar angle θ** (incidence angle relative to surface normal): doable in 2D r–z mesh by rotating the track vector while keeping the mesh fixed. This is what most "angular response" papers actually mean.
- **True azimuthal angle φ** (rotation around the surface normal): for a circular cross-section this is trivial (symmetry), but for the square 100×100 µm and 300×300 µm Petringa SVs, the response varies with φ because corner-to-edge distance differs. True φ-dependence on a square SV requires 3D.

**Recommended interpretation for v4.0:** Polar-angle sweep θ ∈ [0°, 80°] on the existing 2D mesh, with explicit caveat that true azimuthal (φ) variation over the square SV is deferred to v4.0's full-3D feature (item 8 in the project plan). Document this limit prominently.

### Physics

For a charged particle entering at angle θ:

1. **Track length in SV** scales as ℓ = t_SV / cos(θ) for a track that fully crosses the 10 µm epi (until θ is so large the track exits a side face).
2. **Lineal energy** y = ε_imparted / ℓ_chord, where ℓ_chord is averaged over the particle direction distribution. For a parallel beam at θ, ℓ_chord = t/cos(θ) for transit tracks.
3. **CCE per track** depends on which part of the 2D CCE map the track samples. Glancing-angle tracks spend more length in the dead layer → lower CCE per unit length.
4. **Edge events** at large θ: track enters top face but exits a side face → partial energy deposition in SV, additional low-y events.

For microdosimetry, the _spectrum_ y(θ) is the observable. Specifically:

- **y_F(θ)** (frequency-mean lineal energy): expected to scale as 1/cos(θ) at small θ.
- **y_D(θ)** (dose-mean): less sensitive to θ because high-y events dominate the integral.
- **RBE inferred via MKM**: depends on y_D — angular stability of y_D is a clinical figure of merit.

### Clinically Relevant Angles

Proton therapy at INFN-LNS uses pencil-beam scanning with deliberate angulation; clinical relevance:

- **0° (normal incidence):** baseline reference
- **15°–30°:** typical "gantry off-axis" beam in passive scattering
- **45°–60°:** lateral beam in patient frame, common in head-and-neck
- **80°–90°:** edge of clinical relevance, but useful for stress-testing the simulation

Sweep recommendation: θ ∈ {0, 15, 30, 45, 60, 75}° at minimum.

### Figures of Merit

| Quantity            | What it tells you                                 |
| ------------------- | ------------------------------------------------- |
| ⟨ℓ_chord(θ)⟩        | Mean chord length vs angle; needed to normalize y |
| y_F(θ) / y_F(0)     | Frequency-mean angular response                   |
| y_D(θ) / y_D(0)     | Dose-mean angular response (more stable)          |
| Q_collected(θ)      | Charge per event vs angle (couples to noise)      |
| Edge-event fraction | Fraction of events that don't fully cross the SV  |

### Complexity: MEDIUM

- 2D mesh stays the same. The change is in `single_particle.py`: track injection direction becomes a (θ, φ_fixed) parameter.
- LET integration along the angled path: existing `mc_coupling.py` walks the Geant4 step list, just rotate the step coordinates before depositing.
- For 2D r–z, only θ is meaningful; the mesh has implicit cylindrical symmetry. Document this clearly.

### Dependencies on v3.0

- **`single_particle.py`** → add direction vector (currently assumes normal incidence)
- **`mc_coupling.py`** → rotate Geant4 steps before lookup into 2D CCE map
- **`microdosimetry.py`** → existing y-spectrum routines work unchanged; just call them at each θ
- Implicit dependency: build-up 2D (feature 2) — angular response is most interesting where the CCE map is non-trivial

### Anti-features

- Don't simulate the full pencil-beam angular distribution (Gaussian-divergent beam, multiple scattering). Use a δ-function angle for v4.0; beam divergence is v5+.
- Don't claim "azimuthal" sweep on the square SV without 3D — the 2D r–z mesh literally cannot resolve corner-vs-edge differences.

---

## 4. Anisotropic Mobility in 4H-SiC (Tensor Model)

**Classification:** Differentiator. Most published SiC TCAD treats mobility as scalar — including the tensor is a clear novelty and matters quantitatively for thin SVs where drift time competes with diffusion time.

### Physics (Get the Direction Right)

4H-SiC is hexagonal (6mm point group). The mobility tensor is uniaxial: distinct values along c-axis (⟨0001⟩, the depth direction in standard (0001)-cut wafers like Petringa's) and perpendicular to it (basal plane).

**Electrons (from Ishikawa 2023, Schadt & Pensl 1994):**

- μ_∥c (along c-axis, vertical drift in PIN) ≈ **1140–1160 cm²/V·s** at N_D = 2 × 10¹⁵ cm⁻³, T = 300 K
- μ_⊥c (perpendicular, lateral) ≈ **947 cm²/V·s** at same conditions
- Ratio μ*⊥ / μ*∥ ≈ **0.83** for electrons
- **Note:** 4H-SiC is OPPOSITE to 6H and 15R polytypes (where μ*⊥ > μ*∥). Vertical electron transport in 4H is _favored_ — this is why (0001) vertical devices have an advantage.

**Holes (from Ishikawa 2024 JAP, Kagamihara 2018):**

- μ*⊥c > μ*∥c by 20–50% (opposite to electrons!)
- Channel mobility 28 cm²/V·s perpendicular (highest reported for SiC p-MOSFETs)
- Bulk hole mobility values much lower than electron (~120 cm²/V·s perpendicular at low N_A)

**Implications for the Petringa device (p+/n−epi/n+ substrate on (0001)):**

- After ion generates an e-h pair, electrons drift up (toward n+ substrate, along +c-axis) with HIGH mobility μ_∥c.
- Holes drift down (toward p+ contact, along −c-axis) with LOW mobility μ_∥c (hole).
- Hole drift time τ*drift,h ~ t_epi / (μ*∥c,h · E) is the rate-limiting step.
- Lateral diffusion uses different mobilities for e and h — affects edge response and 2D CCE.

### Is the Effect Experimentally Distinguishable?

Yes, but subtle:

- For normal incidence (axial drift), only μ_∥ matters → no anisotropy visible.
- For tilted tracks (feature 3) and edge events (feature 2), the lateral component matters and a tensor model differs from isotropic by ~15–20% in transit time.
- For transient pulse shape (induced current i(t)), tensor model shifts the peak timing by a few hundred ps — possibly detectable with a fast preamp.
- Most published SiC TCAD assumes isotropic ≈ μ*∥. Going to a tensor with μ*⊥ = 0.83·μ_∥ is a small numeric change but a meaningful physics correction.

### Figures of Merit

| Quantity                       | Isotropic vs Tensor difference                 |
| ------------------------------ | ---------------------------------------------- |
| Hole transit time              | ~5–10% (axis-aligned tracks); ~15–20% (tilted) |
| Lateral diffusion length L_n,⊥ | √0.83 ≈ 9% reduction in horizontal direction   |
| 2D CCE map at edges            | Few % CCE change near perimeter                |
| Pulse rise time                | tens to hundreds of ps shift                   |

### Complexity: LOW-MEDIUM

- devsim supports anisotropic mobility via tensor parameters in the mobility model (`Mobility` node-edge model can take per-direction values, or you provide `mu_xx, mu_yy, mu_zz` separately).
- The hard part is consistent application across ALL existing transport simulations — every v3.0 notebook that uses 2D transport will rerun with slightly different numbers.
- Doping dependence: Caughey-Thomas-style fits exist for both μ*∥ and μ*⊥ (Ishikawa 2023 provides empirical equations). Both must be implemented or the model is incomplete.
- Temperature dependence: T^(-α) with α slightly different for ∥ and ⊥ (Ishikawa: α*∥ ≈ 2.4, α*⊥ ≈ 2.0).

### Dependencies on v3.0

- **`sic_material.py`** → replace scalar μ*n(N_D, T), μ_p(N_A, T) with tensor (μ*∥, μ_⊥)(N, T)
- **`device.py` and `device2d.py`** → wire tensor into devsim mobility model
- **All transport simulations** → re-run; expect small but non-zero diffs in regression tests
- **`temperature_sweep.py`** → update T-scaling

### Anti-features

- Don't add anisotropic effective mass for incomplete-ionization calc — second-order on top of second-order.
- Don't model field-dependent anisotropy in the saturated-velocity regime (negligible at the modest fields here, ~10⁴ V/cm).

---

## 5. Graded Epi Doping in 2D

**Classification:** Table stakes. The known v3.0 tech debt (`uniform N_D fails at reverse bias`) directly affects the 2D electrostatics; without graded doping in 2D, the bias-dependent CCE maps will be wrong above the breakdown of the uniform-doping approximation (around −15 V in current model).

### Physics — What Profile Is Physically Correct?

4H-SiC homoepitaxy by CVD (typical for the Petringa device) produces a doping profile shaped by:

1. **Buffer layer:** ~0.5–1 µm thick, doped intermediate between substrate (~10¹⁸–10¹⁹ cm⁻³ p+) and active epi (~10¹⁴ cm⁻³). Buffer doping typically 10¹⁶–10¹⁷ cm⁻³.
2. **Active epi region:** Nominally uniform, but in practice has:
   - **"W-shape" radial profile** across the wafer (gas-flow and temperature-field artifact). Within a single ~mm-scale device this is negligible — assume laterally uniform.
   - **Slight ramp in z** due to memory effects (residual dopant from previous growth) → typically 10–30% increase from top to bottom of active layer.
3. **Substrate interface:** Sharp transition (30–50 nm) from buffer to substrate. The buffer is the engineered "graded" zone.

**Reasonable parameterization for v4.0:**

```
N_D(z) = N_D,sub                          for z > t_epi + t_buf  (substrate)
N_D(z) = N_D,buf · exp(-(z-t_epi)/L_grade) for t_epi < z < t_epi + t_buf  (buffer)
N_D(z) = N_D,epi · (1 + δ · (z/t_epi))     for 0 < z < t_epi   (active)
```

with δ ~ 0.1–0.3 (slight increase toward substrate) and L_grade ~ 100 nm.

The existing 1D v1.0 graded-doping fit (R²=0.998 to C-V) provides anchor points for the 1D z-profile.

### How Does 2D Differ from 1D Graded?

In 1D: N_D = N_D(z) only. Lateral uniformity is assumed and is correct for the active region of a planar device far from edges.

In 2D: For the Petringa planar device, N_D(x, z) ≈ N_D(z) almost everywhere — lateral non-uniformity is real but small (few %). The 2D-specific issue is at edges:

- **Guard ring region** (if simulated, v3.0 alternative-structures): doping under the guard implant differs from active region.
- **Mesa-etched SV** (v3.0 alt structure): sidewalls expose epi at different z-depths → laterally varying _effective_ doping at the surface.

**Practical recommendation:** Keep N_D(z) only (1D profile in 2D mesh), unless modeling a mesa or implanted guard. The big win over v3.0 is fixing the uniform-N_D-fails-at-reverse-bias issue, not adding lateral profile.

### Figures of Merit

| Quantity                     | v3.0 (uniform)               | v4.0 (graded) target                   |
| ---------------------------- | ---------------------------- | -------------------------------------- |
| Depletion width vs V         | Saturates incorrectly        | Matches Hecht/analytical at all V      |
| C-V curve in 2D              | Diverges from 1D above −15 V | Matches 1D within 5% across full range |
| CCE vs V curve               | Has a kink near −15 V        | Smooth, monotonic to 100%              |
| Doping-bias optimization map | Has unphysical region        | Physically valid everywhere            |

### Complexity: LOW-MEDIUM

- 1D graded profile is already in `device.py` (v1.0). Pattern: assign N_D as a function of node z-coordinate.
- 2D port: same pattern with `device2d.py` mesh, using node coordinates. Devsim API supports this directly via `node_solution` or position-dependent material parameters.
- Re-run all v3.0 notebooks that used 2D transport — regression check needed. Expect changes at reverse bias > −15 V.
- The kappa-flat issue (in MEMORY) is unrelated — graded doping doesn't fix the stopping-power table issue.

### Dependencies on v3.0

- **`device2d.py`** → mesh contains nodes with (x, z) coords; assign N_D(z) by z-coordinate
- **`device.py`** → 1D version already done, reuse the same N_D(z) function
- **All v3.0 2D simulations** → regression test for changes
- **`optimization.py`** → existing parametric studies (geometry × doping × bias) interpret "doping" as the active-epi value; semantics unchanged

### Anti-features

- Don't try to fit a profile from C-V data alone — degenerate. Use the published Petringa profile + small calibration.
- Don't model dopant compensation by carbon vacancies (Z1/2 acceptor state) here — that's already in v2.0 radiation damage path.
- Don't model 2D lateral non-uniformity beyond the buffer/substrate transition unless modeling a non-planar structure (mesa, guard).

---

## Feature Dependency Graph

```
[5 Graded epi 2D]  ──► fixes V_bias > 15V regime; precondition for everything below
        │
        ▼
[2 Build-up 2D]    ──► provides CCE_2D(x, z) map
        │
        ├──► [3 Azimuthal/angular response]  (samples the CCE map along tilted tracks)
        │
        └──► [1 Complete noise analysis]     (uses C_d, I_leak, traps to get ENC)
                       │
                       ▼
            [Updated y-spectra + minimum detectable y]
                       │
[4 Anisotropic mobility]  ──── (couples to ALL transport; do early or all results shift)
```

**Recommended phase ordering for the roadmap:**

1. Graded epi 2D (foundation, fixes known tech debt)
2. Anisotropic mobility (affects all downstream — do before generating new "correct" baselines)
3. Build-up 2D CCE map (needs the corrected 2D electrostatics)
4. Complete noise analysis (needs C_d from corrected 2D; can run in parallel with item 5)
5. Azimuthal response (consumes the corrected CCE map + correct mobility tensor)

Items 4 and 5 can be developed in parallel after items 1–3 ship.

---

## MVP Recommendation for v4.0

If milestone scope must contract, prioritize in this order:

1. **Graded epi 2D** — closes a known correctness gap; non-negotiable for paper-quality bias sweeps.
2. **Complete noise (Johnson + 1/f + G-R)** — biggest scientific delta; the "minimum detectable y" claim is currently shot-noise-only and reviewers will catch it.
3. **Build-up 2D** — required for honest comparison with experimental y-spectra; relatively contained scope.
4. **Anisotropic mobility** — physics novelty, small code change, but touches everything (regression burden).
5. **Azimuthal (polar-angle θ in 2D)** — clinically motivating but scientifically the most defer-able if a stretch goal.

Defer to v5: full 3D azimuthal φ-sweep, beam divergence, anisotropic mobility in full 3D geometry, real Geant4 ROOT integration with sample file from LNS.

---

## Anti-Features (v4.0-wide)

| Anti-feature                                                               | Why avoid                                           | What to do instead                                         |
| -------------------------------------------------------------------------- | --------------------------------------------------- | ---------------------------------------------------------- |
| Full 3D mesh for azimuthal φ-sweep                                         | devsim 3D performance unproven; massive memory cost | Polar-angle θ in 2D + caveat about square-SV symmetry      |
| Real Geant4 ROOT file integration                                          | Waiting on sample file from INFN-LNS                | Synthetic ROOT-format CSV fixture as in v3.0               |
| Modeling preamp circuit in SPICE                                           | Out of scope; not device physics                    | ENC parameters from datasheet as inputs                    |
| Self-consistent kappa from stopping-power tables WITHIN the noise pipeline | Tangential to noise analysis                        | Keep kappa fix in its own feature (separate v4.0 item)     |
| Surface-state DLTS calibration from scratch                                | Beyond the 1–2 phase scope per feature              | Treat S_eff as a single calibration knob (already in v1.1) |
| Microphonic / EMC noise                                                    | Not device physics                                  | Out of scope                                               |
| Beam divergence + multiple scattering inside the SV                        | Geant4's job, not ours                              | Use Geant4 step list as-is from input                      |

---

## Confidence Levels

| Feature                           | Physics confidence | Numerical-value confidence                                                          | Sources                                                          |
| --------------------------------- | ------------------ | ----------------------------------------------------------------------------------- | ---------------------------------------------------------------- |
| Noise: shot + Johnson + 1/f + G-R | HIGH (textbook)    | MEDIUM (α_H for 4H-SiC PIN poorly tabulated)                                        | Spieler LBL notes; Levinshtein on 4H-SiC Schottky                |
| Build-up 2D                       | HIGH               | MEDIUM (SRV for SiC poorly characterized, already a v1.1 calibration knob)          | Frontiers 2025 SiC dosimetry review; CMRP silicon SOI literature |
| Azimuthal angular response        | HIGH               | HIGH (only geometry)                                                                | Rosenfeld et al. SOI µdos; CMRP 3D cylindrical                   |
| Anisotropic mobility              | HIGH               | HIGH                                                                                | Ishikawa 2023 phys stat sol; Schadt & Pensl 1994                 |
| Graded epi 2D                     | HIGH               | MEDIUM (Petringa-specific profile from published C-V fit; buffer thickness assumed) | 4H-SiC CVD epi growth literature; existing v1.0 1D fit           |

---

## Sources

- [Experimental and Theoretical Study on Anisotropic Electron Mobility in 4H-SiC (Ishikawa 2023, phys. stat. sol. b)](https://onlinelibrary.wiley.com/doi/10.1002/pssb.202300275)
- [Anisotropy of the electron Hall mobility in 4H, 6H, and 15R silicon carbide (Schadt & Pensl)](https://www.semanticscholar.org/paper/Anisotropy-of-the-electron-Hall-mobility-in-4H,-6H,-Schadt-Pensl/a52d8bffd3888f3a3deb0cc219c9dc514c37d1da)
- [Origin of hole mobility anisotropy in 4H-SiC (Ishikawa 2024, J. Appl. Phys.)](https://pubs.aip.org/aip/jap/article/135/7/075704/3265791/Origin-of-hole-mobility-anisotropy-in-4H-SiC)
- [Electronic Noise — Helmuth Spieler (LBL ICFA notes)](https://www-physics.lbl.gov/~spieler/ICFA_Istanbul/pdf/III_Electronic_Noise.pdf)
- [Noise in Semiconductor Devices (Auburn)](https://www.eng.auburn.edu/~wilambm/pap/2011/K10147_C011.pdf)
- [Low frequency and 1/f noise in wide-gap semiconductors: SiC and GaN](https://www.researchgate.net/publication/3349335_Low_frequency_and_1f_noise_in_wide-gap_semiconductors_Silicon_carbide_and_gallium_nitride)
- [1/f noise in forward biased high voltage 4H-SiC Schottky diodes (Levinshtein et al.)](https://www.sciencedirect.com/science/article/abs/pii/S0038110114000458)
- [Silicon carbide sensors in radiotherapy dosimetry — Frontiers in Sensors 2025](https://www.frontiersin.org/journals/sensors/articles/10.3389/fsens.2025.1622153/full)
- [Microdosimetry for hadron therapy: a state of the art (Frontiers in Physics 2022)](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2022.1035956/full)
- [Cylindrical SOI microdosimeter: charge collection characteristics (Rosenfeld et al., IEEE)](https://ieeexplore.ieee.org/document/4436607/)
- [Effect of different lower detection thresholds in microdosimetric spectra (Bianchi 2021)](https://www.sciencedirect.com/science/article/abs/pii/S1350448721001384)
- [Correction factors to convert microdosimetry measurements in silicon to tissue (Bolst et al.)](https://pubmed.ncbi.nlm.nih.gov/28151733/)
- [Proton microbeam studies of CCE in large area SiC detectors (NIM A 2025)](https://www.sciencedirect.com/science/article/pii/S0168900225010563)
- [Position-resolved charge collection of SiC detectors with epitaxial graphene layer (Sci Rep 2024)](https://www.nature.com/articles/s41598-024-60535-3)
- [Refining SiC epi-growth for high-volume production (Compound Semiconductor)](https://compoundsemiconductor.net/article/106637/Refining_SiC_epi-growth_for_high-volume_production)
- [Estimation of Electron Drift Mobility along the c-Axis in 4H-SiC (vertical Schottky barrier diodes)](https://www.researchgate.net/publication/383382162_Estimation_of_Electron_Drift_Mobility_along_the_c-Axis_in_4H-SiC_by_Using_Vertical_Schottky_Barrier_Diodes)
- [Microdosimetry with a 3D SOI detector in a low energy proton beamline](https://www.sciencedirect.com/science/article/pii/S0969806X20302103)
