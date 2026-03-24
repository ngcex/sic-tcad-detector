# Domain Pitfalls

**Domain:** Adding radiation damage modeling (v2.0) to a validated 1D 4H-SiC TCAD simulator
**Researched:** 2026-03-24
**Confidence:** HIGH for integration/regression pitfalls (verified against codebase), MEDIUM for physics model pitfalls (literature-sourced with cross-validation), LOW for some numerical edge cases (inferred from commercial TCAD reports, not verified in devsim)

---

## Critical Pitfalls

Mistakes that cause silent physics errors, loss of the validated baseline, or wasted implementation effort.

### Pitfall 1: Fluence=0 Does Not Reproduce Pristine Results (Regression Failure)

**What goes wrong:** The radiation damage module introduces new terms into the SRH recombination model (modified lifetimes, additional trap levels, carrier removal). At fluence=0, ALL of these must reduce to exactly the v1.1 physics. But there are multiple ways this silently breaks:

1. **Additive trap levels that don't vanish.** If you add Z1/2 and EH6/7 as explicit SRH trap levels to replace the existing midgap-trap SRH model, the recombination rate changes even at zero defect concentration because the n1/p1 parameters shift from midgap (n1=p1=n_i) to Z1/2-level values (n1=n_i*exp((E_trap-E_i)/kT), p1=n_i*exp((E_i-E_trap)/kT)). These are different even when N_t=0 because the SRH rate expression U_SRH = (n*p - n_i^2) / [tau_p*(n+n1) + tau_n\*(p+p1)] depends on n1 and p1 in the denominator.

2. **Lifetime model replacement.** The v1.1 code sets `taun` and `taup` as scalar region parameters via `srh_lifetime(T, carrier, params)`. If you replace these with a fluence-dependent formula tau(Phi) = 1/(1/tau_0 + K_tau\*Phi), at Phi=0 this returns tau_0. But tau_0 must EXACTLY equal the v1.1 values (tau_n=1e-9 s, tau_p=6e-7 s at 300K, with T-scaling). If tau_0 comes from a different source or is parameterized differently, you get a silent baseline shift.

3. **Carrier removal modifying effective doping.** If N_eff(Phi) = N_D - c\*Phi is implemented by modifying the Donors node model, at Phi=0 it should return N_D. But if the implementation creates a new node model that replaces the existing graded doping profile, the graded profile calibration (N_D_junction=2.9e15, N_D_bulk=8.5e13, L_transition=1e-4) may be lost.

4. **Dark current model interaction.** The v1.1 Hurkx TAT model in `dark_current.py` already uses Z1/2 trap parameters (E_t=0.65 eV, N_t=2.2e13). Radiation damage introduces ADDITIONAL Z1/2 centers. If the damage module overwrites the existing N_t rather than adding to it, the calibrated 18 pA dark current breaks at fluence=0.

**Why it happens:** The damage module touches the same physics (SRH recombination, doping, traps) that the v1.1 baseline depends on. Any non-additive modification breaks the baseline.

**Consequences:** Loss of ALL v1.1 validated results: C-V R^2=0.998, CCE=100% at V>-40V, dark current calibration. You cannot publish v2.0 radiation damage results if the Phi=0 limit no longer matches experiment.

**Prevention:**

1. Design the damage API as purely additive: `apply_radiation_damage(device, region, fluence)` modifies parameters on top of the existing v1.1 state. When fluence=0, this function does nothing (returns immediately).
2. Lifetime degradation: do NOT replace the lifetime model. Instead, after the v1.1 device is fully created and validated, UPDATE the `taun`/`taup` parameters: `new_taun = 1/(1/old_taun + K_tau_n * fluence)`. This preserves v1.1 values at fluence=0 by construction.
3. Carrier removal: create a NEW node model `Donors_eff` = `Donors - c_removal * fluence` and update the Poisson equation to use `Donors_eff`. At fluence=0, `Donors_eff` = `Donors` identically.
4. **Mandatory regression test:** After implementing damage, run the FULL v1.1 test suite at fluence=0 and verify BIT-IDENTICAL results for C-V, I-V, CCE. Any deviation means the damage module is not properly decoupled.

**Detection:** Run C-V at fluence=0 vs the v1.1 golden values. If R^2 drops below 0.998 or dark current changes by more than 1%, the damage module is leaking into the baseline.

**Phase:** Must be enforced from the very first implementation phase. Build the regression test BEFORE writing any damage code.

---

### Pitfall 2: Using Wrong Damage Constants (Orders of Magnitude Scatter in Literature)

**What goes wrong:** The literature on 4H-SiC radiation damage constants shows enormous scatter depending on irradiation conditions, and picking the wrong values can make your predictions wrong by 10x or more. Specific examples:

**Z1/2 introduction rates:**

- Burin et al. (2024, neutron): g_Z1/2 = 5.0 cm^-1
- Luo et al. (2025, 80 MeV proton): g_EH3 = 1.48 cm^-1 (note: they found Z1/2 did NOT increase significantly under their conditions, contradicting other studies)
- Total point defect introduction: ~19 cm^-1 (only 4% of SRIM-predicted primary vacancies -- massive discrepancy)

**Carrier removal rates:**

- Moscatelli et al. (2006, proton): ~240-260 cm^-1 for protons in p-type
- Recent SiC PIN diodes (2025, 253 MeV proton): 4.2-6.4 cm^-1 (linear donor removal)
- The factor ~50x difference between these reflects different device types, proton energies, and whether carrier removal is measured from C-V (effective doping) or Hall effect (free carriers)

**Lifetime damage constants:**

- Luo et al. (2025): 1/tau = a\*ln(Phi) + b (LOGARITHMIC, not linear) with a=2.4e4, b=1.9e6
- Standard Si model: 1/tau = 1/tau_0 + K_tau\*Phi (LINEAR)
- The functional form itself is disputed for SiC

**Why it happens:** Unlike silicon (50+ years of NIEL-scaling validation), 4H-SiC damage studies span only ~15 years, use different crystal qualities, different irradiation facilities, and different measurement techniques. The NIEL-scaling hypothesis (all displacement damage is equivalent when scaled by NIEL) is NOT validated for SiC to the same degree as for Si. Different defect types dominate at different fluence ranges and temperatures.

**Consequences:** If you use g_Z1/2 = 5.0 cm^-1 (from neutron data) for proton irradiation of the Petringa detector, you may overestimate Z1/2 production. If you use the linear 1/tau model when the actual dependence is logarithmic, you overpredict degradation at high fluence and underpredict at low fluence.

**Prevention:**

1. Make ALL damage constants explicit, named parameters with clear provenance:
   ```python
   @dataclass
   class DamageParameters:
       g_Z12: float = 5.0       # cm^-1, Z1/2 introduction rate (Burin 2024, neutron)
       g_EH67: float = 1.6      # cm^-1, EH6/7 introduction rate (Burin 2024)
       c_removal: float = 5.3   # cm^-1, carrier removal rate (recent proton data)
       K_tau_n: float = ...      # cm^2/s, electron lifetime damage constant
       K_tau_p: float = ...      # cm^2/s, hole lifetime damage constant
       source: str = "Burin2024_neutron"  # provenance tag
   ```
2. Implement BOTH linear and logarithmic lifetime degradation models behind a flag. Compare predictions.
3. Document that damage constants are effective parameters tied to specific irradiation conditions. The Petringa detector uses 62 MeV protons; constants from neutron or high-energy proton studies are approximations.
4. Show sensitivity analysis: vary each damage constant by 2x and plot the effect on CCE. This tells you which constants matter most and which are safely uncertain.

**Detection:** If your CCE vs fluence predictions disagree with published data by more than a factor of 2-3, suspect wrong damage constants before suspecting wrong physics. Check: are you using constants from the right particle type, energy range, and device technology?

**Phase:** Must be addressed in the first implementation phase. Create the DamageParameters dataclass with provenance tracking from day one.

---

### Pitfall 3: Treating SiC Defect Chemistry Like Silicon (Wrong Defects, Wrong Physics)

**What goes wrong:** The silicon radiation damage community uses a well-established "Hamburg model" with two effective trap levels (one donor-like, one acceptor-like) plus stable damage. Attempting to transplant this approach to 4H-SiC fails because:

1. **Different dominant defects.** In Si: V-O (A-center), V-V (divacancy), V-P (E-center). In 4H-SiC: Z1/2 (carbon vacancy at Ec-0.67 eV), EH6/7 (Ec-1.6 eV, possibly also carbon vacancy in different charge state), EH4 (cluster defect at Ec-1.03 eV). The Si defects have NO analogues in SiC.

2. **Different compensation mechanism.** In Si, acceptor-like defects compensate n-type doping straightforwardly. In 4H-SiC, Z1/2 is a negative-U center (captures two electrons, with the second capture requiring MORE energy than the first). This means the standard Shockley-Read-Hall occupancy statistics (single-level, detailed balance) are WRONG for Z1/2. The defect has two charge states (0/2-), not the standard (0/-) of a simple acceptor.

3. **EH6/7 donor vs acceptor controversy.** Burin et al. (2024) model EH6/7 as a donor at Ec-1.6 eV. Other groups treat it as acceptor-like. The Burin model uses an electron capture cross-section of 9e-12 cm^2 for EH6/7, which is extraordinarily large (compare Z1/2 at 2e-14 cm^2). This 450x difference means EH6/7 can dominate recombination despite lower introduction rate. But if EH6/7 is actually donor-type, it ADDS to effective doping rather than compensating it. Getting the type wrong reverses the sign of the doping change.

4. **Annealing is completely different.** Si defects anneal at 60-80C (room temperature annealing is significant). SiC defects require much higher temperatures: Z1/2 is stable to ~1500C in some studies, while other radiation-induced defects (S1, S2) anneal at 300-425C. Room-temperature annealing in SiC is negligible for Z1/2 (unlike Si where it is a major effect).

**Why it happens:** Most TCAD radiation damage tutorials, textbooks, and example scripts assume silicon. The Si community has decades of standardized models. SiC models are still actively debated in the literature (2024-2025 papers still disagree on defect assignments and introduction rates).

**Consequences:** Using Si-derived defect models produces qualitatively wrong predictions: wrong voltage at which CCE degrades, wrong fluence threshold for device failure, wrong annealing behavior.

**Prevention:**

1. Do NOT use the Hamburg model or any Si-specific defect parameterization.
2. Use SiC-specific defect parameters from recent TCAD studies (Burin 2024, Luo 2025). Start with the Burin three-defect model (Z1/2 + EH6/7 + EH4) as the most complete published TCAD parameterization.
3. For the negative-U behavior of Z1/2: as a first approximation, treat Z1/2 as a standard single-level acceptor at Ec-0.67 eV. The negative-U correction changes occupancy at intermediate Fermi levels but is secondary for carrier removal in fully depleted detectors where the Fermi level is pinned near midgap.
4. Set annealing activation energies from SiC literature (not Si). For the Petringa detector at room temperature, assume NO annealing of Z1/2 unless explicitly studying thermal recovery.

**Detection:** If your simulation predicts significant CCE recovery at room temperature over hours-to-days timescales, you are using Si annealing parameters. SiC Z1/2 does not anneal at room temperature.

**Phase:** Must be settled before ANY defect model implementation. Define the defect zoo and their parameters as the FIRST task.

---

### Pitfall 4: Carrier Removal Changing the Device Operating Point (Depletion Width Shift)

**What goes wrong:** Carrier removal reduces effective doping: N_eff(Phi) = N_D - c\*Phi. For the Petringa detector with N_D_bulk = 8.5e13 cm^-3 and c = 5.3 cm^-1:

- At Phi = 1e13 p/cm^2: N_eff = 8.5e13 - 5.3e13 = 3.2e13 (62% reduction)
- At Phi = 1.6e13 p/cm^2: N_eff = 8.5e13 - 8.5e13 = 0 (FULL DEPLETION AT ANY BIAS)
- At Phi > 1.6e13 p/cm^2: N_eff goes negative (TYPE INVERSION from n-type to p-type)

This is a dramatic device-level change that cascades through everything:

- Depletion width at -30V changes from ~10 um (pristine) to the full epi thickness (after partial compensation) to a new junction forming from the cathode side (after type inversion)
- C-V characteristic shifts: capacitance flattens as full depletion is reached at lower voltage
- Electric field profile changes shape: from triangular (one-sided junction) to trapezoidal to inverted
- The Hecht equation validation (CCE vs bias) no longer applies because the field profile assumption changes

The critical fluence for the Petringa detector is remarkably LOW: ~1.6e13 p/cm^2. This is because the epi doping (8.5e13) is very low by SiC standards.

**Why it happens:** The Petringa detector was designed for maximum sensitivity (low doping = wide depletion = more signal), not radiation hardness. Low initial doping means carrier removal reaches full compensation at lower fluence. Users who expect SiC to be "radiation hard" may not check whether their specific device geometry hits full compensation within the fluence range of interest.

**Consequences:**

1. **Numerical:** The solver may fail when N_eff crosses zero because the built-in potential and depletion approximations break down. The graded doping profile makes this worse: N_D_junction=2.9e15 is barely affected while N_D_bulk=8.5e13 is fully compensated, creating a non-monotonic effective doping profile.
2. **Physical:** CCE predictions become unreliable if the simulation does not properly handle the transition from partially depleted to fully depleted to type-inverted regimes.
3. **Validation:** Published CCE vs fluence curves for OTHER SiC detectors (different doping) cannot be directly compared.

**Prevention:**

1. Compute the critical fluence Phi_crit = N_D_bulk / c_removal = 8.5e13 / 5.3 = 1.6e13 p/cm^2 BEFORE running any simulations. This sets the fluence scale for your study.
2. Handle the N_eff = 0 crossing explicitly. Below Phi_crit: standard depletion from junction. Above Phi_crit: the epi bulk is intrinsic or p-type, requiring different solver initialization.
3. For the graded profile: compute N_eff(x, Phi) at each position. The junction-side doping (2.9e15) survives to much higher fluence (Phi_crit_junction = 2.9e15/5.3 = 5.5e14), so the profile becomes increasingly non-uniform.
4. Monitor the Newton solver: if iterations increase sharply at a particular fluence, you are likely near the compensation point. Add finer fluence steps around Phi_crit.

**Detection:** Plot N_eff(x) at each fluence. When the minimum of N_eff(x) approaches zero, you are near the critical regime. If the solver diverges, this is almost certainly the cause.

**Phase:** Must be addressed in the carrier removal implementation phase. Compute Phi_crit for the Petringa device geometry as the FIRST step.

---

### Pitfall 5: Dark Current Model Breaks Under Radiation Damage

**What goes wrong:** The v1.1 dark current model in `dark_current.py` was calibrated to match 18 pA at pristine conditions using an effective generation rate N_t = 2.2e13 cm^-3/s with Hurkx field enhancement. Under radiation damage, dark current INCREASES dramatically because:

1. Radiation-introduced defects (Z1/2, EH6/7) create additional generation-recombination centers in the depletion region
2. The depletion width changes (Pitfall 4), changing the volume of the generation region
3. The electric field profile changes, modifying the Hurkx field-enhancement factor

The calibrated N_t = 2.2e13 is an EFFECTIVE parameter that lumps together all generation mechanisms. When you add radiation-induced defects as SEPARATE SRH centers, you are double-counting: the pristine dark current already includes generation from the pre-existing Z1/2 centers (typical as-grown density ~10^12 cm^-3 in good epi).

If you add radiation-induced Z1/2 with density g_Z1/2 \* Phi ON TOP of the calibrated effective N_t, the fluence=0 dark current will be too high (double-counting the pre-existing Z1/2).

**Why it happens:** The v1.1 dark current model was designed as an effective (phenomenological) model, not a physically decomposed model. It cannot cleanly separate pre-existing defects from radiation-induced defects.

**Consequences:** Either the fluence=0 dark current is wrong (if you add radiation defects on top), or the fluence-dependent dark current increase is wrong (if you suppress the radiation contribution to avoid double-counting).

**Prevention:**

1. **Option A (recommended for v2.0):** Keep the v1.1 dark current model untouched at fluence=0. For the radiation-induced dark current INCREASE, add only the ADDITIONAL generation rate from radiation-introduced defects: delta_J_dark(Phi) = q _ integral[sigma_n _ v_th _ N_Z12_rad(Phi) _ n_i \* W(Phi)]. This gives J_dark(Phi) = J_dark(0) + delta_J_dark(Phi), preserving the calibrated baseline.
2. **Option B (for later):** Rebuild the dark current model from scratch with explicit defect populations: N_Z12_total = N_Z12_pristine + g_Z12 \* Phi. This requires re-calibrating the pristine model with a physical (not effective) Z1/2 density, which is a bigger refactoring.
3. Do NOT mix the effective N_t model with explicit radiation-induced defect populations in the same SRH expression.

**Detection:** After adding radiation damage, check J_dark at fluence=0. It must still be ~18 pA (within the tolerance of the v1.1 calibration). If it doubled, you are double-counting.

**Phase:** Must be designed in the architecture phase, before implementing either carrier removal or dark current increase features.

---

## Moderate Pitfalls

### Pitfall 6: Newton Solver Convergence Failure at High Damage Levels

**What goes wrong:** As fluence increases, multiple parameters change simultaneously (lifetime decreases, doping decreases, trap density increases). The Newton-Raphson solver in devsim can fail to converge because:

1. **Lifetime approaching zero:** If 1/tau = 1/tau_0 + K_tau\*Phi, at high Phi the lifetime becomes very short. When tau << transit time, the SRH recombination term dominates the continuity equations, making them stiff. The Jacobian condition number worsens.

2. **Doping compensation near zero:** When N_eff approaches zero, the device transitions from a pn junction to a pi junction to a pp+ junction. The equilibrium carrier concentrations change by orders of magnitude. The initial guess (from the pristine or previous-fluence solution) may be far from the new solution.

3. **Multiple defect levels:** Each SRH trap level adds terms to the generation-recombination rate. With 3 defect species (Z1/2 + EH6/7 + EH4), the recombination model has 3 sets of (n1, p1, tau_n, tau_p) parameters. The Jacobian derivatives become complex, and any error in the analytical derivatives causes Newton to diverge.

4. **devsim-specific:** The GTS TCAD team (Burin 2024) explicitly reported "multiple convergence issues" with 4H-SiC radiation damage simulations, requiring "disabling incomplete ionization" and "fine-tuning simulation grids separately for forward and reverse bias." These same issues will appear in devsim.

**Why it happens:** The pristine SiC device is numerically well-behaved because the physics is dominated by a simple pn junction with well-separated carrier concentrations. Radiation damage destroys this nice structure: low effective doping means p and n concentrations can be comparable throughout the device.

**Prevention:**

1. **Fluence ramping:** Do NOT jump from fluence=0 to fluence=1e14. Use the previous-fluence solution as the initial condition for the next fluence. Step in factors of 2-3x: 0, 1e11, 3e11, 1e12, 3e12, 1e13, 3e13, 1e14.
2. **Tight convergence at low fluence, relaxed at high fluence:** Set abs_error and rel_error tolerance tighter for low fluence (where you need accuracy) and relax slightly at high fluence (where you need convergence).
3. **Monitor Newton iterations:** If iterations exceed 50 at any fluence step, halve the fluence increment. If iterations exceed 100, the Jacobian derivatives are likely wrong.
4. **Keep devsim's incomplete ionization active:** The GTS workaround of disabling incomplete ionization is specific to their solver. In devsim, incomplete ionization is handled at device creation time (pre-computed ionized fractions), not in the solver loop, so it should not cause convergence issues.
5. **Start with a single effective defect level** before adding multiple trap levels. Debug convergence with the simplest model first.

**Detection:** Log Newton iteration count at each fluence step. A sudden increase (e.g., from 10 to 80 iterations) signals approaching a convergence boundary. Plot iteration count vs fluence to identify problematic regimes.

**Phase:** Must be handled during implementation. Build iteration monitoring into the damage simulation loop from the start.

---

### Pitfall 7: Confusing Defect Introduction Rate with Carrier Removal Rate

**What goes wrong:** Two different quantities are used in the literature, often with the same symbol or similar names:

1. **Defect introduction rate** g [cm^-1]: Number of defects per unit fluence. N_defect = g \* Phi. Measured by DLTS. Example: g_Z1/2 = 5.0 cm^-1 means 5 Z1/2 centers per neutron equivalent per cm^3.

2. **Carrier removal rate** c [cm^-1]: Change in effective free carrier concentration per unit fluence. delta_N_eff = c \* Phi. Measured by C-V. Example: c = 5.3 cm^-1 means each proton/cm^2 removes 5.3 free electrons/cm^3.

These are NOT the same quantity. A single Z1/2 center (doubly-charged acceptor, negative-U) removes TWO electrons from the conduction band. Other defects may remove 0, 1, or 2 carriers depending on their charge state and position in the bandgap. The carrier removal rate is the NET effect of ALL defects combined:

c = sum_i(charge_state_i \* g_i)

If you use the Z1/2 introduction rate (5.0) as the carrier removal rate, you may underestimate carrier removal (because you miss contributions from other defects) OR overestimate it (because not all Z1/2 centers are doubly charged).

**Why it happens:** Papers sometimes use "removal rate" loosely to mean either quantity. The symbols g, c, k, beta, and alpha are all used inconsistently across the literature.

**Consequences:** Using g instead of c (or vice versa) in the effective doping formula gives the wrong critical fluence for full compensation, the wrong depletion width vs fluence, and the wrong CCE prediction.

**Prevention:**

1. In the code, name variables explicitly: `g_Z12_introduction` vs `c_carrier_removal`. Never use ambiguous names like `damage_rate`.
2. Use the carrier removal rate c (from C-V measurements) for effective doping: N_eff = N_D - c*Phi. Use introduction rates g (from DLTS) for trap densities: N_Z12 = g_Z12 * Phi.
3. Do NOT derive c from g unless you have a validated model for the charge states of all defects. Use directly measured c values from recent proton irradiation studies (4.2-6.4 cm^-1 for ~250 MeV protons).

**Detection:** If your simulated C-V shift with fluence disagrees with published C-V data by a factor of 2-5, check whether you confused g and c.

**Phase:** Must be clarified in the damage parameters definition phase. Add comments in the DamageParameters dataclass distinguishing the two concepts.

---

### Pitfall 8: Assuming Linear Lifetime Degradation (1/tau = 1/tau_0 + K\*Phi)

**What goes wrong:** The standard model for silicon (and most TCAD textbooks) uses a linear relationship between reciprocal lifetime and fluence:

1/tau(Phi) = 1/tau_0 + K_tau \* Phi

Recent 4H-SiC data (Luo et al. 2025, 80 MeV protons) finds a LOGARITHMIC dependence instead:

1/tau = a \* ln(Phi) + b

The physical reason is not fully understood but may relate to defect clustering: at high fluence, new vacancies form clusters with existing vacancies rather than creating independent recombination centers, giving diminishing returns per additional displacement. Hazdra et al. (2021) similarly found non-linear lifetime behavior in 4H-SiC.

The practical difference is large at high fluence:

- Linear model at Phi=1e14: 1/tau = 1/tau_0 + K\*1e14 (very short lifetime)
- Logarithmic model at Phi=1e14: 1/tau = a\*32.2 + b (much less degraded)

The logarithmic model predicts the device remains functional to higher fluence than the linear model.

**Why it happens:** The linear model is ingrained in the Si community and all TCAD textbooks. It is the default assumption.

**Consequences:** Overpredicting CCE degradation at high fluence (linear model predicts more damage than actually occurs). This matters for radiation hardness claims -- you might conclude the device fails at 1e13 p/cm^2 when it actually survives to 1e14.

**Prevention:**

1. Implement both models behind a parameter flag:
   ```python
   class DamageParameters:
       lifetime_model: str = "linear"  # "linear" or "logarithmic"
       K_tau_n: float = ...  # for linear model
       a_tau: float = 2.4e4  # for logarithmic model (Luo 2025)
       b_tau: float = 1.9e6  # for logarithmic model
   ```
2. Show both predictions in the CCE vs fluence plot with a band showing the difference.
3. For the Petringa detector (62 MeV protons), neither model is directly validated. State this uncertainty explicitly.

**Detection:** If your predicted lifetime at Phi=1e14 is < 1 ps, suspect the linear model is overshooting. Physical SiC lifetimes at 1e14 n_eq/cm^2 are typically ~1-10 ns based on experimental data.

**Phase:** Implement in the lifetime degradation phase. LOW confidence on which model is correct for the Petringa conditions.

---

### Pitfall 9: Ignoring Position-Dependence of Damage in Thin Epi Layers

**What goes wrong:** The standard radiation damage model assumes UNIFORM damage throughout the device: defect density = g \* Phi everywhere. This is valid when the irradiating particles pass through the device without significant energy loss (as for high-energy protons through a 10-um epi layer). But there are two cases where position-dependence matters:

1. **Low-energy protons stopping in the device.** The Petringa detector uses 62 MeV protons for dosimetry. At this energy, the Bragg peak is at ~30 mm in water (~15 mm in SiC), far beyond the 10-um device. So for the primary beam, uniform damage is correct. BUT if the detector is used at the Bragg peak (for microdosimetry), the damage profile is highly non-uniform.

2. **Cumulative damage from the clinical beam itself.** The detector accumulates damage as it measures dose. At FLASH dose rates (20-230 Gy/s), the fluence rate is ~10^9-10^10 p/cm^2/s. Over a typical 100-session clinical use, total accumulated fluence may reach ~10^12-10^13 p/cm^2. This is within the range where carrier removal becomes significant.

3. **The graded doping profile makes position matter even for uniform damage.** With N_D(x) varying from 2.9e15 to 8.5e13, uniform carrier removal (c\*Phi everywhere) has a vastly different fractional effect: 0.2% reduction at the junction vs 60% in the bulk at Phi=1e13.

**Why it happens:** The uniform damage assumption simplifies the model. But the graded doping profile of the Petringa detector means even uniform damage has strongly position-dependent EFFECTS.

**Prevention:**

1. For 62 MeV protons passing through the 10-um epi: uniform damage is justified. Document this assumption with a SRIM/TRIM calculation showing the NIEL is approximately constant across 10 um at this energy.
2. When computing N_eff(x, Phi) = N_D(x) - c\*Phi, this is already position-dependent through N_D(x). Ensure the implementation uses the node model (position-dependent) not a scalar parameter (position-independent).
3. For future extension to Bragg-peak irradiation: the damage profile must be computed from the NIEL(x) distribution, not assumed uniform.

**Detection:** Plot N_eff(x) at several fluences. If N_eff goes negative in the bulk while remaining positive near the junction, you have correctly captured the position-dependent effect of damage on the graded profile.

**Phase:** Address in the carrier removal implementation. Use node models, not scalar parameters.

---

### Pitfall 10: Validating Against Wrong Experimental Data

**What goes wrong:** CCE vs fluence data for SiC detectors exists from several groups, but the experimental conditions vary enormously:

| Group             | Device    | Doping     | Epi thickness | Particle | Fluence range      |
| ----------------- | --------- | ---------- | ------------- | -------- | ------------------ |
| Burin (CERN RD50) | 50 um pad | ~10^14     | 50 um         | neutron  | 5e14-1e16 neq/cm^2 |
| Luo (CSNS)        | PIN       | 5.2e13     | 100 um        | 80 MeV p | 1e11-1e14 neq/cm^2 |
| Petringa (LNS)    | PIN       | 8.5e13     | 10 um         | 62 MeV p | ??? (unpublished)  |
| SiC LGAD (2025)   | LGAD      | gain layer | thin          | proton   | 1e13-3e14 p/cm^2   |

Each of these has different critical fluence, different CCE degradation curves, and different damage constants. You CANNOT validate the Petringa 10-um detector simulation against the Burin 50-um detector data without scaling corrections for:

- Different epi thickness (affects charge collection geometry)
- Different doping (affects critical fluence and depletion width)
- Different irradiation particle/energy (affects NIEL and defect spectrum)
- Different operating bias (affects field profile and CCE)

**Why it happens:** Temptation to validate against the most accessible published data, which may not match the Petringa device.

**Consequences:** Apparent agreement with the wrong experimental data gives false confidence. Disagreement with the wrong data triggers unnecessary model tuning.

**Prevention:**

1. Clearly document which published data your model is compared against and why (or why not).
2. When comparing to data from different devices, normalize CCE to the device-specific parameters: CCE(V/V_fd, Phi/Phi_crit) where V_fd is full depletion voltage and Phi_crit is the critical fluence for that specific doping.
3. The IDEAL validation is against Petringa group experimental data (if/when available). Until then, present model predictions as PREDICTIONS, not validated results.
4. The Luo et al. data (N_D=5.2e13, 100 um, proton) is the closest match to the Petringa device in terms of doping level. Use this for qualitative comparison but note the 10x thickness difference.

**Detection:** If your model matches one dataset perfectly but disagrees with another by 2x, you are likely over-fitting to the first dataset's specific conditions.

**Phase:** Must be addressed when writing the validation notebooks. Define the comparison strategy upfront.

---

## Minor Pitfalls

### Pitfall 11: Annealing Model Using Si Activation Energies

**What goes wrong:** Si radiation damage anneals significantly at room temperature (beneficial annealing timescale ~days at 20C, activation energy ~0.3-0.5 eV). Implementing the same model for SiC gives wrong annealing predictions because:

- Z1/2 in SiC is thermally stable to ~1500C (activation energy for migration ~3-4 eV)
- Radiation-induced S1, S2 centers anneal at 300-425C (still far above room temperature)
- Thermally unstable defects (E1, E2, E3) transform into stable defects (Z1/2, S1, S2) at temperatures up to ~200C

At the Petringa detector's operating temperature (30-40C), effectively NO annealing of the dominant Z1/2 defect occurs. Room-temperature annealing applies ONLY to the unstable intermediate defects, and even that saturates quickly.

**Prevention:**

1. For room-temperature operation: do NOT implement Z1/2 annealing. It does not occur.
2. For thermal annealing studies (elevated temperature): use activation energies specific to each SiC defect species, not Si-derived values. Z1/2 migration: Ea > 3 eV. EH6/7 similarly stable. Only shallow defects and interstitial-related defects anneal at moderate temperatures.
3. If implementing annealing at all, make it optional and off by default for room-temperature simulations.

**Phase:** Annealing implementation phase. Flag as LOW priority for the initial damage model.

---

### Pitfall 12: Mesh Resolution Insufficient for Damage-Induced Field Spikes

**What goes wrong:** When carrier removal creates non-uniform effective doping (especially with the graded profile), the electric field can develop sharp features at the boundary between compensated and uncompensated regions. The existing mesh was designed for the pristine graded profile and may not resolve these features.

**Prevention:**

1. Before running damage simulations, verify the mesh resolves the region where N_eff transitions from positive to near-zero. This may require adding mesh points in the epi bulk region.
2. A simple check: run at a fluence just below Phi_crit, extract E(x), and verify it is smooth. If E(x) shows oscillations or discontinuities, the mesh is too coarse.

**Phase:** Address when implementing carrier removal. May need a `refine_mesh_for_damage()` utility.

---

### Pitfall 13: CCE Computation Method Incompatible with Damaged Device

**What goes wrong:** The v1.1 CCE computation in `charge_collection.py` uses the Hecht equation for analytical benchmarking (CCE = mu*tau*E/d^2 * [1-exp(-d^2/(mu*tau\*E))]). Under radiation damage, this equation requires updating:

- tau becomes tau(Phi), much shorter
- E becomes E(x, Phi), non-uniform due to doping change
- d (collection distance) changes because depletion width changes

But the Hecht equation assumes a UNIFORM electric field and a SINGLE carrier type. These assumptions break down more severely in the damaged device than in the pristine device.

**Prevention:**

1. Use the drift-diffusion CCE (from the devsim solver) as the primary result. The Hecht equation becomes a rough cross-check, not a validation target.
2. Update the Hecht equation parameters (tau, mu, E, d) consistently with the damage model, but present it clearly as an approximation.
3. For severely damaged devices (Phi > Phi_crit), the Hecht equation is qualitatively wrong. Do not use it.

**Phase:** Address when implementing CCE vs fluence notebooks.

---

### Pitfall 14: Temperature-Damage Cross-Coupling Ignored

**What goes wrong:** The v1.1 temperature model and the v2.0 damage model are developed independently. But there are cross-couplings:

- SRH lifetime depends on BOTH T and Phi: tau(T, Phi), not tau(T) \* f(Phi)
- The damage constant K_tau itself may depend on temperature (Arrhenius-type, as found by the IEEE Access 2024 study)
- Carrier removal may be temperature-dependent (defect charge states change with T)

If you implement damage at T=300K only and then try to predict CCE(T, Phi), you may miss the T-dependence of the damage itself.

**Prevention:**

1. For the initial v2.0 implementation: damage at T=300K only is acceptable. Document the limitation.
2. For later extension: the lifetime degradation should be parameterized as tau(T, Phi) = 1/[1/tau_0(T) + K_tau(T)*Phi] where K_tau(T) has its own T-dependence.
3. This cross-coupling is a research question, not a solved problem. Flag for future phases.

**Phase:** Defer to after the basic damage model is validated at 300K.

---

## Phase-Specific Warnings

| Phase Topic                      | Likely Pitfall                                                                        | Mitigation                                                               | Severity |
| -------------------------------- | ------------------------------------------------------------------------------------- | ------------------------------------------------------------------------ | -------- |
| Damage parameter definition      | Pitfall 2 (wrong constants), Pitfall 7 (g vs c confusion)                             | Explicit DamageParameters dataclass with provenance, separate g and c    | CRITICAL |
| Defect model selection           | Pitfall 3 (SiC != Si), Pitfall 8 (linear vs log lifetime)                             | Use SiC-specific defects from Burin 2024; implement both lifetime models | CRITICAL |
| Carrier removal implementation   | Pitfall 4 (operating point shift), Pitfall 9 (position-dependence), Pitfall 12 (mesh) | Compute Phi_crit first; use node models; check mesh resolution           | CRITICAL |
| Fluence=0 regression             | Pitfall 1 (baseline break), Pitfall 5 (dark current double-counting)                  | Additive damage API; mandatory regression test; separate delta_J_dark    | CRITICAL |
| Solver stability at high fluence | Pitfall 6 (Newton convergence)                                                        | Fluence ramping; iteration monitoring; start with single defect          | MODERATE |
| Dark current vs fluence          | Pitfall 5 (double-counting)                                                           | Additive delta_J model; preserve v1.1 calibration at Phi=0               | MODERATE |
| CCE validation                   | Pitfall 10 (wrong comparison data), Pitfall 13 (Hecht breakdown)                      | Use device-matched data; rely on DD solver not Hecht                     | MODERATE |
| Annealing modeling               | Pitfall 11 (Si activation energies)                                                   | SiC-specific Ea; Z1/2 does not anneal at room T                          | MINOR    |
| T-Phi cross-coupling             | Pitfall 14 (ignored coupling)                                                         | Defer; document limitation; 300K only for initial model                  | MINOR    |

---

## Sources

- [Burin et al. (2024) - TCAD modeling of radiation-induced defects in 4H-SiC diodes](https://arxiv.org/abs/2407.11776) -- Three-defect model (Z1/2 + EH6/7 + EH4), introduction rates, convergence issues with GTS TCAD. MEDIUM confidence (neutron irradiation, not proton; 50 um device, not 10 um).
- [Burin et al. (2025) - TCAD Simulations of Radiation Damage in 4H-SiC](https://arxiv.org/html/2407.16710v1) -- GTS framework convergence problems, negative capacitance at low fluence, literature parameter spread. MEDIUM confidence.
- [Luo et al. (2025) - Mechanisms of proton irradiation-induced defects in 4H-SiC PIN detectors](https://arxiv.org/html/2503.09016) -- Logarithmic lifetime model, EH3 introduction rate 1.48 cm^-1, Z1/2 NOT dominant under 80 MeV proton. MEDIUM confidence (different proton energy from Petringa).
- [In-situ radiation damage study of SiC detectors with clinical proton beams (2025)](https://arxiv.org/abs/2510.11304v2) -- Carrier removal rate 4.2-6.4 cm^-1 for ~253 MeV protons. MEDIUM confidence (closest to Petringa conditions but different energy).
- [IEEE Access (2024) 10538275 - Carrier Lifetime Dependence on Temperature and Proton Irradiation in 4H-SiC](https://ieeexplore.ieee.org/document/10538275) -- Arrhenius T-dependence of damage constant K_tau, power-law lifetime-temperature relation with fluence-dependent exponent. MEDIUM confidence.
- [Carrier removal rates in 4H-SiC power diodes (2023)](https://www.sciencedirect.com/science/article/abs/pii/S136980012300464X) -- Predictive analytical model for carrier removal using NIEL. MEDIUM confidence (power diodes, not detectors).
- [Frontiers in Physics (2022) - SiC detectors review](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2022.898833/full) -- Z1/2 as dominant recombination center, general SiC detector properties. HIGH confidence (review article).
- [Hazdra et al. (2021) - Radiation Defects and Carrier Lifetime in 4H-SiC Bipolar Devices](https://onlinelibrary.wiley.com/doi/abs/10.1002/pssa.202100218) -- Non-linear lifetime degradation in SiC. MEDIUM confidence.
- [Ioffe NSM 4H-SiC archive](https://www.ioffe.ru/SVA/NSM/Semicond/SiC/) -- Reference material parameters. HIGH confidence.
- [Frontiers in Physics (2021) - TCAD Modeling of Surface Radiation Damage](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2021.617322/full) -- General TCAD radiation damage modeling practices, effective trap models. MEDIUM confidence (Si-focused but methodology applicable).

---

_Pitfalls research for: v2.0 milestone -- adding radiation damage modeling to validated 1D 4H-SiC TCAD simulator_
_Researched: 2026-03-24_
