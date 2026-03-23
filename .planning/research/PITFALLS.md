# Domain Pitfalls

**Domain:** Adding temperature dependence, surface physics, and transient dynamics to a validated 1D 4H-SiC TCAD simulator (v1.1 milestone)
**Researched:** 2026-03-23
**Confidence:** HIGH for codebase-specific pitfalls (verified against actual source), MEDIUM for physics model pitfalls (literature-sourced)

## Critical Pitfalls

Mistakes that cause rewrites, silent physics errors, or loss of validated results.

### Pitfall 1: Hardcoded 300K Parameters Scattered Throughout the Codebase

**What goes wrong:** The existing codebase has T=300K baked into multiple locations that are easy to miss when threading temperature through. `SiC4H_Parameters` stores `n_i_300`, `NC_300`, `NV_300`, `mu_n_max`, `mu_p_max` as fixed 300K constants. `device.py` uses `params.n_i_300` directly when setting the devsim `n_i` region parameter and sets `params.tau_n`/`params.tau_p` as fixed SRH lifetimes. The `charge_collection.py` module instantiates `_params = SiC4H_Parameters()` at module level and uses fixed 300K mobilities as function defaults (e.g., `mu_e=_params.mu_n_max` in the Hecht equation). The clamped Boltzmann formulation in `poisson.py` uses `V_t` computed once at device creation.

If you change T in `create_sic_device()` without updating every downstream consumer, you get a device where some physics runs at the new T and some silently remains at 300K.

**Why it happens:** The v1.0 code was correctly designed for single-temperature operation. Temperature was a parameter passed to `create_sic_device()` but never exercised at non-300K values. The `compute_ni(T)` function exists in `sic_material.py` but is explicitly documented as "Not used in the v1.0 pipeline."

**Consequences:** Silent physics errors. The simulator appears to work at T!=300K but produces wrong results because n_i, mobility, SRH lifetimes, or Hecht equation benchmarks use stale 300K values. These errors are subtle -- the simulator converges fine, and results look plausible but are quantitatively wrong. Worst case: you validate against 300K data, conclude T-dependence works, but the T-dependent models were never actually active.

**Prevention:**

1. Create a `compute_all_parameters(T, N_D)` function that returns ALL T-dependent quantities: n_i(T), NC(T), NV(T), mu_n(T, N_D), mu_p(T, N_A), E_g(T), tau_n(T), tau_p(T), V_t(T). Wire this into `create_sic_device()`.
2. Audit every usage of `params.n_i_300`, `params.NC_300`, `params.NV_300`, `params.mu_n_max`, `params.mu_p_max`, `params.tau_n`, `params.tau_p`. Specific locations to fix:
   - `device.py` line 153: `value=params.n_i_300` -- must become T-dependent
   - `device.py` lines 196-199: `n1=params.n_i_300`, `p1=params.n_i_300` -- must track n_i(T)
   - `device.py` lines 176-177: `mobility_caughey_thomas(N_D)` -- must accept T
   - `charge_collection.py` line 31: `_params = SiC4H_Parameters()` -- module-level instantiation bakes in 300K at import time
   - `charge_collection.py` lines 37-42: Hecht equation defaults use fixed 300K mobilities
3. After implementing, run the full test suite at T=300K and verify results are BIT-IDENTICAL to v1.0 output. Any deviation means you broke the parameter threading.

**Detection:** Run at T=300K and T=310K. If C-V, CCE, and I-V results are identical, T-dependence is not wired in. For SiC in the 300-313K clinical range, n_i changes by roughly an order of magnitude per 10K (because E_g/kT is so large at 3.26 eV), so even small T changes must produce measurable shifts.

---

### Pitfall 2: Using Silicon Temperature Exponents for 4H-SiC Mobility

**What goes wrong:** Applying Si-derived temperature exponents for mobility to 4H-SiC. The physics is different and the existing `mobility_caughey_thomas()` function has no T-dependence at all -- it returns constant values regardless of T.

The correct 4H-SiC temperature-dependent Caughey-Thomas model (TU Wien Ayalew thesis, Table 3.5):

- mu_max(T) = mu_max(300K) \* (T/300)^gamma, where gamma_n = -2.40 (electrons), gamma_p = -2.15 (holes)
- mu_min(T) = mu_min(300K) \* (T/300)^beta, where beta = -0.5 for both carriers

For Si, the gamma exponents are different (-2.42 for electrons, -2.20 for holes in the Arora model). The difference for holes (gamma_p = -2.15 vs -2.20) matters less than getting the model right from the start.

More importantly: several published SiC TCAD papers use gamma values from the wrong polytype (6H-SiC has gamma_n = -2.70, quite different from 4H-SiC's -2.40). This polytype confusion is well-documented in the literature.

**Why it happens:** Most TCAD examples and online resources are Si-centric. The current `mobility_caughey_thomas()` function already uses correct 4H-SiC values at 300K but simply has no T parameter. Adding T-dependence by searching "Caughey-Thomas temperature model" will return Si exponents.

**Consequences:** At T=313K (40C clinical): mu_n should decrease by ~3% relative to 300K. Using Si exponents gives ~3.1% -- small difference. But using 6H-SiC exponents gives ~4% decrease. The error compounds through drift-diffusion current, transit time, and CCE computation.

**Prevention:**

1. Use ONLY the TU Wien Ayalew thesis values (already cited in `sic_material.py`):
   - gamma_n = -2.40, gamma_p = -2.15 (mu_max exponents)
   - beta_n = beta_p = -0.5 (mu_min exponents)
2. Add these as constants to `SiC4H_Parameters`:
   ```python
   gamma_mu_n: float = -2.40
   gamma_mu_p: float = -2.15
   beta_mu: float = -0.5
   ```
3. Extend `mobility_caughey_thomas()` to accept T:
   ```python
   def mobility_caughey_thomas(N_total, carrier="electron", T=300):
       mu_max_T = mu_max_300 * (T/300)**gamma
       mu_min_T = mu_min_300 * (T/300)**beta
       ...
   ```

**Detection:** Verify mu_n(300K, 1e14) = 949 cm^2/Vs (should match v1.0). Verify mu_n(400K, 1e14) ~ 655 cm^2/Vs (decreases ~31%). If mu_n(400K) is within 5% of mu_n(300K), the T-dependence is not working.

---

### Pitfall 3: Breaking the Validated 300K Baseline (Regression)

**What goes wrong:** Adding new physics (temperature models, surface terms, transient capability) inadvertently changes the 300K steady-state results that are currently validated (C-V R^2=0.998, CCE=100% at V>-40V, flat CCE across FLASH dose rates). Common ways this breaks:

1. Replacing fixed `n_i_300 = 5e-9` with `compute_ni(300)` which returns a slightly different value. Current `compute_ni(300)` returns ~6.5e-9 (verified by reading the code: it uses m_e=0.77*m0, m_h=1.0*m0, M_c=3, Varshni with the same parameters). The 30% discrepancy from the literature value of 5e-9 would shift all I-V and C-V results.
2. Adding surface recombination boundary conditions that are active even at 300K, changing the I-V dark current from the validated ideal-SRH floor.
3. Modifying `ElectronGeneration`/`HoleGeneration` model expressions to include new terms, changing the devsim equation registration.
4. Changing the mesh to accommodate surface physics or transient resolution.

**Why it happens:** The v1.0 results depend on exact numerical values and model registrations. The code has no automated regression tests -- validation is done manually via 5 Jupyter notebooks.

**Consequences:** Loss of the validated baseline. You can no longer claim R^2=0.998 for C-V or flat CCE across FLASH dose rates. You may not know whether deviations come from new physics or from accidentally changed old physics.

**Prevention:**

1. **Before ANY code changes**, extract golden reference values from v1.0 notebooks and store them:
   - C-V capacitances at V = 0, -10, -30V
   - I-V dark current at V = -30, -60V
   - CCE at V = -10, -30, -40V for alpha particles
   - CCE vs dose rate at 20, 100, 230 Gy/s
2. Implement a `test_300K_regression()` function that creates a device at T=300K with v1.0 parameters and verifies results match golden values within tolerance.
3. New physics must be **additive and toggleable**: all new models off by default, explicitly enabled via parameters. When all new physics is disabled, results MUST match v1.0.
4. The `compute_ni(300)` discrepancy (6.5e-9 vs 5e-9) must be resolved BEFORE using it: either calibrate the DOS masses to reproduce 5e-9, or accept the new value and re-validate. Do NOT silently switch.

**Detection:** Run the 300K C-V simulation after every code change. If R^2 drops below 0.995 or CCE changes by more than 0.1%, something broke.

---

### Pitfall 4: Surface Recombination Cannot Explain 18 pA Dark Current in 1D

**What goes wrong:** The goal is to match the experimental 18 pA dark current. The natural approach is to add surface recombination. But a quantitative check shows surface SRH generation cannot produce 18 pA:

Surface generation current: J_surf = q _ n_i _ S

- At 300K with n_i = 5e-9 cm^-3 and S = 1000 cm/s (typical SiO2-passivated SiC):
- J_surf = 1.6e-19 _ 5e-9 _ 1000 = 8e-25 A/cm^2
- For 4 mm^2 area: I_surf = 3.2e-25 A -- 16 orders of magnitude below 18 pA

Space-charge-region generation: J_gen = q _ n_i _ W / (2 \* tau_eff)

- With W = 10e-4 cm, tau_eff = max(tau_n, tau_p) = 6e-7 s:
- J_gen = 1.6e-19 _ 5e-9 _ 10e-4 / (2 \* 6e-7) = 6.7e-28 A/cm^2
- Also negligibly small

The experimental 18 pA over 4 mm^2 = 4.5e-9 A/cm^2. This is ~19 orders of magnitude above the intrinsic SRH generation. The dark current CANNOT be explained by any bulk or surface SRH mechanism at 300K with n_i = 5e-9.

**Why it happens:** SiC's extremely low n_i makes all n_i-proportional currents negligible. The experimental dark current must come from non-n_i mechanisms: perimeter/edge leakage (2D/3D geometry effect), trap-assisted tunneling (field-enhanced generation through deep levels), or measurement artifacts (cable leakage, probe station noise).

**Consequences:** If you spend time implementing a sophisticated surface recombination model to match 18 pA, you will fail. No physically reasonable SRV value can bridge 16 orders of magnitude. You will end up using unphysically large SRV values that have no physical meaning.

**Prevention:**

1. Accept that matching 18 pA in a 1D simulator requires an **effective generation mechanism**, not a physical surface recombination model.
2. The most likely physical mechanisms for 18 pA are:
   - Perimeter leakage: current flows along the SiO2/SiC interface around the device perimeter. This is proportional to perimeter, not area. For a 2x2 mm device, perimeter/area = 8mm/4mm^2 = 2 cm^-1. This is inherently a 2D effect.
   - Trap-assisted tunneling (TAT) through the Z1/2 center (E_C - 0.65 eV): field-enhanced generation in the depletion region. This CAN be modeled in 1D as a field-dependent generation rate.
   - Generation-recombination through deep levels with temperature-activated capture cross-sections.
3. For the 1D model, implement TAT as the primary dark current mechanism. Use the Hurkx model or a simplified field-enhanced generation rate. The Z1/2 trap level at E_C - 0.65 eV is well-documented for 4H-SiC.
4. Treat the TAT parameters (trap density, tunneling mass) as fitting parameters to match 18 pA, and document this clearly.

**Detection:** Before implementing ANY dark current model, compute the theoretical maximum current from each mechanism at the relevant n_i value. If the mechanism cannot produce current within 3 orders of magnitude of the target, it is not the dominant mechanism.

---

### Pitfall 5: Trap-Assisted Tunneling Model Complexity Explosion

**What goes wrong:** 4H-SiC has multiple well-known defect levels that contribute to TAT:

- Z1/2 center at E_C - 0.65 eV (carbon vacancy, dominant in n-type epi, density ~10^12-10^14 cm^-3)
- EH6/7 at E_C - 1.55 eV
- Shallow N donors at E_C - 0.05 eV and E_C - 0.09 eV
- Interface traps at SiC/SiO2 with continuous Dit distribution

Recent research (2024-2025) shows that choosing wrong trap combinations either underestimates or overestimates reverse I-V by orders of magnitude. The full Hurkx TAT model has 5+ free parameters per trap level. With 3-4 trap levels, that is 15-20 parameters to fit against a single I-V curve -- massively under-constrained.

**Why it happens:** The natural instinct is to be physically complete: include all known trap levels, implement full phonon-assisted tunneling with WKB approximation, add field-enhanced emission. But devsim has no built-in TAT model -- everything must be implemented from scratch using custom equation strings, and each additional mechanism multiplies the Jacobian complexity and convergence difficulty.

**Consequences:** Weeks spent implementing a multi-trap TAT model that either (a) does not converge because the Jacobian derivatives are wrong or incomplete, (b) converges but has too many free parameters to provide meaningful fits, or (c) matches 18 pA but only because of parameter overfitting, not physics.

**Prevention:**

1. Start with the **simplest possible model**: a single effective generation rate in the depletion region proportional to exp(E/E0), where E is the local electric field and E0 is a fitting parameter. This captures the essential field-enhancement of TAT without the full Hurkx complexity.
2. The model has exactly 2 fitting parameters: G0 (prefactor) and E0 (field scale). Fit these to match 18 pA at -60V and the voltage dependence of the dark current.
3. Only add complexity (specific trap levels, proper WKB tunneling) if the simple model fails to reproduce the I-V shape across the full 0 to -60V range.
4. All Jacobian derivatives must be analytically correct for Newton convergence. For a field-dependent generation rate, this means derivatives w.r.t. Potential (through the E-field), Electrons, and Holes.

**Detection:** If the Newton solver for the dark current simulation takes > 100 iterations or fails to converge, suspect missing or incorrect Jacobian derivatives. Verify by comparing analytical derivatives against numerical finite-difference derivatives.

---

## Moderate Pitfalls

### Pitfall 6: Transient Timestep Selection for Multi-Scale Dynamics

**What goes wrong:** The FLASH pulse dynamics span 6 orders of magnitude in timescale:

- Dielectric relaxation: tau_d ~ epsilon/(q*mu_n*N_D) ~ 7 ns (bulk epi at N_D=8.5e13)
- Carrier transit: t_tr ~ W/(mu\*E) ~ 35 ps (electrons at 30 kV/cm)
- SRH recombination: tau_p = 600 ns
- Beam pulse duration: 10-200 ms

devsim supports BDF1, BDF2, and TRBDF transient methods but has NO built-in adaptive timestepping. The user must implement timestep control in the Python driver loop. Using a single fixed timestep either misses fast dynamics (too large) or makes ms-scale simulations impossibly slow (too small).

**Why it happens:** Users unfamiliar with stiff PDE systems set tdelta based on the "interesting" physics timescale and miss the fast modes. BDF1 is L-stable (will not oscillate with large dt) but only first-order accurate, introducing numerical diffusion that smears carrier dynamics. BDF2 is second-order but can oscillate with timesteps larger than ~5x the fastest mode.

**Prevention:**

1. Implement adaptive timestepping in the Python driver: start with tdelta ~ 0.1 ns, increase by factor 1.5 when convergence is easy (< 10 Newton iterations), decrease by factor 2 when convergence is hard (> 30 iterations).
2. Use BDF2 for accuracy, with TRBDF for the first step from DC initial condition.
3. For the FLASH problem: consider whether true transient simulation is actually needed. The v1.0 quasi-static approach (DC solve at each dose rate) was successful. True transient adds value only for: (a) intra-pulse onset dynamics (first ~100 ns), (b) inter-pulse carrier decay, (c) plasma buildup over multiple pulses.
4. Monitor charge_error (devsim parameter) at each timestep. If it exceeds tolerance, the timestep is too large.

**Detection:** Run at tdelta and tdelta/2. If CCE or current changes by > 1%, the timestep is too large. Total charge conservation (generated - recombined - collected = stored) should balance within 0.1%.

---

### Pitfall 7: Incomplete Ionization Not Updated with Temperature

**What goes wrong:** `device.py` computes `N_A_ionized = ionized_acceptor_concentration(N_A, T)` at device creation and uses it as a fixed doping value. For donors in the n-epi, the code assumes FULL ionization (`N_D` is used directly, never passed through an ionization model). When T changes:

For nitrogen donors in the n-type epi (most critical for device behavior):

- E_D(hex) = 50 meV: ~85% ionized at 300K, ~82% at 313K (clinical range)
- E_D(cub) = 92 meV: ~55% ionized at 300K, ~52% at 313K

The bulk epi doping N_D_bulk = 8.5e13 is low enough that freeze-out effects are modest at clinical temperatures. But N_D_junction = 2.9e15 has higher ionization (closer to Mott transition) so the graded profile shape subtly changes with T.

For aluminum acceptors in p+ substrate:

- E_A = 220 meV: ~10% ionized at 300K, ~11% at 313K
- N_A = 1e19 means even 10% gives 1e18 -- still degenerately doped. Change is marginal.

**Why it happens:** Incomplete ionization is a self-consistent problem (ionized fraction depends on Fermi level). The v1.0 code solves it once at device creation. For T-sweeps, the device must be recreated at each T anyway (because devsim parameters are set at creation), so incomplete ionization naturally gets the new T.

**Prevention:**

1. For clinical range (303-313K), the change in N_D ionization is < 3%. Pre-computing at each T and treating as fixed is acceptable.
2. Add a `ionized_donor_concentration(N_D, T)` function (currently only acceptors are modeled) for completeness, even if the effect is small.
3. Do NOT implement self-consistent incomplete ionization inside the devsim solve loop for this milestone -- it adds complexity for < 3% effect.

**Detection:** Compute ionized donor fraction at N_D = 8.5e13 and N_D = 2.9e15 for T=300K and T=313K. If the fractional change is < 2%, safe to treat as parametric (recompute at each T, don't iterate).

---

### Pitfall 8: Devsim Contact Equations for Surface Physics in 1D

**What goes wrong:** The existing code uses `simple_physics.CreateSiliconDriftDiffusionAtContact()` for Ohmic contact boundary conditions. Adding surface recombination requires modifying contact equations to include a surface recombination current term. In devsim, contact equations use a different API (`contact_node_model`, `contact_equation`) than region equations (`node_model`, `equation`).

If you accidentally implement surface recombination as a region node_model instead of a contact_node_model, it gets applied at EVERY mesh node, producing completely wrong results -- effectively adding recombination everywhere in the bulk.

**Why it happens:** devsim documentation for custom contact equations is sparse. The `simple_physics` module handles Ohmic contacts but does not include surface recombination. The distinction between contact and region models is not obvious from the API.

**Consequences:** If applied as a region model, SRH+surface recombination rate increases everywhere, drastically reducing carrier concentrations and producing unphysically low CCE. The error may look like "the surface recombination velocity is too high" when the real problem is that it is being applied in the wrong place.

**Prevention:**

1. In a 1D simulation, "surface" means the contact nodes only (x=0 for anode, x=L for cathode). Surface recombination must be added via `devsim.contact_equation()` as an additional `contact_node_model`.
2. Implement and test at one contact at a time. After adding surface recombination at the anode, verify that carrier concentrations change ONLY near x=0, not in the bulk.
3. Study the devsim diode examples for how contact equations are structured before writing custom ones.
4. Given Pitfall 4 (surface recombination alone cannot explain 18 pA), this is lower priority than implementing field-enhanced generation.

**Detection:** After adding any surface term, plot n(x) and p(x) across the entire device. The effect should be localized to within a few diffusion lengths of the contact. If carrier concentrations change uniformly across the device, the model is applied in the wrong scope.

---

### Pitfall 9: SRH Lifetime Temperature Dependence Model for SiC

**What goes wrong:** The current code uses fixed SRH lifetimes: tau_n = 1e-9 s, tau_p = 6e-7 s. These are 300K values. The temperature dependence of SRH lifetime in SiC is NOT a simple power law like Si. The dominant recombination center in 4H-SiC (Z1/2 center, carbon vacancy) has:

- Capture cross-section sigma_n(T) that can be either thermally activated or weakly T-dependent depending on the charge state
- The SRH lifetime tau = 1/(sigma _ v_th _ N_t) where v_th ~ T^(1/2) (thermal velocity)

If the capture cross-section is T-independent, tau ~ T^(-1/2) (lifetime DECREASES with T, opposite to some Si models). If the capture cross-section is thermally activated (sigma ~ exp(-E_b/kT)), lifetime can increase or decrease depending on the barrier E_b.

**Why it happens:** Users apply the Si approximation tau_SRH ~ T^(+1/2) or tau_SRH = constant. For the Z1/2 center in SiC, the experimental data (IEEE Access 2023) shows lifetime increases with T in some samples but decreases in others, depending on the dominant defect.

**Consequences:** Wrong T-dependence of recombination rate affects transient carrier decay, steady-state CCE at different temperatures, and the predicted temperature coefficient of dark current.

**Prevention:**

1. For the clinical range (303-313K, only 13K span), the SRH lifetime change is < 5% regardless of which model is used. Keeping tau constant is acceptable for v1.1.
2. If extending beyond clinical range later, parameterize tau(T) using the thermal velocity scaling: tau(T) = tau_300 \* sqrt(300/T) as the baseline model. This assumes T-independent capture cross-section, which is the simplest defensible assumption.
3. Document the assumption clearly: "SRH lifetimes assumed T-independent for the 303-313K range studied."

**Detection:** Run CCE at T=300K and T=313K with constant lifetime and with tau(T) ~ T^(-1/2). If results differ by < 0.5%, the T-dependence of lifetime is negligible for this application.

---

### Pitfall 10: Mobility Set Once at Device Creation, Not Updated

**What goes wrong:** `device.py` computes mobility via `mobility_caughey_thomas(N_D)` at device creation and sets it as a scalar region parameter (`mu_n`, `mu_p`). This means:

1. Mobility is position-independent (uses a single N_D, not the graded profile)
2. Mobility is not updated if T changes after device creation
3. For transient simulations, mobility cannot vary with local conditions

For the graded doping profile (N_D from 2.9e15 at junction to 8.5e13 in bulk), the mobility variation is modest at 300K: mu_n(2.9e15) ~ 935 vs mu_n(8.5e13) ~ 949, a ~1.5% difference. But at T!=300K, the doping-dependent part and T-dependent part compound.

**Why it happens:** devsim allows both scalar parameters and node models for mobility. The current code uses a scalar parameter because it was simpler and the position variation was small. The `simple_dd.CreateElectronCurrent(device, region, "mu_n")` function references `mu_n` by name and works with either a parameter or a node model.

**Prevention:**

1. For T-dependent parametric sweeps: recreate the device at each T with updated mobility. This is the simplest approach and requires no structural code changes.
2. For position-dependent mobility: create `mu_n` as a node model using the local Donors value. This is a one-line change:
   ```python
   CreateNodeModel(device, region, "mu_n",
       f"{mu_min_T} + ({mu_max_T} - {mu_min_T}) / (1 + (Donors / {N_ref})^{alpha})")
   ```
3. Do NOT add high-field mobility saturation (v_sat model) unless E-fields exceed 100 kV/cm. At -60V across 10 um, E_max ~ 60 kV/cm -- below saturation for SiC.

**Detection:** At 300K, the error from using a single mobility value is < 2%. Monitor this: if the epi doping gradient steepens or bias increases, position-dependent mobility becomes more important.

---

## Minor Pitfalls

### Pitfall 11: SRH n1/p1 Parameters Not Updated with Temperature

**What goes wrong:** The existing code sets `n1 = n_i` and `p1 = n_i` in `device.py` (midgap trap assumption). When T changes, n1 and p1 should update to n_i(T). Currently these are set once at device creation using `params.n_i_300`. If n_i is updated to a T-dependent value in the region parameter but n1/p1 are not, the SRH rate formula becomes inconsistent.

**Prevention:** When updating n_i to n_i(T), also update n1 and p1 to match. For midgap traps, n1 = p1 = n_i is correct at any T. If switching to non-midgap traps (e.g., Z1/2 at E_C - 0.65 eV), compute n1 = n_i _ exp((E_trap - E_i)/kT) and p1 = n_i _ exp((E_i - E_trap)/kT).

---

### Pitfall 12: Auger Coefficients Temperature Dependence

**What goes wrong:** Current Auger coefficients (C_n=5e-31, C_p=2e-31 cm^6/s) are 300K values. Auger coefficients in SiC are poorly characterized as a function of T.

**Prevention:** Keep Auger coefficients constant for v1.1. The v1.0 finding was that Auger is negligible at FLASH dose rates (CCE is flat across 20-230 Gy/s), so even a 2x change in Auger coefficient has no practical impact. Document as assumption.

---

### Pitfall 13: Transient Initial Condition Discontinuity

**What goes wrong:** When switching from DC to transient mode, the initial condition must be a self-consistent DC solution. If you start transient simulation from a non-converged state, the first BDF step sees a large residual and diverges. devsim provides `set_initial_condition()` but the documentation is minimal.

**Prevention:**

1. Always run a full DC solve to convergence before starting any transient simulation.
2. Use `type="transient_dc"` for the first solve to establish the time-derivative state vectors.
3. Start with small tdelta (0.1 \* dielectric_relaxation_time ~ 0.7 ns) for the first few steps, then increase.

---

### Pitfall 14: Charge Conservation Drift in Long Transient Simulations

**What goes wrong:** Over thousands of transient timesteps, numerical errors accumulate. Without monitoring, the total charge drifts, and CCE computed from integrated current becomes unreliable. This is especially dangerous for FLASH simulations where you care about the ratio of collected to generated charge.

**Prevention:**

1. After each transient step, compute total electron and hole charge via integration over the mesh.
2. Track: charge_generated(t) - charge_recombined(t) - charge_collected(t) should equal delta_stored_charge(t).
3. If imbalance exceeds 1%, reduce timestep or tighten convergence (reduce `charge_error` parameter in devsim solve).

---

### Pitfall 15: compute_ni(300) Disagrees with Literature n_i_300

**What goes wrong:** The `compute_ni(T)` function in `sic_material.py` uses specific values for DOS effective masses (m_e=0.77*m0, m_h=1.0*m0) and M_c=3 to compute NC, NV, and then n_i. At 300K, this may return a value different from the literature value `n_i_300 = 5e-9` stored as a constant. Different sources report n_i(300K) for 4H-SiC ranging from 5e-9 to 1.6e-8 cm^-3 depending on which DOS effective masses are used.

If `compute_ni(300)` returns 6.5e-9 while the validated code uses 5e-9, switching to the computed value changes ALL equilibrium carrier concentrations, built-in potential, and dark current by ~30%.

**Prevention:**

1. Before switching to `compute_ni(T)`, verify its 300K output against the literature value.
2. If they disagree, adjust the DOS masses to reproduce the accepted n_i(300K) = 5e-9. This is standard practice in TCAD -- the DOS masses are effective parameters.
3. Alternatively, use the Varshni-based E_g(T) with the constraint that n_i(300K) = 5e-9 to determine an effective NC\*NV product, then use the T-scaling (NC ~ T^1.5, NV ~ T^1.5) for other temperatures.

---

## Phase-Specific Warnings

| Phase Topic                      | Likely Pitfall                                                                                        | Mitigation                                                                                  |
| -------------------------------- | ----------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| Temperature-dependent parameters | Pitfall 1 (scattered 300K constants), Pitfall 2 (Si/6H-SiC exponents), Pitfall 15 (n_i mismatch)      | Full codebase audit, use 4H-SiC Ayalew values, reconcile compute_ni with n_i_300            |
| Regression testing               | Pitfall 3 (breaking 300K baseline)                                                                    | Extract golden values FIRST, build regression test, make new physics toggleable             |
| Dark current / surface physics   | Pitfall 4 (surface recombination too weak), Pitfall 5 (TAT complexity), Pitfall 8 (contact equations) | Do the math first -- surface SRH cannot explain 18 pA; use simple field-enhanced generation |
| Transient FLASH dynamics         | Pitfall 6 (timestep), Pitfall 13 (initial conditions), Pitfall 14 (charge conservation)               | Adaptive timestepping in Python, DC init, conservation monitoring                           |
| Incomplete ionization            | Pitfall 7 (T coupling)                                                                                | Pre-compute at target T, skip self-consistent for clinical range                            |
| Mobility update                  | Pitfall 10 (fixed mobility)                                                                           | Recreate device at each T for parametric sweeps                                             |
| SRH model                        | Pitfall 9 (lifetime T-dependence), Pitfall 11 (n1/p1)                                                 | Keep constant for clinical range, update n1=p1=n_i(T)                                       |

## Sources

- [TU Wien Ayalew thesis - Low-Field Carrier Mobility](https://www.iue.tuwien.ac.at/phd/ayalew/node65.html) -- 4H-SiC Caughey-Thomas parameters: gamma_n=-2.40, gamma_p=-2.15, beta=-0.5. HIGH confidence.
- [TU Wien Ayalew thesis - Incomplete Ionization](https://www.iue.tuwien.ac.at/phd/ayalew/node75.html) -- ionization model for N donors and Al acceptors. HIGH confidence.
- [TCAD models of T and doping dependence in 4H-SiC (Potbhare et al.)](https://www.researchgate.net/publication/259510798) -- validated TCAD T-dependent mobility and bandgap models. MEDIUM confidence.
- [TCAD modeling of radiation-induced defects in 4H-SiC (2024)](https://arxiv.org/html/2407.11776v1) -- TAT modeling pitfalls, Z1/2 center, wrong trap selection effects. MEDIUM confidence.
- [Surface recombination velocities for 4H-SiC (2023)](https://www.sciencedirect.com/science/article/pii/S136980012300673X) -- SRV=150-700 cm/s for SiO2 passivated Si-face. MEDIUM confidence.
- [SiC detectors review (Frontiers in Physics, 2022)](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2022.898833/full) -- Z1/2 center as dominant recombination center. HIGH confidence.
- [DEVSIM solver documentation](https://devsim.net/solver.html) -- BDF1/BDF2/TRBDF transient methods. HIGH confidence.
- [DEVSIM command reference](https://devsim.net/CommandReference.html) -- transient solve API: tdelta, charge_error, type options. HIGH confidence.
- [Carrier lifetime T-dependence in 4H-SiC (IEEE Access)](https://ieeeaccess.ieee.org/featured-articles/lifetimedependence/) -- experimental tau(T) showing sample-dependent behavior. MEDIUM confidence.
- [Modified TCAD for 4H-SiC JBS - electron trapping effects (2025)](https://www.sciencedirect.com/science/article/abs/pii/S0026271425002896) -- TAT with specific traps can overestimate leakage. MEDIUM confidence.
- [Ioffe NSM 4H-SiC archive](https://www.ioffe.ru/SVA/NSM/Semicond/SiC/) -- reference n_i(T), bandgap, mobility values. HIGH confidence.

---

_Pitfalls research for: v1.1 milestone -- adding temperature dependence, surface physics, and transient dynamics to validated 1D 4H-SiC TCAD simulator_
_Researched: 2026-03-23_
