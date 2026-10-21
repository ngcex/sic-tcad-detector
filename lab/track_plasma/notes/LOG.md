# lab/track_plasma running log (issue #9)

These notes cover the axisymmetric single-ion track engine `etna/core/track/` (2026-10-21). All commands run from the repository root:

`uv run python lab/track_plasma/bench/track_bench.py SECTION [--mesh fast|default] [--devsim]`

The outputs are in `notes/SECTION.out`. Wall times were measured on an Apple M4 (10 cores) under a load average of 40-60 from other jobs, so they are pessimistic.

Meshes:

- **FAST:** 21 x 110 nodes, 6930 unknowns. r runs from 1 nm to 20 um (x1.6); z has h_min 5 nm and h_max 0.5 um.
- **DEFAULT:** 31 x 147 nodes, 13671 unknowns.

Unless stated otherwise:

- the stack is the default stack;
- 50 V reverse bias;
- the track is in the n- epi only (the `build_cce_let_table` convention);
- Auger is the constant Galeckas value;
- the time tolerance is rtol 1e-4.

## 1. Parity with the batched 1D engine (`parity`)

Same z mesh, single column (Auger off).

| V   | max \|psi - psi_batched\| | n rel   | steady CCE (batched J_cathode) | track (pair balance) | dCCE    |
| --- | ------------------------- | ------- | ------------------------------ | -------------------- | ------- |
| 0   | 6.7e-16 V                 | 1.3e-14 | 0.872199                       | 0.872199             | 3.3e-11 |
| 10  | 1.8e-15                   | 8.7e-14 | 0.999463                       | 0.999463             | 9e-16   |
| 50  | 7.1e-15                   | 3.3e-13 | 0.999879                       | 0.999879             | 1.8e-13 |
| 100 | 1.4e-14                   | 5.6e-13 | 0.999938                       | 0.999938             | 2.2e-14 |

The direct sum of the cathode contact fluxes is about 3e-5 off from this value. That sum is a small difference of 1e19-cm^-3 majority fluxes. The pair-balance estimator is therefore what `TrackResult.cce` reports:

1 - (recombined + electrons out at the anode + holes out at the cathode + stored) / N_gen

**Low injection.** The axisymmetric track equals the 1D column to within 1e-5 (`tests/test_track_engine.py`). It does NOT equal the steady CCE of the same depth profile:

- with n_i ~ 5e-9 cm^-3, SRH of the excess in the depletion region needs both carriers at the same point at the same time;
- a transient separates them, so it recombines less: 5.9e-5 vs 2.0e-4 of the pairs at 30 V on MESH_TINY.

This means a steady-illumination CCE (`run_cce`) is not the single-particle CCE even at low LET.

## 2. CCE vs LET, Gaussian w = 50 nm (the commercial-TCAD style track; `letscan`)

Recombined fraction:

| LET keV/um | 0.1     | 1       | 10      | 100     | 300    | 1000    | 3000    |
| ---------- | ------- | ------- | ------- | ------- | ------ | ------- | ------- |
| 50 V       | 2.7e-5  | 2.8e-5  | 2.7e-5  | 2.5e-5  | 4.0e-5 | 1.32e-4 | 6.9e-4  |
| 10 V       | 1.32e-4 | 1.38e-4 | 1.46e-4 | 2.62e-4 | 5.1e-4 | 1.15e-3 | 2.67e-3 |
| 0 V, CCE   | 0.8722  | 0.8736  | 0.8805  | 0.8990  | 0.9110 | 0.9237  | 0.9330  |

Peak density at t = 0 is 1.6e15 x (LET / 0.1). At 0 V the CCE RISES with LET. The undepleted epi's diffusion loss (12.8 %) shrinks as the dense track self-screens and the plasma pushes charge out ("funnelling"-like). The table's linear-response assumption therefore fails at 0 V by 7 % between LET 0.1 and 3000.

Wall time per event (FAST): 4-15 s under load, 2-5 s unloaded (`timing`).

## 3. The core decides (`core`, 50 V)

| radial model              | r50     | LET 300: recombined | LET 1000: recombined | peak density (1000) |
| ------------------------- | ------- | ------------------- | -------------------- | ------------------- |
| Gaussian 50 nm            | 42 nm   | 4.0e-5              | 1.3e-4               | 1.6e19              |
| Gaussian 10 nm            | 8.3 nm  | 1.7e-4              | 1.6e-3               | 4.1e20              |
| Gaussian 2 nm             | 1.7 nm  | 3.4e-3              | 3.5e-2               | 9.9e21              |
| WHK 1 MeV/u (no spread)   | 0.50 nm | 3.0e-2              | 0.184                | 8.2e22 (first cell) |
| WHK thermalized s = 1 nm  | 1.2 nm  | 4.7e-3              | 4.7e-2               | 2.2e22              |
| WHK thermalized s = 10 nm | 9.0 nm  | 1.5e-4              | 1.3e-3               | 3.6e20              |
| WHK thermalized s = 50 nm | 42 nm   | 3.9e-5              | 1.3e-4               | 1.6e19              |

Between a sub-nm WHK core and a 50 nm core, the recombination loss spans three decades. The bare WHK core puts mean densities of order 1e22-1e23 cm^-3 in the first cell. At those densities drift-diffusion (and a constant Auger C) is not valid. The thermalization radius s is unknown for 4H-SiC; it is a parameter of the result, NOT a fit.

## 4. Opt-in physics (`physics`, WHK thermalized s = 10 nm)

Recombined fraction:

| model                                                        | LET 300 | LET 1000 |
| ------------------------------------------------------------ | ------- | -------- |
| Auger constant (Galeckas 1997)                               | 1.46e-4 | 1.31e-3  |
| Auger injection-dependent (Tanaka 2023 via arXiv:2410.06798) | 4.4e-5  | 1.24e-4  |
| Auger off (SRH only)                                         | 3.5e-5  | 7.5e-5   |
| constant + carrier-carrier scattering (Lades fit)            | 1.49e-3 | 3.17e-2  |
| injection + carrier-carrier scattering                       | 1.31e-4 | 8.66e-4  |

- **Injection-dependent Auger** (the measured fall of C above ~5e18) cuts the Auger part of the loss by about 10x.
- **Carrier-carrier scattering** lowers the ambipolar diffusivity in the plasma. The plasma lives longer, so Auger acts longer: loss x10-25.
- Both are opt-in. The defaults stay constant Auger and no carrier-carrier scattering.

## 5. Ions at normal incidence (`ions`)

Setup: NIST Bragg profile, WHK at the local speed, thermalized by s; `region="all"` (pairs in the p+ and the n+ too).

| ion                        | V   | s     | CCE    | recombined | wrong contact |
| -------------------------- | --- | ----- | ------ | ---------- | ------------- |
| p 1 MeV (stops in the epi) | 50  | 10 nm | 0.9811 | 1.28e-3    | 1.76e-2       |
| p 3 MeV                    | 50  | 10 nm | 0.9425 | 2.2e-3     | 5.5e-2        |
| alpha 5.486 MeV            | 50  | 10 nm | 0.9405 | 2.5e-3     | 5.7e-2        |
| alpha 8 MeV                | 50  | 10 nm | 0.9451 | 2.2e-3     | 5.3e-2        |
| alpha 5.486 MeV            | 10  | 10 nm | 0.9379 | 3.1e-3     | 5.9e-2        |

The "wrong contact" column is the dead p+ (0.3 um): electrons generated there exit at the anode. For tracks that cross the epi, it also includes the n+, whose holes exit at the cathode. At these LETs recombination is 0.1-0.3 % and nearly bias independent from 10 to 50 V.

## 6. Fig. 4b proxy (`fig4b`, microdosimeter preset, N_D 1e14)

**This is a proxy, not the measurement geometry.** The real track is tilted 65 deg and has no axial symmetry. The proxy is a vertical column with:

- the tilted track's per-PATH density (the true local plasma density);
- the projected depth (7.6 um);
- WHK at the local speed, thermalized by s.

Results (shape = CCE / CCE(-50 V)):

| s     | 0 V     | -5      | -10     | -20     | -30     | -40     | -50 | rms vs data |
| ----- | ------- | ------- | ------- | ------- | ------- | ------- | --- | ----------- |
| 1 nm  | 0.98992 | 0.99913 | 0.99967 | 0.99990 | 0.99996 | 0.99999 | 1   | 0.034       |
| 10 nm | 0.98987 | 0.99912 | 0.99967 | 0.99990 | 0.99996 | 0.99999 | 1   | 0.034       |
| 50 nm | 0.98987 | 0.99912 | 0.99967 | 0.99990 | 0.99996 | 0.99999 | 1   | 0.034       |
| data  | 0.9083  | 0.9662  | 0.9806  | 0.9904  | 0.9932  | 0.9976  | 1   |             |

- The plasma loss is 0.7 % at s = 1 nm, 0.02 % at 10 nm and 0.004 % at 50 nm. In each case it is nearly FLAT in bias above 5 V, so it cancels in the normalised shape.
- The drift-diffusion + SRH + Auger track does not reproduce the data's 9 % deficit at 0 V or its 1/V-like tail, for any core radius tried.
- This confirms lab/microdosimeter LOG §9 with track physics included. The field dependence the data need is not produced by this recombination physics at these LETs.
- **Open candidates:**
  - velocity saturation and field-dependent erosion of a TILTED column;
  - trapping;
  - carrier-carrier scattering at this geometry (not run here);
  - the 3 us shaping (ballistic deficit).
- Nothing was tuned.

## 7. Numerics

**Time tolerance** (`rtol`, LET 1000, s = 10 nm). CCE converges as rtol^~1:

| rtol         | 1e-2    | 1e-3    | 1e-4    | 1e-5   | 1e-6 |
| ------------ | ------- | ------- | ------- | ------ | ---- |
| dCCE vs 1e-6 | -1.4e-5 | -2.6e-6 | -1.9e-7 | 8.7e-8 | 0    |
| steps        | 57      | 101     | 199     | 412    | 874  |

The default rtol is 1e-4.

**Mesh** (`mesh`, LET 1000, s = 10 nm):

| mesh    | M x K    | PHD       | wall   |
| ------- | -------- | --------- | ------ |
| fast    | 21 x 110 | 1.3133e-3 | 6.6 s  |
| default | 31 x 147 | 1.3708e-3 | 15.7 s |
| fine    | 47 x 197 | 1.4113e-3 | 40 s   |

The PHD is converging. The remaining mesh error is about 5-8 % of the PHD itself, which is 1e-4 in CCE. Use DEFAULT or finer for PHD numbers.

**Ledgers:**

- The pair ledger closes to 1e-9-1e-10 at LET >= 100. At LET 0.1-1 it sits at the float64 floor of the 1e19 majority carriers, ~1e-6 (`ledger_floor` is reported).
- Before the potential was solved relative to the dark state (`kernels` docstring), the 0 V, LET 1 ledger error was 13x larger (6.8e-5 vs 5.1e-6) and the step count doubled (412 vs 214: round-off noise triggered LTE rejections).
- Before negative excursions were allowed in the transient Newton (the positivity clamp broke conservation), the ledger error was 5e-5 at LET 1000.

## 8. Electrothermal (`thermal`, opt-in, coupled mu(T), tau(T))

| case                    | dT_max                    | dCCE (coupled - isothermal) | thermal ledger                  | Poynting residual |
| ----------------------- | ------------------------- | --------------------------- | ------------------------------- | ----------------- |
| LET 1000, WHK s = 10 nm | 118 K (core, ps)          | -2.5e-4                     | 3e-8 (adiabatic) / 2e-13 (sink) | 2.3e-5            |
| LET 1000, bare WHK      | 3.0e4 K (unphysical core) | -0.29                       | 7e-8                            | 4.0e-5            |
| LET 300, Gaussian 50 nm | 1.8 K                     | -1.9e-8                     | 8e-8                            | 3.1e-5            |

- **Heat budget** (LET 1000): Joule 1.09e-11 J >> prompt (E_pair - E_g) 9.3e-13 J >> recombination 1e-15 J. Joule per collected charge = V + ~3 V, the built-in-like drop (tested).
- **When heating matters:** for a single ion it is negligible unless the core is sub-10 nm. There the drift-diffusion model is outside its validity anyway. The 3e4 K bare-WHK result is a symptom of that, not a prediction.

## 9. Speed

- **This engine:**
  - Per event: 2-5 s on FAST and 4-13 s on DEFAULT (`timing`, loaded machine).
  - The cost is dominated by SuperLU factorizations (13 ms FAST, 32 ms DEFAULT). Chord Newton reuses each factorization across stages and steps: 27-150 per event instead of ~900.
  - Chord reuse across stages and steps cut the LU count of a LET-1000 event from 424 to 120 (same CCE to 1e-9).
- **GPU:** not attempted. About 1e4 unknowns is a sparse-direct-solve problem; the batched-1D GPU result in lab/README applies.
- **devsim comparison** (`timing --devsim`, `timing_devsim.out`): `build_cce_let_table(hw 50, 50 V, LET 1 / 31.6 / 1000, method="full")` on the 2D Cartesian sheet took 140 s for 3 LETs, including the ~25 s device build.
  - CCE: 0.999967 / 0.999970 / 0.999970.
  - Track engine (FAST, Gaussian cylinder of the same 1 um sigma): 1.1 / 1.2 / 1.9 s per LET. CCE 0.9999686 / 0.9999684 / 0.9999490.
  - At low LET the two agree to 1.6e-6. That is the expected 1D-limit agreement, from different meshes and time steppers.
  - At LET 1000 the cylinder (a denser track than the sheet) already loses 2e-5 more.
  - About 40x faster per LET point, but these are not the same physics: the sheet stays linear by construction.
