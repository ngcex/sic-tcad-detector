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

FAST mesh: PHD numbers carry ~5-8 % mesh error (section 7).

Recombined fraction:

| LET keV/um | 0.1     | 1       | 10      | 100     | 300    | 1000    | 3000    |
| ---------- | ------- | ------- | ------- | ------- | ------ | ------- | ------- |
| 50 V       | 2.7e-5  | 2.8e-5  | 2.7e-5  | 2.5e-5  | 4.0e-5 | 1.32e-4 | 6.9e-4  |
| 10 V       | 1.32e-4 | 1.38e-4 | 1.46e-4 | 2.62e-4 | 5.1e-4 | 1.15e-3 | 2.67e-3 |
| 0 V, CCE   | 0.8722  | 0.8736  | 0.8805  | 0.8990  | 0.9110 | 0.9237  | 0.9330  |

Peak density at t = 0 is 1.6e15 x (LET / 0.1). At 0 V the CCE RISES with LET. The undepleted epi's diffusion loss (12.8 %) shrinks as the dense track self-screens and the plasma pushes charge out ("funnelling"-like). The table's linear-response assumption therefore fails at 0 V by 7 % between LET 0.1 and 3000.

Wall time per event (FAST): 4-15 s under load, 2-5 s unloaded (`timing`).

## 3. The core decides (`core`, 50 V)

FAST mesh: PHD numbers carry ~5-8 % mesh error (section 7).

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

The sub-nm cores are FIRST-CELL LIMITED. At LET 1000 the recombined fraction depends on the axis cell r_min (FAST, bench section `rmin`):

| core | r_min 1 nm | r_min 0.5 nm | r_min 0.25 nm |
|---|---|---|---|
| bare WHK | 0.184 | 0.207 | 0.219 |
| WHK + s = 1 nm | 0.0470 | 0.0504 | 0.0514 |

- **Bare WHK** does not converge. Its first-cell density reaches 6e23 cm^-3, above the atomic density of SiC, so its numbers are not physics.
- **s = 1 nm** converges to about 0.052.
- **s >= 10 nm** is converged in r_min (section 7 is for the mesh).

## 4. Opt-in physics (`physics`, WHK thermalized s = 10 nm)

FAST mesh: PHD numbers carry ~5-8 % mesh error (section 7).

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

FAST mesh: PHD numbers carry ~5-8 % mesh error (section 7).

Setup: NIST Bragg profile, WHK at the local speed, thermalized by s; `region="all"` (pairs in the p+ and the n+ too).

| ion                        | V   | s     | CCE    | recombined | wrong contact |
| -------------------------- | --- | ----- | ------ | ---------- | ------------- |
| p 1 MeV (stops in the epi) | 50  | 10 nm | 0.9811 | 1.28e-3    | 1.76e-2       |
| p 3 MeV                    | 50  | 10 nm | 0.9425 | 2.2e-3     | 5.5e-2        |
| alpha 5.486 MeV            | 50  | 10 nm | 0.9405 | 2.5e-3     | 5.7e-2        |
| alpha 8 MeV                | 50  | 10 nm | 0.9451 | 2.2e-3     | 5.3e-2        |
| alpha 5.486 MeV            | 10  | 10 nm | 0.9379 | 3.1e-3     | 5.9e-2        |

CCE is relative to the pairs generated INSIDE the 11.3 um domain: the 3 MeV proton and both alphas run past the n+ end (alpha 5.486 MeV range 18.5 um).

The "wrong contact" column is the dead p+ (0.3 um): electrons generated there exit at the anode. For tracks that cross the epi, it also includes the n+, whose holes exit at the cathode. At these LETs recombination is 0.1-0.3 % and nearly bias independent from 10 to 50 V.

## 6. Fig. 4b proxy (`fig4b`, microdosimeter preset, N_D 1e14)

FAST mesh: PHD numbers carry ~5-8 % mesh error (section 7).

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
| LET 1000, WHK s = 10 nm | 118 K (core, ps)          | -2.5e-4                     | 3e-8 (adiabatic) / 3e-8 (sink) | 2.3e-5            |
| LET 1000, bare WHK      | 3.0e4 K (unphysical core) | -0.29                       | 7e-8                            | 4.0e-5            |
| LET 300, Gaussian 50 nm | 1.8 K                     | -1.9e-8                     | 8e-8                            | 3.1e-5            |

- **Heat budget** (LET 1000): Joule 1.09e-11 J >> prompt (E_pair - E_g) 9.3e-13 J >> recombination 1e-15 J. Joule per collected charge = V + ~3 V, the built-in-like drop (tested).
- **When heating matters:**
  - 50 nm core: negligible (1.8 K, dCCE -2e-8).
  - 10 nm core at LET 1000: a 118 K ps-scale core rise and dCCE = -2.5e-4, about 19 % of the 1.31e-3 PHD. The rise is stable in time step: dCCE -2.517e-4 at rtol 1e-4 and -2.511e-4 at 1e-5.
  - Below that, drift-diffusion is outside its validity. The 3e4 K bare-WHK result is a symptom of that, not a prediction.
- **Sink ledger.** Since the fix that integrates the conductive flux into the fixed nodes with the step weights, the sink-mode ledger is an independent check: 3e-8. Before the fix it was sources minus stored, i.e. an identity, and the 2e-13 in `thermal.out` comes from that version.
- **Poynting residual vs rtol** (LET 300, Gaussian 20 nm, MESH_TINY):

  | rtol | 1e-3 | 1e-4 | 1e-5 | 1e-6 |
  |---|---|---|---|---|
  | residual | 1.4e-4 | 3.5e-5 | 8.0e-6 | 1.8e-6 |

  The residual falls with the time tolerance, so it is time-discretisation error, not a missing term.
- **0 V and the domain radius.** The CCE rise with LET at 0 V holds when R grows: at R = 20 / 50 / 100 um, CCE(0.1) = 0.87218 / 0.87218 / 0.87217 and CCE(1000) = 0.92372 / 0.92253 / 0.92253. With R = 20 um the high-LET value is 1e-3 high.

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
  - Excluding the devsim build: about 20-35x faster per LET point on a loaded machine. They are not the same physics: the sheet stays linear by construction.

# Follow-up (2026-10-21): tilted tracks, Fig. 4b calibration, float64 floor

Bench: `uv run python lab/track_plasma/bench/tilted_bench.py SECTION [--workers 8]` with SECTION one of `slices`, `scan`, `relin`, `calibrate`.

Outputs:

- `notes/tilted_slices.out`;
- `tilted_scan.out` / `tilted_scan.json`;
- `tilted_relin.out`;
- `tilted_calibrate.out` / `tilted_calibration.json`;
- the plasma tables, in `notes/tables/*.json` (`PlasmaTable.to_dict`).

Code:

- `etna/core/track/cross_section.py`: the transverse slice.
- `etna/core/track/tilted.py`: the 2.5D composition and `calibrate_tau_scale`.
- `kernels.assemble_excess` and `bern_diff`: the excess unknowns.
- `etna/api/track.py`: the facades.

Wall times are on the M4, under a load average of 10-30 from other jobs.

## 10. Was the proxy wrong? Transverse slices (`slices`)

A track tilted 65° sees the field mostly ACROSS it: E⊥ = 0.906 E, E∥ = 0.423 E. The axisymmetric engine used for the §6 proxy can only hold a track parallel to the field.

The slice solves the other limit: the plane perpendicular to a long track, with E⊥ imposed and the engine's physics. It is `cross_section.slice_loss`, a 2D Cartesian half plane with absorbing edges at R.

Checks (`tests/test_track_tilted.py`):

- The linear ψ with n = p = n_i is an exact discrete steady state.
- Sampled radial lumping conserves pairs to 1e-9 and matches the exact erf-product Gaussian.
- At E = 0 with SRH off, the slice reproduces the axisymmetric radial problem on a disc mesh to within 3 %.
- The slice is linear in LET below 0.01 keV/µm: within 0.5 % at E⊥ from 250 to 6.4e4 V/cm.
  - At LET 1 the column's self-field already raises f by 40-50 % at 1e4 V/cm.
  - The first scan's f_lin used LET 1; `relin` recomputed it at 0.01.
- Ledgers close to 1e-8.

Loss of one slice: LET 300 keV/µm, WHK 0.4 MeV/u thermalized by s = 10 nm. The two fields are E⊥ at the track end at −20 / −50 V (N_D 1e14).

| variant | E⊥ 1.45e4 V/cm | E⊥ 4.2e4 V/cm | difference |
|---|---|---|---|
| R 5 / 10 / 20 µm | 6.2e-4 / 8.1e-4 / 9.7e-4 | 2.27e-4 / 2.45e-4 / 2.55e-4 | 4-7e-4 |
| in-plane field gradient (N_D sin²θ), R 10 / 20 | 9.1e-4 / 1.28e-3 | 2.48e-4 / 2.62e-4 | 7e-4 / 1.0e-3 |
| mesh r_min 0.25 nm, x1.25 | 8.36e-4 | 2.51e-4 | 5.9e-4 |
| WHK s 1 nm / 50 nm | 5.6e-3 / 6.8e-4 | 5.0e-3 / 1.1e-4 | 6e-4 / 5.6e-4 |
| Gaussian 300 nm | 6.3e-4 | 8.7e-5 | 5.4e-4 |
| Auger off / injection-dependent | 6.7e-4 / 6.8e-4 | 1.08e-4 / 1.19e-4 | 5.6e-4 |
| carrier-carrier scattering | 2.67e-3 | 2.10e-3 | 5.7e-4 |
| frozen-field v_sat (Lades / Nilsson, ∥c) | 1.05e-3 | 4.7e-4 | 5.8e-4 |
| **τ_n, τ_p x 1/10** | 6.8e-3 | 1.21e-3 | 5.6e-3 |
| **τ_n, τ_p x 1/30** | 1.98e-2 | 3.36e-3 | 1.6e-2 |

**Verdict.** The proxy's conclusion holds. With the library lifetimes (τ_n 1 ns, τ_p 600 ns), the transverse geometry changes the plasma loss by at most 1e-3 between the fields of −20 and −50 V. The data rise by 1e-2 there.

What moves the loss:

- **Lifetime.** The bias-dependent part is the SRH loss during the plasma time. It scales as about 1/E⊥ and as 1/(τ_n + τ_p), the high-injection SRH lifetime.
- **Core size and Auger.** These move the loss nearly flat in bias (s = 1 nm adds 5e-3 at both fields). The free overall scale absorbs that offset.
- **Velocity saturation.** It raises the high-field loss and so FLATTENS the field dependence, which is the wrong direction.
- **Box size (main systematic).** At 1.45e4 V/cm the plasma grows to r ~ N_L q / (2π ε E⊥) ~ 8 µm before erosion. The R and field-gradient dependence (±25-60 %) is the main systematic of the local model.

## 11. 2.5D composition vs Fig. 4b (`scan`)

The composition is `tilted.tilted_alpha_cce`:

CCE(V) = CCE_column(V) − Σ_k w_k f_ex(s_k, E⊥(z_k; V))

- **CCE_column** is the single-particle linear response of the tilted NIST-Bragg generation on the 1D column (track engine, M = 1, a transient at a peak excess of 1e10).
- **f_ex** = f(local LET, local WHK speed, s) − f(LET 0.01).
  - It is tabulated on 6 path nodes x 11 fields (77 slices, 70-270 s with 8 workers) and interpolated.
  - E⊥ = |E| sin θ, from the column's dark state.
- **Device:** the microdosimeter preset (N_D 1e14, p+ 0.3 µm), with epi-only pairs.

Shape = CCE / CCE(−50 V). χ² is computed with a free overall scale (the §9 convention):

| τ scale (τ_p) | 0 V | −5 | −10 | −20 | −30 | −40 | χ² 7 pts (σ 1e-3) | χ² −20..−50 |
|---|---|---|---|---|---|---|---|---|
| 1 (600 ns) | 0.980 | 0.9974 | 0.9989 | 0.9996 | 0.9999 | 1.0000 | 3844 | 51.9 |
| 1/3 | 0.937 | 0.9921 | 0.9968 | 0.9989 | 0.9996 | 0.9999 | 814 | 44.7 |
| 1/10 | 0.819 | 0.9747 | 0.9893 | 0.9964 | 0.9986 | 0.9995 | 7960 | 24.5 |
| 1/20 | 0.710 | 0.9516 | 0.9791 | 0.9928 | 0.9972 | 0.9991 | 36269 | 8.7 |
| 1/30 | 0.664 | 0.9304 | 0.9692 | 0.9893 | 0.9959 | 0.9986 | 54307 | 7.7 |
| 1/50 | 0.580 | 0.8926 | 0.9505 | 0.9824 | 0.9932 | 0.9977 | 99886 | 48.2 |
| data | 0.9083 | 0.9662 | 0.9806 | 0.9904 | 0.9932 | 0.9976 | | |

At τ x1 the 2.5D shape has rms 0.033 from the data, the same as the normal-incidence proxy (0.034). The 50 V plasma loss is 3.1e-4, against 2e-4 for the proxy at s = 10 nm.

## 12. Calibration (`calibrate`): what the data constrain

`calibrate_tau_scale`: flat prior in log τ_scale, PCHIP interpolation of the shapes between the scanned τ, free scale. **This is a CALIBRATION, not a prediction.**

**(a) Base model** (s 10 nm, R 10 µm, N_D 1e14, constant Auger):

| points | σ | τ_scale mode (68 % / 95 %) | τ_p | χ² / dof |
|---|---|---|---|---|
| −20..−50 V | 1e-3 | **0.039** (0.035-0.047 / 0.031-0.057) | 23 ns | 6.4 / 2 |
| −20..−50 V | 5e-4 (digitisation) | 0.039 (0.037-0.042) | 23 ns | 25.4 / 2 |
| −20..−50 V | 2e-3 | 0.039 (0.032-0.066 / 0.026-0.15) | 23 ns | 1.6 / 2 |
| all 7 | 1e-3 | 0.24 | 144 ns | 371 / 5 |

- **Tail:** reproduced only at σ ≥ 1e-3 (the §9 convention), with χ² 6.4 for 2 dof. At the stated read-off σ of 5e-4 it is not reproduced (χ² 25).
- **All seven points:** no single τ fits them (χ² 371 for 5 dof). At the tail mode the model gives 0.681 / 0.938 / 0.972 at 0 / −5 / −10 V, against data 0.908 / 0.966 / 0.981. That is far too much loss at low bias.

**(b) Systematic band of the tail-best τ_scale.** Each variant's shift from the base is interpolated in log τ between its 1/10 and 1/30 runs:

| variant | tail-best τ_scale |
|---|---|
| s 1 / 50 nm | 0.0385 / 0.0389 |
| R 5 / 20 µm | 0.029 / 0.051 |
| field gradient (R 20) | 0.070 |
| v_sat bracket | 0.046 |
| injection-dependent Auger | 0.039 |
| carrier-carrier scattering | 0.039 |
| T 295 K | 0.039 |
| N_D 0.5 / 1.5 / 2e14 | 0.040 / 0.039 / 0.041 |

Band: **τ_scale 0.03-0.07, i.e. τ_p ≈ 17-42 ns.** This is also the band of τ_n + τ_p, the high-injection SRH lifetime. The box and gradient systematics dominate, and they are larger than the statistical 68 % interval.

**(c) Where the low-bias excess comes from.** Plasma loss split by the local E⊥:

| τ | V | < 2e3 V/cm | 2e3-1e4 | > 1e4 | column / column(50) |
|---|---|---|---|---|---|
| 1/30 | 0 | 0.173 | 0.045 | 1e-4 | 0.881 |
| 1/30 | −5 | 0 | 0.062 | 0.008 | 0.997 |
| 1/30 | −10 | 0 | 0.014 | 0.018 | 0.998 |
| 1/30 | −20..−50 | 0 | 0 | 0.0135 → 0.0036 | |
| 1/20 | 0 | 0.169 | 0.032 | 8e-5 | **0.9085** |

- **0 to −10 V lie OUTSIDE the validity of the 2.5D model, and it cannot decide them.** These points are dominated by path elements below 1e4 V/cm:
  - there the plasma grows to ~10 µm, past the box and past the scale over which the real field varies;
  - the slice background is intrinsic-depleted, while the real material near the depletion edge is neutral n-type;
  - the R 5 vs 20 µm shapes at 0 V differ (0.76 vs 0.64 at τ/30).

  The low-field treatment was not adjusted to fit them.
- **The column alone at τ/20 gives 0.9085 at 0 V** (data 0.9083). This is a coincidence worth noting, not a fit.

**(d) Identifiability.** Singular values of the scale-projected shape Jacobian in (log10 τ, log10 s, N_D/1e14), σ 1e-3:

| τ | points | singular values | column norms (log τ / log s / N_D) |
|---|---|---|---|
| 1/30 | tail | 20.1, 0.35, 0.002 | 20.1 / 0.05 / 0.38 |
| 1/30 | all | 352, 31, 0.7 | 315 / 2.7 / 159 |
| 1/10 | tail | 5.7, 0.12, 7e-4 | 5.7 / 0.02 / 0.13 |

- The tail constrains ONE combination: the plasma SRH lifetime.
- **The core radius s is NOT constrained.** The tail χ² is 7.6 / 7.7 / 7.7 for s = 1 / 10 / 50 nm at τ/30, and 24.7 / 24.5 / 24.5 at τ/10. The posterior of s equals its prior. **The Fig. 4b data say nothing about the 4H-SiC track core radius (thermalization length).**
- N_D enters only through the low-bias points, which the model cannot represent.

**Plausibility.**

- **Literature.** τ_p ≈ 23 ns (band 17-42 ns) is about 100x below the high-injection (ambipolar) lifetimes of good 4H-SiC epi. Chowdhury et al., arXiv:1011.1039, give "in the range of 2 µs" (only the abstract was opened). It is also below the 500 ns hole lifetime of Verzellesi 2002 (an abstract, cited in lab/microdosimeter LOG §9).
- **The device's own dark current does not test it.** Its `N_t` (1.36e12) is an effective generation RATE in cm⁻³ s⁻¹ (`etna/core/dark_current.py`), not a trap density. Read as midgap SRH generation it gives n_i / G ≈ 4e-21 s, which is not a lifetime: SiC leakage is not thermal SRH generation. The dark current therefore cannot confirm or exclude a short τ.
- **No independent second test.** A 23 ns lifetime would also act in the column (diffusion collection), where on its own it reproduces the 0 V point at τ/20. The two readings are not independent tests.

**Decision.**

- **No preset.** The plasma-SRH reading of the tail needs a TUNED lifetime, 30-100x below good-epi values. It misses 0 to −10 V, and it holds only at σ ≥ 1e-3. That is not the defensible, non-tuned parameter set the task requires.
- **What is exposed.** `run_tilted_alpha_track(tau_scale=...)` exposes the knob. The calibration is recorded here and on the issue.
- **What the tail data constrain**, if the mechanism is plasma SRH: τ_n + τ_p ≈ 23 ns. The statistical 68 % interval is 21-28 ns at σ 1e-3; the systematic band is 17-42 ns.
- **What they do not constrain:** the core radius, the Auger model and carrier-carrier scattering, all of which are flat in bias.
- **Other mechanisms not excluded.** Other mechanisms that give a 1/E-like tail are not excluded by this work. Examples are trapping within the 3 µs shaping (§9: needs implausible trap densities) or a calibration artefact.

## 13. Float64 floor: excess unknowns (now the default)

`run_transient(excess=None)` solves for (dψ, dn, dp) about the dark state. `kernels.assemble_excess` evaluates the SG flux difference with `bern_diff`, which is accurate to 1.1e-14 relative against a 50-digit Decimal reference (x in ±80, d from 1e-14 to 2). It also evaluates the recombination excess without cancellation.

The absolute form stays selectable (`excess=False`). It is used automatically for thermal coupling and `energy=True`.

MESH_FAST results:

| case | CCE absolute | CCE excess | ledger abs → excess | floor abs → excess |
|---|---|---|---|---|
| default stack, 0 V, LET 1 | 0.873557230 | 0.873558014 | 5.1e-6 → 2.5e-8 | 9.6e-7 → 9.2e-14 |
| default stack, 0 V, LET 1000 | 0.923716151 | 0.923716157 | 1.0e-8 → 9.8e-9 | 8.9e-10 → 8.7e-14 |
| default stack, 50 V, LET 1 | 0.999969306 | 0.999969308 | 1.4e-7 → 7.0e-9 | 5.7e-7 → 5.5e-14 |
| default stack, 50 V, LET 1000 | 0.9998684349 | 0.9998684349 | 1.5e-10 → 1.6e-10 | 7.4e-10 → 7.2e-14 |
| micro, 0 V, LET 1 | 0.986214753 | 0.986214960 | 1.5e-6 → 1.9e-8 | 8.2e-7 → 8.0e-14 |
| micro, 50 V, LET 1000 | 0.9998695579 | 0.9998695579 | 1.6e-10 → 7.9e-11 | 7.0e-10 → 6.8e-14 |

- **Agreement.** At LET ≥ 100 the two forms agree to ≤ 6e-9. At LET 1 the absolute CCE was floor-limited, off by up to 8e-7.
- **Ledger.** The excess ledger is now set by the Newton tolerance (1e-9 per step), not by round-off.
- **Cost** (MESH_FAST, unloaded): +20-28 % per event. LET 1 goes from 2.2 to 2.9 s and LET 1000 from 2.8 to 3.4 s. The extra cost is the `bern_diff` terms; the step and LU counts are unchanged.
