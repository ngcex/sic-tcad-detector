# c-axis low-field mobility default (issue #6 follow-up): log

Started 2026-10-21.

> **Final state (section 12): strict reading.** The default `c_axis` changes the
> electrons only; the holes stay the legacy (isotropic) set because the primaries disagree
> on them. The holes x 0.90 variant described in sections 3-11 (the first version of this
> branch) is the opt-in set `c_axis_holes_0p9`; its numbers there remain valid for it.

Code:

- `etna/core/physics_models.py`: `LowFieldMobility` (two forms, a direction) and the sets
  `c_axis` (new default), `legacy_perp` (previous default), `ishikawa2025` (opt-in), plus
  the `ishikawa2025` high-field set.
- `etna/core/sic_material.py`: `DEFAULT_MOBILITY`, `SiC4H_Parameters.mobility`,
  `mobility_set`, the dispatch in `mobility_caughey_thomas_T`.
- `etna/core/device.py` (`mobility_ct_expression`, `create_sic_device(mobility=)`; no mesh
  code touched), `etna/core/device2d.py`, `etna/core/batched/{material,sweeps,dark}.py`,
  `etna/core/engine2d/device.py`, `etna/core/track/grid.py`, `etna/api/device.py`
  (`DeviceConfig.mobility`), `etna/api/simulation.py` (engine groups), `etna/api/models.py`
  and `etna/core/devsim_models.py` (direction-aware high-field v_sat).

Tests: `tests/test_mobility_caxis.py`. Bench: `lab/mobility_caxis/bench/caxis_bench.py`
(1D facades, about 2 min with the track section), `lab/mobility_caxis/bench/
rebuild_let_tables.py` (the CCE(LET) tables).

## 1. Claims of the first pass, checked against what was opened

| Claim                                               | Verdict                                             | Source opened                                                                                                                                                                                                                                                                                                                                                                                                                  |
| --------------------------------------------------- | --------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| The default set is the PERPENDICULAR-to-c set       | **confirmed, primary**                              | Ayalew thesis node65: Table 3.5 headers are mu_perp,300 (950 / 40 / 1.94e17 / 0.61 / beta -0.5 / gamma -2.40 for electrons; 125 / 15.9 / 1.76e19 / 0.34 / -0.5 / -2.15 holes); eq. 3.81 mu_n,par = 1.2 mu_n, mu_p,par = mu_p                                                                                                                                                                                                   |
| The electron gamma -2.40 is the parallel-to-c value | **confirmed, secondary** (Schaffer 1994 not opened) | review arXiv:2410.06798 Table XXVIII footnote c ("T > 250 K in <11-20> (gamma_max = -2.4 towards <0001>)") and text (Schaffer: T^-2.4 parallel, T^-2.15 perpendicular above 200 K). Roschke & Schwierz 2001 (row Rosc01) label -2.4 perpendicular, so Ayalew's table mixes the two. Measured (Ishikawa thesis Fig. 2.17, Hall, 300-600 K): beta_par 2.58-2.78 and beta_perp 2.40-2.58 at N_D < 1e16; so -2.40 is close to both |
| The hole gamma -2.15 is borrowed from the electrons | **confirmed, secondary**                            | review text: "Schaffer et al. stated gamma_max = -2.15 only for electrons". Numerically close to measurement: hole Hall beta 2.18-2.36 (Ishikawa thesis Fig. 4.7, both directions)                                                                                                                                                                                                                                             |
| Ayalew's (T/300)^-0.5 on mu_min is missing          | **confirmed, primary**                              | Ayalew eq. 3.79 / Table 3.5 beta -0.5. Not an agreed law: Ishikawa's drift fit has +0.8 (electrons) and -2.2 (holes) on mu_min. Irrelevant at the epi doping (mu_min matters above ~1e17)                                                                                                                                                                                                                                      |

## 2. Sources

**Opened, primary:**

- R. Ishikawa, _Low-Field and High-Field Anisotropic Carrier Transport in SiC_, PhD thesis,
  Kyoto University (2025), doi:10.14989/doctor.k26017, open access since 2026-03-23
  (http://hdl.handle.net/2433/296110). Hall bars and pulsed-conductance bridges along
  <0001> and <1-100> on (11-20) epilayers, N_D 5e14-4e18, N_A 8e15-7e18, 120-600 K.
  - Eq. 2.16, electron DRIFT mobility parallel to c (Hall / Hall factor):
    mu = (1190 - 10) t^-3.4 / (1 + (N_D / 1.0e17)^0.54 t^-2.2) + 10 t^0.8;
    Eq. 2.17 perpendicular: (980 - 30) t^-3.2 / (1 + (N_D / 1.2e17)^0.55 t^-1.6) + 30 t^-0.2.
  - Eq. 2.11 / 2.12, the Hall fits (parallel mu_max 1260, perpendicular 1040).
  - Eq. 4.6 / 4.7, HOLE Hall mobility parallel (83 / 20 / 6.4e18 / 0.83) and perpendicular
    (94 / 40 / 6.2e18 / 0.72); the hole drift ratio mu_par / mu_perp is about 0.9 (Sec.
    4.3.1, resistivity ratio, independent of N_A and T).
  - Ch. 3: electron v_sat 1.7e7 cm/s in both directions (Table 3.1, N_D <= 1e16), beta
    0.62-0.87 parallel. Ch. 5: holes do not saturate up to 200 kV/cm (Figs. 5.7, 5.9).
- R. Ishikawa, M. Kaneko, T. Kimoto, Solid State Phenom. 360, 205 (2024),
  doi:10.4028/p-MAe9pZ (CC BY): vertical-SBD drift mobility along c 1070 +- 290 cm2/Vs at
  N_D 9.7e14 (Eq. 2.16 gives 1101).
- T. Kimoto, Jpn. J. Appl. Phys. 58, 018002 (2019): electron mu_par 16 % above mu_perp
  (Hall, r_H = 1; 912 vs 784 at 2.4e16).
- A. Koizumi, J. Suda, T. Kimoto, J. Appl. Phys. 106, 013716 (2009): hole Hall mobility,
  Al-doped, perpendicular to c (CT fit with mu_min = 0; read for context).
- T. Ayalew thesis (TU Wien 2004), nodes 64-65 (sections 3.3, 3.3.1).

**Secondary only** (via the review, arXiv:2410.06798 v8, Tables XXVII, XXVIII, XXX, Figs.
37, 39): Schaffer 1994 (MRS 339, 595), Roschke & Schwierz 2001, Hatakeyama 2003 (hole
ratio 1.15), Mnatsakanov 2001, Lades 2000. Not opened (paywalled): all of these, and
Ishikawa's APEX 14, 061005 (2021), pssb 260, 2300275 (2023), JAP 135, 075704 (2024),
whose numbers are in the thesis. The Kimoto & Cooper textbook values were not opened
(the review's Fig. 37 lists [Kimo14a] at mu_perp/mu_par 0.85).

**Discrepancy recorded:** the review's [Ishi24] perpendicular hole row has mu_min 20
(mu_0 74 = 94 - 20); thesis Eq. 4.7 has mu_min 40. The parallel row agrees (20 / 63).

## 3. Do the primaries agree?

**Electrons, parallel to c, 300 K: yes.**

| N_D (cm^-3)                              | 1e14 | 1e15 | 2.3e15 | 1.2e16 | 1.1e17 |
| ---------------------------------------- | ---- | ---- | ------ | ------ | ------ |
| Ishikawa drift par (Eq. 2.16) = `c_axis` | 1162 | 1099 | 1054   | 905    | 585    |
| Ishikawa drift perp (Eq. 2.17)           | 961  | 916  | 883    | 771    | 516    |
| legacy (Schaffer Hall perp)              | 941  | 915  | 893    | 809    | 573    |
| Ishikawa Hall par (Eq. 2.11)             | 1251 | 1220 | 1193   | 1084   | 757    |
| drift par / drift perp                   | 1.21 | 1.20 | 1.19   | 1.17   | 1.13   |

- The thesis's perpendicular drift fit reproduces the legacy set within 3.5 % up to 5e15:
  the legacy set IS the perpendicular drift mobility, and the c-axis set is that times the
  anisotropy 1.2 (Ayalew eq. 3.81), 1.16 (Kimoto 2019), 1/0.83 (review consensus, Table
  XXVII), all consistent.
- The vertical-SBD drift value (1070 +- 290 at 9.7e14) agrees (1101).
- The thesis's pulsed-conductance mobilities (Table 3.1: 1190 / 1200 / 1170 at 2e15 / 5e15
  / 1e16) used the Hall electron density, so they are Hall-like and match Eq. 2.11, not
  2.16: not an independent check of the drift value.
- **Uncertainty:** the drift values rest on the thesis's electron Hall factor (1.08 at
  1e14 to 1.30 at 1e17). Spread of the opened values at low doping: 1070-1174 (drift),
  i.e. about +-8 %.

**Electron temperature law: no.** Ishikawa's drift fit has t^-3.4 on (mu_max - mu_min),
fitted 300-600 K (the electron Hall factor rises with T); the Hall fits have -3.2 / -2.9
and the single-exponent Hall betas 2.4-2.8; Schaffer's parallel exponent is -2.4. At 450 K
the drift law gives 309.5 cm2/Vs at the epi bulk vs 443.9 with -2.40 (-30 %); below 300 K
-3.4 is an extrapolation that grows fast (2159 at 250 K vs 1819). **Decision:** `c_axis`
keeps the legacy T law (gamma -2.40 on mu_max, mu_min T-independent), which the review
attributes to the parallel direction; the full thesis law is the opt-in `ishikawa2025`.

**Holes: direction ratio yes, magnitude no.** mu_par / mu_perp: 0.9 (Ishikawa drift
ratio, primary), 1/1.15 = 0.87 (Hatakeyama 2003, secondary), 1.0 (Ayalew; outdated per the
review). Magnitude at low N_A, 300 K: Schaffer 124 (Hall perp), Hatakeyama 2003 120,
Koizumi 2009 114 (Hall perp), Ishikawa 94 (Hall perp) / 83 (Hall par) / about 78
(conductance par, Fig. 5.8). **Decision:** `c_axis` holes = legacy x 0.90 (only the agreed
part changes): 111.4 cm2/Vs at the epi bulk, uncertainty -30 / +5 %. The Ishikawa Hall
holes are the opt-in `ishikawa2025` (83.0).

**High field: kept opt-in.** Electrons v_sat 1.7e7 isotropic (Ishikawa) vs lades2000's
2.2e7 perp / 1.64e7 par (Nilsson MC); holes not saturated at 200 kV/cm. The lades2000
parallel hole parameters (v_sat 6.5e6, beta 1.2) with mu0 = 85 give 4.1e6 / 5.2e6 cm/s
at 100 / 200 kV/cm against about 4-5e6 / 5e6 read off thesis Fig. 5.7a (N_A 1e16): they
are consistent with the measurement. With the c-axis default, `high_field_mobility` now
uses the parallel v_sat (`PhysicsModels.vertical_is_parallel`).

## 4. Implementation

- `LowFieldMobility` has `direction` ("c_axis" / "perp"), `form` ("ct" = Ayalew eqs.
  3.77-3.79, "arora" = thesis Eqs. 2.10 / 4.5) and `T_range_K`. `coefficients(T, carrier)`
  gives the T-resolved constants shared by the numpy evaluation and `devsim_expression`
  (same order of operations; `tests/test_mobility_caxis.py` evaluates the devsim string
  against numpy to 2e-15).
- One dispatch point: `mobility_caughey_thomas_T(N, T, carrier, params)` and
  `device.mobility_ct_expression(...)` read `params.mobility`; "legacy_perp" runs the old
  code path on the `SiC4H_Parameters` fields (bit-identical). Every builder that creates
  the parameters takes `mobility=` (`create_sic_device`, `create_sic_2d_device`,
  `make_batch` -> `pb.mobility`, `edge_mobilities`, `dark_sweep`, `generation_cce`, the
  track `TrackDevice.mobility`); devsim device infos carry it in `info["params"]`, and the
  engine twins read it from there (fluence, trap kinetics, dark, FLASH, engine2d).
- `DeviceConfig.mobility` (validated) is emitted by `device_kwargs` only when not the
  default, the `incomplete_ionization` pattern; the engine groups key on it.
- The historical presets (`legacy_config()`, `nplus_calibrated_config()`) carry
  `mobility="legacy_perp"`: they reproduce the numbers of their period.
  `measured_cv_config()` (= the default) and `microdosimeter_config()` take the default.
- `PhysicsModels`: "etna_default" = the device's own set; an `anisotropy` is refused on a
  c-axis set; `high_field_params(T, device_mobility)` takes the parallel v_sat whenever
  the vertical transport is parallel to c. The devsim twin (`devsim_models`) and the
  facades pass the device's set.
- The app sidebar has a "Low-field mobility" selectbox.
- fv2d / alternative structures take the default through `SiC4H_Parameters()` (fv2d
  `build_mesh(params=)` selects another); their devsim baselines
  (`lab/geometry2d/notes/devsim_baselines.json`, `devsim_psi.npz`) were regenerated with
  `lab/geometry2d/bench/devsim_baselines.py`: fv2d parity holds (tests/test_fv2d.py).

## 5. Before / after (default dosimeter; `caxis_bench.py`, 2026-10-21, M4)

before = `legacy_perp`, after = `c_axis`.

**(0) Mobility (cm2/Vs).** Epi bulk 3.38e13 / junction 8.19e15 / 1e19, 300 K; then the
bulk at 250 / 350 / 450 K.

|                | 3.38e13 | 8.19e15 | 1e19  | 250 K  | 350 K | 450 K |
| -------------- | ------- | ------- | ----- | ------ | ----- | ----- |
| e before       | 945.4   | 834.7   | 115.3 | 1464.2 | 653.1 | 357.4 |
| e after        | 1174.4  | 947.4   | 100.6 | 1819.1 | 811.3 | 443.9 |
| e ishikawa2025 | 1174.4  | 947.4   | 100.6 | 2159.0 | 703.4 | 309.5 |
| h before       | 123.8   | 117.5   | 75.7  | 183.1  | 88.9  | 51.9  |
| h after        | 111.4   | 105.8   | 68.1  | 164.8  | 80.0  | 46.7  |
| h ishikawa2025 | 83.0    | 82.8    | 45.7  | 125.7  | 58.4  | 33.0  |

**(1) C-V:** identical (max |C_after / C_before - 1| = 0). Every C-V fit, the Photons Fig.
4b W labels (residuals -0.0253 / -0.0015 / +0.0203 um both ways) and the Bruzzi C-V check
are unchanged.

**(2) alpha CCE (`run_cce`, batched; devsim agrees to 1e-6):**

| V            | -2       | -5       | -10      | -20      | -40      | -70      | -100     |
| ------------ | -------- | -------- | -------- | -------- | -------- | -------- | -------- |
| before       | 0.993707 | 0.998831 | 0.999429 | 0.999703 | 0.999845 | 0.999910 | 0.999936 |
| after        | 0.993031 | 0.998715 | 0.999374 | 0.999676 | 0.999831 | 0.999901 | 0.999930 |
| ishikawa2025 | 0.990714 | 0.998305 | 0.999179 | 0.999577 | 0.999781 | 0.999872 | 0.999910 |

Decomposition (bench-only hybrid sets, CCE minus before at -2 / -10 / -40 / -100 V):
c-axis electrons only +8.5e-6 / +4.3e-6 / +1.5e-6 / +6.8e-7; c-axis holes only -6.9e-4 /
-6.0e-5 / -1.6e-5 / -6.6e-6; both -6.8e-4 / -5.5e-5 / -1.4e-5 / -5.7e-6; with the legacy
electron mu_min 40 (p+ minority electrons) identical to 1e-8. **The CCE loss is the hole
diffusion out of the undepleted epi tail; the electron change is irrelevant there.**

**(3) Full-width drift transit time (low field):**

| V                | -10    | -40   | -70   | -100  |
| ---------------- | ------ | ----- | ----- | ----- |
| t_e before (ps)  | 125.1  | 34.0  | 17.6  | 11.3  |
| t_e after        | 103.6  | 28.6  | 14.7  | 9.4   |
| t_h before       | 935.6  | 249.5 | 130.3 | 84.6  |
| t_h after        | 1039.5 | 277.2 | 144.7 | 94.0  |
| t_h ishikawa2025 | 1376.2 | 364.4 | 191.0 | 124.6 |

t_e x0.826-0.841 (the doping-weighted ratio of the depleted region), t_h x1.1111.
With `lades2000` velocity saturation the electron transit becomes 149.0 / 82.1 / 70.6 /
66.4 ps (x1.44 / 2.87 / 4.80 / 7.10 over low field; parallel v_sat 1.64e7 now used), the
hole 1130 / 391 / 268 / 223 ps; the previous opt-in route (legacy_perp + burin2024 +
lades2000) gave 149.2 / 81.7 / 70.5 / 66.4 ps and 1166 / 400 / 273 / 226 ps: the c-axis
default reproduces it within 1 % (electrons) and 3 % (holes). With the `ishikawa2025`
high-field set (v_sat 1.7e7, beta 0.76): 197.5 / 105.0 / 86.2 / 78.3 ps.

**(4) Dark current (4 mm2):** 7.637e-12 / 1.478e-11 / 1.800e-11 / 2.487e-11 A at -20 /
-40 / -60 / -100 V both ways (relative change < 2e-6; -3.3e-3 at 0 V, 1.16 pA). The
current is TAT-generation dominated, so the 18 pA N_t calibration holds; no refit.

**(5) CCE vs T, -30 V:** dd 0.999813 / 0.999797 / 0.999783 / 0.999757 -> 0.999795 /
0.999778 / 0.999763 / 0.999734 at 250 / 300 / 350 / 450 K; hecht 0.995046 / 0.994405 /
0.993800 / 0.992675 -> 0.995986 / 0.995467 / 0.994977 / 0.994065.

**(6) Tilted-alpha CCE, microdosimeter preset (Petringa Micro Fig. 4b):** raw 0.87545 /
0.98042 / 0.99726 / 0.99928 / 0.99961 at 0 / -1 / -2 / -5 / -10 V -> 0.86410 / 0.97833 /
0.99697 / 0.99921 / 0.99958. Normalised at -50 V against the digitised data: rms 0.0354
-> 0.0351, residual at 0 V +0.0845 -> +0.0837 (data 0.9083, model 0.9928 -> 0.9920),
residuals at >= -5 V unchanged to 1e-4. **Marginally better, not meaningfully:** the
model stays flat where the data drop by 9 % (plasma / columnar recombination not in the
DD model; LOG of lab/microdosimeter).

**(7) FLASH steady-state CCE, -30 V:** 0.999809 -> 0.999791 at every dose rate.

**(8) Proton damage CCE (`run_radiation_damage`, -40 V, legacy damage model):**
0.99982 / 0.99979 / 0.99973 / 0.99959 / 0.99930 / 0.99874 -> 0.99981 / 0.99977 /
0.99970 / 0.99956 / 0.99926 / 0.99870 at 1e10 / 2.3e10 / 5.2e10 / 1.2e11 / 2.6e11 /
6e11 p/cm2 (NaN beyond compensation both ways). With the charged traps (opt-in):
0.996109 / 0.972203 / 0.354958 / 0.152176 -> 0.996191 / 0.973500 / 0.353762 / 0.151621
at 1e13 / 1e14 / 1e15 / 5e15: the faster electrons win where trapping limits the
collection (1e13-1e14), the slower holes lose elsewhere; all within 1.3e-3.

**(9) FLASH transient CCE (`run_transient`, devsim, -30 V):** 0.999725 / 0.999868 /
0.999869 / 0.999870 / 0.999801 / 0.999822 -> 1.000000 (raw 1.0000022) / 0.999758 /
0.999759 / 0.999842 / 0.999822 / 0.999792 at 20 / 50 / 100 / 150 / 200 / 230 Gy/s. The
moves (+2.8e-4 to -1.1e-4, not monotonic in dose rate) are of the order of the
direct-anode transient estimator's known round-off (~1e-4, tests/test_default_stack.py),
not resolved physics; the steady-state equivalent moves by -1.8e-5.

**(10) Single-ion track engine (opt-in, MESH_FAST, gaussian 50 nm, 50 V):** LET 10: CCE
0.999972 -> 0.999976, t95 of the cathode charge 121.4 -> 132.6 ps; LET 1000: CCE 0.999868
-> 0.999860, recombined 1.317e-4 -> 1.396e-4, t95 232.3 -> 268.1 ps (+15 %, the slower
holes). Ledger residuals 1e-8 / 2e-10.

**Microdosimetric spectra** (`run_microdosimetry`): independent of the mobility (pure MC
pipeline). The CCE(LET) tables: section 6.

## 6. CCE(LET) tables (`rebuild_let_tables.py`, notebook-16 keywords, 6 / 4 workers)

| LET (keV/um)                      | 0.5       | 4.43      | 39.2      | 500       |
| --------------------------------- | --------- | --------- | --------- | --------- |
| tracked table before (older code) | 0.9999682 | 0.9999718 | 0.9999723 | 0.9999723 |
| legacy_perp, rebuilt on this code | 0.9999625 | 0.9999690 | 0.9999701 | 0.9999703 |
| c_axis (now tracked in data/)     | 0.9999664 | 0.9999722 | 0.9999735 | 0.9999737 |

- hw 50 um (100 um table), 50 V. t_collection (95 % of the terminal charge): 0.1168 ns
  (legacy_perp) -> 0.1298 ns (c_axis, +11.1 % = 1 / 0.90: the holes set it); the old
  tracked column was the stale 0.208 ns.
- 300 um table (hw 150 um, 10 LETs): 0.9999568 / 0.9999707 / 0.9999732 / 0.9999737 at 0.5
  / 5 / 50 / 500 keV/um (tracked before 0.9999594 / 0.9999704 / 0.9999718 / 0.9999719).
- Build times 782 s (100 um) and 1293 s (300 um) with 6 workers; legacy 100 um 1214 s with 4. All points converged. The legacy rebuild is kept in `notes/let_legacy_perp/`.
- Every CCE moves by < 4e-6: the tables are low injection, fully depleted at 50 V, and the
  2D sheet track has no plasma model, so the mobility only sets the timing.

## 7. Validation comparisons: better or worse?

- **Bruzzi 2023 C-V** and **Photons Fig. 4b W labels**: C-V only, identical.
- **Petringa Micro Fig. 4b (tilted alpha)**: rms 0.0354 -> 0.0351, 0 V residual +0.0845 ->
  +0.0837. Marginally closer, far below the 9 % discrepancy that the DD model cannot
  produce. Not evidence for or against the mobility.
- **18 pA dark-current bound** (N_t calibration): unchanged to 2e-6, no refit needed.
- **Calibration circularity check:** the knobs fitted on the legacy mobility are the C-V
  doping (independent of the mobility) and N_t (dark current insensitive to it). The
  damage constants (K_tau, removal) are literature values, not refitted. So no knob
  tuned on the legacy mobility moves a validation number here.
- **Breakdown (Hatakeyama 2004)**: unchanged (impact ionization depends on the field only).

## 8. Open / not done

- The hole magnitude parallel to c: 83 (Ishikawa Hall) to 112.5 (default); a minority-
  carrier hole mobility in n- SiC (the dosimeter's transport) is not measured by Hall.
  The opt-in `ishikawa2025` brackets the low side: CCE at -2 V 0.990714 vs 0.993031.
- The electron T law above 300 K: Ishikawa's drift fit (-3.4) gives 30 % less electron
  mobility at 450 K than the kept -2.40; temperature outputs carry that uncertainty.
- No mobility tensor in 2D (lateral electrons ~20 % too mobile with the c-axis scalar).
- High-field: `ishikawa2025` (v_sat 1.7e7, isotropic) registered, not the default (no
  high-field model is the default); the parallel `lades2000` electron v_sat (1.64e7, MC)
  is within 4 % of it.
- The FLASH transient estimator noise (~1e-4) hides the mobility effect there.
- fv2d / alternative structures (notebook 19) take the default; no `DeviceConfig` route.

## 9. Re-pinned tests, regenerated files

Re-pinned to the c-axis numbers (the legacy value stays in a comment or a second assert):
`tests/test_nplus.py` (1D spike CCE, 2D CCE with/without n+ generation),
`tests/test_default_stack.py` (default CCE record, transient anode), `tests/test_charge_collection.py`
(Hecht 300 K), `tests/test_material.py` (legacy formula tests now take legacy params;
device mobility 1096.28 / 109.01), `tests/test_microdosimeter_preset.py` (0 V 0.9920, old
30 V doping 0.9583), `tests/test_viz_figures.py` (Fig. 4b model 0.9920),
`tests/test_radiation_traps.py` (charged-trap pins), `tests/test_physics_models.py`
(ayalew2004 vs legacy). `tests/test_displacement_current.py` uses `legacy_perp` on its
uniform-field device (the c-axis mu_n / mu_p ratio leaves 3 time steps in its electron-
only window; the q N v / d check is set-agnostic). `tests/test_app_device_sidebar.py`:
the new field. Regenerated: `data/cce_let_table_{100,300}um.json`,
`examples/measurements/data/it_synthetic.csv` (make_synthetic.py; I(T) moves by < 5e-6),
`lab/geometry2d/notes/devsim_baselines.json` / `devsim_psi.npz` (devsim_baselines.py).
Notebooks carry no committed outputs; nb14's first-order transit times (78 / 285 ps, high-
field velocities) are not ETNA mobility numbers and were left.

## 10. Mesh convergence of the change (`section_convergence`, batched, h / h/2 / h/4)

|                   | -2 V       | -10 V      | -40 V      |
| ----------------- | ---------- | ---------- | ---------- |
| after - before, h | -6.7608e-4 | -5.4935e-5 | -1.4082e-5 |
| h/2               | -6.7760e-4 | -5.4906e-5 | -1.4054e-5 |
| h/4               | -6.7838e-4 | -5.4901e-5 | -1.4044e-5 |

The change is mesh-converged to 0.3 % (-2 V) and 0.06 % (-10 V); the CCE itself moves by
2e-5 (-2 V) and 1e-7 (-10 V) from h to h/4, both sets alike.

## 11. Decision points for the reviewer

- **Every moved CCE pin comes from the hole x 0.90**, not from the electrons
  (decomposition, section 5 (2)): c-axis electrons alone move the CCE by at most +8.5e-6,
  so an electrons-only variant would keep every legacy pin within 1e-5. The hole ratio
  rests on one opened primary (Ishikawa thesis: resistivity ratio about 0.9; its Hall
  83 / 94 = 0.88 and conductance about 0.82 agree in sign) and one secondary
  (Hatakeyama 2003, 1/1.15); Ayalew eq. 3.81 (primary, older) says isotropic. If the
  stricter reading of "only if the primaries agree" is wanted, set `mu_max_p` /
  `mu_min_p` of "c_axis" back to the legacy 125 / 15.9 (one line); the electron change
  alone then reproduces the legacy CCE pins to < 1e-5 while transit times keep t_e x0.83.
- **Merge hazard:** `lab/geometry2d/notes/devsim_baselines.json` / `devsim_psi.npz` and
  `data/cce_let_table_*.json` depend on the 2D meshes that a parallel branch changes.
  Whoever merges second reruns `lab/geometry2d/bench/devsim_baselines.py` (minutes) and
  `lab/mobility_caxis/bench/rebuild_let_tables.py data_out c_axis 6` (13 + 22 min) and
  copies the tables to data/. The npz cannot be hand-merged.
- **ishikawa2025 below 300 K** warns (`MobilityExtrapolationWarning`), not refused.
- **lab scripts** that read `SiC4H_Parameters` mobility fields directly
  (`lab/solver2d/solver2d/engine.py`, `lab/microdosimeter/mdlab.py`, bench/hypotheses*)
  still compute the perpendicular (legacy) set; they are frozen lab records.

## 12. Final state: strict reading (coordinator decision, 2026-10-21)

`c_axis` = Ishikawa Eq. 2.16 electrons (T law gamma -2.40) + the legacy holes (isotropic,
as Ayalew eq. 3.81): the primaries agree on the electrons only. The holes x 0.90 set is
the opt-in `c_axis_holes_0p9`; `ishikawa2025` (Hall holes) stays opt-in. Merged with main
(contact-cell graded meshes for alternative structures, trap-kinetics follow-up); the
devsim baselines (`depth_mesh="legacy"`, as tests/test_fv2d.py pins) and the CCE(LET)
tables were regenerated once on the merged code.

| Output (default dosimeter)                  | before (`legacy_perp`)                    | after (`c_axis`)                          |
| ------------------------------------------- | ----------------------------------------- | ----------------------------------------- |
| mu_n / mu_p at the epi bulk, 300 K (cm2/Vs) | 945.4 / 123.8                             | 1174.4 / 123.8                            |
| C-V, V_BD, Bruzzi, Photons W labels         |                                           | identical                                 |
| alpha CCE -2 / -10 / -40 / -100 V           | 0.993707 / 0.999429 / 0.999845 / 0.999936 | 0.993716 / 0.999433 / 0.999847 / 0.999937 |
| CCE change, mesh h / h/2 / h/4, -2 V        |                                           | +8.519e-6 / +8.543e-6 / +8.547e-6         |
| t_e -10 / -40 / -100 V (ps)                 | 125.1 / 34.0 / 11.3                       | 103.6 / 28.6 / 9.4                        |
| t_h (ps)                                    | 935.6 / 249.5 / 84.6                      | unchanged                                 |
| t_e with lades2000, -100 V                  | 52.6 ps (perp v_sat)                      | 66.4 ps (parallel v_sat 1.64e7)           |
| dark current -20 / -60 / -100 V             | 7.637e-12 / 1.800e-11 / 2.487e-11 A       | change < 2e-6                             |
| CCE vs T, dd, -30 V, 250 / 300 / 450 K      | 0.999813 / 0.999797 / 0.999757            | 0.999814 / 0.999799 / 0.999759            |
| Hecht CCE vs T, 300 K                       | 0.994405                                  | 0.995475                                  |
| FLASH steady state, -30 V                   | 0.999809                                  | 0.999811                                  |
| FLASH transient 20 / 100 Gy/s               | 0.999725 / 0.999869                       | 0.999725 / 0.999869                       |
| proton damage (legacy model) 2.6e11 / 6e11  | 0.99930 / 0.99874                         | 0.99932 / 0.99880                         |
| charged traps 1e13 / 1e14 / 1e15 / 5e15     | 0.996109 / 0.972203 / 0.354958 / 0.152176 | 0.996417 / 0.974836 / 0.357614 / 0.152972 |
| tilted alpha, microdosimeter, 0 V (raw)     | 0.87545                                   | 0.87545                                   |
| Petringa Fig. 4b rms                        | 0.0354                                    | 0.0354                                    |
| track engine t95, LET 10 / 1000             | 121.4 / 232.3 ps                          | 120.7 / 237.2 ps                          |
| CCE(LET) 100 um, 0.5 / 500 keV/um           | 0.9999625 / 0.9999703 (pre-merge)         | 0.9999662 / 0.9999742                     |
| CCE(LET) t_collection                       | 0.1168 ns                                 | 0.1168 ns                                 |

Every legacy CCE pin holds within 1e-5 (largest +8.5e-6 at -2 V). What moves: electron
transit (x0.83), the trapping-limited CCE (charged traps up to +2.6e-3 at 1e14; Hecht
+1.1e-3), the velocity-saturated electron transit (parallel v_sat). Re-pinned again
(legacy values in comments): tests/test_nplus.py, test_default_stack.py,
test_charge_collection.py, test_material.py, test_microdosimeter_preset.py,
test_viz_figures.py, test_radiation_traps.py, test_mobility_caxis.py; examples
it_synthetic.csv regenerated.

After the merge also: tests/test_contact_mesh.py (mesa references re-derived with
`lab/mobility_caxis/bench/contact_mesh_refs.py`: 0.963068316 / 0.999091802 converged,
0.96924012 / 0.99987739 legacy mesh; legacy_perp reproduces the previous pins to 5e-9 /
4e-8) and tests/test_balanced_current.py (legacy geometry now with `legacy_perp`).
Full suite on the merged branch: 1731 passed, 8 skipped.
