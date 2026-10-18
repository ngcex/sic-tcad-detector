# PERF-05: follow-ups after the default-stack flip and core cleanup (2026-10-17)

From the independent review of 72beb9b (no blockers):

- `calibrate_graded_doping` (device.py) and `calibrate_graded_doping_2d` (device2d.py) default
  to p+ 1 um and do not pass nplus, so they now build p+ 1 um + default n+ (neither model).
  Callers: notebook 02, scripts/run_calibration_2d.py. Decide: calibrate against legacy
  geometry explicitly, or against the stack.
- `alternative_structures` builders and low-level batched helpers (`node_positions`,
  `make_batch`) still default to legacy geometry: "every core default reads DEFAULT_STACK"
  is overstated.
- `test_cce_increases_with_bias` is pinned to legacy until the single-particle estimator
  (Option A, advisor) lands; `test_v3_cce_center_100um_preserved` now guards the DEFAULT
  model (0.95346904), name is misleading.
- `microdosimetric_sweep(n_jobs>1)` spawn workers do not pin numba/BLAS threads
  (oversubscription); no app path uses n_jobs today.
- `cce_heatmap_2d` counts the n+ layer as active epi (assigned to the estimator agent).

## n+ mesh coarsening (lab, 2026-10-17): gate failed, production unchanged

The 2x node count is in the EPI (devsim grades 1 nm cells at the epi/n+ line back to the 50 nm
epi spacing, ~330 extra rows on a 10 um epi), not in the n+. Candidates move W by 2.4e-5
(k3, -14 % nodes, 0.80x devsim) to 3.5e-4 (e05n02, -39 %, 0.51x); the fine mesh itself is only
converged to 5e-6..3.5e-5 in W. Decision (engine2d already serves the app in ~1 s per 2D
point; devsim is the fallback/validation path): keep the fine mesh. Revisit if devsim 2D
builds become a bottleneck. Data: lab/solver2d/notes/nplus_mesh_gate.csv, LOG.md.

## -30 V C-V gap (+2.56 %): fringe test (lab/finite_implant, 2026-10-17)

2D finite p+ implant slice (devsim): fringe capacitance is ~1.5 % of the plane C at -30 V
(2.2 % peak at reach-through), not 0.1-1 %. Refitting 0/-10 V with it closes 60 % of the gap
with no charge; it closes it exactly only with oxide charge Qf ~ +2e12 cm^-2 (high for
SiO2/4H-SiC), plates make -10 V worse. NOT ruled out. Four explanations each reproduce the
three paper points with one knob: fringe+Qf, 530 nm buffer, 9.75 um epi, ~2.5 % larger area;
they separate only with the full measured C-V (knee -5..-15 V; eps*A/C at -5 V: buffer 6.55,
9.75 um epi 6.36, fringe 7.08 um) and C-V on two die sizes. Needs data from the owner.
Also open: temperature_sweep.py:116-119 and charge_collection.py:135 still use N_ref*0.001 for
mobility doping; figures/cv_validation_bruzzi2023.png stale.

## -30 V C-V gap: the full measured curve (lab/cv_gap, 2026-10-18)

Photons Fig. 4b digitised at every 250 mV (0.02 um; Microdosimetry Fig. 3b is the same data
set, agrees to 0.02 um). The knee is at -1..-4 V, not -8..-10: W(-2 V) = 7.18 um (production
2.99), W(-5 V) = 9.27 (5.93). All four 3-point candidates and the production stack are excluded
(chi2 1e5-2.7e5 for 230 points); after refitting the epi profile none of the single knobs
(fringe, buffer-to-1e19, area, epi thickness) fits the -10..-60 V creep either (chi2 >= 1900).
What fits (chi2 198/224): junction layer ~8e15 / L 0.46 um, bulk ~3.4e13, and a graded n-/n+
transition still at only ~1e16 at 10 um. The production doping FORM cannot represent the data;
proposal (not applied) in lab/cv_gap/notes/LOG.md section 6. Bruzzi & Verroi: different wafer,
no constraint on N_D_bulk; its +3.1 % is the Schottky-vs-pn V_bi.

## Mott-Schottky cross-check of the digitised C-V (2026-10-18)

1/C^2 slope of the digitised Photons Fig. 4b points (C = eps*A/W, A = 4 mm2): -3..-0.5 V: N = 3.9e13;
-6..-2 V: 1.2e14 (this window gives the authors' "order 1e14"); -10..-3 V: 4.8e14; -30..-10 V: 5.7e15.
So the authors' "order 1e14, consistent with nominal" is a loose slope over the knee window; the
digitisation is consistent with their statement and the profile is far from uniform (doping rises
steeply toward the epi/n+ end). Bruzzi & Verroi is a Schottky on a 30 um epi at 2.46e15: no
constraint on this device (see VALIDATION_LITERATURE.md, corrected).
Decision pending (user): adopt an opt-in graded n-/n+ transition (NJ 8.2e15, L 0.46, NB 3.4e13,
Lb 1.14 um, N_end 4.5e16) with gates; flipping the default puts N_D_bulk below every documented
doping statement and the epi/n+ position is no longer fixed by the data.
Cheapest discriminating measurements: C-V to -100...-200 V (models predict eps*A/C 9.96-9.98 um
at -100 V, 10.09-10.11 um at -200 V) or SIMS/SRP through the epi/n+ interface.

## Opt-in measured-C-V preset implemented (feat/buffer-profile, 2026-10-18)

`etna.measured_cv_config()` / `DeviceConfig(buffer_thickness_um, buffer_end_doping)`: graded
n-/n+ transition in the devsim builders (1D/2D), the batched engine and the 2D engine
(parity 1D dW 1.9e-15, dpsi 1.4e-14 V, dCCE 1.2e-10; 2D dpsi 7.1e-14 V, dCCE 1.4e-14), engine
boxes widened to N_D_bulk 2e13 and the buffer keys, sidebar checkbox, export labels, N_t
refit 2.849024e12. Plateau rms 0.013 um vs the digitised curve; the knee stays at 0.24 um rms
(the exponential form; a tabulated profile would be needed). DEFAULT_STACK unchanged; impact
table and open items in lab/buffer_profile/notes/LOG.md.

## DEFAULT FLIPPED to the measured-C-V profile (feat/default-measured-cv, 2026-10-18)

Owner-authorised (Fable review). `DeviceConfig()` / `DEFAULT_STACK` = the former
`MEASURED_CV_STACK` (now an alias; `measured_cv_config()` == `DeviceConfig()`): p+ 0.3 um,
NJ 8.185880e15, NB 3.382506e13, L 0.4638630 um, graded n-/n+ transition 1.136634 um to
4.516421e16, n+ 1 um 1e19, N_t 2.849024e12 (`SiC4H_Parameters.N_t`), E_t 0.67. Previous
default kept bit-identical as `etna.nplus_calibrated_config()` (`NPLUS_CALIBRATION`, now its
own record); `legacy_config()` unchanged. Core builders' `buffer_*` defaults, the engine2d
gate (missing key = default transition), `microdosimetric_sweep` and `run_all` read the
stack; `device_kwargs` always emits the buffer keys (None = abrupt epi end).
Reproducibility (lab/default_flip/bench/snapshot.py + compare.py, 4 mm2): new default ==
old preset, nplus_calibrated == old default, legacy == legacy, all worst rel diff 0.0
(donors, C-V, CCE, dark, peak E, compensation fluence, sweep gate); direct core calls
(`create_sic_device(graded)`, `cce_vs_bias()`) == `DeviceConfig()` (donors to 9.3e-16 rel, CCE
2.6e-11 devsim vs batched).

| quantity (4 mm2, 300 K)                    | old default (nplus_calibrated)       | new default (measured C-V)          |
| ------------------------------------------ | ------------------------------------ | ----------------------------------- |
| V_fd analytic / 1.2 V_fd gate              | 9.93 V / 11.92 V                     | 2.74 V / 3.28 V                     |
| reach-through (VFD_DEFINITION)             | 12.30 V                              | 3.10 V                              |
| eps*A/C 0 / -2 / -5 / -10 / -30 / -60 V um | 1.700 2.992 5.931 9.500 9.979 9.996  | 1.655 6.939 9.252 9.519 9.750 9.872 |
| CCE 0 / -1 / -2 / -3 / -5 / -10 V          | 0.8407 0.8931 0.9331 0.9612 0.9893 0.9993 | 0.8703 0.9743 0.9944 0.9978 0.9989 0.9995 |
| CCE -50 V                                  | 0.999887                             | 0.999883                            |
| dark current -10 / -30 / -60 V (pA)        | 6.08 / 14.72 / 18.00                 | 6.30 / 15.76 / 18.00                |
| peak E excl. step spikes, -30 / -60 V      | 8.31e4 / 1.13e5 V/cm (x 0.40 um)     | 9.06e4 / 1.21e5 V/cm (x 0.40 um)    |
| compensation fluence 5.6 / 62 MeV p        | 3.30e13 / 4.72e13 p/cm2              | 1.35e13 / 1.93e13 p/cm2             |
| run_all sweep grid (hw 25/50 x 20/40/60 V) | 6/6 valid, engine2d                  | 6/6 valid, engine2d                 |
| microdosimetric gate by bias, hw 25        | 5, 10 V shallow; valid from 15 V     | valid from 5 V (2 of 7 rows change) |
| 2D engine, one point hw 50 / 40 V          | 0.21 s, CCE 0.999860                 | 0.22 s, CCE 0.999855                |
| 2D default mesh hw 50 (devsim build)       | 17784 nodes, 0.05 s                  | 17784 nodes, 0.05 s                 |
| run_all, all 8 sections                    | 11.0 s                               | 11.1 s                              |
| run_cv 241 pts / sweep grid                | 0.12 s / 0.91 s (batched / engine2d) | 0.17 s / 0.91 s (batched / engine2d)|

DESIGN-1 operating bias: gate >= 3.3 V; the paper's -10 V is 3x the gate (3.6x V_fd), -30 V
9x (11x V_fd); recommendation kept at >= -10 V (margin for the fit's doping uncertainty:
V_fd 3.9-7.6 V if the bulk is at the documented 0.5-1e14). Export labels state: fit to the
digitised Fig. 4b, plateau rms 0.013 um, knee rms 0.24 um, N_D_bulk 3.4e13 below the
documented 0.5-1e14, epi/n+ position set by the fit, not by SIMS.
Decisions taken here: `calibrate_graded_doping(_2d)` now build the LEGACY geometry
explicitly (item 1 above); `alternative_structures` graded fallbacks pinned to the
previous default's doping (legacy geometry; the guard-ring diode does not converge past
3 V with NB 3.4e13 / L 0.46 um and no transition); notebook 19 pins the same.
Fast paths still serve the default: `validated_domain_reason(DeviceConfig())` None (run_cv /
run_cce batched), engine2d gate None for the default 2D kwargs (run_all sweep, field).
Parity re-run on the new default (lab/buffer_profile/bench, no-buffer cases now explicit):
1D 11 cases to -100 V worst dW 1.9e-15, dpsi 1.4e-14 V, dCCE 1.2e-10, donors 7.3e-15;
2D 9 cases worst dpsi 1.5e-13 V (20 um epi, 100 V; 7.1e-14 before), dCCE 1.5e-14,
engine 0.1-0.3 s vs devsim 4.5-27.5 s per case. v3 cce_center re-pinned 0.95377414 ->
0.95394391 (+0.018 %, record in tests/baselines/v3_frozen.json).
## 2D CCE contact-node normalisation (2026-10-18)

compute_cce_2d divided by the box-rule integral of G over ALL nodes, but contact nodes are
Dirichlet (never injected). Now: `injected_generation` (contact nodes zeroed) is the
normaliser of compute_cce_2d and contact_current_sharing; alpha_stripe_generation /
cce_lateral_scan default zero_on_contacts=True; engine2d `_generation` zeroes the anode and
cathode rows. Shifts (lab/default_flip/bench/cce2d_contacts.py, devsim == engine):
default stack alpha stripe (hw 25, 20 / 50 V) 0.9997234 / 0.9998827: unchanged (contact
share 0: the stripe is confined to the epi); default + nplus_generation 50 V 0.962262 ->
0.962708 (share 4.64e-4); legacy planar hw 25 50 V 0.998545 -> 0.999804, hw 50 20 V 0.996396
-> 0.997653 (share 1.26e-3); uniform generation through the n+ (v3 recipe) 0.954383 ->
0.954963 at 50 V (share 6.07e-4); v3 cce_center 0.95394391 -> 0.95452362; test_nplus 2D
abrupt-n+ nplus_generation 0.9621205 -> 0.9625666. Engine vs devsim after the change: max
|dCCE| 7.9e-15 over these cases (<= 1e-11 gate).
Notebooks re-run into scratch (base 68af67e vs this branch; tracked notebooks untouched):
nb06 75 s: I_reverse(300 K) 3.680e-10 -> 3.941e-10 A/cm2 (+7.1 %: N_t refit + profile),
W(0 V, 300 K) 1.689 -> 1.649 um (edge), CCE(300 K) 0.9942 -> 0.9944, dI/dT clinical
4.758e-10 -> 5.095e-10 A/cm2/K. nb07 55 s: N_t 2.81e12 -> 2.85e12; its own -30 V device
check 18.4 -> 19.7 pA (ratio to 18 pA 1.02 -> 1.09; the notebook compares a -30 V point
with the -60 V bound). nb12 31 s: compensation now at 1.94e13 (was 5e13 in its grid):
the 1.94e13 grid point is reported as fully compensated (not solved); the 2.78e13 /
3.73e13 / 5e13 points were already. nb20 259 s: sweep valid 28 -> 32 of 36 (V_fd column
6.96 / 10.16 -> 3.91 / 6.47 V at N_D_bulk 5e13 / 8.5e13), baseline half-width points 7 ->
8, center CCE 20 V 0.999726 -> 0.999723, dark at 50 V (hw 50) 4.500e-14 -> 4.441e-14 A.
Full suite: 913 passed in 290 s (-n auto, loadfile).
Single-particle (etna.core.single_particle) aligned to the same convention: the ion track
zeroes the contact nodes and Q_generated is the box-rule integral (was the triangle rule
over all nodes). Legacy hw 25, 10 V, LET 1: 0.961613 -> 0.962940 vs steady-state 0.962560
(gap -9.5e-4 -> +3.8e-4); default hw 25, 50 V (anode): 0.999970 -> 0.999966 (box vs
triangle only; the track touches no contact with the n+ layer), so nb16's default-stack
CCE(LET) table moves by ~4e-6 (lab/default_flip/bench/sp_contacts.py).
