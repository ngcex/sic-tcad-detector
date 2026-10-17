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
