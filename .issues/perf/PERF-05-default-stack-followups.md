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
