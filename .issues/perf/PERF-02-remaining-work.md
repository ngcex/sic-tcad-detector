---
id: PERF-02
title: Remaining work after the batched engine promotion (handoff checklist)
severity: improvement
category: performance
status: open
labels: [performance, handoff, tests]
---

State at hand-off (2026-10-16): branch `fix/audit-v7-top` is pushed to origin; nothing merged to main.
Read PERF-01 first. Audit issues A-*, B-* are still open (D-01, D-04, C-01, C-04 fixed; A-01 partial).

## 1. Finish "complete the suite" (was in progress, interrupted)
- [ ] Slow-test baseline never captured: run `pytest -n auto --dist loadfile -m slow` on HEAD and on
      commit 8428e62 (a git worktree), diff the outcomes. Only the run on HEAD had started (killed).
      Slow-only files: test_api_device, test_api_field, test_alternative_structures (+ slow tests in others).
- [ ] Verify xdist + numba: cold `etna/core/batched/__pycache__` race (first -n auto run was cold and
      passed: 527 passed), warm run time, and whether 10 workers x 10 numba threads oversubscribe
      (try NUMBA_NUM_THREADS=1 for workers). Do not force `-n auto` in addopts (pdb).
- [ ] README "devsim note" (README.md ~line 90) still says run DD modules one file at a time; update it
      with the measured xdist command once the slow run is verified.
- [ ] test_charge_collection_2d still costs ~100 s (one 2D device build). Could be cut by a lighter mesh.
- [x] Leftover agent worktrees under `.claude/worktrees/` removed (2026-10-16).

## 2. Extend the engine (each needs a parity run vs devsim before widening the box)
- [ ] run_field (psi/E at one bias) on the engine: `reverse_sweep(keep_psi=True)` already gives psi.
- [ ] Widen `_ENGINE_BOX` in etna/api/simulation.py (N_A below 1e18, T outside 250-400, forward bias).
- [ ] Batched facades for the other sweeps that rebuild a devsim device per point:
      cce_vs_fluence (needs lifetime/doping damage in the engine), temperature sweep, dark current
      (TAT/SRV not modelled), flash/transient (transient is time-dependent).
- [ ] CCE batch memory: ContinuityBatch holds ~30 (N,B) temporaries (1.2 GB at B=10000); chunk the batch
      if larger sweeps become routine.
- [ ] 2D: no engine path yet (PERF-01 option D). The 2D tests dominate the remaining suite time.
- [x] Fix devsim cv_sweep 0.5 V lattice (A-01) so the fallback path is also exact (tests/test_cv_sweep_lattice.py).

## 3. Audit issues still open (see .issues/audit-v7/issues/*-INDEX.md)
Highest value next: B-01 (stale TAT field in sweep_iv_vs_temperature, leakage 7.6x low), B-02
(_compute_node_efield garbage on 2D), B-03 (microdosimetric_sweep fails with N_D_junction=None),
D-02/B-06 (DESIGN-1 "18 pA" not produced by the model: measured 141 pA), A-02..A-15.

## 4. Housekeeping
- [ ] Open a PR from fix/audit-v7-top into main when ready (not done; user has not asked).
- [ ] lab/gpu_solver stays as the research record (MLX/numpy backends live only there).

## 5. Session 2026-10-16/17 findings (adaptive ramp, slow-suite profile)
- Slow suite before: 19 min (73 tests, only 4 touched the batched engine; MLX is lab-only). The 2D build
  (`create_2d_dd_device`) dominates: 72% of its time is UMFPACK `factor_numeric` (41 ms x ~1800 Newton
  iterations for a 50 um / -50 V build). Full suite now 606 passed in 653 s.
- Done: `ramp_bias(V_step_max=)` adaptive ramp (2 V steps, coarse tolerance on intermediate points, strict
  on the target, snapshot/restore + halving on failure). 2D build 2.2x (50 um 102->47 s, 150 um 390->177 s),
  parity ~1e-9. Enabled in create_2d_dd_device, run_field, flash (2), cce_vs_fluence, cce_vs_bias_at_fluence,
  cce_post_anneal, cce_vs_epi_thickness (1D gains only 1.3-1.9x: the ramp is a smaller share than guessed).
  NOT enabled: sweep_iv_vs_temperature (1.0x), cce_vs_bias and cv_sweep (ramp failure = reported truncation
  voltage), TAT dark-current loops (lagged E-field/gamma make the final state path dependent), transient.
- cce_vs_fluence at 1e14 and cce_post_anneal gave NaN with fixed steps (ramp fails at the first 0.1 V step,
  effective doping < 1e10); adaptive converges. Value 0.9907 has no independent reference yet.
- Stale-device hazard (open, pre-existing): `devsim.solve` is global; a live converged 2D device sits at a
  ~1e-9 noise floor (strict rel 1e-10 never met) and can make other builds fail to converge. Tests worked
  around by deleting devices early and running fixture tests last. Not yet checked: do the facades that do
  not reset devsim (radiation_damage, temperature_sweep, flash, transient) leave devices alive?
- Fixed A-01 in cv_sweep (off-lattice target solved at the lattice point below it): tests/test_cv_sweep_lattice.py.
- Next: SuperLU (or other) replacement of devsim's UMFPACK shim callback (time one factorization first);
  run_field / temperature sweep onto the batched engine; device_list check after facade calls.

## 6. Predictor + split (same session, later)
- Newton start for Potential extrapolated from the two previous converged states in `_ramp_bias_adaptive`
  (carriers: no gain, left alone). Measured factorizations per 2 V step 28 -> ~8. 2D build 50 um 102 -> 11.5 s,
  150 um 390 -> 74 s (V_step_max=5 in create_2d_dd_device); 1D safe sites 3.3-6.6x. Parity unchanged.
- TestReverseBiasConvergence moved to tests/test_device2d_convergence.py. Full suite: 606 passed in 305 s
  (was ~19 min for the slow part alone).
- Facades leave no devsim devices alive (checked 1D and 2D, repeated calls): stale-device item closed.
- Strict rel 1e-10 sits below the ~1e-9 noise floor of 2D solves: strict solves at the target sometimes hit
  the 40-iteration cap and fall back to the relaxed retry (~15% of factorizations). NOT changed: the
  tolerance is a numerical contract (test_api_engine_switch compares at rtol 1e-10); user decision.
- cce_vs_fluence: 5e13 still NaN in both modes (new open failure); 7e13 and 1e14 now converge (no independent
  reference). Grid where fixed steps converge: identical values.
- B-01 plan (read-only agent): make sweep_iv_vs_temperature use dark_current_sweep (per-step TAT refresh,
  fixed steps; raise on empty result). Fixing moves leakage at 300 K from ~18 pA to ~141 pA, so N_t=2.2e13
  ("calibrated ~18 pA") was likely tuned on the stale path; notebook 06 outputs change. Needs a decision on
  recalibrating N_t before touching it.
- Next candidates: run_field and temperature sweep onto the batched engine; SuperLU timing for one
  factorization; remaining slow tests: test_2d_vs_1d_cv_centerline (153 s), test_center_matches_1d (132 s).

## 7. Tolerance, B-01, N_t (2026-10-17, commits 1f6a6bd, 69df110, 0af94c1)
- DC rel tolerance 1e-10 -> 1e-8 via `etna/core/solver_tol.DC_REL_ERROR` at the measured steady-state sites
  (transient, flash_recombination, single_particle, device.py still 1e-10). 2D build 50 um 7.7 s, 150 um 28.9 s
  (from 102/390 s); zero strict failures at 150 um; a live converged 2D device no longer breaks later builds.
  Full suite 606 passed in 225 s.
- B-01 fixed: sweep_iv_vs_temperature -> dark_current_sweep (test vs dark_current_sweep, was 0.13). Notebook 06 now shows
  I(350K)/I(280K) = 1.39e6 ~ n_i ratio 1.37e6 (was 1.5x with the stale field).
- N_t recalibrated 2.2e13 -> 2.8e12 (calibration input: bound "<18 pA up to -60 V", 4 mm2; linear in N_t). Closes D-02/B-06.
  run_dark_current / fluence / sensitivity leakage drops ~7.9x. Notebooks 06/07 re-executed; notebook 20 (feasibility)
  and scripts/create_notebook_20.py NOT regenerated: stale, check before quoting its leakage.
- Open: B-02 (TAT E-field node indexing is 1D-only; add a dimension guard in setup_tat_model after grepping 2D users),
  cce_vs_fluence NaN at 5e13, remaining slow tests (test_2d_vs_1d_cv_centerline 157 s: check device2d._ramp_to),
  ParametricSweep/batch page coverage on the engine, SuperLU timing (lower priority now).
