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
      commit 98a14b4 (a git worktree), diff the outcomes. Only the run on HEAD had started (killed).
      Slow-only files: test_api_device, test_api_field, test_alternative_structures (+ slow tests in others).
- [ ] Verify xdist + numba: cold `etna/core/batched/__pycache__` race (first -n auto run was cold and
      passed: 527 passed), warm run time, and whether 10 workers x 10 numba threads oversubscribe
      (try NUMBA_NUM_THREADS=1 for workers). Do not force `-n auto` in addopts (pdb).
- [ ] README "devsim note" (README.md ~line 90) still says run DD modules one file at a time; update it
      with the measured xdist command once the slow run is verified.
- [ ] test_charge_collection_2d still costs ~100 s (one 2D device build). Could be cut by a lighter mesh.
- [ ] Remove the leftover merged worktrees under `.claude/worktrees/` (testpaths now keeps pytest out).

## 2. Extend the engine (each needs a parity run vs devsim before widening the box)
- [ ] run_field (psi/E at one bias) on the engine: `reverse_sweep(keep_psi=True)` already gives psi.
- [ ] Widen `_ENGINE_BOX` in etna/api/simulation.py (N_A below 1e18, T outside 250-400, forward bias).
- [ ] Batched facades for the other sweeps that rebuild a devsim device per point:
      cce_vs_fluence (needs lifetime/doping damage in the engine), temperature sweep, dark current
      (TAT/SRV not modelled), flash/transient (transient is time-dependent).
- [ ] CCE batch memory: ContinuityBatch holds ~30 (N,B) temporaries (1.2 GB at B=10000); chunk the batch
      if larger sweeps become routine.
- [ ] 2D: no engine path yet (PERF-01 option D). The 2D tests dominate the remaining suite time.
- [ ] Fix devsim cv_sweep 0.5 V lattice (A-01) so the fallback path is also exact.

## 3. Audit issues still open (see .issues/audit-v7/issues/*-INDEX.md)
Highest value next: B-01 (stale TAT field in sweep_iv_vs_temperature, leakage 7.6x low), B-02
(_compute_node_efield garbage on 2D), B-03 (microdosimetric_sweep fails with N_D_junction=None),
D-02/B-06 (DESIGN-1 "18 pA" not produced by the model: measured 141 pA), A-02..A-15.

## 4. Housekeeping
- [ ] Open a PR from fix/audit-v7-top into main when ready (not done; user has not asked).
- [ ] lab/gpu_solver stays as the research record (MLX/numpy backends live only there).
