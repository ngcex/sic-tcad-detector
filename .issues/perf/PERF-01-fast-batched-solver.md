---
id: PERF-01
title: Make ETNA simulations much faster — batched CPU float64 solver first, GPU only where measured to win
severity: improvement
category: performance
status: open
location:
  - etna/core/{drift_diffusion,poisson,device,cv_analysis,charge_collection,dark_current,temperature_sweep,flash_recombination,transient}.py
  - etna/api/sweep.py (ParametricSweep)
  - lab/gpu_solver/ (Fable research worktree, in progress)
labels: [performance, solver, numba, gpu, apple-silicon, research]
---

## Summary

Sweeps (bias x fluence x temperature x geometry x doping) are the dominant workload and
are slow. Every sweep point builds a fresh devsim device (mesh, models, equilibrium solve,
DD init) and solves on CPU through devsim's global-state, single-threaded, float64 sparse
direct solver. Goal: reduce wall time of typical sweeps by one to two orders of
magnitude without losing accuracy versus the current results.

## Where the time goes (hypotheses to MEASURE first, not assume)

1. Per-point overhead: device creation + equilibrium solve + model setup is repeated for
   every point of cce_vs_fluence, temperature sweeps, flash/transient sweeps.
2. Solver cost: a 1D device has only ~10^2-10^3 nodes; the Newton step is a
   block-tridiagonal solve (microseconds in compiled code), so devsim's generic sparse
   machinery and Python-level model evaluation dominate, not linear algebra.
3. No parallelism: devsim state is process-global and not thread-safe (now guarded by a
   lock, audit C-04), so independent sweep points run strictly serially.
4. Redundant work: continuation ramps restart from equilibrium for each point instead of
   warm-starting from the neighbouring converged point.

## Options, in order of expected return per unit of risk

| # | Option | Expected gain | Risk |
|---|---|---|---|
| A | Cheap wins inside devsim path: warm-start along sweeps, reuse device across bias points, `multiprocessing` with one devsim per process (state is per-process) | 2-8x on a multi-core M4, near zero numerical risk | low |
| B | Batched 1D solver in numpy/numba, float64 CPU: Scharfetter-Gummel + Newton/Gummel with block-tridiagonal (Thomas) solves, `prange` over independent points | 10-100x on sweeps (hypothesis) | medium: must re-implement and validate physics |
| C | Same batched solver on Apple GPU (MLX / torch-MPS / Metal kernels) | extra gain only for very large batches | high: Metal has no float64; needs scaling, log/quasi-Fermi variables, or mixed precision with float64 refinement |
| D | 2D: structured-mesh Newton with GPU or multithreaded CG / multigrid | the only place a big matrix exists, GPU plausible | high |

Key point: a GPU does not help one small 1D system (launch + transfer overhead exceeds the
compute). The win is batching many independent points, and most of that win is already
available on CPU in float64 (option B). Decide GPU (C/D) only on measured data.

## Plan

1. Benchmarks (committed under lab/gpu_solver/bench/): single 1D bias ramp; batched sweep
   of ~1000 points; 2D mesh. Record devsim baseline times and a profile (cProfile / py-spy)
   to confirm where time goes.
2. Option A quick wins in the main package behind identical outputs; measure.
3. Option B as a separate module, validated against devsim before anything replaces it.
4. Spike C/D with MLX, torch-MPS, Accelerate; keep only if measured faster at equal
   accuracy.

## Acceptance criteria

- Every claimed speedup is measured, same problem, same accuracy, with the reproducing
  command; negative results are documented.
- New solver agrees with devsim on equilibrium potential, I-V, C-V, depletion width and
  CCE within stated tolerances (target <=1e-6 relative in float64 CPU, tolerance
  documented separately for any float32 path).
- No change to existing outputs until a new path passes validation; old path stays the
  reference until then.
- No extra dependency in the root pyproject without a measured reason; experiments live
  in lab/ with their own pyproject until promoted.
- Tests: tolerance-based comparison suite vs devsim (fast subset in CI, full set marked
  slow).

## Notes

- Related: audit C-04 (devsim lock) serializes devsim access; option A multiprocessing
  sidesteps it because state is per-process.
- Fable's lab worktree is exploring B/C/D now; its README with measured numbers will feed
  step 4. Numba was not named explicitly in its brief and should be a mandatory baseline.
- Effort: A = S/M, B = L, C/D = L+ (research).
