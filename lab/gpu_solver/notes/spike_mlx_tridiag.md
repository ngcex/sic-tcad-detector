# Spike: MLX/Metal float32 vs numpy float64 for a batched 1D tridiagonal Newton inner loop

Date: 2026-10-15. Hardware: Apple M4, 16 GB unified memory, macOS 24.6 (Darwin).
Software: mlx 0.32.4, numpy 2.5.4 (Accelerate), Python 3.13 via uv.
MLX GPU has no float64 (`float64 is not supported on the GPU`); all MLX numbers below are float32 on the GPU stream.

## Reproduce

```
cd lab/gpu_solver && uv run python notes/spike_mlx_tridiag.py
```

(`--quick` runs a reduced matrix. Full run takes ~25 s.)

Timing protocol: wall time with `time.perf_counter`, median of 5 runs after 2 warmups. For MLX, `mx.eval()` is called on every output
before the clock stops. Inputs are already resident (numpy arrays or MLX arrays); host<->device conversion is NOT included.
`mx.compile` trace/compile time is excluded from the timed runs (first-call cost was 20-70 ms per function).

## Test matrices

- Nonuniform 1D mesh: N points over 11e-4 cm, junction at the middle, spacing h geometrically graded from 3e-8 cm at the junction
  up to (capped at) 1e-5 cm at the edges; grading ratio chosen by bisection so the intervals sum exactly to the half-length
  (for N=327 this gives 163+163 intervals; coarse edge region is uniform at 1e-5).
- Poisson-like tridiagonal: sub[i] = 1/h[i-1], super[i] = 1/h[i], diag[i] = -(1/h[i-1] + 1/h[i]) - d_i,
  with d_i = 10^U(-2, 8) random per row ("spec" regime), Dirichlet rows (diag 1) at both ends, rhs ~ N(0,1).
  1/h ranges 1e5 .. 3.3e7. Dense 2-norm condition number of system 0 (N=327): 1.8e8.
- Reference: numpy float64 Thomas. All float32 variants are compared against it.

## Methods

1. `np64 Thomas (B,N)` / `(N,B)`: numpy float64 Thomas, Python loop over N, each step vectorised over B. (B,N) touches a strided
   column per step, (N,B) a contiguous row.
2. `mx32 Thomas loop`: same algorithm in MLX float32, Python loop over N, layout (N,B); eager and wrapped in `mx.compile`.
3. `mx32 PCR`: parallel cyclic reduction, N padded to the next power of two (327->512, 1000->1024) with identity rows, log2(N) sweeps
   fully vectorised over (B,N) with concatenate-based shifts; eager and `mx.compile`d.
4. `mx32 Metal kernel`: `mx.fast.metal_kernel`, one GPU thread per system, sequential Thomas in float32 over its N unknowns,
   grid=(B,1,1), threadgroup=(min(256,B),1,1), scratch `cp` as a second output. Two layouts: (B,N) as specified (thread-strided
   reads) and transposed (N,B) (coalesced reads across threads).
5. Newton assembly F,J (2 exp + ~10 elementwise ops) on (B,N) arrays, numpy float64 vs MLX float32 eager vs compiled.
6. Full Newton iteration = assembly + Thomas solve + update, numpy float64 (N,B) vs MLX float32 (compiled assembly + Metal (N,B) kernel).

## Results, run 1 (primary): tridiagonal solve wall time (ms)

| N    | B     | np64 Thomas (B,N) | np64 Thomas (N,B) | mx32 Thomas loop eager | mx32 Thomas loop compiled | mx32 PCR eager | mx32 PCR compiled | mx32 Metal kernel (B,N) | mx32 Metal kernel (N,B) |
| ---- | ----- | ----------------- | ----------------- | ---------------------- | ------------------------- | -------------- | ----------------- | ----------------------- | ----------------------- |
| 327  | 1     | 0.829             | 0.936             | 16.4                   | 6.947                     | 1.156          | 3.117             | 1.125                   | 0.344                   |
| 327  | 100   | 0.867             | 0.758             | 9.805                  | 5.380                     | 1.071          | 0.632             | 0.320                   | 0.293                   |
| 327  | 1000  | 5.760             | 2.192             | 12.2                   | 7.203                     | 15.3           | 9.100             | 0.510                   | 0.467                   |
| 327  | 10000 | 67.0              | 12.0              | 11.3                   | 6.783                     | 163.2          | 91.7              | 2.892                   | 2.044                   |
| 1000 | 1     | 2.507             | 2.221             | 38.5                   | 20.6                      | 1.082          | 0.507             | 0.405                   | 0.476                   |
| 1000 | 100   | 4.010             | 2.688             | 35.9                   | 20.4                      | 2.662          | 1.490             | 0.535                   | 0.634                   |
| 1000 | 1000  | 12.7              | 6.741             | 36.0                   | 21.8                      | 29.0           | 17.4              | 0.928                   | 1.053                   |
| 1000 | 10000 | 240.2             | 83.0              | 36.9                   | 23.2                      | 285.8          | 175.1             | 8.047                   | 4.411                   |

## Results, run 2 (repeat, same command, a few minutes later): tridiagonal solve wall time (ms)

| N    | B     | np64 Thomas (B,N) | np64 Thomas (N,B) | mx32 Thomas loop eager | mx32 Thomas loop compiled | mx32 PCR eager | mx32 PCR compiled | mx32 Metal kernel (B,N) | mx32 Metal kernel (N,B) |
| ---- | ----- | ----------------- | ----------------- | ---------------------- | ------------------------- | -------------- | ----------------- | ----------------------- | ----------------------- |
| 327  | 1     | 0.827             | 0.743             | 22.7                   | 6.993                     | 0.922          | 0.500             | 0.453                   | 0.457                   |
| 327  | 100   | 0.955             | 0.839             | 13.5                   | 7.063                     | 2.381          | 1.540             | 0.836                   | 0.900                   |
| 327  | 1000  | 6.446             | 3.986             | 13.5                   | 6.315                     | 13.4           | 9.310             | 0.527                   | 0.643                   |
| 327  | 10000 | 61.4              | 21.0              | 11.2                   | 8.673                     | 144.7          | 91.5              | 2.746                   | 1.698                   |
| 1000 | 1     | 2.790             | 2.225             | 46.1                   | 29.9                      | 1.520          | 0.569             | 0.410                   | 0.399                   |
| 1000 | 100   | 3.886             | 3.316             | 49.4                   | 41.5                      | 2.931          | 1.695             | 0.577                   | 0.631                   |
| 1000 | 1000  | 16.8              | 8.338             | 66.1                   | 29.5                      | 67.0           | 20.5              | 0.940                   | 1.050                   |
| 1000 | 10000 | 264.5             | 60.9              | 38.2                   | 22.9                      | 315.7          | 190.8             | 8.390                   | 4.463                   |

Run-to-run variability: the Metal kernel and the B=10000 cells are stable to ~10-20 %; mid-size cells (B=100..1000) vary by up to
~2x between runs (e.g. np64 (N,B) 327x1000: 2.19 vs 3.99 ms; mx loop eager 1000x1000: 36 vs 66 ms). The N=327, B=1 row of run 1
(Metal (B,N) 1.125, PCR compiled 3.117) is the first config measured and caught one-off JIT/clock-ramp noise; run 2 is representative
(0.45 / 0.50 ms). Conclusions below use ratios that hold in both runs.

Observations on the solve:

- numpy layout: (N,B) (contiguous row per step) is always faster than (B,N); 1.1x at small B, 2.6-5.6x at B=10000. Use (N,B) in numpy.
- numpy float64 has a floor of ~0.75 ms (N=327) / ~2.2 ms (N=1000) from the Python loop over N, independent of B up to B~100.
- The MLX Python-loop Thomas is dispatch-bound: ~10-20 ms (N=327) / 20-45 ms (N=1000) regardless of B. `mx.compile` roughly halves
  it but it is still 5-30x slower than numpy at B<=1000. Only at B=10000 does it beat numpy (N,B), and then it is still 3-5x slower
  than the Metal kernel. Not a viable path.
- PCR: fine at B=1..100 (~0.5-1.5 ms compiled, latency-bound) but it does log2(N)*~20 full-array passes, so it scales as O(B*N*log N)
  memory traffic and loses badly at B>=1000 (9-190 ms). Loses to numpy at B=1000 and B=10000. Not a viable path at scale.
- Metal kernel (one thread per system, sequential Thomas): 0.3-0.5 ms floor (MLX dispatch+sync latency), then ~0.2 us per system at
  N=327 and ~0.45-0.8 us per system at N=1000. The transposed (N,B) coalesced layout is 1.4-1.8x faster than (B,N) at B=10000 and
  equal within noise below that. This is the only MLX path that wins everywhere.

## float32 accuracy (spec regime, d_i = 10^U(-2,8)): max over systems of ||x32 - x64||_inf / ||x64||_inf

| N    | B     | method                  | raw float32 | row-scaled (in f32) | row-scaled (in f64, then cast) |
| ---- | ----- | ----------------------- | ----------- | ------------------- | ------------------------------ |
| 327  | 1     | mx32 Thomas loop eager  | 4.27e-08    | 1.02e-08            | 4.30e-08                       |
| 327  | 1     | mx32 PCR eager          | 7.55e-08    | 4.27e-08            | 1.02e-08                       |
| 327  | 1     | mx32 Metal kernel (B,N) | 4.27e-08    | 1.02e-08            | 7.55e-08                       |
| 327  | 1     | mx32 Metal kernel (N,B) | 4.27e-08    | 1.02e-08            | 7.55e-08                       |
| 327  | 100   | mx32 Thomas loop eager  | 2.87e-07    | 2.87e-07            | 2.08e-07                       |
| 327  | 100   | mx32 PCR eager          | 1.66e-07    | 3.52e-07            | 1.60e-07                       |
| 327  | 100   | mx32 Metal kernel (B,N) | 1.80e-07    | 2.52e-07            | 1.32e-07                       |
| 327  | 100   | mx32 Metal kernel (N,B) | 1.80e-07    | 2.52e-07            | 1.32e-07                       |
| 327  | 1000  | mx32 Thomas loop eager  | 3.42e-07    | 3.97e-07            | 3.18e-07                       |
| 327  | 1000  | mx32 PCR eager          | 4.32e-07    | 3.60e-07            | 5.57e-07                       |
| 327  | 1000  | mx32 Metal kernel (B,N) | 3.43e-07    | 5.75e-07            | 2.90e-07                       |
| 327  | 1000  | mx32 Metal kernel (N,B) | 3.43e-07    | 5.75e-07            | 2.90e-07                       |
| 327  | 10000 | mx32 Thomas loop eager  | 6.46e-07    | 5.83e-07            | 5.65e-07                       |
| 327  | 10000 | mx32 PCR eager          | 9.38e-07    | 5.18e-07            | 6.70e-07                       |
| 327  | 10000 | mx32 Metal kernel (B,N) | 6.92e-07    | 6.57e-07            | 6.69e-07                       |
| 327  | 10000 | mx32 Metal kernel (N,B) | 6.92e-07    | 6.57e-07            | 6.69e-07                       |
| 1000 | 1     | mx32 Thomas loop eager  | 2.32e-08    | 1.19e-07            | 9.10e-08                       |
| 1000 | 1     | mx32 PCR eager          | 9.10e-08    | 9.10e-08            | 9.10e-08                       |
| 1000 | 1     | mx32 Metal kernel (B,N) | 2.32e-08    | 1.19e-07            | 9.10e-08                       |
| 1000 | 1     | mx32 Metal kernel (N,B) | 2.32e-08    | 1.19e-07            | 9.10e-08                       |
| 1000 | 100   | mx32 Thomas loop eager  | 3.97e-07    | 4.01e-07            | 2.42e-07                       |
| 1000 | 100   | mx32 PCR eager          | 3.81e-07    | 3.25e-07            | 2.50e-07                       |
| 1000 | 100   | mx32 Metal kernel (B,N) | 4.23e-07    | 5.23e-07            | 1.51e-07                       |
| 1000 | 100   | mx32 Metal kernel (N,B) | 4.23e-07    | 5.23e-07            | 1.51e-07                       |
| 1000 | 1000  | mx32 Thomas loop eager  | 6.11e-07    | 7.55e-07            | 6.77e-07                       |
| 1000 | 1000  | mx32 PCR eager          | 5.72e-07    | 6.74e-07            | 6.45e-07                       |
| 1000 | 1000  | mx32 Metal kernel (B,N) | 4.13e-07    | 7.04e-07            | 2.67e-07                       |
| 1000 | 1000  | mx32 Metal kernel (N,B) | 4.13e-07    | 7.04e-07            | 2.67e-07                       |
| 1000 | 10000 | mx32 Thomas loop eager  | 8.49e-07    | 9.13e-07            | 5.36e-07                       |
| 1000 | 10000 | mx32 PCR eager          | 9.48e-07    | 9.36e-07            | 6.26e-07                       |
| 1000 | 10000 | mx32 Metal kernel (B,N) | 7.62e-07    | 6.94e-07            | 4.65e-07                       |
| 1000 | 10000 | mx32 Metal kernel (N,B) | 7.62e-07    | 6.94e-07            | 4.65e-07                       |

In the spec regime every float32 method is accurate to <1e-6 relative (worst 9.5e-7), with or without row scaling; row scaling
makes no systematic difference. That is because random d_i up to 1e8 makes most rows strongly diagonally dominant.

## float32 accuracy vs conditioning (N=327, B=100, Metal kernel (B,N); numpy float32 Thomas as cross-check)

| d_i regime                        | cond_2 (system 0, dense) | mx32 Metal raw | mx32 Metal row-scaled | np32 Thomas raw | np32 Thomas row-scaled |
| --------------------------------- | ------------------------ | -------------- | --------------------- | --------------- | ---------------------- |
| laplacian (d_i = 1e-2 everywhere) | 1.6e+09                  | 5.31e-04       | 1.26e-02              | 2.60e-04        | 2.21e-02               |
| spec (d_i = 10^U(-2,8))           | 1.8e+08                  | 2.67e-07       | 2.75e-07              | 2.92e-07        | 2.33e-07               |
| dominant (d_i = 1e8 everywhere)   | 2.3e+08                  | 5.25e-08       | 5.25e-08              | 5.25e-08        | 5.25e-08               |

This is the honest caveat. When the Jacobian is a nearly pure nonuniform Laplacian (depleted region: carrier terms negligible against
1/h ~ 1e5..3e7), float32 Thomas is only accurate to 3e-4 .. 5e-4 relative, and row scaling by the diagonal makes it WORSE
(1e-2 .. 2e-2). numpy float32 reproduces the same numbers, so this is float32 arithmetic, not a kernel bug. The same effect shows in
the full Newton iteration below (psi update discrepancy 2e-2 with s=1e-3, cond 1.6e9; with s=1e5 the discrepancy drops to 3e-6 and with
s=1e8 to 4e-7; measured separately, N=327, B=8).

**Answer to the accuracy question:** float32 solves of these ill-scaled systems ARE accurate to ~1e-5 relative (in fact <1e-6) as long
as rows have an O(1/h)-or-larger diagonal shift d_i (the spec'd regime). They are NOT accurate to 1e-5 when the system degenerates to a
pure Laplacian over the fine mesh (3e-4 .. 5e-4, and 1e-2 if row-scaled). A production float32 Poisson solver would need either
mixed-precision iterative refinement (float32 solve, float64 residual on the CPU stream) or acceptance of ~1e-3 relative accuracy in
depleted regions.

## Newton assembly (F and J, 2 exp + ~10 elementwise ops), wall time (ms), run 1

| N    | B     | np64  | mx32 eager | mx32 compiled | F rel err f32 vs f64 |
| ---- | ----- | ----- | ---------- | ------------- | -------------------- |
| 327  | 1     | 0.009 | 0.578      | 0.200         | 7.9e-08              |
| 327  | 100   | 0.430 | 0.349      | 0.316         | 7.9e-08              |
| 327  | 1000  | 6.190 | 1.029      | 0.662         | 1.2e-07              |
| 327  | 10000 | 57.3  | 5.552      | 2.219         | 1.2e-07              |
| 1000 | 1     | 0.013 | 0.496      | 0.173         | 9.5e-08              |
| 1000 | 100   | 1.445 | 0.585      | 0.527         | 1.1e-07              |
| 1000 | 1000  | 17.5  | 2.517      | 1.380         | 1.2e-07              |
| 1000 | 10000 | 179.1 | 15.2       | 5.966         | 1.4e-07              |

(run 2: np64 6.25 / 69.9 ms and mx32 compiled 0.35 / 2.80 ms at N=327, B=1000 / 10000.)
Assembly is where MLX shines: numpy float64 exp on 3.3M doubles costs ~57 ms, MLX compiled float32 2.2 ms (26x). At B=1 numpy wins by
20-50x because MLX has a ~0.2-0.5 ms dispatch+sync floor. Crossover at B ~ 100.

## One full Newton iteration (assembly + tridiagonal solve + update), N=327, wall time (ms)

| N   | B     | np64 (N,B) Thomas | mx32 Metal (N,B) eager | mx32 Metal (N,B) compiled | speedup (compiled) | rel diff psi (f32 vs f64) |
| --- | ----- | ----------------- | ---------------------- | ------------------------- | ------------------ | ------------------------- |
| 327 | 1000  | 8.920             | 1.590                  | 1.178                     | 7.6x               | 2.1e-02                   |
| 327 | 10000 | 77.3              | 8.090                  | 4.278                     | 18x                | 1.9e-02                   |

(run 2: np64 11.1 / 101.6 ms, mx32 compiled 1.052 / 7.598 ms -> 10.6x / 13.4x.)
`mx.compile` of a function containing the `mx.fast.metal_kernel` call works. The psi discrepancy of 2e-2 is the Laplacian-conditioning
effect described above (s=1e-3 makes the test Jacobian cond 1.6e9); it is not a speed artefact.

## What worked / what did not

- `mx.fast.metal_kernel` worked on the first attempt with the documented signature (`inputs=`, `template=`, `grid=`, `threadgroup=`,
  `output_shapes=`, `output_dtypes=`); scalar B, N passed as 1-element int32 arrays. One bug was found and fixed in the transposed
  kernel's back-substitution index (`x[k+1]` -> `x[k+B]`); the accuracy table caught it (error 1.1 -> 3e-7).
- `mx.compile` of the 327-step Python-loop Thomas traced in ~70 ms and ran, but is still dispatch-bound.
- PCR is correct (same accuracy as Thomas) but not competitive beyond B~100.

## Conclusion

1. The only MLX path that beats numpy float64 is the custom Metal kernel (one thread per system, sequential Thomas, transposed (N,B)
   layout): solve 0.47 vs 2.2 ms at B=1000 (4.7x) and 2.0 vs 12 ms at B=10000 (5.9x) for N=327; 1.05 vs 6.7 ms (6.4x) and 4.4 vs 83 ms
   (19x) for N=1000. The MLX Python-loop Thomas and PCR both lose to numpy at B>=1000 and should be dropped.
2. The crossover B for the Metal kernel is at or below B=1 for the solve alone (numpy's Python loop over N costs 0.75 ms regardless of B;
   the kernel's dispatch floor is 0.3-0.5 ms) and B ~ 100 for the assembly; for a full Newton iteration MLX wins by 7.6-10.6x at B=1000
   and 13-18x at B=10000 (1.2 vs 8.9 ms and 4.3 vs 77 ms, N=327).
3. Float32 accuracy on the spec'd ill-scaled systems (d_i spanning 1e-2..1e8) is <1e-6 relative for every method, with or without row
   scaling, so "~1e-5 relative" is met there.
4. Float32 is NOT accurate to 1e-5 when the system is a near-pure Laplacian on the fine mesh (cond ~1e9): 3e-4..5e-4 raw, and row scaling
   makes it worse (1e-2). A real Poisson Newton step in a depleted region is in this regime (measured psi discrepancy 2e-2 at s=1e-3).
5. Recommendation: MLX float32 is worth it for B>=1000 only with mixed precision (float32 Metal Thomas as the inner solve + float64 residual
   / iterative refinement, e.g. on the mx.cpu stream or numpy), not as a drop-in replacement; below B~100 stay with numpy float64 (N,B).
