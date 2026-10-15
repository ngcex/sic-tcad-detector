# Spike: PyTorch MPS (float32) vs numpy float64 for a batched 1D tridiagonal Newton inner loop

Machine: Apple M4, 16 GB, macOS 15.7.1 (Darwin 24.6), Python 3.13.5, torch 2.14.1, numpy 2.5.4, torch CPU threads = 4.
All numbers below are measured (median of 5 after 2 warmups; `torch.mps.synchronize()` before stopping the clock on MPS).

## Reproduce

```
uv venv /private/tmp/claude-501/-Users-ngcex-projects-physics-petringa/df2eeea8-54e5-4d29-83f3-c7a1539b1919/scratchpad/torch_env --python 3.13
uv pip install --python /private/tmp/claude-501/-Users-ngcex-projects-physics-petringa/df2eeea8-54e5-4d29-83f3-c7a1539b1919/scratchpad/torch_env/bin/python torch numpy
/private/tmp/claude-501/-Users-ngcex-projects-physics-petringa/df2eeea8-54e5-4d29-83f3-c7a1539b1919/scratchpad/torch_env/bin/python \
  lab/gpu_solver/notes/spike_torch_mps.py --json spike_results.json
```

Install took ~3 s (wheels cached). Full run takes a few minutes. Test problem: N+1 spacings geometric 3e-8 .. 1e-5 cm,
sub/super = 1/h, diag = -(1/h_l+1/h_r) - d_i, d_i log-uniform [1e-2, 1e8], random normal rhs.

## 1. MPS availability / float64

`torch.backends.mps.is_available()` = True (built = True). float64 on MPS fails:

```
TypeError: Cannot convert a MPS Tensor to float64 dtype as the MPS framework doesn't support float64. Please use float32 instead.
```

`torch.linalg.solve` (batched, float32) and `torch.compile` (with MPS tensors) both work.

## 2-5. Tridiagonal solve, wall time in ms (build+solve for dense)

`np64` = numpy float64. `np64 PCR` = extra reference. `MPS` columns are float32. `PCR-c` = torch.compile'd PCR.
`dense` = diag_embed build + torch.linalg.solve. `*` none. skip = skipped by the memory rule (B=10000, or matrix >= 4 GB; an
earlier attempt at dense B=1000,N=1000 (4.0 GB matrix) hit an MPS out-of-memory error, so the rule was made strict).

| N    | B     | np64 Thomas | np64 PCR | MPS Thomas loop | MPS PCR eager | MPS PCR-c | MPS dense | CPU f64 dense |
| ---- | ----- | ----------- | -------- | --------------- | ------------- | --------- | --------- | ------------- |
| 327  | 1     | 0.903       | 0.317    | 7.16            | 1.88          | 0.273     | 2.03      | 0.456         |
| 327  | 100   | 0.871       | 2.05     | 9.14            | 1.37          | 0.781     | 21.9      | 48.1          |
| 327  | 1000  | 3.91        | 30.1     | 11.4            | 6.52          | 1.08      | 212       | 465           |
| 327  | 10000 | 46.9        | 406      | 18.1            | 63.1          | 9.36      | skip      | skip          |
| 1000 | 1     | 3.15        | 0.467    | 25.4            | 1.5           | 0.327     | 5.54      | 6.13          |
| 1000 | 100   | 3.1         | 8.49     | 28.7            | 2.07          | 1.38      | 301       | 624           |
| 1000 | 1000  | 13.8        | 140      | 36.2            | 21.2          | 3.54      | skip      | skip          |
| 1000 | 10000 | 289         | 1740     | 57.9            | 246           | 39.6      | skip      | skip          |

Speedup of best MPS method over numpy f64 Thomas (best MPS = compiled PCR except where noted):

| N    | B=1  | B=100 | B=1000 | B=10000 |
| ---- | ---- | ----- | ------ | ------- |
| 327  | 3.3x | 1.1x  | 3.6x   | 5.0x    |
| 1000 | 9.6x | 2.2x  | 3.9x   | 7.3x    |

Eager (uncompiled) MPS PCR loses to numpy Thomas at B=1000 and B=10000 for N=327 (6.52 vs 3.91, 63.1 vs 46.9). The MPS Thomas loop
beats numpy only at B=10000 (18.1 vs 46.9 at N=327; 57.9 vs 289 at N=1000). Note numpy PCR at B=1 (0.317 ms, N=327) is
faster than numpy Thomas (0.903 ms), so the B=1 "3.3x" is against the weaker numpy baseline (vs numpy PCR it is 1.2x).
Dense batched solve is far worse than everything else on both MPS and CPU.

torch.compile compile time (first call, excluded from the timings above): up to 2.19 s; some calls ~0 s (cache hits).

## Accuracy (float32 vs numpy f64 Thomas; max over systems of ||x32-x64||_inf/||x64||_inf)

Worst case over all (B, N) configs run:

| method                                         | raw      | row-scaled by diagonal |
| ---------------------------------------------- | -------- | ---------------------- |
| numpy f32 Thomas (same arithmetic as MPS loop) | 3.00e-06 | 2.05e-06               |
| MPS f32 Thomas loop                            | 3.00e-06 | 2.05e-06               |
| MPS f32 PCR                                    | 4.03e-06 | 2.25e-06               |
| MPS f32 dense                                  | 1.86e-06 | 1.19e-06               |
| numpy f64 PCR (check)                          | 5.28e-15 | 5.00e-15               |

Error grows mildly with B (more systems, worse worst case): e.g. N=327 PCR raw 5.3e-07 (B=1) to 3.4e-06 (B=10000).
Row scaling gives a small gain (roughly 1.3x to 2x), not an order of magnitude. The test matrices are diagonally dominant by construction
(|diag| = 1/h_l + 1/h_r + d_i >= sum of off-diagonals), which is benign for Thomas; real Jacobians near punch-through or deep reverse bias were NOT tested.

## 6. Newton assembly (elementwise F and J on (B, N)), ms

| N    | B     | numpy f64 | MPS f32 eager | MPS f32 torch.compile | torch CPU f64 eager |
| ---- | ----- | --------- | ------------- | --------------------- | ------------------- |
| 327  | 1     | 0.00896   | 0.459         | 0.386                 | 0.017               |
| 327  | 100   | 0.459     | 0.613         | 0.414                 | 0.264               |
| 327  | 1000  | 6.43      | 1.85          | 0.518                 | 3.22                |
| 327  | 10000 | 67.7      | 7.75          | 0.885                 | 28.2                |
| 1000 | 1     | 0.0135    | 0.317         | 0.383                 | 0.020               |
| 1000 | 100   | 1.66      | 0.846         | 0.422                 | 0.926               |
| 1000 | 1000  | 20.4      | 4.85          | 2.47                  | 9.88                |
| 1000 | 10000 | 283       | 21.9          | 5.95                  | 88.5                |

torch.compile works on MPS (it ran and gave speedups). Below B~100 numpy wins (MPS has ~0.3-0.4 ms floor per call).

## 7. Full Newton iteration (assembly + tridiagonal solve + update), N=327, ms

| B     | numpy f64 (assembly + Thomas) | MPS f32 PCR eager | MPS f32 PCR compiled | MPS f32 dense |
| ----- | ----------------------------- | ----------------- | -------------------- | ------------- |
| 1000  | 13.9                          | 8.55              | 2.25                 | 239           |
| 10000 | 164                           | 83.3              | 14.4                 | skip          |

Speedup of compiled MPS over numpy: 6.2x at B=1000, 11.4x at B=10000 (eager: 1.6x and 2.0x).
One-step accuracy of the MPS step vs numpy f64 step (same psi, relative inf-norm of dx): 1.77e-05 (B=1000), 2.18e-05 (B=10000).
That test uses random psi in [-20, 20] and random doping, not a converged device state.

## Conclusion

1. MPS has no float64; everything on the GPU is float32 (linear-solve error ~1e-6 relative, one-step dx error ~2e-5 on this synthetic test).
2. Best path: fully vectorized parallel cyclic reduction under torch.compile (not eager, not dense, not a Python Thomas loop).
3. Versus numpy f64 (Thomas, single process) at N=327: crossover at about B~100 (0.78 vs 0.87 ms); 3.6x at B=1000, 5.0x at B=10000 for the solve alone.
4. Full Newton iteration at N=327: 6.2x (B=1000) and 11.4x (B=10000) faster compiled-MPS than numpy; below B~100 numpy is as fast or faster.
5. Caveats: float32 convergence of the real nonlinear problem (exp(psi) range, ill-conditioned Jacobians) is untested, so a speedup does not by itself justify float32 Newton; numpy baseline is single-threaded plain Python-loop Thomas.
