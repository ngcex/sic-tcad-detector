# lab/gpu_solver — batched 1D SiC diode solver on Apple Silicon

Research spike (2026-10-15, Apple M4 10-core, 16 GB, macOS 24.6, Python 3.13,
numpy 2.5 / numba 0.68 / mlx 0.32.4). Question: can a batched solver on the
Apple GPU make ETNA's sweeps (bias x doping/fluence x temperature) much faster
than the devsim path, at the same accuracy? Everything below is measured, with
the reproducing command; nothing is extrapolated unless marked.

**Answer in one line:** yes, 100-10 000x — but the win comes from batching plus
compiled float64 CPU code (numba), not from the GPU. The MLX/Metal float32 path
works and is validated, but never beats numba at any batch size on this machine.

## What exists

```
gpusolver/physics.py       4H-SiC numbers, mesh, graded doping, alpha generation (ETNA-identical)
gpusolver/poisson_np.py    batched quasi-Fermi Poisson, numpy float64 (reference; base class)
gpusolver/poisson_numba.py same, numba float64, per-device Newton loop in prange   <- primary
gpusolver/poisson_mlx.py   same, MLX/Metal float32, custom Thomas kernel + mx.compile
gpusolver/continuity.py    SG electron/hole continuity + SRH + generation, float64 CPU, Gummel
gpusolver/extract.py       depletion width (ETNA's n/N_D >= 0.5 rule), C, E-field
gpusolver/sweep.py         make_batch / reverse_sweep (W, C, psi vs V) / cce_sweep
gpusolver/ref_devsim.py    devsim reference builds + sweeps (needs `uv sync --extra ref`)
bench/                     problems.py, bench_devsim.py, bench_batched.py, devsim mesh
tests/test_vs_devsim.py    9 tolerance tests vs devsim (float64 1e-10, float32 1e-4)
notes/                     spike results (MLX, torch-MPS), float32 ablation
```

Physics scope: ETNA's 1D p+/n- graded diode (etna/core/device.py mesh anchors,
incomplete Al ionisation, calibrated n_i(T), devsim ohmic contacts). One batch
= B devices sharing one mesh (same geometry), any per-device doping / T / bias.

### The key simplification (measured, not assumed)

Under reverse bias the devsim drift-diffusion solution of this SiC diode is
numerically identical to a single nonlinear Poisson equation with flat
quasi-Fermi levels (phi_n = V_cathode, phi_p = V_anode): max |delta psi| =
1.4e-14 V, W and E to 1e-14 relative at every point 0..-60 V (currents are
~1e-24 A/cm^2, so they carry no electrostatic information). W(V), C-V, E(x),
full-depletion voltage therefore need one unknown per node, not three, and the
batched Newton converges with direct 5 V bias jumps (devsim ramps in 0.5 V).
CCE adds the continuity equations as a float64 Gummel stage (3-4 iterations).

## Measured results

Same mesh (devsim's 327 nodes), same bias points, same device parameters.
Errors are vs the devsim DD result. devsim timings exclude Python import.
Reproduce (run one at a time; other load on the machine distorts timings):

```
cd lab/gpu_solver && uv sync --extra ref
uv run python bench/bench_devsim.py                     # ~4 min, writes results/devsim.json
uv run python bench/bench_batched.py --backend numba    # float64, tol 1e-8 (devsim-grade)
uv run python bench/bench_batched.py --backend numpy
uv run python bench/bench_batched.py --backend mlx      # float32, tol 1e-4
uv run python bench/bench_batched.py --backend numba --tol 1e-4   # same tol as mlx
uv run pytest -q                                        # ~6 s
```

| Problem                                      | devsim               | numba f64 | numpy f64 | MLX f32  | max err vs devsim (numba / mlx) |
| -------------------------------------------- | -------------------- | --------- | --------- | -------- | ------------------------------- |
| P1: 1 device, 121 biases 0..-60 V, W(V)      | 1.44 s (12 ms/solve) | **14 ms** | 270 ms    | 97 ms    | W 9e-15 / 2.8e-6                |
| P2: 100 devices x 10 biases, W(V) (1000 pts) | 102.5 s              | **10 ms** | 123 ms    | 43-57 ms | W 5e-15 / 1.5e-6                |
| P3: CCE, 100 devices x 5 biases (500 pts)    | 98.8 s               | **69 ms** | 191 ms    | 111 ms   | CCE 2.9e-7 / 8.3e-7 (abs)       |

Speed-ups (numba vs devsim): P1 ~100x, P2 ~10 000x, P3 ~1400x. At tol 1e-4
(the float32 setting) numba is 13 / 9 / 96 ms, i.e. the tolerance is not what
separates the backends.

Batch scaling, one 0..-45 V sweep (10 biases), random dopings, wall time:

| B      | numba f64 | numpy f64 | MLX f32    | devsim (extrapolated from 1.02 s/device) |
| ------ | --------- | --------- | ---------- | ---------------------------------------- |
| 1      | 1.3 ms    | 58 ms     | 21 ms      | 1.0 s                                    |
| 10     | 2 ms      | 76 ms     | 29 ms      | 10 s                                     |
| 100    | 10 ms     | 120 ms    | 70-96 ms   | 102 s (measured, P2)                     |
| 1 000  | 109 ms    | 675 ms    | 214-256 ms | ~17 min                                  |
| 10 000 | 1.1 s     | 9.7 s     | 1.3 s      | ~2.8 h                                   |

Numba scales linearly at ~0.1 ms per device per 10-bias sweep (10 cores). MLX
only approaches it at B = 10 000 (1.3 vs 1.1 s) and is float32.

Float32 accuracy of the MLX path (tests/test_vs_devsim.py, default device,
0..-60 V): psi 3.2e-7 V, W 6e-7 relative, E 1e-5 of E_max, CCE 8e-7. It needs
the potential carried as a float32 hi/lo pair: without it the Newton update
stalls between 1e-4 and 1e-3 V_t and tol = 1e-4 never converges
(notes/ablate_f32_split.py).

## Hardware facts verified

- MLX: `float64 is not supported on the GPU` (CPU stream only). torch-MPS:
  same (`MPS framework doesn't support float64`).
- Carrier densities in SiC span 1e-35..1e19 cm^-3: not representable in
  float32 at all. Only the potential (with quasi-Fermi-referenced exponents
  that are <= 0 for physical states) can live on the GPU; the continuity
  equations stay float64 on the CPU (mixed precision).
- Per-call floor: MLX ~0.2-0.5 ms per compiled kernel + sync, torch-MPS
  ~0.3-0.4 ms. A 327-node system is ~30 us of numba work, so the GPU can only
  win when B >~ 1000 systems are launched together — and then numba's 10 cores
  are already there.
- Tridiagonal solves on the GPU (notes/spike_mlx_tridiag.md,
  notes/spike_torch_mps.md): Python-loop Thomas is dispatch-bound (10-45 ms);
  PCR loses to numpy at B >= 1000; the only competitive MLX path is a custom
  Metal kernel running one sequential Thomas per thread on an (N, B) layout
  (4.7x numpy at B=1000 for the solve alone). torch-MPS: compiled PCR is 6-11x
  numpy per Newton iteration at B >= 1000, dense batched solve is hopeless.
  float32 Thomas on a near-pure Laplacian Jacobian (depleted region) is only
  3e-4 accurate per solve; Newton tolerates this because the residual is what
  drives convergence.

## What does NOT work / negative results

- **GPU for 1D is not worth it.** numba float64 beats MLX float32 at every
  batch size measured (1..10 000) at equal or better accuracy. torch-MPS was
  measured on synthetic Newton iterations only (6-11x vs numpy at B >= 1000);
  it was not built into the solver because the MLX kernel path, which is
  faster than torch's PCR per the spikes, already loses to numba.
- numpy float64 vectorised over the batch is 10-30x slower than numba: the
  Python-level Thomas loop (327 steps) costs ~0.75 ms per call regardless of
  B, and at large B the exp()-dominated assembly is memory-bound.
- The MLX Gummel/CCE loop only converges if the cathode current is evaluated
  on the continuity-consistent state (before the Poisson update); the current
  is a ~2e-7 relative difference of two SG terms and float32 potential noise
  swamps it otherwise (first attempt: 121 iterations, not converged).
- Not done: 2D. The data gives no reason to try a GPU 2D path before a
  compiled CPU one (see next steps). Not done: forward-bias I-V, transient,
  temperature-dependent incomplete ionisation of donors (ETNA has none either),
  anything with devsim's 2D/mesa/guard-ring structures.
- Dependence on devsim's mesh: validation uses devsim's node positions
  (bench/mesh_devsim_default.npy); `physics.make_mesh` reproduces the anchors
  but not devsim's exact grading (150 vs 327 nodes), so results on it differ
  from devsim's at the discretisation-error level, not the solver level.

## Recommended next steps

1. Promote the numba path (poisson_numba + continuity + sweep) as ETNA's
   batched 1D engine for C-V, W(V), E(x), CCE-vs-bias/fluence sweeps; keep the
   devsim path as the reference and the tolerance tests as the gate. Fluence
   and temperature enter only through doping/lifetime/n_i arrays, which are
   already per-device batch inputs.
2. Replace devsim's 0.5 V ramp with direct jumps + the secant predictor
   (both backends converge at 5 V steps; 1.5x fewer Newton iterations).
3. Add a batched mesh generator that reproduces devsim's grading exactly, or
   accept ours and re-baseline the calibration once.
4. For forward bias / high injection (FLASH), a coupled Newton (3 unknowns per
   node, block-tridiagonal) in numba is the natural extension; the continuity
   stage here is Gummel and will slow down when the coupling is strong.
5. 2D: structured-mesh Newton in numba (banded LU or multigrid) first; revisit
   the GPU only if a profile shows the linear solve dominating at >= 1e5 nodes.
6. Drop the MLX backend from a promoted version unless B >= 10 000 sweeps
   become routine; it is kept here as a validated, documented experiment.
