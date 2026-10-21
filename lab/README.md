# lab/

Research spikes. Each directory keeps its own `notes/LOG.md` (commands and numbers).
Nothing in `etna/`, `app/` or `tests/` imports code from `lab/`, but `lab/` is NOT
disposable: `tests/test_measured_cv_preset.py` reads the digitised C-V data in
`lab/cv_gap/notes/` (as does `lab/buffer_profile`), and comments in `etna/` and `tests/`
cite the `lab/*/notes/LOG.md` files as the provenance of calibrated defaults.

## Archived: `lab/gpu_solver` (tag `archive/gpu-solver-lab`)

Removed from the tree on 2026-10-18; recover with
`git checkout archive/gpu-solver-lab -- lab/gpu_solver` (or `git show archive/gpu-solver-lab:lab/gpu_solver/README.md`).

Batched 1D SiC diode solver spike (2026-10-15, Apple M4): can a batched solver on the
Apple GPU make ETNA's bias x doping x temperature sweeps much faster than devsim at the
same accuracy? Yes, 100-10 000x (1 device x 121 biases 1.44 s -> 14 ms; 100 devices x
10 biases 102 s -> 10 ms; CCE 100 x 5 99 s -> 69 ms; errors vs devsim 1e-14 in W,
3e-7 in CCE), but the win came from batching plus compiled float64 CPU code (numba),
not from the GPU: the MLX/Metal float32 path worked and was validated (1e-6 in W) yet
never beat numba at any batch size (B = 1..10 000). Closed negative result for the GPU;
the numba path was promoted to `etna/core/batched/` (run_cv / run_cce and the batch sweeps).

## `lab/solver2d`

2D numba engine spike (2026-10-16), promoted to `etna/core/engine2d/`; production parity
against devsim lives in `tests/test_engine2d*.py`. See `lab/solver2d/README.md` for the
devsim reference dumps (`notes/*.npz`), which are regenerated on demand, not committed.

## `lab/geometry2d`

The general 2D engine `etna/core/fv2d` (issue #4, 2026-10-21) handles laterally
non-uniform and axisymmetric devices. These benches cover:

- parity with devsim on the notebook-19 structures;
- mesh convergence and adaptive refinement;
- the choice of linear solver;
- the validity of the axisymmetric proxy;
- whether 3D is tractable.

`tests/test_fv2d.py` reads `notes/devsim_baselines.json` and `notes/devsim_psi.npz`;
regenerate them with `bench/devsim_baselines.py`.
