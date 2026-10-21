# fv3d: 3D Poisson on tensor grids (issue #4 follow-up)

Package `etna/core/fv3d/`, facade `etna.run_field_3d`, docs `docs/USAGE.md` section 11.2,
tests `tests/test_fv3d.py`. Benches in `lab/fv3d/bench/`. All numbers below were measured
on the shared 10-core M4 (16 GB) under load averages of 2-40 (other agents running), so
wall times are upper bounds with +-2x noise; iteration counts and accuracies are not
load-dependent.

## 1. What is solved

Nonlinear Poisson with Boltzmann carriers and FLAT quasi-Fermi levels (phi_n = the
n-contact voltage, phi_p = the p-contact voltage): the reverse-bias dark state of 4H-SiC
(n_i ~ 1e-9 cm^-3), exactly the dark-ramp model of `fv2d.Solver2D.set_bias`. Electrodes
(n+ / p+ columns, trench walls, implanted planes) are ohmic Dirichlet bodies at
V + -V_t ln(N / n_i). Box discretisation (devsim / fv2d) on tensor node lines. Mirror
(zero-flux) faces: the quarter cell of a square / rectangular lattice IS the periodic
lattice (`test_quarter_cell_equals_the_mirrored_lattice`: quarter vs explicit half cell
agree to < 1e-7 V). Not representable: the diagonal (eighth-cell) mirror and true
periodic wrap for lattices without two mirror planes.

Also: `carriers=False` = fully-depleted charge, one linear solve (field above full
depletion); Ramo weighting potentials (Laplace, readout at 1, other contacts 0);
`ramo_signal` (drift line of one charge in E with Caughey-Thomas mu(E), induced charge
and current from the weighting potential).

## 2. Solver

- Newton on psi, globally convergent: the residual is the gradient of a strictly convex
  energy, Armijo backtracking on it; inexact Newton (Eisenstat-Walker forcing).
  Continuation guess (extrapolation / quasi-Fermi shift, clipped to the electrode
  potential range) and bias bisection on failure. A direct jump 0 -> 60 V converges
  (test).
- Linear solves: CG preconditioned by ONE geometric multigrid V-cycle (matrix-free,
  numba prange): every other line per axis, rediscretised box operator on the coarse
  lines, Newton diagonal restricted with P^T, trilinear P masked on Dirichlet nodes,
  red-black Gauss-Seidel (2 pre / 2 post, colours reversed after the coarse correction:
  symmetric preconditioner, checked to 1e-16 in the test), SuperLU at <= 3000 nodes.
- Dirichlet elimination keeps the operator SPD (fv2d's identity rows would make it
  nonsymmetric).
- Cut cells: electrodes given as LEVEL SETS place the Dirichlet value on the true
  surface (an edge cut at fraction theta gets coupling g / theta; after elimination only
  the free node's diagonal and RHS change, so the matrix stays symmetric). Staircase
  circles converge at order 0.5-0.9 (coax max error 0.49 V of 100 V at h = 0.125 um);
  cut cells at order 1.8-1.9 (2.5e-3 V at h = 0.125 um).
- Field: second-order differences, normal component zeroed on the box faces (the
  discrete natural boundary). With np.gradient's one-sided difference at a mirror face
  the symmetric NULL point between two anode columns showed |E| = 400-1200 V/cm
  (decreasing with h): an O(h) artefact (see section 6).

Multigrid robustness (iterations to 1e-10, `scratch` tests and the benches): uniform grids
9-13 iterations independent of N (32^3 .. 200^3). Graded column-cell grids (aspect
ratios up to ~20 between the 0.1-0.25 um column lines and 2 um bulk lines): 20-65
iterations (23 / 43 / 65 at h_fine 0.5 / 0.2 / 0.1 um, 33k-237k nodes). A Galerkin
(RAP) hierarchy with symmetric GS was not better (21 / 36 / 53, slower per iteration):
the loss is the point smoother on anisotropic cells, not the coarse operator. Extreme
grading (geometric ratio 1.3 over a whole axis, aspect ratio > 1000) stalls PCG: keep
the line grading moderate (the builders cap spacing growth at 1.15 per line).

## 3. Validation (`lab/fv3d/bench/validate.py`, 6-16 s)

| check                                                                                                             | result                                                                                                                                                                                                                                                                                                                       |
| ----------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| z-extruded 2D column cell (`axisym_validity.square_cell`, 109x109x3) vs fv2d Solver2D on the same lines, 10-100 V | max                                                                                                                                                                                                                                                                                                                          | dpsi | 1.8e-15 - 1.4e-14 V; depleted 0.99081 at 75 V (fv2d bench: 75 V first sweep point >= 99 %) |
| over-depleted slab (carriers=False, 70 graded lines) vs the exact parabola                                        | 1.6e-12 V                                                                                                                                                                                                                                                                                                                    |
| slab with Boltzmann carriers vs abrupt one-sided junction W = sqrt(2 eps (Vbi + V - 2kT/q) / qN) (0-20 V)         | W(n = N/2) is 0.08-0.19 um (0.2-0.45 L_D) shorter: the depletion-approximation edge vs the n = N/2 point, shrinking with bias                                                                                                                                                                                                |
| coax r 5 / 20 um, 100 V, carriers=False vs analytic                                                               | cut cell: max error 0.110 / 0.032 / 0.0092 / 0.0025 V at h = 1 / 0.5 / 0.25 / 0.125 um (order 1.77, 1.81, 1.89); staircase (cut_cell=False) 2.09 / 1.44 / 0.89 / 0.49 V (order 0.5-0.9)                                                                                                                                      |
| coax with carriers vs fv2d (r, z) (Boltzmann, 5-50 V)                                                             | 0.06-0.11 V, NOT convergent: the nm-scale carrier spill-over layer at ohmic electrode surfaces (L_D(1e19) ~ 1 nm) is resolved by neither grid; fv2d itself moves 0.045 V between 0.05 and 0.01 um surface lines. Irrelevant at the depletion-voltage level (0.1 V of 50-100 V); the carrier-free coax converges at 2nd order |
| manufactured Laplace on graded lines (mirror sides)                                                               | max nodal error order 1.93, 1.95, 2.01                                                                                                                                                                                                                                                                                       |

devsim 3D was not run (not cheap: devsim's 3D needs an external tetrahedral mesh and a
direct solver, the issue #4 extrapolation is 31 GB per factorisation at 0.9 M nodes);
the references are fv2d parity (exact), analytic solutions and Richardson.

## 4. Grid convergence of the column cell (`lab/fv3d/bench/convergence.py`)

Lines scaled by s (h_fine 0.25 s um, h_bulk 2 s um); the line families are not nested
(graded_lines with scaled spacing), so the sequences are partly pre-asymptotic: the
error estimates below are the Richardson estimate where the observed order is in
(0.3, 4), else |q(s_2) - q(s_3)|.

A, through columns (z-invariant, the issue #4 case), pitch 50 um, r 5 um:

| s    | nodes | V99 (n < N/2 nodes) | V99 (charge criterion) |
| ---- | ----- | ------------------- | ---------------------- |
| 2    | 5476  | 74.796              | 75.403                 |
| 1.41 | 13520 | 74.620              | 75.069                 |
| 1    | 28672 | 74.569              | 74.860                 |
| 0.71 | 61952 | 74.411              | 74.724                 |

V99 = 74.5 +- 0.25 V (Richardson 74.47-74.55). The 2D fv2d bench's "75 V" is the first
sweep point past 99 %, consistent.

B, single-sided 30 um deep columns + planar p+ back contact, 50 um thick (true 3D):

| s    | nodes  | V99    | V99 (charge) |
| ---- | ------ | ------ | ------------ |
| 2    | 58867  | 62.108 | 62.241       |
| 1.41 | 135200 | 62.083 | 62.232       |
| 1    | 237568 | 61.958 | 62.067       |

V99 = 62.0 +- 0.15 V (pre-asymptotic: the estimate is the last difference). The bulk
fraction with |E| < 1 kV/cm at 150 V is 0.03-0.2 % in both cases (+-0.05 % absolute,
not asymptotic: a small region counted box by box).

## 5. Scaling and tractability (`lab/fv3d/bench/scaling.py --big`)

n^3 box Laplacian with a Dirichlet column and back plane (c = 0, the hardest case),
PCG to 1e-10; peak = tracemalloc of numpy allocations (setup + solve):

| N             | MG levels | PCG it | PCG s | peak MB | B/node | SuperLU factor s / fill MB |
| ------------- | --------- | ------ | ----- | ------- | ------ | -------------------------- |
| 32^3 = 32768  | 3         | 9      | 0.26  | 3.0     | 92     | 3.8 / 169                  |
| 40^3 = 64000  | 3         | 11     | 0.46  | 5.8     | 91     | 14.6 / 484                 |
| 64^3 = 262144 | 4         | 10     | 0.44  | 24      | 91     | -                          |
| 96^3 = 884736 | 4         | 10     | 1.11  | 80      | 91     | -                          |
| 128^3 = 2.1e6 | 5         | 10     | 1.39  | 190     | 91     | -                          |
| 160^3 = 4.1e6 | 5         | 13     | 4.13  | 371     | 91     | -                          |
| 200^3 = 8.0e6 | 5         | 13     | 9.19  | 724     | 91     | -                          |

MG and SuperLU solutions agree to 1e-11. Memory is exactly linear (91 B/node for a
linear solve); time ~ N^1 for N >= 2.5e5 (iterations bounded). At the issue #4 named
point (0.9 M nodes) one linear solve is ~1.1 s and 80 MB, against the extrapolated
780 s / 31 GB per SuperLU factorisation.

Full nonlinear solve of the true-3D column cell (30 um columns + p+ back, 50 um) at
870 000 nodes (87 x 100 x 100): setup + equilibrium 10.8 s, 50 V 14.9 s (22 Newton,
129 CG), 100 V 10.5 s (11 Newton, 105 CG); peak 0.24 GB = 271 B/node.

Largest tractable cell: memory is the limit only far out. At 271 B/node an 8 GB budget
holds ~3e7 nodes (e.g. 300 x 300 x 330 lines); one bias then costs ~100-130 CG
iterations at ~2.6 s each (8 M nodes: 0.7 s per iteration measured, linear in N) =
roughly 5 min per bias on this machine, a full depletion scan ~1 h. Practical sweet
spot: 1-5 M nodes, 10-60 s per bias.

## 6. Physics findings

1. **The "low-field saddles" of the square column cell are null points.** The edge
   midpoint between two anode columns (x = a, y = 0 of the quarter cell) lies on two
   mirror planes of the lattice, so Ex = Ey = 0 there exactly; with through columns
   Ez = 0 too: |E| = 0, at every bias. The 90-360 V/cm of the issue #4 2D bench were
   the one-sided np.gradient at the mirror face (O(h); 3D reproduces 1185 / 790 / 594 /
   407 V/cm at s = 2 / 1.41 / 1 / 0.71 with that estimator, 0 with the face-normal
   component set to 0). What matters for charge collection is the SIZE of the low-field
   region around the null point (and that carriers generated exactly on the symmetry
   lines drift slowly): the bulk fraction below 1 kV/cm at 150 V is ~0.1 %.
2. See section 7 for the depletion maps (filled from `design_maps.py`).

## 7. Design maps (`lab/fv3d/bench/design_maps.py`)

(filled below)

## 8. GPU (`uv run --with mlx python lab/fv3d/bench/gpu_smoother.py`)

MLX 0.32.4 on the M4 GPU (Metal) works when added ad hoc (`--with mlx`, not a
dependency). Red-black GS sweep in float32 (masked full-array updates) vs the numba f64
prange kernel:

| grid  | numba f64 sweep | MLX f32 sweep | ratio                           |
| ----- | --------------- | ------------- | ------------------------------- |
| 64^3  | 25 ms           | 3.2 ms        | 7.9x (CPU timing load-inflated) |
| 128^3 | 24.8 ms         | 25.7 ms       | 0.97x                           |
| 192^3 | 23.1 ms         | 47.9 ms       | 0.48x                           |
| 256^3 | 88.8 ms         | 90.9 ms       | 0.98x                           |

f32 vs f64 single-sweep difference 1e-7 relative. Used as the finest-level smoother
inside the f64 PCG: SAME iteration count (10 at 128^3 and 192^3), true residual 4e-11 /
9e-11, solution identical to 8e-15 relative: the f32 smoother costs no accuracy because
the outer Krylov iteration is f64. But the end-to-end PCG is 3-4x SLOWER (10.9 vs 2.6 s,
27 vs 9.4 s) from the per-V-cycle host/device conversions. The stencil is memory-bound
and the M4 CPU and GPU share the same unified memory bandwidth; the best a fully
GPU-resident f32 multigrid could gain is ~1-2x on the smoother. Result: CPU-only (numba
f64); the GPU is not worth a dependency here (consistent with the 1D GPU spike).
