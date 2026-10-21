# Fast single-ion transient for the CCE(LET) tables (2026-10-21)

Goal: replace devsim's transient in `build_cce_let_table` (22+ min table rebuilds, the
suite's longest test unit) with a numba transient on the same device, with parity
against devsim, adaptive time stepping and error estimates. Code:
`etna/core/track/planar.py` (engine), `etna/core/single_particle_fast.py` (table,
`build_cce_let_table(engine="fast")`). Docs: docs/USAGE.md 19.4. Issues #4 (fv2d), #9
(track engine).

All wall times below were taken on a shared 10-core M4 under a load average of 30-55
(about a dozen agents); they are upper bounds, the ratios are what matter. The earliest
runs (load ~3) are marked.

## 1. Engine

- Mesh: a mesh-only devsim build of `create_sic_2d_device` (same keywords), imported by
  `etna.core.fv2d.mesh_from_devsim` (node volumes, EdgeCouple / EdgeLength, doping) plus
  the node mobilities `mu_n_node` / `mu_p_node`; their arithmetic edge average equals
  devsim's edge `mu_n` / `mu_p` to < 1e-12 (asserted at import). 0.1-0.4 s.
- Dark state: fv2d (equilibrium, flat-quasi-Fermi ramp, G = 0 Gummel check). Dark
  potential vs devsim's at 10 V, hw 2: max 4.0e-11 V.
- Transient: `TrackSystem(geometry=...)` (the track engine's coupled Newton, numba
  `assemble_excess`, SuperLU) on the 5-point neighbour lists, excess unknowns about the
  dark state (dark contact current exactly 0). Physics = devsim's 2D DD (static Poisson,
  SG, edge-averaged local mobility, midgap SRH; no Auger / carrier scattering: C_n =
  C_p = 0). Subclass `PlanarSystem` only swaps the SuperLU ordering (MMD on A^T + A: at
  hw 50, 53k unknowns, fill 5.1M vs 6.3M with COLAMD, solve 10.9 vs 19.0 ms, factor 0.35
  vs 0.44 s; hw 150: 28M vs 35M, 44 vs 64 ms).
- Contacts: conduction I_c = q sum_c (Phi_p - Phi_n) (unconstrained excess fluxes of the
  contact nodes), contact charge Q_c = -q sum_c rpsi (devsim's get_contact_charge: the
  displacement series agree to 1e-7 of the peak, below), displacement dQ_c/dt. The
  anode + cathode terminal currents cancel to < 1e-9 of Q_collected per step (test).
- `transient_devsim_grid`: devsim's scheme (1 ps BDF1 injection with G = gen / dt,
  dt = max(1e-13, min(1e-9, 0.1 t)), stop when |dQ| < 1e-7 Q_run for 3 steps, CCE =
  sum |I| dt with I_dark = 0).
- `transient_adaptive`: TR-BDF2 with the charge-weighted LTE control of
  `etna.core.track.solver` (rtol 1e-4 of N_gen per step), instantaneous deposition,
  stop at stored excess < 1e-9 N_gen; step-size hysteresis (keep h unless the controller
  asks for >= 1.5x or a cut): Newton then reuses the exact factorisation, ~3 solves per
  stage instead of ~10 (hw 2: 6.5 -> 1.7 s, load ~3). CCE = |sum I dt| (signed; the
  anode current changes sign briefly: sum |I| dt is 1e-7 higher), equal to the pair
  balance 1 - (recombined + wrong contact + stored) / N_gen to <= 1e-9; ledger <= 1e-9.

## 2. Parity with devsim (bench/parity.py; references bench/devsim_reference.py -> data/)

notes/parity.out (hw 150 block: first run; hw 2 / hw 50: notes/parity_hw2_hw50.out).
Same mesh, generation (Q_generated identical to 1e-16), same time grid (identical step
times; the stop may differ by one step).

| case (50 V unless noted) | series dI / peak | disp. / peak | CCE twin - devsim | twin - devsim sign-aware | 2 I_dark t_end / Q_gen | t95 twin / devsim - 1 |
| ------------------------ | ---------------- | ------------ | ----------------- | ------------------------ | ---------------------- | --------------------- |
| hw 2, LET 500            | 7.9e-11          | 7.5e-10      | +4.6e-10          | -2.2e-10                 | 6.7e-10                | 0                     |
| hw 2, LET 15.8           | 3.6e-9           | 2.3e-8       | +1.1e-8           | -1.0e-8                  | 2.1e-8                 | 0                     |
| hw 2, LET 5              | 1.1e-8           | 1.2e-7       | +4.8e-8           | -2.0e-8                  | 6.7e-8                 | 0                     |
| hw 2, LET 0.5            | 5.1e-8           | 6.5e-7       | +3.6e-7           | -2.6e-7                  | 6.7e-7                 | -2e-6                 |
| hw 2, 10 V, LET 10       | 1.2e-8           | 1.4e-7       | +6.2e-8           | +1.9e-8                  | 4.3e-8                 | +2e-7                 |
| hw 150, LET 0.5          | 4.0e-6           | 3.6e-5       | +1.72e-5          | -1.5e-6                  | 5.3e-5                 | -6e-6                 |
| hw 50, LET 500           | 7.2e-10          | 3.5e-9       | +1.1e-8           | -4.8e-9                  | 1.6e-8                 | 0                     |
| hw 50, LET 0.5           | 6.6e-7           | 8.7e-6       | +8.0e-6           | -1.1e-6                  | 1.7e-5                 | -5e-6                 |

Reading: the twin IS devsim's discrete solution (LET 500: 5e-10 in CCE, 8e-11 in the
series). The low-LET differences are devsim's: (a) the estimator sum (|I| - |I_dark|) dt
with I_dark = round-off of the 1e18-1e19 contact layers (8.0e-15, 2.0e-13, 6.0e-13 A/cm
at hw 2, 50, 150: it scales with the width) and of the opposite sign to the signal
current, an offset ~2 I_dark t_end / Q_gen (column 6; the twin has I_dark = 0);
(b) round-off noise of the same origin around the small LET-0.5 signal (the sign-aware
devsim sum still misses by 1.5e-6 at hw 150).

## 3. Time: devsim's grid is 8e-6 low in CCE (bench/time_convergence.py, notes/time.out)

hw 2, LET 15.8, 50 V:

| run                            | CCE          | t_collection (ns) | steps |
| ------------------------------ | ------------ | ----------------- | ----- |
| devsim grid (10 %, 1 ps), twin | 0.9999743504 | 0.116788          | 94    |
| grid 5 %, 0.5 ps               | 0.9999777355 | 0.113172          | 186   |
| grid 2.5 %, 0.25 ps            | 0.9999796871 | 0.111417          | 368   |
| grid 1.25 %, 0.125 ps          | 0.9999808027 | 0.110576          | 726   |
| adaptive rtol 1e-3             | 0.9999816950 | 0.109821          | 123   |
| adaptive rtol 1e-4 (default)   | 0.9999823233 | 0.109779          | 234   |
| adaptive rtol 1e-5             | 0.9999821875 | 0.109772          | 500   |
| adaptive rtol 1e-6             | 0.9999821642 | 0.109769          | 1075  |

(The grid rows: devsim's estimator, stop rtol 1e-9.)
The grid increments 3.4e-6, 2.0e-6, 1.1e-6 halve (first order) towards the adaptive
value; the default rtol 1e-4 is within 1.6e-7 (CCE) and 1e-4 (t_collection, relative) of
rtol 1e-6. devsim's grid: CCE -8.15e-6 (the same at every LET 0.5-500 and at hw 2 / 50 /
150), t_collection +6.4 % (0.1168 vs 0.1098 ns); at 10 V -2.9e-5 and +7 % (0.624 vs 0.583 ns).

Origin (bench/time_error_origin.py, notes/time_error_origin.out): for a LINEAR transient
the backward-Euler collected charge is grid-independent (A sum_k h_k x_k = x_end - x_0),
so the error needs a nonlinearity. With the SRH lifetimes x 1e6 the grid and adaptive
CCE agree to 3.4e-10 (vs 8.0e-6 with SRH): the depleted-epi SRH of the excess pairs,
U ~ dn dp / (tau_p dn + tau_n dp), homogeneous of degree 1 (LET-independent CCE) but
not linear, so the coarse 1 ps injection step / 10 % grid misjudge the electron-hole
overlap while the pairs separate.

## 4. Lateral cut (bench/lateral_cut.py, notes/lateral.out; hw 50, LET 15.8)

Cut at W um beyond the track (reflecting side, half outer cell = box-exact):
W = 3 / 5 / 10 / 20 / full: CCE 0.9999823218 / ..217 / ..209 / ..210 / ..209,
t_collection 0.1097680 / 0.1097784 / 0.1097782 / 0.1097790 / 0.1097790 ns. At hw 150
the 10 um cut is -7.5e-10 off the full width (parity.out). Default 10 um: 8208 nodes at
hw 50 (vs 17784), 10944 at hw 150 (vs 49248).

## 5. Mesh (bench/mesh_convergence.py, notes/mesh.out; hw 50, cut 10 um, LET 15.8)

devsim's mesher with every mesh-line spacing / f (`_mesh_line_scale`):

| f                 | nodes (nx x ny) | CCE          | t_collection (ns) |
| ----------------- | --------------- | ------------ | ----------------- |
| 1 (devsim's mesh) | 12 x 684        | 0.9999823209 | 0.1097782         |
| 1.41              | 16 x 967        | 0.9999828384 | 0.1096225         |
| 2                 | 22 x 1369       | 0.9999832705 | 0.1095365         |
| 4                 | 42 x 2737       | 0.9999838146 | 0.1094415         |

Increments 9.5e-7 (1 -> 2), 5.4e-7 (2 -> 4): roughly first order; extrapolated limit
~0.999984-0.999985, so devsim's mesh is ~2e-6 low in CCE and +0.3 % in t_collection.
Both engines share this error (parity uses devsim's mesh); not changed here.

## 6. Timing

Per CCE(LET) point, single thread (notes/timing.out; NUMBA / OMP / BLAS threads 1, as the
n_jobs workers run; CPU seconds, load 32-36 during the run):

| point (LET 15.8, 50 V) | nodes | fast adaptive | fast devsim-grid twin | devsim |
|---|---|---|---|---|
| hw 50, 10 um cut | 8208 | 14.7 s (238 steps) | 7.4 s (81 steps) | |
| hw 150, 10 um cut | 10944 | 21.1 s | 11.6 s | |
| hw 50, full width | 17784 | 35.5 s | 18.0 s | 62.9 s (80 steps) + 3.6 s build |

At load ~7 the same adaptive hw 50 cut point takes 4.3-5.1 s (newton_tol 1e-10, now the
default: 10 % fewer solves than 1e-11, CCE moved 5e-12). Profile: residual evaluations
~50 %, SuperLU solves ~25 %, factorisations ~15 %. Build (mesh import + fv2d dark
state): 0.4-0.9 s CPU at the cut, 2.3 s at hw 50 full.

Whole committed tables (bench/rebuild_tables.py, n_jobs 6, notes/rebuild_*.out):

| table (notebook-16 keywords) | adaptive, load ~7 | adaptive, load ~45 | devsim-grid twin, load ~7 | devsim |
|---|---|---|---|---|
| 100 um (hw 50, 20 LETs) | 41.0 s | 239 s | 10.1 s | ~15 min (docs/USAGE.md 12, 2026-10-18, not re-measured) |
| 300 um (hw 150, 10 LETs) | 43.7 s | 141 s | 15.6 s | ~26 min |

The twin reproduces the committed tables to <= 1e-6 in CCE at LET >= 13 keV/um (1.1e-8
at LET 500 for 100 um, 4.2e-8 for 300 um; same synthesised rows, method "linear") and by
the dark-current offset below that (up to 8.0e-6 / 1.72e-5 at LET 0.5).

## 7. Decision on data/cce_let_table_*.json: NOT regenerated

The rule was: regenerate only if the fast engine reproduces the committed CCE to 1e-6.
It does not, and the difference is in devsim's numbers, not the fast engine's:
the dark-current offset (up to 1.7e-5 at LET 0.5, hw 150) and the 10 % grid time error
(-8.15e-6 at every LET). The devsim-grid twin reproduces devsim's discrete solution
(5e-10 where the offset is negligible), the adaptive scheme the converged one. A rebuild
with the adaptive engine would RAISE every committed CCE by 8e-6 (more at low LET); that
is a physics/data decision for the maintainers (bench/rebuild_tables.py writes both
variants to a scratch directory). The t_collection_ns column of the committed files is
already stale (single_particle docstring).

## 8. Tests and the suite critical path

Default suite (`-m "not exhaustive"`), serial wall times at load 4-8 (2026-10-21):

| before (devsim) | s | after | s |
|---|---|---|---|
| TestBuildCCELetTableSmallMesh: fixtures (full + linear tables) | 58.6 | tests/test_single_particle_fast.py (whole file, module fixtures, 7 tests, each <= 6.4 s) | 18.1 |
| ... test_build_once_and_linear_parity | 6.6 | test_single_particle_fast_parity.py::test_devsim_grid_twin_matches_devsim_and_build_once | 21.6 |
| ... test_parallel_table_equals_serial[linear-3] | 16.9 | test_single_particle_fast_parity.py::test_devsim_grid_twin_series_at_10V | 8.1 |
| TestCCEVsBias | 14.4 | | |
| total | 96.5 | total | 47.8 |

Full suite, `uv run pytest -n 6 --durations`, same machine, load 5-11 (2026-10-21):

| | before (devsim classes in the default suite, main's schedule) | after |
|---|---|---|
| suite wall | 197.2 s (1983 passed, 8 skipped) | 186.8 s (1989 passed, 8 skipped) |
| longest unit | TestBuildCCELetTableSmallMesh ~140 s (fixtures 101.6 + parallel 35.8 + build-once call) | test_engine2d failing_config 65.7 s |
| single-particle units | + TestCCEVsBias 25.1 s, legacy 22.7 s | fast file ~45 s (fixtures 19 + error-estimate test 20 + ...), parity 39.4 s and 16.1 s, legacy 20.5 s |

The suite wall is bound by the total work over 6 workers, so it drops by only 5 %; the
longest unit (the critical path once the other units are spread) halves. Under xdist load
the devsim parity test (3 devsim transients: the 2-LET shared-device table and the fresh
device of the build-once check) takes 39 s, above the 10 s per-test guideline: it is the
one devsim check kept on purpose. `test_single_particle_legacy.py` (devsim, hw 25,
legacy stack, cathode contact) is at 20-23 s and was not switched: the legacy stack has
no devsim parity test in the fast files, and adding one would cost about what it saves.

Both devsim classes stay, marked `exhaustive` (run with `-m exhaustive`; the [full-2]
parametrization was already exhaustive: 20.8 s). Under xdist the largest of these units
drops from ~82 s (HEAVY: 160 s under load) to 21.6 s; `tests/_xdist_schedule.py` updated
(the parity file is split per test). The suite's critical path is now set by other units:
`test_engine2d.py::test_sweep_row_failing_config_hw25` (38.4 s serial at load 4; HEAVY 150)
is a devsim DC ramp to 80 V on the legacy (no n+) stack that fails at 67.4 V after 1873
Newton iterations in 25 solves (70 % UMFPACK numeric factorisation, mostly the failing
tiers past 57.5 V; the engine2d row takes 0.2 s). Not a transient: not accelerated here;
running the devsim reference at hw 2 (laterally uniform failure) would be the cheap
option, not done (it changes what the test pins).

## 9. Gaps

- Mesh adaptivity is an estimate (devsim's mesher refined uniformly), not a refinement
  loop: the parity mesh is devsim's, whose ~2e-6 CCE error both engines share.
- Planar `create_sic_2d_device` only (anode + cathode, Cartesian, laterally uniform for
  the cut); the alternative structures (fv2d builders) and the axisymmetric devices are
  not wired, nor the incomplete-ionization option (fv2d has none).
- `test_engine2d.py::test_sweep_row_failing_config_hw25` (~80-150 s) is a devsim DC
  ramp failure, not a transient: this engine does not touch it.
