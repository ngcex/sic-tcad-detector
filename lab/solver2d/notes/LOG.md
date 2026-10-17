# lab/solver2d running log (2D numba engine spike)

Machine: shared M4 10-core; load average 6-25 during this work. All commands from repo root.

## STEP 0 -- devsim 2D mesh structure (2026-10-16)

`uv run python lab/solver2d/bench/dump_mesh.py 25` / `... 50` -> notes/mesh_hw{25,50}.npz
`uv run python lab/solver2d/bench/check_mesh.py lab/solver2d/notes/mesh_hw50.npz`

| hw (um) | nodes | nx x ny | edges (hor/ver/diag) | triangles |
|---|---|---|---|---|
| 25 | 3924 | 12 x 327 | 11095 (3597/3912/3586) | 7172 |
| 50 | 8502 | 26 x 327 | 24801 (8175/8476/8150) | 16300 |
| 150 | 23544 (from ref npz) | 72 x 327 (assumed same y grid) | | |

- Tensor-product grid: every node at (xs[i], ys[j]); devsim node order is NOT
  j*nx+i nor i*ny+j (devsim inserts row by row) -> map by coordinates.
- One diagonal per rectangle, all the same orientation ((i,j+1)-(i+1,j)).
  Diagonal |EdgeCouple| <= 2.8e-33 cm (vs >= 6.6e-8 for axis edges) -> exact 5-point FV stencil.
- NodeVolume == rectangular dual cell (dx_dual*dy_dual) to 2.2e-16 rel;
  EdgeCouple(horizontal) == dy_dual, EdgeCouple(vertical) == dx_dual to 2e-16 rel.
- Region "sic" spans y = 0..11 um; air_top / air_bot are separate 24-node regions with
  NO equations and no interface -> inert. Contacts: anode = full row y=0 (nx nodes),
  cathode = full row y=11 um (nx nodes). Lateral boundaries x=0, x=hw: natural Neumann.
- The 327 y-lines are bit-identical to the 1D devsim mesh (lab/gpu_solver/bench/mesh_devsim_default.npy).
- ref_50_50.npz (dark, V=50): lateral spread over the 26 columns: Potential 7e-15 V,
  log10 n 4.8e-10, log10 p 3.9e-10 -> the dark 2D state is the 1D solution on the y grid
  (doping depends on y only, full-width contacts, Neumann sides).
- Extended dump (coordinator request): notes/mesh_hw{25,50}.npz hold per node x, y, NodeVolume,
  Donors/Acceptors/NetDoping, Potential/Electrons/Holes (V=0); per edge endpoints, EdgeCouple, EdgeLength;
  triangle list; contact node lists (anode/cathode rows). notes/scan_hw25_V50.npz and scan_hw50_V50.npz
  (`bench/devsim_scan_dump.py HW 50 NPTS`) hold the dark V=50 state, edge currents + generation node
  models at the first/last stripe, get_contact_current (e, h) per stripe, the triangle I_gen and CCE.
- Column-to-column spread at V=50 (ref_50_50, dark): Potential 7.1e-15 V, log10 Electrons 4.8e-10,
  log10 Holes 3.9e-10. Dark state is x-independent -> supports the spec's inference that the
  edge-vs-centre CCE difference (0.998369 vs 0.998495 at hw=50, 1.3e-4) is numerical (mesh / quadrature
  / uncollected contact-row generation), not physics. (Inference; not separately proven.)

## devsim conventions pinned (hw=25, V=50, stripe at x=0)

- get_contact_current(cathode) == sum over cathode-row nodes of EdgeCouple-weighted edge flux, with
  NO NodeVolume x node-model term: e 4.0588184715343e-08 (devsim) vs 4.0588184715343813e-08 (edge sum);
  edge+node would be 4.0639e-08. So generation inside the cathode-row half cells (G up to 9.1e17 there,
  alpha range 15 um > 10 um epi) is in I_gen but never collected: ~1.3e-3 systematic CCE deficit.
  Reproduced, not fixed.
- I_gen: integrate_over_mesh_2d triangle-vertex mean; engine reproduces it bit-for-bit
  (250838050012.9449 both) using devsim's element list mapped to k.
- psi moves by 1.4e-6 V (centre stripe) / 2.0e-6 V (edge) between dark and illuminated -> NOT < 1e-8,
  so psi cannot be frozen; the scan is a full Gummel (n, p, quasi-Fermi Poisson) per stripe.

## STEP 1 -- equilibrium parity

`uv run python -m pytest lab/solver2d/tests -q -p no:cacheprovider` (5 tests, ~1.2 s)
- NetDoping vs devsim: 4.4e-16 rel.
- Equilibrium Potential, full 2D Newton from neutral guess (12 it): max |dpsi| 6.1e-16 V (hw=25), 7.2e-16 V (hw=50).
- Same via 1-column solve + 2D polish: 5.6e-16 / 7.2e-16 V; log10 n, p 1e-14.

## STEP 2 -- reverse bias + lateral scan

Engine design: banded LU without pivoting (numba, half-bandwidth nx, x fastest), one scalar equation at a time.
Dark ramp on a 1-column grid (exact for this laterally uniform problem) + one 2D Newton polish whose LU is
kept; electron and hole Jacobians factored once at the dark state; each stripe = chord-Gummel (exact residuals,
frozen dark LUs), 4 iterations to |du|<1e-10 V_t and |dI/I|<1e-9.

Dark V=50 vs ref_*_50.npz: hw=50 max|dpsi| 1.4e-14 V, log10 n 3.6e-11, log10 p 1.4e-13 (where > 1 cm^-3);
hw=150 (72 x 327): 2.1e-14 V, 3.2e-10, 1.3e-13.
Illuminated (stripe x=0, hw=50) vs devsim: psi 2.8e-14 V, log10 n 4.3e-10, log10 p 5.1e-14.

### GO/NO-GO timing (2026-10-16 ~17:15, load avg 2.7-2.9, no pytest running)

`uv run python lab/solver2d/bench/bench_go_nogo.py 50 150` -> notes/bench_go_nogo.out
Engine single-threaded (no prange), JIT cache warm (warm-up incl. a hw=25 run 0.20 s; fresh compile not measured
separately, est. a few s -- guess). Engine time INCLUDES building devsim's mesh-only device.

| hw | devsim build+ramp | devsim scan | devsim total | engine mesh | dark | factor n,p | scan | engine total | speedup | max dCCE |
|---|---|---|---|---|---|---|---|---|---|---|
| 50 | 10.10 s | 2.85 s | 12.95 s | 0.030 | 0.004 | 0.005 | 0.055 | 0.093 s | 139x | 1.45e-7 |
| 150 | 38.14 s | 10.32 s | 48.46 s | 0.088 | 0.056 | 0.102 | 0.460 | 0.706 s | 69x | 1.11e-7 |

Engine re-run right after devsim: 0.095 s / 0.711 s (no drift). => GO (>= 5x at both widths, parity met).

### Re-run with the geometry-general variant (17:18, load 1.6-1.8)

bench_go_nogo.py now also runs `dark(V, full2d=True)`: the whole ramp (equilibrium + 5 V secant steps)
with 2D Newton + 2D banded LU each iteration (81 2D LUs), i.e. without the lateral-uniformity shortcut.
CCE arrays (engine, engine full2d, devsim cce_lateral_scan itself) saved to notes/cce_ref_lateral_scan.npz
(this, not scan_hw50_V50.npz -- which is my own loop with re-solved baselines -- is the CCE reference;
the two devsim CCE sets agree to the printed 8 digits).

| hw | devsim total | engine (1-column ramp) | speedup | engine (full 2D ramp) | speedup | max dCCE |
|---|---|---|---|---|---|---|
| 50 | 12.95 s | 0.093 s | 139x | 0.295 s (dark 0.204) | 44x | 1.45e-7 |
| 150 | 48.01 s | 0.706 s | 68x | 4.92 s (dark 4.27) | 9.8x | 1.11e-7 |

2D LU cost: ~2.5 ms (hw=50, N=8502, w=26), ~52 ms (hw=150, N=23544, w=72) -> N w^2 scaling.
Dark cathode current at V=50, hw=50: devsim 1.33e-15 A/cm, engine 6.29e-15 A/cm (engine dark state is
flat quasi-Fermi, not DD-consistent; irrelevant to CCE, which does no dark subtraction and whose
illuminated state is solved with full Gummel). Systematic CCE offset -0.5..-1.5e-7 vs devsim: unexplained
(guess: devsim's rel 1e-8 stopping on a current that is a small difference of large SG terms).

## Follow-up: is edge-vs-centre CCE a mesh artefact? (17:24, load ~4.7)

`uv run python lab/solver2d/bench/mesh_artefact.py` -> notes/mesh_artefact.out. Engine only, V=50, sigma 2 um,
same y lines/doping/BCs. Uniform-x grids built directly as Grid2D (same diagonal orientation); rebuilding the
production grid this way reproduces its scan to 2.2e-16.

| mesh (hw=50) | centre CCE | mid (x=hw/2) | edge CCE | edge/centre-1 |
|---|---|---|---|---|
| production graded | 0.99849483 | -- | 0.99836869 | -1.26e-4 |
| uniform 0.5 um | 0.99848689 | 0.99847356 | 0.99846023 | -2.67e-5 |
| uniform 0.25 um | 0.99848023 | 0.99847356 | 0.99846690 | -1.34e-5 |
| uniform 0.125 um | 0.99847689 | 0.99847357 | 0.99847023 | -6.68e-6 |
(hw=25: production -1.62e-4; uniform 0.5/0.25/0.125 um: -2.67e-5/-1.34e-5/-6.68e-6, same values.)

On uniform meshes the remainder is O(h) and antisymmetric (centre +h, edge -h around the mid-stripe value
0.99847356). That points at I_gen's triangle quadrature (integrate_over_mesh_2d): with one diagonal orientation,
area/6*(v00+2v10+2v01+v11) differs from the box rule by A/12*(v10+v01-v00-v11) ~ -h_x h_y dG/dx dG/dy. The sign of
dG/dx flips between a stripe at x=0 and one at x=hw, and the term cancels at mid-domain.
Check, I_gen = sum(NodeVolume*G) (the box rule the solver uses for the source):
  production mesh: centre 0.998473564, mid 0.998473560, edge 0.998473548 -> edge/centre-1 = -1.6e-8
  uniform 0.5 um: -2e-10; uniform 0.25 um: +2e-10; centre/mid-1 about 3-5e-9 everywhere.
Conclusion: PROVEN numerical. Using the box rule for I_gen collapses edge/centre to 1 within 2e-8 even on the production
mesh. The whole 1.3e-4 edge/centre difference comes from the I_gen triangle quadrature (mesh-dependent, O(h)). It is not
the solver's lateral resolution, not the Neumann mirror (a half-Gaussian at x=0 or x=hw = mid stripe to 5e-9), and
not the cathode-row loss (laterally uniform; it gives the common ~1.5e-3 deficit, 0.998474 everywhere).

## Promotion to etna/core/engine2d (2026-10-16 evening)

Package: etna/core/engine2d/{grid,banded,kernels,device,scan}.py; entry points
`engine2d_unsupported_reason(...)`, `lateral_cce_scan_2d(...)`; integration only in
etna/core/optimization.py (`microdosimetric_sweep(..., solver="devsim"|"auto")`, records get a
"solver" column). Tests: tests/test_engine2d.py (16 fast + 3 slow, all pass; slow 21 s).

Changes vs the lab spike:
- I_generated = sum(NodeVolume*G) (box rule), matching compute_cce_2d after fd631e3.
- Mesh: still devsim's mesh-only build (create_sic_2d_device, no equations, ~30-90 ms); not
  regenerated in numpy. Every engine assumption is ASSERTED on that device (grid_from_arrays):
  full tensor grid, NodeVolume == rectangular dual, Donors/Acceptors identical across columns,
  anode/cathode == full rows j=0 / j=ny-1. Doping is taken from devsim's evaluated node models.
- Cathode current evaluated in a balanced (well-conditioned) form: electrons from the summed
  interior balance sum V(G-U) + electron flux at the anode (minority there), holes by the direct
  cathode edge flux (minority there). Same quantity as get_contact_current at convergence. The
  direct edge form (devsim's) is a small difference of huge SG terms when the epi next to the
  cathode is quasi-neutral: e.g. hw=25, epi 20 um, N_D 5e14, 80 V: point-to-point scatter of the
  10 CCE values 1.4e-6 (direct, also seen in devsim) vs 2.6e-8 (balanced); default 50 V: 3.1e-8 vs 4.1e-9.
- Punch-through: devsim walls measured (`bench/devsim_wall.py 25 5e13 8.5e13 5e14`, 0.5 V steps
  from 50 V): last converged 54.5 V (5e13), 56.5 V (8.5e13), 80 V reached (5e14). The quantity that
  lines up at the wall is the dark field on the last cell before the cathode:
  4.745e4 / 4.627e4 V/cm at the last converged points, 4.795e4 / 4.677e4 at the failing ones;
  converged at 4.30e4 (5e13, 50 V) and 3.20e4 (5e14, 80 V). (Not n/N_D: flat-qF n underflows in
  the depleted epi, and "n near the cathode" is mesh-position dependent.)
  Rule: field <= 4.45e4 -> engine answers; > 5.0e4 -> Engine2DPunchThrough (caller writes the
  same NaN row devsim produces); in between -> Engine2DFallback (devsim decides). Calibrated at
  hw=25, epi 10 um only (guess: devsim's wall may move with hw/epi; the gate checks hw 25-150).
- Gummel non-convergence or a Newton failure -> Engine2DFallback (devsim decides).
- Spot checks via _sweep_one, auto vs devsim (hw=25): epi 5 um 20 V dC 2e-9; epi 5 um 50 V both NaN
  (engine punch-through 0.02 s vs devsim failure 11.2 s); epi 20 um 50 V dC 2.8e-8; epi 20 um 5e14 80 V
  dC 5.2e-7, dE 2.1e-6, dRatio 1.6e-6 (devsim's direct-current noise; see above); epi 10 um 5e14 80 V dC 5.2e-7.

## GATE: notebook 20 grid, 36 configs (hw 25/50/100/150, epi 10 um, N_D 5e13/8.5e13/5e14, V 20/50/80)

`uv run python lab/solver2d/bench/gate_nb20.py --jobs 4` (load 2.8 at start, no pytest), then
`uv run python lab/solver2d/bench/gate_analyse.py`. CSVs: notes/gate_devsim.csv, gate_auto.csv
(gate_auto_final.csv = rerun with final code: identical to 1e-16).

- Wall time, n_jobs=1: devsim 1187.4 s, auto 8.81 s (135x). auto n_jobs=4: 4.57 s, identical results.
- Per config (median by hw, devsim / engine): 25: 5.15 / 0.042 s; 50: 13.5 / 0.104; 100: 28.5 / 0.288;
  150: 47.9 / 0.663 s. The 8 failing configs: devsim 10.6-155 s each, engine 0.03-0.15 s.
- Cold JIT (empty NUMBA_CACHE_DIR, fresh process): first config 0.94 s (warm cache: 0.12 s);
  `import etna.core.optimization` ~1 s either way. Spawned workers pay that import, which is why
  n_jobs=4 gives only 1.9x on a 9 s run.
- All 36 rows were served by engine2d (no fallbacks). NaN rows identical (8: V=80, N_D 5e13 / 8.5e13,
  every hw). devsim ramp failures at 54.25 V (5e13, all hw) and 58.0 / 58.0 / 56.6 / 56.6 V (8.5e13,
  hw 25/50/100/150); engine: Engine2DPunchThrough, cathode-cell field 5.29e4 / 5.48e4 V/cm.
- is_valid identical (also is_fully_depleted, passes_cce_floor).
- max |dCCE| center 1.17e-5 (> 1e-5 target), edge 9.4e-6, ratio 2.7e-6, cce_std 2.5e-6; signed mean
  dCCE -1.1e-6 (engine lower).
- Strict ranking NOT identical; tie-aware ranking (ratio rounded to 5 decimals, grid-order tiebreak)
  identical. The valid-group ratios span 3.0e-6 in devsim vs 1.5e-8 in the engine: after the box-rule
  fix every valid ratio is 1 and the strict order is devsim's numerical noise.
- The one row beyond 1e-5: hw=150, N_D 5e14, 20 V (partially depleted, CCE 0.816, is_valid False).
  The differences grow with hw: 2.7e-6 (50), 6.8e-6 (100), 1.17e-5 (150). Engine centre CCE for this
  config is 0.81611930 / 0.81611930 / 0.81611929 at hw 25 / 50 / 150 (hw-independent to 1e-9);
  devsim gives 0.81612016 (hw 25) and 0.81613094 (hw 150). Three extra devsim Newton iterations at
  rel 1e-14 leave the devsim value unchanged (stagnated). The engine's own direct edge-flux
  evaluation also drifts with hw (-3.6e-6 at 150) while the balanced form does not. Conclusion: devsim
  has ~1e-5 roundoff error in this partially depleted (quasi-neutral cathode side) regime; the engine
  is the more accurate side.

## Box widening (TASK 1, 2026-10-16 ~19:00, load 1.6-4)

Walls: `bench/devsim_wall2.py HW EPI_UM V0 V1 ND...` (devsim 0.5 V steps from a converged bias, hw=25):
| epi um | N_D_bulk | last converged V (field V/cm) | first failed V (field) |
|---|---|---|---|
| 5 | 5e13 | 26.5 (4.632e4) | 27.0 (4.732e4) |
| 5 | 8.5e13 | 27.0 (4.585e4) | 27.5 (4.685e4) |
| 5 | 5e14 | 37.0 (4.821e4) | 37.5 (4.921e4) |
| 7.5 | 5e13 | 40.0 (4.670e4) | 40.5 (4.736e4) |
| 7.5 | 5e14 | 62.5 (4.661e4) | 63.0 (4.728e4) |
| 10 | 5e13 | 54.5 (4.745e4) | 55.0 (4.795e4) |
| 10 | 8.5e13 | 56.5 (4.627e4) | 57.0 (4.677e4) |
| 15 | 5e13 | 80 reached (4.476e4) | - |
| 20 | 5e13 | 80 reached (2.952e4) | - |
| 20 | 8.5e13 | 80 reached (2.305e4) | - |

Same field rule separates converged/failed at every epi (failures >= 4.685e4, converged <= 4.821e4);
thresholds unchanged (engine <= 4.45e4 < devsim decides < 5.0e4 < punch-through NaN).

Parity grid: `bench/parity_box.py START STOP` (3 chunks, each < 3 min) -> notes/parity_box.csv,
_sweep_one solver='auto' vs 'devsim', 24 points:

```
 hw  epi_um           nd  V solver_a  center_d            dC            dE            dR  valid_d  valid_a   t_d    t_a
 25     5.0 5.000000e+13  5 engine2d  0.983523 -4.922353e-08 -3.341402e-08  1.607436e-08     True     True  2.31  0.092
 25     5.0 8.823059e+13 20 engine2d  0.994375  2.096925e-09  6.354144e-10 -1.469777e-09     True     True  2.67  0.092
 25     5.0 5.000000e+14 35 engine2d  0.994435  2.020600e-09  2.257342e-09  2.380667e-10     True     True  3.14  0.029
 25     5.0 5.000000e+14 37   devsim       NaN           NaN           NaN           NaN    False    False  4.87  4.873
 50     5.0 5.000000e+13 10 engine2d  0.993400 -5.057928e-08 -4.178771e-08  8.849982e-09     True     True  6.80  0.074
 25     5.0 8.500000e+13 50 engine2d       NaN           NaN           NaN           NaN    False    False 12.67  0.020
 25     7.5 5.000000e+13 20 engine2d  0.995550 -2.418130e-09 -1.463040e-09  9.593597e-10     True     True  3.25  0.106
 25     7.5 5.000000e+14 65   devsim       NaN           NaN           NaN           NaN    False    False 10.42 10.450
 25     7.5 8.823059e+13 35 engine2d  0.995736 -1.681085e-10 -1.767736e-10 -8.702372e-12     True     True  3.59  0.035
 50     7.5 5.000000e+14 50 engine2d  0.995773 -2.563675e-08 -2.916150e-08 -3.539715e-09     True     True 12.42  0.091
 25    10.0 5.000000e+13  5 engine2d  0.864047  5.994777e-08  1.521038e-08 -5.177656e-08    False    False  3.37  0.043
 25    10.0 8.823059e+13 10 engine2d  0.963119  1.747448e-07  1.476699e-07 -2.811174e-08    False    False  3.91  0.039
 25    10.0 5.000000e+14 35 engine2d  0.924837 -5.401995e-07 -5.414102e-07 -1.309047e-09    False    False  5.24  0.040
 50    10.0 8.500000e+13 65 engine2d       NaN           NaN           NaN           NaN    False    False 28.15  0.053
 25    10.0 5.000000e+14 75 engine2d  0.998299 -4.673510e-07  5.020491e-07  9.710516e-07     True     True  6.31  0.041
 25    15.0 5.000000e+13 80   devsim  0.999738  0.000000e+00  0.000000e+00  0.000000e+00     True     True  6.25  6.332
 25    15.0 8.823059e+13 50 engine2d  0.999635  4.663177e-09 -1.092284e-08 -1.559171e-08     True     True  6.11  0.051
 25    15.0 5.000000e+14  5 engine2d  0.666305  1.985336e-07  2.820717e-07  1.253749e-07    False    False  3.75  0.054
 50    15.0 5.000000e+13 35 engine2d  0.999497 -3.059433e-08 -2.814684e-08  2.448714e-09     True     True 15.61  0.123
 25    20.0 5.000000e+13 80 engine2d  0.999746 -7.084600e-11 -5.440359e-11  1.644662e-11     True     True  7.98  0.056
 25    20.0 8.823059e+13 20 engine2d  0.998949 -3.312713e-08 -4.884073e-08 -1.573014e-08    False    False  6.45  0.057
 25    20.0 5.000000e+14 10 engine2d  0.783241 -3.756942e-07 -2.997107e-07  9.701117e-08    False    False  5.22  0.061
 50    20.0 5.000000e+14 65 engine2d  0.990542  4.103059e-06  4.246861e-06  1.451746e-07    False    False 25.50  0.145
 25    20.0 5.000000e+13  5 engine2d  0.899261  1.689414e-09  3.856723e-09  2.410107e-09    False    False  4.79  0.063
```

max|dC| 4.10e-06, max|dE| 4.25e-06, max|dRatio| 9.71e-07; NaN rows identical True; is_valid identical True; solver counts {'engine2d': 21, 'devsim': 3}
devsim-decides points: (25,5um,5e14,37 V) devsim NaN (its adaptive ramp fails below the 0.5 V-step wall);
(25,7.5um,5e14,65 V) devsim NaN; (25,15um,5e13,80 V) devsim converges (field 4.476e4, just above ENGINE_MAX).
Largest dC 4.1e-6 at (50,20um,5e14,65 V, partially depleted, valid False); max ratio diff ~1e-6 from devsim's
edge noise at (25,10um,5e14,75 V): devsim ratio 0.99999903, engine 1.0000000003.

Final ENGINE2D_BOX: hw 25-150, epi 5-20 um, N_D_bulk 5e13-5e14, V 5-80. hw 100/150 evidence only from the
36-config gate (epi 10); other epi only at hw 25/50 (lateral uniformity is asserted per device).

## TASK 2: run_field 2D -> engine (2026-10-16 ~19:30)

etna/api/simulation.py: run_field tries the engine first for 2D configs eligible under
`_field_engine_unsupported_reason` (new helper -> etna.core.engine2d.engine2d_dark_unsupported_reason:
graded, N_A 1e19, T 300, substrate 1 um, N_D None, hw 25-150, epi 5-20 um, N_D_bulk 5e13-5e14,
N_D_junction 1e15-1e16, L_transition 0.3-3 um, reverse bias 5-80 V). Engine: `dark_state_2d` (new, scan.py) =
mesh-only devsim build + 1-column ramp + 2D polish, returned in devsim's node order. Above
FIELD_ENGINE_MAX (4.45e4 V/cm) or on any engine error -> devsim path (unchanged; reset_devsim_fully only
there). The node |E| code was factored into `_node_field_2d` and is shared by both paths.
metadata["solver"] = "engine2d" | "devsim".

Parity (DeviceConfig(half_width_um=hw), engine vs forced devsim path, same process):
| hw | bias | engine | devsim | x/y order | NetDoping | max dPotential | max rel dE (floor 1e-3 Emax) |
|---|---|---|---|---|---|---|---|
| 25 | -20 | 0.190 s | 2.50 s | identical | identical | 5.3e-15 V | 4.5e-13 |
| 25 | -50 | 0.128 s | 3.67 s | identical | identical | 7.1e-15 V | 2.8e-13 |
| 50 | -20 | 0.177 s | 6.66 s | identical | identical | 5.3e-15 V | 4.5e-13 |
| 50 | -50 | 0.176 s | 9.78 s | identical | identical | 1.4e-14 V | 3.4e-13 |
Engine time is mostly the LinearNDInterpolator |E| post-processing (same code as devsim path).
Tests: tests/test_api_field_engine.py (10 fast + 4 slow), tests/test_api_field.py, tests/test_app_field_page.py,
tests/test_engine2d.py fast: all pass (slow field files 11.6 s).
Note: tests/test_api_field.py's 2D test (hw=50, -10 V) now runs the engine path; the devsim 2D path is covered
by the forced-devsim parity and fallback tests in test_api_field_engine.py.

## n+ substrate spike (2026-10-16 ~20:00): is the devsim wall a contact artefact?

Scripts: lab/solver2d/bench/nplus.py (1D; monkeypatches etna.core.device.build_diode_mesh and
set_graded_doping_profile in-process), nplus2d.py (2D; patched copy of create_sic_2d_device installed into
etna.core.charge_collection_2d). Variants: base = production; mesh = epi end refined (ps 1 nm at the epi end,
+1 um region) but NO n+ doping; np19 / np18 = 1 um n+ at 1e19 / 1e18 under the epi, cathode on the n+ end.
1D device = create_dd_device graded defaults (N_D_junction 2.9e15, L 1 um), N_D_bulk 8.5e13 (= default) and 5e13.

### Walls (last converged reverse bias)
1D, devsim adaptive ramp (ramp_bias V_step_max=5, to 150 V): base, np19, np18 all reach 150 V (both N_D) --
the adaptive secant ramp steps over the wall even without n+.
1D, fixed 0.5 V steps from 0 V: base 67.5 V (8.5e13) / 54.5 V (5e13); mesh-only 77.0 / 60.0 V;
np19 and np18: 150 V reached (both N_D).
2D hw=25, fixed 0.5 V steps from 0 V:
| epi | N_D_bulk | base | mesh-only | n+ 1e19 |
|---|---|---|---|---|
| 5 um | 5e13 | 26.5 V (LOG table) | - | 150 V reached (112 s) |
| 5 um | 5e14 | 37.0 V (LOG table) | - | 150 V reached (112 s) |
| 10 um | 5e13 | 54.5 V (LOG table) | 60.0 V | 150 V reached (160 s) |
| 10 um | 5e14 | 94.0 V (fails 94.5) | - | 150 V reached (163 s) |
=> the wall is a contact/mesh artefact: refining the mesh alone moves it (54.5 -> 60 V 2D, 54.5 -> 60 V 1D),
an n+ layer removes it (>= 150 V everywhere tested).

### Calibrated quantities, n+ vs base (1D)
W: etna rule (first n-side node with n/Donors >= 0.5, clipped to the epi thickness); C = eps*0.05/W.
| N_D_bulk | V | W base (um) | W np19 | W np18 | dC/C np19 |
|---|---|---|---|---|---|
| 8.5e13 | 0 | 1.69708 | 1.69708 | 1.69708 | 0 |
| 8.5e13 | -5 | 5.775706 | 5.775733 | 5.775735 | -5e-6 |
| 8.5e13 | -10 | 9.592595 | 9.213316 | 9.216782 | +4.1e-2 |
| 8.5e13 | -30 | 9.984307 | 9.895584 | 9.901486 | +9.0e-3 |
| 8.5e13 | -50 | 9.986025 | 9.940143 | 9.946114 | +4.6e-3 |
| 5e13 | -5 | 7.291363 | 7.269091 | 7.269334 | +3.1e-3 |
| 5e13 | -10 | 9.926384 | 9.588434 | 9.59383 | +3.5e-2 |
| 5e13 | -50 | 9.986084 | 9.940636 | 9.946658 | +4.6e-3 |
(With n+ the electron tail of the n+/n- high-low junction (Debye length in n- ~0.4 um) keeps n/N_D >= 0.5 a
few hundred nm before the epi end, so W saturates at ~9.94 um instead of ~9.99.)

CCE (cce_vs_bias replica, G = alpha profile 15 um range, 1e18 peak): devsim's cathode contact current is
unusable with n+ (a small difference of huge SG terms in the 1e19 region: CCE_direct jumps by 1-2 %, >1).
Evaluated with the balanced form Ie = sum_interior V*ElectronGeneration - Jn(anode edge), Ih direct
(base: balanced == direct to 5e-8):
| N_D | G domain | V=-10 | -20 | -30 | -40 | -50 |
|---|---|---|---|---|---|---|
| 8.5e13 base | (epi = device) | 0.966579 | 0.996169 | 0.997735 | 0.998272 | 0.998475 |
| 8.5e13 np19 | epi only | 0.998867 | 0.999562 | 0.999708 | 0.999778 | 0.999820 |
| 8.5e13 np18 | epi only | 0.998857 | 0.999560 | 0.999706 | 0.999777 | 0.999820 |
| 8.5e13 np19 | epi + n+ (alpha reaches it) | 0.961164 | 0.961942 | 0.962095 | 0.962170 | 0.962216 |
| 5e13 base | | 0.990597 | 0.996812 | 0.997936 | 0.998340 | 0.998498 |
| 5e13 np19 | epi only | 0.999065 | 0.999563 | 0.999705 | 0.999776 | 0.999819 |
n+ acts as a back-surface field (holes generated near the cathode are no longer lost into the ohmic contact):
+3.3 % at -10 V, +1.3e-3 at -50 V (8.5e13). If the alpha deposits in the n+ too, CCE drops by ~3.7 %.

Dark current (dark_current_sweep, N_t 2.8e12, area 0.04 cm^2; direct extraction also broken with n+:
4.7e-9 A flat for np19, so balanced form shown; base balanced vs direct differ by 2e-3 at -10 V, 3e-4 at
-30 V, 2e-6 at -60 V -- the SRV contact term):
| N_D | -10 V | -30 V | -60 V |
|---|---|---|---|
| 8.5e13 base (direct, official) | 5.948897e-12 | 1.457471e-11 | 1.795679e-11 |
| 8.5e13 base (balanced) | 5.960995e-12 | 1.457048e-11 | 1.795675e-11 |
| 8.5e13 np19 (balanced) | 6.110770e-12 (+2.5 %) | 1.471496e-11 (+1.0 %) | 1.801384e-11 (+0.32 %) |
| 8.5e13 np18 (balanced) | 6.104578e-12 | 1.471070e-11 | 1.802023e-11 (+0.35 %) |
| 5e13 base (balanced) | 6.056268e-12 | 1.485636e-11 | 1.795674e-11 |
| 5e13 np19 (balanced) | 6.222897e-12 (+2.8 %) | 1.500854e-11 (+1.0 %) | 1.801366e-11 (+0.32 %) |

### Decision rule
Wall moves beyond 100 V: YES (>= 150 V, 1D and 2D). Calibrated range changes < 1e-3: NO (C at -10 V +4.1 %,
-30 V +0.9 %; CCE -10 V +3.3 %, -50 V +1.3e-3; dark current -60 V +0.32 %, -10 V +2.5 %).
=> NOT a pure validity extension. Redo: graded-doping C-V calibration (calibrate_graded_doping, its 2D twin;
W targets at -10/-30 V move 4 % / 0.9 %), N_t fit of the dark current (-60 V +0.3 %, more at low bias),
and decide whether the alpha deposits in the n+ (CCE -3.7 %).

## n+ substrate in production (feat/nplus-substrate, 2026-10-16)

Phase 1 (balanced extraction, no n+): see commit message; lab/solver2d/bench/balanced_check.py.
Direct-vs-balanced dark-current gap is an N_t-independent ~3e-13 A/cm^2 round-off floor of the
direct form (balanced current exactly linear in N_t), not the SRV term as the spike note said.

Phase 2 (option `nplus_thickness_cm` / `N_D_nplus` / `nplus_generation` in create_sic_device,
create_sic_2d_device; DeviceConfig.nplus_thickness_um / N_D_nplus). lab/solver2d/bench/nplus_check.py:
- default (None) bit-for-bit vs phase 1 (bitforbit.py: mesh, W, C, CCE, dark, run_cv/run_cce/run_dark_current).
- reproduces the spike exactly: W(-10/-30/-50 V) 9.213316/9.895584/9.940143 um; CCE (epi-only G)
  0.998867 0.999562 0.999708 0.999778 0.999820 (raw collected/generated, < 1).
- wall: 1D fixed 0.5 V steps 67.5 V -> 150 V; 2D hw 25, epi 10 um, N_D_bulk 5e13: 150 V (161 s).
- dark current depended on the n+ thickness (f_depl floor of the effective-generation proxy
  generates throughout the n+: -60 V 1.8014e-11 A at 1 um, 1.8114e-11 at 5 um). The proxy is now
  confined to x <= epi end when an n+ layer exists: 1/5/20 um -> 1.7982941/1.7982946/1.7982946e-11 A
  (-60 V, 0.04 cm2, N_t 2.8e12; no n+: 1.7956753e-11, +0.15 %).
- alpha deposit in the n+ (nplus_generation=True; 15 um range, 10 um epi, I_gen 1.9224e-4 vs
  1.5919e-4 epi-only), CCE at -10 / -50 V by n+ thickness (n+ inherits epi tau_p 0.6 us, mu_p 75.7,
  L_p ~11 um): 1 um 0.9612/0.9622 (track truncated, I_gen 1.7259e-4), 5 um 0.9466/0.9479,
  10 um 0.9680/0.9694, 20 um 0.9758/0.9772, 50 um 0.9770/0.9784. Not converged in thickness below
  ~50 um and set by the assumed n+ hole lifetime; a realistic short n+ lifetime drives it to the
  no-collection bound CCE_epi * 1.5919/1.9224 = 0.828. Epi-only: 0.998868/0.999820 (any thickness).

Phase 3 (batched 1D engine with n+): lab/solver2d/bench/nplus_parity1d.py, 10 configs to -100 V:
max rel dW 6.1e-15, |dpsi| 1.4e-14 V, |dCCE| 3.2e-11. p+ 0.3/0.5 um without n+: dW 1.6e-15, dCCE 3.2e-12.
Engine box: p+ 0.3-3 um, n+ 1-5 um, N_D_nplus 1e18-1e19. 2D engine: rejects n+ (devsim).
Reach-through with n+ = edge within one Debye length of the epi end (12.85 V default device;
the mesh rule gives 17.0 V without n+ and never fires with n+).

Phase 4 (recalibration, record only: etna.core.device.NPLUS_CALIBRATION; lab/solver2d/bench/nplus_calibrate.py):
- C-V targets W = 1.7 / 9.5 / 9.73 um at 0 / -10 / -30 V (A = 4 mm2): NO exact solution for
  (N_D_junction, N_D_bulk, L) with n+ (best: N_D_bulk at the 1e12 bound, -30 V +1.5 %) nor without
  (best: N_D_junction at the 1e17 bound, L 0.11 um, -30 V +2.6 %; the shipped defaults also miss
  -30 V by +2.61 %). Past full depletion W saturates at the epi end (n+: 9.90 um at -30 V).
  p+ thickness has no effect on W (1 vs 0.3 um identical).
- 2-point exact solve (0, -10 V) at held L (n+, p+ 0.3 um), residual -30 V and deliverables:
  | L um | N_D_junction | N_D_bulk | res -30 V | CCE -10 / -50 V | I(N_t 2.8e12) -10/-30/-60 V (A, 4 mm2) |
  | 0.500 | 6.770553e15 | 9.228039e13 | +1.77 % | 0.999104 / 0.999819 | 6.140e-12 / 1.5216e-11 / 1.79667e-11 |
  | 0.750 | 3.948229e15 | 7.985437e13 | +1.76 % | 0.999070 / 0.999819 | 6.155e-12 / 1.5028e-11 / 1.79659e-11 |
  | 0.987 | 2.974271e15 | 6.335936e13 | +1.74 % | 0.999021 / 0.999819 | 6.121e-12 / 1.4845e-11 / 1.79657e-11 |
  | 1.250 | 2.451669e15 | 3.920880e13 | +1.71 % | 0.998954 / 0.999819 | 6.066e-12 / 1.4663e-11 / 1.79656e-11 |
  | 1.500 | 2.179422e15 | 4.732929e12 | +1.7 %  | 0.998897 / 0.999818 | 6.029e-12 / 1.4560e-11 / 1.79655e-11 |
  (L = 2 um: no 2-point solution, N_D_bulk at the 1e12 bound.) Record: L = 0.987 um (held).
- N_t: 2.8e12 -> 2.805349e12 (I(-60 V) 17.9657 -> 18.0000 pA; -10 V 6.132e-12, -30 V 1.4873e-11 A).
- E_t 0.65 -> 0.67 (radiation_damage Z1/2, Burin Table I): leakage <= 1.8e-6 (n+ model), but on the
  default no-n+ model run_dark_current then stops at -89.5 V instead of -100 V (wall), so the default
  SiC4H_Parameters.E_t stays 0.65 until the flip. run_dark_current reports no "truncated" flag (open item).
- p+ 1 -> 0.3 um (n+ model): dark -0.28 % at -10 V, -0.10 % at -60 V (f_depl floor in the p+);
  CCE <= 8.5e-7.
- Old default (no n+, p+ 1 um, 2.93e15/8.82e13/0.987, N_t 2.8e12) -> record: CCE -10/-30/-50 V
  0.963266/0.997719/0.998474 -> 0.999021/0.999706/0.999819; dark -10/-30/-60 V
  5.962/14.564/17.957 pA -> 6.132/14.873/18.000 pA; W 0/-10/-30 V 1.697/9.515/9.984 -> 1.700/9.500/9.899 um.

Post-phase-4 checks:
- 2D devsim with n+ (hw 25, n+ 1 um): balanced raw CCE 0.9995626 (-20 V) / 0.9998488 (-60 V),
  = 1D, laterally flat; the direct contact current gives 1.0136 / 0.9891 (unusable). Alpha in the
  n+ (1 um): 0.9619 / 0.9622. Fixed: cce_lateral_scan's in_mesh_generation_fraction now counts the
  n+ depth when nplus_generation=True (0.8978 vs 0.8280 epi-only).
- Alpha-in-n+ bracket on the record (p+ 0.3 um, record doping), raw CCE at -10 / -30 / -50 V:
  epi-only (any n+ thickness) 0.999022 / 0.999706 / 0.999819;
  in n+, 5 um 0.946929 / 0.947792 / 0.947922; 20 um 0.976105 / 0.977040 / 0.977178;
  50 um 0.977322 / 0.978260 / 0.978398 (n+ inherits epi tau_p 0.6 us -> L_p ~11 um);
  no collection from the n+: 0.999819 * 1.591911/1.922426 = 0.8279 (-50 V).

## 2D engine with n+ (feat/engine2d-nplus, 2026-10-17)

### Experiment: engine internals on the n+ mesh vs devsim
`uv run python lab/solver2d/bench/nplus2d_parity.py START STOP [csv]` (3 chunks in parallel, load ~6)
-> notes/nplus2d_parity.csv. Engine side = Device2D on grid_from_devsim(..., nplus_thickness_cm,
N_D_nplus, nplus_generation) (public gate bypassed), generation as cce_lateral_scan
(confine_generation_to_epi unless nplus_generation); devsim side = create_2d_dd_device +
dark dump + cce_lateral_scan (balanced extract_steady_state_current), 10 lateral points.
No engine physics change was needed: doping and mesh come from devsim's evaluated node
models, the ohmic BCs from NetDoping at the contact rows (n+ at the cathode), and the
engine's cathode current was already the balanced form (electrons = interior volume
balance + anode minority flux, holes = direct cathode minority flux) -- the same identity
as extract_steady_state_current.

21 configs: hw 25/50/100/150; epi 5/7.5/10/15/20 um; N_D_bulk 5e13/8.8e13/5e14; n+ 1/2/3/5 um;
N_D_nplus 1e18/3e18/1e19; V 5-100; nplus_generation False (16) / True (5).
| quantity | max over the 21 configs |
|---|---|
| W (column x=0, n/Donors >= 0.5 rule) | 2.8e-14 um (dW 0 in 17/21) |
| Potential | 1.4e-13 V |
| log10 n, log10 p (nodes with n, p > 1 cm^-3 in both) | 7.5e-13, 8.3e-13 |
| CCE, every lateral point | 1.2e-12 (partially depleted 20 um/5e14/10 V); <= 3.4e-13 elsewhere |
| Gummel iterations per stripe | <= 4 |
devsim converged at every point (no wall with n+ up to 100 V). Engine times in the CSV are
from a loaded machine incl. per-process JIT and are NOT speed-up evidence.

### Promotion: gate, sweep and run_field wiring
- etna.core.engine2d: ENGINE2D_NPLUS_BOX = n+ 1-5 um, N_D_nplus 1e18-1e19, V_bias 5-100 V
  (with n+; without n+ the box and the field rule are unchanged: 5-80 V). Outside it ->
  Engine2DFallback (devsim). With n+ the cathode-cell field rule (FIELD_ENGINE_MAX /
  FIELD_NAN_MIN, calibrated on the contact-on-epi wall) is not applied and
  Engine2DOutsideValidity is never raised: the last cell is inside the n+ and devsim
  has no wall there. dark_state_2d / engine2d_dark_unsupported_reason: same n+ box,
  V_cathode 5-100, nplus_generation rejected (meaningless for the dark state).
  Lateral uniformity is still asserted on every built mesh (grid_from_arrays).
- _sweep_one task tuple takes an optional 9th element (extra create_sic_2d_device
  keywords, forwarded to the gate, the engine and create_2d_dd_device);
  microdosimetric_sweep(nplus_thickness_cm=None, N_D_nplus=None, nplus_generation=False).
  A devsim ramp failure on an n+ device is reported as a plain error, not as the
  contact-on-epi "outside model validity" row.
- run_all sweep: no longer skips n+ configs (engine box gate as before; summary states the
  n+ layer or "no n+ substrate"); run_field with n+ goes to the engine via device_kwargs.

Production-path parity, `bench/nplus2d_sweep_parity.py sweep 0 5` / `5 10` (2 in parallel)
-> notes/nplus2d_sweep_parity.csv, _sweep_one auto vs devsim, 10 corner configs (hw 25-150,
epi 5-20, 5-100 V, n+ 1/2/3/5 um, 1e18/1e19, 5 with nplus_generation): every row solver
engine2d; max |dCCE| centre 4.2e-14, edge 9.0e-14, |dratio| 5.2e-14; is_valid identical
(10/10); no NaN / failure rows on either side.
`bench/nplus2d_sweep_parity.py field` (run_field engine vs forced devsim, same process):
| hw | n+ | bias | engine | devsim | x/y order, NetDoping | max dPotential | max rel dE (floor 1e-3 Emax) |
|---|---|---|---|---|---|---|---|
| 25 | 1 um 1e19 | -20 | 0.29 s | 5.65 s | identical | 5.3e-15 V | 4.0e-12 |
| 25 | 1 um 1e19 | -100 | 0.21 s | 13.6 s | identical | 1.4e-14 V | 7.7e-12 |
| 50 | 5 um 1e18 | -60 | 0.46 s | 37.5 s | identical | 2.8e-14 V | 1.9e-11 |
| 150 | 2 um 1e19 | -100 | 1.14 s | 160.7 s | identical | 5.7e-14 V | 1.4e-11 |
(rel dE is 10-50x the no-n+ value 3-4e-13: |E| is the LinearNDInterpolator post-processing,
which sees the 1 nm cells at the n-/n+ step; potential agrees to 1e-14 V.)

Timing (`bench/nplus2d_bench.py`, serial, JIT warm, load 1.7-2.3; epi 10 um, 8.5e13,
n+ 1 um 1e19, 60 V, 10 lateral points; engine best of 3, includes the mesh-only devsim build):
| hw | devsim build+ramp | devsim total | engine | speed-up | max dCCE |
|---|---|---|---|---|---|
| 50 | 22.3 s | 28.8 s | 0.215 s | 134x | 8.1e-15 |
| 150 | 84.8 s | 108.6 s | 1.41 s | 77x | 5.0e-15 |

Default (no n+) unchanged: notebook-20 36-config auto sweep with HEAD code (git archive) vs
this branch: DataFrame.equals True (bit-identical).

run_field with n+ AND non-default graded doping (the dark gate admits N_D_junction 1e15-1e16,
L 0.3-3 um), `bench/nplus2d_sweep_parity.py field doping`:
| hw | doping (N_D_junction / L um / N_D_bulk / epi um) | n+ | bias | engine | devsim | dPotential | rel dE |
|---|---|---|---|---|---|---|---|
| 25 | 1e16 / 0.3 / 5e13 / 10 | 5 um 1e18 | -100 | 0.30 s | 14.4 s | 2.8e-14 V | 9.9e-12 |
| 25 | 1e15 / 2.9 / 5e14 / 10 | 1 um 1e19 | -100 | 0.19 s | 14.4 s | 1.4e-14 V | 9.0e-12 |
| 50 | 5e15 / 2.0 / 1e14 / 20 | 3 um 3e18 | -5 | 0.65 s | 16.0 s | 1.8e-15 V | 1.1e-12 |
| 25 | 1e16 / 1.0 / 2e14 / 5 | 2 um 1e19 | -100 | 0.16 s | 10.6 s | 2.8e-14 V | 1.0e-11 |
x/y order and NetDoping identical in every case. (Pre-existing edge, not n+: L_transition_um=3.0
is rejected by the dark gate because 3.0*1e-4 > 3e-4 in floating point -> devsim.)

## Stack box: p+ 0.3-1 um and graded doping in the 2D engine gate (feat/nplus-default, 2026-10-17)

Probe: the recalibrated stack (NPLUS_CALIBRATION: p+ 0.3 um, n+ 1 um 1e19, 2.974271e15 /
6.335936e13 / L 0.987 um) was rejected by both 2D gates: `_FIXED_KWARGS` / `_DARK_FIXED` pinned
substrate_thickness_cm = 1e-4 and (scan gate) N_D_junction / L_transition to device2d's
2.931417e15 / 9.872267e-5 by exact equality, so even the old DeviceConfig doping (2.93e15 /
0.987 um) went to devsim. No engine change needed (mesh and doping come from devsim's mesh-only
build; the junction position is info["substrate_thickness_cm"]).

`uv run python lab/solver2d/bench/stack2d_parity.py parity A B` (3 chunks in parallel, load ~4.5)
-> notes/stack2d_parity.csv. 16 configs: hw 25/50/100/150; epi 5/7.5/10/15/20 um; p+ 0.3/0.5/0.7/1 um;
N_D_junction 1e15-1e16, L 0.3-3 um, N_D_bulk 5e13-5e14 (calibrated, legacy 2.93e15/8.82e13 and
corners); with n+ (10, 2 with nplus_generation) and without (6); 10-100 V.
| quantity | max over the 16 configs |
|---|---|
| W (column x=0) | 2.1e-14 um |
| Potential | 1.4e-13 V |
| log10 n / log10 p | 1.7e-11 (no-n+ rows near the cathode, n ~ 1e2) / 7.2e-13 |
| CCE, every lateral point | 6.3e-12 (partially depleted 20 um / 5e13 / 10 V); <= 3.6e-14 elsewhere |
| Gummel iterations per stripe | <= 5 |

Field rule off the old doping (`stack2d_parity.py wall`, devsim 0.5 V ramp, no n+, hw 25, epi 10):
- p+ 0.3, calibrated doping: last converged 55.0 V field 4.666e4, failed 55.5 V field 4.716e4.
- p+ 1, legacy 2.93e15 / 8.82e13 / 0.987: last converged 56.0 V 4.549e4, failed 56.5 V 4.598e4.
The second failure is below the old "every failure above 4.685e4" but above FIELD_ENGINE_MAX
(4.45e4): devsim decides in that band, so the rule (engine below 4.45e4, NaN above 5e4) holds.

Gate: ENGINE2D_BOX += substrate_thickness_cm 0.3e-4-1e-4, N_D_junction 1e15-1e16, L_transition
0.3e-4-3e-4 (= the dark box ranges); ENGINE2D_DARK_BOX += substrate_thickness_cm 0.3e-4-1e-4;
None = create_sic_2d_device default. microdosimetric_sweep takes substrate_thickness_cm /
N_D_junction / L_transition (forwarded with the n+ keys); run_all's sweep section forwards the
config's p+ / doping (it used to drop them: the 2D sweep ran on device2d's doping whatever the
sidebar said) and gates on the config's doping_profile / N_A / T (a 310 K config used to be
swept at 300 K under its label).

Production path (`stack2d_parity.py sweep`, `_sweep_one` auto vs devsim, 6 configs incl. the
calibrated stack with/without n+, legacy, two corners): every row engine2d; |dCCE| centre
<= 1.4e-14, edge <= 1.2e-14, |dratio| <= 1.6e-14; is_valid / is_fully_depleted / V_fd identical.
Calibrated stack V_fd (analytic graded) = 8.208 V (legacy doping: 10.40 V).
`stack2d_parity.py field` (run_field engine vs forced devsim):
| hw | stack | bias | engine | devsim | x/y, NetDoping | dPotential | rel dE |
|---|---|---|---|---|---|---|---|
| 25 | calibrated + n+ | -20 | 0.53 s | 5.7 s | identical | 4.4e-15 V | 5.1e-12 |
| 50 | calibrated + n+ | -50 | 0.34 s | 23.5 s | identical | 2.1e-14 V | 1.4e-11 |
| 25 | calibrated + n+ | -100 | 0.21 s | 13.3 s | identical | 1.4e-14 V | 9.7e-12 |
| 25 | calibrated, no n+ | -30 | 0.15 s | 3.1 s | identical | 5.3e-15 V | 3.8e-14 |
| 50 | legacy, no n+ | -20 | 0.21 s | 7.0 s | identical | 5.3e-15 V | 4.5e-13 |
