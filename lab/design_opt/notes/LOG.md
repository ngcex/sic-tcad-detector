# Design optimisation on the fast engines (issue #10)

Bench: `uv run python lab/design_opt/bench/design_bench.py [throughput analytic dosimeter robust micro]`.
Date 2026-10-21, Apple M4 (10 cores) shared with ~12 agents. **The load average was 38-53
during every run below**, so absolute times are 2-5x pessimistic. Earlier quiet-ish
measurements are given in brackets.

## 1. Evaluator parity (tests/test_design_device.py)

`etna.design.evaluate_configs` drives the core kernels directly:

- `reverse_sweep`, `cce_sweep` twice (the alpha profile and a uniform generation);
- `dark_sweep` with the batch doping and `tat_constants`;
- `breakdown_voltage`;
- `run_radiation_damage` per config.

Agreement with the facades on 3 configs at -10 V:

| quantity                 | facade                 | difference          |
| ------------------------ | ---------------------- | ------------------- |
| C, W                     | `run_cv`               | 0                   |
| CCE                      | `run_cce`              | <= 1.6e-15          |
| dark current             | `run_dark_current`     | <= 6.9e-15 relative |
| transit time (lades2000) | `run_transit_time`     | 0                   |
| V_BD                     | `run_breakdown`        | <= 7.7e-12 V        |
| CCE at fluence           | `run_radiation_damage` | 0                   |

Outside `_ENGINE_BOX` the point is returned `ok=False` ("no devsim"). The test asserts this with the devsim paths monkeypatched to raise.

## 2. Throughput (cheap metric set, warm numba cache)

**One shared mesh** (epi 10 um):

| B   | per device                                                            |
| --- | --------------------------------------------------------------------- |
| 1   | 58.7 ms [41 ms at lower load]                                         |
| 8   | 8.3 ms [13 ms]                                                        |
| 64  | 10.9 ms [7.1 ms; profile 5.8 ms, of which the engine is about 3.2 ms] |

**Mesh grouping.** Of 64 designs with epi uniform in 8-12 um:

| epi values                     | mesh groups | per device |
| ------------------------------ | ----------- | ---------- |
| continuous                     | 64          | 56 ms      |
| snapped to the 0.25 um lattice | 17          | 27 ms      |

The lattice is what makes batching pay.

**Costly metrics**, extra per device at B = 8:

- `cce_fluence`: +322 ms under load [+55 ms]. It goes through `run_radiation_damage`, which builds a devsim mesh-and-doping reference per call (no solve) and is not batched across configs. This is the bottleneck of the default problem.
- `V_bd`: +55 ms. Breakdown is a constraint-free post-check of the front by default; it never binds in this box (> 1.3 kV against 10-100 V).

**Default problem, 60 designs:**

- At -10 V: 44 s of evaluation plus 24 s of GP/acquisition time, so 1.35 designs/s under load.
- Robust worst case (33 tolerance samples per design for the cheap metrics, 2 for `cce_fluence`): 0.44 designs/s. p5: 0.36 designs/s.

**Fix found on the way.** On a cache hit, `throughput_designs_per_s` was the hit count divided by about 0 s; it now counts cache misses only.

## 3. Sample efficiency on analytic problems (5 seeds, median [min, max])

| problem (budget)                        | BO                      | Sobol, same budget    |
| --------------------------------------- | ----------------------- | --------------------- |
| Branin, gap to 0.397887 (36)            | 1.2e-4 [1.3e-5, 3.0e-4] | 0.26 [0.058, 0.92]    |
| constrained quadratic, gap to 0.18 (24) | 2.8e-5 [8e-7, 9.8e-5]   | 0.054 [0.0049, 0.093] |
| ZDT1 d=3, HV / true HV (50)             | 0.91 [0.87, 0.94]       | 0.10 [0.025, 0.32]    |

GP and acquisition cost about 3 s per run under load (0.4-1.9 s quiet).

## 4. Default dosimeter problem (-10 V; max sensitivity vs max CCE after 1e14 p/cm2 at 62 MeV)

**BO against random search.** At 60 designs, with HV taken with reference (S = 0, CCE = 0.8):

| seed | BO   | Sobol |
| ---- | ---- | ----- |
| 0    | 4.83 | 3.40  |
| 1    | 4.88 | 3.31  |
| 2    | 4.95 | 3.05  |

**The front is set almost entirely by the epi thickness.** Sensitivity is close to proportional to the epi: 7.45 nC/Gy at 5.5 um, 12.9 at 9.5 um, 34.1 at 25.25 um. CCE(1e14) goes the other way: 0.992 at 5.5 um, 0.975 at 9.5 um, 0.877 at 25.25 um.

**The doping is driven by the gate.** The binding constraint is the full-depletion gate (1.2 V_fd <= 10 V):

- thick designs push N_D_bulk to the box floor, 2e13;
- thin designs sit at about 8.5e13-1e14 for radiation hardness.

**Dark current** is 3.7-4.6 pA, against a budget of 18 pA per 4 mm2: it never binds at -10 V.

**At -50 V** the gate is loose:

- the front moves to epi 26-30 um (S 36-41 nC/Gy, CCE(1e14) 0.967-0.971);
- leakage is 13.5-17.9 pA on 4 mm2, so it approaches the per-area budget;
- the HV trace still rises at 60 designs (0.497 -> 0.671): not converged.

## 5. Robust mode (DESIGN-1 tolerances, 16 Sobol + 16 corners + nominal)

| mode    | max epi on the front | S range (nC/Gy) | CCE(phi) range | gate                             |
| ------- | -------------------- | --------------- | -------------- | -------------------------------- |
| nominal | 25.25 um             | 7.45 - 34.1     | 0.877 - 0.992  |                                  |
| worst   | 14.0 um              | 6.78 - 17.6     | 0.9445 - 0.991 | many sit at 0.97-0.99 worst case |
| p5      | 17.5 um              | 7.05 - 22.6     |                |                                  |

The nominal optimum is **not robust**. Its thick-epi end fails the gate once the epi grows by 1 um and the doping rises by 20 %. The worst-case front keeps N_D_junction near the robust floor (1.25e15) and L at its robust floor (0.6 um).

## 6. How far the optimum rests on unvalidated physics (front of section 4, re-evaluated)

**Legacy removal model with the si-neq convention.** 0 % of the front is solvable: 1e14 p/cm2 is beyond full compensation for every design (Phi_neq 1.94e14). The radiation-hardness axis therefore exists ONLY in the calibrated kinetic-trap model, in its semi-insulating regime, which has not been checked against a measurement. A design rule taken from this axis is a model prediction.

**si-neq convention (trap model).** CCE(1e14) drops by 0.02-0.03 on thick designs (0.877 -> 0.855) and by about 0.002 on thin ones. 100 % of the front stays feasible and non-dominated, so the ranking is unchanged.

**`legacy_perp` mobility.** CCE moves by at most 0.015. The front stays 100 % non-dominated at -10 V; at -50 V, 88 % does (it is flat there).

**Sensitivity** depends on NIST stopping powers and the low-injection uniform-generation CCE. Its proportionality to the epi is geometric and robust. Its absolute value (13.5 nC/Gy for the default 4 mm2 stack at 62 MeV) is not checked against a measurement.

**Graded transition.** Its thickness and end doping were fitted at 10 um and are applied at every epi: an extrapolation, flagged in the result notes.

## 7. Microdosimeter cell (fv2d, -30 V, 24 designs, 99 s under load: 0.25 designs/s)

**The front.** The efficient (CCE >= 0.95) area fraction rises with the cell radius: 0.57 at 15 um and 0.77 at 42.5 um.

**What sets the minimum.** The minimum CCE, 0.81-0.83, is the n+ column region in every design: generation inside the column is lost. The column radius and the epi barely matter.

**Resolution.** The efficient-area boundary is resolved only to the scan spacing (2-6 um with 8 points).

**Validity.** Annular-trench model only, not a column lattice (USAGE s.11.1).
