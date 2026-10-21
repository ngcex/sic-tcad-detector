# Physical-model library (issue #6): log

Started 2026-10-21. Code:

- `etna/core/physics_models.py`: the registry.
- `etna/core/avalanche.py`: ionization integrals and breakdown.
- `etna/core/batched/models_dd.py`: the coupled avalanche / high-field DD.
- `etna/core/devsim_models.py`: the devsim twin.
- `etna/api/models.py`: the facades.

Tests: `tests/test_physics_models.py` (24 tests, ~7 s) and `tests/test_physics_models_devsim.py` (4 tests, ~10 s). Docs: `docs/USAGE.md` §18.

## Reproduce

```bash
uv run python lab/physics_models/bench/models_bench.py --devsim   # ~21 s on an M4 (10 cores)
```

## 1. Research

Parameters are only entered from sources that were opened.

**Primary:**

- Hatakeyama et al., APL 85, 1380 (2004), Table I. Open repository copy: https://tsukuba.repo.nii.ac.jp/record/18553/files/APL_85-8.pdf
- Hatakeyama, SISPAD 2005, Eq. 1. http://in4.iue.tuwien.ac.at/pdfs/sispad2005/1562052.pdf
- Ayalew thesis (TU Wien 2004), nodes 61-75. https://www.iue.tuwien.ac.at/phd/ayalew/

**Secondary:**

- Burin et al., arXiv:2410.06798 v8: Tables XIII, XXVII, XXVIII, XXXIII and XXXIV.
- Kyuregyan, arXiv:1506.00980. Agrees with the review on Loh 2008, Ng 2003 and Niwa 2014.

**Not opened** (paywalled): Konstantinov 1997, Niwa 2015, Zhao 2019, Loh 2008, Roschke & Schwierz 2001, Lindefelt 1998, Khan & Cooper 2000. Their numbers carry status `verified-secondary`, with no DOI recorded.

**Findings on the existing ETNA physics.** These are documented in the issue; the defaults are unchanged.

- **Caughey-Thomas direction.** The set is the Schaffer 1994 basal-plane set, perpendicular to c, but it is used for vertical (c-axis) transport.
- **T exponents.**
  - The electron gamma -2.40 is Schaffer's <0001> value.
  - The hole gamma -2.15 is borrowed from the electrons.
  - Ayalew also has (T/300)^-0.5 on mu_min. ETNA has none; the `ayalew2004` low-field set adds it.
- **Ayalew's 4H Chynoweth set** carries the 6H b_n with a factor-10 change (review App. A4). It is registered as `ayalew2004_legacy`, and it overestimates V_BD by +48 to +94 %.

## 2. Results (bench, 2026-10-21, M4)

### (1) V_BD vs Hatakeyama 2004 measured, 1940 (1e16/N_D)^0.8 V

Modelled device: p+ 0.3 um 1e19 / 10 um n / n+ 1 um 1e19.

| N_D                                                  | 3e16    | 5e16    | 1e17    | 2e17    |
| ---------------------------------------------------- | ------- | ------- | ------- | ------- |
| measured fit (V)                                     | 805.6   | 535.3   | 307.5   | 176.6   |
| hatakeyama2004 (consistency: fitted on these diodes) | -0.6 %  | -0.3 %  | +1.8 %  | +7.0 %  |
| niwa2015                                             | -6.4 %  | -2.7 %  | +5.2 %  | +18.1 % |
| zhao2019                                             | -10.8 % | -8.1 %  | -2.1 %  | +8.0 %  |
| loh2008                                              | -17.8 % | -14.6 % | -8.2 %  | +2.0 %  |
| stefanakis2021                                       | -2.1 %  | +4.3 %  | +16.2 % | +33.5 % |
| ayalew2004_legacy                                    | +47.8 % | +55.7 % | +70.6 % | +93.9 % |

- **Peak field vs F_C = 2.70 (1e16/N_D)^-0.1 MV/cm:** hatakeyama2004 -1.6 / -2.1 / -2.6 / -3.2 %.
- **Timing:** about 0.1-0.2 s per set for 4 devices batched, about 80-130 batched Poisson solves (the first call includes compilation).
- **Open:** the +7 % at 2e17 (hatakeyama2004). Candidate causes:
  - the 0.8 power-law fit itself;
  - the p+ depletion (N_A- 1.3e18 with the legacy ionized fraction);
  - the 2-5 MV/cm fit range edge (F_max 3.53 MV/cm).

### (2) Default dosimeter (DeviceConfig()), 1D parallel-plane V_BD, mesh refinement (h, h/2, h/4)

| Set               | V_BD (V) | h/2     | h/4     | F_max (MV/cm) |
| ----------------- | -------- | ------- | ------- | ------------- |
| hatakeyama2004    | 2225.78  | 2225.89 | 2225.93 | 2.290         |
| niwa2015          | 1955.44  | 1955.52 | 1955.57 | 2.022         |
| zhao2019          | 1928.18  | 1928.25 | 1928.29 | 1.995         |
| loh2008           | 1825.64  | 1825.69 | 1825.72 | 1.893         |
| stefanakis2021    | 1916.93  | 1916.99 | 1917.04 | 1.984         |
| ayalew2004_legacy | 2489.87  | 2489.94 | 2489.98 | 2.551         |

- **Convergence** is about first order: the |dV| ratio over h -> h/2 -> h/4 is 1.5-2.3.
- **The 2h quadrature-only estimate** (0.001-0.003 V) underestimates the true mesh error, because it ignores the Poisson discretisation. The facade therefore reports the real re-solve, `mesh_error_V`.
- **Facade cost:** `run_breakdown(DeviceConfig())` takes 0.80 s, including the mesh-check re-solve.

### (3) Temperature, default dosimeter V_BD (V)

| Set            | 250 K   | 300 K  | 350 K  | 400 K  | 448 K  |
| -------------- | ------- | ------ | ------ | ------ | ------ |
| niwa2015       | 1953.3* | 1955.4 | 1976.0 | 2005.5 | 2038.3 |
| zhao2019       | 1929.0  | 1928.2 | 1944.8 | 1970.0 | 1998.7 |
| stefanakis2021 | 1903.5* | 1916.9 | 1931.3 | 1946.2 | 1961.1 |

\* outside the set's T range (`allow_T_extrapolation`). zhao2019's dip from 250 to 300 K is the hole c, d issue the review flags.

### (4) Coupled block Newton vs lagged avalanche term

Default dosimeter, alpha generation, hatakeyama2004.

| V/V_BD | M McIntyre | gain DD | (1-1/M) rel | coupled     | lagged (sweeps, time, J diff) |
| ------ | ---------- | ------- | ----------- | ----------- | ----------------------------- |
| 0.900  | 2.151      | 2.151   | +6.6e-7     | 3x2, 3.5 ms | 15, 14 ms, -2.4e-9            |
| 0.970  | 9.809      | 9.809   | +3.9e-6     | 3x2, 21 ms  | 48, 55 ms, -2.2e-8            |
| 0.990  | 34.486     | 34.491  | +4.3e-6     | 3x2, 12 ms  | 136, 186 ms, -7.5e-8          |
| 0.997  | 121.953    | 122.017 | +4.3e-6     | 4x2, 6 ms   | 421, 514 ms, -2.6e-7          |

The lagged scheme's sweep count grows like M (ratio 1 - 1/M); the coupled scheme's does not.

### (5) Effects on the default dosimeter

**CCE (alpha):**

|                                     | -10 V   | -40 V   | -70 V   | -100 V  |
| ----------------------------------- | ------- | ------- | ------- | ------- |
| lades2000 minus default             | -5.0e-5 | -5.9e-5 | -6.2e-5 | -6.5e-5 |
| lades2000 + burin2024 minus default | -1.7e-4 | -1.2e-4 | -1.2e-4 | -1.2e-4 |

**Electron transit time** (full depleted width; ratio to low field in brackets):

|                       | -10 V            | -40 V           | -70 V           | -100 V          |
| --------------------- | ---------------- | --------------- | --------------- | --------------- |
| lades2000             | 156.5 ps (1.25x) | 72.0 ps (2.12x) | 57.8 ps (3.28x) | 52.6 ps (4.65x) |
| lades2000 + burin2024 | 149.2 ps (1.19x) | 81.7 ps (2.40x) | 70.5 ps (4.00x) | 66.4 ps (5.87x) |

**Hole transit time:**

|                       | -10 V           | -40 V          | -70 V          | -100 V         |
| --------------------- | --------------- | -------------- | -------------- | -------------- |
| lades2000             | 985 ps (1.05x)  | 312 ps (1.25x) | 199 ps (1.53x) | 157 ps (1.85x) |
| lades2000 + burin2024 | 1166 ps (1.25x) | 400 ps (1.60x) | 273 ps (2.09x) | 226 ps (2.67x) |

**With avalanche** (hatakeyama2004):

|                     | -100 V   | -575 V   | -1050 V | -1525 V | -2000 V |
| ------------------- | -------- | -------- | ------- | ------- | ------- |
| dark multiplication | 1        | 1        | 1.00001 | 1.01368 | 2.03299 |
| alpha CCE           | 0.999936 | 0.999989 | 1       | 1.01453 | 2.11515 |

### (6) devsim twin vs engine

Device: p+ / 3 um 1e17 / n+, G 1e18 in the epi, 300 V, hatakeyama2004 + lades2000.

- J: devsim -1.6390581572e-3, engine -1.6390581572e-3, relative difference -3.1e-12.
- Time: devsim 1.97 s (build, ramp and solve) against 13 ms for the engine (Poisson ramp and solve).

## 3. Design notes and decisions

**devsim twin.**

- devsim's symdiff treats a referenced model as a constant. The field-dependent mobility must therefore enter the current with its derivative.
  - Fully inlined, it took about 6 s per current in symdiff (CreateEdgeModelDerivatives).
  - The chain rule written out (J = F(E) J_low) takes about 0.1 s.
- The impact term is an `edge_volume_model` with explicit derivative models, d|J| = sgn(J) dJ.
- |E| in alpha is regularised as sqrt(E^2 + 1) so that the derivative does not overflow at E = 0.

**Current extraction.**

- The batched engine's balanced form, and devsim's `extract_steady_state_current`, sum node generation models. Neither sees the edge-volume impact term.
- The new solvers therefore read Jn + Jp on the edge with the smallest |SG terms|. With the models off this equals the balanced form to round-off: 0.9994291317960069 both ways, default device, -10 V.
- The peak-field edge is NOT suitable. It sits on the p+ junction, where holes are majority, and gave a 1.9e-4 cancellation error.

**Dark current with avalanche.** It runs without the Hurkx factor. The closed form exceeds its validity limit (1.1e6 V/cm) long before avalanche fields.

**Not coupled (open in issue #6):**

- BGN in Poisson.
- 2D anisotropy and anisotropic impact ionization.
- Field-dependent trapping.
- Transient runs with the models.
