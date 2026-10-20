# ETNA usage guide

This is the single usage document for ETNA (Epitaxial Tcad Numerical Analysis). It is
written for a physicist or engineer who knows semiconductor radiation detectors but has
not read this code. Every Python block below was run against this repository (version
5.0.0, 2026-10-18); the numbers quoted come from those runs on an Apple M4 (10 cores) and
will differ slightly on other machines. File paths are relative to the repository root.

Contents

1. [What ETNA is and is not](#1-what-etna-is-and-is-not)
2. [What to trust (one-page table)](#2-what-to-trust)
3. [Installation](#3-installation)
4. [Quick start](#4-quick-start)
5. [The device stack, presets, units and sign conventions](#5-the-device-stack-presets-units-and-sign-conventions)
6. [The `run_*` facades](#6-the-run_-facades)
7. [Solver provenance and fast paths](#7-solver-provenance-and-fast-paths)
8. [Validity and limits](#8-validity-and-limits)
9. [The Streamlit app](#9-the-streamlit-app)
10. [Exports (Python)](#10-exports-python)
11. [Microdosimetric sweep](#11-microdosimetric-sweep)
12. [Notebooks](#12-notebooks)
13. [`lab/` and how provenance is recorded](#13-lab-and-how-provenance-is-recorded)
14. [Tests](#14-tests)
15. [Development workflow](#15-development-workflow)
16. [Troubleshooting](#16-troubleshooting)

---

## 1. What ETNA is and is not

ETNA is a 1D/2D drift-diffusion TCAD library for **4H-SiC p+/n-/n+ radiation detectors**
(the INFN-LNS Petringa group's 4 mm² proton dosimeter and its 25 mm² microdosimeter
sibling). It solves Poisson plus electron/hole continuity on a finite-volume mesh with
4H-SiC material models, and derives from it C-V, depletion and field profiles,
charge-collection efficiency (CCE) vs bias, an effective dark-current budget, CCE vs proton
fluence, a temperature sweep and a FLASH dose-rate sensitivity bound. On top of the library
sit a Streamlit UI (`app/`), engineering exports (`etna/export/`), optional manim
animations (`etna/viz/`) and 21 analysis notebooks.

Three solvers sit behind one API: [devsim](https://devsim.org/) (the reference finite-volume
solver), a batched numba 1D engine (`etna/core/batched/`) and a numba 2D engine
(`etna/core/engine2d/`). The engines reproduce devsim inside parameter boxes where they were
checked against it. Section 7 explains how each result tells you which solver produced it.

What ETNA is **not**:

- It is not a commercial TCAD package. There is no process simulation, no 3D, no edge or
  guard-ring physics on the default path, no oxide charge, and no impact ionisation or
  breakdown model.
- It is not validated against measurement in general. The distinction matters throughout
  this guide:
  - **Verified against devsim** (code-to-code). Inside documented boxes, the numba
    engines agree with devsim to dC/C ~3e-8 (devsim's own central-difference
    capacitance), dCCE <= ~1e-9 and dpsi ~1e-13 V. This shows the engines solve the same
    equations. It says nothing about whether the equations describe the real device.
  - **Calibrated**. The default doping profile is a fit to one digitised measured C-V
    curve (Photons paper, Fig. 4b). The dark-current trap rate `N_t` is fitted to the
    measured bound of < 18 pA at -60 V. Agreement with the data used for a fit is fit
    quality, not validation.
  - **Literature check**. A C-V comparison with Bruzzi & Verroi 2023 (a Schottky diode on
    a different 30 µm wafer) checks the solver scale only: +3.1 % at -40 V, < 1 % with
    V_bi matched (`.planning/VALIDATION_LITERATURE.md`). It does not validate this device.
  - **Not reproduced / open**. The microdosimeter's measured alpha CCE-vs-bias shape is
    not reproduced by any drift-diffusion variant (section 8). FLASH plasma recombination
    is not implemented. The tissue-equivalence (kappa) table is a placeholder; the NIEL
    hardness factors are SR-NIEL data in the Si 1-MeV-neq convention, with the SiC
    neutron normalisation unavailable (`data/niel/README.md`).
    SIMS is unavailable, so the epi/n+ position comes from the C-V fit.

---

## 2. What to trust

| Quantity                                                  | Status                                                      | Basis                                                                                                                                                                   | Where                                                                                                                       |
| --------------------------------------------------------- | ----------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| Engine results (C, W, CCE, potential) vs devsim           | **verified** (code-to-code)                                 | 1D: dC/C <= ~3e-8, dW <= 1e-17 cm, dCCE <= ~1e-9 (9e-10 at the tabulated-box corners) to -100 V; 2D: dpsi ~1e-13 V                                                      | `tests/test_batched_engine.py`, `tests/test_api_engine_switch.py`, `tests/test_engine2d*.py`, `lab/engine_box/notes/LOG.md` |
| C-V / eps*A/C of the default device                       | **calibration** (fit)                                       | default doping fitted to the digitised Photons Fig. 4b: plateau rms 0.013 µm, knee rms 0.24 µm; 1.655 / 9.519 / 9.750 µm at 0 / -10 / -30 V vs paper 1.68 / 9.52 / 9.73 | `lab/cv_gap/notes/LOG.md`, `etna/core/device.py` (`DEFAULT_STACK`)                                                          |
| Doping profile N_D(x)                                     | **inferred**, not measured                                  | inversion of one C-V curve; bulk 3.4e13 is below the documented 0.5-1e14; no SIMS                                                                                       | section 8                                                                                                                   |
| Full-depletion voltage                                    | **definition-dependent**                                    | analytic V_fd 2.74 V, reach-through 3.10 V, sweep gate 1.2 V_fd = 3.28 V (default stack)                                                                                | `etna.export.labels.VFD_DEFINITION`                                                                                         |
| CCE vs bias (1D, alpha in epi)                            | **model** (low injection, steady state)                     | engine verified vs devsim; no measured CCE reproduced for the dosimeter                                                                                                 | `run_cce` docstring                                                                                                         |
| Microdosimeter alpha CCE shape (Fig. 4b, Phys. Med. 2025) | **not reproduced**                                          | flat to 0.5 % from 0 V in the model; data drop to 0.908 at 0 V                                                                                                          | `run_cce_tilted_alpha` docstring, `lab/microdosimeter/notes/LOG.md` §9                                                      |
| Dark current                                              | **calibration** (budget estimate)                           | one N_t fitted to < 18 pA at -60 V, 4 mm²; bias shape f_depl heuristic x Hurkx (1+Gamma), temperature shape assumed                                                                                         | `run_dark_current`, DESIGN-1                                                                                                |
| CCE vs proton fluence                                     | **relative shape only**                                     | NIEL hardness factors: SR-NIEL data, Si-neq convention (SiC normalisation open); points at or beyond the compensation fluence are NaN                                                                            | `run_radiation_damage` docstring                                                                                            |
| FLASH dose-rate CCE                                       | **numerical sensitivity bound**                             | no plasma or high-injection physics                                                                                                                                     | `run_flash_recombination`, `run_transient` docstrings                                                                       |
| Tissue-equivalence kappa(E)                               | **data-blocked**                                            | `data/srim/*.csv` are placeholders                                                                                                                                      | `data/srim/README.md`                                                                                                       |
| 2D field maps                                             | **verified** vs devsim inside the engine2d box              | laterally uniform device; no edge termination                                                                                                                           | `tests/test_engine2d*.py`, `tests/test_field_engine_default_100V.py`                                                        |
| Reference layout (GDS/SVG)                                | **reference only, not TCAD-validated, not for fabrication** | dimensions from `data/foundry/layout_params.toml`; ring, pad, dicing, oxide and 5x5 die are PLACEHOLDER drawing values until the mask data is entered                   | `deliverables/DIMENSIONS-provenance.md`, section 10                                                                         |

Summary: the numerics are well verified. The physics is calibrated on one C-V curve and
one leakage bound. Everything else is a model prediction that has not been checked
against a measurement of this device.

---

## 3. Installation

Requirements: Python >= 3.13 and [uv](https://docs.astral.sh/uv/). devsim ships as a
wheel; it looks for BLAS/LAPACK at import (on macOS it uses the system Accelerate/LAPACK;
on Linux install OpenBLAS/LAPACK, e.g. `apt-get install libopenblas-dev liblapack-dev`).

```bash
git clone git@github.com:ngcex/sic-tcad-detector.git
cd sic-tcad-detector
uv sync --extra dev --extra viz     # library + dev tools + optional animation/GDS
```

- `dev` is the one development dependency set (`pyproject.toml`): pytest,
  pytest-xdist, jupyter, nbstripout, ruff and meshio (independent VTK/VTU reader in the
  export tests).
- `viz` (optional) adds `manim` (Export page animations) and `gdstk` (GDSII download).
  On macOS manim's pycairo needs `brew install pkgconf cairo pango` first. Without
  `viz`, the Export page hides the animation section and the GDS button, and the
  matching tests skip themselves.
- **Always sync with `--extra dev`.** A plain `uv sync` installs no extras and removes
  pytest from `.venv`. `uv run pytest` can then fall through to a system pytest with the
  wrong interpreter (audit C-15).
- After cloning, run `uv run nbstripout --install --attributes .gitattributes` once, so
  that git strips notebook outputs on commit (section 12).

Check the install:

```python
import importlib.metadata
import etna

print(etna.__version__, importlib.metadata.version("etna"))  # 5.0.0 5.0.0
```

The version is defined once, in `etna/_version.py`; `pyproject.toml` reads it through
hatch (`dynamic = ["version"]`).

At import devsim prints a few lines about which BLAS/LAPACK it loaded (e.g.
`Loading "liblapack.dylib": ALL BLAS/LAPACK LOADED`). That is normal. Solver logs are
silenced (section 7).

---

## 4. Quick start

```python
from etna import DeviceConfig, run_cv, run_all

cfg = DeviceConfig()                      # the default 4 mm2 dosimeter (section 5)
cv = run_cv(cfg, v_start=0.0, v_stop=-30.0, n_points=31)

print(cv.sim_type, cv.x[:3], cv.y[:3])    # 'cv', bias in V, capacitance in F (for area_cm2)
print(cv.metadata["solver"])              # 'batched' (the 1D engine; section 7)
print(cv.metadata["validated_domain"])    # True: inside the engine/devsim cross-check box
print(cv.metadata["reach_through_V"])     # 3.1 (V, a positive magnitude)
print(cv.metadata["W_from_C"][[0, 10, 30]] * 1e4)  # eps*A/C in um at 0/-10/-30 V: ~1.655 9.519 9.750

report = run_all(cfg)                     # every section in one call, ~3 s on the default device
for row in report.rows():
    print(row["section"], row["status"], row["solver"], row["key figures"])
```

Every facade returns a `SimResult` (`etna/api/results.py`):

| field      | meaning                                                                                                                                                                                                          |
| ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `config`   | the `DeviceConfig` used                                                                                                                                                                                          |
| `sim_type` | `"cv"`, `"cce"`, `"field"`, `"damage"`, `"dark_current"` (x = bias), `"temperature"`, `"flash"`, `"transient"`, `"microdosimetry"`; the app's Dark Current page builds a `"dark_current_vs_T"` aggregate (x = T) |
| `x`, `y`   | primary axis and output (numpy arrays; see each facade in section 6)                                                                                                                                             |
| `metadata` | dict: solver provenance (section 7) plus facade-specific arrays and scalars                                                                                                                                      |
| `mesh`     | `MeshData` (node coordinates in cm, node values, regions, contacts) for `run_field`; `None` otherwise                                                                                                            |

`run_all` returns a `Report` (`etna/api/run_all.py`). It provides `report.ok`,
`report.rows()` (one summary dict per section), `report.get("cv").result` (the
`SimResult`) and `report.to_json()`. A failing section becomes a labelled row; it never
aborts the run.

---

## 5. The device stack, presets, units and sign conventions

### 5.1 The default stack

`DeviceConfig()` (`etna/api/device.py`) is the measured-C-V model, default since
2026-10-18 (`etna.core.device.DEFAULT_STACK`):

| layer                   | default                                                                                                                     | `DeviceConfig` field                                                 |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------- |
| p+ anode                | 0.3 µm, N_A 1e19 cm^-3 (anode contact at x = 0)                                                                             | `substrate_thickness_um` (**the p+ anode, despite the name**), `N_A` |
| n- epi, graded          | 10 µm; N_D(x) = N_D_bulk + (N_D_junction - N_D_bulk) exp(-x/L): 8.18588e15 at the junction, L 0.463863 µm, bulk 3.382506e13 | `epi_thickness_um`, `N_D_junction`, `N_D_bulk`, `L_transition_um`    |
| graded n-/n+ transition | the last 1.136634 µm of the epi; N_D rises log-linearly to 4.516421e16 at the epi end                                       | `buffer_thickness_um`, `buffer_end_doping`                           |
| n+                      | 1 µm modelled (the conductive substrate truncated), 1e19, cathode on its far end                                            | `nplus_thickness_um`, `N_D_nplus`                                    |
| area                    | **0.04 cm² = 4 mm²** (the dosimeter); C and I are absolute for this area                                                    | `area_cm2`                                                           |
| dark-current fit        | N_t 1.355575e12 cm^-3 s^-1 (Hurkx refit 2026-10-20), E_t 0.67 eV                                                            | `N_t`, `E_t`                                                         |
| temperature             | 300 K                                                                                                                       | `T`                                                                  |

`DeviceConfig` validates on construction: it rejects non-positive sizes, N_D >= N_A,
T outside 100-600 K and an inconsistent transition, with a `ValueError` that names the
field. Set `half_width_um` (µm) for a 2D device (lateral half-width; depth is the 2D
`y`).

### 5.2 Presets

| preset                            | call                                            | what it is                                                                                                                                                            | when to use                                                                                                                                                                                  |
| --------------------------------- | ----------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| default (measured-C-V parametric) | `DeviceConfig()` == `etna.measured_cv_config()` | the stack above: an exponential epi plus a graded n-/n+ transition, fitted to the digitised Photons Fig. 4b eps*A/C curve                                             | default for all dosimeter work                                                                                                                                                               |
| tabulated measured-C-V (opt-in)   | `etna.measured_cv_tabulated_config()`           | 19-knot log-linear N_D table from a model-independent inversion of the same curve (knee rms 0.014 µm vs 0.24 µm for the parametric form), N_t refitted to 1.356120e12 | sensitivity check of results that depend on the knee (low bias). It is still one C-V curve, not SIMS. The first 0.5 µm and 8.5-9.25 µm are unresolved (`lab/tabulated_profile/notes/LOG.md`) |
| previous default                  | `etna.nplus_calibrated_config()`                | 2026-10-17 model: exponential epi solved on 2 C-V points, abrupt 1e19 n+ step, N_t 1.343970e12 (Hurkx refit)                                                          | reproducing numbers of that period bit for bit. Its C-V knee (-8..-10 V) is excluded by the full digitised curve                                                                             |
| legacy                            | `etna.legacy_config()`                          | pre-n+ model: cathode directly on the n- epi, p+ 1 µm, 2.93e15 / 8.82e13 / 0.987 µm, N_t 2.8e12, E_t 0.65                                                             | reproducing old notebooks only. Its deep-bias devsim wall (~ -64 V) is a contact artefact                                                                                                    |
| microdosimeter                    | `etna.microdosimeter_config(p_plus_um=0.3)`     | Petringa 2025 Phys. Med. 140: 25 mm² (`area_cm2=0.25`), uniform N_D 1e14, abrupt n+, V_fd 9.0 V                                                                       | with `run_cce_tilted_alpha` for the 65° 241Am test; layout/geometry studies                                                                                                                  |

Each preset is also exposed as a dict (`etna.MEASURED_CV_PRESET`,
`etna.MEASURED_CV_TABULATED_PRESET`, `etna.NPLUS_CALIBRATED_PRESET`, `etna.LEGACY_PRESET`,
`etna.MICRODOSIMETER_PRESET`). Every preset function accepts overrides:

```python
import etna
from etna import DeviceConfig

assert etna.measured_cv_config() == DeviceConfig()
thin = DeviceConfig(epi_thickness_um=8.0)                      # one field changed
flat = DeviceConfig(doping_profile="uniform", N_D=1e14)       # uniform epi: abrupt epi/n+ unless set
old = etna.legacy_config(T=320.0)
md = etna.microdosimeter_config(p_plus_um=1.0)
twod = DeviceConfig(half_width_um=50.0)                        # 2D slice, 50 um half-width
print(md.area_cm2, md.N_D, flat.buffer_thickness_um)           # 0.25 1e14 None
```

### 5.3 Incomplete ionization (opt-in, `incomplete_ionization=True`)

By default (`incomplete_ionization=False`) ETNA keeps the **legacy** ionization: the p+
Al carries one scalar ionized fraction (0.132 at 1e19, 300 K,
`etna.core.incomplete_ionization.ionized_acceptor_concentration`) and the donors are
fully ionized. Every baseline and fit in this guide uses it. The scalar's low-doping
branch is known to be unphysical below ~1e18 (audit-v7 A-02: it rises with N_A);
`DeviceConfig` rejects N_A < 1e17.

`DeviceConfig(incomplete_ionization=True)` switches on the **carrier-dependent** model
(any preset accepts it as an override, e.g. `etna.measured_cv_tabulated_config(
incomplete_ionization=True)`). The ionized densities follow the local carriers inside
the Poisson equation (Burin et al., arXiv:2410.06798, Eq. 82; Pearson-Bardeen
dE(N) = dE0 - alpha N^(1/3), Eq. 85, clamped at 0):

| dopant        | levels (dE0, alpha)                                                           | g | source                                                |
| ------------- | ----------------------------------------------------------------------------- | - | ----------------------------------------------------- |
| N (donor)     | hexagonal 70.9 meV, 3.38e-5 meV cm; cubic 123.7 meV, 4.65e-5 meV cm (1/2 each) | 2 | Kagamihara et al., JAP 96, 5601 (2004); Burin Tab. XXV |
| Al (acceptor) | 220 meV, 2.32e-5 meV cm                                                       | 4 | Achatz et al., APL 92, 072103 (2008); Burin Tab. XXIV  |

N_C(T), N_V(T) come from `etna.core.sic_material.effective_dos`. Depleted dopants are
fully ionized and only neutral regions freeze out. At 300 K the neutral fractions are:
p+ Al 1e19 2.9 % (legacy scalar 13.2 %), n+ 1e19 51 %, the 4.5e16 transition end 87 %,
the 8.2e15 junction layer 96 %, the 3.4e13 bulk 99.97 %. **Assumed, not from the
sources:** no T dependence of g, N in Eq. 85 is the local density of the same dopant
(no compensation term), no Al impurity band or deep level, and the dopant occupancy is
quasi-static (C is the low-frequency C). Custom parameters:
`DeviceConfig(incomplete_ionization=IonizationModel(acceptor_levels=(DopantLevel(1.0,
0.205, 1.7e-8),)))` (`etna.core.incomplete_ionization`; alpha in eV cm).

Effect on the default device (300 K, legacy -> carrier-dependent):

| quantity                         | legacy                | carrier-dependent     |
| -------------------------------- | --------------------- | --------------------- |
| N_A_ionized (p+)                 | 1.322e18              | 2.913e17              |
| V_bi (psi(cathode) - psi(anode)) | 3.198 V               | 3.142 V               |
| eps*A/C at 0 / -2 / -10 / -30 V  | 1.655 / 6.939 / 9.519 / 9.750 µm | 1.619 / 6.867 / 9.519 / 9.750 µm |
| CCE at -10 V                     | 0.999473              | 0.999466              |
| reach-through                    | 3.10 V                | 3.15 V                |

**Why it is not the default.** The two measured-C-V presets were fitted with the legacy
model. With the carrier-dependent model the parametric default keeps its gates (plateau
rms 0.0137 µm, knee rms 0.2445 µm against the stated < 0.3 µm), but the tabulated preset
fails its knee gate: knee rms 0.0764 µm against <= 0.05 µm (13/13 markers still within
1 sigma; plateau rms 0.0076 µm). Its W(0) moves 1.679 -> 1.654 µm. Making it the default
needs a re-inversion of the C-V with the model on, then re-pinned baselines (and the
N_t refits). Solvers: devsim (1D, 2D), the batched 1D C-V/CCE engine (`run_cv`,
`run_cce`, batches, `run_temperature_sweep`, `run_cce_tilted_alpha`) and the 2D engine
(`run_field` 2D, lateral CCE) carry the model, with parity to devsim at the existing
levels (C <= 3.2e-8, W <= 1e-13 cm, CCE <= 1e-9, 2D psi <= 1e-9 V;
`tests/test_incomplete_ionization_engines.py`). The dark-current, FLASH, fluence and
C-V-at-fluence engines do not carry it: with the model on those paths run on devsim.
The Streamlit sidebar has no toggle yet.

### 5.4 Units and signs

- **Inputs (`DeviceConfig`)**: lengths in µm, doping in cm^-3, T in K, area in cm²,
  N_t in cm^-3 s^-1, E_t in eV below Ec.
- **Bias is negative for reverse bias** in every facade (`v_stop=-50.0` means 50 V
  reverse). Internally the cathode is driven positive (`V_target = -bias`). Forward bias
  is out of scope for the engines.
- **Voltages derived from bias** (`reach_through_V`, `cross_checked_to_V`, V_fd) are
  **positive magnitudes**.
- **Outputs**:

  | output                                  | unit                       |
  | --------------------------------------- | -------------------------- |
  | `run_cv` `y` (C)                        | F, absolute for `area_cm2` |
  | `W_from_C`, `depletion_widths`          | cm                         |
  | `run_field` 1D `x`                      | µm (depth)                 |
  | `run_field` 1D `y`                      | V/cm                       |
  | `mesh.x_coords` / `mesh.y_coords`       | cm                         |
  | `run_cce` `I_collected` / `I_generated` | A/cm²                      |
  | `run_dark_current` `y`                  | A, absolute for `area_cm2` |
  | fluences                                | p/cm²                      |
  | dose rates                              | Gy/s                       |

- **Area.** The default 4 mm² is the area at which N_t was fitted to the 18 pA bound.
  Until 2026-10-18 the default was 1e-4 cm², so older notebooks and scripts show C and I
  400x smaller.

---

## 6. The `run_*` facades

All facades are importable from `etna` (`etna/__init__.py`) and live in
`etna/api/simulation.py`; their docstrings are the reference. Each takes a
`DeviceConfig` first. Timings are single calls on the default device, after the numba
cache is warm.

### `run_cv(config, v_start=0.0, v_stop=-50.0, n_points=40)`

Returns small-signal C-V. `x` is bias (V), `y` is C = dQ/dV (F), the observable of a
measured C-V. Metadata:

- `W_from_C` (eps*A/C, cm) and `depletion_widths` (the n/N_D = 0.5 depletion edge, cm;
  not eps*A/C near full depletion)
- `one_over_C_squared`, `area_cm2`
- `truncated`, `requested_v_stop`
- `reach_through_V`, `reach_through_known`, `reach_through_method`
- provenance (section 7)

1D only. Takes about 0.1 s on the batched engine.

```python
import numpy as np
from etna import DeviceConfig, run_cv

r = run_cv(DeviceConfig(), 0.0, -30.0, 31)
mott = r.metadata["one_over_C_squared"]          # 1/C^2, for a Mott-Schottky plot
print(r.metadata["truncated"], r.metadata["requested_v_stop"], r.x.min())  # False -30.0 -30.0
print(sorted(r.metadata))
```

### `run_cce(config, v_start=-10.0, v_stop=-40.0, n_points=30)`

Returns CCE vs reverse bias for alpha-like generation confined to the n- epi (CCE
relative to the in-epi generation). `x` is bias (V), `y` is CCE in [0, 1]. Metadata:

- `I_collected` (per bias point, A/cm²) and `I_generated` (A/cm²)
- `truncated`, `requested_v_stop`
- `reach_through_V`, `reach_through_known`
- provenance

1D only (a 2D config raises `NotImplementedError`). Takes about 0.1 s on the engine.

```python
from etna import DeviceConfig, run_cce

r = run_cce(DeviceConfig(), v_start=-1.0, v_stop=-10.0, n_points=10)
print(r.metadata["solver"], [round(float(v), 4) for v in r.y[[0, 1, -1]]])  # batched [0.9743 0.9944 0.9995]
```

### `run_cce_tilted_alpha(config, v_start=0, v_stop=-60, n_points=61, theta_deg=65, alpha_range_um=18, dead_layer_um=0, shape="bragg", generation_in_p_plus=False)`

Models the microdosimeter's tilted-alpha test. It runs on the batched engine only and
raises `etna.api.simulation.OutsideModelValidity` outside the engine box. Metadata:

- `cce_vs_E0` (collected over E0/E_pair, which includes dead-layer and p+ losses)
- `projected_end_um`, `dead_layer_energy_fraction`, `generation_shape`
- `reach_through_V`

The Bragg shape is qualitative, not ASTAR/SRIM data. The model **does not** reproduce the
paper's CCE-vs-bias shape (section 8).

```python
from etna import microdosimeter_config, run_cce_tilted_alpha

r = run_cce_tilted_alpha(microdosimeter_config())
print(r.metadata["projected_end_um"], r.metadata["reach_through_V"])  # ~7.61 um, 9.0 V
```

### `run_field(config, bias_V=-50.0)`

Returns the field and potential at one reverse bias.

- **1D**: `x` is depth (µm), `y` is |E| (V/cm). Metadata holds `potential` (V) and
  `net_doping` (cm^-3). 1D always runs on **devsim** (about 0.3 s).
- **2D** (`half_width_um` set): `x` and `y` are **empty**. Use `mesh.x_coords`
  (lateral, cm), `mesh.y_coords` (depth, cm) and `mesh.node_values["ElectricField"]`.
  The 2D engine serves configs inside its box (about 0.5 s); otherwise devsim runs.
- A legacy-model 2D bias past the contact wall raises `OutsideModelValidity`
  (section 8).

```python
from etna import DeviceConfig, run_field
from etna.api.results import peak_field_1d

cfg = DeviceConfig()
r = run_field(cfg, bias_V=-30.0)
peak, depth_um, excluded = peak_field_1d(r.x, r.y, cfg)   # skips the nm spikes at abrupt steps
print(r.metadata["solver"], f"{peak:.3g} V/cm at {depth_um:.2f} um", excluded)  # devsim 9.06e+04 V/cm at 0.40 um True

r2 = run_field(DeviceConfig(half_width_um=50.0), bias_V=-30.0)
print(r2.metadata["solver"], r2.x.size, r2.mesh.x_coords.size, r2.mesh.y_coords.max() * 1e4)  # engine2d 0 17784 ~11.3 um
```

### `run_radiation_damage(config, fluences=None, V_bias=-40.0, proton_energy_MeV=5.6)`

Returns CCE vs proton fluence. `x` is fluence (p/cm²), `y` is CCE. Metadata:

- `valid` (bool per point) and `reason` (str per point)
- `compensation_fluence` (p/cm²)
- `solver` (the solver that ran the points: `batched` on the default stack, `devsim`
  outside the engine box, joined when mixed) and `solver_per_point`
- `V_bias`, `energy_MeV`

The donor-removal model fully compensates the low-doped epi at a fluence that depends on
the profile: 1.01e12 p/cm² at 5.6 MeV (3.48e12 at 62 MeV) on the default stack, 3.88e11
(1.33e12) on the tabulated preset. Points at or beyond it are **NaN** with a reason, not solved. Choose the grid
from `compensation_fluence`, and always check `valid`. With `fluences=None` the grid
adapts to the device (`etna.api.simulation.default_fluence_grid`): 8 log-spaced points
from 0.01x to 3x the compensation fluence: 6 are solved, and the last 2 are normally NaN
with a reason (the 3x point always; `compensation_fluence` is a lower bound on the mesh
threshold), showing where the model ends. Explicit grids are used as given.

NIEL hardness factors (audit C-5 / B-07, since 2026-10-20): kappa(E) =
NIEL_p,Si(E) / NIEL_n,Si(1 MeV) from SR-NIEL-7 proton NIEL (1-1000 MeV, 59 energies,
log-log interpolation) over the ASTM 95 MeV mb normalisation: 6.67 at 5.6 MeV, 2.75 at
30, 1.94 at 62, 1.82 at 70, 1.12 at 150 MeV (the old placeholders were 0.50 / 0.35 /
0.33 / 0.22 at 30-150 MeV and 5.6 MeV was silently clamped to 0.5). Energies outside
1-1000 MeV raise `ValueError`. This is the 1-MeV-neq convention of the defect
introduction rates (Burin et al., fitted to reactor data labelled in Si neq). Open: the
SiC/Si NIEL ratio of the reactor neutrons is not available (SR-NIEL's SiC neutron
calculator needs a login), so the SiC-specific normalisation is unvalidated; if that
ratio exceeds the proton one (0.94 at 5.6 MeV, 0.85 at 62 MeV:
`radiation_damage.niel_ratio_sic_si`), kappa overstates the SiC damage per proton. At
5.6 MeV NIEL scaling is also known to overpredict measured Si damage by tens of percent.
Data and provenance: `data/niel/` (raw SR-NIEL pages, CSVs, fetch script, README). The
Radiation Damage page offers 30, 62, 70 and 150 MeV only.

```python
import numpy as np
from etna import DeviceConfig, run_radiation_damage
from etna.api.simulation import compensation_fluence

cfg = DeviceConfig()
phi_c = compensation_fluence(cfg, 5.6)                       # 1.01e12 p/cm2
r = run_radiation_damage(cfg, fluences=np.geomspace(1e9, 0.9 * phi_c, 5), V_bias=-40.0)
print(f"{phi_c:.3g}", r.metadata["valid"].all(), np.round(r.y, 4))
```

### `run_dark_current(config, v_start=0.0, v_stop=-100.0, n_points=20, N_t=None, S_n=None, S_p=None, solver="auto")`

Returns the total reverse leakage. `x` is bias (V), `y` is the total current (A, for
`area_cm2`). Metadata:

- `I_SRH`, `I_TAT`, `I_SRV` (a model partition, not a measured decomposition; I_SRV is
  zero at the ohmic cathode). They carry the sign of I_total (generation positive):
  I_SRH + I_TAT + I_SRV ≈ I_total, with I_SRH a tiny negative net-recombination term
- `N_t`, `E_t`, `area_cm2`
- `truncated`, `reached_v_stop`
- `reach_through_V`, `past_reach_through`, `validity_note` (filled for the legacy model
  only)

Runs on the batched 1D engine inside its box (section 7; the same TAT/SRH/SRV model
and frozen-field fixed point as devsim, parity <= ~1e-11 relative, about 20 ms for 20
points), otherwise on devsim (about 2-3 s). `solver="devsim"` forces the devsim reference,
`solver="engine"` raises outside the box. Quote the result as a budget estimate
(section 2).

```python
from etna import DeviceConfig, run_dark_current

r = run_dark_current(DeviceConfig(), 0.0, -60.0, 7)
print([f"{i * 1e12:.2f}" for i in r.y[[1, 3, 6]]], "pA at -10/-30/-60 V")  # 3.99 11.88 18.00
```

Model (since 2026-10-20, audit B-08): Z1/2 SRH plus the fitted effective generation
N_t f_depl, both with the Hurkx (1992) trap-assisted-tunnelling enhancement tau ->
tau / (1 + Gamma), Gamma = 2 sqrt(3 pi) (F/F_Gamma) exp[(F/F_Gamma)^2], F_Gamma =
sqrt(24 m_t (kT)^3) / (q hbar) = 3.69e5 V/cm at 300 K (m_t = 0.25 m0, assumed). Gamma is
~1.8 at 1e5 V/cm, so the bias shape is f_depl x (1 + Gamma(E)); before, a Kt branch kept
Gamma = 1 below ~4.5 MV/cm. N_t was refitted so that I(300 K, -60 V, 4 mm²) is still
18.000 pA: 2.849024e12 -> 1.355575e12 (default), 2.849219e12 -> 1.356120e12 (tabulated),
2.805336e12 -> 1.343970e12 (previous default); the legacy preset keeps its historical
2.8e12 (37.4 pA at -60 V now). Other biases move: -10 / -30 / -100 V give 3.99 / 11.88 /
24.87 pA (were 6.30 / 15.76 / ~18.2). The closed form is the saddle-point value of
Hurkx's integral and is ~30 % above the integral at 1e5 V/cm; fields above its validity
limit (~1.1e6 V/cm at 300 K) are logged (devsim) or sent to devsim (engine).

### `run_temperature_sweep(config, temperatures=None, voltage=-30.0, method="hecht")`

Returns CCE vs T at one bias (default 250-350 K, 5 points); `config.T` is ignored.
`"hecht"` is an analytical Hecht CCE on the engine's depletion width, `"dd"` a
drift-diffusion CCE. Temperatures inside the engine T box (150-450 K) run on the engine
(milliseconds); others run on devsim. `metadata["solver_per_T"]` records the path per
point.

```python
import numpy as np
from etna import DeviceConfig, run_temperature_sweep

r = run_temperature_sweep(DeviceConfig(), np.array([250.0, 300.0, 350.0]), voltage=-30.0, method="dd")
print(r.metadata["solver"], np.round(r.y, 4))
```

### `run_flash_recombination(config, dose_rates_Gy_s=None, V_bias=-30.0, E_MeV=62, solver="auto")` and `run_transient(config, dose_rates=None, V_bias=-30.0, epi_thickness_cm=None)`

FLASH dose-rate sweeps. `x` is dose rate (Gy/s), `y` is CCE. `run_transient` also
returns `transient_cce_raw` (unclipped). `run_flash_recombination` runs on the batched
engine inside its scope (section 7; about 40 ms, parity <= ~7.5e-12 with Auger on),
otherwise on devsim (about 1-3 s); `run_transient` runs on devsim (about 2 s). Any
variation with dose rate is a **numerical sensitivity bound**. No plasma or high-injection physics is modelled. Both
solvers of the steady-state sweep solve **SRH + Auger** (C_n 5e-31, C_p 2e-31 cm^6/s,
Galeckas et al. 1997, 300 K constants; override with
`etna.core.flash_recombination.cce_vs_dose_rate(..., C_n=, C_p=)`). Auger lowers CCE by
a dose-rate-independent offset: 3.4e-7 on the default stack at -30 V, up to 5.6e-5 on
notebook 05b's grid (5 um epi, N_D_bulk 5e14, -10 V), from minority holes recombining in
the 1e19 n+ substrate; it scales ~linearly with C_n, so the literature spread (up to
~2e-30) allows up to ~2e-4 there. It adds no dose-rate dependence (<= 1.5e-8). Until
2026-10-20 the Auger term was silently dropped (both solvers SRH-only). `run_transient`
still is SRH-only (open item, `.issues/TECH-DEBT-2026-10-18.md`).

```python
from etna import DeviceConfig, run_flash_recombination, run_transient

f = run_flash_recombination(DeviceConfig())
t = run_transient(DeviceConfig())
print(f.y.round(6), t.metadata["transient_cce_raw"].round(6))
```

#### Transient terminal current (conduction + displacement)

The transient drivers (`TransientSolver.simulate_pulse`, `simulate_pulse_train`,
`single_particle.simulate_single_particle`) record three series per BDF1 step:

| key                     | quantity                                                                                       |
| ----------------------- | ---------------------------------------------------------------------------------------------- |
| `currents`              | conduction current, `extract_contact_current` (continuity equations only)                      |
| `currents_displacement` | displacement current, `displacement_current(Q_k, Q_k-1, dt_k)` of the Poisson contact charge `extract_contact_charge` |
| `currents_total`        | terminal current = conduction + displacement: the physical I(t)                                 |

- **CCE and Q are unchanged.** They still integrate `currents`. The displacement term
  telescopes to `Q_contact_end - Q_contact_start`. So once the device is back in its dark
  state, the total and the conduction current carry the same charge (`Q_induced` vs
  `Q_collected`).
- **The pulse shape comes from the total current.** `analyze_current_pulse(..., currents_total=...)`
  takes `I_peak`, `t_peak` and `t_collection` from the total current. The conduction-only
  values are kept as `*_conduction`.
- **`t_collection` definition.** It is the time from the injection start at which the
  cumulative charge first reaches 95 %, interpolated linearly inside the step.
- **Validation.** The tests are in `tests/test_displacement_current.py`.
  - Parallel plate: a 1e9 cm^-2 pair sheet in a fully depleted 10 µm n- gap (no n+,
    -30 V) gives I = q N v / d within 3 %, both while both carriers drift and once only
    the electrons do. The induced charge is q N to 2e-3.
  - 1D Shockley-Ramo identity: at every step, total(anode) = -total(cathode) =
    (1/L) ∫ (J_n + J_p) dx (`ramo_current_1d`).
- **How large the correction is depends on the contact.**
  - _Contact on depleted material_ (legacy no-n+ cathode at full depletion, the
    parallel-plate test): the conduction current is ~0 until the carriers arrive, then
    it spikes. Its peak is 18.36 A/cm^2 at 30.3 ps, against the true 5.69 A/cm^2 at
    2.5 ps. In 2D (legacy stack, hw 5 µm, 50 V) the conduction currents of the two
    contacts miss Kirchhoff by up to 1.5x the anode current; the totals sum to zero
    to 3e-8.
  - _Contact on a heavily doped ohmic layer_ (every current preset): conduction already
    carries nearly all of the terminal current. The displacement share is only
    significant on steps not long compared with the layer's dielectric relaxation time
    (~0.1 ps for the p+). Default 2D device (hw 50 µm, 50 V, LET 10, anode), measured
    2026-10-20:

    | quantity          | displacement share or effect                                  |
    | ----------------- | ------------------------------------------------------------- |
    | 1 ps injection step | 5 % of the current                                          |
    | later steps       | <= 1.6 % of the peak                                          |
    | I_peak            | 5.49e-5 total vs 5.42e-5 conduction A/cm                      |
    | t_collection      | 0.11679 ns (new definition) vs 0.11739 ns (old step-end value) |
    | CCE               | 0.999970 either way                                           |

- **Old definition and the stale table column.** Before 2026-10-20, `t_collection` was the
  end time of the BDF1 step that crossed 95 % of the conduction charge (`searchsorted` on
  a cumulative sum that assumes it is monotonic). That made it quantised to the 10 %
  step growth. The tracked `data/cce_let_table_100um.json` / `_300um.json` predate this
  change.
  - Their `t_collection_ns` column (0.208 ns for every LET) is stale: the 2026-10-18
    rebuild gave 0.117 ns.
  - No code reads that column: `load_cce_let_table` interpolates CCE only, and CCE is
    unaffected.
  - The tables were not regenerated (~40 min, `ETNA_REBUILD_CCE_TABLES=1` in notebook 16).

### `run_microdosimetry(config, mc_csv_path, sv_thickness_um=10.0, sv_width_um=100.0, sv_depth_um=None, chord_model="cauchy")`

Pure data pipeline (no devsim). It reads Monte Carlo per-step depositions from a CSV
(columns `event_id, x, y, z, edep`; edep in keV, positions in cm) and returns the lineal
energy spectrum: `x` is y (keV/µm), `y` is y·d(y). Metadata holds `y_F`, `y_D`, `f_y`,
`l_bar_um` (mean chord length), `chord_model` and `sv_dimensions_um` (width, depth,
thickness). `config` is carried for provenance only. `mc_csv_path` may be
a path or a file-like object (e.g. `io.BytesIO`). It raises `ValueError` for non-positive or
non-finite sensitive-volume dimensions, an unknown `chord_model`, a CSV with no events,
non-finite or negative `edep`, or no event depositing energy.

Mean chord length (audit D-15 follow-up, 2026-10-20): `chord_model="cauchy"` (default) is
Cauchy's l̄ = 4V/S of the finite box w × d × t (`sv_depth_um=None` means a square footprint,
d = w), so the width matters: 16.67 µm for DESIGN-2's 100 × 100 × 10 µm SV (the new default
width; it was 150 µm), 17.65 µm for 150 × 150 × 10 µm. `chord_model="slab"` gives the
infinite-slab limit l̄ = 2t (20 µm, the value used before 2026-10-20). Both assume isotropic
irradiation; for a normally incident beam the chord is t, which neither gives. The Cauchy
value raises y (and y_F, y_D) by 2t/l̄ = 1.2 for the default SV relative to the slab value.

```python
from etna import DeviceConfig, run_microdosimetry

r = run_microdosimetry(DeviceConfig(), "data/synthetic_mc_events.csv")
print(round(float(r.metadata["y_F"]), 2), round(float(r.metadata["y_D"]), 2), round(r.metadata["l_bar_um"], 3))  # 20.67 63.8 16.667
s = run_microdosimetry(DeviceConfig(), "data/synthetic_mc_events.csv", chord_model="slab")
print(round(float(s.metadata["y_F"]), 2), round(float(s.metadata["y_D"]), 2), s.metadata["l_bar_um"])  # 17.23 53.22 20.0
```

#### Lateral (edge) CCE weighting of MC hits (core API, audit B-11)

`etna.core.mc_coupling.process_mc_ensemble(events, cce_interp, cce_lateral=scan)`
multiplies each event's CCE(LET) by f(r) = CCE_2D(r) / CCE_2D(0) of a 2D lateral scan
(`etna.core.charge_collection_2d.cce_lateral_scan`), r the distance of the event's
edep-weighted lateral centroid (x, y; z is depth) from `sv_center_cm`
(`lateral_metric="chebyshev"` for a square SV, `"radial"`, `"x"`). Events beyond the scan
are clamped to the edge value (`lateral_outside="edge"`, counted in `n_outside_lateral`
and logged), set to 0 (`"zero"`) or rejected (`"raise"`). Without `cce_lateral` the CCE is
position-independent, as before. Assumptions: CCE(LET, r) = CCE(LET) f(r) (separable, f
from the scan's alpha-stripe generation), the centroid stands for the whole track. The
default 2D device is laterally uniform, so its scan gives f = 1 (no edge physics); a
structure with an edge does: the guard-ring device (`alternative_structures.
create_guard_ring_device`, SV half-width 50 µm, 3 µm gap, 10 V, anode read-out) gives
f = 0.9959 at 41 µm, 0.977 at 46 µm and 0.49 at 51.5 µm (middle of the gap;
`lab/edge_cce/notes/LOG.md`).
`run_microdosimetry` does not apply it (it has no CCE step).

### Batches and sweeps

- `run_cv_batch(configs, ...)` and `run_cce_batch(configs, ...)` solve many configs
  together on the engine, grouped by mesh. The rest go through devsim one by one.
  A failing config yields an empty `SimResult` with `metadata["error"]` and
  `solver == "failed"`; the other configs are kept.
- `ParametricSweep(base_config, param, values, sim_fn, sim_kwargs)` clones the config
  once per value (`dataclasses.replace`) and calls the facade.

```python
from etna import DeviceConfig, ParametricSweep, run_cce, run_cce_batch, run_cv_batch

rs = run_cv_batch([DeviceConfig(epi_thickness_um=e) for e in (8.0, 10.0, 12.0)], 0.0, -20.0, 21)
print([r.metadata["solver"] for r in rs])                      # ['batched', 'batched', 'batched']
cs = run_cce_batch([DeviceConfig(N_D_bulk=n) for n in (3e13, 1e14)], -5.0, -30.0, 6)
print([round(float(r.y[-1]), 5) for r in cs])                  # CCE at -30 V per config

sweep = ParametricSweep(DeviceConfig(), "epi_thickness_um", [8.0, 10.0, 12.0], run_cce,
                        {"v_start": -5.0, "v_stop": -30.0, "n_points": 6})
for r in sweep.run():
    print(r.config.epi_thickness_um, r.metadata["solver"], round(float(r.y[-1]), 5))
```

### `run_all(config, sections=None)`

Runs, in one process: `cv`, `cce`, `field`, `damage`, `dark`, `temperature`, `flash`,
and, on request, `sweep` (a small 2D-engine microdosimetric grid; section 11). It asserts
between sections that no devsim device leaked. On the default device the whole run takes
about 0.5 s with every fast path active (all sections but the 1D `field` on the engines;
legacy preset 2.1 s: its 0..-100 V dark sweep runs on devsim). A section's `status` is `ok`, `warning` (it ran
but returned no valid point; `reason` says why), `error` or `skipped`; `rep.ok` is
True when no section is `error` or `warning` (skipped sections do not count). A curve with invalid points starts its key
figures with "k of n points valid", and the damage row adds the compensation fluence.

```python
from etna import DeviceConfig, run_all

rep = run_all(DeviceConfig(), sections=("cv", "cce", "field"))
print(rep.ok, [(s.name, s.status, s.solver) for s in rep.sections])
print(rep.get("field").result.metadata["solver"])
```

---

## 7. Solver provenance and fast paths

Every `SimResult.metadata` carries:

| key                             | meaning                                                                                                                                                                                                                               |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `solver`                        | what produced the numbers: `"batched"` (1D numba engine), `"engine2d"` (2D numba engine), `"devsim"`, `"failed"` (batch slot that raised), or a descriptive string (temperature sweep: `"hecht (analytical; W from batched engine)"`) |
| `validated_domain`              | True if the device lies in the box where the engine was cross-checked against devsim                                                                                                                                                  |
| `domain_note`                   | why not (None inside)                                                                                                                                                                                                                 |
| `cross_checked_to_V`            | engine results only: the deepest reverse bias pinned against devsim by the parity tests (100 V with an n+ layer, 40 V for the legacy model; dark current: 100 V, legacy 60 V; FLASH: 50 V)                                            |
| `fallback_reason`               | present when an engine was tried and devsim took over (e.g. `"batched engine did not converge"`, or `run_field`'s `"2D engine not applicable: ..."`)                                                                                  |
| `truncated`, `requested_v_stop` | a devsim sweep stopped before `v_stop` (convergence wall); `x` holds only the reached points                                                                                                                                          |

Which solver runs:

- **`run_cv`, `run_cce`, `run_cv_batch`, `run_cce_batch`, `run_temperature_sweep`,
  `run_cce_tilted_alpha`**: the batched 1D engine (`etna/core/batched/`) inside
  `_ENGINE_BOX` (`etna/api/simulation.py`):

  | parameter    | box                     |
  | ------------ | ----------------------- |
  | T            | 150-450 K               |
  | N_A          | 1e18-1e19               |
  | N_D_junction | 1e15-1e16               |
  | N_D_bulk     | 2e13-2e15               |
  | L_transition | 0.3-3 µm                |
  | uniform N_D  | 1e14-1e16               |
  | epi          | 5-30 µm                 |
  | p+           | 0.3-3 µm                |
  | n+           | 1-5 µm at 1e18-1e19     |
  | transition   | 0.5-2.5 µm to 1e16-1e19 |
  | table N_D    | 5e12-5e16               |

  Outside the box, or if the engine does not converge, devsim runs. `run_cce_tilted_alpha`
  has no devsim path and raises instead. Ask before running:

  ```python
  from etna import DeviceConfig
  from etna.api.simulation import validated_domain_reason

  print(validated_domain_reason(DeviceConfig()))        # None: inside
  print(validated_domain_reason(DeviceConfig(T=500.0))) # 'T=500 outside [150, 450]' -> devsim
  ```

- **2D `run_field`** and the microdosimetric sweep use the 2D engine
  (`etna/core/engine2d/`) inside its boxes (`ENGINE2D_BOX` / `ENGINE2D_DARK_BOX`, epi
  5-30 µm, n+ and transition boxes). Otherwise devsim runs.
- **Dark current and FLASH** (since 2026-10-18): `run_dark_current`,
  `run_flash_recombination` and the core `dark_current_vs_fluence`,
  `sweep_iv_vs_temperature`, `sensitivity_sweep`, `cce_vs_dose_rate` and
  `parametric_cce_sweep` run on the batched engine (`etna/core/batched/dark.py`) when
  the device is in the box above (graded/uniform, 150-450 K, reverse bias). Limits on
  top of it: the legacy model only to -60 V for the dark current (deeper sweeps go to
  devsim, whose legacy truncation record is the documented one) and -40 V for FLASH;
  FLASH widens N_D_junction to 1e15-2e17 (notebook 05b's scaled junction doping). Each
  takes `solver="auto" | "engine" | "devsim"`. Parity: `tests/test_dark_current_engine.py`
  (I_total, I_TAT <= 1e-9 relative, measured <= 1.4e-11), `tests/test_flash_engine.py`
  (CCE <= 1e-9, measured <= 7.5e-12 with Auger on).
- **2D dark current** (`etna.core.optimization.get_dark_current_2d`): a device the 2D
  engine accepts is laterally uniform, so its current is the 1D engine's J times the SV
  area (vs the devsim 2D solve: 2.3e-5 / 5.0e-6 since the Hurkx Gamma, 2.4e-7 / 3.3e-7
  before; the node field discretisation entering f_depl and Gamma(E);
  `tests/test_dark_current_2d_engine.py`); otherwise devsim 2D.
- **Incomplete ionization** (`incomplete_ionization=True`, §5.3): the 1D C-V/CCE engine
  and the 2D engine solve it; dark current, FLASH, fluence and C-V-at-fluence go to
  devsim (`metadata["domain_note"]` / the engine reason names `incomplete_ionization`).
- **devsim only**: 1D `run_field`, `run_transient`, `dark_current_sweep(device_info)` (a
  devsim device you built), and any engine fallback.
- **Fluence sweeps**: `run_radiation_damage` (through the core `cce_vs_fluence`) uses the
  batched engine where its box allows (`tests/test_fluence_engine.py`) and devsim
  otherwise. `metadata["solver"]` names the solver that ran the points and
  `metadata["solver_per_point"]` records it per point (empty for the NaN points beyond
  compensation).

The **box** is a statement about where the two solvers agree. It is not a statement
about physical validity: a device inside the box is solved the same way by both codes,
and that says nothing about agreement with a measurement.

Logging and environment:

- devsim's per-solve stdout log (about 65 lines per solve) is silenced by
  `etna.core.devsim_lock.quiet_devsim()`, a context manager the facades already use. Set
  `ETNA_DEVSIM_VERBOSE=1` to see it again. Etna's own diagnostics go through `logging`
  (logger names `etna.*`).
- `ETNA_FIGURE_DIR` sets where `etna.core.plotting.save_figure` writes (default:
  `./figures` of the current working directory). Run notebooks from `notebooks/` (giving
  `notebooks/figures/`, which is gitignored) or set the variable. Running from the repo
  root writes over the tracked `figures/*.png`.
- numba parallel kernels use all cores. Worker pools (`n_jobs > 1` in the optimisation
  sweeps) pin `NUMBA_NUM_THREADS` / `OMP_NUM_THREADS` to 1 per worker.

---

## 8. Validity and limits

**Full depletion has two definitions.** Both are stated in `etna.export.labels.VFD_DEFINITION`.

- _Analytic V_fd_: the Poisson double integral of the graded profile up to the **start
  of the n-/n+ transition**, minus V_bi. Default stack: **2.74 V**.
- _Reach-through_ (`metadata["reach_through_V"]`): the smallest |V| at which the
  depletion edge is within one Debye length of the transition start (engine dark sweep,
  0.05 V resolution). Default stack: **3.10 V**.
- The 2D sweep gate is |V| >= **1.2 V_fd = 3.28 V**.
- Depleting the transition itself is the slow -10..-60 V creep of eps*A/C, not full
  depletion.

| preset                                       | V_fd                 | reach-through |
| -------------------------------------------- | -------------------- | ------------- |
| previous default (`nplus_calibrated_config`) | 9.93 V               | 12.3 V        |
| legacy                                       |                      | 17.0 V        |
| tabulated                                    | 3.37 V (decade rule) | 4.85 V        |
| microdosimeter                               |                      | 9.0 V         |

DESIGN-1 still recommends >= -10 V operation, as margin for the doping uncertainty below.

**The doping profile is an inversion of one C-V curve, not SIMS.** The default profile
reproduces the digitised Photons Fig. 4b (230 points, 0..-60 V) by construction, so C-V
agreement is not validation. Its bulk value, 3.4e13 cm^-3, lies **below** the documented
0.5-1e14 cm^-3. The epi/n+ position is set by the fit. If the bulk is at the documented
value, V_fd would be 3.9-7.6 V. The -30 V residual of the old model had four one-knob
explanations (fringe capacitance plus oxide charge, buffer, thinner epi, larger area);
they are discussed in `.issues/perf/PERF-05-default-stack-followups.md` and
`lab/finite_implant/notes/LOG.md`.

**What a measurement would settle:**

| measurement                                | what it settles                                                     |
| ------------------------------------------ | ------------------------------------------------------------------- |
| SIMS or SRP through the epi/n+ interface   | the profile directly                                                |
| C-V to -100..-200 V                        | the transition form                                                 |
| C-V on two die sizes                       | fringe capacitance vs area                                          |
| leakage vs bias and T                      | the shape of the dark-current model; today only one point is fitted |
| a CCE-vs-bias measurement on the dosimeter | a check of the transport model                                      |

**Peak field.** Abrupt doping steps (the p+/n- junction, the epi/n+ step) give nm-scale
field spikes that depend on the mesh. `peak_field_1d` (`etna/api/results.py`) excludes
±0.1 µm around each step (`INTERFACE_SPIKE_WINDOW_UM`). `run_all`, the exports and the
viz all use it, so they report one number: 9.06e4 V/cm at -30 V and 1.21e5 V/cm at
-60 V, at x = 0.40 µm, on the default stack.

**Convergence walls and the outside-validity wall.** With the n+ layer (all current
presets except legacy) devsim has no contact wall: on the default stack it was checked
to converge to -150 V, and the engines are cross-checked against it to -100 V. The **legacy** model puts the ohmic cathode on the n- epi, so
devsim stops converging around -64 V at 300 K (where it stops depends on the
linear-algebra build: a 2 V-step dark sweep reaches -100 V on macOS but -78.9 V on the ubuntu
CI runner). This is a contact artefact, not punch-through.

- `run_cv` / `run_cce` / `run_dark_current` on the devsim path return the reached points
  with `truncated=True`.
- The batched engine does not reproduce the wall, and returns the full grid with
  `truncated=False`.
- 2D `run_field` on the legacy model raises `OutsideModelValidity` _before_ solving when
  the cathode-cell field would exceed 5e4 V/cm, with a message meant for the user.

**Compensation fluence.** See `run_radiation_damage`. Beyond `compensation_fluence(config,
E)` the epi is fully compensated in this model, and the point is NaN with a reason.

**Microdosimeter alpha test not reproduced.** With the preset (N_D 1e14), the model CCE is
flat to 0.5 % from 0 V. The data drop to 0.908 at 0 V and rise 1 % from -20 to -50 V. No
drift-diffusion variant fits (chi2 >= 228). A field-dependent columnar-recombination loss
does fit, but it is lab-only (`lab/microdosimeter/notes/LOG.md` §9).

**Incomplete ionization** is the legacy scalar by default; the carrier-dependent model
is opt-in (§5.3) and not yet calibrated against the C-V.

**Not modelled:**

- high injection or plasma (FLASH)
- edge termination and guard rings on the default path
- oxide charge
- impact ionisation and breakdown
- trapping within a shaping time and ballistic deficit
- 3D effects

---

## 9. The Streamlit app

```bash
uv run streamlit run app/main.py --server.address localhost
```

The committed `.streamlit/config.toml` binds the app to localhost when it is launched from
the repo root (Streamlit reads the config of the launch directory); the flag makes that
explicit from any directory. The app has no authentication: for shared use pass
`--server.address 0.0.0.0` behind an authenticating reverse proxy (audit C-16). One solve at a time:
devsim state is process-global, and the facades serialise access through a lock.

**Sidebar** (`app/components/device_sidebar.py`). It is rendered on every page; the
configuration persists across pages.

- _Dimensionality_ (1D/2D). 2D shows _Half-width (µm)_.
- _Doping profile_ (graded/uniform).
- _Epi thickness_, _Substrate thickness_ (the p+ anode), _N_A_, the graded triplet
  (_N_D junction_, _N_D bulk_, _Transition length_) or _N_D uniform_.
- _graded n-/n+ transition (C-V fit)_ checkbox, with its thickness and end doping.
- _Temperature_, _Area (cm²)_ (default 4.000e-02).
- Three preset checkboxes. They override the fields; if more than one is checked, the
  first wins:
  - _historical model (cathode directly on n- epi)_: `legacy_config`
  - _previous calibrated model (n+, no buffer)_: `nplus_calibrated_config`
  - _tabulated measured-C-V profile (19 knots, opt-in)_: `measured_cv_tabulated_config`
- Invalid values give a red error that names the field, and the pages refuse to run.
- A caption says whether a 1D device is inside the engine/devsim cross-check box.

**Pages** (`app/main.py`, `app/workflows/*.py`). Each page has a _Run simulation_
button, caches its result against the device and warns when the sidebar has changed
since. Each shows the solver and validity notes, and offers a CSV download whose `#`
header carries solver, domain, config and version.

| page                    | purpose                                                                         | inputs                                                                                                                                                             | outputs                                                                                                                                                                                                                                                                                                                                                                                                          |
| ----------------------- | ------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Home                    | orientation, current config summary                                             | none                                                                                                                                                               |                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Run All                 | every section for the sidebar device (`etna.run_all`)                           | section multiselect (the 2D-engine sweep is opt-in)                                                                                                                | summary table (status, solver, key figures, reason, wall time), a curve per section, JSON/CSV                                                                                                                                                                                                                                                                                                                    |
| C-V Analysis            | `run_cv`                                                                        | max reverse bias (default -50 V), points (40)                                                                                                                      | C-V, Mott-Schottky, eps*A/C and depletion edge vs V, full-depletion marker, CSV                                                                                                                                                                                                                                                                                                                                  |
| Charge Collection (CCE) | `run_cce`                                                                       | V_start (-10), V_stop (-40), points (30); reverse bias only                                                                                                        | CCE vs bias, CSV                                                                                                                                                                                                                                                                                                                                                                                                 |
| Field Map               | `run_field`, 1D or 2D                                                           | reverse bias (default -50 V, down to -200)                                                                                                                         | 1D: E(x), potential, net doping, CSV. 2D: heatmap viewer (field, potential, doping)                                                                                                                                                                                                                                                                                                                              |
| Radiation Damage        | `run_radiation_damage`                                                          | proton energy (30/62/70/150 MeV), fluence range (default 0.7e-3 to 0.7 x the compensation fluence), points (6), bias (-20 V)                                       | CCE vs fluence with NaN gaps explained. A persistent banner says kappa/NIEL are placeholders                                                                                                                                                                                                                                                                                                                     |
| Dark Current            | total leakage vs T at one bias (`ParametricSweep` over T of `run_dark_current`) | T range (250-400 K), points (6), bias (-20 V), advanced N_t                                                                                                        | I_total vs T, CSV. A calibration banner is shown; only the total is shown                                                                                                                                                                                                                                                                                                                                        |
| Microdosimetry          | `run_microdosimetry`                                                            | MC events CSV upload, SV thickness (10 µm) and width                                                                                                               | y·d(y), y_F / y_D, CSV. Independent of the sidebar                                                                                                                                                                                                                                                                                                                                                               |
| Batch Sweep             | `ParametricSweep`                                                               | one of epi thickness, p+ thickness, N_A, N_D junction, N_D bulk, transition length, T, area; comma-separated values; facade (CCE vs bias, C-V, CCE vs temperature) | overlaid curves, bulk CSV, per-value failures listed                                                                                                                                                                                                                                                                                                                                                             |
| Export                  | engineering files (section 10)                                                  | spec sheet (reuses the Run All report when it matches); profile bias (default -30 V); animation button                                                             | layout parameter file in use (upload an alternative TOML) with its PLACEHOLDER list; spec sheet JSON/CSV/Markdown/PDF; cross-section SVG/PNG; mesh/profile CSV, VTK and VTU (ParaView/VisIt/meshio); two 480p manim MP4s (needs `viz`); reference layout SVG and GDSII (needs `gdstk`) with the disclaimer "reference layout, not TCAD-validated; not a fabrication layout while any dimension is a PLACEHOLDER" |

Page details not in the table:

- The Field Map viewer draws the anode/cathode contacts (`MeshData.contacts`) as lines at
  their depth and, in 2D, outlines the region (`MeshData.regions`).
- The Dark Current CSV has columns `T_K, I_total_A`.
- Batch Sweep values must be finite and positive; sweeping the area with a CCE facade is
  disabled (CCE does not depend on area).
- The shared CSV serializer (`app/components/results.py::to_csv_bytes`) covers every facade:
  `run_temperature_sweep`, `run_flash_recombination` and `run_transient` included, and a native
  `run_dark_current` result exports `bias_V, I_total_A` with the truncation record. The C-V
  CSV column for 1/C² is `one_over_C2_per_F2`.

As of the v5.0 milestone audit (`.planning/milestones/v5.0-MILESTONE-AUDIT.md`), 7 of the
21 notebook workflows had no page. The Run All page has since added the FLASH
steady-state section (notebook 04). Notebooks 08, 12, 14, 16, 19 and 20 remain
notebook-only.

---

## 10. Exports (Python)

`etna.export` holds pure builders: they take a config plus finished results and never
simulate. Each one carries `etna.export.VALIDATION_LABEL`: "simulation output, validated
against devsim reference; not validated against measurement except where stated".

```python
from etna import DeviceConfig, run_all, run_field
from etna import export as ex

cfg = DeviceConfig()
rep = run_all(cfg, ex.SPEC_SECTIONS)                 # ('cv', 'cce', 'field', 'dark')
spec = ex.build_spec(cfg, rep)
md, js, csv = ex.spec_to_markdown(spec), ex.spec_to_json(spec), ex.spec_to_csv(spec)
pdf = ex.spec_to_pdf(spec, cfg)                      # bytes
svg = ex.cross_section_svg(cfg)                      # depth to scale, lateral schematic
prof = run_field(cfg, bias_V=-30.0)
vtk = ex.profile_to_vtk(prof)                        # legacy ASCII VTK; profile_to_csv too
vtu = ex.profile_to_vtu(prof)                        # XML UnstructuredGrid (lines / quads)
layout = ex.layout_svg("dosimeter")                  # REFERENCE-NOT-FAB; layout_gds needs gdstk
print(len(md) > 0, pdf[:4], ex.gds_available()[0])
```

### Mesh / profile files

`profile_to_csv` / `profile_to_vtk` / `profile_to_vtu` write, per mesh node of a
`run_field` result: `NetDoping_cm-3` (N_D - N_A, ionized acceptors), `Potential_V`,
`ElectricField_V_per_cm` (|E|) and `Electrons_cm-3` / `Holes_cm-3`. Coordinates in um:
1D x = depth from the anode; 2D x = lateral, y = depth. The VTU uses VTK_LINE cells
in 1D, VTK_QUAD cells on the 2D tensor grid and VTK_VERTEX otherwise; `bias_V` and
`validated_domain` are FieldData, the label / solver / units an XML comment. Carrier
provenance: devsim path = node solutions of the solved device; 2D engine path = the
engine dark state, equal to devsim's to rtol 1e-9 wherever the density exceeds
1e10 cm^-3 (default stack, hw 25 um, -20 V, `tests/test_export_mesh.py`); in the
depleted region both are below 1 cm^-3 (the engine underflows towards 0), i.e. zero
for any comparison. Results built without carriers export the three base columns.

### Layout / mask dimension file

Every lateral and passivation dimension used by the reference layout, the spec sheet
(`spec["layout"]`) and the cross-section is read from **`data/foundry/layout_params.toml`**
(nothing is hard-coded). One table per dimension and device:

```toml
[dosimeter.guard_ring_gap]
value = 150.0            # > 0; counts: integer >= 0
unit = "um"              # mm | um | nm, or "count"
class = "PLACEHOLDER"    # DOCUMENTED | MEASURED_FROM_FIGURE | PLACEHOLDER
source = "PLACEHOLDER: ..."
uncertainty = "unknown"  # >= 0 in `unit`, or "unknown"
note = "p+ edge to inner edge of the first ring"   # optional (also doc_item)
```

Devices `dosimeter` and `microdosimeter`; dimensions `active_side`, `die_side`,
`guard_ring_count`, `guard_ring_gap` (p+ edge to the first ring), `guard_ring_width`,
`guard_ring_spacing` (used only with two or more rings), `pad_size` (square corner pad,
placed between the outermost ring and the die edge), `dicing_street` (half drawn
around the die), `oxide_thickness`, optional `passivation_thickness`. The file is
seeded from `deliverables/DIMENSIONS-provenance.md`: the 2 x 2 / 5 x 5 mm active areas
are DOCUMENTED, the 2.7 +- 0.15 mm dosimeter die MEASURED_FROM_FIGURE; ring count,
widths and gaps, pad, dicing street, oxide thickness and the 5 x 5 die are PLACEHOLDER
drawing values, **not foundry dimensions**. Replace them with the mask data and set
`class = "DOCUMENTED"` plus a `source` naming the mask set.

- Other file: `ex.load_layout_params(path)`, the `layout_params=` argument of
  `build_spec` / `spec_to_pdf` / `cross_section_svg` / `cross_section_png` /
  `layout_svg` / `layout_gds` (path, mapping or `LayoutParams`), the environment
  variable `ETNA_LAYOUT_PARAMS`, or the upload box on the Export page.
- Validation (`LayoutParamsError`, the message names the field, e.g.
  `dosimeter.guard_ring_width.value: must be > 0`): schema version, required keys, units,
  positive values, integer counts, known classes, die larger than the active area, rings
  inside the die, pad between the outermost ring and the die edge.
- Placeholders: while a device's table has any PLACEHOLDER dimension (even an unused
  one, e.g. the ring spacing with one ring), its GDSII top
  cell ends in `REFERENCE-NOT-FAB` (`REFERENCE` once none is left), the placeholder names
  are written as a layout label, quoted in the spec sheet and the cross-section, shown as
  a warning on the Export page, and every export raises a
  `PlaceholderDimensionWarning`. Every layout carries "reference layout, not
  TCAD-validated".
- The spec sheet / cross-section pick the device whose active square equals
  `config.area_cm2` (0.04 cm^2 -> dosimeter, 0.25 cm^2 -> microdosimeter); other areas
  get no layout dimensions, stated as such.

Animations (`etna.viz`, needs `viz`): `etna.viz.animation_data(cfg)`, then
`render_animation(data, scene, out_dir)` for `scene in etna.viz.SCENES`
(`StackScene`, `DepletionScene`). Check `etna.viz.manim_available()` first.

---

## 11. Microdosimetric sweep

`etna.core.optimization.microdosimetric_sweep(...)` runs a combinatorial grid of
half-width × epi thickness × N_D_bulk × bias. Each point gets a lateral CCE scan, a
validity gate (|V| >= 1.2 V_fd, and the legacy contact wall) and a ranking. It uses the 2D
engine where eligible (`solver="auto"`) and devsim otherwise; failing configs become NaN
rows with `failure_reason`. `run_all(cfg, sections=("sweep",))` runs a small,
engine-only version: half-widths 25/50 µm × biases 20/40/60 V on the sidebar epi and
doping. If any point falls outside the engine box, the section is skipped rather than
switched to a slow devsim run.

Caveats:

- **Extrapolation.** The graded transition (Lb 1.14 µm, N_end 4.5e16) was fitted on the
  10 µm epi. The sweep applies it unchanged at other epi thicknesses, and labels such rows
  `transition Lb/N_end fitted on the 10 um epi, applied unchanged at this epi thickness:
model extrapolation` (`TRANSITION_EXTRAPOLATION_NOTE`, `etna/api/run_all.py`).
- The modelled device is laterally uniform, so `edge_center_ratio` is 1 by construction.
  It is not a uniformity result.
- CCE is relative to in-epi generation.

```python
from etna import DeviceConfig, run_all

rep = run_all(DeviceConfig(), sections=("sweep",))
sec = rep.get("sweep")
print(sec.status, sec.solver, sec.result.metadata["summary"])
```

---

## 12. Notebooks

`notebooks/` holds 21 notebooks, committed **without outputs**. `.gitattributes` routes
`*.ipynb` through nbstripout, and CI checks `nbstripout --verify notebooks/*.ipynb`. Run
them from `notebooks/` (`uv run jupyter lab`), so that their figures land in
`notebooks/figures/` (gitignored). Every notebook figure, whether written by `save_figure`
or by a direct `savefig`, goes to `$ETNA_FIGURE_DIR` or `./figures` of the directory the
kernel started in; notebooks 15-20 resolve it before their `os.chdir('..')`.

Notebooks written before 2026-10-18 may state numbers for older defaults; the model they
run is whatever `DeviceConfig()` / the core defaults are now, unless they pin a preset.

| notebook                       | topic                                          | notes                                                                                                                                                   |
| ------------------------------ | ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 01_phase1_validation           | material parameters, analytical electrostatics |                                                                                                                                                         |
| 02_electrical_characterization | C-V, Mott-Schottky, doping profile, I-V        |                                                                                                                                                         |
| 03_charge_collection           | CCE vs bias, generation profiles, Hecht        | generator `scripts/create_notebook_03.py`                                                                                                               |
| 04_flash_recombination         | FLASH Auger/SRH CCE vs dose rate               | sensitivity bound only                                                                                                                                  |
| 05_dark_current_vs_fluence     | dark current vs proton fluence                 | TAT on the batched engine                                                                                                                               |
| 05b_flash_parametric_studies   | FLASH parametric sweeps                        | runs its 60-condition sweep when the untracked cache is missing                                                                                         |
| 06_temperature_dependence      | material properties and CCE vs T               |                                                                                                                                                         |
| 07_dark_current                | TAT / SRH / surface terms                      |                                                                                                                                                         |
| 08_transient_flash             | FLASH pulse transients                         | generator `create_notebook_08.py`                                                                                                                       |
| 09_radiation_damage            | defect introduction, lifetimes, NIEL scaling   | markdown still says NIEL placeholders / kappa 0.35 (SR-NIEL since 2026-10-20)                                                                                                                                      |
| 10_cce_vs_fluence              | CCE degradation vs fluence                     |                                                                                                                                                         |
| 11_dark_current_cv_evolution   | leakage and C-V vs fluence                     |                                                                                                                                                         |
| 12_multi_defect_comparison     | 3-defect vs single-defect models               |                                                                                                                                                         |
| 13_parametric_optimization     | radiation-hardness sweep                       |                                                                                                                                                         |
| 14_validation                  | comparison with published irradiation data     |                                                                                                                                                         |
| 15_2d_electrostatics_cce       | 2D field and CCE of the microdosimeter         | generator `create_notebook_15.py`                                                                                                                       |
| 16_single_particle_cce         | single-particle 2D transient CCE, CCE(LET)     | generator `create_notebook_16.py`; loads `data/cce_let_table_*.json`; `ETNA_REBUILD_CCE_TABLES=1` rebuilds them (~40 min) into `$ETNA_FIGURE_DIR/data/` |
| 17_mc_coupling                 | MC energy deposition to pulse height           | generator executes it in place                                                                                                                          |
| 18_microdosimetric_spectra     | lineal energy spectra, y_F / y_D               | kappa placeholders; generator executes it                                                                                                               |
| 19_alternative_structures      | mesa, 3D electrode, ΔE-E, guard ring           | pins the previous-default doping on the legacy geometry (PERF-06)                                                                                       |
| 20_feasibility_report          | microdosimeter design study                    | outside-validity rows reported                                                                                                                          |

Every notebook calls `set_devsim_quiet(True)` (`etna.core.devsim_lock`, exported from
`etna`) in its first code cell: devsim's Newton log of the direct core calls is dropped
(it used to stream ~1e5 lines per notebook into the outputs); `ETNA_DEVSIM_VERBOSE=1`
brings it back. Use it in your own scripts too.

Run times measured on 2026-10-18 (Apple M4, executed into a scratch copy, the machine shared
with other jobs): 01 3 s, 02 17 s, 03 15 s, 04 5 s, 05 6.5 s, 05b 5.3 s (full 60-condition sweep, on
the engine; it still prints its old "~80 s" estimate), 06 29 s, 07 5.6 s, 08 21 s, 09 4 s, 10 11 s, 11 34 s, 12 27 s, 13 17 s, 14 18 s, 15 60 s, 16 71 s (tables
loaded; one 2D ion-track transient remains), 17 4 s, 18 3 s, 19 193 s, 20 123 s. Notebooks 05, 05b, 06, 07, 20
were re-timed on 2026-10-19 after the dark-current / FLASH engine port, back to back with the
previous code (scratch copies, same session): 05 35.9 -> 6.5 s, 05b 91.2 -> 5.3 s, 06 50.6 ->
28.5 s (its C-V / CCE-vs-T sweeps are still devsim), 07 27.2 -> 5.6 s (its
dark_current_sweep on a user-built devsim device stays devsim), 20 271.6 -> 123.1 s (the
2D dark current; the rest is the 2D CCE scans). Rebuilding the
notebook-16 tables (`ETNA_REBUILD_CCE_TABLES=1`) takes about 40 min (100 um ~15 min, 300 um ~26 min; one device per table, linear-response synthesis of the non-simulated LETs) and writes them to
`$ETNA_FIGURE_DIR/data/`, never over `data/`. Notebook 17 rewrites `data/synthetic_mc_events.csv` (bit-identical, seeded) relative to the repo root (they `chdir('..')`):
run them from a scratch copy unless you mean to regenerate those tracked files. The test suite and CI do not
execute notebooks (only the no-outputs and generator checks run);
`scripts/regression_sweep_v3_notebooks.py` is a historical manual tool for the v3
deliverables.

**Regeneration.** Where a generator exists it is the source of truth; otherwise the
notebook is. See the table in `scripts/README.md`. After editing either side run
`uv run python scripts/check_notebook_generators.py` (exit 1 on drift; CI runs it).
Generators run from `scripts/` and write `../notebooks/`. Generators 17 and 18 also
execute the notebook, so strip the outputs before committing.

**Rule:** never commit outputs (`uv run nbstripout notebooks/X.ipynb`), and never execute
notebooks from the repo root into the tracked `figures/`.

---

## 13. `lab/` and how provenance is recorded

`lab/` holds research spikes. Each has `notes/LOG.md`, with the commands run and the
numbers obtained. Nothing in `etna/`, `app/` or `tests/` imports from `lab/`, but
`tests/test_measured_cv_preset.py` reads the digitised C-V in `lab/cv_gap/notes/`, and
calibrated constants in `etna/` cite these logs (`lab/README.md`).

| directory               | question                                                                                             |
| ----------------------- | ---------------------------------------------------------------------------------------------------- |
| `lab/cv_gap`            | digitisation of the measured C-V (Photons Fig. 4b); which doping forms fit it; the 19-knot inversion |
| `lab/buffer_profile`    | the graded n-/n+ transition as a model extension; engine parity                                      |
| `lab/default_flip`      | bit-for-bit reproducibility of the presets across the 2026-10-18 default flip                        |
| `lab/tabulated_profile` | the tabulated preset (WP-E)                                                                          |
| `lab/finite_implant`    | is the -30 V residual lateral fringe capacitance?                                                    |
| `lab/microdosimeter`    | the microdosimeter preset; why its alpha CCE shape is not reproduced                                 |
| `lab/engine_box`        | widening the engine T and epi boxes; devsim ramp steps                                               |
| `lab/solver2d`          | the 2D numba engine spike (promoted to `etna/core/engine2d/`)                                        |
| GPU spike (archived)    | tag `archive/gpu-solver-lab`: numba beat the GPU for 1D batches                                      |

Where decisions and evidence live:

- `.issues/`: issues and work packages. `.issues/perf/PERF-0x-*.md` covers performance
  and the default stack, `.issues/audit-v7/` the latest audit (one file per finding,
  indexed in `*-INDEX.md`), and `.issues/TECH-DEBT-2026-10-18.md` the open-debt inventory.
- `.planning/`: physics audits (`PHYSICS_AUDIT_v*.md`, `PHYSICS_REVIEW_v6.md`), the
  literature check (`VALIDATION_LITERATURE.md`) and milestone audits.
- `deliverables/`: foundry-facing specs (DESIGN-1 dosimetry p-n, DESIGN-2 microdosimeter
  geometry; Markdown is the source and the PDFs were regenerated from it on 2026-10-18).
  `DIMENSIONS-provenance.md` classes every layout dimension as DOCUMENTED,
  MEASURED-FROM-FIGURE or PLACEHOLDER.
- The data behind a figure or number is regenerated by the script named next to it, e.g.
  `scripts/validation/cv_vs_bruzzi2023.py` produces `figures/cv_validation_bruzzi2023.png`.

---

## 14. Tests

```bash
uv run pytest -n auto -q                          # default run (the local gate)
uv run pytest -n auto -q -m exhaustive            # the deselected deep variants
uv run pytest -q -m "" tests/test_api_cv.py       # everything in one file
uv run pytest -n0 -x tests/test_api_cce.py        # serial, e.g. with --pdb
uv run pytest -n auto -q --durations=20           # slowest tests
```

- 92 modules and 1092 tests. `pytest.ini` deselects the `exhaustive` marker (7 redundant
  deep variants, each implied by a kept test), so the default run collects 1085.
  The `slow` marker labels the devsim integration tests but does not deselect them.
- `pytest.ini` adds `--dist loadgroup --no-loadscope-reorder` (inert without `-n`), so
  `-n auto` alone uses the tuned scheduler. Do not pass `--dist loadfile` any more (it
  overrides the tuned one), and use `-n0` rather than `-p no:xdist` (the latter makes
  the `--dist` options unknown).
- **Scheduling** (`tests/_xdist_schedule.py`): each test is tagged with an
  `xdist_group`: its whole file if it uses a module-scoped fixture, alone if it is on
  the `SPLIT_PER_TEST` list (heavy independent tests), else its class, else its file.
  The heaviest groups (`HEAVY`, rough CPU seconds) are dispatched first. A stale
  `HEAVY` entry only costs balance. Add a new 30 s+ test there if it becomes the tail.
- **devsim leak guard** (`tests/conftest.py`): `devsim.solve` solves _every_ device alive
  in the process, so one device left behind by a test made each later devsim solve in
  that worker ~20x slower (1D `cce_vs_bias` 10 s -> 208 s with one leaked 2D device).
  This showed up as random 100-400 s tests depending on the order. An autouse fixture
  now deletes devices created during a test; devices owned by module- or class-scoped
  fixtures survive. About 45 tests still leave devices (mostly `test_poisson.py`,
  `test_device2d.py`, `test_material.py`, `test_dark_current.py` helpers). List them
  with `ETNA_REPORT_DEVSIM_LEAKS=/tmp/leaks.txt uv run pytest -n auto -q`.
- **Timing** (10-core M4, 4 performance + 6 efficiency cores; the machine was shared
  with other jobs for every number here). 2026-10-18, before: `--dist loadfile`,
  502 s wall, 2173 s user CPU. After: 248 s wall, 1682 s user CPU (load average 4-11).
  The floor is the slowest single tests: run serially they take 87 s
  (`test_device2d_convergence.py` -50 V at 300 µm), 69 s (`test_devsim_quiet` direct
  core calls), 62 s (single-particle CCE vs bias), 46 s and 39 s. Inside a 10-worker
  run they take ~2x longer, on an efficiency core or with contended memory, so
  ~3-4 minutes is the practical wall time.
- numba kernels use `cache=True` and cache to `etna/**/__pycache__`. A cold cache costs
  only ~1-2 s per worker (`test_batched_engine.py`: 7.9 s cold vs 5.8 s warm), so
  there is no warm-up fixture.
- Threads: `tests/conftest.py` pins numba/BLAS to one thread per xdist worker (PERF-02).
- Families: `test_api_*` (facades), `test_app_*` (Streamlit `AppTest`, mocked and
  un-mocked smoke tests per page), engine parity (`test_batched_engine.py`,
  `test_api_engine_switch.py`, `test_engine2d*.py`, `test_nplus_engine*.py`,
  `test_engine_T_box.py`, `test_temperature_engine.py`, `test_fluence_engine.py`,
  `test_dark_current_engine.py`, `test_dark_current_2d_engine.py`, `test_flash_engine.py`),
  presets (`test_default_stack.py`, `test_nplus_calibration.py`,
  `test_measured_cv_preset.py`, `test_measured_cv_tabulated.py`,
  `test_microdosimeter_preset.py`) and the core physics modules.
- CI (`.github/workflows/ci.yml`) runs the lock check, ruff (`ruff check .` with
  `ruff.toml`: pyflakes F, E9, W6 over etna/, app/, tests/, scripts/), notebook hygiene, the generator check and the default run with
  `-n 4`. It runs only on pushes to `main`, pull requests and manual dispatch, not on
  every branch push, because the local run is the gate. The numba cache is kept in
  `actions/cache`, keyed on the exact hash of `etna/**/*.py` + `uv.lock`, with source
  mtimes pinned so numba accepts it. The test timeout is 45 min and has not been
  re-measured since the leak fix (it was 37-52 min before). A nightly job runs
  `-m exhaustive` and a `pip-audit` of the locked set.

---

## 15. Development workflow

- **Changing an engine or a solver path:** run the parity tests listed above before and
  after. Any widening of a box (`_ENGINE_BOX`, `ENGINE2D_*_BOX`) needs a parity bench
  against devsim at the new corners, recorded in a `lab/*/notes/LOG.md` (see
  `lab/engine_box/notes/LOG.md` for the pattern) and pinned by a test.
- **Never change a calibrated constant without before/after numbers.** This covers
  `DEFAULT_STACK`, `N_t`, mobility and lifetime parameters and the presets. Record the
  numbers (C-V points, V_fd, CCE, dark current, compensation fluence) in the issue or lab
  log. The previous defaults stay available as presets, and `tests/test_default_stack.py`
  / `tests/test_nplus_calibration.py` pin them bit for bit.
- **Dependencies:** edit `pyproject.toml`, then `uv lock`, and commit `uv.lock`.
  `uv lock --check` must pass (CI). Audit with
  `uv export --frozen --no-hashes --all-extras --no-emit-project > /tmp/r.txt && uvx pip-audit --no-deps --disable-pip -r /tmp/r.txt`.
- **Version:** bump `etna/_version.py` only.
- **Notebooks:** see section 12. Outputs are stripped, and generators must stay in sync.
- **Provenance:** a new finding gets an `.issues/` file or a `lab/` log entry with the
  command that reproduces it. Comments in code cite those files.

---

## 16. Troubleshooting

| symptom                                                            | cause and fix                                                                                                                                                                                                                                                      |
| ------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| A devsim device you built with `etna.core.*` disappeared           | devsim-path facades (`run_cv` fallback, 1D `run_field`, `run_dark_current`) call `reset_devsim_fully()` (`etna/core/devsim_reset.py`), which deletes **every** devsim device in the process. Build your own devices after the facade calls, or in another process. |
| results that depend on a previous cell (after an interrupted call) | a stale device from an interrupted call. Restart the kernel, or call `etna.core.devsim_reset.reset_devsim_fully()`. `run_all` asserts that no device leaks between sections.                                                                                       |
| `metadata["truncated"] == True`, shorter `x`                       | a devsim convergence wall at deep reverse bias. For the legacy model it is the cathode-on-epi artefact (about -64 V): use the default stack, or stay above it.                                                                                                     |
| `OutsideModelValidity`                                             | the requested point is outside the model's validity (legacy 2D past the contact wall; `run_cce_tilted_alpha` outside the engine box). The message says what to change.                                                                                             |
| NaN points in a damage sweep                                       | the fluence is at or past `compensation_fluence`; see `metadata["reason"]`.                                                                                                                                                                                        |
| First call is slow (about 1 s), later calls about 0.1 s            | numba compilation. The batched-engine kernels use `cache=True`, so later processes load the compiled code from `__pycache__`. Changing numba or the source recompiles.                                                                                             |
| A C-V or CCE call takes seconds instead of milliseconds            | the config is outside the engine box (`validated_domain_reason(cfg)`), or the engine fell back (`metadata["fallback_reason"]`).                                                                                                                                    |
| `ModuleNotFoundError` (plotly, streamlit) under `uv run pytest`    | `.venv` was synced without `--extra dev`, so a system pytest ran. Run `uv sync --extra dev`.                                                                                                                                                                       |
| Export page: "Animations are unavailable" / no GDSII button        | install `--extra viz` (macOS: `brew install pkgconf cairo pango` first). Note that a later plain `uv sync --extra dev` removes it again.                                                                                                                           |
| A notebook overwrote `figures/*.png`                               | it ran from the repo root. Run from `notebooks/`, or set `ETNA_FIGURE_DIR`, then `git checkout figures/`.                                                                                                                                                          |
| You want devsim's solver output                                    | `ETNA_DEVSIM_VERBOSE=1`.                                                                                                                                                                                                                                           |
