# ETNA usage guide

This is the single usage document for ETNA (Epitaxial Tcad Numerical Analysis). It is
written for a physicist or engineer who knows semiconductor radiation detectors but has
not read this code. Every Python block below was run against this repository (version
5.0.0; re-run 2026-10-21); the numbers quoted come from those runs on an Apple M4 (10
cores) and will differ slightly on other machines; timings are dated where given. File
paths are relative to the repository root.

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
17. [Measurement inputs (`etna.measurements`)](#17-measurement-inputs-etnameasurements)
18. [Physical-model library (opt-in): high-field mobility, anisotropy, avalanche and breakdown](#18-physical-model-library-opt-in-high-field-mobility-anisotropy-avalanche-and-breakdown)
19. [Single-ion track engine (opt-in, `etna.core.track`)](#19-single-ion-track-engine-opt-in-etnacoretrack)

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
    a different 30 µm wafer) checks the solver scale only: W +3.1 % at -40 V, +1.0 % with
    V_bi matched (paper W to 2 significant figures) (`.planning/VALIDATION_LITERATURE.md`). It does not validate this device.
  - **Not reproduced / open**. The microdosimeter's measured alpha CCE-vs-bias shape is
    not reproduced by any drift-diffusion variant (section 8). FLASH plasma recombination
    is not implemented. The NIEL hardness factors are SR-NIEL data, by default (since
    2026-10-21) in the `sic-damage` convention built on the NJOY SiC neutron NIEL
    (r_n = 1.069; `data/niel/README.md`); the literature introduction rates remain
    unvalidated and the default damage model's removal is calibrated on other devices'
    C-V at one temperature. The tissue-equivalence kappa(E) uses real NIST
    PSTAR/ASTAR data (SiC by Bragg additivity, `data/nist_star/`).
    SIMS is unavailable, so the epi/n+ position comes from the C-V fit.

---

## 2. What to trust

| Quantity                                                  | Status                                                      | Basis                                                                                                                                                                      | Where                                                                                                                       |
| --------------------------------------------------------- | ----------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| Engine results (C, W, CCE, potential) vs devsim           | **verified** (code-to-code)                                 | 1D: dC/C <= ~3e-8, dW <= 1e-17 cm, dCCE <= ~1e-9 (9e-10 at the tabulated-box corners) to -100 V; 2D: dpsi ~1e-13 V                                                         | `tests/test_batched_engine.py`, `tests/test_api_engine_switch.py`, `tests/test_engine2d*.py`, `lab/engine_box/notes/LOG.md` |
| C-V / eps*A/C of the default device                       | **calibration** (fit)                                       | default doping fitted to the digitised Photons Fig. 4b: plateau rms 0.013 µm, knee rms 0.24 µm; 1.655 / 9.519 / 9.750 µm at 0 / -10 / -30 V vs paper 1.68 / 9.52 / 9.73    | `lab/cv_gap/notes/LOG.md`, `etna/core/device.py` (`DEFAULT_STACK`)                                                          |
| Doping profile N_D(x)                                     | **inferred**, not measured                                  | inversion of one C-V curve; bulk 3.4e13 is below the documented 0.5-1e14; no SIMS                                                                                          | section 8                                                                                                                   |
| Full-depletion voltage                                    | **definition-dependent**                                    | analytic V_fd 2.74 V, reach-through 3.10 V, sweep gate 1.2 V_fd = 3.28 V (default stack)                                                                                   | `etna.export.labels.VFD_DEFINITION`                                                                                         |
| CCE vs bias (1D, alpha in epi)                            | **model** (low injection, steady state)                     | engine verified vs devsim; no measured CCE reproduced for the dosimeter                                                                                                    | `run_cce` docstring                                                                                                         |
| Microdosimeter alpha CCE shape (Fig. 4b, Phys. Med. 2025) | **not reproduced** (calibration only) | 1D: flat to 0.7 % from 0 V; tilted track + plasma (2.5D) with library lifetimes: rms 0.033. The −20..−50 V tail needs τ_n+τ_p ≈ 23 ns (band 17-42), ~100x below good epi; 0..−10 V not fitted | `run_cce_tilted_alpha` docstring, `lab/microdosimeter/notes/LOG.md` §9, `lab/track_plasma/notes/LOG.md` §§10-12 |
| Dark current                                              | **calibration** (budget estimate)                           | one N_t fitted to < 18 pA at -60 V, 4 mm²; bias shape f_depl heuristic x Hurkx (1+Gamma), temperature shape assumed                                                        | `run_dark_current`, DESIGN-1                                                                                                |
| CCE vs proton fluence (default model, since 2026-10-21)   | **model**: removal **calibrated**, CCE **predicted**        | kinetic trap preset: literature levels (Burin/Capan) + direct removal d **calibrated** to 3 published 252.7 MeV C-V slopes (other devices, one temperature: removal mechanism NOT identified, a slow deep acceptor fits equally); NIEL `sic-damage` (computed r_n, spectrum-dependent 1.02-1.20); CCE vs fluence and beyond compensation **predicted**, not measured | `run_radiation_damage` docstring, section 6 "Defaults since 2026-10-21" |
| Compensation fluence                                      | **model** (calibrated removal x NIEL)                       | 1.80e12 p/cm² at 62 MeV (default stack); 4.36e12 legacy scalar, 3.48e12 legacy si-neq                                                                                   | `compensation_fluence`                                                                                                      |
| FLASH dose-rate CCE                                       | **numerical sensitivity bound**                             | no plasma or high-injection physics                                                                                                                                        | `run_flash_recombination`, `run_transient` docstrings                                                                       |
| Tissue-equivalence kappa(E)                               | **data** (NIST) + approximation                             | NIST PSTAR/ASTAR (2026-10-20); SiC by Bragg additivity (error ~few % for p > 1 MeV, unverified estimate); LET-ratio convention, 0.43 @1 MeV -> 0.375 @150 MeV (p -> water) | `data/nist_star/README.md`, `compute_kappa_table`                                                                           |
| 2D field maps                                             | **verified** vs devsim inside the engine2d box              | laterally uniform device; no edge termination                                                                                                                              | `tests/test_engine2d*.py`, `tests/test_field_engine_default_100V.py`                                                        |
| Reference layout (GDS/SVG)                                | **reference only, not TCAD-validated, not for fabrication** | dimensions from `data/foundry/layout_params.toml`; ring, pad, dicing, oxide and 5x5 die are PLACEHOLDER drawing values until the mask data is entered                      | `deliverables/DIMENSIONS-provenance.md`, section 10                                                                         |

Summary: the numerics are well verified. The physics is calibrated on one C-V curve and
one leakage bound. Everything else is a model prediction that has not been checked
against a measurement of this device.

---

## 3. Installation

Requirements: Python >= 3.13 and [uv](https://docs.astral.sh/uv/). devsim ships as a
wheel; it looks for BLAS/LAPACK at import (on macOS it uses the system Accelerate/LAPACK;
on Linux install OpenBLAS/LAPACK, e.g. `apt-get install libopenblas-dev liblapack-dev`).

```bash
git clone https://github.com/ngcex/sic-tcad-detector.git
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

### Installing without a checkout

ETNA is distributed from GitHub only (no PyPI release; the name `etna` is taken there).
Into an existing uv-managed environment:

```bash
uv pip install "git+https://github.com/ngcex/sic-tcad-detector"                  # library
uv pip install "etna[viz] @ git+https://github.com/ngcex/sic-tcad-detector"      # + gdstk / manim
uv pip install "git+https://github.com/ngcex/sic-tcad-detector@<commit-or-tag>"  # pinned
```

Or build a wheel from a checkout with `uv build` (writes `dist/etna-*.whl` and the
sdist), then `uv pip install dist/etna-*.whl`.

- **The runtime data is in the wheel.** The hatch `force-include` table in
  `pyproject.toml` copies the runtime subset of `data/` to `etna/_data/`: the NIST
  PSTAR/ASTAR CSVs, the SR-NIEL CSVs, the legacy stopping-power tables, the CCE(LET)
  tables, the digitised validation points, `foundry/layout_params.toml`, the measurement
  templates and the provenance READMEs. Raw provenance (`data/nist_star/raw/`, the
  SR-NIEL `.html` pages and fetch script) and notebook inputs (`synthetic_mc_events.csv`)
  stay in the sdist and the repository. `etna/_paths.py` resolves every data file in
  this order: `$ETNA_DATA_DIR` if set, then the packaged `etna/_data/`, then the
  checkout's `data/` (source checkout, editable install). Scripts, notebooks and `lab/`
  therefore keep reading `data/` in the repository. `tests/test_packaging.py` fails when
  a new runtime CSV or JSON in `data/` is missing from the table.
- **The layout parameter file stays user-editable.** In a wheel the default is package
  data. Copy it out and point `ETNA_LAYOUT_PARAMS` at the copy (section 10):

  ```bash
  cp "$(python -c 'import etna.export as e; print(e.LAYOUT_PARAMS_PATH)')" my_layout.toml
  export ETNA_LAYOUT_PARAMS=$PWD/my_layout.toml
  ```

- **Not in the wheel:** the Streamlit app (`app/`), notebooks, scripts, tests and `lab/`.
  Use a checkout for those.
- The CI job `wheel-smoke` (`.github/workflows/ci.yml`) does the following:
  - builds the sdist and the wheel, and checks the metadata with `twine check`
  - installs the wheel into a clean venv outside the checkout
  - runs `scripts/wheel_smoke.py --app` with `python -I`

  The script checks that `etna` is imported from site-packages with the packaged data,
  then runs `run_cv`, `run_cce`, `run_all`, the SiC->water kappa, `w_creep()`, and the
  layout params with SVG/GDS and the spec PDF. It also runs the app, copied next to the
  installed package: the entry script, the measurement templates, the LET tables and
  `/_stcore/health`. To run it locally:

  ```bash
  uv build -o /tmp/dist && uv venv /tmp/v -p 3.13
  VIRTUAL_ENV=/tmp/v uv pip install /tmp/dist/etna-*.whl gdstk
  cd /tmp && /tmp/v/bin/python -I "$OLDPWD/scripts/wheel_smoke.py" --app
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

report = run_all(cfg)                     # every section in one call, < 1 s with a warm numba cache
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

`DeviceConfig` validates on construction: it rejects non-positive or non-finite sizes,
doping and area, N_D >= N_A, T outside 100-600 K and an inconsistent transition, with a
`ValueError` that names the field. The same check runs again when a facade (`run_*`,
`build_device`) receives the config, so a field changed in place after construction
(`cfg.substrate_thickness_um = 0.0`) gives that `ValueError` and never a raw devsim
error. Set `half_width_um` (µm) for a 2D device (lateral half-width; depth is the 2D
`y`).

The epi doping profile is a strict contract. The field is typed
`DopingProfile = Literal["graded", "uniform"]`; the allowed values are
`etna.api.device.DOPING_PROFILES`.

- `doping_profile` is exactly `"graded"` or `"uniform"`. Matching is case-sensitive:
  `"Graded"`, `"step"` or a typo raises a `ValueError` that lists the allowed values.
  The config is written into every CSV and figure header, so each physical model has
  one spelling only.
- `"uniform"` requires `N_D`. There is no hidden default, and leaving it `None` raises.
- `"graded"` requires `N_D=None`. The graded epi is set by `N_D_junction`, `N_D_bulk`
  and `L_transition_um` (or by `donor_table`), and the builders never read `N_D`, so an
  explicit value raises instead of being ignored silently.

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

| dopant        | levels (dE0, alpha)                                                            | g   | source                                                 |
| ------------- | ------------------------------------------------------------------------------ | --- | ------------------------------------------------------ |
| N (donor)     | hexagonal 70.9 meV, 3.38e-5 meV cm; cubic 123.7 meV, 4.65e-5 meV cm (1/2 each) | 2   | Kagamihara et al., JAP 96, 5601 (2004); Burin Tab. XXV |
| Al (acceptor) | 220 meV, 2.32e-5 meV cm                                                        | 4   | Achatz et al., APL 92, 072103 (2008); Burin Tab. XXIV  |

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

| quantity                         | legacy                           | carrier-dependent                |
| -------------------------------- | -------------------------------- | -------------------------------- |
| N_A_ionized (p+)                 | 1.322e18                         | 2.913e17                         |
| V_bi (psi(cathode) - psi(anode)) | 3.198 V                          | 3.142 V                          |
| eps*A/C at 0 / -2 / -10 / -30 V  | 1.655 / 6.939 / 9.519 / 9.750 µm | 1.619 / 6.867 / 9.519 / 9.750 µm |
| CCE at -10 V                     | 0.999473                         | 0.999466                         |
| reach-through                    | 3.10 V                           | 3.15 V                           |

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
`DeviceConfig` first. Timings (measured 2026-10-21, Apple M4, load average 1-4) are
the second call of the same facade on the default device in one process, numba disk
cache loaded. The first call in a fresh process adds the cache load (run_cv 0.13 s
instead of 7 ms); with an empty numba cache it compiles instead (about 0.4-0.6 s per
kernel family: run_cv 0.56 s, run_cce 0.41 s, 2D field 0.67 s).

### `run_cv(config, v_start=0.0, v_stop=-50.0, n_points=40)`

Returns small-signal C-V. `x` is bias (V), `y` is C = dQ/dV (F), the observable of a
measured C-V. Metadata:

- `W_from_C` (eps*A/C, cm) and `depletion_widths` (the n/N_D = 0.5 depletion edge, cm;
  not eps*A/C near full depletion)
- `one_over_C_squared`, `area_cm2`
- `truncated`, `requested_v_stop`
- `reach_through_V`, `reach_through_known`, `reach_through_method`
- provenance (section 7)

1D only. Takes about 10 ms (40 points) on the batched engine.

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

1D only (a 2D config raises `NotImplementedError`). Takes about 30-40 ms (30 points) on
the engine.

```python
from etna import DeviceConfig, run_cce

r = run_cce(DeviceConfig(), v_start=-1.0, v_stop=-10.0, n_points=10)
print(r.metadata["solver"], [round(float(v), 4) for v in r.y[[0, 1, -1]]])  # batched [0.9719, 0.9937, 0.9994]
```

### `run_cce_tilted_alpha(config, v_start=0, v_stop=-60, n_points=61, theta_deg=65, alpha_range_um=18, dead_layer_um=0, shape="bragg", generation_in_p_plus=False)`

Models the microdosimeter's tilted-alpha test. It runs on the batched engine only and
raises `etna.api.simulation.OutsideModelValidity` outside the engine box. Metadata:

- `cce_vs_E0` (collected over E0/E_pair, which includes dead-layer and p+ losses)
- `projected_end_um`, `dead_layer_energy_fraction`, `generation_shape`
- `reach_through_V`

Since 2026-10-20 `shape="bragg"` is the NIST ASTAR curve in SiC (Si + C by Bragg
additivity, Bohr straggling, peak/entrance 2.49), depth-scaled to `alpha_range_um`
(NIST CSDA range: 18.50 um at 5.5 MeV); `"bragg_powerlaw"` is the former qualitative
curve. The 1D/2D alpha profiles (`alpha_depth_shape`) use the same data by default
(range 18.43 um at 5.486 MeV; the old 15 um flat box is `bragg_peak_ratio=None`).
The model **does not** reproduce the
paper's CCE-vs-bias shape (section 8).

```python
from etna import microdosimeter_config, run_cce_tilted_alpha

r = run_cce_tilted_alpha(microdosimeter_config())
print(r.metadata["projected_end_um"], r.metadata["reach_through_V"])  # ~7.61 um, 9.0 V
```

### `run_field(config, bias_V=-50.0, solver="auto")`

Returns the field and potential at one reverse bias.

- **1D**: `x` is depth (µm), `y` is |E| (V/cm). Metadata holds `potential` (V) and
  `net_doping` (cm^-3). 1D runs on the **batched engine** inside its box (the run_cv
  box, |bias| up to `cross_checked_to_V` = 40 V without n+ / 100 V with it, static
  ionization; about 1.5 ms, `metadata["solver"] == "batched"`): same mesh and keys as devsim,
  potential within 1.5e-14 V, |E| within 4.5e-12 of max E, Electrons / Holes within 2e-9
  relative (`tests/test_field_engine_1d.py`). Otherwise, or with `solver="devsim"`, devsim
  (about 50 ms). `solver="engine"` raises where no engine serves the call.
- **2D** (`half_width_um` set): `x` and `y` are **empty**. Use `mesh.x_coords`
  (lateral, cm), `mesh.y_coords` (depth, cm) and `mesh.node_values["ElectricField"]`.
  The 2D engine serves configs inside its box (about 0.3 s for a new geometry, mostly
  the mesh; about 10 ms for another bias on the same geometry); otherwise devsim runs.
- A legacy-model 2D bias past the contact wall raises `OutsideModelValidity`
  (section 8).

```python
from etna import DeviceConfig, run_field
from etna.api.results import peak_field_1d

cfg = DeviceConfig()
r = run_field(cfg, bias_V=-30.0)
peak, depth_um, excluded = peak_field_1d(r.x, r.y, cfg)   # skips the nm spikes at abrupt steps
print(r.metadata["solver"], f"{peak:.3g} V/cm at {depth_um:.2f} um", excluded)  # batched 9.06e+04 V/cm at 0.40 um True

r2 = run_field(DeviceConfig(half_width_um=50.0), bias_V=-30.0)
print(r2.metadata["solver"], r2.x.size, r2.mesh.x_coords.size, r2.mesh.y_coords.max() * 1e4)  # engine2d 0 17784 ~11.3 um
```

### `run_radiation_damage(config, fluences=None, V_bias=-40.0, proton_energy_MeV=5.6, *, particle="proton", niel_convention=None)`

Returns CCE vs proton fluence. `x` is fluence (p/cm²), `y` is CCE. Metadata:

- `valid` (bool per point) and `reason` (str per point)
- `damage_model` (`"charged_traps"` or `"carrier_removal"`), `damage_model_label` (the
  trap model's provenance label) and `damage_model_note` (why `"auto"` ran the legacy
  model, else empty)
- `compensation_fluence` (p/cm², of the model that ran), `legacy_compensation_fluence`
  (the scalar-removal reference) and, with traps, `trap_compensation_fluence`
- `particle`, `niel_convention`, `kappa`
- `solver` (the solver that ran the points: `batched` / `batched+traps` on the default
  stack, `devsim` outside the engine box, joined when mixed) and `solver_per_point`
- `V_bias`, `energy_MeV`

#### Defaults since 2026-10-21 (issues #5, #8) and how to get the old numbers

The maintainer changed two defaults; both old behaviours stay selectable by name:

| | default since 2026-10-21 | previous default (selectable) |
| --- | --- | --- |
| NIEL convention | `niel_convention=None` = `"sic-damage"` (`radiation_damage.DEFAULT_NIEL_CONVENTION`) | `niel_convention="si-neq"` |
| damage model | `DeviceConfig(radiation_traps="auto")`: the calibrated kinetic preset `trap_presets.radmanovac_calibrated` (direct donor removal) | `radiation_traps=False` or `"legacy"`: scalar donor removal, eta 5 /neq, points beyond compensation refused |

```python
import dataclasses
from etna import DeviceConfig, run_radiation_damage
old = run_radiation_damage(DeviceConfig(radiation_traps=False), fluences, -40.0, 62.0,
                           niel_convention="si-neq")   # bit-identical to the old default
new = run_radiation_damage(DeviceConfig(), fluences, -40.0, 62.0)
```

- **`radiation_traps` values (strict).** `"auto"` (default), `"legacy"` / `False`,
  `True` (Burin Table I, steady occupancy: the opt-in model of 2026-10-18), a preset
  name of `etna.core.trap_presets.PRESETS`, or a `RadiationTrapModel`. Anything else
  raises (`1`, `"yes"`, unknown names).
- **`"auto"` scope.** The charged-trap engine is batched-1D only. Where it cannot run
  (`incomplete_ionization=True`, T outside 250-400 K, legacy geometry beyond 40 V,
  `solver="devsim"` in the core function) `"auto"` runs the legacy model, logs a
  warning and says why in `damage_model_note`. Explicit trap models never fall back:
  they raise. `"auto"` is therefore not exclusive with `incomplete_ionization`; an
  explicit trap model is.
- **Zero fluence** stays on the legacy engine path: the pristine point is bit-identical
  to the legacy sweep. Non-irradiated simulations (`run_cv`, `run_cce`, ...) do not read
  `radiation_traps` and are unchanged.
- **Which paths use which model.** The trap model runs in `run_radiation_damage`,
  `run_all` ("damage"), the Radiation Damage page and the core
  `cce_vs_fluence(radiation_traps=...)`. The core `cce_vs_fluence` keeps
  `radiation_traps=None` (legacy) as its own default (it is the devsim-parity primitive).
  `cce_vs_bias_at_fluence`, the C-V / dark current at fluence (notebooks 05, 11, 12,
  13) and `RadiationDamageParams`-based helpers stay on the scalar model; they follow
  the new NIEL default.
- **Fluence grid.** With `fluences=None` (`default_fluence_grid`): trap model, 10
  log-spaced points from 0.01x its compensation fluence to max(3x it, 1e16 neq), all
  solved; legacy model, 8 points from 0.01x to 3x its compensation fluence (the last 2
  normally NaN with a reason).
- **Speed.** Default sweep (10 points, 62 MeV, -40 V, default stack): 1.2 s in a fresh
  process with the numba cache on disk, 0.8 s warm (`lab/radiation_traps/bench/
  default_flip.py`; the kinetic dark state is one batched implicit run for all points).
  The first run after a clean checkout compiles the kernels (~15 s, once).

Compensation fluence (`compensation_fluence(config, E)`: min N_D / (rate x kappa) of the
model that runs; trap model: rate = neutral-bulk removal = every acceptor level filled +
d = 7.4 + 4.69 = 12.09 /neq):

| device, energy | default (traps, sic-damage) | legacy, sic-damage | legacy, si-neq (old default) |
| --- | --- | --- | --- |
| default stack, 62 MeV | 1.80e12 | 4.36e12 | 3.48e12 |
| default stack, 5.6 MeV | 4.75e11 | 1.15e12 | 1.01e12 |
| tabulated preset, 62 MeV | 6.90e11 | 1.67e12 | 1.33e12 |

CCE at -40 V, 62 MeV, default stack, at fractions of each model's own compensation
fluence (`lab/radiation_traps/bench/default_flip.py`):

| x Phi_comp | 0.1 | 0.5 | 0.9 | 1.5 | 10 |
| --- | --- | --- | --- | --- | --- |
| old default (legacy, si-neq; Phi_comp 3.48e12) | 0.999618 | 0.998882 | 0.998294 | NaN | NaN |
| legacy, sic-damage (4.36e12) | 0.999618 | 0.998882 | 0.998294 | NaN | NaN |
| calibrated traps, si-neq, d 2.87 (1.70e12) | 0.999787 | 0.999572 | 0.999377 | 0.999113 | 0.996910 |
| **new default** (calibrated traps, sic-damage, d 4.69; 1.80e12) | 0.999798 | 0.999619 | 0.999456 | 0.999237 | 0.997448 |

Only the **new default** row is on the current c-axis electron mobility (re-measured
2026-10-21 after the mobility default changed); the other rows were computed on
`mobility="legacy_perp"` and move by up to about 1e-4 on the current default.

At the same absolute fluences (0.1 ... 10 x 3.48e12 p/cm²) the new default gives
0.999751 / 0.999414 / 0.999122 / 0.998745 / 0.995838 (old: 0.999618 / 0.998882 /
0.998294 / NaN / NaN). The new default is LESS damaging below compensation (the legacy
lifetime K_tau counts every level as a trap, the occupancy lifetimes only free ones) and
keeps solving beyond it. These CCE numbers are predictions (no measured CCE vs fluence
of this device).

#### NIEL convention and the calibrated removal (consistency algebra)

The two kinds of per-neq numbers in the chain transform differently:

1. **Neutron-anchored literature rates** (Burin et al. arXiv:2407.16710 Table I: Z1/2
   5.0, EH4 2.4, EH6/7 1.6 /neq; footnote b cites Hazdra et al., "Point defects in
   4H-SiC epilayers introduced by neutron irradiation" (title read on the arXiv page;
   that paper not opened). The rates are labelled in reactor 1-MeV neq. `sic-damage`
   is normalised by r_n exactly so that a fission-spectrum neutron has the SAME kappa
   in `si-neq` and `sic-damage` (`tests/test_neutron_niel.py`): these rates carry over
   unchanged. Per PROTON they now give r_p(E)/r_n fewer defects (x0.80 at 62 MeV, x0.75
   at 252.7 MeV). That is the intended physics correction of the decision, not
   bookkeeping.
2. **Proton-calibrated numbers** (the slow removal d, fitted to Radmanovac et al.
   arXiv:2510.11304 Table 2 slopes per p/cm² at 252.7 MeV). The data are per proton and
   convention independent; the model reads them as
   g_p = kappa_c(252.7) x [G_lit(kappa_c Phi_p) + d_c],
   with G_lit the per-neq removal that the (convention-independent) literature levels
   produce in the depleted window. Hence d_c is convention DEPENDENT and is RE-FITTED in
   each convention (`lab/radiation_traps/bench/convention_refit.py`, ~6 min):

   | convention | kappa(252.7) | d (/neq) | chi2 (2 dof) | model g per p (onsemi, W2, W4) | literature-level share per p |
   | --- | --- | --- | --- | --- | --- |
   | si-neq | 0.8514 | 2.873 +- 0.221 | 1.13 | 4.95, 4.81, 4.41 | 2.41, 2.26, 1.91 |
   | sic-damage | 0.6390 | 4.694 +- 0.303 | 1.03 | 4.86, 4.74, 4.47 | 1.80, 1.67, 1.43 |
   | sic-1mev | 0.6051 | 5.097 +- 0.324 | 1.04 | 4.84, 4.73, 4.48 | 1.70, 1.57, 1.35 |

   measured: 5.60 +- 0.81, 4.60 +- 0.33, 4.48 +- 0.26. The fitted per-proton slopes agree
   between conventions to 0.09-0.11 /p, well inside the data sigmas: the per-proton
   calibration is invariant, as it must be. Rescaling instead (d_si x kappa_si /
   kappa_c = 3.83 /neq in sic-damage) would under-predict g by ~0.6 /p, because G_lit
   does not scale.
3. **No double counting.** The preset records its convention
   (`RadiationTrapModel.niel_convention`); `resolve_traps("auto", convention)` builds it
   in the run's convention, and every fluence sweep (`cce_vs_fluence`,
   `KineticTrapSolver.from_config`, `calibrate_cv`) refuses a calibrated model under
   another convention (`check_trap_convention`, ValueError). The literature levels
   carry no convention (they are neutron-anchored).
4. **What is NOT invariant** (and is reported, not hidden): the share of the removal
   carried by the literature levels (above), the energy transfer from 252.7 MeV to other
   energies (kappa_c(62) / kappa_c(252.7) = 2.43 in sic-damage vs 2.28 in si-neq: the
   SiC vs Si NIEL energy dependence), hence the 62 MeV neutral-bulk removal per proton
   (18.8 /p sic-damage vs 20.0 /p si-neq) and the compensation fluence (1.80e12 vs
   1.70e12 p/cm², +6 %), and the CCE beyond compensation (fewer literature traps per
   proton, more non-trapping removal).
5. **Unconverted checks.** The flat-C-V brackets of Kroll (24 GeV/c, authors' Si kappa
   0.62) and Li (~80 MeV, kappa unstated) are labelled in Si neq by the authors and
   cannot be converted here (the proton tables stop at 1 GeV; Li's kappa unknown).
   They are checks, not fit inputs (`bench/mechanisms.py`), and remain in Si neq.

#### NIEL hardness factors

NIEL hardness factors (audit C-5 / B-07, since 2026-10-20): SR-NIEL-7 proton NIEL
(1-1000 MeV, 59 energies, log-log interpolation). Default `sic-damage`: kappa(E) =
NIEL_p,SiC(E) / (r_n NIEL_n,Si(1 MeV)): 5.89 at 5.6 MeV, 1.55 at 62, 0.859 at 150, 0.639
at 252.7 MeV. `si-neq` (the previous default): NIEL_p,Si(E) / NIEL_n,Si(1 MeV) over the
ASTM 95 MeV mb normalisation: 6.67 at 5.6 MeV, 2.75 at 30, 1.94 at 62, 1.82 at 70, 1.12
at 150 MeV (the old placeholders were 0.50 / 0.35 / 0.33 / 0.22 at 30-150 MeV and 5.6
MeV was silently clamped to 0.5). Energies outside 1-1000 MeV raise `ValueError`. The
SiC/Si NIEL ratio of protons r_p (0.94 at 5.6 MeV, 0.85 at 62 MeV:
`radiation_damage.niel_ratio_sic_si`) is below the fission-neutron r_n = 1.069, so
`si-neq` overstates the SiC damage per proton by r_n/r_p. At 5.6 MeV NIEL scaling is
also known to overpredict measured Si damage by tens of percent. Data and provenance:
`data/niel/` (raw SR-NIEL pages, CSVs, fetch script, README).

```python
import numpy as np
from etna import DeviceConfig, run_radiation_damage
from etna.api.simulation import compensation_fluence

cfg = DeviceConfig()
phi_c = compensation_fluence(cfg, 5.6)                       # 4.75e11 p/cm2
r = run_radiation_damage(cfg, fluences=np.geomspace(1e9, 10 * phi_c, 5), V_bias=-40.0)
print(f"{phi_c:.3g}", r.metadata["valid"].all(), np.round(r.y, 4))
```

#### Neutron NIEL and NIEL conventions (opt-in, `particle=`, `niel_convention=`)

The SiC neutron NIEL is computed (2026-10-21) with NJOY2016 HEATR on ENDF/B-VIII.0
(Lindhard/Robinson partition, E_d 21/21 eV; `scripts/build_neutron_niel.py`, data and
checks in `data/niel/README.md`, code `etna.core.neutron_niel`). Checks: Si at 1 MeV is
0.998 x the ASTM 95 MeV mb; the Si fission average is +0.7 % vs the Griffin
(SAND92-0094) damage function; an independent Python elastic integration matches HEATR
to 2e-4. Results: NIEL_n,SiC(1 MeV) = 2.299e-3 MeV cm^2/g; SiC/Si for a Watt fission
spectrum r_n = 1.069. No published SiC neutron NIEL value was found to compare with, and
r_n depends on the (unknown) calibration reactor spectrum: 1.02-1.20 for plausible
shapes.

`etna.core.radiation_damage.hardness_table(particle, convention)` gives the factor
(neq per particle); `run_radiation_damage`, `compensation_fluence`,
`trap_compensation_fluence` and `default_fluence_grid` take `particle=` and
`niel_convention=` (keywords), and the core functions take `niel_table=`:

| convention | kappa | 62 MeV proton |
| ---------- | ----- | ------------- |
| `si-neq` (default until 2026-10-21) | NIEL_x,Si(E) / NIEL_n,Si(1 MeV, ASTM) | 1.942 |
| `sic-damage` (**default**) | NIEL_x,SiC(E) / (r_n NIEL_n,Si(1 MeV)): same SiC damage as the Si-labelled reactor neq of the rates | 1.552 |
| `sic-1mev` | NIEL_x,SiC(E) / NIEL_n,SiC(1 MeV): SiC-native neq | 1.470 |

`particle="neutron"` reads the fluences as neutrons and the energy as the neutron energy
(bin-averaged NJOY table, 1/20 decade bins, ~1e-5 to 18.9 MeV, no extrapolation; 1 MeV:
0.953 in `si-neq`, because the resonant damage function is averaged over the bin; 14 MeV:
1.89, +4-5 % above Griffin). With the fission r_n, `si-neq` overstates the SiC damage
per proton by r_n/r_p (x1.13 at 5.6 MeV, x1.25 at 62, x1.30 at 150 MeV), which is the
factor `sic-damage` removes; `sic-damage` is the default since 2026-10-21. Bench:
`uv run python lab/neutron_niel/bench.py` (0.05 s).

```python
from etna import DeviceConfig, run_radiation_damage
from etna.api.simulation import compensation_fluence

cfg = DeviceConfig(radiation_traps="legacy")
compensation_fluence(cfg, 62.0, niel_convention="si-neq")        # 3.48e12 p/cm2
compensation_fluence(cfg, 62.0)                                  # 4.36e12 p/cm2 (sic-damage)
r = run_radiation_damage(cfg, proton_energy_MeV=14.0, particle="neutron")
r.metadata["kappa"], r.metadata["particle"]                      # 1.535, 'neutron'
```

#### Charged deep traps beyond compensation (`radiation_traps=True`: Burin Table I, steady)

`DeviceConfig(radiation_traps=True)` (or a `RadiationTrapModel`) replaces the
donor-removal scalar with the radiation-induced deep levels as charged species in
Poisson, as in the Burin et al. TCAD approach (arXiv:2407.16710). Real irradiated 4H-SiC
turns semi-insulating and keeps detecting. The legacy model refuses those points; this
one solves them. Since 2026-10-21 the DEFAULT (`"auto"`) is the calibrated kinetic
preset of the next subsections; `True` stays this steady Burin model. The numbers in this
subsection are in the `si-neq` convention (pass `niel_convention="si-neq"` to reproduce
them). Code: `etna.core.radiation_traps` and `PoissonBatch(traps=...)`.

- **Traps.** N_t,k = eta_k * Phi_neq on the epi nodes, from Burin Table I via
  `RadiationDamageParams`:
  - Z1/2: acceptor, Ec-0.67, eta 5.0.
  - EH4: acceptor, Ec-1.03, eta 2.4.
  - EH6/7: donor, Ec-1.60, eta 1.6, sigma_n 1.4e-14 (audit C-4, not Burin's 9e-12).
- **Donors.** N_D is not reduced. Removal emerges from the trap charge.
- **Occupancy.** Steady-state SRH, f = (c_n n + c_p p1) / (c_n (n + n1) + c_p (p + p1)),
  with p1 = n_i^2 / n1, so f is Fermi-Dirac at equilibrium. Degeneracy g = 1 (assumed).
  Z1/2 is a single level (Burin Table I).
- **Which state sets f.** f comes from the dark state at the applied bias, solved
  self-consistently by Newton. It is then frozen as fixed charge in the illuminated
  Gummel solve.
  - Reason: alpha pulses are ns long, while the trap time constants run from ms to
    1e14 s.
  - Check: CCE changes by < 1e-4 between G = 1e16 and 1e20.
- **Lifetimes.** Unchanged (`compute_K_tau`).
- **Mobility.** Taken at the pristine donors.
- **Solver.** Batched 1D engine only. `solver="devsim"` and `incomplete_ionization`
  are refused. An unconverged point is NaN with a reason: there is no silent legacy
  fallback.
- **Fluence grid.** The default grid reaches 1e16 neq. Metadata gives
  `damage_model="charged_traps"`.

Numbers (default stack, 62 MeV, kappa 1.942; `lab/radiation_traps/bench/validate.py`,
log in `lab/radiation_traps/notes/LOG.md`; computed with the perpendicular-to-c mobility
of that date, now `mobility="legacy_perp"`. With the c-axis electron default of
2026-10-21 the traps row at -40 V reads 0.999780 / 0.999273 / - / 0.99642 / 0.97484 / - /
0.35761 / 0.15297: the faster electrons collect more where trapping limits the collection,
by up to 2.6e-3; `tests/test_radiation_traps.py`):

| Phi_p (p/cm^2)    | 1e11     | 1e12     | 3e12     | 1e13  | 1e14  | 3e14  | 1e15  | 5e15  |
| ----------------- | -------- | -------- | -------- | ----- | ----- | ----- | ----- | ----- |
| CCE legacy, -40 V | 0.999777 | 0.999247 | 0.998347 | NaN   | NaN   | NaN   | NaN   | NaN   |
| CCE traps, -40 V  | 0.999777 | 0.999248 | 0.998345 | 0.996 | 0.972 | 0.747 | 0.355 | 0.152 |
| CCE traps, -200 V |          |          |          |       | 0.995 |       | 0.872 | 0.340 |

Bulk physics (N_D 3.38e13):

- **Emergent neutral-region removal: 7.4 /neq** (both acceptors filled). For comparison,
  the legacy scalar is 5, and Radmanovac et al. (arXiv:2510.11304) measured 5-7.5 at
  252.7 MeV.
- **This model's own compensation fluence** is therefore min N_D / (7.4 kappa):
  **2.35e12 p/cm^2 at 62 MeV** on the default stack, below the legacy 3.48e12.
  Metadata gives it as `trap_compensation_fluence` (and, since 2026-10-21, as
  `compensation_fluence` of the run), and `etna.api.simulation.trap_compensation_fluence`
  computes it; `legacy_compensation_fluence` is the scalar reference. The 3e12 column in the table above is already past it:
  CCE still matches legacy because the epi is over-depleted at -40 V.
- **N_eff below compensation does not match legacy.** Legacy N_eff falls with fluence,
  while the trap model's depleted N_eff rises. Only CCE agrees there.
- **Semi-insulating beyond compensation.** E_F reaches Ec-1.08 eV at 1e14 neq and
  Ec-1.20 eV at 1e16. The resistivity proxy is above 1e14 Ohm cm.
- **No type inversion.** E_F stays above midgap.

**Calibrated vs predicted.** No parameter was fitted to CCE or C-V data.

- The trap table is literature (Burin Table I).
- CCE below compensation is predicted, and it equals the legacy result (difference
  < 1e-5).
- CCE beyond compensation is a model prediction that depends on one stated assumption.
- Qualitative literature check, on different devices and biases:
  - **Gsponer-like device** (50 um, 1.5e14, 1e15 neq): CCE is 0.36 / 0.70 / 0.97 at
    -100 / -300 / -600 V. Measured: above 80 %.
  - **De Napoli-like device** (20 um, 1e14, at 1.3e14 / 7.3e14 cm^-2): losses are
    0.8 % / 11 % at -300 V and 2.8 % / 56 % at -100 V. Measured: 7 % and 22.5 %.

Same decade, bias-dependent, not a fit.

**Main assumption (read before using numbers beyond compensation).** In a depleted
region, steady-state occupancy empties both acceptors (emission times: Z1/2 about 40 ms,
EH4 about 30 min). EH6/7 stays about 90 % ionised. The depleted space charge is therefore
N_D + 1.43 Phi_neq, which is positive and grows with fluence. Two consequences:

- Full depletion needs more bias at high fluence. That bias need is what drives the CCE
  drop at -40 V.
- The model does not reproduce the falling 1/C^2 slope of measured C-V removal. In the
  depleted region it has the opposite sign.

EH6/7 (emission time about 1e14 s) never reaches steady state on lab time scales.

- If the occupancy is frozen at the neutral-bulk value, the compensated epi stays
  charge-neutral and acts as a dielectric with a uniform field. CCE would then be
  trapping-limited, which is higher.
- `trap_diagnostics` reports both cases. The engine uses steady state (the TCAD DC
  choice).

Small-signal C with traps is not implemented (it raises). At 0 to -1 V a semi-insulating
epi under the DC alpha proxy is photoconductive: a few points do not converge, and CCE is
not monotonic there. For C-V with traps, and for occupancies that depend on the
measurement history, see the next subsection.

#### Kinetic (time-dependent) trap occupancy and C-V (opt-in, `etna.core.trap_kinetics`)

This is a library module only: no `DeviceConfig` field, no facade, and every default is
unchanged. It integrates the SRH rate equation of every trap at every node,
df/dt = c_n n (1-f) - e_n f - c_p p f + e_p (1-f), along a bias / illumination
history. Poisson is solved with dark flat quasi-Fermi levels at each instant. Context and
sources: GitHub issue #5.

**Numerics**

- Backward Euler is folded into a numba Newton kernel, f_new = (f_old + dt a) / (1 + dt lam). It is
  implicit and unconditionally stable for the ~1e8 /s capture rates of a neutral epi, with
  an analytic Jacobian, parallel over the batch.
- The two limits are tested: dt = inf is the `radiation_traps` steady-state kernel (psi
  identical to 1e-14 V), and dt = 0 is a fixed-charge Poisson.
- Adaptive dt uses step doubling. The error is weighted by each trap's share of the
  local space charge (`norm="charge"`, the default) or is plain occupancy
  (`norm="occupancy"`).
- Bias ramps use a secant predictor and steps of at most `dV_max`.
- Convergence on a 20 um, 1e14 device, C-V from 0 to -40 V at 1 V/s, max |dC/C| against
  tol = 1e-5:

  | tol  | steps | max abs(dC/C) |
  | ---- | ----- | ------------- |
  | 1e-2 | 56    | 1e-4          |
  | 1e-3 | 150   | 2.3e-5        |
  | 1e-4 | 784   | 4.2e-6        |

- The analytic single-trap relaxation in a fully depleted epi is reproduced, and the
  error falls with tol (`tests/test_trap_kinetics.py`).

**API**

```python
from etna import DeviceConfig
from etna.core.trap_kinetics import KineticTrapSolver, hold, simulate_cv, neff_from_cv

cfg = DeviceConfig(epi_thickness_um=50.0, doping_profile="uniform", N_D=4e13,
                   radiation_traps=True)   # Burin table (the default "auto" = calibrated preset)
s = KineticTrapSolver.from_config(cfg, [0.0, 1e12, 2e12], energy_MeV=252.7,
                                  niel_convention="si-neq")  # p/cm2; numbers below: si-neq
cv = simulate_cv(s, V_points, freq_Hz=1e6, occupancy="kinetic", initial="equilibrium",
                 storage=(hold(600.0, 0.0),), rate_V_per_s=1.0)
cv.C_ac, cv.C_qs, neff_from_cv(V_points, cv.C_ac)   # F/cm2 (area 1), cm^-3
```

- **`occupancy`** takes one of three values:
  - `"kinetic"`: the history.
  - `"steady"`: the steady state at each bias (the `radiation_traps` choice).
  - `"frozen"`: the `initial` occupancy held fixed.
- **`initial`** (the occupancy left by the irradiation and storage) takes one of four
  forms:
  - `"equilibrium"`: Fermi-Dirac at 0 V.
  - `"in_beam"`: high injection, f = c_n / (c_n + c_p).
  - `"neutral"`: local neutral bulk.
  - An array.

  EH6/7 never relaxes, so this is a stated assumption: report results across policies.

- **`solver.run(state, [Segment(...), ...])`** runs general histories: ramps, holds,
  uniform illumination `G`.
  - Illumination uses a drift-limited photo-carrier estimate, capped at G * 1 us. It is
    crude.
  - Temperature is set per solver: build one solver per T.
- **C is computed with the trap charge frozen at AC.**
  - `C_qs` is the quasi-static small-signal C.
  - `C_ac` adds the dielectric-relaxation impedance of the undepleted region in series:
    1/Y = 1/(i w C_qs) + sum h / (sigma + i w eps).
  - On conducting devices `C_ac` = `C_qs` (tested to 0.5 %). On a semi-insulating epi
    C -> eps A / d.
- **Other observables:**
  - `generation_current` gives the steady-state SRH generation through the traps, and
    the net hole emission of the actual occupancy.
  - `flatband_resistance` gives the specific resistance of the compensated epi (the
    ohmic bound of the forward current).

**Results** (`lab/radiation_traps/bench/kinetics.py`, ~40 s; log in
`lab/radiation_traps/notes/KINETICS.md`). PREDICTED unless marked CALIBRATED.

- **Sign of the C-V removal is fixed.** Radmanovac-like device (50 um, 4e13,
  252.7 MeV, 1 MHz, N_eff at W 10-40 um), removal rate g per p/cm2:

  | Model                           | g (/p)     |
  | ------------------------------- | ---------- |
  | Steady state                    | 0.2        |
  | Kinetic                         | 2.7        |
  | Legacy scalar                   | 4.3        |
  | EH6/7 as acceptor (alternative) | 4.2        |
  | Measured (arXiv:2510.11304)     | 5.6 +- 0.8 |
  - Steady state is effectively no removal: N_eff rises, then falls only near the edge.
  - Kinetic is EH4-dominated. Z1/2 empties within ~40 ms in a depleted region.
  - The kinetic g is the same for the equilibrium start and for an in-beam start plus
    10 min or 1 day at 0 V.

- **The table's magnitude is short by ~2x, and scaling cannot fix it.**
  - CALIBRATED to that single g, the global rate scale is 2.06 +- 0.28. That predicts
    compensation at 3.1e12 p/cm2, against N0/g = 7.1e12.
  - Scaling EH4 alone gives 2.50 +- 0.39.
  - The g sensitivities per unit scale (3-fluence subset, g = 2.51 there) are
    EH4 2.06, Z1/2 0.69 and all 2.92: Z1/2 near the depletion edge still adds to g.
  - A two-number reading solves the C-V slope and the compensation fluence N0/g
    together. It assumes the measured linear N_eff(Phi) extrapolates to compensation, as
    it does for CNM W4. Result: s_EH4 = 2.86 +- 0.39 and s_Z1/2 = -0.06, which is about
    zero.
    - The sigma of s_Z1/2 from g alone vanishes, because both rows scale with g. Its
      real uncertainty is set by that assumption.
    - So the data are consistent with a slow acceptor at about 7 /neq carrying the
      removal, and Z1/2 carrying none at room temperature.
    - This is an assumption-laden reading, not a fit. Li et al. (arXiv:2503.09016)
      likewise report that Z1/2 does "not appear to significantly contribute".
  - The open question is in issue #5.
  - **Superseded** by the next subsection and `lab/radiation_traps/notes/MECHANISMS.md`. With
    onsemi's own flat-C-V bracket, the Z1/2 scale is 0.8 +- 0.65 (consistent with nominal),
    and the slow component (EH4 plus the calibrated permanent removal) is about 5.3 /neq.
    The "Z1/2 carries none" reading came from an extrapolated compensation fluence.
- **Temperature (PREDICTED, testable).** The kinetic g falls from 2.51 /p at 300 K to
  1.95 /p at 325 K and 0.77 /p at 350 K: EH4 emission speeds up with T. No measurement
  of this was found.
- **History.** After 1 h at -150 V, a down-sweep shows N_eff within 1 % of the
  up-sweep. Edge traps re-capture in microseconds, so no significant C-V hysteresis is
  predicted. No opened source reports any.
- **C-V flattening.** Gsponer-like device: 50 um, 1.5e14, 3x3 mm2, 10 kHz.
  - At >= 3e13 neq, C is flat to 0.9 %, at 15.46-15.59 pF. Measured: flat at 5e14 to
    1e16, about 15.5 pF.
  - The model flattening onset is about 1.6e13 neq for 6.2e13 doping. Kroll et al.
    measured flattening from about 3.1e13 neq (Si-kappa conversion).
- **Loss of rectification (ohmic bound).** At 1 kV forward the bound is
  I_fwd <= 1.5e-11 A at 1e15 neq and 1.4e-12 A at 1e16. Measured: < 1 nA, and equal to
  the reverse floor (< 10 pA). Qualitative only; double injection is not modelled.
- **Leakage.** Bulk trap generation rises with fluence but stays at 1e-25 to 1e-21 A,
  far below the ~1e-11 A floor of every measurement opened, so it cannot be tested. The
  measurements show reverse leakage unchanged (Radmanovac et al.) or falling with
  fluence (Kroll et al., Li et al.). This module does not address that: the dominant
  measured leakage mechanism is not in it.

**Calibration** (`etna.core.trap_calibration.calibrate_cv`). It fits only multiplicative
scales on the unmeasured introduction rates, globally (`"all"`) or per trap name, to
measured C-V curves.

- Inputs are `(fluence, CVData)` pairs from `etna.measurements.load_cv_csv`. Each
  file's frequency and sigma are used.
- It returns:
  - the scales with their Gauss-Newton 1-sigma, also scaled by the reduced chi^2;
  - the correlation matrix and the singular values;
  - identifiability flags and warnings for bounds, correlation > 0.95 and missing
    frequency.
- The device itself (geometry, pristine N_D) is input from pristine C-V or SIMS and is
  not fitted. The initial occupancy and the history are stated inputs.
- Tested by recovering a SYNTHETIC truth through the CSV loader: 1.7025 +- 0.0037 for
  1.7, within 1 sigma; the test bound is 3 sigma. The same model, mesh and protocol
  generate and fit those data, so this validates the plumbing and the covariance, not
  whether the model is adequate.
- No measured irradiated C-V curve ships with the repository, so no calibrated 62 MeV
  number exists yet. The INFN-LNS C-V vs fluence data
  (`.planning/CARRIER_REMOVAL_REVIEW.md` §6) would provide it.

**Not modelled**

- AC trap response (the lambda effect).
- Transport-limited quasi-Fermi splitting in a semi-insulating epi under bias.
- Poole-Frenkel emission.
- Negative-U two-level Z1/2 in the default table (available as an opt-in two-level
  trap, next subsection).
- Temperature changes within a run.

#### Missing C-V removal: two-level traps, direct removal, presets, calibration (opt-in)

Follow-up of issue #5. The kinetic Burin table gives about half of the measured
depleted-region removal (g 2.6 vs 5.6 +- 0.8 per p/cm2, Radmanovac et al.,
arXiv:2510.11304). Evidence, every number and its source:
`lab/radiation_traps/notes/MECHANISMS.md`; benches `bench/mechanisms.py` (about
4 min) and `bench/cce_kinetic.py` (about 50 s). Since 2026-10-21 the
`radmanovac_calibrated` preset below is the DEFAULT damage model
(`DeviceConfig().radiation_traps == "auto"`); `radiation_traps=True` is still the steady
Burin table. The numbers of this subsection were computed in the `si-neq` convention
(d = 2.87 /neq); in the default `sic-damage` convention the re-fitted d is 4.69 /neq and
the per-proton slopes are the same within 0.1 /p (see "NIEL convention and the
calibrated removal").

**New model pieces** (fields of `RadiationTrapModel` / `DeepTrap`):

- **Two-level traps.** `DeepTrap(..., E2_below_Ec=, sigma_n2=, sigma_p2=)` adds the
  level of the second electron (negative U: E2 deeper). Example: Z1/2 = V_C(=/-/0).
  - Occupancy is carried per level row as g1 = P(>= 1 e) and g2 = P(2 e).
  - Steady state: the zero-flux chain. Kinetics: a 2x2 backward-Euler system in the
    numba kernel, with an analytic Jacobian.
  - Tests: the steady state equals the analytic chain and Fermi-Dirac; the relaxation
    converges to the matrix exponential.
  - The steady Poisson kernel refuses two-level traps. The CCE path routes them through
    `trap_kinetics`, whose steady limit is identical.
- **`donor_removal`** (cm^-1 per neq): direct, occupancy-independent donor loss on the
  epi.
- **`occupancy="kinetic"`**: CCE uses the kinetic dark state. That state is the
  `initial` state at 0 V, a ramp at `ramp_V_per_s`, then a hold of `hold_s`; its trap
  charge is then frozen. The default is `"steady"`.
- **`lifetimes="occupancy"`**: CCE lifetimes count only levels able to capture. A
  filled level cannot capture an electron. The default is `"legacy"`, which uses
  `compute_K_tau`.

**Presets** (`etna.core.trap_presets.trap_preset(name)`; pass the returned model as
`DeviceConfig(radiation_traps=...)`):

| Preset | Content | Label |
| --- | --- | --- |
| `burin_steady` | the default trap table | literature |
| `burin_kinetic` | same table, kinetic occupancy | literature |
| `literature_kinetic` | Z1/Z2 two-level (Capan et al. JAP 124 (2018) 245701, Table I: energies and Arrhenius emission, [Z2]:[Z1] = 4.4, 2.5 V_C/neq assumed) + EH4 + EH6/7 (sigma_n 1e-13) + kinetic + occupancy lifetimes | PREDICTED |
| `radmanovac_calibrated` (**the default damage model since 2026-10-21**, `"auto"`) | `literature_kinetic` + direct removal d per convention: 4.69 +- 0.30 /neq (sic-damage, default), 2.87 +- 0.22 (si-neq), 5.10 +- 0.32 (sic-1mev) | CALIBRATED to Radmanovac Table 2 (per proton; re-fitted per convention, see "NIEL convention and the calibrated removal") |
| `radmanovac_calibrated_deep_acceptor` | same rate on a midgap acceptor (Ec-1.9) | CALIBRATED; identical C-V |

**Calibration on published numbers** (`trap_calibration.calibrate_removal`):

- `RemovalObservation` takes either a removal slope g read at fixed depletion width
  (`kind="slope"`) or the fluence at which the C-V turns flat (`kind="flat_fluence"`,
  `semi_insulating_fluence`).
- The fit is linear in the parameters: scales of trap rates, or the absolute
  `donor_removal`. It uses bounded least squares with re-linearisation.
- It reports sigma, correlation, singular values and the weakest parameter
  direction.

**Results** (PREDICTED unless CALIBRATED):

- **Negative-U Z1/2 does not close the gap.** Z1/Z2 empty within about 0.1 s in
  depletion, so the slope does not see them. g is 2.41 /p with 2.5 V_C/neq and 3.13 /p
  with 5.0 V_C/neq.
- **EH4 as a donor collapses g** to 0.35 /p. This follows the CAV(+/0) assignment of
  Karsthof et al. (arXiv:2007.03985).
- **Hazdra's EH4 pair** changes g by less than 0.2.
- **EH6/7 as an acceptor** gives g = 3.8 /p.
- **Identifiability.**
  - On the 3 Radmanovac slopes, a slow EH4-like acceptor and a permanent removal
    (direct, or a midgap acceptor) are degenerate, with correlation -0.999.
  - Direct removal and a midgap acceptor give identical C-V below about 500 K.
  - Fast acceptors (Z1/2) are seen only through the flat-C-V fluence: scale 0.8 +- 0.65.
- **Calibrated fits.**
  - Direct removal: d = 2.87 +- 0.22 /neq, chi2 = 1.13 for 2 dof. Model g:
    onsemi 4.94, W2 4.81, W4 4.40.
  - EH4 scale instead: 2.28 +- 0.10, chi2 = 2.06. The chi2 difference does not
    discriminate (the sigmas are width spreads). The preset uses direct removal
    because the EH4-only fit needs EH4 at 2.3x its DLTS rate, and the removal
    literature points to a permanent compensating species.
- **Predicted flat-C-V fluences**, all inside the measured brackets:
  - onsemi: 4.95e12 p/cm2;
  - W4: above 3.5e13;
  - Kroll et al.: 1.2e13 neq;
  - Gsponer et al.: below 5e14 neq.
- **Discriminating prediction: temperature.** g(onsemi) at 300 / 325 / 350 K:
  - direct removal: 4.94 / 4.56 / 3.22;
  - EH4-scale alternative: 5.19 / 4.32 / 1.32.
- **EH6/7 sigma_n.** 9e-12 (Gaggl/Burin) emits 40-160x faster than measured EH6/7 DLTS
  peaks in this code's prefactor, and its cited source has no EH6/7 value. 1.4e-14
  (audit C-4) has no primary source and emits 4-15x too slowly. About 1e-13 is
  emission-consistent at 1.60 eV.
  - C-V removal is insensitive to sigma_n.
  - The default stays 1.4e-14 (MECHANISMS.md section 4).
- **CCE beyond compensation** (default stack, 62 MeV), at -40 V, for 1e15 / 5e15
  p/cm2:

  | Model                                     | CCE at 1e15 | CCE at 5e15 |
  | ----------------------------------------- | ----------- | ----------- |
  | steady traps                              | 0.355       | 0.152       |
  | kinetic, legacy lifetimes                 | 0.831       | 0.472       |
  | kinetic, occupancy lifetimes (calibrated) | 0.919       | 0.686       |
  - The steady-to-kinetic step dominates.
  - The calibrated removal itself changes CCE by less than 0.01 up to 1e15.
  - The unknown post-irradiation occupancy of EH6/7 dominates beyond about 3e14. With
    an in-beam start and a 1 h hold, CCE is 0.477 at 1e15. Report CCE as a range.

**Not modelled.** AC trap response, Poole-Frenkel emission, temperature changes within
a run, degeneracy factors (all 1), and trap-assisted leakage beyond SRH. The removal
species' identity is still open: the C-V data cannot fix it, but a C-V at 350 K or
EH4 DLTS at the same beam could.

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
and frozen-field fixed point as devsim, parity <= ~1e-11 relative, about 8 ms for 20
points), otherwise on devsim (about 0.26 s for the same 20 points). `solver="devsim"`
forces the devsim reference, `solver="engine"` raises outside the box. Quote the result as a budget estimate
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
engine inside its scope (section 7; about 5 ms, parity <= ~7.5e-12 with Auger on),
otherwise on devsim (about 0.5 s with `solver="devsim"`); `run_transient` runs on devsim
(about 1 s). Any variation with dose rate is a **numerical sensitivity bound**. No plasma or high-injection physics is modelled. Both
solvers of the steady-state sweep solve **SRH + Auger** (C_n 5e-31, C_p 2e-31 cm^6/s,
Galeckas et al. 1997, 300 K constants; override with
`etna.core.flash_recombination.cce_vs_dose_rate(..., C_n=, C_p=)`). Auger lowers CCE by
a dose-rate-independent offset: 3.4e-7 on the default stack at -30 V, up to 5.6e-5 on
notebook 05b's grid (5 um epi, N_D_bulk 5e14, -10 V), from minority holes recombining in
the 1e19 n+ substrate; it scales ~linearly with C_n, so the literature spread (up to
~2e-30) allows up to ~2e-4 there. It adds no dose-rate dependence (<= 1.5e-8). Until
2026-10-20 the Auger term was silently dropped (both solvers SRH-only); `run_transient`
dropped it until 2026-10-21 (now SRH + Auger too). Its direct-anode estimator does not
resolve the offset (round-off ~1e-4 of q*integral(G)): on the default stack Auger moves
the transient CCE by ~1e-12, at the 05b corner by 4e-8..9e-5 depending on the dose rate;
use the steady-state numbers above for the Auger magnitude.

```python
from etna import DeviceConfig, run_flash_recombination, run_transient

f = run_flash_recombination(DeviceConfig())
t = run_transient(DeviceConfig())
print(f.y.round(6), t.metadata["transient_cce_raw"].round(6))
```

#### Transient terminal current (conduction + displacement)

The transient drivers (`TransientSolver.simulate_pulse`, `simulate_pulse_train`,
`single_particle.simulate_single_particle`) record three series per BDF1 step:

| key                     | quantity                                                                                                              |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `currents`              | conduction current, `extract_contact_current` (continuity equations only)                                             |
| `currents_displacement` | displacement current, `displacement_current(Q_k, Q_k-1, dt_k)` of the Poisson contact charge `extract_contact_charge` |
| `currents_total`        | terminal current = conduction + displacement: the physical I(t)                                                       |

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

    | quantity            | displacement share or effect                                   |
    | ------------------- | -------------------------------------------------------------- |
    | 1 ps injection step | 5 % of the current                                             |
    | later steps         | <= 1.6 % of the peak                                           |
    | I_peak              | 5.49e-5 total vs 5.42e-5 conduction A/cm                       |
    | t_collection        | 0.11679 ns (new definition) vs 0.11739 ns (old step-end value) |
    | CCE                 | 0.999970 either way                                            |

- **Old definition and the stale table column.** Before 2026-10-20, `t_collection` was the
  end time of the BDF1 step that crossed 95 % of the conduction charge (`searchsorted` on
  a cumulative sum that assumes it is monotonic). That made it quantised to the 10 %
  step growth. The tracked `data/cce_let_table_100um.json` / `_300um.json` were rebuilt on
  2026-10-21 with the c-axis mobility default (`lab/mobility_caxis/bench/
rebuild_let_tables.py`, notebook-16 keywords, about 12 min each with 6 workers, after the
  merge of the contact-cell meshes), so their `t_collection_ns` column is the current
  definition: 0.1168 ns at every LET (the holes set it; the earlier tracked tables had the
  stale 0.208 ns). CCE 0.9999662-0.9999742 (100 um) and 0.9999569-0.9999741 (300 um); a
  `legacy_perp` rebuild (before that merge) gives 0.9999625-0.9999703 (100 um), the
  previously tracked table 0.9999682-0.9999723. The JSONs carry `"mobility": "c_axis"`.
  `load_cce_let_table` interpolates CCE only.

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
f = 0.9948 at 41 µm, 0.9726 at 46 µm and 0.488 at 51.5 µm (middle of the gap; graded depth
mesh, NIST Bragg alpha; `lab/edge_cce/notes/LOG.md`). The f ratios move by <= 4e-4 with
the depth mesh; the absolute CCE moved by -5.7e-3 (2026-10-21, `lab/mesh_contact`).
`run_microdosimetry` does not apply it (it has no CCE step).

#### Energy-deposition sampler, event-by-event pulse heights, Allpix Squared exchange (core API, issue #7)

These are new opt-in modules. Nothing in the `run_*` facades calls them, and no defaults
changed. Before them ETNA had no transport of its own: an event's deposit was whatever an
external Geant4 file said, or the deterministic LET x chord.

**`etna.core.mc_transport.sample_track_events(particle, energy_MeV, sv_um=(100, 100, 10), ...)`**

A numba float64 parallel condensed-history sampler for `"proton"`, `"alpha"` and `"c12"` in
`"sic"`, `"si"` or `"water"`, through a box SV. What goes into it:

- **Irradiation.** `irradiation="normal"`, an angle in degrees (parallel beam), or
  `"isotropic"` (mu-randomness; the mean chord is Cauchy's 4V/S).
- **Mean loss.** NIST PSTAR/ASTAR electronic stopping. SiC is Bragg-composed. C-12 is a
  MODEL: z_eff² x PSTAR at the same velocity, with z_eff from Ziegler's formula as in
  Geant4.
- **Fluctuations.** A line-by-line port of the Geant4 11.3 models: `"urban"` for protons
  and `"ion"` for alpha/C-12 (G4IonFluctuations' Gaussian branch, without Yang's factor;
  `"universal"` at high velocity), as in G4EmStandardPhysics_option4. `"bohr"` and
  `"none"` are also available.
- **Delta rays.** Explicit above `tcut_keV` (default: electron range = 5 % of the
  smallest SV side, 6.6 keV for 10 um), with a restricted loss below. They travel
  straight with the NIST ESTAR CSDA range (`data/nist_estar/`) times
  `delta_range_factor`, so they can leave the SV. `surround_um` adds SiC around the SV:
  delta rays made there can enter, and the beam is degraded by it.
- **Stepping.** Adaptive (`step_fraction`, the energy-loss fraction per step), so an
  alpha stopping in the SV ends at the CSDA range.
- **Random numbers.** A SplitMix64 stream per (seed, event), started at a hash of the
  pair, so events are independent (issue #17: before, event k+1 reused event k's draws
  shifted by one and the standard errors were 1.6-2x too small). Results are
  independent of the thread count, and `event_offset` continues a stream.

It returns `TrackEvents`: `edep_keV` (= primary + delta part), `chord_um`,
`E_sv_entry_MeV`, entry point and direction, `stop_path_um` and optional `hits`
(`max_hits>0`; recorded primary deposits are spread over points `hit_step_um` = 1 um
apart, recording only). The deterministic baseline is
`fluctuations="none", delta_rays=False`, which gives exactly the CSDA mean along the chord
(`mean_csda_deposit_keV`).

**`etna.core.mc_response`**

- `lineal_moments(edep, l_bar)` gives event-level y_F and y_D with delta-method standard
  errors.
- `adaptive_lineal_energy(particle, E, rel_tol=0.01, ...)` samples batches until both
  relative standard errors are below `rel_tol`.
- `default_l_bar_um(sv_um, irradiation)` (the default l_bar of `adaptive_lineal_energy`)
  is the sampler's mean chord: t for a normal beam, V / A_proj with
  A_proj = wx wy cos(theta) + wy t sin(theta) for a beam tilted in x-z, 4V/S for
  isotropic. It is not the slab t / cos(theta) (60 deg, 100x100x10 um: 17.05 um, not
  20; issue #18).
- `pulse_heights(edep, chord, cce_interp, fano=None, enc_electrons=0, threshold_keV=0)`
  runs, per event: N = epsilon/W, then optional Gaussian Fano smearing, then
  CCE(LET = epsilon/chord) from `load_cce_let_table`, then ENC noise, then the
  threshold.
  - W defaults to ETNA's 7.8 eV.
  - Allpix's 4H-SiC default is 7.6 eV with F = 0.1. ETNA sets no default F.

**`etna.core.mc_exchange`** (file formats only; no Allpix, Geant4 or ROOT install)

- `write_deposition_csv(path, events)` writes and `read_deposition_csv(path, unit_length="mm",
unit_energy="MeV")` reads the Allpix `DepositionReader` CSV (`Event: N` headers from
  0, `pdg, t, E, x, y, z, detector, track, parent`).
  - The read DataFrame goes straight into `process_mc_ensemble`.
  - `tracks_to_dataframe(events)` gives one row per event: position at the chord
    midpoint, plus `ekin_MeV` (SV-entry energy) for the kappa(E) lookup.
- `write_init` / `read_init` handle the Allpix INIT field format, laid out as
  `FieldParser::parse_init_file` reads it: 1-based i j k, z fastest, header sizes
  (thickness, x, y) in um, values on cells, z from the back to the implant side.
- `export_allpix_fields(run_field(...), outdir)` writes three files and returns the
  matching `ElectricFieldReader` / `DopingProfileReader` / `WeightingPotentialReader`
  config snippet:
  - `etna_efield.init`: V/cm, E = -grad psi of the solved potential, so it points to +z
    (the p+ side) in the depleted epi.
  - `etna_doping.init`: 1/cm³, n positive.
  - `etna_weighting.init`: planar pad z/t.

  1D results give a laterally uniform map. 2D half-cell results are mirrored to a full
  plane (`PIXEL_FULL`).

  The snippet sets `field_depth` / `doping_depth` / `potential_depth` to the exported
  thickness (the whole stack, 11.3 um for the default device; `depth_range_um` picks a
  sub-range). Allpix otherwise maps each file onto the whole sensor thickness and
  silently stretches it (issue #19: 47 V instead of 53.2 V across a 10 um sensor at
  -50 V). The Allpix sensor must be at least that thick (Allpix rejects deeper maps);
  pass `sensor_thickness_um=` to get a warning on any mismatch. The maps sit at the top
  (implant side) of the sensor, while `write_deposition_csv` hits are SV-centred.

- These files are checked against the quoted parser layout and against each other, not in
  an Allpix run.
  - Allpix truncates E/W per deposit to an integer in `DepositionReader`, so avoid many
    sub-W deposits.
  - Allpix has no SiC mobility model (use `constant`/`custom`).

```python
from etna.core import mc_transport as mt, mc_response as mr

ev = mt.sample_track_events("proton", 62.0, sv_um=(100, 100, 10), n_events=200_000)
m = mr.lineal_moments(ev.edep_keV, 10.0)
print(round(m["y_F"], 2), round(m["y_D"], 2))  # 2.53 2.85 (bare SV; deterministic LET x t: 2.77 2.77)
```

Benchmarks: `uv run python lab/mc_transport/bench/bench.py` (25-100 s, depending on machine load; numbers in
`lab/mc_transport/notes/LOG.md`).

- **Consistency (not validation; same tables).** The mean deposit equals the NIST CSDA
  integral within 0.1 % (8e6 events).
- **Variance.** 0.83-1.04 x the Vavilov/Bohr variance xi T_max (1 - beta²/2).
- **Vavilov regime (kappa 0.1-0.5).** The most probable value is within 0.2-6 % of the
  numerical Vavilov density. The FWHM is 10-12 % narrower.
- **Landau regime (kappa < 0.03).** The MPV is 2-12 % above PDG Eq. 34.12 and the FWHM is
  1.2-1.8 x 4 xi. The PDG and Bichsel 1988 report that thin-Si distributions are wider
  than Landau's. There is no Geant4 run yet, so how close the port is to Geant4 itself is
  unchecked.
- **Alpha 5.486 MeV in SiC.** Range 18.41 um vs CSDA 18.43 um. Straggling 0.156 um vs
  the Bohr integral 0.165 um.
- **Delta rays.** A bare 10 um SV keeps 0.914 of the mean (62 MeV p). A 100-300 um SiC
  surround restores 1.01.
- **Effects separated (62 MeV p, normal beam).**
  - Straggling alone keeps y_F (2.77) and raises y_D from 2.77 to 3.61.
  - Delta-ray escape from a bare SV lowers y_F to 2.53 and y_D to 2.85.
  - A 300 um SiC surround brings them back to 2.76 / 3.22.
  - The bare-SV default is therefore not a detector prediction for an SV inside a chip:
    set `surround_um` to the real overlayer/substrate.
- **Speed.** 10 M events/s on 10 threads and 2 M events/s on 1 (Apple M4, light load),
  falling to 2-4 M events/s under load average 45-110. A GPU port was not attempted: one
  spectrum takes under 1 s on the CPU.

Not modelled:

- nuclear interactions: about 1e-4 per 10 um traversal for C-12, from the computed
  Tripathi cross sections;
- multiple scattering of the primary;
- electron scattering (straight delta rays give the maximal escape distance);
- track structure below T_cut = 1 keV, which makes it invalid for sub-um SVs.

The Geant4 cross-check is planned (issue #7).

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
between sections that no devsim device leaked. On the default device every section runs on
the engines and the whole run takes about 0.06 s with a warm numba cache (legacy preset
about 0.9 s: its 0..-100 V dark sweep runs on devsim; measured 2026-10-21). A section's `status` is `ok`, `warning` (it ran
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
- **devsim only**: `run_transient`, `dark_current_sweep(device_info)` (a
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
  `notebooks/figures/`) or set the variable. Both `notebooks/figures/` and the root
  `figures/` are gitignored; the only tracked figure is
  `figures/cv_validation_bruzzi2023.png` (section 13).
- `ETNA_DEVSIM_STOCK_SOLVER=1` keeps devsim's stock UMFPACK callback. By default
  `etna.core.devsim_fastsolve` installs the same callback with its per-solve library
  reload cached (bit-identical results, `tests/test_devsim_fastsolve.py`).
- devsim bias points start from the engine solution (`etna.core.drift_diffusion.seeded_bias`):
  the engine's dark state at the target bias plus one strict devsim solve replaces the
  0.5-2 V ramp from 0 V (with Auger or the TAT/Hurkx model the engine state is only the
  Newton start). It is never used with incomplete ionization, illumination, forward bias
  or laterally non-uniform 2D devices; there, or if that solve fails, the old ramp runs
  (`tests/test_seeded_bias.py`).
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
E)` the epi is fully compensated. The default (calibrated trap) model solves it as a
semi-insulating epi (a prediction); the legacy scalar model (`radiation_traps=False`)
returns NaN with a reason.

**Microdosimeter alpha test not reproduced.** With the preset (N_D 1e14), the model CCE is
flat to 0.5 % from 0 V. The data drop to 0.908 at 0 V and rise 1 % from -20 to -50 V. No
drift-diffusion variant fits (chi2 >= 228). A field-dependent columnar-recombination loss
does fit, but it is lab-only (`lab/microdosimeter/notes/LOG.md` §9).

**Incomplete ionization** is the legacy scalar by default; the carrier-dependent model
is opt-in (§5.3) and not yet calibrated against the C-V.

**Validation variants as figures** (`etna.viz.validation`, run 2026-10-20, pinned in
`tests/test_viz_figures.py`). Both use existing solvers only; the reference points are in
`data/validation/` with their source and uncertainty in the file header. "Micro Fig. 4b"
(alpha CCE, microdosimeter paper) and "Photons Fig. 4b" (C-V, dosimeter paper) are different
figures.

- `tilted_alpha_shape()`: `run_cce_tilted_alpha(microdosimeter_config())` at its defaults,
  normalised at -50 V, against the digitised Micro Fig. 4b points (read-off estimates,
  ±0.0005). Shape only, nothing tuned. Result: model 0.9928 at 0 V vs data 0.9083 (the
  same with the perpendicular-to-c mobility used before 2026-10-21); rise -20 -> -50 V
  +0.0001 (model) vs +0.0096 (data); rms 0.0354. Not reproduced, as stated above.
- `w_creep(config=DeviceConfig())`: `run_cv` (batched) W(V), depletion edge and eps*A/C,
  against the Photons Fig. 4b W labels 1.68 / 9.52 / 9.73 µm (±0.005 rounding). Result on
  the default stack (graded epi + 1.14 µm n-/n+ transition): eps*A/C 1.655 / 9.519 / 9.750
  µm; creep -10 -> -30 V +0.232 µm (eps*A/C) and +0.256 µm (edge) vs +0.21 µm in the labels.
  The default doping is a fit to this C-V curve, so this is a regression of the fit, not an
  independent validation. (The lab numbers in `lab/microdosimeter/notes/LOG.md` §4 were for
  the older default without the transition.)

**Not modelled:**

- high injection or plasma (FLASH). A single dense ion track (plasma / columnar
  recombination) is handled only by the opt-in track engine (§19), at normal incidence.
- edge termination and guard rings on the default path
- oxide charge
- impact ionisation and breakdown on the default facades. They are available opt-in, 1D
  only, through `run_breakdown` / `run_avalanche` ([§18](#18-physical-model-library-opt-in-high-field-mobility-anisotropy-avalanche-and-breakdown)).
  The breakdown voltage there is the parallel-plane value, with no edge termination.
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

| page                    | purpose                                                                         | inputs                                                                                                                                                                                                                                                                                                                                 | outputs                                                                                                                                                                                                                                                                                                                                                                                                          |
| ----------------------- | ------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Home                    | orientation, current config summary                                             | none                                                                                                                                                                                                                                                                                                                                   |                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Simulate → Figures      | one workflow: device -> observables -> figures, data, animation (`etna.viz`)    | device (sidebar or a preset: default, microdosimeter, tabulated, previous calibrated, legacy); observables (C-V, field, CCE, dark current, FLASH, temperature, transient (slow, devsim), tilted alpha, CCE(LET) table, two validation variants; 2D: maps incl. carriers); sweep V_stop / points, single bias; figure format, DPI, size | engine-box badge, progress bar, session cache, per-observable Plotly + solver / validity / honesty notes, publication figure preview, PNG/SVG/PDF per figure, zip of all figures + CSVs, animation (manim or matplotlib, mp4/gif) with st.video                                                                                                                                                                  |
| Run All                 | every section for the sidebar device (`etna.run_all`)                           | section multiselect (the 2D-engine sweep is opt-in)                                                                                                                                                                                                                                                                                    | summary table (status, solver, key figures, reason, wall time), a curve per section, JSON/CSV                                                                                                                                                                                                                                                                                                                    |
| C-V Analysis            | `run_cv`                                                                        | max reverse bias (default -50 V), points (40)                                                                                                                                                                                                                                                                                          | C-V, Mott-Schottky, eps*A/C and depletion edge vs V, full-depletion marker, CSV                                                                                                                                                                                                                                                                                                                                  |
| Charge Collection (CCE) | `run_cce`                                                                       | V_start (-10), V_stop (-40), points (30); reverse bias only                                                                                                                                                                                                                                                                            | CCE vs bias, CSV                                                                                                                                                                                                                                                                                                                                                                                                 |
| Field Map               | `run_field`, 1D or 2D                                                           | reverse bias (default -50 V, down to -200)                                                                                                                                                                                                                                                                                             | 1D: E(x), potential, net doping, CSV. 2D: heatmap viewer (field, potential, doping)                                                                                                                                                                                                                                                                                                                              |
| Radiation Damage        | `run_radiation_damage` | damage model selector (default: calibrated kinetic traps = `radiation_traps="auto"`; legacy carrier removal; Burin Table I steady traps), particle (proton/neutron), NIEL convention (default sic-damage), energy (1-1000 MeV protons), fluence range (traps: 0.01 x the model's compensation fluence to 1e16 neq; legacy: 0.7e-3 to 0.7 x), points (6), bias (-20 V) | CCE vs fluence with NaN gaps explained (legacy). A persistent banner says kappa is SR-NIEL in the sic-damage convention (r_n = 1.069) and what is calibrated; an info box gives the compensation fluence of the chosen model and the legacy reference |
| Dark Current            | total leakage vs T at one bias (`ParametricSweep` over T of `run_dark_current`) | T range (250-400 K), points (6), bias (-20 V), advanced N_t                                                                                                                                                                                                                                                                            | I_total vs T, CSV. A calibration banner is shown; only the total is shown                                                                                                                                                                                                                                                                                                                                        |
| Microdosimetry          | `run_microdosimetry`                                                            | MC events CSV upload, SV thickness (10 µm) and width (100 µm, depth = width)                                                                                                                                                                                                                                                           | y·d(y), y_F / y_D, CSV. Independent of the sidebar                                                                                                                                                                                                                                                                                                                                                               |
| Batch Sweep             | `ParametricSweep`                                                               | one of epi thickness, p+ thickness, N_A, N_D junction, N_D bulk, transition length, T, area; comma-separated values; facade (CCE vs bias, C-V, CCE vs temperature)                                                                                                                                                                     | overlaid curves, bulk CSV, per-value failures listed                                                                                                                                                                                                                                                                                                                                                             |
| Export                  | engineering files (section 10)                                                  | spec sheet (reuses the Run All report when it matches); profile bias (default -30 V); animation button                                                                                                                                                                                                                                 | layout parameter file in use (upload an alternative TOML) with its PLACEHOLDER list; spec sheet JSON/CSV/Markdown/PDF; cross-section SVG/PNG; mesh/profile CSV, VTK and VTU (ParaView/VisIt/meshio); two 480p manim MP4s (needs `viz`); reference layout SVG and GDSII (needs `gdstk`) with the disclaimer "reference layout, not TCAD-validated; not a fabrication layout while any dimension is a PLACEHOLDER" |
| Measurements            | owner-supplied C-V, SIMS/SRP, I(T) files vs the model (section 17)              | CSV upload per tab (header-only templates downloadable); C-V smoothing lambda; SIMS junction depth, knot spacing, resampling; I(T) prefactor exponent                                                                                                                                                                                  | validation messages (file line + column); measured-vs-model figure with residuals; donor-table JSON (C-V inversion or SIMS resampling), comparison / E_a JSON, Markdown report, PNG                                                                                                                                                                                                                              |

Page details not in the table:

- Simulate → Figures caches every result in the session under (config hash, observable,
  inputs): pressing _Run_ again with the same device and inputs recomputes nothing, and a
  changed V_stop only reruns the sweeps. Figure bytes are built when a download is clicked
  (and cached); the preview toggle renders a 150 dpi PNG. Results for another device are
  flagged as stale. A 2D device the engine does not serve falls back to `run_field`, without
  carrier maps (the page says why).

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
21 notebook workflows had no page. Since then Run All has added the FLASH steady-state
section (notebook 04) and Simulate → Figures the FLASH transient (part of notebook 08) and
the committed CCE(LET) table (the output of notebook 16). Notebooks 12, 14, 19 and 20
remain notebook-only.

---

## 10. Exports (Python)

`etna.export` holds pure builders: they take a config plus finished results and never
simulate. Each one carries `etna.export.VALIDATION_LABEL`: "simulation output; fast
engines cross-checked against the devsim reference solver; not validated against
measurement except where stated". The devsim comparison is a solver-to-solver
cross-check of the same model, not a validation of the model.

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
  variable `ETNA_LAYOUT_PARAMS`, or the upload box on the Export page. In an installed
  wheel the default (`ex.LAYOUT_PARAMS_PATH`, under `etna/_data/foundry/`) is package
  data: copy it out, edit the copy and set `ETNA_LAYOUT_PARAMS` (section 3).
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

### Publication figures and animations (`etna.viz`)

`figure_for(obj, kind=None, size="single")` turns a finished result into a matplotlib figure
in one house style (Okabe-Ito colours, inward ticks, 8 pt text, units on every axis: pF,
µm, kV/cm, A, Gy/s) and `export_figure(fig, "png" | "svg" | "pdf", dpi=300, size=None)`
returns the bytes (`save_figure(fig, path)` writes them). Sizes: `"single"` 89 x 70 mm,
`"double"` 183 x 105 mm, `"slide"` 254 x 143 mm, or `(w_mm, h_mm)`. The canvas is exact (no
tight bbox): a PNG is int(w·dpi) x int(h·dpi) pixels; SVG/PDF keep the text as text. The
style is applied with `rc_context`, so global matplotlib settings are untouched.

```python
import etna
from etna.viz import figure_for, export_figure, kinds_for, save_figure

cfg = etna.DeviceConfig()
cv = etna.run_cv(cfg, 0.0, -60.0, 61)
print(kinds_for(cv))                        # ['cv', 'depletion', 'mott_schottky', 'stack']
save_figure(figure_for(cv, "depletion"), "W_vs_V.pdf")
png = export_figure(figure_for(etna.run_field(cfg, -30.0)), "png", dpi=600, size="double")
```

Kinds (`etna.viz.KINDS`): `cv`, `mott_schottky`, `depletion`, `field` (E and psi vs depth
over the layer bands), `cce` (also the tilted-alpha result), `cce_let` (from
`etna.viz.load_let_table("data/cce_let_table_100um.json")`), `dark_current`, `flash`,
`transient`, `temperature`, `damage`, `microdosimetry`, `map_potential` / `map_field` /
`map_doping` / `map_electrons` / `map_holes` (2D, filled contours plus labelled contour
lines), `stack` (any result or a `DeviceConfig`), `tilted_alpha_shape`, `w_creep`.

Every figure carries a footer: etna version, solver, the 12-hex config hash
(`etna.viz.config_hash(cfg)`, sha256 of the sorted config), UTC date, the validation label,
and the notes that apply. No solver returns an uncertainty and the engine/devsim parity is
far too small to draw, so curves have **validity bands, not error bands**: |V| beyond
`cross_checked_to_V` hatched, the reach-through voltage dotted, outside-box / truncation /
fallback notes, dark-current points past reach-through open, unsolved damage points marked
with their reason, and the honesty notes (dark current = budget estimate, FLASH = sensitivity
bound, kappa Si-neq with the SiC normalisation open, CCE(LET) table ~6e-6 off fresh builds). Error bars appear only on
digitised data (section 8).

2D carrier maps: `etna.viz.maps.carrier_state_2d(cfg, bias_V)` runs the 2D engine and also
returns n and p (`run_field` does not). Against devsim Electrons / Holes on the same device
the relative difference is <= 1.5e-13 wherever the density is >= 1e2 cm^-3 (half-width 50 µm
at -30 V and 25 µm at -10 V); below that both are at their numerical floor and the map clips.
Configurations outside the engine (or the legacy model) raise `CarrierMapUnavailable`.

Animations: `data = etna.viz.animation_data(cfg)` (real run_cv W and run_field E(x) at 0 to
-60 V), then `etna.viz.animate(data, scene, fmt="mp4" | "gif", engine="auto")` for `scene in
etna.viz.SCENES` (`StackScene`, `DepletionScene`). `auto` uses manim (`render_animation`,
needs `viz`) for mp4 when it imports, otherwise `render_frames`, the matplotlib fallback:
one frame per solved bias, no interpolation; mp4 through ffmpeg (H.264, yuv420p), gif through
Pillow (no ffmpeg needed); 2.4-2.8 s per scene measured on a heavily loaded M4.

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

### 11.1 General 2D engine for laterally non-uniform and axisymmetric structures (`etna.core.fv2d`, opt-in)

`etna.core.engine2d` handles only laterally uniform diodes. `etna.core.fv2d` solves any 2D
tensor-grid device:

- per-node doping and mobility;
- any number of contacts on any nodes (top segments, side walls, whole electrode columns);
- Cartesian or axisymmetric (r, z) geometry.

The physics and discretisation are devsim's: box method, Scharfetter-Gummel currents,
SRH, ohmic contacts and static doping. Incomplete ionization is not supported in fv2d.

How it solves:

- The code is numba float64. Newton uses SuperLU; the Gummel chord solves use a banded LU
  or SuperLU, with the shorter grid dimension numbered fastest.
- The dark state at reverse bias comes from an adaptive Newton ramp on the flat
  quasi-Fermi Poisson. `dark_check()` runs a G = 0 Gummel pass to verify it is the DD
  steady state: 2 iterations, |du| < 1e-14 V on the notebook-19 structures.
- The illuminated solves are Gummel with chord LUs, run in parallel threads.

Nothing in the default paths uses fv2d.

```python
from etna.core.alternative_structures import create_guard_ring_device
from etna.core.fv2d import (mesh_from_builder, prepared_solver, cce_lateral_scan_fv2d,
                            contact_sharing_fv2d, cut_node_sets, alpha_stripe_generation_fv2d)

mesh, _ = mesh_from_builder(create_guard_ring_device)        # devsim's mesh, verbatim
s = prepared_solver(mesh, {"cathode": 10.0})                  # anode and guard at 0 V
info = mesh.info
ns = cut_node_sets(mesh, info["sv_half_width_cm"] + 0.5 * info["gap_cm"])
scan = cce_lateral_scan_fv2d(s, n_points=11, contact="anode", node_sets=ns,
                             x_max_cm=info["sv_half_width_cm"] + 1.5e-4)
frac, residual = contact_sharing_fv2d(s, alpha_stripe_generation_fv2d(mesh, 50e-4), ns)
```

The scan result has the keys of `charge_collection_2d.cce_lateral_scan`, so
`mc_coupling.process_mc_ensemble(cce_lateral=scan)` accepts it unchanged. It also carries
`residuals` (current conservation) and `direct_mismatch` (balanced vs direct
extraction).

`cce_lateral_scan(device_info, ..., engine="fv2d")` re-solves a devsim-biased device on
its own mesh with fv2d, at the device's current contact biases. `mc_coupling.lateral_cce_factor`
needs a scan that starts at x = 0. For the 3D electrode, scan from r = 0 so that the
low-CCE n+ column region is included, rather than starting at the top-gap edge as
notebook 19 does.

**Native meshes and adaptive refinement.**

- `etna.core.fv2d.structures` holds `guard_ring_structure`, `mesa_structure` and
  `electrode3d_structure`. They take the same keywords as the devsim builders (the
  ALT_STACK pin, graded epi). Their doping is a numpy twin of the devsim models, equal to
  round-off on devsim's mesh.
- The native box geometry is exact: 2π ∫r dr ring volumes. devsim's per-triangle
  `CylindricalNodeVolume` is about 1e-4 larger in the bulk and 5/3 larger on the axis.
- `refine.adaptive_cce_scan(structure, voltages, scan_kwargs, tol)` inserts tensor lines
  where the field, the space charge or a contact cell call for them. It re-solves the
  scan on every level and stops on a geometric-extrapolation error estimate.

**Accuracy and speed.** Measured in `lab/geometry2d/notes/LOG.md`:

- On devsim's mesh, parity with devsim is 1e-16 to 2.4e-12 in CCE and 1e-14 V in
  potential, for the mesa at 50 V, the guard ring at 10 V (including three-contact
  sharing) and the 3D electrode at 50 V.
- The full workflow (import or build, ramp, dark check, scan) is faster than devsim by
  5x (mesa), 13x (guard ring) and 27x (3D electrode) in one shared load window. Comparing
  each tool's best observed time, the factors are about 15x, 33x and 78x.

**Depth mesh at a contact on the epi (`depth_mesh`, since 2026-10-21).** The ALT_STACK
cathode is an ohmic contact directly on the n- epi, a perfect sink for holes: holes
generated within its boundary layer (the undepleted epi at partial depletion, a layer of
~kT/(qE) once depleted) are lost, and the CCE converges at FIRST order in the size of the
contact cell. The builders `create_mesa_device`, `create_guard_ring_device`,
`create_delta_e_e_device` and their fv2d twins therefore take
`depth_mesh="graded"` (default; `etna.core.contact_mesh`): a 0.1 nm cathode cell growing
by 1.25x to 50 nm, 50 nm over the last 1 um of the epi, and 1 nm junction lines.
`depth_mesh="legacy"` gives the meshes before 2026-10-21 (100 nm cathode cell, 10 nm
junction) for reproducing earlier numbers; a dict overrides single grading keys.
`create_sic_2d_device` takes the same keyword with default `"legacy"` (every pinned number
and the engine2d grid); use `"graded"` for the legacy-stack planar diode (`nplus_thickness_cm=None`).

| mesa centre CCE (ALT_STACK)                                                                                                                      | 10 V      | 50 V      |
| ------------------------------------------------------------------------------------------------------------------------------------------------ | --------- | --------- |
| converged (uniform-bisection limit, two starting meshes agree to 2.5e-8; `legacy_perp` mobility; the c-axis default gives 0.9630683 / 0.9990918) | 0.9630651 | 0.9990906 |
| graded mesh (default)                                                                                                                            | +5.4e-5   | -2.9e-6   |
| legacy mesh                                                                                                                                      | +6.17e-3  | +7.86e-4  |
| legacy-stack planar `create_sic_2d_device`, legacy mesh (27.5 nm cathode cell)                                                                   | +1.4e-3   | +6.8e-4   |

At 50 V the converged 1 - CCE (9.1e-4) is mostly this cathode sink; the default n+ stack,
which has no such contact, gives 1.25e-4. Cost: the guard ring goes from 11.5k to 20.4k
nodes (devsim 10 V ramp and scan about 110 s, fv2d 2.4 s). The default n+ stack is not
affected: bisecting the contact cells of the 1D engine moves the CCE by < 2.3e-10, a 0.1 nm
2D cathode cell by 4e-10 (the cathode sits on the 1e19 n+, no generation reaches it), so
the 1D/2D default meshes and the CCE(LET) tables are unchanged. Details, the 1D general
mesh error (C 2.7e-3 at -1 V, CCE 1.7e-4 at -1 V, < 1e-6 beyond -5 V) and every number
that moved: `lab/mesh_contact/notes/LOG.md`.

The AMR values quoted here before (0.96294 mesa, 0.96281 guard ring at 10 V, estimate
5e-6) were 1.2e-4 / 2.6e-4 off: `refine.adaptive_cce_scan`'s geometric-extrapolation
estimate is a heuristic and underestimated the error 15-20x on this pre-asymptotic
sequence. Use uniform bisection with Richardson (`lab/mesh_contact/bench/alt_cathode_study.py`)
for a reference value.

**Validity of the axisymmetric proxy.** The (r, z) 3D electrode with a p+ wall is a
faithful model of annular-trench cells. It is not a model of column-lattice 3D
detectors. For a 50 µm square cell with r = 5 µm columns:

- 99 % depletion is reached at 75 V, against 30 V for the equal-area cylinder;
- the square cell keeps low-field saddles (about 100-400 V/cm at 50-100 V) that the
  cylinder cannot have.

For such layouts, solve the Cartesian cross-section (`lab/geometry2d/bench/axisym_validity.py`).

**3D.** True 3D DD with direct solvers is not tractable here. Extrapolated to a 0.9 M-node
3D cell, it would take about 780 s and 31 GB per factorisation. A 3D engine needs
Krylov + AMG, see the LOG.

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

| notebook                       | topic                                          | notes                                                                                                                                                                                     |
| ------------------------------ | ---------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 01_phase1_validation           | material parameters, analytical electrostatics |                                                                                                                                                                                           |
| 02_electrical_characterization | C-V, Mott-Schottky, doping profile, I-V        |                                                                                                                                                                                           |
| 03_charge_collection           | CCE vs bias, generation profiles, Hecht        | generator `scripts/create_notebook_03.py`                                                                                                                                                 |
| 04_flash_recombination         | FLASH Auger/SRH CCE vs dose rate               | sensitivity bound only                                                                                                                                                                    |
| 05_dark_current_vs_fluence     | dark current vs proton fluence                 | TAT on the batched engine; fluences are fractions of the legacy `legacy_compensation_fluence(DeviceConfig(), 62)` (4.36e12 p/cm², sic-damage; 3.48e12 si-neq before 2026-10-21) since 2026-10-21                                                          |
| 05b_flash_parametric_studies   | FLASH parametric sweeps                        | runs its 60-condition sweep when the untracked cache is missing                                                                                                                           |
| 06_temperature_dependence      | material properties and CCE vs T               |                                                                                                                                                                                           |
| 07_dark_current                | TAT / SRH / surface terms                      |                                                                                                                                                                                           |
| 08_transient_flash             | FLASH pulse transients                         | generator `create_notebook_08.py`                                                                                                                                                         |
| 09_radiation_damage            | defect introduction, lifetimes, NIEL scaling   | NIEL plot in both conventions (sic-damage default 1.552, si-neq 1.942 at 62 MeV); summary table at 0.1/0.5/0.9 x the legacy compensation fluence; note cell on the 2026-10-21 defaults                                           |
| 10_cce_vs_fluence              | CCE degradation vs fluence                     | Figure 1b: the default calibrated kinetic trap model vs the legacy removal, 0.01-100 x its compensation fluence (1.80e12 p/cm²); the other figures study the legacy model; fluences are fractions of the legacy `legacy_compensation_fluence(DeviceConfig(), 62)` (4.36e12 p/cm², sic-damage; 3.48e12 si-neq before 2026-10-21) since 2026-10-21                                                                                     |
| 11_dark_current_cv_evolution   | leakage and C-V vs fluence                     | fluences are fractions of the legacy `legacy_compensation_fluence(DeviceConfig(), 62)` (4.36e12 p/cm², sic-damage; 3.48e12 si-neq before 2026-10-21) since 2026-10-21; C-V at 0.1/0.3/0.6/0.9 x                                                           |
| 12_multi_defect_comparison     | 3-defect vs single-defect models               | fluences are fractions of the legacy `legacy_compensation_fluence(DeviceConfig(), 62)` (4.36e12 p/cm², sic-damage; 3.48e12 si-neq before 2026-10-21) since 2026-10-21                                                                                     |
| 13_parametric_optimization     | radiation-hardness sweep                       | 2 x 3 x 2 grid at 0.6 x the compensation fluence of the lowest-doped config (was 4 x 4 x 4 at 1e13, 15-20 min)                                                                            |
| 14_validation                  | model vs published irradiation trends          | fractions of the legacy `legacy_compensation_fluence(DeviceConfig(), 62)` (4.36e12 p/cm², sic-damage; 3.48e12 si-neq before 2026-10-21) + 1 refused point; trapping-only loss below it; unsourced 0.95/0.70/0.40 curve and onset ~1e12 removed 2026-10-21 |
| 15_2d_electrostatics_cce       | 2D field and CCE of the microdosimeter         | generator `create_notebook_15.py`                                                                                                                                                         |
| 16_single_particle_cce         | single-particle 2D transient CCE, CCE(LET)     | generator `create_notebook_16.py`; loads `data/cce_let_table_*.json`; `ETNA_REBUILD_CCE_TABLES=1` rebuilds them (~40 min) into `$ETNA_FIGURE_DIR/data/`                                   |
| 17_mc_coupling                 | MC energy deposition to pulse height           | generator executes it in place                                                                                                                                                            |
| 18_microdosimetric_spectra     | lineal energy spectra, y_F / y_D               | real NIST kappa, 62 MeV assumed; generator executes it                                                                                                                                    |
| 19_alternative_structures      | mesa, 3D electrode, ΔE-E, guard ring           | pins the previous-default doping on the legacy geometry (PERF-06)                                                                                                                         |
| 20_feasibility_report          | microdosimeter design study                    | outside-validity rows reported                                                                                                                                                            |

Every notebook calls `set_devsim_quiet(True)` (`etna.core.devsim_lock`, exported from
`etna`) in its first code cell: devsim's Newton log of the direct core calls is dropped
(it used to stream ~1e5 lines per notebook into the outputs); `ETNA_DEVSIM_VERBOSE=1`
brings it back. Use it in your own scripts too.

Run times measured on 2026-10-21 (Apple M4, `jupyter nbconvert --execute` from `notebooks/`
into a scratch directory with `ETNA_FIGURE_DIR` in scratch; wall clock including the kernel
start, numba cache warm, load average 2-3): 01 3.7 s, 02 12 s, 03 4.8 s, 04 3.7 s, 05 3.4 s,
05b 4.2 s (full 60-condition sweep on the engine, empty cache; it still prints its old "~80 s"
estimate), 06 32 s, 07 3.5 s, 08 9.1 s, 09 3.2 s, 10 4.5 s, 11 4.6 s, 12 3.6 s, 13 3.2 s,
14 3.3 s, 15 34 s, 16 53 s (tables loaded; one 2D ion-track transient remains), 17 2.9 s,
18 2.9 s. Notebooks 19 and 20 were not re-measured (over 2 min each: 193 s and 123 s on
2026-10-18). Notebooks 05, 05b, 06, 07, 20
were re-timed on 2026-10-19 after the dark-current / FLASH engine port, back to back with the
previous code (scratch copies, same session): 05 35.9 -> 6.5 s, 05b 91.2 -> 5.3 s, 06 50.6 ->
28.5 s (its C-V / CCE-vs-T sweeps are still devsim), 07 27.2 -> 5.6 s (its
dark_current_sweep on a user-built devsim device stays devsim), 20 271.6 -> 123.1 s (the
2D dark current; the rest is the 2D CCE scans). Rebuilding the
notebook-16 tables (`ETNA_REBUILD_CCE_TABLES=1`) takes about 40 min (2026-10-18, not re-measured; 100 um ~15 min, 300 um ~26 min; one device per table, linear-response synthesis of the non-simulated LETs) and writes them to
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

**Rule:** never commit outputs (`uv run nbstripout notebooks/X.ipynb`). Run notebooks from
`notebooks/` (or a scratch copy) so their figures stay out of the repo root.

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
| `lab/geometry2d`        | the general 2D engine `etna/core/fv2d`: parity, AMR, axisymmetric proxy, 3D tractability             |
| `lab/mesh_contact`      | contact-cell mesh convergence: graded depth mesh of the contact-on-epi builders, what moved          |
| `lab/edge_cce`          | lateral (edge) CCE map of the guard ring for the MC weighting                                        |
| GPU spike (archived)    | tag `archive/gpu-solver-lab`: numba beat the GPU for 1D batches                                      |

Where decisions and evidence live:

- `.issues/`: issues and work packages. `.issues/perf/PERF-0x-*.md` covers performance
  and the default stack, `.issues/audit-v7/` the latest audit (one file per finding,
  indexed in `*-INDEX.md`), and `.issues/TECH-DEBT-2026-10-18.md` the open-debt inventory.
- `.planning/`: physics audits (`PHYSICS_AUDIT_v*.md`, `PHYSICS_REVIEW_v6.md`), the
  literature check (`VALIDATION_LITERATURE.md`) and milestone audits.
- `deliverables/`: foundry-facing specs (DESIGN-1 dosimetry p-n, DESIGN-2 microdosimeter
  geometry; Markdown is the source and the PDFs were regenerated from it on 2026-10-21).
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

- 136 modules and 1766 tests (2026-10-21). `pytest.ini` deselects the `exhaustive` marker
  (6 redundant deep variants, each implied by a kept test), so the default run collects
  1760: all pass with `--extra dev --extra viz`; with `--extra dev` only, 1752 pass and 8
  skip (the manim and gdstk tests). The `slow` marker labels the devsim integration tests
  but does not deselect them.
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
- **Timing** (10-core M4, 4 performance + 6 efficiency cores; the machine was shared with
  other jobs). 2026-10-21, 1508 default tests: 169 s wall, 902 s user CPU (load average
  3 at the start); 205 s wall, 948 s user CPU with load average up to 25. Expect about
  3 minutes. History: 2026-10-18 before the
  scheduler and leak fixes, `--dist loadfile`, 502 s wall, 2173 s user CPU; after,
  248 s wall, 1682 s user CPU. The engine-seeded devsim start (`seeded_bias`) and the
  test-speed work of 2026-10-19..21 cut the CPU time further. The floor is set by the
  slowest single tests, which run ~2x slower inside a loaded 10-worker run.
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
- CI (`.github/workflows/ci.yml`) installs `--extra dev --extra viz` (apt cairo/pango so
  pycairo and manimpango build) and runs the lock check, ruff (`ruff check .` with
  `ruff.toml`: pyflakes F, E9, W6 over etna/, app/, tests/, scripts/), notebook hygiene, the generator check and the default run with
  `-n 4`. It runs only on pushes to `main`, pull requests and manual dispatch, not on
  every branch push, because the local run is the gate. The numba cache is kept in
  `actions/cache`, keyed on the exact hash of `etna/**/*.py` + `uv.lock`, with source
  mtimes pinned so numba accepts it. The test timeout is 45 min; the test job took 14 min on 2026-10-21
  (37-52 min before the leak fix). A nightly job runs
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

| symptom                                                            | cause and fix                                                                                                                                                                                                                                                                           |
| ------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| A devsim device you built with `etna.core.*` disappeared           | devsim-path facades (`run_cv` fallback, `run_field` / `run_dark_current` outside the engine box) call `reset_devsim_fully()` (`etna/core/devsim_reset.py`), which deletes **every** devsim device in the process. Build your own devices after the facade calls, or in another process. |
| results that depend on a previous cell (after an interrupted call) | a stale device from an interrupted call. Restart the kernel, or call `etna.core.devsim_reset.reset_devsim_fully()`. `run_all` asserts that no device leaks between sections.                                                                                                            |
| `metadata["truncated"] == True`, shorter `x`                       | a devsim convergence wall at deep reverse bias. For the legacy model it is the cathode-on-epi artefact (about -64 V): use the default stack, or stay above it.                                                                                                                          |
| `OutsideModelValidity`                                             | the requested point is outside the model's validity (legacy 2D past the contact wall; `run_cce_tilted_alpha` outside the engine box). The message says what to change.                                                                                                                  |
| NaN points in a damage sweep                                       | legacy model (`radiation_traps=False`, or "auto" outside the engine: `metadata["damage_model_note"]`): the fluence is at or past `compensation_fluence`; see `metadata["reason"]`.                                                                                                                                                                                                        |
| First call is slow (about 0.5 s), later calls about 10 ms          | numba compilation (measured 2026-10-21, run_cv; section 6). The batched-engine kernels use `cache=True`, so later processes load the compiled code from `__pycache__`. Changing numba or the source recompiles.                                                                         |
| A C-V or CCE call takes seconds instead of milliseconds            | the config is outside the engine box (`validated_domain_reason(cfg)`), or the engine fell back (`metadata["fallback_reason"]`).                                                                                                                                                         |
| `ModuleNotFoundError` (plotly, streamlit) under `uv run pytest`    | `.venv` was synced without `--extra dev`, so a system pytest ran. Run `uv sync --extra dev`.                                                                                                                                                                                            |
| Export page: "Animations are unavailable" / no GDSII button        | install `--extra viz` (macOS: `brew install pkgconf cairo pango` first). Note that a later plain `uv sync --extra dev` removes it again.                                                                                                                                                |
| Notebook figures appeared in the repo-root `figures/`              | it ran from the repo root. Run from `notebooks/`, or set `ETNA_FIGURE_DIR`. `figures/` is gitignored; if the tracked `figures/cv_validation_bruzzi2023.png` changed, restore it with `git checkout figures/`.                                                                           |
| You want devsim's solver output                                    | `ETNA_DEVSIM_VERBOSE=1`.                                                                                                                                                                                                                                                                |

---

## 17. Measurement inputs (`etna.measurements`)

The owners (INFN-LNS, the foundry) will supply C-V to -100..-200 V, SIMS/SRP through the
epi/n+ interface and I(T) (`.issues/TECH-DEBT-2026-10-18.md`, "data requests"). This module
reads them as editable CSV inputs, validates them, compares them with the model and turns
them into model inputs. It does not change the physics: the model side is `run_cv`,
`run_cv_batch` and `run_dark_current` on a `DeviceConfig`. **The repository ships no
measurement data**: `data/measurements/templates/` holds header-only templates and
`data/measurements/README.md` lists every column, unit and metadata key. Every test uses
synthetic data generated in the test from the model itself.

File rules (all three kinds): units are in the column names (`C_pF`, `depth_nm`, `I_pA`,
`T_K`); unknown names are refused with the accepted list; `# key: value` lines carry
metadata (an empty value = not given); bias is in the ETNA sign (reverse negative) unless
`# bias_convention: reverse_positive`. A C-V with more forward than reverse points, or
forward bias above 3 V, and an I(T) with any forward bias, are refused as probable sign
errors. Every problem of a file is reported at once by one
`MeasurementFileError` (a `ValueError`), each line naming the file line (1-based, counting
comments and header) and the column, e.g. `line 4, column 'C_pF': 'abc' is not a finite
number`. Every result has `to_json()` (with the data file's sha256 and the model's config
hash, solver and ETNA version), `notes` and `limits`; `report_markdown(...)` writes one
report; `etna.measurements.plots` draws the figures in the `etna.viz` style (export with
`etna.viz.export_figure`).

### C-V

```python
import etna, etna.measurements as M

d = M.load_cv_csv("my_cv.csv")              # area: column > '# area_cm2:' > area_cm2=...
cmp = M.compare_cv(d, etna.DeviceConfig())  # W = eps*A/C and C/A on the measured biases
cmp.rms_W_um, cmp.rms_C_rel, cmp.notes
dep = M.depletion_profile(d)                # differential N(W), depletion approximation
inv = M.invert_cv(d, etna.DeviceConfig())   # -> inv.donor_table, inv.config
open("donor_table.json", "w").write(inv.to_json())
cfg = inv.config                            # the fitted device itself
print(M.device_rebuild_code(inv.config, "donor_table.json"))  # code rebuilding it later
```

- `DeviceConfig(donor_table=table)` alone is NOT the fitted device: it keeps the
  defaults (epi, substrate, n+, mobility, traps) and the parametric transition. Use
  `inv.config`, or `device_rebuild_code(inv.config, json_name)`, which writes every
  non-default field out (issue #20; the Measurements page shows that code).
- Capacitance is stored per area (F/cm²); `W = EPS / (C/A)` with the model's permittivity,
  so data and model are converted identically. The model area does not enter the
  comparison (per-area C); a data area different from `config.area_cm2` is noted.
- `compare_cv` runs `run_cv` from 0 V to the deepest measured bias. When the measured
  biases lie on a uniform grid from 0 V the model is solved **at** them (exact); otherwise
  on a 0.05 V grid and interpolated linearly (max interpolation error < 0.01 µm at the knee
  of the default device, `step_V` to refine). The file's temperature replaces `config.T`
  (`use_data_temperature=False` to keep it). Forward-bias points are kept in the data but
  not compared.
- `depletion_profile`: `N(W) = 2 / (q eps d(1/C'^2)/dV_rev)`, slope of a local linear fit
  over `window` points. Model-free and approximate: it smears N over a Debye length,
  is wrong where N changes fast (junction layer, n-/n+ rise) and diverges at full
  depletion. It is the start of the inversion.
- `invert_cv` is the `lab/cv_gap` knot method (`fitlib.fit_knots`, which produced the
  19-knot `MEASURED_CV_TABLE`) on the public API: log10 N_D at fixed knots (default: the lab
  knots scaled to the epi, dense at the junction and the epi end; `knots_um=` to choose),
  log-linear between them, least squares on `(W_data - W_model)/sigma_W` plus
  `lam * diff2(log10 N) / sqrt(knot spacing)` (`lam` = 1). Every trial is a full Poisson +
  small-signal C solve on the batched engine; the knots are bounded to its validated table
  range (5e12..5e16 cm⁻³), a device outside the engine box is refused. Without a sigma
  column every point gets `default_sigma_W_um` = 0.02 µm and chi2 is relative.
- Runtime scales with the model grid. Data on a uniform grid from 0 V are solved at their
  own biases: 241 points 0..-60 V 36 s, 401 points 0..-200 V 30 s (synthetic, M4).
  Non-uniform data use a 0.05 V grid capped at `max_grid_points` = 801 (noted); an
  uncapped 4001-point grid (0..-200 V) took 600 s.
- Solver path (a deliberate deviation from "public API only"): `run_cv_batch` solves one
  donor table per group, so the K+1 Jacobian columns would be K+1 separate solves. The
  inversion stacks them into one `etna.core.batched` PoissonBatch built with the same calls
  as the facade (`device_kwargs` -> `make_batch` -> `reverse_sweep`), used only after it
  reproduced `run_cv` on the start profile to 1e-9 µm (otherwise the public path is used,
  noted); the reported fit is a `run_cv` call. Pinned by
  `test_stacked_batch_matches_run_cv_per_table`.
- Round trip (tests/test_measurements_cv.py, noise-free synthetic curve from the 19-knot
  preset, 121 points 0..-60 V, inverted on the default geometry): converged, rms W < 5 nm,
  every knot inside the probed depth within 0.1 dex (measured: 0.03 dex at 241 points).
  This is an inverse crime (same model, no noise, the true knots): it checks the
  machinery, not the accuracy on a real device.

Limits of the C-V route (also in every result's `limits`):

- The C-V constrains N_D only between W at the least-reverse bias used (about W(0 V) =
  1.7 µm on the dosimeter) and W at the deepest bias. Knots outside are set by the
  regularisation (`inv.constrained`, `inv.probed_depth_um`). Minimum bias range: the
  deepest probed depth is W at the deepest bias. On the 19-knot preset (synthetic)
  eps*A/C = 9.26 / 9.53 / 9.74 / 9.89 / 9.99 / 10.01 µm at -5 / -10 / -30 / -60 / -100 /
  -200 V: the bulk (2-9 µm) needs only -5 V, the n-/n+ rise (9.25-10 µm) is probed between
  -10 and -100 V, and -100..-200 V adds about 20 nm at the epi end (where the rise reaches
  the 5e16 table bound). A curve that stops at the knee says nothing about the epi end.
- Depth resolution is about one Debye length (`inv.debye_um`: 0.37 µm at 1e14, 1.2 µm at
  1e13 cm⁻³); finer knots are regularisation-driven (noted).
- Area error and perimeter (fringe) capacitance are degenerate with doping: a 2.5 % area
  change or a fringe charge reproduce the same -30 V shift (`lab/cv_gap/notes/LOG.md`).
  Give the measured area and, ideally, guard-ring-connected C-V.
- Quasi-static model: no frequency dispersion, no trap response, no series resistance. The
  measurement frequency is recorded only.
- Beyond `cross_checked_to_V` (100 V with the n+ layer) the engine is not cross-checked
  against devsim; the -100..-200 V part of a new C-V is compared but flagged.

### SIMS / SRP

```python
prof = M.load_sims_csv("my_sims.csv")        # depth from the SURFACE, N and/or Al in cm^-3
tab = M.resample_profile(prof, etna.DeviceConfig(), spacing_um=0.25)
tab.donor_table, tab.rel_uncertainty, tab.censored, tab.junction_source, tab.N_A_plateau
cfg = tab.to_device_config(etna.DeviceConfig())   # table + p+ thickness (+ N_A plateau)
M.profile_vs_model(tab, etna.DeviceConfig())["rms_log10"]
```

- Depth origin: `donor_table` depth is below the p+/n- junction. The junction is the
  `junction_depth_um` argument, else the file's `# junction_depth_um:`, else the Al/N
  crossing, else the model's p+ thickness (`tab.junction_source` says which). The table
  stops at the config's epi end; the n+ substrate stays `N_D_nplus`.
- Epi end as measured: `tab.epi_end_from_data_um` is the depth below the junction where
  N reaches the n+ level (0.1 x `N_D_nplus`). When it differs from the config's epi by
  more than a knot spacing a note says so (the table then holds n+-level values, which
  `DeviceConfig` refuses above N_A or runs on devsim above 5e16): resample with
  `epi_thickness_um` set to it.
- `to_device_config` checks the result against the engine box (`tab.domain_reason`,
  noted): e.g. a junction shallower than 0.3 µm or an Al plateau above 1e19 runs on devsim.
- Resampling: `method="bin_mean"` (default; mean of log10 N in the bin around each knot,
  noise-robust, rounds structure narrower than a bin) or `"interp"` (log-linear
  interpolation, exact on a noise-free profile). Per-knot relative uncertainty: the mean
  measurement sigma in the bin (not divided by sqrt(n): SIMS calibration errors are
  correlated) combined with the bin scatter / sqrt(n).
- Detection limit (column or metadata): knots at or below it are `censored` and
  `to_device_config` refuses them unless `allow_censored=True`, so a SIMS floor never
  becomes the epi doping silently. Knots outside the measured depth are `extrapolated`
  (edge value held).
- Acceptors: `DeviceConfig` has a scalar `N_A` and a p+ thickness only, so the Al profile
  is reported as data, the plateau (median Al above the junction) is offered as `N_A`, and
  the Al tail into the epi is not modelled.
- SIMS measures total atoms, SRP the net free-carrier density; the technique is recorded
  and both are used as N_D.

### I(T) and the activation energy

```python
it = M.load_it_csv("my_it.csv")              # one fit per bias (V_V column or '# bias_V:')
an = M.analyze_it(it, etna.DeviceConfig(), T_power=0.0)
for c in an.comparisons:
    print(c.bias_V, c.data_fit.E_a_eV, c.data_fit.sigma_E_a_eV, c.model_fit.E_a_eV, c.z)
```

- Fit: `ln(|I|/T^p) = ln A - E_a/(k_B T)`, weighted least squares with `sigma_ln =
sigma_I/|I|` (absolute covariance and chi2/dof with file sigmas; residual-scaled
  covariance and >= 3 temperatures without). `p` = `T_power` is a choice (0 = plain
  Arrhenius, 1.5 = generation current, 2 or 3 as in parts of the SiC literature); E_a moves
  by about k_B T per unit of p, so data and model are fitted with the same p.
- Model E_a: `run_dark_current` at the data temperatures and bias, the same fit on
  |I_model(T)|. Its T dependence is an assumption (`N_t(T) ~ n_i(T)`, i.e. E_a = E_g/2,
  plus the T^1.5 of n_i, the Varshni E_g(T) and the Hurkx Gamma(T)): on the default device
  at -30 V, 280-340 K, p = 0 the fitted model E_a is 1.69 eV against E_g/2 = 1.63 eV. The
  comparison (`delta_eV`, `z` = delta / sigma of the data fit) tests that assumption; a
  curved Arrhenius plot (see the residuals) means more than one mechanism.

### Unit and plausibility checks, examples, data request

On top of the schema checks, every loader runs `etna.measurements.plausibility` on the
converted values. **Refused** (added to the `MeasurementFileError` with a unit hint): a
median depletion width eps/(C/A) below 1 nm or above 1 cm ("values look like F in a pF
column (x1e12)", "pF in an F column"), an area outside 1e-6-100 cm^2, |V| > 5 kV, a SIMS
concentration above the SiC atomic density (m^-3 in a cm^-3 column) or all values <= 30
(log10 given), depths beyond 5 mm (nm in `depth_um`), T outside 4-2000 K, |I| > 0.1 A.
**Warned** (`MeasurementPlausibilityWarning`, also `plausibility_report(data)`): W
outside 0.05-500 um, |V| > 1.5 kV, test frequency outside 10 Hz-100 MHz, T outside
77-700 K (with a degC hint for a `T_K` column below 250), C rising by > 5 % with reverse
bias, SIMS values below 1e12 cm^-3 or a profile ending above 0.2 um, |I| above 0.1 mA or
below 10 aA. The bounds are sanity limits, not model validity ranges.

`examples/measurements/` holds one SYNTHETIC file of each kind generated from the model
(`make_synthetic.py`, truth = `measured_cv_tabulated_config()`: C-V 0 to -200 V, SIMS
N/Al from the surface, I(T) at -10/-50 V, 300-380 K) and `roundtrip.py`, which loads them
and recovers the truth (C-V rms W 8e-10 um, donor table 0.000 dex, E_a to 1e-11 eV) in
under a second; `tests/test_measurement_examples.py` checks the files are reproducible
and the round trip. With a realistic SIMS N floor above the epi doping the epi knots are
censored and refused (the test shows it). What to ask INFN-LNS and the foundry for, and
how each file is used: `docs/DATA_REQUEST_CHECKLIST.md`.

### Streamlit page

**Measurements** (`app/workflows/measurements.py`): three tabs (C-V, SIMS / SRP, I(T)),
each with the header-only template download, a file upload, the validation messages (line
and column), and buttons that run the comparison / inversion / resampling / E_a fit against
the sidebar device. A C-V file without an area uses the sidebar area, and the page says so.
Outputs: figures (PNG), the comparison and E_a JSON, the donor-table JSON
(`etna.donor_table/1`, readable by `load_donor_table_json`) and a Markdown report. Results
are cached in the session and flagged when the file or the sidebar device changed.

---

## 18. Physical-model library (opt-in): high-field mobility, anisotropy, avalanche and breakdown

**What this adds.** `etna.core.physics_models` adds models on top of the default physics (issue #6):

- velocity saturation;
- the c-axis mobility anisotropy;
- impact ionization, with breakdown and multiplication;
- registered bandgap-narrowing parameters.

Each parameter set records its provenance. **Nothing changes unless you ask** apart from the low-field mobility set, which is a `DeviceConfig` field read by every engine (default `c_axis` since 2026-10-21, next subsection; `mobility="legacy_perp"` restores the earlier numbers bit for bit). The other models are reached only through the three facades below. They run on the batched 1D engine, which is checked against a devsim twin.

```python
from etna import DeviceConfig, PhysicsModels, run_breakdown, run_avalanche, run_transit_time

cfg = DeviceConfig()                                   # default 4 mm2 dosimeter
bd = run_breakdown(cfg)                                # default set: hatakeyama2004 along c
print(round(bd.metadata["V_BD"], 2), round(bd.metadata["mesh_error_V"], 3))   # 2225.78 0.105

m = PhysicsModels(impact_ionization="hatakeyama2004", high_field_mobility="lades2000")
a = run_avalanche(cfg, m, v_start=-100.0, n_points=6, approach_breakdown=True)
print(a.y.round(3))   # CCE with gain (not clipped): [1. 1.074 4.053 21.827 95.504 388.867]
print(a.metadata["M_mcintyre"].round(3))  # McIntyre limit, same field: [... 95.57 389.781]

t = run_transit_time(cfg, PhysicsModels(high_field_mobility="lades2000"), -10, -100, 4)
print((t.metadata["t_e"] / t.metadata["t_e_low_field"]).round(2))   # [1.44 2.87 4.8  7.1 ]
```

### Default low-field mobility: parallel to the c-axis (since 2026-10-21)

The current of the vertical (0001) diode flows along c. Until 2026-10-21 every solver used the Schaffer/Ayalew Caughey-Thomas set, which was measured **perpendicular** to c (Ayalew Table 3.5 is headed mu_perp; its eq. 3.81 gives mu_n,par = 1.2 mu_n, mu_p,par = mu_p). `DeviceConfig.mobility` now selects the set for every engine (devsim 1D/2D, batched 1D, engine2d, track, dark, FLASH, fluence):

| `mobility=`                   | Electrons                                                                                                                    | Holes                                                                   | Use                                                      |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------- | -------------------------------------------------------- |
| `"c_axis"` (default)          | Ishikawa 2025 thesis Eq. 2.16, drift mobility parallel to c, at 300 K (1190 / 10 / 1e17 / 0.54); T law gamma -2.40 on mu_max | unchanged (legacy set, isotropic as Ayalew eq. 3.81)                    | vertical transport                                       |
| `"legacy_perp"`               | Schaffer/Ayalew perpendicular (950 / 40 / 1.94e17 / 0.61, gamma -2.40)                                                       | Schaffer/Ayalew (125 / 15.9 / 1.76e19 / 0.34, gamma -2.15)              | the earlier numbers, bit for bit; the historical presets |
| `"c_axis_holes_0p9"` (opt-in) | as `c_axis`                                                                                                                  | legacy x 0.90 (thesis drift ratio mu_par / mu_perp, one opened primary) | the hole-anisotropy alternative                          |
| `"ishikawa2025"` (opt-in)     | Eq. 2.16 with its full T law (fitted 300-600 K, warns outside)                                                               | thesis Eq. 4.6 (Hall, parallel to c)                                    | sensitivity to the T law and the hole magnitude          |

Mobility at the epi bulk (3.38e13), 300 K: electrons 945.4 -> 1174.4 cm2/Vs (+24.2 %, uncertainty about +-8 %, the Hall-factor conversion); holes 123.8, unchanged. The primaries agree on the electrons (the thesis's perpendicular drift fit equals the legacy set within 3.5 % up to 5e15; vertical-SBD drift 1070 +- 290 at 1e15; mu_par / mu_perp 1.16-1.21). They do not agree on the electron T law (drift fit -3.4 vs Hall -2.4 to -3.2: -2.40 kept) nor on the holes (magnitude 83-124 cm2/Vs Hall at low N_A; anisotropy 0.9 Ishikawa, 1/1.15 Hatakeyama 2003, isotropic Ayalew), so the holes stay as they were (uncertainty along c about -35 / +0 %). Sources, cross-checks and decisions: `lab/mobility_caxis/notes/LOG.md`.

Effect on the default dosimeter (`lab/mobility_caxis/bench/caxis_bench.py`, before = `legacy_perp`, after = `c_axis`):

| Output                                                                  | Before                                      | After                                                                                   |
| ----------------------------------------------------------------------- | ------------------------------------------- | --------------------------------------------------------------------------------------- |
| C-V, W(V), every C-V fit (Photons Fig. 4b labels, Bruzzi), V_BD         |                                             | identical (electrostatics)                                                              |
| alpha CCE (`run_cce`) at -2 / -10 / -40 / -100 V                        | 0.993707 / 0.999429 / 0.999845 / 0.999936   | 0.993716 / 0.999433 / 0.999847 / 0.999937 (+8.5e-6 to +6.8e-7, mesh-converged to 0.3 %) |
| electron transit, full depleted width, -10 / -40 / -100 V               | 125.1 / 34.0 / 11.3 ps                      | 103.6 / 28.6 / 9.4 ps (x0.83)                                                           |
| hole transit, same                                                      | 935.6 / 249.5 / 84.6 ps                     | unchanged                                                                               |
| dark current at -20 / -60 / -100 V                                      | 7.637e-12 / 1.800e-11 / 2.487e-11 A         | same to < 2e-6 (generation-dominated; the 18 pA N_t calibration holds)                  |
| CCE vs T at -30 V, dd, 250 / 300 / 450 K                                | 0.999813 / 0.999797 / 0.999757              | 0.999814 / 0.999799 / 0.999759                                                          |
| Hecht CCE vs T at -30 V, 300 K                                          | 0.994405                                    | 0.995475                                                                                |
| tilted-alpha CCE (microdosimeter preset) at 0 V; Petringa Fig. 4b rms   | 0.87545; 0.0354                             | 0.87545; 0.0354                                                                         |
| proton-damage CCE at 2.6e11 / 6e11 cm^-2, -40 V                         | 0.99930 / 0.99874                           | 0.99932 / 0.99880                                                                       |
| charged-trap CCE at 1e13 / 1e14 / 1e15 p/cm2, -40 V                     | 0.996109 / 0.972203 / 0.354958              | 0.996417 / 0.974836 / 0.357614                                                          |
| FLASH steady-state CCE, -30 V                                           | 0.999809                                    | 0.999811                                                                                |
| FLASH transient CCE, 20 / 100 Gy/s                                      | 0.999725 / 0.999869                         | 0.999725 / 0.999869                                                                     |
| track engine (opt-in), LET 10 / 1000, 50 V: t95 of the cathode charge   | 121.4 / 232.3 ps                            | 120.7 / 237.2 ps                                                                        |
| CCE(LET) table, 2D sheet track, hw 50 um, 50 V (tracked `data/` tables) | 0.9999625-0.9999703, t_collection 0.1168 ns | 0.9999662-0.9999742 (+3.7e-6 to +3.9e-6), t_collection 0.1168 ns                        |

The legacy CCE pins hold within 1e-5; only the transit times and the trapping-limited collection (charged traps, Hecht) move. `"c_axis_holes_0p9"` gives CCE 0.993031 at -2 V and a hole transit x1.111 (the hole diffusion out of the undepleted tail sets the CCE loss). 2D devices take the scalar c-axis mobility in both directions: vertical transport is right, lateral electrons are about 20 % too mobile (no mobility tensor yet). Tests: `tests/test_mobility_caxis.py` (registry, cross-checks against the opened primaries, devsim-expression parity, legacy reproduction, engine/devsim parity of the selection).

### The selection: `PhysicsModels`

`PhysicsModels` is frozen and hashable. Its fields take a registry name or a parameter-set instance:

| Field                 | Registry names                                                                                                                                                                      | Default                       |
| --------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------- |
| `low_field_mobility`  | `etna_default` (the device's own set, `DeviceConfig.mobility`), `c_axis`, `legacy_perp`, `ishikawa2025`, `ayalew2004` (legacy_perp plus Ayalew's (T/300)^-0.5 on mu_min)            | `etna_default`                |
| `high_field_mobility` | `lades2000`: mu = mu0 / (1 + (mu0 E / v_sat)^beta)^(1/beta), v_sat ~ T^-0.44, beta ~ T^1; `ishikawa2025`: electrons v_sat 1.7e7 isotropic, beta 0.76, no T law (holes as lades2000) | None                          |
| `anisotropy`          | `ayalew2004` (mu_n,par = 1.2 mu_n), `burin2024` (electrons x 1/0.83, holes x 1/1.15); only with a perpendicular low-field set (`legacy_perp`, `ayalew2004`)                         | None                          |
| `impact_ionization`   | `hatakeyama2004` (default of `run_breakdown`), `niwa2015`, `zhao2019`, `loh2008`, `stefanakis2021`, `ayalew2004_legacy`; all are alpha = a (1 + c dT) exp(-[b (1 + d dT) / E]^m)    | None                          |
| `bandgap_narrowing`   | `lades2000` (Slotboom: Delta E_g = C (F + sqrt(F^2 + 0.5)), F = ln(N / 1e17))                                                                                                       | None, registered only (below) |

**Validation is strict.** These are all errors:

- an unknown or miscapitalised name;
- a wrong type;
- an impact-ionization set without parallel-to-c coefficients, because the 1D diode's field is along c;
- a temperature outside the set's measured range, unless `allow_T_extrapolation=True`. Hatakeyama 2004 is 300 K only; when extrapolated, the results carry `T_extrapolated`.

`describe()` gives the provenance record that every result stores in `metadata["physics_models"]`.

**Direction.** Since 2026-10-21 the solvers' default low-field set is `c_axis`, the mobility PARALLEL to c, the direction of the vertical current (see "Default low-field mobility" above). `legacy_perp` (`DeviceConfig(mobility="legacy_perp")`) is the previous default, the Caughey-Thomas set measured perpendicular to c (Schaffer 1994 via Ayalew Table 3.5), which ETNA used for vertical transport until then. An `anisotropy` model rotates a perpendicular set to the c-axis and is refused with a c-axis set (it would double-count). The high-field model uses `v_sat,par` whenever the vertical transport is parallel to c (`PhysicsModels.vertical_is_parallel`: a c-axis device, or an anisotropy), `v_sat,perp` on a `legacy_perp` device without anisotropy (the previous behaviour, bit for bit).

### Provenance (what was opened)

Status **primary** means the numbers were read in the primary source. Status **secondary** means they were read in a compilation that restates it:

- J. Burin et al., "TCAD parameters for 4H-SiC: a review", arXiv:2410.06798, Tables XIII, XXVII, XXXIII and XXXIV;
- the TU Wien thesis of T. Ayalew (2004).

| Set                                                 | Status                                                                                                          | Notes                                                                                                                                                                                                                                                           |
| --------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `hatakeyama2004`                                    | primary: APL 85, 1380 (2004), doi:10.1063/1.1784520, Table I                                                    | <0001> and <11-20>; 300 K; fitted 2-5 MV/cm. The SISPAD-2005 anisotropic interpolation is `hatakeyama_anisotropic_alpha` (tested, not coupled into 2D)                                                                                                          |
| `niwa2015`, `zhao2019`, `loh2008`, `stefanakis2021` | secondary (review Table XIII)                                                                                   | Primaries not opened. The review flags Zhao 2019's hole c, d. The `niwa2015` 300-448 K and `stefanakis2021` 300-573 K ranges are ETNA assumptions. The `loh2008` direction is not stated in the review                                                          |
| `ayalew2004_legacy`                                 | primary, **legacy**                                                                                             | The review finds b_n carried over from 6H. It is off by +48 to +94 % against Hatakeyama's measured V_BD below. Do not use it for predictions                                                                                                                    |
| `c_axis` mobility (default)                         | electrons primary (Ishikawa 2025 Kyoto thesis Eq. 2.16, drift, parallel to c); holes: `legacy_perp` (unchanged) | 300 K values exactly Eq. 2.16; T law: gamma -2.40 on mu_max (Schaffer's <0001> exponent, secondary) instead of the thesis's -3.4 (a disagreement between drift and Hall fits). Holes unchanged because the primaries disagree (lab/mobility_caxis/notes/LOG.md) |
| `c_axis_holes_0p9` mobility (opt-in)                | as `c_axis`; holes x 0.90 (thesis drift ratio, primary; Hatakeyama 2003 1/1.15 secondary)                       | Ayalew eq. 3.81 (primary, older) says isotropic: not the default                                                                                                                                                                                                |
| `ishikawa2025` mobility                             | primary (thesis Eqs. 2.16 and 4.6)                                                                              | Full thesis T laws, fitted 300-600 K; holes are HALL mobilities                                                                                                                                                                                                 |
| `legacy_perp` / `ayalew2004` mobility               | primary (Ayalew Table 3.5)                                                                                      | Perpendicular to c. The electron gamma -2.40 is Schaffer's <0001> value; the hole -2.15 is borrowed from electrons. `legacy_perp` keeps mu_min T-independent                                                                                                    |
| `lades2000` high-field                              | electrons primary (Ayalew Table 3.6); parallel-to-c and holes secondary (Nilsson 2000 MC)                       | **Hole beta and T laws are an assumption** (no published fit)                                                                                                                                                                                                   |
| `burin2024` / `ayalew2004` anisotropy               | secondary / primary                                                                                             |                                                                                                                                                                                                                                                                 |
| `lades2000` BGN                                     | primary (Ayalew Table 3.4)                                                                                      | Lindefelt's own four-term form is restated inconsistently by two secondaries, so it is not offered                                                                                                                                                              |

### Numerics

**Breakdown** (`etna.core.avalanche.breakdown_voltage`). The McIntyre ionization integrals are computed with alpha constant on each mesh edge (the edge field) and integrated in closed form per edge. 1 - I_n = exp(-phi(W)) (1 - I_p) holds to round-off. The search is batched over devices: a secant march on ln I_p, then Illinois to 1e-3 V, with Poisson sub-steps of at most 10 V (scaled by T/300 below 300 K).

- `mesh_error_V` is |V_BD(h) - V_BD(h/2)| from a full re-solve on the bisected mesh. Convergence is about first order (ratio 1.5-2.3 over h, h/2, h/4), so the base value is within about 2 x that number of the converged one.
- `quadrature_error_V` is only the 2h-vs-h quadrature part and is much smaller.

**Drift-diffusion with avalanche** (`etna.core.batched.models_dd.ModelContinuityBatch`). Edge impact generation g = alpha_n |Jn| + alpha_p |Jp| is integrated over the edge volume, as devsim's `edge_volume_model`. At fixed potential the system is linear in (n, p) apart from SRH. Each outer iteration therefore runs a full Newton on (n, p) with the 2x2-block tridiagonal Jacobian, then a Poisson update.

A lagged (Gummel) avalanche term converges with ratio 1 - 1/M. On the default dosimeter (`lab/physics_models/bench/models_bench.py`, 2026-10-21, M4):

| M    | lagged sweeps / time | coupled outer x Newton / time |
| ---- | -------------------- | ----------------------------- |
| 2.2  | 15 / 14 ms           | 3 x 2 / 3.5 ms                |
| 34.5 | 136 / 186 ms         | 3 x 2 / 12 ms                 |
| 122  | 421 / 514 ms         | 4 x 2 / 6 ms                  |

Past V_BD the system is singular, and the solve reports `converged=False`. `run_avalanche` marks those points NaN ("at or beyond V_BD").

**Current extraction.** The current is Jn + Jp on the best-conditioned edge, where the Scharfetter-Gummel terms are smallest. The discrete equations conserve Jn + Jp exactly. With every model off, the coupled solver reproduces `ContinuityBatch` to <= 1e-15.

**High-field mobility.** mu(E) on each edge is re-evaluated from the potential at every outer iteration. In devsim the field factor is written into the currents through the chain rule.

**Parity** (`tests/test_physics_models_devsim.py`). Device: p+ / 3 um 1e17 / n+, with uniform generation 1e18 cm^-3 s^-1 in the epi, at 300 V (M ~ 34).

- Engine vs devsim: 3e-12 relative, avalanche and high-field together.
- The devsim run takes 2.0 s (build, ramp and solve) against 13 ms for the engine, ramp included.
- DD gain vs the generation-weighted McIntyre M: below 5e-6 relative in (1 - 1/M) up to M = 122.

### Validation and calibration labels

**Breakdown vs measurement.** Hatakeyama 2004 measured V_BD = 1940 (1e16/N_D)^0.8 V on (0001) p+n mesas (10 um epi, N_D 3e16-2e17). The table compares a modelled p+ / 10 um n / n+ diode. These coefficients were fitted on these same diodes, so for `hatakeyama2004` this is **consistency, not validation**. The other sets are independent comparisons and show the model spread:

| N_D (cm^-3)       | 3e16    | 5e16    | 1e17    | 2e17    |
| ----------------- | ------- | ------- | ------- | ------- |
| measured fit (V)  | 805.6   | 535.3   | 307.5   | 176.6   |
| hatakeyama2004    | -0.6 %  | -0.3 %  | +1.8 %  | +7.0 %  |
| niwa2015          | -6.4 %  | -2.7 %  | +5.2 %  | +18.1 % |
| zhao2019          | -10.8 % | -8.1 %  | -2.1 %  | +8.0 %  |
| loh2008           | -17.8 % | -14.6 % | -8.2 %  | +2.0 %  |
| stefanakis2021    | -2.1 %  | +4.3 %  | +16.2 % | +33.5 % |
| ayalew2004_legacy | +47.8 % | +55.7 % | +70.6 % | +93.9 % |

The peak field at breakdown, against Hatakeyama's critical field F_C = 2.70 (1e16/N_D)^-0.1 MV/cm (SISPAD 2005), is between -1.6 and -3.2 % for `hatakeyama2004`.

**Default dosimeter V_BD**: a 1D parallel-plane prediction and an **upper bound**.

- The values are 2225.8 V (hatakeyama2004), 1955.4 V (niwa2015), 1928.2 V (zhao2019), 1825.6 V (loh2008) and 1916.9 V (stefanakis2021). Peak field 1.9-2.3 MV/cm, mesh error at most 0.11 V.
- The spread between sets (about 20 %) is the real uncertainty.
- The real device breaks down at its edge termination, which is not modelled.
- The epi is low doped and in punch-through at about 3 V. That extrapolates beyond the measured N_D range of the sets: their field range still covers 1.9-2.3 MV/cm, but the doping regime differs.

**Temperature.** `niwa2015`, `zhao2019` and `stefanakis2021` give a positive temperature coefficient of V_BD from 300 to 448 K on the default dosimeter: 1955 -> 2038 V (+4.2 %), 1928 -> 1999 V (+3.7 %) and 1917 -> 1961 V (+2.3 %).

### Effects on the outputs (default dosimeter, model vs default physics)

| Output                                                      | Effect                                                                                                                                                                    | Why                                                                                                                                                                 |
| ----------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| C-V                                                         | none                                                                                                                                                                      | Pure electrostatics; none of these models enters Poisson                                                                                                            |
| CCE (alpha, steady state) with `lades2000`                  | -9.1e-5 at -10 V to -1.1e-4 at -100 V (c-axis default; -5.0e-5 / -6.5e-5 on a `legacy_perp` device)                                                                       | Slower carriers recombine slightly more; lifetime >> transit                                                                                                        |
| `legacy_perp` device + `burin2024` anisotropy + `lades2000` | -1.7e-4 / -1.2e-4                                                                                                                                                         | The previous opt-in route to c-axis transport (an anisotropy is refused on the c-axis default)                                                                      |
| Electron transit time (full depleted width), `lades2000`    | x1.44 at -10 V, x7.10 at -100 V (66.4 ps instead of 9.4 ps; parallel v_sat 1.64e7 since the default is c-axis). `legacy_perp`: x1.25 / x4.65 (52.6 ps instead of 11.3 ps) | Velocity saturation is strong at 1e5 V/cm. Holes x1.10 to x2.54 (1028 / 215 ps, parallel hole v_sat 6.5e6). This is what a pulse-shape / timing analysis should use |
| CCE with avalanche, `hatakeyama2004`                        | 1.00001 at -1050 V, 1.0145 at -1525 V, 2.12 at -2000 V                                                                                                                    | APD-like gain only within ~25 % of V_BD; none at the operating -10 to -100 V (M - 1 < 1e-5)                                                                         |
| Dark current with avalanche                                 | x1.0137 at -1525 V, x2.03 at -2000 V                                                                                                                                      | Multiplication of the TAT effective generation; see the note below                                                                                                  |

**The dark current with avalanche** uses the TAT effective generation and the Z1/2 SRH level **without** the Hurkx (1 + Gamma) factor, whose closed form is invalid above about 1.1e6 V/cm. Its absolute value is therefore below `run_dark_current`'s at the same bias. Read it as a multiplication factor (`metadata["gain"]`).

### Facades

- **`run_breakdown(config, models=None, *, n_points=30, M_max=1e4, v_max=1e4, v_tol=1e-3, mesh_check=True)`**
  - x: bias on a geometric schedule toward V_BD.
  - y: M_p.
  - metadata: `V_BD`, `V_BD_refined_mesh`, `mesh_error_V`, `quadrature_error_V`, `E_max_at_breakdown`, `field_extrapolated`, `M_n`, `M_p`, `I_n`, `I_p`, `E_max`, `note`, provenance.
- **`run_avalanche(config, models=None, v_start=-100, v_stop=-1000, n_points=10, *, source="alpha"|"dark", approach_breakdown=False, M_max=1e3, generation_rate=1e18)`**
  - y: unclipped CCE (alpha) or dark current (A).
  - metadata: `J`, `J_baseline`, `gain`, `M_mcintyre`, `V_BD`, `converged`, `reasons`.
  - `approach_breakdown=True` replaces the linear grid with a geometric schedule toward V_BD: the adaptive near-breakdown sweep. Its last point is V_BD - (V_BD - |v_start|) / `M_max`. `M_max` is a bias depth, not a target multiplication: 1e3 ends near M ~ 400 on the default device.
- **`run_transit_time(config, models=None, v_start=-5, v_stop=-100, n_points=20)`**
  - Drift-only, full-width crossing times of the depleted epi (field > 1e3 V/cm).
  - Metadata: `t_e`, `t_h`, `t_e_low_field`, `t_h_low_field`, `W_dep`.

All three are 1D only. They refuse 2D configs, `radiation_traps`, and `bandgap_narrowing`.

**devsim twin.** `etna.core.devsim_models.apply_physics_models(device_info, models)` installs the same models on a 1D devsim DD device. Call it after `add_generation_to_dd`: any later equation re-registration drops the impact term, and `total_current` then raises.

### Not done (open in issue #6)

- **Bandgap narrowing** is registered with provenance but not coupled into Poisson. Contact potentials and n_i(x) touch the shared Poisson code. The expected effect is about -0.04 V on V_bi at the 1e19 p+, below the C-V fit's resolution.
- **2D**: the anisotropic mobility and anisotropic impact ionization (Hatakeyama interpolation) are not coupled into engine2d or devsim 2D.
- **Field-dependent trapping** beyond the Hurkx TAT factor.
- **Transient runs with the models.** `run_transient` and FLASH do not take models.

---

## 19. Single-ion track engine (opt-in, `etna.core.track`)

The track engine solves one ion track in the default stack (or any `DeviceConfig`), axisymmetric in (r, z) at normal incidence. It models the dense track core, plasma (columnar) recombination by SRH + Auger, the pulse-height defect (PHD) and, opt-in, the lattice temperature.

It is a separate, opt-in engine with its own facades (`run_track_event`, `run_track_cce_let`, `run_tilted_alpha_track`, section 19.2). No other `run_*` facade calls it, and `build_cce_let_table` (2D Cartesian sheet track, linear) is unchanged. Issue #9; numbers and commands are in `lab/track_plasma/notes/LOG.md`.

```python
from etna import DeviceConfig
from etna.core.track.phd import MESH_FAST, cce_vs_let, simulate_track, ion_let_profile, whk_along
from etna.core.track.radial import gaussian, whk, thermalized, beta_of

r = simulate_track(DeviceConfig(), LET_keV_um=1000.0, V_bias=50.0,
                   radial=thermalized(whk(beta_of(1.0)), 10e-7), mesh=MESH_FAST)
print(r.cce, r.recombined_fraction, r.ledger_residual)   # 0.998687 0.00131 2e-9

df = cce_vs_let(DeviceConfig(), (1, 100, 1000), 10.0, gaussian(50e-7), mesh=MESH_FAST)

let, epu = ion_let_profile("alpha", 5.486)               # NIST Bragg profile in SiC
r = simulate_track(DeviceConfig(), LET_keV_um=let, V_bias=50.0, region="all",
                   radial=lambda z: thermalized(whk_along(epu)(z), 10e-7), mesh=MESH_FAST)
```

**Model**

- **Discretisation.** Cylindrical finite volumes; the first radial cell is 0.5-2 nm and spacing grows geometrically to R = 20-30 µm. The equations are fully coupled Newton on (ψ, n, p) with Scharfetter-Gummel fluxes, devsim's own discretisation. Both contact planes are ohmic over the full area, and the outer cylinder is reflecting.
- **Time integration.** TR-BDF2 with an embedded local-truncation-error estimate, on a charge-weighted error per step (`rtol`, default 1e-4). Steps run from 1 fs to µs. The modified (chord) Newton reuses one SuperLU factorisation across stages and steps.
- **Dark state.** Identical to the batched 1D engine (ψ to 1e-14 V; steady CCE to 1e-11).
- **Generation.** `radial` gives the fraction of pairs per annulus exactly, so a core finer than the mesh is a conserved sub-grid quantity. Available models:
  - `gaussian(w)`: commercial-TCAD style.
  - `core_penumbra`.
  - `whk(beta)`: Waligórski-Hamm-Katz, transferred from water by mass thickness (an approximation).
  - `thermalized(model, s)`: a planar Gaussian spread s.
- **Instantaneous or pulsed deposition.** Deposition is instantaneous by default; `pulse_sigma_s` gives a Gaussian pulse instead.
- **Opt-in physics** (`etna.core.track.grid.device_from_config`):
  - `auger_model="injection"`: Auger C falls above ~5e18 cm⁻³ (Tanaka 2023 via arXiv:2410.06798).
  - `carrier_scattering=True`: Conwell-Weisskopf, Lades 4H-SiC fit.
- **CCE.** `TrackResult.cce` is the induced charge from the pair balance: 1 − (recombined + wrong-contact minority + stored)/N_gen. The direct majority-carrier contact sums (`cce_cathode`) carry about 1e-5 of float64 round-off in the 1e19 layers.

**Checks** (`tests/test_track_engine.py`):

- Each run reports `ledger_residual` (pairs) and `extra["quad_err"]`.
- With `energy=True`, the discrete Poynting theorem ∫J·E = Σ_c ψ_c ∫I_c − ΔW_field is checked; it closes to 3.5e-5 at rtol 1e-4 and falls with the tolerance (1.8e-6 at 1e-6), so the residual is time-discretisation error.
- The heat solver matches the analytic line source to within 3 %.

**Results to know** (50 V unless stated, details in the LOG):

- **The core radius decides the plasma loss.** At LET 1000 keV/µm the recombined fraction is:
  - 1.3e-4 with a 50 nm core;
  - 1.3e-3 for WHK spread by 10 nm;
  - about 5e-2 for 1 nm;
  - 0.18-0.22 for the bare sub-nm WHK core. That value is not mesh-converged: the first-cell density exceeds the atomic density, and drift-diffusion is outside its validity there.

  These are `MESH_FAST` numbers: the PHD carries 5-8 % mesh error, and the sub-nm cores depend on the axis cell.

  The thermalisation radius for 4H-SiC is not known. Every PHD number from this engine is conditional on `radial`.

- **Normal-incidence ions at 10–50 V.** Protons and 5.5 / 8 MeV alphas recombine 0.1–0.3 % of their pairs. The p+ dead layer (wrong-contact loss) is larger, at 2–6 %.
- **0 V.** The CCE _rises_ with LET (0.872 → 0.933 from 0.1 to 3000 keV/µm), so the linear-table assumption fails without depletion.
- **Steady-state CCE is not the single-particle CCE** even at low injection. SRH in the depletion region needs both carriers at once (n_i ≈ 0).
- **The 241Am tilted-alpha data (Fig. 4b) are not reproduced** with the library lifetimes. This holds for the normal-incidence proxy and for the proper tilted (2.5D) model of section 19.2, and it confirms `lab/microdosimeter/notes/LOG.md` §9. See 19.3 for what the data do constrain.

**Electrothermal (opt-in, `thermal={"couple": True, "sink": "none"|"cathode"}`)**

- **Physics.** Anisotropic κ (4.71 W/cm/K radial, 3.24 W/cm/K along c; arXiv:1712.00830) and c_v = 0.69 J/g/K × 3.21 g/cm³. Heat sources are Joule J·E (SG edge power), E_g per recombination, and prompt (E_pair − E_g) heat. Mobility and SRH τ follow T; transport stays isothermal SG (no thermodiffusion).
- **Energy ledger.** The thermal ledger closes to ≤ 1e-7, adiabatic or with the cathode heat sink. The sink heat is the integrated conductive flux, so the ledger is an independent check.
- **Result.** The effect of heating depends on the core size:
  - 50 nm core: negligible (1.8 K).
  - 10 nm core at LET 1000: ΔT_max 118 K for picoseconds and ΔCCE −2.5e-4, about 19 % of the PHD (stable in time step).
  - Smaller cores: the result is not meaningful.

**Cost and limits**

- **Cost.** 2–5 s per event on `MESH_FAST` (6930 unknowns) and 4–15 s on `MESH_DEFAULT`, dominated by sparse LU. Use DEFAULT or finer for PHD numbers: the mesh error is about 5–8 % of the PHD.
- **Not modelled:**
  - a full 3D tilted track (section 19.2 composes 2D transverse slices instead);
  - self-consistent velocity saturation (only a frozen-field bracket exists, `mu_field`);
  - band-gap narrowing;
  - hot carriers;
  - nuclear-stopping or window defects.

### 19.1 Excess unknowns (the float64 floor, default since the issue #9 follow-up)

`run_transient` now solves for the CHANGE from the dark state, (dψ, dn, dp), by default (`kernels.assemble_excess`). Every term is the exact difference from the dark state, evaluated without cancellation:

- SG flux: a[dn_k B(D) − dn_l B(−D)] + a[n_dk ΔB(D0, δ) − n_dl ΔB(−D0, −δ)], with ΔB(x, d) = B(x + d) − B(x) from `bern_diff` (series / expm1 forms; 1e-14 relative against a 50-digit reference, tested);
- recombination: (n_d dp + p_d dn + dn dp) g(n, p) + (n_d p_d − n_i²)(g − g_d);
- Poisson: V (dp − dn) + the Laplacian of dψ.

In absolute unknowns the excess in a 1e19 cm⁻³ layer is resolved only to ~eps × 1e19 per node and stage: the `ledger_floor` was ~1e-6 of N_gen at LET 1 keV/µm. In excess unknowns it is ~n_steps × 4e-16 (MESH_FAST, default stack, 0 V, LET 1: ledger 5.1e-6 → 2.5e-8, floor 9.6e-7 → 9e-14). The CCE agrees with the absolute form to 6e-9 (LET 1000) and 8e-7 (LET 1, where the absolute form was floor-limited).

`excess=False` restores the absolute form (first release). The thermal coupling and `energy=True` use absolute states, so `excess=None` (default) switches to the absolute form for them; `excess=True` with either raises.

### 19.2 Tilted tracks: transverse slices and the 2.5D composition

A tilted track is not axisymmetric. Two new modules handle it:

- **`etna.core.track.cross_section`**: the plane perpendicular to a long, locally uniform track, with the transverse field E⊥ = E sin θ imposed far away. It is a 2D Cartesian half plane (mirror plane through the track and the field), with unit length along the track and absorbing Dirichlet edges at |u| = R, v = R. The physics and numerics are the engine's (`TrackSystem(dev, geometry=...)` now accepts any finite-volume mesh). Backgrounds:
  - `"depleted"`: n = p = n_i, linear ψ, exact discrete steady state (tested); `N_space` adds an in-plane space charge for the field gradient;
  - `"neutral"`: n = N_D, E⊥ = 0.

  `slice_loss(E_perp, LET, radial, ...)` returns the recombined fraction of the slice. `mu_field=True` is a frozen-far-field velocity-saturation bracket (Lades 2000 / Nilsson parallel-to-c v_sat).

- **`etna.core.track.tilted`**: `tilted_alpha_cce(config, V_values, tau_scale=..., s_nm=...)`. CCE(V) = CCE_column(V) − Σ_k w_k f_ex(s_k, E⊥(z_k; V)), where:
  - CCE_column is the single-particle linear response of the tilted generation on the 1D column (track engine, M = 1);
  - f_ex = f(LET(s)) − f(LET → 0) is the plasma excess of the slice, tabulated in a `PlasmaTable` (6 path nodes × 11 fields, ~77 slices) and interpolated.

  Approximations, reported in `notes`:
  - the parallel component E cos θ is dropped (it only adds erosion, so the slice overestimates the loss);
  - the box R matters at low field, where the plasma grows to several µm before erosion;
  - undepleted path segments use the E → 0 depleted-slice excess.

Facades (`etna.api.track`, exported from `etna`):

| function                                                                                                        | returns                                                                 |
| --------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------- |
| `run_track_event(config, LET_keV_um, V_bias, radial="whk", s_nm=10, mesh="fast")`                               | x = t (s), y = cathode current (e/s); metadata cce, recombined, ledgers |
| `run_track_cce_let(config, LET_values, V_bias, radial="gaussian", s_nm=50)`                                     | x = LET, y = CCE; PHD, recombined                                       |
| `run_tilted_alpha_track(config=microdosimeter, V_values, tau_scale=1, s_nm=10, R_um=10, workers=1, table=None)` | x = −V, y = CCE; cce_column, plasma_loss, shape, table, notes           |

### 19.3 Fig. 4b with the tilted model: a calibration, not a reproduction

Details are in `lab/track_plasma/notes/LOG.md` §§10-12. The bench is `lab/track_plasma/bench/tilted_bench.py`. Device: the microdosimeter preset with N_D 1e14.

**With the library lifetimes (τ_n 1 ns, τ_p 600 ns) the tilted track does not reproduce Fig. 4b.**
- The plasma loss changes by at most 1e-3 between the fields of −20 and −50 V; the data rise by 1e-2 there.
- The shape misses by rms 0.033, the same as the normal-incidence proxy (0.034). The tilt was not the missing ingredient.

**What moves the bias-dependent loss.**
- The bias-dependent loss is SRH recombination during the plasma time, which scales as about 1/E⊥ and 1/(τ_n + τ_p).
- The following move the loss nearly flat in bias, so the free overall scale absorbs them:
  - the core radius (thermalization length) s;
  - the Auger model;
  - carrier-carrier scattering.
- Velocity saturation flattens the field dependence.

**Calibration** (free scale, the −20 to −50 V tail):

| quantity | value |
|---|---|
| τ_scale | 0.039, 68 % interval 0.035-0.047 at σ 1e-3 |
| τ_n + τ_p ≈ τ_p | ≈ 23 ns |
| systematic band (box R, field gradient, v_sat, N_D) | 17-42 ns |
| χ² | 6.4 for 2 dof at σ 1e-3; 25 at the digitisation σ of 5e-4 |

This lifetime is about 100x below good-epi ambipolar lifetimes (~2 µs, arXiv:1011.1039, abstract).

**What the fit cannot do.**
- No single τ fits all seven points (χ² 371 for 5 dof).
- The 0 to −10 V points depend on path elements below 1e4 V/cm, outside the validity of the 2.5D model.

**What the data constrain.**
- The singular values of the scale-projected Jacobian show ONE constrained direction: the plasma SRH lifetime.
- **The core radius is not constrained:** Δχ² < 0.2 between s = 1 and 50 nm.
- Hence no preset: `tau_scale` is a user knob of `run_tilted_alpha_track`.
