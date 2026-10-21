# ETNA — Epitaxial Tcad Numerical Analysis (for 4H-SiC detectors)

ETNA is a Python TCAD library with a Streamlit UI for **4H-SiC p⁺/n⁻/n⁺ radiation
detectors**: the INFN-LNS Petringa group's 4 mm² proton dosimeter and its 25 mm²
microdosimeter sibling. It solves Poisson plus electron/hole drift-diffusion in 1D and 2D
and derives C-V, depletion and field profiles, charge-collection efficiency (CCE), a
dark-current budget, CCE vs proton fluence, temperature and FLASH dose-rate sweeps, and
microdosimetric spectra. It is built on the open-source [devsim](https://devsim.org/)
solver plus two numba engines (1D batched, 2D) that reproduce devsim inside checked
parameter boxes.

ETNA is a research code, not a commercial TCAD package. Its numerics are verified
against devsim; its physics is **calibrated on one measured C-V curve and one leakage
bound**, and most other outputs are model predictions not yet checked against a
measurement of this device (see [What to trust](#what-to-trust)).

The full guide is **[docs/USAGE.md](docs/USAGE.md)**: every facade with a runnable
example, presets, units, solver provenance, validity limits, the app, exports,
notebooks, tests and troubleshooting. This README is the overview.

---

## Quick start

Requires Python ≥ 3.13 and [uv](https://docs.astral.sh/uv/) (the only supported
installer; do not use pip/venv). devsim ships as a wheel; on Linux install
OpenBLAS/LAPACK first (`apt-get install libopenblas-dev liblapack-dev`).

```bash
git clone https://github.com/ngcex/sic-tcad-detector.git
cd sic-tcad-detector
uv sync --extra dev                 # library + pytest, xdist, jupyter, nbstripout, ruff, meshio
# optional: uv sync --extra dev --extra viz   (manim animations + gdstk GDSII export;
#           on macOS run `brew install pkgconf cairo pango` first)

uv run python examples/cv_example.py                                 # C-V of the default device
uv run streamlit run app/main.py --server.address localhost         # the interactive app
uv run pytest -n auto -q                                             # the test suite
```

Always sync with `--extra dev`: a plain `uv sync` removes the dev tools from `.venv`, and
`uv run pytest` may then run a system pytest with the wrong interpreter. A later
`uv sync --extra dev` without `--extra viz` also removes manim and gdstk again.

Without a checkout (library only; ETNA is distributed from GitHub, not PyPI):

```bash
uv pip install "git+https://github.com/ngcex/sic-tcad-detector"               # into the active venv
uv pip install "etna[viz] @ git+https://github.com/ngcex/sic-tcad-detector"   # + gdstk / manim
# or from a checkout: uv build, then uv pip install dist/etna-*.whl
```

The wheel carries the runtime data (NIST tables, SR-NIEL CSVs, validation points, layout
defaults, measurement templates) as `etna/_data/`. Raw provenance files (scraped HTML)
are only in the sdist and the repository. The Streamlit app, notebooks, scripts and tests
are not in the wheel: run them from a checkout. Details:
[USAGE §3](docs/USAGE.md#3-installation).

Python API (there is no command-line tool; `etna` is an importable package):

```python
from etna import DeviceConfig, run_cv, run_cce, run_all

cfg = DeviceConfig()                       # default 4 mm2 dosimeter (measured-C-V stack)
cv = run_cv(cfg, 0.0, -30.0, 31)           # SimResult: x = bias (V), y = C (F, for area_cm2)
print(cv.metadata["solver"], cv.metadata["reach_through_V"])        # batched 3.1
cce = run_cce(cfg, v_start=-1.0, v_stop=-10.0, n_points=10)        # y = CCE in [0, 1]
report = run_all(cfg)                      # every section in one call
for row in report.rows():
    print(row["section"], row["status"], row["solver"])
```

Reverse bias is negative. Every result records which solver produced it
(`metadata["solver"]`: `batched`, `engine2d` or `devsim`) and whether the device lies in
the engine/devsim cross-check box (`metadata["validated_domain"]`). The public API is
listed in [`etna/__init__.py`](etna/__init__.py); each `run_*` facade is documented in
[USAGE §6](docs/USAGE.md#6-the-run_-facades).

---

## What to trust

| Output                                                           | Status                            | Basis                                                                                                                                                                                                                                                            |
| ---------------------------------------------------------------- | --------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Engine results (C, W, CCE, field, dark current, FLASH) vs devsim | **verified**, code to code        | 1D dC/C ≲ 3e-8, dCCE ≲ 1e-9 to −100 V; 2D potential ~1e-13 V. Says nothing about the real device                                                                                                                                                                 |
| C-V of the default device                                        | **calibration** (fit)             | the doping profile is fitted to the digitised Photons Fig. 4b C-V (plateau rms 0.013 µm, knee rms 0.24 µm). Agreement is fit quality, not validation                                                                                                             |
| Doping profile N_D(x)                                            | **inferred**                      | inversion of one C-V curve, no SIMS; bulk 3.4e13 cm⁻³ is below the documented 0.5–1e14                                                                                                                                                                           |
| CCE vs bias                                                      | **model**                         | low injection, steady state; no measured CCE of the dosimeter reproduced. The microdosimeter's measured alpha CCE-vs-bias shape is **not** reproduced                                                                                                            |
| Dark current                                                     | **calibration** (budget estimate) | one trap rate N_t fitted to < 18 pA at −60 V on 4 mm²; bias and temperature shapes are model assumptions                                                                                                                                                         |
| CCE vs proton fluence                                            | **relative shape only**           | SR-NIEL proton hardness factors (Si 1-MeV-neq convention). The **SiC neutron NIEL normalisation is open**, so the absolute SiC fluence scale is unvalidated. Points at or beyond the compensation fluence (3.48e12 p/cm² at 62 MeV on the default stack) are NaN |
| FLASH dose-rate CCE                                              | **numerical sensitivity bound**   | SRH + Auger only; high-injection plasma recombination is not implemented                                                                                                                                                                                         |
| Tissue-equivalence κ(E)                                          | **data + approximation**          | real NIST PSTAR/ASTAR tables; SiC by Bragg additivity of Si + C. κ is the SiC → water LET ratio, 0.43 at 1 MeV → 0.375 at 150 MeV for protons                                                                                                                    |
| Reference layout (GDSII/SVG)                                     | **REFERENCE-NOT-FAB**             | drawn from `data/foundry/layout_params.toml`; guard rings, pad, dicing street, oxide and the microdosimeter die size are PLACEHOLDER values until mask data is entered                                                                                           |

The only literature check (`scripts/validation/cv_vs_bruzzi2023.py`, Bruzzi & Verroi 2023) uses a Schottky diode on a different 30 µm wafer: depletion width +3.1 % at −40 V, +1.0 %
with V_bi matched (the paper quotes W to 2 significant figures). It checks the solver
scale, not this device.

**Opt-in models** (default `False`, so every baseline stays as pinned):

- `DeviceConfig(incomplete_ionization=True)`: carrier-dependent N/Al ionization inside
  Poisson (Burin et al.). The default keeps a fixed ionized fraction for the p⁺ Al. Not
  the default because the measured-C-V presets were fitted without it
  ([USAGE §5.3](docs/USAGE.md#53-incomplete-ionization-opt-in-incomplete_ionizationtrue)).
- `DeviceConfig(radiation_traps=True)`: radiation-induced deep traps (Z1/2, EH4, EH6/7)
  as charged species in Poisson, so `run_radiation_damage` can solve past the
  compensation fluence (also a checkbox on the app's Radiation Damage page).
  Predictions there rest on a stated occupancy assumption
  ([USAGE §6](docs/USAGE.md#charged-deep-traps-beyond-compensation-opt-in-radiation_trapstrue)).

- `etna.PhysicsModels` with `run_breakdown` / `run_avalanche` / `run_transit_time`: a
  model library with provenance, 1D only (issue #6). It covers velocity saturation, c-axis
  anisotropy, and impact ionization (Hatakeyama 2004 and four other sets) with avalanche
  breakdown and gain. The default facades do not use it
  ([USAGE §18](docs/USAGE.md#18-physical-model-library-opt-in-high-field-mobility-anisotropy-avalanche-and-breakdown)).

**Not modelled:** high injection / plasma, edge termination and guard rings on the
default path (the opt-in breakdown is the parallel-plane value), oxide charge, 3D. Full list and the
full-depletion definitions in [USAGE §8](docs/USAGE.md#8-validity-and-limits).

---

## Ways to use it

| Interface                       | Entry point                                                                                                                                                                                                                                            | Details                                                               |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------- |
| Streamlit app                   | `uv run streamlit run app/main.py --server.address localhost`                                                                                                                                                                                          | [USAGE §9](docs/USAGE.md#9-the-streamlit-app)                         |
| Python API                      | `import etna` (`run_cv`, `run_cce`, `run_field`, `run_radiation_damage`, `run_dark_current`, `run_temperature_sweep`, `run_flash_recombination`, `run_transient`, `run_microdosimetry`, `run_cce_tilted_alpha`, batches, `ParametricSweep`, `run_all`) | [USAGE §4–7](docs/USAGE.md#4-quick-start)                             |
| Presets                         | `DeviceConfig()` (= `etna.measured_cv_config()`), `measured_cv_tabulated_config()`, `nplus_calibrated_config()`, `legacy_config()`, `microdosimeter_config()`                                                                                          | [USAGE §5.2](docs/USAGE.md#52-presets)                                |
| Publication figures, animations | `etna.viz` (`figure_for`, `export_figure`, `animate`)                                                                                                                                                                                                  | [USAGE §10](docs/USAGE.md#publication-figures-and-animations-etnaviz) |
| Engineering exports             | `etna.export` (spec sheet JSON/CSV/MD/PDF, cross-section, mesh CSV/VTK/VTU, reference GDSII/SVG)                                                                                                                                                       | [USAGE §10](docs/USAGE.md#10-exports-python)                          |
| Measurement inputs              | `etna.measurements` (C-V, SIMS/SRP and I(T) CSV loaders, model comparison, C-V inversion)                                                                                                                                                              | [USAGE §17](docs/USAGE.md#17-measurement-inputs-etnameasurements)     |
| Notebooks                       | 21 notebooks in `notebooks/`, committed without outputs                                                                                                                                                                                                | [USAGE §12](docs/USAGE.md#12-notebooks)                               |

**The app** (12 pages): Home, Simulate → Figures, Run All, C-V Analysis, Charge
Collection (CCE), Field Map, Radiation Damage, Dark Current, Microdosimetry, Batch Sweep,
Export, Measurements. The device is set once in the sidebar (geometry, doping,
temperature, area, presets) and shared by every page. Each result shows its solver and
validity notes, and its CSV download carries the same provenance in a `#` header. The
app has no authentication and binds to localhost (`.streamlit/config.toml`); for shared
use put it behind an authenticating reverse proxy.

**Notebooks:** run them from `notebooks/` so their figures land in the gitignored
`notebooks/figures/`:

```bash
cd notebooks && uv run jupyter lab 01_phase1_validation.ipynb
```

Notebooks 12 (multi-defect comparison), 14 (validation vs literature trends), 19
(alternative structures) and 20 (feasibility report) have no app page; 08 (FLASH
transient) and 16 (CCE vs LET) are covered in part by Simulate → Figures.

---

## Where outputs go

| Producer                                                | Output location                                                             |
| ------------------------------------------------------- | --------------------------------------------------------------------------- |
| Python API                                              | in memory (`SimResult`, `Report`); nothing is written unless you save it    |
| `etna.core.plotting.save_figure` (notebooks)            | `$ETNA_FIGURE_DIR`, else `./figures` of the working directory (gitignored)  |
| `etna.viz.save_figure(fig, path)` / `export_figure`     | the path you give / bytes                                                   |
| `etna.export` builders                                  | bytes or strings; you write them                                            |
| Streamlit app                                           | browser downloads (CSV, JSON, PNG/SVG/PDF, zip, MP4/GIF, VTK/VTU, GDSII)    |
| Notebook 16 table rebuild (`ETNA_REBUILD_CCE_TABLES=1`) | `$ETNA_FIGURE_DIR/data/`, never over `data/`                                |
| Notebook 17                                             | rewrites the tracked `data/synthetic_mc_events.csv` (seeded, bit-identical) |

Generated figures are not tracked: `figures/` is gitignored except the Bruzzi validation
plot `figures/cv_validation_bruzzi2023.png`.

---

## Reproducing the validation numbers

```bash
uv run python scripts/validation/cv_vs_bruzzi2023.py     # rewrites figures/cv_validation_bruzzi2023.png
uv run pytest -n auto -q tests/test_batched_engine.py tests/test_api_engine_switch.py  # engine vs devsim parity
uv run pytest -q tests/test_measured_cv_preset.py tests/test_viz_figures.py          # C-V fit and validation variants
```

- Engine/devsim parity: `tests/test_batched_engine.py`, `tests/test_engine2d*.py`,
  `tests/test_dark_current_engine.py`, `tests/test_flash_engine.py`,
  `tests/test_fluence_engine.py` (list in [USAGE §14](docs/USAGE.md#14-tests)).
- C-V calibration: the digitised curve and fit live in `lab/cv_gap/` (`notes/LOG.md`);
  `etna.viz.validation.w_creep()` and `tilted_alpha_shape()` draw the model against the
  digitised points in `data/validation/`
  ([USAGE §8](docs/USAGE.md#8-validity-and-limits)).
- Every calibrated constant cites the `lab/*/notes/LOG.md` with the commands and numbers
  behind it ([`lab/README.md`](lab/README.md)).

---

## Foundry layout parameters

Every lateral and passivation dimension of the reference layout, the spec sheet and the
cross-section comes from **[`data/foundry/layout_params.toml`](data/foundry/layout_params.toml)**
(one table per device and dimension: `value`, `unit`, `class`, `source`, `uncertainty`).
To use the mask data, edit that file (or a copy) and set `class = "DOCUMENTED"` with a
`source` naming the mask set. Point ETNA at a copy with

- `ETNA_LAYOUT_PARAMS=/path/to/layout_params.toml`, or
- `etna.export.load_layout_params(path)` / the `layout_params=` argument of the export
  builders, or
- the upload box on the app's Export page.

In an installed wheel the default file is read-only package data. Copy it out first with
`cp "$(python -c 'import etna.export as e; print(e.LAYOUT_PARAMS_PATH)')" my_layout.toml`,
then set `ETNA_LAYOUT_PARAMS=$PWD/my_layout.toml`.

The file is validated on load (`LayoutParamsError` names the field). While a device has
any `PLACEHOLDER` dimension, its GDSII top cell ends in `REFERENCE-NOT-FAB`, the exports
list the placeholder names and raise a `PlaceholderDimensionWarning`. Details:
[USAGE §10](docs/USAGE.md#layout--mask-dimension-file).

`deliverables/` holds the two foundry-facing physics specifications (process-agnostic
feasibility specs, not mask layouts): **DESIGN-1** dosimetry p-n diode and **DESIGN-2**
microdosimeter sensitive-volume geometry, as Markdown (the source) and PDF.
`deliverables/DIMENSIONS-provenance.md` classes every layout dimension as DOCUMENTED,
MEASURED-FROM-FIGURE or PLACEHOLDER.

---

## Environment variables

| Variable                          | Effect                                                                                                                                      |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `ETNA_FIGURE_DIR`                 | where notebook figures (`etna.core.plotting.save_figure`) and rebuilt CCE(LET) tables go; default `./figures` of the working directory      |
| `ETNA_LAYOUT_PARAMS`              | alternative layout parameter TOML for `etna.export`                                                                                         |
| `ETNA_DATA_DIR`                   | read runtime data from this directory (laid out like the repository `data/`) instead of the packaged `etna/_data` or the checkout's `data/` |
| `ETNA_DEVSIM_VERBOSE=1`           | show devsim's per-solve log (silenced by default)                                                                                           |
| `ETNA_DEVSIM_STOCK_SOLVER=1`      | keep devsim's stock UMFPACK callback instead of ETNA's cached one (results are bit-identical either way)                                    |
| `ETNA_REBUILD_CCE_TABLES=1`       | notebook 16 rebuilds the CCE(LET) tables (~40 min) instead of loading `data/cce_let_table_*.json`                                           |
| `ETNA_CCE_TABLE_JOBS`             | worker processes for that rebuild (default min(8, CPU count))                                                                               |
| `ETNA_REPORT_DEVSIM_LEAKS=<file>` | tests: list the tests that leave devsim devices behind                                                                                      |

---

## Tests

```bash
uv run pytest -n auto -q                   # default run: the local gate
uv run pytest -n auto -q -m exhaustive     # the 6 deselected deep variants
uv run pytest -n0 -x tests/test_api_cce.py # serial, e.g. with --pdb
```

1653 tests in 128 modules. `pytest.ini` deselects the 6 `exhaustive` ones, so the default
run is 1647 tests: all pass with `--extra dev --extra viz`, and with `--extra dev` only
1639 pass and 8 skip (manim/gdstk). On a 10-core Apple M4 shared with other jobs, `-n auto`
took 169–205 s on 2026-10-21 (600 s on a heavily loaded machine): plan on about 3 minutes. `pytest.ini` supplies the tuned
xdist scheduler; use `-n0` for pdb.

CI (`.github/workflows/ci.yml`) runs on pushes to `main`, pull requests and manual
dispatch: lock check, ruff, notebook hygiene, the generator drift check and the default
suite with `-n 4`. A nightly job runs the exhaustive tests and `pip-audit`. Scheduling,
the devsim leak guard and test families: [USAGE §14](docs/USAGE.md#14-tests).

---

## Repository layout

| Path                     | Contents                                                                                                                                                                                                                  |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `etna/`                  | the package: `api/` (public facades, `DeviceConfig`), `core/` (solvers and physics: devsim builders, `batched/` 1D engine, `engine2d/`, damage, dark current, microdosimetry), `export/`, `viz/`, `measurements/`         |
| `app/`                   | Streamlit app: `main.py`, `workflows/` (one module per page), `components/`                                                                                                                                               |
| `notebooks/`             | 21 analysis notebooks, without outputs                                                                                                                                                                                    |
| `tests/`                 | pytest suite (API, app `AppTest`, engine/devsim parity, presets, core physics)                                                                                                                                            |
| `docs/USAGE.md`          | the usage guide                                                                                                                                                                                                           |
| `data/`                  | NIST stopping powers (`nist_star/`), SR-NIEL tables (`niel/`), layout parameters (`foundry/`), measurement templates (`measurements/`), digitised validation points (`validation/`), CCE(LET) tables, synthetic MC events |
| `deliverables/`          | DESIGN-1 / DESIGN-2 foundry specs and dimension provenance                                                                                                                                                                |
| `lab/`                   | research spikes; each `notes/LOG.md` records the commands and numbers behind a calibrated default                                                                                                                         |
| `scripts/`               | notebook generators and drift check, the Bruzzi validation script ([`scripts/README.md`](scripts/README.md))                                                                                                              |
| `examples/`              | `cv_example.py`                                                                                                                                                                                                           |
| `.issues/`, `.planning/` | audit findings, tech-debt inventory, physics audits and reviews, literature check                                                                                                                                         |

---

## Development notes

- After cloning, run `uv run nbstripout --install --attributes .gitattributes` once:
  notebook outputs are stripped on commit, and CI rejects notebooks with outputs and
  generators that drift from their notebooks (`scripts/check_notebook_generators.py`).
- Dependencies: edit `pyproject.toml`, run `uv lock`, commit `uv.lock`.
- The version is defined once, in `etna/_version.py`.
- Engine changes need parity runs; calibrated constants need before/after numbers
  ([USAGE §15](docs/USAGE.md#15-development-workflow)).

## License

MIT, see [LICENSE](LICENSE).

## Acknowledgements

Built on [devsim](https://devsim.org/) (J. E. Sanchez). The radiation-damage model follows
Burin et al., arXiv:2407.16710 (2024). Developed for, and calibrated against data from,
the Petringa group's 4H-SiC detector programme at INFN-LNS Catania.
