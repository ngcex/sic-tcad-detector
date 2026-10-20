# ETNA — Epitaxial Tcad Numerical Analysis (for 4H-SiC detectors)

A Python TCAD (Technology Computer-Aided Design) library and Streamlit UI for
modelling the electrical, transient, and radiation-damage behaviour of
**4H-SiC p-n junction radiation detectors**, aimed at proton dosimetry and
microdosimetry applications.

It is **not** a commercial TCAD package: it is a self-contained simulator built
on the open-source [`devsim`](https://devsim.org/) finite-volume semiconductor
solver, with calibrated 4H-SiC material physics, an installable Python API
(`etna`), a Streamlit UI, and a set of analysis notebooks — developed for the
Petringa group's 4H-SiC detector work at INFN-LNS Catania.

**Start here: [`docs/USAGE.md`](docs/USAGE.md)**, the single usage guide (scope and a
"what to trust" table, installation, presets and units, every `run_*` facade with a
runnable example, solver provenance, validity limits, the app, notebooks, tests,
troubleshooting). This README is the overview.

---

## What it does

- **Electrostatics & C–V** — built-in potential, depletion width, electric
  field, capacitance–voltage extraction. The default doping profile is a FIT to the
  digitised measured C-V of this device (Photons Fig. 4b: plateau rms 0.013 µm, knee
  rms 0.24 µm), so C-V agreement is fit quality, not validation; the only literature check
  (Bruzzi & Verroi 2023, a Schottky diode on a different, 30 µm epi wafer at 2.46e15) is a
  solver-scale check, not a validation of this device: the +3.1 % at −40 V is the Schottky-vs-pn
  built-in-voltage difference (<1 % with V_bi matched). The measured C-V of this device
  (Photons Fig. 4b, digitised in `lab/cv_gap/`) puts the knee at −1…−4 V; the previous
  abrupt-n⁺ default had it at −8…−10 V, hence the 2026-10-18 default flip — see
  `.issues/perf/PERF-05-default-stack-followups.md` and `.planning/VALIDATION_LITERATURE.md`.
- **Charge collection (CCE)** — drift-diffusion transport, Hecht analysis,
  CCE vs bias, single-particle and transient response.
- **Dark current** — total leakage vs temperature from an effective
  trap-assisted generation model (calibration input: one trap density N_t fitted
  to the experimental bound < 18 pA up to −60 V, 4 mm²; bias/temperature
  dependence heuristic; a budget estimate, not a prediction).
- **Radiation damage** — defect introduction (Z1/2, EH4, EH6/7), carrier-lifetime
  degradation, CCE vs fluence, annealing kinetics (Burin et al. 2024 model).
- **Microdosimetry** — lineal-energy spectra, y_F / y_D, tissue-equivalence
  scaling, Monte-Carlo (Geant4) coupling.
- **2D devices & alternative structures** — graded-doping 2D solver, mesa,
  3D-electrode, ΔE-E telescope (guard-ring hand-scored only; edge physics not
  modelled); operating-point and batch sweeps.

21 Jupyter notebooks under `notebooks/` reproduce the underlying figures and
results; the Streamlit UI (`app/`) exposes the same workflows interactively
without needing to read or write any Python.

---

## Streamlit UI

```bash
uv run streamlit run app/main.py --server.address localhost
```

The committed `.streamlit/config.toml` binds the app to localhost when it is launched from
the repo root (Streamlit reads the config of the launch directory); `--server.address
localhost` makes that explicit from any directory. The app has no authentication: for shared
use pass `--server.address 0.0.0.0` behind an authenticating reverse proxy (audit C-16).

A 10-page interactive app built on the `etna` public API. Configure a device once in the
sidebar: geometry, doping, temperature, area, and preset checkboxes for the historical,
previous-calibrated and tabulated models. The configuration persists across pages:

| Page                        | Workflow                                                                                                 |
| --------------------------- | -------------------------------------------------------------------------------------------------------- |
| **Home**                    | Orientation and current device summary                                                                   |
| **Run All**                 | Every simulation for the sidebar device (`etna.run_all`), JSON/CSV                                       |
| **C-V Analysis**            | C–V curve, Mott-Schottky, depletion width, full-depletion marker                                         |
| **Charge Collection (CCE)** | CCE vs reverse bias, configurable min/max bias sweep range                                               |
| **Field Map**               | 1D/2D field & doping at a chosen reverse bias, heatmap viewer                                            |
| **Radiation Damage**        | CCE vs proton fluence, energy/bias sweep, κ(E) disclosure banner                                         |
| **Dark Current**            | Total leakage vs temperature (calibrated budget estimate)                                                |
| **Microdosimetry**          | MC-events CSV upload → y·d(y) spectrum, y_F / y_D readouts                                               |
| **Batch Sweep**             | Parametric sweep over a curated device field, overlaid curves, bulk CSV                                  |
| **Export**                  | Spec sheet (JSON/CSV/MD/PDF), cross-section, mesh CSV/VTK, animations, schematic GDS (REFERENCE-NOT-FAB) |

Every result page shows which solver produced it (batched engine, 2D engine
or devsim) and whether the device lies inside the engine/devsim cross-check domain; the CSV
download carries the same provenance (solver, cross-check domain, device config,
software version) in its `#` header. Page-by-page inputs and outputs are in
[`docs/USAGE.md` section 9](docs/USAGE.md#9-the-streamlit-app).

---

## Quick start

Requires Python 3.13+ and [`uv`](https://docs.astral.sh/uv/). The `devsim` solver comes
as a wheel; on Linux it needs OpenBLAS/LAPACK (e.g.
`apt-get install libopenblas-dev liblapack-dev`).

```bash
# clone
git clone https://github.com/ngcex/sic-tcad-detector.git
cd sic-tcad-detector

# install the package and the dev extra (pytest, pytest-xdist, jupyter, nbstripout, ruff);
# add --extra viz for the Export page's animations (manim) and GDSII (gdstk)
uv sync --extra dev

# use the library directly
uv run python -c "from etna import DeviceConfig, run_cv; print(run_cv(DeviceConfig()).metadata['solver'])"

# or launch the interactive UI
uv run streamlit run app/main.py --server.address localhost

# run the default test suite in parallel (the scheduler comes from pytest.ini)
uv run pytest -n auto -q

# open a notebook (start from notebooks/ so saved figures land in notebooks/figures/)
cd notebooks && uv run jupyter lab 01_phase1_validation.ipynb
```

> **Always sync with `--extra dev`.** A plain `uv sync` installs no extras and removes
> pytest from `.venv`. `uv run pytest` may then fall through to a system pytest
> (audit C-15).
>
> **Tests:** there are 1092 tests in 92 modules. `pytest.ini` deselects the 7 `exhaustive`
> tests (run them with `-m exhaustive`), so the default run is 1085 tests, and it is the
> real gate. On a 10-core machine `-n auto` takes about 4 minutes (248 s, down from
> 502 s with `--dist loadfile` on the same shared machine). `pytest.ini` supplies
> `--dist loadgroup`; `tests/_xdist_schedule.py` balances the work and
> `tests/conftest.py` deletes devsim devices that tests leak. Use `-n0` when
> debugging with pdb. See [docs/USAGE.md §14](docs/USAGE.md#14-tests).
>
> **CI:** `.github/workflows/ci.yml` runs on pushes to `main`, on pull requests and on
> manual dispatch. It checks the lock, ruff, notebook hygiene and the generators, then
> runs the default suite with `-n 4`. The exhaustive tests and `pip-audit` run nightly.

### Python API

```python
from etna import DeviceConfig, run_cv, run_cce, run_field, run_radiation_damage

cfg = DeviceConfig(epi_thickness_um=10.0, N_A=1e19)  # default area: the 4 mm2 dosimeter
result = run_cv(cfg)          # SimResult: bias (x), capacitance in F for area_cm2 (y), metadata
result = run_cce(cfg, v_start=-1.0, v_stop=-40.0)
result = run_field(cfg)       # result.mesh (MeshData); 1D: x = depth (um), y = |E| (V/cm)
result.metadata["solver"]     # "batched" | "engine2d" | "devsim"
```

Each facade, its return values and metadata keys, with an example that runs, is in
[`docs/USAGE.md` section 6](docs/USAGE.md#6-the-run_-facades).

**Solvers.** `run_cv` / `run_cce` (and their `run_*_batch` twins used by
`ParametricSweep`) run on a batched numba 1D engine (`etna.core.batched`) when
the device lies in the cross-check parameter box, where it was cross-checked
against devsim, and fall back to devsim otherwise or on non-convergence. The same
engine serves `run_dark_current` (TAT/SRH/SRV model), `run_flash_recombination`, the
temperature and fluence sweeps and the 2D dark current of laterally uniform devices
(milliseconds instead of seconds; parity with devsim <= ~1e-11). 2D
`run_field` uses a numba 2D engine (`etna.core.engine2d`) where eligible.
`metadata["solver"]`, `metadata["validated_domain"]` / `"domain_note"` record
what ran; C-V/CCE results also carry `metadata["reach_through_V"]`, the
full-depletion voltage. `DeviceConfig` rejects out-of-range values
(`ValueError`), and batch sweeps keep the completed points when one value fails
(its slot is an empty result with `metadata["error"]`).

See `etna/__init__.py` for the full public API surface (`DeviceConfig`,
`SimResult`, `MeshData`, `run_cv`, `run_field`, `run_cce`,
`run_radiation_damage`, `run_dark_current`, `run_temperature_sweep`,
`run_flash_recombination`, `run_transient`, `run_microdosimetry`,
`run_cv_batch`, `run_cce_batch`, `run_cce_tilted_alpha`, `run_all`, `ParametricSweep`,
and the presets `measured_cv_config`, `measured_cv_tabulated_config`,
`nplus_calibrated_config`, `legacy_config`, `microdosimeter_config`).

---

## Repository layout

| Path            | Contents                                                                                                                                                                                                                                                     |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `etna/`         | Installable Python package: `etna/api/` (public facades — `DeviceConfig`, `run_cv`, `run_cce`, …) and `etna/core/` (internal solver modules — material params, Poisson, drift-diffusion, CCE, dark current, radiation damage, microdosimetry, 2D devices, …) |
| `app/`          | Streamlit UI: `app/main.py` (navigation shell), `app/workflows/` (one page per simulation type), `app/components/` (sidebar, plotting, CSV export)                                                                                                           |
| `notebooks/`    | 21 analysis notebooks (validation, characterisation, FLASH, radiation damage, microdosimetry, feasibility); committed without outputs                                                                                                                        |
| `tests/`        | 92 pytest modules (1092 tests; 1085 in the default run): `test_api_*.py` (library facades), `test_app_*.py` (Streamlit AppTest, incl. un-mocked per-page smoke tests), engine/devsim parity, presets, plus the core-physics modules                          |
| `docs/`         | `docs/USAGE.md`, the usage guide                                                                                                                                                                                                                             |
| `lab/`          | Research spikes, each with `notes/LOG.md` (commands and numbers behind calibrated defaults; see `lab/README.md`)                                                                                                                                             |
| `scripts/`      | Notebook generators, the generator drift check, validation scripts (`scripts/README.md`)                                                                                                                                                                     |
| `.issues/`      | Audit findings, performance issues and the tech-debt inventory                                                                                                                                                                                               |
| `data/`         | Material/stopping-power tables (`data/srim/` holds **placeholders** for κ — see below)                                                                                                                                                                       |
| `deliverables/` | Foundry-facing detector design specs (PDF + Markdown)                                                                                                                                                                                                        |
| `figures/`      | Generated publication figures                                                                                                                                                                                                                                |
| `examples/`     | Standalone usage scripts (e.g. `examples/cv_example.py`)                                                                                                                                                                                                     |

---

## Detector design specs (foundry hand-off)

`deliverables/` contains two **process-agnostic physics specifications** intended
to be handed to a foundry (ST / FBK) as a feasibility study — not mask layouts:

- **DESIGN-1** — dosimetry p-n diode (zero-bias capable)
- **DESIGN-2** — microdosimeter sensitive-volume geometry (junction type left to fab)

Each states explicitly what is a validated model target vs. what the foundry must
supply. The Markdown files are the source; the PDFs were regenerated from them on
2026-10-18, after the default-stack flip. `deliverables/DIMENSIONS-provenance.md`
classifies every layout dimension as DOCUMENTED, MEASURED-FROM-FIGURE or PLACEHOLDER.

---

## Scientific status & known limitations

This code has undergone deep physics audits (`.planning/PHYSICS_AUDIT_v*.md`),
a full v5.0 integration audit (`.planning/milestones/v5.0-MILESTONE-AUDIT.md`,
21/25 requirements satisfied, 4 partial with disclosed upstream cause), and a
literature validation pass (`.planning/VALIDATION_LITERATURE.md`) comparing
simulated C–V against an open-access 4H-SiC Schottky device on a different wafer
(Bruzzi & Verroi 2023): a solver-scale check only (+3.1 % at −40 V is the Schottky-vs-pn
built-in-voltage difference; <1 % with V_bi matched), not a validation of this device. In the interest of honesty, the current limitations are:

- ✅ **Verified numerics (code-to-code, not against measurement):** the batched 1D and
  2D numba engines reproduce devsim inside their parameter boxes (1D dC/C ≲ 3e-8,
  dCCE ≲ 1e-10 to −100 V, dark current ≲ 1.4e-11 relative, FLASH CCE ≲ 3e-12; 2D
  potential ~1e-13 V), and the 2D graded model's C–V
  matches its 1D twin (R² ≥ 0.99 over 0…−50 V). Electrostatics, C–V and CCE-vs-bias are
  therefore solved consistently, but their agreement with this device rests on the
  C–V calibration below: the C–V is a fit, and no measured CCE of the dosimeter has been
  reproduced. See the "what to trust" table in
  [`docs/USAGE.md`](docs/USAGE.md#2-what-to-trust).
- ⚠️ **Calibration, not prediction:** dark current is a calibration input: one
  trap density N_t fitted to the experimental bound < 18 pA up to −60 V (4 mm²;
  N_t = 2.849024e12 on the default stack, E_t = 0.67 eV); quote it as a budget
  estimate only.
- ⚠️ **FLASH dose-rate:** the high-injection _plasma-recombination_ physics is
  **not implemented** — FLASH dose-rate outputs are exploratory sensitivity
  bounds, not a validated mechanistic prediction. The steady-state sweep solves
  SRH + Auger (active since 2026-10-20; before, the Auger term was silently
  dropped): Auger only adds a dose-rate-independent CCE offset (≤ 5.6e-5, n⁺ layer).
- ⚠️ **Device stack = calibration (default since 2026-10-18: the measured-C-V profile):** p⁺ anode
  0.3 µm / 10 µm graded n⁻ epi (N_D_junction 8.19e15, L 0.46 µm, N_D_bulk 3.38e13) with
  a graded n⁻/n⁺ transition over its last 1.14 µm (log-linear to 4.5e16 at 10 µm) /
  n⁺ substrate (1 µm modelled, 10¹⁹ cm⁻³, cathode on its far end)
  (`etna.core.device.DEFAULT_STACK`; `etna.measured_cv_config()` is an alias). The
  profile is a **fit to the digitised Photons Fig. 4b** (plateau rms 0.013 µm over
  −5…−59 V, knee rms 0.24 µm): the whole curve is a calibration input. **N_D_bulk
  3.4e13 is below the documented 0.5–1e14** statements, and the epi/n⁺ position is set
  by the fit, not by SIMS. `run_cv`
  returns the small-signal C = dQ/dV (the observable of a measured C-V; W = ε·A/C in
  `metadata["W_from_C"]`); the n/N_D = 0.5 depletion edge is metadata
  (`depletion_widths`). There is no deep-bias convergence wall (devsim converges to
  ≥ −150 V; the engine is cross-checked to −100 V). Full depletion, by _definition_
  (`etna.export.labels.VFD_DEFINITION`): analytic V_fd **2.74 V** (to the transition
  start), reach-through **3.10 V** (edge within one Debye length of the transition
  start, `metadata["reach_through_V"]`), sweep gate 1.2 × V_fd = **3.28 V**. The paper's
  −10 V is 3× that gate (3.6× V_fd); DESIGN-1 keeps ≥ −10 V as the practical operating
  bias (margin for the fit's doping uncertainty: V_fd 3.9–7.6 V if the bulk is at the
  documented 0.5–1e14).
- ⚠️ **Previous default:** `etna.nplus_calibrated_config()` (UI: "previous calibrated
  model") reproduces the 2026-10-17 numbers bit for bit: abrupt epi/n⁺ step, doping
  solved on the 0 / −10 V C-V points (2.958237e15 / 8.252732e13 / L 0.987 µm held),
  N_t 2.805336e12; V_fd 9.93 V, reach-through 12.3 V, ε·A/C(−30 V) 9.98 vs 9.73 µm
  (+2.6 %), C-V knee at −8…−10 V (excluded by the full digitised curve).
- ⚠️ **Legacy model:** `etna.legacy_config()` (UI: "historical model") reproduces the
  pre-flip numbers: cathode directly on the n⁻ epi, p⁺ 1 µm, 2.93e15 / 8.82e13 /
  0.987 µm, N_t 2.8e12, E_t 0.65. Only there does devsim stop converging at deep bias
  (≈ −64 V) — an artefact of the missing n⁺ substrate, **not** a physical limit.
  Direct core calls (`etna.core.*` builders and sweeps) default to the same stack
  (including the graded transition); pass `nplus_thickness_cm=None,
buffer_thickness_cm=None, buffer_end_doping=None, substrate_thickness_cm=1e-4` (or
  `device_kwargs(legacy_config())`) for the legacy geometry.
- 🔧 **Data-blocked (machinery ready, real data needed):**
  - Tissue-equivalence κ(E): the Bragg-additivity machinery is in place
    (`compute_kappa_table(source="bragg")`) but `data/srim/*.csv` are placeholders
    — drop in real NIST PSTAR proton stopping powers (see `data/srim/README.md`).
    Surfaced as a persistent warning banner on the Radiation Damage page.
  - NIEL hardness factors for SiC are placeholders pending SR-NIEL data.

Residual gaps between the notebook workflows and the Streamlit UI are tracked as v6.0
tech debt in the milestone audit above, not left as silent omissions. That audit found
7 of 21 workflows with no page; the Run All page has since added the FLASH steady-state
section (notebook 04). Notebooks 08 (FLASH transient), 12 (multi-defect comparison),
14 (validation vs literature), 16, 19 and 20 remain notebook-only.

---

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgements

Built on [`devsim`](https://devsim.org/) (J. E. Sanchez). Radiation-damage model
follows Burin et al., arXiv:2407.16710 (2024). Developed for and calibrated
against experimental data from the Petringa group's 4H-SiC detector program at
INFN-LNS Catania.

## Development note: notebook outputs

Notebook outputs are stripped on commit (`.gitattributes` + `nbstripout`). After cloning, run
`uv run nbstripout --install --attributes .gitattributes` once so that your local git applies
the filter. CI rejects tracked notebooks with outputs (`nbstripout --verify`) and generators
that drift from their notebooks (`scripts/check_notebook_generators.py`). To regenerate
outputs locally, re-run a notebook from `notebooks/` (`jupyter nbconvert --execute`). Run
from the repo root, the notebook figures would write over the tracked `figures/`; set
`ETNA_FIGURE_DIR` to redirect them. Notebook 17 also rewrites `data/synthetic_mc_events.csv`
relative to the repo root; notebook 16 loads `data/cce_let_table_*.json` and rebuilds them
(~40 min) only with `ETNA_REBUILD_CCE_TABLES=1`, into `$ETNA_FIGURE_DIR/data/`. Notebook run
times are listed in docs/USAGE.md section 12 (the slowest, 19, about 3 min).

## Development environment

The `dev` extra is the one development dependency set: pytest, pytest-xdist, jupyter,
nbstripout and ruff. `uv sync --extra dev --extra viz` installs it plus the optional
animation and layout dependencies. manim needs `pkg-config`, `cairo` and `pango`; on
macOS, `brew install pkgconf cairo pango`. gdstk is needed for the GDSII download.

A later `uv sync --extra dev` removes the `viz` packages again; the Export page then hides
the animation section and the GDSII button. The package version is defined once, in
`etna/_version.py`. The development workflow (parity tests when changing an engine,
before/after numbers for any calibrated constant) is described in
[`docs/USAGE.md` section 15](docs/USAGE.md#15-development-workflow).
