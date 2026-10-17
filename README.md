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

---

## What it does

- **Electrostatics & C–V** — built-in potential, depletion width, electric
  field, capacitance–voltage extraction. R² = 0.998 is fit-recovery (same
  points used to calibrate the doping profile); the only literature check
  (Bruzzi & Verroi 2023, a Schottky diode on a different, 30 µm epi wafer at 2.46e15) is a
  solver-scale check, not a validation of this device: the +3.1 % at −40 V is the Schottky-vs-pn
  built-in-voltage difference (<1 % with V_bi matched). The measured C-V of this device
  (Photons Fig. 4b, digitised in `lab/cv_gap/`) puts the knee at −1…−4 V where the model has it
  at −8…−10 V — see `.issues/perf/PERF-05-default-stack-followups.md` and
  `.planning/VALIDATION_LITERATURE.md`.
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

22 Jupyter notebooks under `notebooks/` reproduce the underlying figures and
results; the Streamlit UI (`app/`) exposes the same workflows interactively
without needing to read or write any Python.

---

## Streamlit UI

```bash
streamlit run app/main.py
```

An 8-page interactive app built on the `etna` public API — configure a device
once in the sidebar (persists across pages) and explore:

| Page                        | Workflow                                                          |
| --------------------------- | ----------------------------------------------------------------- |
| **C-V Analysis**            | C–V curve, Mott-Schottky, depletion width, full-depletion marker  |
| **Charge Collection (CCE)** | CCE vs reverse bias, configurable min/max bias sweep range        |
| **Field Map**               | 1D/2D field & doping at a chosen reverse bias, heatmap viewer     |
| **Radiation Damage**        | CCE vs proton fluence, energy/bias sweep, κ(E) disclosure banner  |
| **Dark Current**            | Total leakage vs temperature (calibrated budget estimate)         |
| **Microdosimetry**          | MC-events CSV upload → y·d(y) spectrum, y_F / y_D readouts        |
| **Batch Sweep**             | Parametric sweep over any device field, overlaid curves, bulk CSV |

Every result page shows which solver produced it (batched engine, 2D engine
or devsim) and whether the device lies inside the engine/devsim cross-check domain; the CSV
download carries the same provenance (solver, cross-check domain, device config,
software version) in its `#` header.

---

## Quick start

Requires Python 3.13+ and the `devsim` solver. Dependency management uses
[`uv`](https://docs.astral.sh/uv/).

```bash
# clone
git clone https://github.com/ngcex/sic-tcad-detector.git
cd sic-tcad-detector

# install the package (devsim, numpy, scipy, numba, matplotlib, plotly, streamlit, pandas + dev extras)
uv sync --extra dev

# use the library directly
uv run python -c "from etna import DeviceConfig, run_cv; print(run_cv(DeviceConfig()))"

# or launch the interactive UI
uv run streamlit run app/main.py

# run the fast tests
uv run pytest tests/test_api_*.py tests/test_app_*.py -q

# open a notebook
uv run jupyter lab notebooks/01_phase1_validation.ipynb
```

> **Parallel tests:** `uv run pytest -n auto --dist loadfile -q` runs the full suite
> (~700 tests) in about 4 minutes; use `-n0` when debugging with pdb.
>
> **devsim note:** the full drift-diffusion test suite is slow, and stacking
> many DD device builds in one interpreter process can exhaust devsim's
> resources. Run DD-heavy modules one file at a time (e.g.
> `uv run pytest tests/test_api_cce.py -m slow -q`) rather than the whole
> suite in a single `pytest -q` invocation.

### Python API

```python
from etna import DeviceConfig, run_cv, run_cce, run_field, run_radiation_damage

cfg = DeviceConfig(epi_thickness_um=10.0, N_A=1e19)
result = run_cv(cfg)          # SimResult: bias (x), capacitance (y), metadata
result = run_cce(cfg, v_start=-1.0, v_stop=-40.0)
result = run_field(cfg)       # includes result.mesh (MeshData) for 2D devices
result.metadata["solver"]     # "batched" | "engine2d" | "devsim"
```

**Solvers.** `run_cv` / `run_cce` (and their `run_*_batch` twins used by
`ParametricSweep`) run on a batched numba 1D engine (`etna.core.batched`) when
the device lies in the cross-check parameter box, where it was cross-checked
against devsim, and fall back to devsim otherwise or on non-convergence. 2D
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
`ParametricSweep`).

---

## Repository layout

| Path            | Contents                                                                                                                                                                                                                                                     |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `etna/`         | Installable Python package: `etna/api/` (public facades — `DeviceConfig`, `run_cv`, `run_cce`, …) and `etna/core/` (internal solver modules — material params, Poisson, drift-diffusion, CCE, dark current, radiation damage, microdosimetry, 2D devices, …) |
| `app/`          | Streamlit UI: `app/main.py` (navigation shell), `app/workflows/` (one page per simulation type), `app/components/` (sidebar, plotting, CSV export)                                                                                                           |
| `notebooks/`    | 22 analysis notebooks (validation, characterisation, FLASH, radiation damage, microdosimetry, feasibility)                                                                                                                                                   |
| `tests/`        | 57 pytest modules (~700 tests) — `test_api_*.py` (library facades), `test_app_*.py` (Streamlit AppTest, incl. un-mocked per-page smoke tests), plus the core-physics modules                                                                                 |
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
supply.

---

## Scientific status & known limitations

This code has undergone deep physics audits (`.planning/PHYSICS_AUDIT_v*.md`),
a full v5.0 integration audit (`.planning/milestones/v5.0-MILESTONE-AUDIT.md`,
21/25 requirements satisfied, 4 partial with disclosed upstream cause), and a
literature validation pass (`.planning/VALIDATION_LITERATURE.md`) comparing
simulated C–V against an open-access 4H-SiC Schottky device on a different wafer
(Bruzzi & Verroi 2023): a solver-scale check only (+3.1 % at −40 V is the Schottky-vs-pn
built-in-voltage difference; <1 % with V_bi matched), not a validation of this device. In the interest of honesty, the current limitations are:

- ✅ **Solid / publishable:** electrostatics, C–V, CCE-vs-bias, the calibrated
  graded-doping 2D model (converges to −50 V; 2D C–V matches the 1D twin,
  R² ≥ 0.99 over 0…−50 V).
- ⚠️ **Calibration, not prediction:** dark current is a calibration input: one
  trap density N_t fitted to the experimental bound < 18 pA up to −60 V (4 mm²;
  N_t = 2.805336e12 on the default stack, E_t = 0.67 eV); quote it as a budget
  estimate only.
- ⚠️ **FLASH dose-rate:** the high-injection _plasma-recombination_ physics is
  **not implemented** — FLASH dose-rate outputs are exploratory sensitivity
  bounds, not a validated mechanistic prediction.
- ✅ **Device stack (default since the n⁺ flip):** p⁺ anode 0.3 µm / 10 µm graded
  n⁻ epi / n⁺ substrate (1 µm modelled, 10¹⁹ cm⁻³, cathode on its far end), doping
  solved exactly on the 0 / −10 V C-V points (`etna.core.device.DEFAULT_STACK`): those
  two points are **calibration inputs**, not reproduced validation points. `run_cv`
  returns the small-signal C = dQ/dV (the observable of a measured C-V; W = ε·A/C in
  `metadata["W_from_C"]`); the n/N_D = 0.5 depletion edge is metadata
  (`depletion_widths`). There is no deep-bias convergence wall (devsim converges to
  ≥ −150 V; the engine is cross-checked to −100 V). Reach-through is **12.3 V** by a
  _definition_: the depletion edge within one Debye length of the n⁻/n⁺ step
  (`metadata["reach_through_V"]`; paper: ~10 V; 17 V under the legacy mesh rule).
- ⚠️ **Known C-V residual:** ε·A/C at −30 V is 9.98 µm vs the measured 9.73 µm
  (+2.6 %; Photons Fig. 4b labels 1.68 / 9.52 / 9.73 µm), the one out-of-fit C-V point.
  A ~0.5 µm graded n⁻/n⁺ buffer, an epi ~2.5 % thinner than nominal and edge/fringe
  capacitance falling with bias each close the gap; three points cannot separate them,
  so none is picked and nothing is refitted (lab/microdosimeter/notes/LOG.md §4, §6).
- ⚠️ **Legacy model:** `etna.legacy_config()` (UI: "historical model") reproduces the
  pre-flip numbers: cathode directly on the n⁻ epi, p⁺ 1 µm, 2.93e15 / 8.82e13 /
  0.987 µm, N_t 2.8e12, E_t 0.65. Only there does devsim stop converging at deep bias
  (≈ −64 V) — an artefact of the missing n⁺ substrate, **not** a physical limit.
  Direct core calls (`etna.core.*` builders and sweeps) default to the same stack; pass
  `nplus_thickness_cm=None, substrate_thickness_cm=1e-4` (or
  `device_kwargs(legacy_config())`) for the legacy geometry.
- 🔧 **Data-blocked (machinery ready, real data needed):**
  - Tissue-equivalence κ(E): the Bragg-additivity machinery is in place
    (`compute_kappa_table(source="bragg")`) but `data/srim/*.csv` are placeholders
    — drop in real NIST PSTAR proton stopping powers (see `data/srim/README.md`).
    Surfaced as a persistent warning banner on the Radiation Damage page.
  - NIEL hardness factors for SiC are placeholders pending SR-NIEL data.

Residual gaps between the notebook workflows and the Streamlit UI (7 of 21
workflows have no dedicated page yet, e.g. FLASH transient, multi-defect
comparison, validation-vs-literature) are tracked as v6.0 tech debt in the
milestone audit above — not silent omissions.

---

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgements

Built on [`devsim`](https://devsim.org/) (J. E. Sanchez). Radiation-damage model
follows Burin et al., arXiv:2407.16710 (2024). Developed for and calibrated
against experimental data from the Petringa group's 4H-SiC detector program at
INFN-LNS Catania.

## Development note: notebook outputs

Notebook outputs are stripped on commit (`.gitattributes` + `nbstripout`). After cloning run
`uv run nbstripout --install --attributes .gitattributes` once so your local git applies the filter.
Re-run a notebook (`jupyter nbconvert --execute`) to regenerate its outputs locally.

## Development environment

`uv sync --extra dev --extra viz` installs the test tools (pytest-xdist: `pytest -n auto --dist loadfile`) and the optional
animation dependency (manim, needs `pkg-config`, `cairo`, `pango`; on macOS `brew install pkgconf cairo pango`). A plain
`uv sync --extra dev` removes manim again; the Export page then simply hides the animation section.
