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
  field, capacitance–voltage extraction. Validated against measured C–V
  (R² = 0.998).
- **Charge collection (CCE)** — drift-diffusion transport, Hecht analysis,
  CCE vs bias, single-particle and transient response.
- **Dark current** — trap-assisted generation (SRH/TAT/SRV decomposition),
  temperature dependence (calibrated, single-point).
- **Radiation damage** — defect introduction (Z1/2, EH4, EH6/7), carrier-lifetime
  degradation, CCE vs fluence, annealing kinetics (Burin et al. 2024 model).
- **Microdosimetry** — lineal-energy spectra, y_F / y_D, tissue-equivalence
  scaling, Monte-Carlo (Geant4) coupling.
- **2D devices & alternative structures** — graded-doping 2D solver, mesa,
  3D-electrode, ΔE-E telescope, guard-ring; parametric optimisation and batch
  sweeps.

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
| **C-V Analysis**            | C–V curve, Mott-Schottky plot, depletion width, CSV download      |
| **Charge Collection (CCE)** | CCE vs reverse bias, configurable min/max bias sweep range        |
| **Field Map**               | 1D/2D electric field & doping, quantity-selector heatmap          |
| **Radiation Damage**        | CCE vs proton fluence, energy/bias sweep, κ(E) disclosure banner  |
| **Dark Current**            | J_SRH / J_TAT / J_SRV decomposition vs temperature                |
| **Microdosimetry**          | MC-events CSV upload → y·d(y) spectrum, y_F / y_D readouts        |
| **Batch Sweep**             | Parametric sweep over any device field, overlaid curves, bulk CSV |

Every page's results are downloadable as CSV with device-config and
software-version metadata embedded in the header.

---

## Quick start

Requires Python 3.13+ and the `devsim` solver. Dependency management uses
[`uv`](https://docs.astral.sh/uv/).

```bash
# clone
git clone https://github.com/ngcex/sic-tcad-detector.git
cd sic-tcad-detector

# install the package (devsim, numpy, scipy, matplotlib, plotly, streamlit, pandas + dev extras)
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
```

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
| `tests/`        | 46 pytest modules — `test_api_*.py` (library facades), `test_app_*.py` (Streamlit AppTest), plus the original core-physics test modules                                                                                                                      |
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

This code has undergone deep physics audits (`.planning/PHYSICS_AUDIT_v*.md`)
and a full v5.0 integration audit (`.planning/milestones/v5.0-MILESTONE-AUDIT.md`,
21/25 requirements satisfied, 4 partial with disclosed upstream cause). In the
interest of honesty, the current limitations are:

- ✅ **Solid / publishable:** electrostatics, C–V, CCE-vs-bias, the calibrated
  graded-doping 2D model (converges to −50 V, matches C–V R² ≥ 0.99 over 0…−50 V).
- ⚠️ **Calibration, not prediction:** dark current is a single-point fit; quote it
  as a budget estimate only.
- ⚠️ **FLASH dose-rate:** the high-injection _plasma-recombination_ physics is
  **not implemented** — FLASH dose-rate outputs are exploratory sensitivity
  bounds, not a validated mechanistic prediction.
- ⚠️ **Deep reverse-bias non-convergence:** the 1D drift-diffusion solver's
  `ramp_bias` does not converge past ≈ −60 to −66 V for the default device —
  a solver-robustness issue distinct from the (fixed) 2D doping calibration.
  CCE and Field Map pages surface this as a graceful error, not a crash.
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
follows Burin et al., arXiv:2407.16710 (2024). Developed for and validated
against experimental data from the Petringa group's 4H-SiC detector program at
INFN-LNS Catania.
