# SiC-TCAD — 4H-SiC Radiation-Detector Simulator

A Python TCAD (Technology Computer-Aided Design) toolkit for modelling the
electrical, transient, and radiation-damage behaviour of **4H-SiC p-n junction
radiation detectors**, aimed at proton dosimetry and microdosimetry
applications.

It is **not** a commercial TCAD package: it is a self-contained simulator (~15k
lines) built on the open-source [`devsim`](https://devsim.org/) finite-volume
semiconductor solver, with calibrated 4H-SiC material physics and a set of
analysis notebooks.

---

## What it does

- **Electrostatics & C–V** — built-in potential, depletion width, electric
  field, capacitance–voltage extraction. Validated against measured C–V
  (R² = 0.998).
- **Charge collection (CCE)** — drift-diffusion transport, Hecht analysis,
  CCE vs bias, single-particle and transient response.
- **Dark current** — trap-assisted generation, temperature dependence
  (calibrated, single-point).
- **Radiation damage** — defect introduction (Z1/2, EH4, EH6/7), carrier-lifetime
  degradation, CCE vs fluence, annealing kinetics (Burin et al. 2024 model).
- **Microdosimetry** — lineal-energy spectra, y_F / y_D, tissue-equivalence
  scaling, Monte-Carlo (Geant4) coupling.
- **2D devices & alternative structures** — graded-doping 2D solver, mesa,
  3D-electrode, ΔE-E telescope, guard-ring; parametric optimisation.

22 Jupyter notebooks under `notebooks/` reproduce the figures and results.

---

## Quick start

Requires Python 3.11+ and the `devsim` solver.

```bash
# clone
git clone https://github.com/ngcex/<repo-name>.git
cd <repo-name>

# install dependencies (devsim, numpy, scipy, matplotlib, pytest)
pip install -r requirements.txt          # or: uv pip install -r requirements.txt

# run the fast tests
pytest -q                                 # devsim DD tests are slow; see note below

# open a notebook
jupyter lab notebooks/01_phase1_validation.ipynb
```

> **devsim note:** the full drift-diffusion test suite is slow and stacking many
> DD device builds in one interpreter can exhaust devsim's process resources.
> Run DD-heavy test classes one at a time (e.g.
> `pytest tests/test_device2d.py::TestCalibrationCV`) rather than the whole file.

---

## Repository layout

| Path            | Contents                                                                                                                          |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `src/`          | Simulator modules (material params, Poisson, drift-diffusion, CCE, dark current, radiation damage, microdosimetry, 2D devices, …) |
| `notebooks/`    | 22 analysis notebooks (validation, characterisation, FLASH, radiation damage, microdosimetry, feasibility)                        |
| `tests/`        | 25 pytest modules                                                                                                                 |
| `data/`         | Material/stopping-power tables (`data/srim/` holds **placeholders** for κ — see below)                                            |
| `deliverables/` | Foundry-facing detector design specs (PDF + Markdown)                                                                             |
| `figures/`      | Generated publication figures                                                                                                     |

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

This code has undergone two deep physics audits (`.planning/PHYSICS_AUDIT_v*.md`).
In the interest of honesty, the current limitations are:

- ✅ **Solid / publishable:** electrostatics, C–V, CCE-vs-bias, the calibrated
  graded-doping 2D model (converges to −50 V, matches C–V R² ≥ 0.99 over 0…−50 V).
- ⚠️ **Calibration, not prediction:** dark current is a single-point fit; quote it
  as a budget estimate only.
- ⚠️ **FLASH dose-rate:** the high-injection _plasma-recombination_ physics is
  **not implemented** — FLASH dose-rate outputs are exploratory sensitivity
  bounds, not a validated mechanistic prediction.
- 🔧 **Data-blocked (machinery ready, real data needed):**
  - Tissue-equivalence κ(E): the Bragg-additivity machinery is in place
    (`compute_kappa_table(source="bragg")`) but `data/srim/*.csv` are placeholders
    — drop in real NIST PSTAR proton stopping powers (see `data/srim/README.md`).
  - NIEL hardness factors for SiC are placeholders pending SR-NIEL data.

---

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgements

Built on [`devsim`](https://devsim.org/) (J. E. Sanchez). Radiation-damage model
follows Burin et al., arXiv:2407.16710 (2024).
