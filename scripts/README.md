# scripts/

Run everything from the repository root with `uv run python scripts/<name>.py`.

## Notebook generators: which file is the source of truth

Decision of 2026-10-18 (audit-v7 D-06, D-07). A generator is kept only if it reproduces
the committed notebook cell for cell. `scripts/check_notebook_generators.py` checks
this without writing or executing anything (exit 1 on drift); run it after editing
either side.

| Notebook                                   | Source of truth                                  | Note                                                                                                                                                                                                           |
| ------------------------------------------ | ------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 03_charge_collection                       | `create_notebook_03.py`                          | in sync; writes `../notebooks/` (run from `scripts/`)                                                                                                                                                          |
| 04_flash_recombination                     | the notebook                                     | generator deleted: it predated the v5 honesty caveats (null result of the SRH+Auger model, no plasma physics)                                                                                                  |
| 05b_flash_parametric_studies               | the notebook                                     | generator deleted: notebook has the cache-completeness check the generator lacked                                                                                                                              |
| 08_transient_flash                         | `create_notebook_08.py`                          | synced (dropped `savefig.bbox_inches`, which the notebook never had)                                                                                                                                           |
| 15_2d_electrostatics_cce                   | `create_notebook_15.py`                          | ex `create_notebook_15_v2.py`; the v1 generator (absolute `/Users/...` path, 29 cells, not the committed notebook) was deleted                                                                                 |
| 16_single_particle_cce                     | `create_notebook_16.py`                          | notebook regenerated from it on 2026-10-18 (backward-Euler estimator text, `etna.core` logger)                                                                                                                 |
| 17_mc_coupling, 18_microdosimetric_spectra | `create_notebook_17.py`, `create_notebook_18.py` | in sync; these two also EXECUTE the notebook in place after writing it                                                                                                                                         |
| 19_alternative_structures                  | `create_notebook_19.py`                          | in sync                                                                                                                                                                                                        |
| 20_feasibility_report                      | the notebook                                     | generator deleted: the notebook carries the 2026-10 corrections (operating-point ranking, edge/center = 1 by construction, outside-model-validity rows, `solver='auto'`) that the 2026-07 generator would undo |
| all others                                 | the notebook                                     | never had a generator                                                                                                                                                                                          |

Generators other than 17/18 only write the notebook source; outputs are stripped on
commit by nbstripout (`.gitattributes`). Never execute a notebook into `notebooks/`
to refresh tracked figures by accident: notebooks write `../figures/*.png` (and,
via `etna.core.plotting.save_figure`, `$ETNA_FIGURE_DIR` or `./figures` relative to the
working directory, i.e. `notebooks/figures/`, gitignored; before 2026-10-18 `etna/figures/`).

## Notebook housekeeping (2026-10-18)

- `notebooks/03_executed.ipynb` removed: the same 12 cells as `03_charge_collection.ipynb`
  with the pre-fix "CCE decreases monotonically with epi" check, no outputs, unreferenced.
- `notebooks/05_parametric_studies.ipynb` renamed `05b_flash_parametric_studies.ipynb`: it
  collided with the deliverable `05_dark_current_vs_fluence.ipynb` and is not in the v3
  deliverable list (`tests/baselines/v3_frozen.json`).

## Other scripts

| Script                             | Status                                                                                                                                                                                                                      |
| ---------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `validation/cv_vs_bruzzi2023.py`   | regenerates `figures/cv_validation_bruzzi2023.png` (~25 s): solver-scale check on a different (Schottky) device, V_bi matched                                                                                               |
| `check_notebook_generators.py`     | generator/notebook drift check (above)                                                                                                                                                                                      |
| `freeze_v3_baselines.py`           | historical (v3.0): froze `tests/baselines/v3_frozen.json` with the v3 graded defaults                                                                                                                                       |
| `regression_sweep_v3_notebooks.py` | historical (v3.0): executes the deliverable notebooks in memory with cwd `notebooks/` (the notebook files are not rewritten, but their `../figures/*` side outputs ARE overwritten) |
| `diagnose_1d_2d_parity.py`         | historical (Phase 26): 1D vs 2D reverse-bias diagnostic with the v3 graded profile; imports fixed to `etna.core`                                                                                                            |
| `run_calibration_2d.py`            | historical (Phase 26): 2D Nelder-Mead refit; no longer patches `etna/core/device2d.py` (defaults are `DEFAULT_STACK`, fitted on measured C-V) and writes `26-CALIBRATION-RERUN.md` instead of overwriting the frozen record |
