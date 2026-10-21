# Contributing to ETNA

Thank you for helping. ETNA is a research TCAD code for 4H-SiC detectors, maintained by
one person on a best-effort basis (see [SUPPORT.md](SUPPORT.md)). Bug reports, physics
questions, measurement data and pull requests are all welcome.

Everyone taking part follows the [Code of Conduct](CODE_OF_CONDUCT.md). Report security
issues privately as described in [SECURITY.md](SECURITY.md), not in public issues.

## Ways to contribute

- **Report a bug.** Use the bug template. Include the ETNA version (`etna.__version__`),
  OS, Python, the `metadata` of the result, and a minimal script.
- **Ask a physics question.** Use the physics-question template. Questions about whether
  a model applies to your device are welcome.
- **Contribute data.** Use the data-contribution template. Measured C-V, SIMS/SRP, I(T)
  or CCE data improve the calibration more than anything else.
  [docs/DATA_REQUEST_CHECKLIST.md](docs/DATA_REQUEST_CHECKLIST.md) lists what is
  needed and the file formats.
- **Send code or docs.** Open an issue first for anything larger than a small fix, so the
  approach can be agreed before you invest time.

## Development setup (uv only)

ETNA uses [uv](https://docs.astral.sh/uv/) for everything. **Do not use pip, venv or
conda** for the development environment.

```bash
git clone https://github.com/ngcex/sic-tcad-detector.git
cd sic-tcad-detector
uv sync --extra dev                    # library + pytest, xdist, jupyter, nbstripout, ruff, meshio
# optional: uv sync --extra dev --extra viz   (manim + gdstk; macOS: brew install pkgconf cairo pango)
```

Always sync with `--extra dev`. A plain `uv sync` removes the dev tools, and `uv run
pytest` may then fall through to a system pytest. On Linux, install OpenBLAS/LAPACK
first (`apt-get install libopenblas-dev liblapack-dev`).

To add a dependency, edit `pyproject.toml`, run `uv lock`, and commit `uv.lock`
together with the change (CI runs `uv lock --check`).

## Tests and benchmarks

```bash
uv run pytest -n auto -q                  # the default run: the local gate, must be green
uv run pytest -n0 -x tests/test_api_cv.py # one file, serial
uv run --no-sync ruff check .             # lint
```

Rules (details in [docs/USAGE.md §14-15](docs/USAGE.md#14-tests)):

- **Every change comes with tests.** Each new test should run in under about 10 s.
  Redundant deep variants get the `exhaustive` marker, which the default run deselects.
- **Numerics changes need evidence.** Changing an engine, a solver path or an engine box
  needs a parity check against devsim, or against an analytic limit, at the affected
  corners. A convergence study is needed where a discretisation changes.
- **Benchmarks and studies go under `lab/<topic>/`.** Put the bench script there and
  record the numbers in `lab/<topic>/notes/LOG.md`, with the command that reproduces
  them. Do not put long runs in the test suite.
- **Speed.** Use numba float64 kernels first. Use a GPU path only when it is measured to
  beat numba on a realistic workload.
- **Calibrated constants.** Never change one without before/after numbers. The previous
  default stays available as a preset and is pinned by a test.
- **Notebooks.** Commit them without outputs (`nbstripout`). Generated notebooks must
  match their `scripts/create_notebook_*.py` generator
  (`uv run python scripts/check_notebook_generators.py`).
- **Runtime data files.** A new file under `data/` that the library reads must be:
  - loaded through `etna._paths.data_path(...)`, never through a path relative to
    `__file__`;
  - added to the `[tool.hatch.build.targets.wheel.force-include]` table in
    `pyproject.toml`.

  `tests/test_packaging.py` enforces both.

## Honesty rules for physics claims

These rules matter more than style:

- **Label every result for what it is.** Use exactly one of:
  - **verified** (code to code, or against an analytic limit);
  - **calibration** (fitted to data, so agreement with that data is fit quality);
  - **validation** (compared with an independent measurement not used in the fit);
  - **prediction** or **model** (not compared with anything).

  Never upgrade a "calibration" or "relative shape only" result to "validation". The
  [README "What to trust" table](README.md#what-to-trust) is the reference.

- **Cite only sources you have opened.** Mark anything taken from a search snippet or
  from memory as **unverified**. Never invent a number, a reference or a DOI. Data
  digitised from figures records the figure, the paper and the digitisation method.
- **Report limitations.** Report truncated sweeps, fallbacks and out-of-box
  configurations through the result `metadata`, never silently.

## Documentation

- [docs/USAGE.md](docs/USAGE.md) is the single usage document. Add a short section or
  paragraph for every new user-facing feature, with a runnable example.
- Add an entry under `## [Unreleased]` in [CHANGELOG.md](CHANGELOG.md).
- The documentation site builds from `docs/`. Check it with
  `uvx --from mkdocs==1.6.1 --with mkdocs-material==9.7.7 mkdocs build --strict`.

## Pull requests

- Branch from `main`. Keep each PR to one topic. Do not reformat unrelated code.
- Fill in the PR template (tests, docs, provenance, label honesty).
- CI (`ci.yml`, `matrix.yml`) must be green. The local `uv run pytest -n auto -q` is
  the real gate.
- By contributing you agree that your contribution is licensed under the
  [MIT licence](LICENSE) of the project.
