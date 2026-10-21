## What and why

<!-- One topic per PR. Link the issue: "Closes #N" / "Part of #N". -->

## Evidence

<!-- Numbers before/after, parity vs devsim or an analytic limit, convergence study, timings (with machine). Paste commands that reproduce them. -->

## Checklist

- [ ] Tests added or updated (each new test < ~10 s; deep variants marked `exhaustive`)
- [ ] `uv run pytest -n auto -q` green locally (list any failures that already exist on `main`)
- [ ] `uv run --no-sync ruff check .` clean; no unrelated reformatting
- [ ] Docs: short section or note in `docs/USAGE.md` for user-facing changes; `CHANGELOG.md` entry under `[Unreleased]`
- [ ] Benchmarks / studies under `lab/<topic>/` with the numbers in its `notes/LOG.md`
- [ ] New runtime data in `data/` is loaded via `etna._paths.data_path` and listed in the pyproject `force-include` table
- [ ] Dependencies changed only through `pyproject.toml` + `uv lock` (no pip)
- [ ] **Provenance**: every number, constant and reference has a source that was actually opened; anything else is marked _unverified_
- [ ] **Label honesty**: no "calibration", "model" or "relative shape only" result is presented as "validation"; changed defaults keep the previous behaviour available as an option
