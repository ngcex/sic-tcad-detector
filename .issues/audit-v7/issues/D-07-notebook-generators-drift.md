---
id: D-07
title: Notebook generators drifted from committed notebooks (stale src.* logger names in nb16/nb20, extra rcParams in gen08, hand-edits in nb05)
severity: low
category: docs
status: fixed  # 2026-10-18 sync: 7acbe10 + 9dec765 (wp-c1): generators kept only where they reproduce the notebook (scripts/check_notebook_generators.py, all in sync), 04/05/20 deleted, decisions in scripts/README.md
location: scripts/create_notebook_05.py, scripts/create_notebook_08.py:72, scripts/create_notebook_16.py, scripts/create_notebook_20.py; notebooks/05_parametric_studies.ipynb, 08_transient_flash.ipynb, 16_single_particle_cce.ipynb, 20_feasibility_report.ipynb
labels: [notebooks, generators, drift, logging]
source: raw-D F7
---

## Summary

The generators and the committed notebooks disagree in both directions, so neither is authoritative:

- nb16 and nb20 call `logging.getLogger('src.single_particle')` / `'src.optimization'`, while their
  generators were updated to `etna.core.*`. In nb20, `setLevel(WARNING)` therefore silences a logger
  that no longer exists, and INFO output from `etna.core.optimization` gets through.
- `create_notebook_08.py:72` adds `'savefig.bbox_inches': 'tight'`, which is not in the committed notebook.
- nb05 contains a cache-completeness block that the generator does not emit (I did not verify this independently).
- Output-path conventions are inconsistent: 03/04 write `../notebooks/` (run from `scripts/`), 05/08 use
  `PROJECT_DIR`, and 15_v2/16-20 write `notebooks/` (run from the repo root).

## Reproduction/Evidence

```
$ grep -o "getLogger([^)]*)" notebooks/16_*.ipynb notebooks/20_*.ipynb scripts/create_notebook_16.py scripts/create_notebook_20.py
scripts/create_notebook_16.py:getLogger('etna.core.single_particle')
scripts/create_notebook_20.py:getLogger('etna.core.optimization')
notebooks/20_feasibility_report.ipynb:getLogger('src.optimization')
notebooks/16_single_particle_cce.ipynb:getLogger('src.single_particle')
$ grep -n "savefig.bbox" scripts/create_notebook_08.py notebooks/08_transient_flash.ipynb
scripts/create_notebook_08.py:72:  'savefig.bbox_inches': 'tight'   (no hit in the notebook)
$ grep -n "\.ipynb" scripts/create_notebook_0{3,4,5,8}.py  -> ../notebooks (03,04), PROJECT_DIR (05,08)
```

## Impact

Notebook verbosity is cosmetic. The larger problem is that there is no single source of truth, so
future edits will keep diverging.

## Proposed fix

Pick one source of truth:
(a) delete `scripts/create_notebook_*.py`, since the notebooks are committed and are the artifact, and fix
the two `getLogger` lines in nb16/nb20 directly (`'src.' -> 'etna.core.'`); or
(b) regenerate 05/08/16/20 from the generators, normalize all output paths to
`Path(__file__).resolve().parents[1] / "notebooks"`, and add a CI check that re-running each generator
reproduces the committed code cells.

## Test plan

`grep -rn "getLogger('src\." notebooks/` returns nothing. For (b), the generator-diff check passes for
every generator.

## Effort

S (a) / M (b)
