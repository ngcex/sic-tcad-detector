---
id: D-06
title: scripts/create_notebook_15.py is a dead superseded generator writing to a nonexistent absolute path
severity: low
category: dead-code
status: fixed  # 2026-10-18 sync: f21d4ff (wp-c1): dead create_notebook_15.py v1 deleted, create_notebook_15_v2.py -> create_notebook_15.py
location: scripts/create_notebook_15.py:477-480; scripts/create_notebook_15_v2.py:609
labels: [scripts, dead-code, notebooks]
source: raw-D F6
---

## Summary

`create_notebook_15.py` writes to `/Users/ngcex/projects/physics/etna/notebooks/15_2d_electrostatics_cce.ipynb`,
a directory that does not exist, so it raises `FileNotFoundError`. It is the pre-rewrite version.
`create_notebook_15_v2.py` writes to the repo-relative `notebooks/15_...` and is the generator that
`.planning/phases/20-*/20-02-SUMMARY.md` records as producing the committed notebook.
Severity is downgraded to low (raw: medium). The script is inert and nothing runs it.

## Reproduction/Evidence

```
$ sed -n 474,481p scripts/create_notebook_15.py
nbformat.write(nb, "/Users/ngcex/projects/physics/etna/notebooks/15_2d_electrostatics_cce.ipynb")
$ grep -n nbformat.write scripts/create_notebook_15_v2.py
609:nbformat.write(nb, "notebooks/15_2d_electrostatics_cce.ipynb")
$ grep -rn create_notebook_15 . (excluding .git/.issues)
-> only .planning/ phase docs (historical); no code references
```

I did not re-run the raw report's cell-diff harness (13/14 cells differing). The stale absolute path and
the planning record are enough to show v1 is superseded.

## Impact

Maintainers get confused about which generator is canonical. Running v1 crashes.

## Proposed fix

`git rm scripts/create_notebook_15.py`, then `git mv scripts/create_notebook_15_v2.py scripts/create_notebook_15.py`.
Optionally update the `.planning` references.

## Test plan

`ls scripts/create_notebook_15*.py` shows a single file. Running it from the repo root regenerates nb15
with identical code cells.

## Effort

S
