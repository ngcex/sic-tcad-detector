---
id: C-15
source: raw-C F15
title: `uv run pytest` silently falls through to a system pytest after a plain `uv sync`, because dev deps are an optional extra rather than a uv dependency group
severity: low
category: improvement
status: open
location:
  - pyproject.toml:20-24 ([project.optional-dependencies] dev = [pytest, jupyter])
  - README.md:76, :85
labels: [tooling, dx, uv]
---

## Summary

pytest lives in the `dev` optional extra. After a plain `uv sync`, which does not install
extras, `.venv/bin/pytest` does not exist, and `uv run pytest` resolves to a pytest elsewhere on
PATH (here `/opt/miniconda3/bin/pytest`, a different interpreter with no plotly). Some app test
modules then fail at collection with `ModuleNotFoundError: plotly` while others pass, which is
easy to misread as a code bug.

Note: the repo's current `.venv` has pytest, because the developer ran `uv sync --extra dev`,
so `uv run which pytest` is fine there today. This is a fresh-clone / re-sync footgun. README:76
does say `uv sync --extra dev`, but a later plain `uv sync` silently removes pytest.

## Reproduction / Evidence

Reproduced in an isolated scratch environment so the repo `.venv` was untouched:

```
S=<scratchpad>
UV_PROJECT_ENVIRONMENT=$S/venv_plain uv sync --frozen -q
ls $S/venv_plain/bin/pytest                 -> No such file or directory
UV_PROJECT_ENVIRONMENT=$S/venv_plain uv run --frozen which pytest
                                            -> /opt/miniconda3/bin/pytest
UV_PROJECT_ENVIRONMENT=$S/venv_plain uv run --frozen pytest tests/test_app_csv_export.py -q
                                            -> ERROR tests/test_app_csv_export.py ... 1 error during collection
```

(The scratch venv was deleted afterwards.)

## Impact

Contributor and CI confusion. Tests run under the wrong interpreter, with partial false failures.

## Proposed fix

```diff
 # pyproject.toml
-[project.optional-dependencies]
-dev = [
-    "pytest>=7.0",
-    "jupyter>=1.0",
-]
+[dependency-groups]
+dev = [
+    "pytest>=7.0",
+    "jupyter>=1.0",
+]
+# (optionally keep [project.optional-dependencies].dev too, for pip users)
```

`uv sync` installs the `dev` group by default. Then `uv lock`. In README, change line 76 to
`uv sync`. Alternatively, always use `uv run --extra dev pytest`.

## Test plan

In a fresh scratch env (`UV_PROJECT_ENVIRONMENT=<tmp> uv sync && uv run which pytest`), pytest
resolves inside the env and `uv run pytest tests/test_app_csv_export.py -q` collects and passes.

## Effort

S
