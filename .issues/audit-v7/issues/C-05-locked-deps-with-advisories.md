---
id: C-05
source: raw-C F5
title: uv.lock pins transitive dependencies with published advisories (runtime: urllib3, anyio, gitpython; dev: tornado, jupyterlab, jupyter-server, notebook, mistune, soupsieve)
severity: medium
category: security
status: open
location:
  - uv.lock:30-31 (anyio 4.14.1), :548-549 (gitpython 3.1.50), :2223-2224 (urllib3 2.7.0)
  - pyproject.toml:10-24 (no audit step; no floor pins for transitive deps)
labels: [security, dependencies, supply-chain]
---

## Summary

pip-audit against the locked dependency set reports advisories in several transitive packages.
None is imported directly by etna, app or scripts
(`grep -rn "import requests|import urllib3|import anyio|import git" etna app scripts` returns nothing).

Corrections to the raw finding, from my own re-audit:

- anyio has 3 advisories (PYSEC-2026-4023/4024/4025), not 1.
- The raw file lists tornado 6.5.7 among "directly declared deps with no findings". That is wrong twice. Tornado is not declared in pyproject.toml, and it has 6 advisories. It is pulled in only via the `dev` extra (jupyter/ipykernel), not by Streamlit 1.58, which uses starlette.
- The dev extra also brings jupyterlab, jupyter-server, notebook, mistune and soupsieve with findings.

## Reproduction / Evidence

Runtime set (no extras):

```
uv export --frozen --no-hashes > req_C_runtime.txt
uvx pip-audit --no-deps --disable-pip -r req_C_runtime.txt
  anyio      4.14.1  x3  (PYSEC-2026-4023/4024/4025, fix 4.14.2)
  gitpython  3.1.50  x32 (PYSEC-2026-3783.., fixes 3.1.51-3.1.59)
  urllib3    2.7.0   x3  (PYSEC-2026-4175/4176/4177, fix 2.8.0)
```

All extras (`uv export --frozen --no-hashes --all-extras`) add:

```
  tornado        6.5.7  x6  (PYSEC-2026-3928, GHSA-chx6-46f5-w4vp, GHSA-c2m8-h5v5-343r,
                             GHSA-3hv7-mjh2-fv65, GHSA-wwv5-g3v4-889x, GHSA-8423-8fgw-73vq; fix 6.5.9)
  jupyterlab     4.6.1  x12 (fix 4.6.4)
  jupyter-server 2.20.0 x1  (fix 2.21.0)
  notebook       7.6.0  x1  (fix 7.6.3)
  mistune        3.3.2  x1  (fix 3.3.3)
  soupsieve      2.8.4  x2  (fix 2.9.0)
```

Dependency paths (`uv tree --frozen --invert --package ...`): urllib3 <- requests <- streamlit
and jupyterlab-server; gitpython and anyio <- streamlit; tornado <- ipykernel / jupyter-client
<- jupyter (etna extra: dev).

## Impact

Low demonstrated exploitability. The app makes no outbound HTTP, and gitpython is used only by
Streamlit's repo-info display. But the lock is what gets installed and deployed. The jupyter
server stack is network-facing when notebooks are run, and the Streamlit server binds all
interfaces (C-16).

## Proposed fix

```
uv lock --upgrade-package urllib3 --upgrade-package anyio --upgrade-package gitpython \
        --upgrade-package tornado --upgrade-package jupyterlab --upgrade-package jupyter-server \
        --upgrade-package notebook --upgrade-package mistune --upgrade-package soupsieve
uv sync --extra dev
```

Add an audit target, for example a CI step or Makefile:

```
audit:
	uv export --frozen --no-hashes --all-extras > /tmp/req.txt && uvx pip-audit --no-deps --disable-pip -r /tmp/req.txt
```

## Test plan

- After the upgrade, the pip-audit command above reports "No known vulnerabilities found" for both the runtime and all-extras exports.
- The full test suite (`uv run --extra dev pytest -q`, slow tests included) still passes, because the upgrades touch Streamlit's transitive stack.

## Effort

S
