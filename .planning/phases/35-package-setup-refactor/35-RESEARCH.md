# Phase 35: Package Setup & Refactor - Research

**Researched:** 2026-06-26
**Domain:** Python packaging (hatchling/pyproject.toml), import-path rename, editable install
**Confidence:** HIGH

---

## Summary

Phase 35 is a pure refactoring phase: move the existing flat `src/` directory to `petringa/core/`, rewrite ~330+ `from src.X` import statements across source modules, tests, scripts, and notebooks, and introduce a `pyproject.toml` with hatchling as the build backend. No physics logic changes. No new dependencies (all eight runtime deps are already installed in the project venv). The acceptance gate — `pytest -q` green on all 25 test modules with `tests/baselines/v3_frozen.json` byte-for-byte unchanged — requires a fully working devsim installation because the slow-marked integration tests run on a bare `pytest -q` invocation.

The single architectural complexity is a scope boundary: Success Criterion #2 (`python -c "from petringa import DeviceConfig"`) requires `DeviceConfig` to exist in the package, but `DeviceConfig` is formally assigned to Phase 36 (LIB-01). Resolution strategy: ship a minimal `DeviceConfig` skeleton dataclass in `petringa/__init__.py` during Phase 35 so the smoke test passes, and hand off the full implementation to Phase 36. The planner must make this boundary explicit.

The rename scope is wider than just `src/` → `petringa/core/`: 17 of 24 source modules also import each other with `from src.X`, scripts contain 59 such imports (some as string literals embedded inside notebook-creation scripts), and 22 notebooks contain 81 more. All must be addressed — failing to touch the string-literal `patch("src.single_particle....")` call in `tests/test_mc_coupling.py` will silently break mock assertions after the rename.

**Primary recommendation:** Mechanical textual rewrite — `from src.` → `from petringa.core.` and `import src.` → `import petringa.core.` — using a reliable sed-or-Python script, then hand-fix the six string/path/logger-name references. Do not switch to relative imports; the diff would be larger and harder to review without any benefit for this phase.

---

## User Constraints (from CONTEXT.md)

> No CONTEXT.md exists for Phase 35. Constraints are drawn from STATE.md and REQUIREMENTS.md.

### Locked Decisions (from STATE.md)

- Package name: `petringa` (installable with `pip install -e .` / `uv pip install -e .`)
- Build backend: hatchling via `pyproject.toml` (replaces `requirements.txt`)
- Public API lives in `petringa/api/`; internal modules in `petringa/core/` are not public contract
- Source rename: `src/` → `petringa/core/`
- Acceptance gate: `pytest -q` green + `v3_frozen.json` baseline byte-for-byte unchanged
- Tool for all installs: **uv** (not pip/venv) — per project memory
- No physics changes in any v5.0 phase — refactor only

### Deferred Ideas (OUT OF SCOPE for Phase 35)

- Full `DeviceConfig`, `SimResult`, `MeshData` dataclass implementation (Phase 36)
- `run_cv()`, `run_field()`, `run_cce()` facade implementation (Phases 36-37)
- Streamlit UI (Phases 38-43)
- Any new physics model

---

<phase_requirements>

## Phase Requirements

| ID     | Description                                                                                                                                                                 | Research Support                                                                                                |
| ------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- |
| PKG-01 | Developer can install with `uv pip install -e .` from a `pyproject.toml` with hatchling build backend, replacing `requirements.txt`                                         | pyproject.toml structure verified from hatch docs; hatchling 1.30.1 available on PyPI                           |
| PKG-02 | Package declares all runtime deps (`devsim`, `numpy`, `scipy`, `matplotlib`, `plotly`, `streamlit`, `pandas`) and optional `[dev]` extras (`pytest`, `jupyter`)             | All deps confirmed present in .venv; pyproject.toml optional-dependencies pattern documented                    |
| PKG-03 | All 25 pytest modules pass unchanged after refactor (only import paths updated `src.X` → `petringa.core.X`); `tests/baselines/v3_frozen.json` regression baseline unchanged | 25 test files confirmed; import-scope audit complete; baseline file confirmed at tests/baselines/v3_frozen.json |

</phase_requirements>

---

## Architectural Responsibility Map

| Capability                                 | Primary Tier                    | Secondary Tier            | Rationale                                                                 |
| ------------------------------------------ | ------------------------------- | ------------------------- | ------------------------------------------------------------------------- |
| Package build metadata                     | Build system (pyproject.toml)   | —                         | Standard Python packaging; hatchling reads pyproject.toml at install time |
| Module discovery (editable install)        | Build system (hatchling)        | —                         | hatchling wheel target's `packages` key controls what lands on sys.path   |
| Import-path rewrite                        | Source files (all .py + .ipynb) | —                         | Mechanical text transform; no runtime tier boundary                       |
| Public API export (`petringa/__init__.py`) | Package root                    | petringa/api/ (Phase 36+) | `__init__.py` re-exports symbols; internal modules stay in core/          |
| Regression guard (v3_frozen.json)          | Test layer                      | devsim runtime            | Slow tests call devsim; devsim must be functional for PKG-03 to pass      |

---

## Standard Stack

### Core

| Library   | Version                                                     | Purpose                                      | Why Standard                                                                |
| --------- | ----------------------------------------------------------- | -------------------------------------------- | --------------------------------------------------------------------------- |
| hatchling | 1.30.1 [VERIFIED: PyPI via `pip3 index versions hatchling`] | Build backend for pyproject.toml             | PEP 517/518 compliant; lightweight; declared in STATE.md as locked decision |
| uv        | 0.9.8 [VERIFIED: `uv --version`]                            | All installs including `uv pip install -e .` | Project memory mandates uv over pip/venv                                    |

### Supporting

| Library    | Version                                                                         | Purpose                                         | When to Use                                           |
| ---------- | ------------------------------------------------------------------------------- | ----------------------------------------------- | ----------------------------------------------------- |
| devsim     | 2.10.0 [VERIFIED: `.venv/lib/python3.13/site-packages/devsim-2.10.0.dist-info`] | TCAD solver; required by slow integration tests | Runtime dep; must be importable for PKG-03 slow tests |
| numpy      | 2.4.3 [VERIFIED: venv]                                                          | Array computation                               | Runtime dep declaration in pyproject.toml             |
| scipy      | 1.17.1 [VERIFIED: venv]                                                         | Scientific routines                             | Runtime dep declaration                               |
| matplotlib | 3.10.8 [VERIFIED: venv]                                                         | Plotting                                        | Runtime dep declaration                               |
| plotly     | not yet in venv [ASSUMED]                                                       | Interactive plots                               | Needed for UI phases; declare now per PKG-02          |
| streamlit  | not yet in venv [ASSUMED]                                                       | UI framework                                    | Needed for UI phases; declare now per PKG-02          |
| pandas     | 3.0.1 [VERIFIED: venv]                                                          | Data frames                                     | Runtime dep declaration                               |
| pytest     | 9.0.2 [VERIFIED: venv]                                                          | Test runner                                     | Optional [dev] extra per PKG-02                       |
| jupyter    | 1.1.1 [VERIFIED: venv]                                                          | Notebook execution                              | Optional [dev] extra per PKG-02                       |

**Installation (for setting up fresh venv):**

```bash
uv pip install -e ".[dev]"
```

**Version verification commands run:**

```bash
# Python version
python3 --version                    # 3.13.5
uv --version                         # 0.9.8
pip3 index versions hatchling        # 1.30.1 current
ls .venv/lib/python3.13/site-packages/ | grep -E "devsim|numpy|scipy|matplotlib|pandas|pytest|jupyter"
```

---

## Package Legitimacy Audit

> All packages in Phase 35 are already-declared dependencies (either in requirements.txt or venv). No new packages are introduced by this phase. slopcheck not run because the install list is identical to the existing runtime environment.

| Package                                      | Registry | Age                           | Downloads                      | Source Repo                    | slopcheck | Disposition                               |
| -------------------------------------------- | -------- | ----------------------------- | ------------------------------ | ------------------------------ | --------- | ----------------------------------------- |
| hatchling                                    | PyPI     | 4+ yrs [ASSUMED: established] | Very high [ASSUMED]            | github.com/pypa/hatch          | not run   | Approved — official PyPA project          |
| devsim                                       | PyPI     | 5+ yrs [ASSUMED]              | Low/niche [ASSUMED]            | github.com/devsim/devsim       | not run   | Approved — already installed and tested   |
| numpy/scipy/matplotlib/pandas/pytest/jupyter | PyPI     | 10+ yrs                       | Hundreds of millions [ASSUMED] | well-known                     | not run   | Approved — tier-1 scientific Python stack |
| plotly                                       | PyPI     | 8+ yrs [ASSUMED]              | Very high [ASSUMED]            | github.com/plotly/plotly.py    | not run   | Approved — widely used                    |
| streamlit                                    | PyPI     | 6+ yrs [ASSUMED]              | High [ASSUMED]                 | github.com/streamlit/streamlit | not run   | Approved — widely used                    |

**Packages removed due to slopcheck [SLOP] verdict:** none
**Packages flagged as suspicious [SUS]:** none

_slopcheck was not run; all packages are either already in the project venv or are blue-chip ecosystem libraries. All tagged [ASSUMED] for install age/download counts; the packages themselves are verified by venv presence._

---

## Architecture Patterns

### System Architecture Diagram

```
Before Phase 35:
  project root
  ├── src/                    ← flat, importable via sys.path (root on path)
  │   ├── __init__.py         ← empty
  │   ├── device.py
  │   ├── device2d.py
  │   └── ... (22 more .py)
  ├── tests/                  ← from src.X import Y
  ├── scripts/                ← from src.X import Y + sys.path.insert(0, root)
  ├── notebooks/              ← from src.X import Y
  └── requirements.txt

After Phase 35:
  project root
  ├── petringa/               ← installable Python package
  │   ├── __init__.py         ← exports DeviceConfig stub + future public API
  │   ├── core/               ← renamed from src/
  │   │   ├── __init__.py
  │   │   ├── device.py
  │   │   ├── device2d.py
  │   │   └── ... (22 more .py, internal imports rewritten)
  │   └── _version.py
  ├── tests/                  ← from petringa.core.X import Y
  ├── scripts/                ← from petringa.core.X import Y
  ├── notebooks/              ← from petringa.core.X import Y
  └── pyproject.toml          ← replaces requirements.txt; hatchling backend
```

Data flow for PKG-01/PKG-03:

```
uv pip install -e .
  → hatchling reads pyproject.toml
  → discovers petringa/ package (packages = ["petringa"])
  → registers petringa/ on sys.path (editable)

pytest -q
  → test_*.py: from petringa.core.X import Y
  → [fast] unit tests: no devsim required
  → [slow] @pytest.mark.slow: devsim called directly; v3_frozen.json compared
```

### Recommended Project Structure

```
petringa/               ← package root (renamed from src/)
├── __init__.py         ← minimal: exports DeviceConfig stub for SC#2
├── _version.py         ← version = "5.0.0"
└── core/               ← all existing src/ modules
    ├── __init__.py     ← empty (internal namespace)
    ├── sic_material.py
    ├── device.py
    ├── device2d.py
    ├── poisson.py
    ├── drift_diffusion.py
    ├── cv_analysis.py
    ├── charge_collection.py
    ├── charge_collection_2d.py
    ├── generation_profiles.py
    ├── single_particle.py
    ├── mc_coupling.py
    ├── microdosimetry.py
    ├── radiation_damage.py
    ├── dark_current.py
    ├── flash_recombination.py
    ├── temperature_sweep.py
    ├── transient.py
    ├── incomplete_ionization.py
    ├── alternative_structures.py
    ├── devsim_reset.py
    ├── optimization.py
    ├── validation.py
    ├── analytical.py
    ├── plotting.py
    └── plotting2d.py
```

### Pattern 1: pyproject.toml with hatchling — Explicit Package Selection

**What:** Declare the `petringa/` package explicitly so hatchling does not accidentally include `app/`, `tests/`, `scripts/`, `notebooks/`, `data/` in the wheel.
**When to use:** Whenever the project root contains sibling directories that are NOT part of the package.
**Example:**

```toml
# Source: https://hatch.pypa.io/latest/config/build/
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "petringa"
version = "5.0.0"
description = "SiC TCAD Simulator Library for the Petringa group"
requires-python = ">=3.13"
dependencies = [
    "devsim>=2.10.0",
    "numpy>=1.24",
    "scipy>=1.11",
    "matplotlib>=3.7",
    "plotly>=5.0",
    "streamlit>=1.30",
    "pandas>=2.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.0",
    "jupyter>=1.0",
]

[tool.hatch.build.targets.wheel]
packages = ["petringa"]
```

**Critical:** Without `packages = ["petringa"]`, hatchling auto-discovers all Python packages in the project root, which may include unintended directories. [CITED: https://hatch.pypa.io/latest/config/build/]

### Pattern 2: Editable Install for Development

**What:** `uv pip install -e .` registers the package in editable mode so changes to `petringa/` take effect immediately without reinstalling.
**Example:**

```bash
# Install package with dev extras
uv pip install -e ".[dev]"

# Verify
python -c "from petringa import DeviceConfig"
python -c "import petringa; print(petringa.__version__)"
```

The editable install makes `petringa/` the importable namespace; `from petringa.core.device import ...` resolves to `./petringa/core/device.py` directly.

### Pattern 3: Mechanical Import Rewrite — sed-based

**What:** Transform all `from src.X` → `from petringa.core.X` and `import src.X` → `import petringa.core.X` across all Python source files.
**Why not relative imports:** Switching from `from src.X import Y` to `from .X import Y` is a larger semantic change, makes the code harder to grep, and offers no benefit for this refactor phase.
**Example script:**

```bash
# Source files, tests, scripts — mechanical substitution
find /path/to/project -name '*.py' \
  ! -path '*/.venv/*' ! -path '*/.git/*' \
  -exec sed -i '' \
    -e 's/from src\.\([^ ]*\)/from petringa.core.\1/g' \
    -e 's/import src\.\([^ ]*\)/import petringa.core.\1/g' \
  {} +
```

**String-literal refs require manual edits (6 total — see Runtime State Inventory).**

**Notebooks (.ipynb) require a separate mechanism — sed skips them silently.** The sed pass above covers all `.py` files, including the string literals embedded in `create_notebook_*.py` scripts. But the 22 `.ipynb` files are JSON, not plain text, so `find -name '*.py'` skips them entirely. A Python json-walker is required:

```python
# Rewrite from src. imports in all .ipynb notebook cells
import json, re, pathlib

for nb_path in pathlib.Path("notebooks").glob("*.ipynb"):
    nb = json.loads(nb_path.read_text())
    changed = False
    for cell in nb.get("cells", []):
        new_src = []
        for line in cell.get("source", []):
            new_line = re.sub(r"from src\.([^ ]+)", r"from petringa.core.\1", line)
            new_line = re.sub(r"import src\.([^ ]+)", r"import petringa.core.\1", new_line)
            if new_line != line:
                changed = True
            new_src.append(new_line)
        cell["source"] = new_src
    if changed:
        nb_path.write_text(json.dumps(nb, indent=1))

# Verify zero residual src. imports in notebooks
import subprocess, sys
result = subprocess.run(
    ["grep", "-r", "from src.", "notebooks/"],
    capture_output=True, text=True
)
if result.stdout:
    print("ERROR: residual src. imports found:", result.stdout)
    sys.exit(1)
```

Post-rewrite verification: `grep -r 'from src\.' notebooks/` must return zero results.

### Pattern 4: Minimal DeviceConfig Stub for SC#2

**What:** A minimal dataclass in `petringa/__init__.py` that satisfies `python -c "from petringa import DeviceConfig"` without implementing the full API (which belongs to Phase 36).
**When to use:** Phase 35 — the full implementation is Phase 36 (LIB-01).
**Example:**

```python
# petringa/__init__.py  — Phase 35 skeleton
"""petringa: SiC TCAD Simulator Library."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

__version__ = "5.0.0"


@dataclass
class DeviceConfig:
    """Device configuration stub — full implementation in Phase 36."""
    epi_thickness_um: float = 10.0
    substrate_thickness_um: float = 1.0
    half_width_um: Optional[float] = None
    N_A: float = 1e19
    doping_profile: str = "graded"
    N_D: Optional[float] = None
    N_D_junction: float = 2.93e15
    N_D_bulk: float = 8.82e13
    L_transition_um: float = 0.987
    T: float = 300.0
    area_cm2: float = 1e-4


__all__ = ["DeviceConfig", "__version__"]
```

This stub is intentionally minimal — it carries the right field names per the design spec (§3.1) so Phase 36 can fill in `__post_init__`, validation, and `build_device()` without changing the public interface.

### Anti-Patterns to Avoid

- **Adding `src/` back to sys.path in conftest.py:** The old approach was that pytest found `src.` imports because the project root was on `sys.path`. After the rename and `uv pip install -e .`, sys.path gets `petringa/` via the editable install. Adding a conftest.py that re-inserts `src/` would mask rename failures silently.
- **Switching to relative imports:** `from .device import X` inside `petringa/core/device.py` works, but changes 91 import lines in source modules to a different semantic form. Unnecessary for a pure rename.
- **`hatchling` auto-discovery without `packages = [...]`:** Without explicit package declaration, hatchling may vacuum up `tests/`, `scripts/`, `notebooks/` as sub-packages if they contain `__init__.py` files. Set `[tool.hatch.build.targets.wheel] packages = ["petringa"]` explicitly.
- **Using `pip install -e .` instead of `uv pip install -e .`:** Project memory mandates uv; all install steps in the plan must use uv.

---

## Don't Hand-Roll

| Problem           | Don't Build            | Use Instead                          | Why                                                                                                                                                         |
| ----------------- | ---------------------- | ------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Import rewriting  | Custom AST transformer | sed/Python string replace            | The rewrite is purely textual — `from src.X` → `from petringa.core.X`; regex handles it exactly; AST transform adds risk of touching non-import occurrences |
| Package discovery | Custom finder          | hatchling `packages` key             | PEP 517/518; hatchling handles editable installs, wheel builds, and metadata automatically                                                                  |
| Version string    | Duplicated literals    | `_version.py` imported by `__init__` | Single source of truth; hatchling can read it dynamically or use static value                                                                               |

**Key insight:** This phase is a mechanical rename, not a build system from scratch. The only hand-work needed is the six string/path/logger references that a text substitution cannot safely handle.

---

## Runtime State Inventory

> Phase 35 is a rename/refactor trigger. Each category answered explicitly.

| Category                                                             | Items Found                                                                                                                                                                    | Action Required                                                                                                                                                                                                                                                               |
| -------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Stored data                                                          | None — no datastore keys on Python module names                                                                                                                                | None                                                                                                                                                                                                                                                                          |
| Live service config                                                  | None — no external services reference `src.` module names                                                                                                                      | None                                                                                                                                                                                                                                                                          |
| OS-registered state                                                  | None — no launchd/cron/pm2 entries reference `src.*`                                                                                                                           | None                                                                                                                                                                                                                                                                          |
| Secrets/env vars                                                     | None — no `.env` files or CI env vars reference `src.*`                                                                                                                        | None                                                                                                                                                                                                                                                                          |
| Build artifacts                                                      | `.venv/lib/python3.13/site-packages/` — no `petringa.egg-info` or `src.egg-info` or `.pth` file for `src/` found; existing venv has devsim+deps but NOT the `petringa` package | After install: `uv pip install -e .` creates editable link for `petringa/`; no stale artifact cleanup needed                                                                                                                                                                  |
| **Code references** — `from src.X` / `import src.X`                  | **326 import statements** across: src/ (91), tests/ (176), scripts/ (59)                                                                                                       | Mechanical textual rewrite: `from src.` → `from petringa.core.`, `import src.` → `import petringa.core.`                                                                                                                                                                      |
| **String-literal src refs** (6 total — not caught by import rewrite) | 1. `tests/test_mc_coupling.py:301` — `patch("src.single_particle.ion_track_generation_2d")`                                                                                    | Change to `patch("petringa.core.single_particle.ion_track_generation_2d")`                                                                                                                                                                                                    |
|                                                                      | 2. `tests/test_mc_coupling.py:330` — same patch target                                                                                                                         | Same fix                                                                                                                                                                                                                                                                      |
|                                                                      | 3. `tests/test_radiation_damage.py:452` — `module_path = "src/radiation_damage.py"` (AST-parsed for structural test)                                                           | Change to `"petringa/core/radiation_damage.py"`                                                                                                                                                                                                                               |
|                                                                      | 4. `scripts/run_calibration_2d.py:51` — `DEVICE2D_PATH = pathlib.Path("src/device2d.py")`                                                                                      | Change to `pathlib.Path("petringa/core/device2d.py")`                                                                                                                                                                                                                         |
|                                                                      | 5. `scripts/create_notebook_16.py:110` — `logging.getLogger('src.single_particle')`                                                                                            | Change to `logging.getLogger('petringa.core.single_particle')`                                                                                                                                                                                                                |
|                                                                      | 6. `scripts/create_notebook_20.py:189` — `logging.getLogger('src.optimization')`                                                                                               | Change to `logging.getLogger('petringa.core.optimization')`                                                                                                                                                                                                                   |
| **String-literal src refs in scripts (notebook-building code)**      | scripts/create_notebook_03.py, 04, 05, 08, 15, 15_v2, 16, 17, 18, 19, 20 embed `from src.X` as string literals inside Python source strings they write to .ipynb cells         | These scripts generate notebook cells; fix the string content so generated notebooks get `from petringa.core.X`                                                                                                                                                               |
| **Notebook imports** (22 notebooks, 81 import statements)            | All 22 .ipynb files have cells with `from src.X import Y`                                                                                                                      | Phase 35 scope decision: rewrite notebook imports now (cheap, mechanical) while noting their full execution is validated in Phase 43 — see Open Questions                                                                                                                     |
| **sys.path manipulation in scripts**                                 | `scripts/freeze_v3_baselines.py`, `scripts/diagnose_1d_2d_parity.py`, `scripts/run_calibration_2d.py` — insert project root into sys.path so `import src.*` resolves           | After rename + editable install, the sys.path insert is still needed for scripts run as `uv run python scripts/...`; the project root stays valid since `petringa/` lives there. No action needed on the sys.path lines themselves, only on the `from src.` lines below them. |

**Canonical question — after all files updated, what runtime systems still have the old string cached/stored/registered?**
Answer: None. All references are in version-controlled Python source files. No external service, OS registry, or datastore keys on `src.*` module names.

---

## Common Pitfalls

### Pitfall 1: mock.patch string target not updated

**What goes wrong:** `patch("src.single_particle.ion_track_generation_2d")` uses a dotted module path string. After the rename, the patch target resolves to the old (now non-existent) path, so the mock silently does not apply — tests pass for the wrong reason (no mock) or fail with `ModuleNotFoundError`.
**Why it happens:** `mock.patch` takes a string module path, not an imported name; text search for `from src.` misses it.
**How to avoid:** Explicitly grep for `patch("src.` and `patch('src.` as part of the rename checklist — already identified (2 occurrences in `tests/test_mc_coupling.py`).
**Warning signs:** Tests that previously exercised mock paths now test real code paths; `assert mock_itg.call_count == 1` assertions fail.

### Pitfall 2: `pytest -q` includes slow devsim integration tests

**What goes wrong:** SC#3 says "all 25 test modules pass." `pytest.ini` registers the `slow` marker but does NOT exclude it by default. Running `pytest -q` without `-m "not slow"` runs the full suite including `@pytest.mark.slow` tests in `test_v3_baseline_regression.py`. These call devsim, build 2D devices, and compare against `v3_frozen.json`. This means PKG-03 acceptance requires devsim **fully functional** (BLAS loaded, solver working), not just importable.
**Why it happens:** Confusing marker registration (suppresses warning) with marker deselection (actually skips tests).
**How to avoid:** Run acceptance suite with bare `pytest -q` (no `-m` flag), verifying the devsim-backed slow tests pass. The devsim BLAS warning (`libopenblas.dylib: MISSING DLL` / `liblapack.dylib: ALL BLAS/LAPACK LOADED`) is expected on macOS arm64 and does not indicate failure as long as LAPACK loads.
**Warning signs:** All tests except slow ones pass — this is NOT a passing suite for PKG-03.

### Pitfall 3: hatchling auto-discovering unintended packages

**What goes wrong:** If `[tool.hatch.build.targets.wheel] packages = ["petringa"]` is omitted, hatchling scans for all Python packages in the project root. If `tests/`, `scripts/`, `notebooks/`, or `data/` directories happen to have `__init__.py` files, they get included in the wheel and pollute the installed package.
**Why it happens:** Hatchling's default auto-discovery is greedy.
**How to avoid:** Always specify `packages = ["petringa"]` explicitly in the wheel target.
**Warning signs:** `pip show petringa` lists unexpected files; `python -c "import tests"` succeeds after install.

### Pitfall 4: Script sys.path manipulation survives rename but masks failures

**What goes wrong:** Scripts like `freeze_v3_baselines.py` do `sys.path.insert(0, project_root)` before `from src.X import Y`. After rename, `project_root` is still valid (petringa/ lives there), BUT if the `from src.X` lines are not updated, Python finds neither `src/` (renamed) nor the new path, throwing `ModuleNotFoundError`. However, if only the sys.path line is removed and lines not updated, failure mode is the same.
**How to avoid:** Update the `from src.X` lines in scripts to `from petringa.core.X`; the sys.path insert becomes redundant after `uv pip install -e .` but can remain for scripts that pre-date the install.

### Pitfall 5: Notebook cell string-embedded `from src.` imports

**What goes wrong:** Scripts like `create_notebook_03.py` embed `"from src.generation_profiles import (\n"` as Python string literals that get written to .ipynb cells. A mechanical import rewrite (`sed` on `from src.`) will correctly rewrite the Python module itself (e.g., the outer `.py` file's imports), but the **embedded string content** also contains `from src.` and must be updated too, or the generated notebooks will have broken imports.
**How to avoid:** The sed command replacing `from src.` must also process string content inside `create_notebook_*.py` files. Since these are Python files, the same sed pass handles both Python code and string literals. But verify after rewrite by grepping for residual `from src.` in scripts.

---

## Code Examples

Verified patterns from official sources:

### Complete pyproject.toml for petringa

```toml
# Source: https://hatch.pypa.io/latest/config/build/ + https://hatch.pypa.io/latest/config/metadata/
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "petringa"
version = "5.0.0"
description = "SiC TCAD Simulator Library"
requires-python = ">=3.13"
dependencies = [
    "devsim>=2.10.0",
    "numpy>=1.24",
    "scipy>=1.11",
    "matplotlib>=3.7",
    "plotly>=5.0",
    "streamlit>=1.30",
    "pandas>=2.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.0",
    "jupyter>=1.0",
]

[tool.hatch.build.targets.wheel]
packages = ["petringa"]
```

### Import rewrite: before and after

```python
# BEFORE (src/device.py and all callers)
from src.poisson import setup_poisson, solve_equilibrium
from src.sic_material import SiC4H_Parameters

# AFTER (petringa/core/device.py and all callers)
from petringa.core.poisson import setup_poisson, solve_equilibrium
from petringa.core.sic_material import SiC4H_Parameters
```

### mock.patch string target: before and after

```python
# BEFORE (tests/test_mc_coupling.py:301)
with patch("src.single_particle.ion_track_generation_2d") as mock_itg:

# AFTER
with patch("petringa.core.single_particle.ion_track_generation_2d") as mock_itg:
```

### AST structural test path: before and after

```python
# BEFORE (tests/test_radiation_damage.py:452)
module_path = "src/radiation_damage.py"

# AFTER
module_path = "petringa/core/radiation_damage.py"
```

### Editable install smoke test sequence

```bash
uv pip install -e ".[dev]"
python -c "from petringa import DeviceConfig; print('OK')"
python -c "from petringa.core.device import create_sic_device; print('core OK')"
pytest -q
```

---

## Open Questions

1. **DeviceConfig scope boundary between Phase 35 and Phase 36**
   - What we know: Success Criterion #2 for Phase 35 is `python -c "from petringa import DeviceConfig"`. REQUIREMENTS.md maps LIB-01 (full `DeviceConfig` dataclass) to Phase 36. The design spec (§4, Phase 1) bundles both the rename and `DeviceConfig` creation together.
   - What's unclear: Does Phase 35 need only a minimal skeleton dataclass (fields but no validation, no `build_device()`) to satisfy the smoke test, with Phase 36 filling in the full implementation? Or should Phase 35 ship nothing (breaking SC#2) and Phase 36 deliver the full dataclass?
   - Recommendation: Ship a minimal `DeviceConfig` skeleton in `petringa/__init__.py` during Phase 35 — just the `@dataclass` with the fields from the design spec §3.1 (no physics, no devsim calls). This satisfies the explicit SC#2 smoke test and is < 20 lines. Phase 36 upgrades it with validation, `build_device()`, and the full API. The planner should make this boundary explicit and coordinate with Phase 36.

2. **Notebook import rewrite: in or out of Phase 35 scope?**
   - What we know: PKG-03 covers test modules only; the 22 notebooks' import correctness is verified in Phase 43. There are 81 `from src.X` imports across 22 notebooks.
   - What's unclear: Rewriting notebook imports in Phase 35 is cheap (same mechanical pass) but leaves notebooks broken if deferred. If left, Phase 43 must do a notebook import rewrite sweep in addition to execution verification.
   - Recommendation: Include notebook import rewrite in Phase 35 as a separate task (it's cheap, mechanical, and doesn't gate the acceptance test). Mark notebook execution as Phase 43 scope. This avoids the 22 notebooks being in a broken-import state for the entire v5.0 window (Phases 35-43).

3. **`devsim` entry in pyproject.toml dependencies — installability on a fresh machine**
   - What we know: devsim is a C extension (prebuilt wheel for macOS arm64). It installs via `pip install devsim` but it is not a PyPI package in the traditional sense — it may require a custom index or direct URL.
   - What's unclear: Whether `uv pip install -e .` on a fresh machine with `devsim>=2.10.0` in dependencies will successfully locate and install the devsim wheel, or whether devsim requires manual handling (BLAS/LAPACK prereqs noted in venv output).
   - Recommendation: The planner should add a note that `devsim` in `dependencies` is a declaration for documentation/tooling purposes, and that `pip install -e .` from a fresh venv may require devsim to be pre-installed separately (devsim install verification is already a project convention). This doesn't affect Phase 35's acceptance test (devsim is already present in the project venv).

---

## Environment Availability

| Dependency | Required By                | Available        | Version                                                | Fallback                                                                      |
| ---------- | -------------------------- | ---------------- | ------------------------------------------------------ | ----------------------------------------------------------------------------- |
| Python     | All                        | ✓                | 3.13.5                                                 | —                                                                             |
| uv         | All installs               | ✓                | 0.9.8                                                  | —                                                                             |
| devsim     | PKG-03 slow tests          | ✓                | 2.10.0 (BLAS: LAPACK ok, libopenblas missing on macOS) | No fallback — slow tests require devsim                                       |
| pytest     | PKG-03                     | ✓                | 9.0.2                                                  | —                                                                             |
| numpy      | tests + slow tests         | ✓                | 2.4.3                                                  | —                                                                             |
| scipy      | tests                      | ✓                | 1.17.1                                                 | —                                                                             |
| matplotlib | tests                      | ✓                | 3.10.8                                                 | —                                                                             |
| pandas     | tests                      | ✓                | 3.0.1                                                  | —                                                                             |
| plotly     | pyproject.toml declaration | ✗ (not in .venv) | —                                                      | Declare in pyproject.toml; install will pull it                               |
| streamlit  | pyproject.toml declaration | ✗ (not in .venv) | —                                                      | Declare in pyproject.toml; install will pull it                               |
| hatchling  | Build backend              | ✗ (not in .venv) | 1.30.1 on PyPI                                         | `uv pip install -e .` pulls hatchling automatically via build-system.requires |

**Missing dependencies with no fallback:**

- None that block Phase 35 execution (devsim is present)

**Missing dependencies with fallback:**

- plotly, streamlit: not yet in .venv; will be pulled on `uv pip install -e ".[dev]"` or later when UI phases need them. Phase 35 only declares them in pyproject.toml; no test exercises them.
- hatchling: auto-pulled as part of the build backend during `uv pip install -e .`; no pre-install needed.

**devsim BLAS note:** Running `python -c "import devsim"` on macOS arm64 emits `Loading "libopenblas.dylib": MISSING DLL` followed by `Loading "liblapack.dylib": ALL BLAS/LAPACK LOADED`. This is normal — LAPACK loads successfully via the fallback; devsim is fully functional. The PKG-03 slow tests will run correctly.

---

## Validation Architecture

### Test Framework

| Property           | Value                                                 |
| ------------------ | ----------------------------------------------------- |
| Framework          | pytest 9.0.2                                          |
| Config file        | `pytest.ini` (project root)                           |
| Quick run command  | `pytest -q -m "not slow"`                             |
| Full suite command | `pytest -q` (includes @pytest.mark.slow devsim tests) |

**Critical note:** The PKG-03 acceptance gate is `pytest -q` (full suite, including slow tests). The quick run command excludes slow tests for day-to-day development. The Phase 35 final verification must run the full suite.

### Phase Requirements → Test Map

| Req ID   | Behavior                                        | Test Type   | Automated Command                                                                       | File Exists?                            |
| -------- | ----------------------------------------------- | ----------- | --------------------------------------------------------------------------------------- | --------------------------------------- |
| PKG-01   | `uv pip install -e .` completes without error   | smoke       | `uv pip install -e . && python -c "import petringa"`                                    | ❌ Wave 0 — add to plan as install step |
| PKG-02   | pyproject.toml declares all deps and dev extras | smoke       | `python -c "import importlib.metadata; print(importlib.metadata.requires('petringa'))"` | ❌ Wave 0                               |
| PKG-03   | All 25 test modules pass                        | integration | `pytest -q`                                                                             | ✅ 25 existing test\_\*.py files        |
| SC#2     | `from petringa import DeviceConfig` works       | smoke       | `python -c "from petringa import DeviceConfig; print('OK')"`                            | ❌ Wave 0 — need petringa/**init**.py   |
| Baseline | v3_frozen.json byte-for-byte unchanged          | regression  | `git diff --exit-code tests/baselines/v3_frozen.json` (VCS check — see note)            | ✅ file in git                          |

**Baseline integrity note:** `test_v3_baseline_regression.py` reads `v3_frozen.json` and asserts that the *current code* reproduces the frozen values within tolerance. It does NOT assert the JSON file itself is unmodified — if `freeze_v3_baselines.py` were re-run (which also gets its imports rewritten in this phase), it would overwrite the baseline and the test would still pass with new numbers. The only correct way to assert byte-for-byte baseline integrity is: `git diff --exit-code tests/baselines/v3_frozen.json`. This must be a distinct phase-gate step, separate from `pytest -q`.

### Sampling Rate

- **Per task commit:** `pytest -q -m "not slow"` (fast tests only; < 30 s)
- **Per wave merge:** `pytest -q` (full suite including slow)
- **Phase gate:** Full suite green (`pytest -q`) + `python -c "from petringa import DeviceConfig"` + `python -c "from petringa.core.device import create_sic_device"` + `git diff --exit-code tests/baselines/v3_frozen.json` (confirms baseline was not regenerated) before verification

### Wave 0 Gaps

- [ ] `petringa/__init__.py` — stub DeviceConfig + `__version__`; covers SC#2 and PKG smoke
- [ ] `petringa/core/__init__.py` — empty; required for package namespace
- [ ] `petringa/_version.py` — version = "5.0.0"
- [ ] `pyproject.toml` — replaces requirements.txt; covers PKG-01, PKG-02
- [ ] `pytest.ini` update — add `testpaths = tests` if needed after rename (currently no testpaths directive; pytest finds tests/ by convention)

_(Existing 25 test\__.py files cover PKG-03 once imports are rewritten — no new test files needed)\*

---

## Security Domain

> This phase has no new attack surface: it is a pure rename + packaging step with no network calls, no new endpoints, and no new user-facing functionality. ASVS categories are not applicable. `security_enforcement` is not set in `.planning/config.json`, but the phase does not warrant a security domain section beyond this note.

---

## State of the Art

| Old Approach                     | Current Approach                           | When Changed                                            | Impact                                                   |
| -------------------------------- | ------------------------------------------ | ------------------------------------------------------- | -------------------------------------------------------- |
| `requirements.txt` + manual venv | `pyproject.toml` + `uv pip install -e .`   | PEP 517/518 standardized ~2018; hatchling emerged ~2022 | Proper editable installs, metadata introspection, extras |
| `setup.py` / `setup.cfg`         | `pyproject.toml`                           | Deprecated in favor of pyproject.toml                   | Simpler, no imperative code at build time                |
| `python setup.py develop`        | `pip install -e .` / `uv pip install -e .` | PEP 660 (2021)                                          | Standards-compliant editable installs without `setup.py` |

**Deprecated/outdated:**

- `requirements.txt` alone: Still valid for pinned CI deps, but not a substitute for package metadata; does not enable `pip install -e .`
- `setup.py develop`: Removed from modern pip; use `pip install -e .` or `uv pip install -e .`

---

## Assumptions Log

| #   | Claim                                                                                          | Section                      | Risk if Wrong                                                                                                                                                 |
| --- | ---------------------------------------------------------------------------------------------- | ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| A1  | plotly, streamlit will install cleanly with `uv pip install -e .` on macOS arm64 Python 3.13   | Standard Stack / Environment | Low risk — both have arm64 wheels; worst case is version incompatibility requiring a pin                                                                      |
| A2  | devsim `>=2.10.0` in pyproject.toml dependencies is installable from PyPI on a fresh machine   | Open Questions #3            | Medium risk — if devsim requires a custom index, the `[project] dependencies` entry may fail on other machines; acceptance test on existing venv not affected |
| A3  | hatchling will be pulled automatically as part of `uv pip install -e .` without pre-installing | Build system                 | LOW — this is standard PEP 517 behavior; uv handles build-system.requires automatically                                                                       |
| A4  | plotly and streamlit package names on PyPI are exactly `plotly` and `streamlit`                | Standard Stack               | LOW — these are unambiguous long-established packages; risk is near-zero                                                                                      |

**If this table is empty:** Not the case — A2 (devsim on fresh machine) is worth planner attention.

---

## Sources

### Primary (HIGH confidence)

- `https://hatch.pypa.io/latest/config/build/` — hatchling wheel target configuration, `packages` key, editable install behavior [CITED]
- `https://hatch.pypa.io/latest/config/metadata/` — `[project]` table, `dependencies`, `optional-dependencies` structure [CITED]
- `pip3 index versions hatchling` — confirmed hatchling 1.30.1 is current on PyPI [VERIFIED: PyPI via pip3 index]
- `.venv/lib/python3.13/site-packages/` directory scan — confirmed versions of devsim, numpy, scipy, matplotlib, pandas, pytest, jupyter [VERIFIED: filesystem]
- `uv --version` — confirmed uv 0.9.8 available [VERIFIED: CLI]
- `grep` scan of all Python files — 326 `from src.` / `import src.` statements; 6 string-literal src refs [VERIFIED: grep on filesystem]
- `.planning/STATE.md` — locked decisions (package name, build backend, rename target) [VERIFIED: project file]
- `docs/superpowers/specs/2026-06-26-simulator-library-ui-design.md` — DeviceConfig field names, project structure, Phase 1 scope [VERIFIED: project file]

### Secondary (MEDIUM confidence)

- `tests/test_v3_baseline_regression.py` — confirms slow tests call devsim; `pytest -q` runs them [VERIFIED: code read]
- `pytest.ini` — confirms `slow` marker registered but NOT auto-excluded [VERIFIED: file read]
- `tests/test_mc_coupling.py`, `tests/test_radiation_damage.py` — confirmed string-literal `src.` refs in mock.patch and AST path [VERIFIED: grep + code read]

### Tertiary (LOW confidence)

- None — all claims in this research are verified or cited.

---

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH — versions confirmed in .venv; hatchling from PyPI
- Architecture: HIGH — design spec, STATE.md, and code structure all agree
- Import scope: HIGH — grep-verified counts (326 import statements; 6 string refs; 81 notebook imports)
- Pitfalls: HIGH — two confirmed from code inspection (mock.patch, slow test inclusion); others from hatchling docs

**Research date:** 2026-06-26
**Valid until:** 2026-09-26 (stable packaging ecosystem; hatchling API is stable)
