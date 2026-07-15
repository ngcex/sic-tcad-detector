# Phase 35: Package Setup & Refactor - Pattern Map

**Mapped:** 2026-06-26
**Files analyzed:** ~55 (3 new + 25 core modules + 25 tests + scripts + notebooks)
**Analogs found:** 2 / 3 new files (pyproject.toml has no codebase analog)

---

## Context: Two-Bucket Structure

Phase 35 is a rename/refactor — not a feature phase. The deliverable splits into:

- **Bucket A — Genuinely new files (3):** `pyproject.toml`, `etna/__init__.py`, `etna/core/__init__.py` / `etna/_version.py`. These need copy-from patterns.
- **Bucket B — Moved/rewritten files (~55):** All `src/*.py` modules, all `tests/*.py`, `scripts/*.py`, `notebooks/*.ipynb`. The "pattern" is the textual transformation rule, not a per-file analog. One worked exemplar covers all.

---

## File Classification

| New/Modified File                  | Role         | Data Flow                    | Closest Analog                              | Match Quality                     |
| ---------------------------------- | ------------ | ---------------------------- | ------------------------------------------- | --------------------------------- |
| `pyproject.toml`                   | config       | N/A                          | `requirements.txt` (same dep list, partial) | no-analog (content from RESEARCH) |
| `etna/__init__.py`             | package-root | N/A                          | `src/sic_material.py` (dataclass style)     | style-match only                  |
| `etna/_version.py`             | config       | N/A                          | none                                        | no-analog (trivial one-liner)     |
| `etna/core/__init__.py`        | namespace    | N/A                          | `src/__init__.py` (empty)                   | exact (both empty)                |
| `etna/core/*.py` (25 modules)  | transform    | import-rewrite               | `src/*.py` (same file, renamed)             | exact (same file)                 |
| `tests/*.py` (25 test files)       | transform    | import-rewrite               | `tests/*.py` (same file)                    | exact (same file)                 |
| `scripts/*.py` (15 scripts)        | transform    | import-rewrite               | `scripts/*.py` (same file)                  | exact (same file)                 |
| `notebooks/*.ipynb` (22 notebooks) | transform    | import-rewrite (json-walker) | `notebooks/*.ipynb` (same file)             | exact (same file)                 |

---

## Pattern Assignments

### Bucket A: New Files

---

### `pyproject.toml` (config)

**Analog:** `requirements.txt` (partial — use as dep-list cross-check only)

**requirements.txt contents** (current, `/Users/ngcex/projects/physics/etna/requirements.txt`):

```
devsim>=2.10.0
numpy>=1.24
scipy>=1.11
matplotlib>=3.7
pytest>=7.0
```

**Delta vs pyproject.toml:** pyproject.toml adds `plotly>=5.0`, `streamlit>=1.30`, `pandas>=2.0` to runtime deps (not in requirements.txt), moves `pytest` to `[dev]` optional extra, and adds `jupyter>=1.0` to `[dev]`. The `pandas` addition is not cosmetic — `test_mc_coupling.py` imports `pandas` at top-level, meaning the runtime dep was undeclared in requirements.txt.

**Full pattern** (from RESEARCH §Pattern 1 + §Code Examples):

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "etna"
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
packages = ["etna"]
```

**Critical:** `packages = ["etna"]` must be explicit. Without it, hatchling auto-discovery may include `tests/`, `scripts/`, `notebooks/` if they contain `__init__.py`. `tests/__init__.py` exists (confirmed by `ls tests/`), so this is a real risk.

---

### `etna/__init__.py` (package-root, DeviceConfig stub)

**Analog:** `src/sic_material.py` — this codebase's dataclass style convention

**Dataclass style pattern** (`/Users/ngcex/projects/physics/etna/src/sic_material.py`, lines 12-30):

```python
from dataclasses import dataclass
import numpy as np


@dataclass
class SiC4H_Parameters:
    """4H-SiC material parameters at 300K.

    All values sourced from literature with citations.
    Units are CGS (cm, cm^-3, eV, s, F/cm) for devsim compatibility.
    """

    # --- Bandgap ---
    E_g: float = 3.26  # eV at 300K, Ioffe NSM Archive
    E_g_0: float = 3.2965625  # eV at 0K, calibrated so Varshni gives E_g=3.26 at 300K
    E_g_alpha: float = 6.5e-4  # eV/K, Varshni parameter (Ioffe NSM)
    E_g_beta: float = 1300.0  # K, Varshni parameter (Ioffe NSM)

    # --- Dielectric constant ---
```

**Pattern conventions to copy:**

- `@dataclass` with field-per-line and inline comment (units/source)
- No `__post_init__` in stub — that is Phase 36 scope
- `from dataclasses import dataclass` (not `dataclasses.dataclass`)

**Full stub content** (from RESEARCH §Pattern 4):

```python
# etna/__init__.py  — Phase 35 skeleton
"""etna: SiC TCAD Simulator Library."""

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

**Phase boundary note:** This stub satisfies SC#2 (`python -c "from etna import DeviceConfig"`). Phase 36 (LIB-01) adds `__post_init__` validation and `build_device()` without changing the public interface. Planner must mark this explicitly in the Phase 35 plan so Phase 36 knows it is upgrading, not creating.

---

### `etna/_version.py` (config)

**Analog:** None — trivial one-liner.

```python
__version__ = "5.0.0"
```

`etna/__init__.py` imports this as: `from etna._version import __version__` (or simply define `__version__` inline in `__init__.py` directly — either approach works; the RESEARCH recommendation is `_version.py` as single source of truth).

---

### `etna/core/__init__.py` (namespace)

**Analog:** `src/__init__.py` — confirmed empty (1-line file, blank).

Copy exactly: empty file (or single-line docstring `"""etna.core: internal simulation modules."""`). No exports — internal namespace only per RESEARCH §Decisions.

---

## Bucket B: Transformation Pattern (Applies to All ~55 Moved/Rewritten Files)

### Cross-Cutting Import Rewrite Rule

**Source of truth:** `/Users/ngcex/projects/physics/etna/src/charge_collection_2d.py`, lines 32-39 (densest cluster of internal cross-module imports in the codebase):

```python
# BEFORE (src/charge_collection_2d.py lines 32-39)
from src.device2d import create_sic_2d_device
from src.poisson import setup_poisson, solve_equilibrium
from src.drift_diffusion import (
    setup_sic_drift_diffusion,
    extract_contact_current,
    ramp_bias,
)
from src.charge_collection import add_generation_to_dd
```

```python
# AFTER (etna/core/charge_collection_2d.py — same lines, rewritten)
from etna.core.device2d import create_sic_2d_device
from etna.core.poisson import setup_poisson, solve_equilibrium
from etna.core.drift_diffusion import (
    setup_sic_drift_diffusion,
    extract_contact_current,
    ramp_bias,
)
from etna.core.charge_collection import add_generation_to_dd
```

**Rule:** `from src.X` → `from etna.core.X` and `import src.X` → `import etna.core.X`. Do NOT switch to relative imports (`from .X import`); the diff is larger and grep becomes harder.

**Mechanical rewrite command (`.py` files only):**

```bash
find /Users/ngcex/projects/physics/etna -name '*.py' \
  ! -path '*/.venv/*' ! -path '*/.git/*' \
  -exec sed -i '' \
    -e 's/from src\.\([^ ]*\)/from etna.core.\1/g' \
    -e 's/import src\.\([^ ]*\)/import etna.core.\1/g' \
  {} +
```

**Notebook rewrite (`.ipynb` files — JSON, sed skips them):**

```python
import json, re, pathlib

for nb_path in pathlib.Path("notebooks").glob("*.ipynb"):
    nb = json.loads(nb_path.read_text())
    changed = False
    for cell in nb.get("cells", []):
        new_src = []
        for line in cell.get("source", []):
            new_line = re.sub(r"from src\.([^ ]+)", r"from etna.core.\1", line)
            new_line = re.sub(r"import src\.([^ ]+)", r"import etna.core.\1", new_line)
            if new_line != line:
                changed = True
            new_src.append(new_line)
        cell["source"] = new_src
    if changed:
        nb_path.write_text(json.dumps(nb, indent=1))
```

**Verification after rewrite:**

```bash
grep -r 'from src\.' /Users/ngcex/projects/physics/etna/src/ \
  /Users/ngcex/projects/physics/etna/tests/ \
  /Users/ngcex/projects/physics/etna/scripts/ \
  /Users/ngcex/projects/physics/etna/notebooks/
# Must return zero results
```

---

## Shared Patterns

### The 6 String-Literal Fixes (NOT caught by mechanical rewrite — highest risk)

These are hard-coded string values that reference `src.*` module paths, not Python import statements. The sed pass and json-walker above will NOT fix them. Each requires a manual targeted edit.

| #   | File                             | Line | Current value                                          | Required value                                                   |
| --- | -------------------------------- | ---- | ------------------------------------------------------ | ---------------------------------------------------------------- |
| 1   | `tests/test_mc_coupling.py`      | 301  | `patch("src.single_particle.ion_track_generation_2d")` | `patch("etna.core.single_particle.ion_track_generation_2d")` |
| 2   | `tests/test_mc_coupling.py`      | 330  | `patch("src.single_particle.ion_track_generation_2d")` | `patch("etna.core.single_particle.ion_track_generation_2d")` |
| 3   | `tests/test_radiation_damage.py` | 452  | `module_path = "src/radiation_damage.py"`              | `module_path = "etna/core/radiation_damage.py"`              |
| 4   | `scripts/run_calibration_2d.py`  | 51   | `DEVICE2D_PATH = pathlib.Path("src/device2d.py")`      | `pathlib.Path("etna/core/device2d.py")`                      |
| 5   | `scripts/create_notebook_16.py`  | 110  | `logging.getLogger('src.single_particle')`             | `logging.getLogger('etna.core.single_particle')`             |
| 6   | `scripts/create_notebook_20.py`  | 189  | `logging.getLogger('src.optimization')`                | `logging.getLogger('etna.core.optimization')`                |

**Why items 1-2 are the highest risk:** `mock.patch` resolves the string as a dotted module path at runtime. If not updated, the mock silently does not apply — the test may pass for the wrong reason (no mock) or fail with `ModuleNotFoundError`. `assert mock_itg.call_count == 1` assertions in surrounding test logic will give false results.

### Notebook-Building Scripts: Embedded `from src.` in String Literals

Scripts `create_notebook_03.py`, `04`, `05`, `08`, `15`, `15_v2`, `16`, `17`, `18`, `19`, `20` embed `from src.X import Y` as Python string literals that get written into notebook cells. The mechanical sed pass handles these correctly because the scripts are `.py` files and sed processes all text in the file, including string content. No special handling required — but verify after the sed pass:

```bash
grep -n 'from src\.' /Users/ngcex/projects/physics/etna/scripts/create_notebook_*.py
# Must return zero results after sed pass
```

Confirmed example before sed (`scripts/create_notebook_03.py` line 38):

```python
"from src.generation_profiles import (\n"
```

After sed, this becomes:

```python
"from etna.core.generation_profiles import (\n"
```

### sys.path Manipulation in Scripts (No Action Required)

Scripts `freeze_v3_baselines.py`, `diagnose_1d_2d_parity.py`, `run_calibration_2d.py` insert the project root into `sys.path` before importing. After rename + `uv pip install -e .`, the project root is still the correct path (because `etna/` lives there). The `sys.path.insert` lines themselves do not need to change — only the `from src.X` import lines below them.

**Pattern from** `scripts/run_calibration_2d.py` lines 39-43:

```python
_PROJECT_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from src.device2d import calibrate_graded_doping_2d  # <- only this line changes
```

After rewrite, the import becomes `from etna.core.device2d import calibrate_graded_doping_2d`. The `sys.path.insert` block is retained as-is (harmless after editable install).

### pytest.ini (No Change Required)

Current `/Users/ngcex/projects/physics/etna/pytest.ini`:

```ini
[pytest]
markers =
    slow: marks tests as slow (devsim integration tests, >10s each)
```

No `testpaths` or `pythonpath` directive is present. pytest discovers `tests/` by convention. After editable install, `from etna.core.X` imports resolve via the installed package — no `pythonpath` addition needed. pytest.ini stays unchanged.

**However:** The planner should confirm that `[tool.pytest.ini_options]` testpaths are not needed. If the planner moves config to `pyproject.toml`, the equivalent section is:

```toml
[tool.pytest.ini_options]
markers = ["slow: marks tests as slow (devsim integration tests, >10s each)"]
```

---

## No Analog Found

| File                   | Role   | Reason                                                                                                  |
| ---------------------- | ------ | ------------------------------------------------------------------------------------------------------- |
| `pyproject.toml`       | config | No existing pyproject.toml in project; content from RESEARCH §Code Examples and official hatchling docs |
| `etna/_version.py` | config | Trivial one-liner; no analog needed                                                                     |

---

## Metadata

**Analog search scope:** `src/`, `tests/`, `scripts/`, `notebooks/`
**Files scanned:** 25 src modules, 25 test files, 15 scripts, 22 notebooks, `requirements.txt`, `pytest.ini`, `src/__init__.py`, `tests/__init__.py`
**Pattern extraction date:** 2026-06-26

**Scope counts from RESEARCH (grep-verified):**

- 326 `from src.` / `import src.` statements across: src/ (91), tests/ (176), scripts/ (59)
- 6 string-literal `src.*` refs requiring manual fixes (identified above)
- 81 `from src.X` imports in 22 notebooks (json-walker handles)
- 11 `create_notebook_*.py` scripts with embedded `from src.` string literals (sed pass handles)
