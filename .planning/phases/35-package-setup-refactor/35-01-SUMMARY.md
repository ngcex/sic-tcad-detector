---
phase: 35-package-setup-refactor
plan: 01
subsystem: infra
tags: [python-packaging, hatchling, pyproject, uv, dataclass]

# Dependency graph
requires:
  - phase: 25
    provides: prior milestone completion baseline (v3.0 shipped, 20 notebooks)
provides:
  - "Installable petringa package registered via editable install"
  - "pyproject.toml single source of truth for build metadata and dependencies"
  - "DeviceConfig stub dataclass importable from petringa"
  - "petringa/_version.py single-source version string"
affects:
  [
    36-core-api-deviceconfig-cv-field,
    37-core-api-cce-facades,
    all-subsequent-v5.0-phases,
  ]

# Tech tracking
tech-stack:
  added: [hatchling, plotly, streamlit, pandas]
  patterns:
    - "petringa/_version.py as single-source version string, re-exported via petringa/__init__.py"
    - "DeviceConfig stub dataclass (field-per-line, inline unit comments) — matches src/sic_material.py convention"
    - "[tool.hatch.build.targets.wheel] packages explicit list to avoid tests/ auto-discovery leak"

key-files:
  created:
    - pyproject.toml
    - petringa/__init__.py
    - petringa/_version.py
  modified:
    - README.md

key-decisions:
  - "petringa/core/ deliberately NOT created in this plan — reserved for Plan 02's git mv src petringa/core to avoid nesting bug"
  - "requirements.txt removed via git rm — pyproject.toml is now the sole dependency manifest per PKG-01"
  - "pandas>=2.0 added to runtime deps (was previously an undeclared dependency used by test_mc_coupling.py)"

patterns-established:
  - "Single-source version: petringa/_version.py defines __version__, __init__.py re-exports it"
  - 'Explicit packages = ["petringa"] in [tool.hatch.build.targets.wheel] to prevent tests/scripts/notebooks leaking into the wheel'

requirements-completed: [PKG-01, PKG-02]

# Metrics
duration: 12min
completed: 2026-07-01
---

# Phase 35 Plan 01: Package Scaffold Summary

**Installable `petringa` package via hatchling with pyproject.toml (7 runtime deps + dev extras), DeviceConfig stub dataclass, and single-source version string — editable install verified working.**

## Performance

- **Duration:** 12 min
- **Started:** 2026-07-01T18:06:45Z (approx, per STATE.md session start)
- **Completed:** 2026-07-01
- **Tasks:** 1
- **Files modified:** 4 (3 created, 1 modified) + 1 removed

## Accomplishments

- `pyproject.toml` created with hatchling build backend, all 7 runtime dependencies (devsim, numpy, scipy, matplotlib, plotly, streamlit, pandas), and `[dev]` optional extras (pytest, jupyter)
- `petringa/__init__.py` created with `DeviceConfig` stub dataclass (11 fields matching design spec §3.1) and `__version__` re-export
- `petringa/_version.py` created as single source of truth for version string (5.0.0)
- `requirements.txt` removed (git rm) — replaced by pyproject.toml per PKG-01
- README.md install instructions updated to `uv pip install -e ".[dev]"`
- Editable install completed successfully (`uv pip install -e ".[dev]"`); `from petringa import DeviceConfig` and `petringa.__version__` both verified working

## Task Commits

Each task was committed atomically:

1. **Task 1: Create pyproject.toml and petringa package scaffold** - `b426987` (feat)

**Plan metadata:** (this commit, to follow)

## Files Created/Modified

- `pyproject.toml` - Build metadata: hatchling backend, project deps, dev extras, explicit wheel packages list
- `petringa/__init__.py` - DeviceConfig stub dataclass (11 fields) + `__version__` re-export + `__all__`
- `petringa/_version.py` - `__version__ = "5.0.0"` single-source version string
- `README.md` - Install instruction updated from `pip install -r requirements.txt` to `uv pip install -e ".[dev]"`
- `requirements.txt` - Removed (git rm), superseded by pyproject.toml

## Decisions Made

- Did not create `petringa/core/` in this plan — deliberately deferred to Plan 02, which relies on `git mv src petringa/core` creating the directory fresh (creating it here first would cause `git mv` to nest `src/` inside `petringa/core/` instead of renaming)
- `pandas>=2.0` included in runtime deps even though not in the old `requirements.txt`, because `tests/test_mc_coupling.py` imports it at top level — this was an undeclared dependency now made explicit
- Used `git rm` (not manual delete) for `requirements.txt` to keep the removal tracked in git history

## Deviations from Plan

None - plan executed exactly as written. A local formatter hook auto-wrapped two long lines in `petringa/__init__.py` (the `Optional[float] = None` field defaults) after the Write tool call; this is a whitespace-only reformatting with no semantic change and was verified not to break any acceptance criteria.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required. `uv pip install -e ".[dev]"` pulled `plotly`, `streamlit`, and their transitive dependencies from PyPI automatically (per threat model T-35-SC, disposition: accept).

## Next Phase Readiness

- Plan 02 (same phase) can now safely run `git mv src petringa/core` since `petringa/core/` does not exist yet
- `petringa` package is installed in editable mode in the active `.venv`; `DeviceConfig` stub and `__version__` are importable
- Phase 36 (Core API) can build on this scaffold — `DeviceConfig` stub is ready for `__post_init__` validation and `build_device()` additions without changing the public interface

---

## Self-Check: PASSED

- FOUND: `/Users/ngcex/projects/physics/petringa/pyproject.toml`
- FOUND: `/Users/ngcex/projects/physics/petringa/petringa/__init__.py`
- FOUND: `/Users/ngcex/projects/physics/petringa/petringa/_version.py`
- MISSING (expected, by design): `/Users/ngcex/projects/physics/petringa/petringa/core/` — deliberately deferred to Plan 02
- MISSING (expected, by design): `/Users/ngcex/projects/physics/petringa/requirements.txt` — removed per PKG-01
- FOUND commit: `b426987` (feat(35-01): scaffold installable petringa package with hatchling)
- Verified: `from petringa import DeviceConfig` exits 0
- Verified: `petringa.__version__` == "5.0.0"
- Verified: `grep 'packages.*petringa' pyproject.toml` matches

---

_Phase: 35-package-setup-refactor_
_Completed: 2026-07-01_
