---
id: D-09
title: 45 ruff F-rule findings (unused imports, dead assignments incl. leftover I_generated, unasserted W_60)
severity: low
category: dead-code
status: partial  # 2026-10-18 sync: 7acbe10 (wp-c1): the two scripts/ findings (diagnose numpy, regression_sweep traceback) fixed; etna/ and tests/ findings are WP-C2
location: etna/core/{device.py:23,140; dark_current.py:42,633; charge_collection.py:964,1353,1397; optimization.py:25; radiation_damage.py:1061; temperature_sweep.py:21,24; alternative_structures.py:60,346; device2d.py:318; single_particle.py:156-157; flash_recombination.py:259}; scripts/diagnose_1d_2d_parity.py:48; scripts/regression_sweep_v3_notebooks.py:24; tests/* (26 findings)
labels: [lint, dead-code, ruff]
source: raw-D F9
---

## Summary

`ruff check --select F401,F811,F841,F821,F632,F541` reports exactly 45 findings, 28 of them auto-fixable.
There are no F821 (undefined-name) errors, so nothing here crashes at runtime. The findings worth noting:

- `charge_collection.py:964,1397`: `I_generated` is computed and never used, left over from an older CCE formula.
- `tests/test_poisson.py:164`: `W_60` is computed and never asserted, so an assertion is missing.
- `tests/test_app_run_mockability.py:18/45`: F811 redefinition of `DeviceConfig`.

Small correction to the raw report: `single_particle.py:156-157` is F841 (unused `device`/`region`), not
an f-string issue. The F541 finding is `flash_recombination.py:259` only.

## Reproduction/Evidence

```
$ uvx ruff check --select F401,F811,F841,F821,F632,F541 --output-format concise etna scripts tests examples app
... (full list matches the location field)
Found 45 errors.
[*] 28 fixable with `--fix` (17 unsafe fixes available with `--fix --unsafe-fixes`).
```

## Impact

Noise. The unused `I_generated` and `W_60` hint at logic that is incomplete or stale.

## Proposed fix

```
uvx ruff check --fix --select F401,F541 etna scripts tests
```

Then remove the F841 assignments by hand (check that each RHS has no side effects first, especially
devsim calls). Add to `tests/test_poisson.py` around line 164: `assert W_60 < epi_thickness` (or
`W_60 <= W_200`). Add a `ruff check --select F` step to CI or pre-commit.

## Test plan

The ruff command returns "All checks passed". `uv run pytest tests/test_poisson.py -q` is green with the new assertion.

## Effort

S
