---
id: D-08
title: Nine etna.core public functions are referenced only by their own tests
severity: low
category: dead-code
status: fixed 2026-10-18 (wp-c2 910e2a6 + 834e1ae: 3 of the 9 removed - compute_cce_from_current, compute_microdosimetric_means, plot_electric_field_multi; the other 6 kept with a one-line reason in their docstrings)
location: etna/core/charge_collection.py:145 compute_cce_from_current; :1225 cce_anneal_vs_temperature; etna/core/cv_analysis.py:66 depletion_width_from_capacitance; etna/core/dark_current.py:920 dark_current_post_anneal; etna/core/mc_coupling.py:166 list_root_trees; :314 events_to_charge_profiles; etna/core/microdosimetry.py:226 compute_microdosimetric_means; etna/core/plotting.py:76 plot_electric_field_multi; etna/core/sic_material.py:308 effective_dos
labels: [dead-code, api-surface]
source: raw-D F8
---

## Summary

A word-boundary grep over `etna/ app/ scripts/ examples/ notebooks/` finds no call site for any of the
nine symbols. The only non-test hits are self-references in log strings or docstrings
(`cce_anneal_vs_temperature`, `dark_current_post_anneal`, `events_to_charge_profiles`). Each symbol is
referenced by exactly one test file. `etna/__init__.py` exports only the facades.
Severity is downgraded to low (raw: medium). Unused library surface has no runtime effect, and some of
these (`effective_dos`, `depletion_width_from_capacitance`) are reasonable utility API.

## Reproduction/Evidence

```
$ for s in ...; do grep -rnw "$s" etna/ app/ scripts/ examples/ notebooks/ | grep -v "def $s" | wc -l; grep -rlw "$s" tests/ | wc -l; done
compute_cce_from_current nontest=0 testfiles=1
cce_anneal_vs_temperature nontest=1 (log f-string, charge_collection.py:1298) testfiles=1
depletion_width_from_capacitance nontest=0 testfiles=1
dark_current_post_anneal nontest=2 (log strings, dark_current.py:1090,1110) testfiles=1
list_root_trees nontest=0 testfiles=1
events_to_charge_profiles nontest=1 (module docstring, mc_coupling.py:13) testfiles=1
compute_microdosimetric_means nontest=0 testfiles=1
plot_electric_field_multi nontest=0 testfiles=1
effective_dos nontest=0 testfiles=1
```

I did not re-run vulture. The raw report's claim that `junction_offset` in `plot_electric_field_multi` is
an unread parameter was not re-checked.

## Impact

Maintenance burden and test time spent on code nothing uses. Some of these functions (for example
`compute_cce_from_current`) may encode older formulas, which risks someone picking them up by mistake.

## Proposed fix

Triage each one. Delete those that are superseded (`compute_cce_from_current`, which CCE now gets from
`compute_cce_from_dd`, plus `events_to_charge_profiles` and `cce_anneal_vs_temperature` if they are
superseded) along with their tests. Keep the genuine utilities and mark them as library API in their
docstrings. Optionally add a vulture CI step with an allowlist.

## Test plan

After deletion, `uv run pytest tests/ -q` is green and grep finds no dangling references.

## Effort

M
