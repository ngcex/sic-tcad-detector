---
id: B-12
title: Dead code and test-only symbols in the physics-model slice (dead assignment, unused import/guard/variable, unwired anneal chain, stale docstring)
severity: low
category: dead-code
status: fixed 2026-10-18 (wp-c2: ce2e236 import warnings; 910e2a6 compute_cce_from_current, compute_microdosimetric_means, plot_electric_field_multi removed, list_root_trees kept; 834e1ae anneal chain / events_to_charge_profiles kept as public API with reasons, doping overridable via device_kwargs and read from graded_stack_kwargs (a78cbdf); dd84a64 _R_sic_mm, None guard, citation, validate_2d_vs_1d reason; the hecht_cce_partial_depletion defaults text was already corrected; mu_e dead assignment already gone)
location: etna/core/{charge_collection,dark_current,radiation_damage,mc_coupling,microdosimetry,plotting,plotting2d,generation_profiles,optimization}.py (per item below)
labels: [dead-code, cleanup, raw-F12]
---

## Summary

This is one finding with several items. Each was checked with `grep -rln <name> --include='*.py' --include='*.ipynb' etna app tests scripts notebooks`. "Test-only" means the only hits are the defining module and its test file.

- [x] **charge_collection.py:105.** `mu_e = mobility_caughey_thomas_T(...)` is overwritten immediately at :108-112, so the assignment is dead. That line is the only use of the import at :31, so the import becomes unused once the line goes.
- [x] **`charge_collection.compute_cce_from_current`** is test-only (charge_collection.py and test_charge_collection.py).
- [x] **`cce_post_anneal`, `cce_anneal_vs_temperature`** (charge_collection) and **`dark_current_post_anneal`** are test-only. As a result, `radiation_damage.compute_annealed_params`, `defect_recovery_fractions`, `annealing_fraction` and `AnnealingParams` have no consumer outside tests (hits limited to radiation_damage, charge_collection, dark_current and test_radiation_damage). `cce_post_anneal` and its sibling also hardcode `N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1.0e-4` in internal calls (charge_collection.py:1090-1092, 1136-1138) instead of taking them as parameters.
- [x] **`mc_coupling.list_root_trees`, `mc_coupling.events_to_charge_profiles`** are test-only.
- [x] **`microdosimetry.compute_microdosimetric_means`** is test-only. It duplicates the y_F/y_D code in `lineal_energy_spectrum`.
- [x] **`plotting.plot_electric_field_multi`** is test-only and also crashes for one voltage (B-09).
- [x] **`plotting2d.validate_2d_vs_1d`** is used only by tests and the superseded `scripts/create_notebook_15.py`. It does not appear in `notebooks/15_2d_electrostatics_cce.ipynb`.
- [x] **generation_profiles.py `_R_sic_mm`** is a dead variable kept with `noqa: F841` (B-10).
- [x] **optimization.py:148-156.** The `raw_center_cce is None` guard cannot be reached. `compute_cce_2d` always returns `float` (`return 0.0` or `return float(np.clip(...))`, charge_collection_2d.py:197,200), and `cce_lateral_scan` collects those values into an array.
- [x] **radiation_damage.py:1061.** `import warnings` inside `radiation_hardness_sweep` is unused. `grep -n "warnings" etna/core/radiation_damage.py` finds only that line.
- [x] **`charge_collection.hecht_cce_partial_depletion` docstring** lists defaults "950 / 1e-9 / 125 / 6e-7", but the code resolves `_UNSET` from `SiC4H_Parameters`. It also cites "Sciencedirect S0168900205006443" without author or title.

## Reproduction / Evidence

- ran: the grep loop over all names above. Output:
  - `compute_cce_from_current`: charge_collection.py, test_charge_collection.py
  - `cce_post_anneal`, `cce_anneal_vs_temperature`: same two files
  - `dark_current_post_anneal`: dark_current.py, test_dark_current.py
  - `defect_recovery_fractions`, `annealing_fraction`: radiation_damage.py, test_radiation_damage.py
  - `list_root_trees`, `events_to_charge_profiles`: mc_coupling.py, test_mc_coupling.py
  - `compute_microdosimetric_means`: microdosimetry.py, test_microdosimetry.py
  - `plot_electric_field_multi`: plotting.py, test_plotting.py
  - `validate_2d_vs_1d`: plotting2d.py, test_plotting2d.py, scripts/create_notebook_15.py
- ran: `sed -n 100,114p etna/core/charge_collection.py` shows the dead assignment. `grep -n mobility_caughey_thomas_T etna/core/charge_collection.py` shows only :31 (import) and :105.
- ran: `grep -n "def hecht_cce_partial_depletion" -A40` shows the stale defaults at docstring lines 199-208.

## Impact

Maintenance burden. The anneal chain looks like supported functionality but is not reachable from the API, app or notebooks, and it hardcodes doping. The stale docstring misleads readers about the values in use.

## Proposed fix

- Delete the dead assignment at charge_collection.py:105 and the now-unused import, the `_R_sic_mm` line, the `import warnings` at radiation_damage.py:1061, and the unreachable `None` guard (keep only the float branch).
- Anneal chain: either wire it into the API or a notebook, with `N_D_junction`, `N_D_bulk` and `L_transition` as parameters mirroring `cce_vs_fluence`, or move it to an `etna/core/experimental/annealing.py` module with a clear "not validated / not wired" docstring.
- Remove `compute_cce_from_current`, `compute_microdosimetric_means` (or make `lineal_energy_spectrum` call it), `list_root_trees` and `validate_2d_vs_1d`, or document them as public utilities.
- Fix the `hecht_cce_partial_depletion` docstring: "Default: from SiC4H_Parameters (T-dependent)", plus a full citation.
- Add ruff `F841`/`F401` (without the noqa) and, optionally, `vulture` to CI.

## Test plan

`uv run ruff check etna --select F401,F841` is clean. The full test suite stays green after deletions, minus the tests of deleted functions. If the anneal chain is kept, add a test that passes non-default doping and observes a different CCE.

## Effort

M
