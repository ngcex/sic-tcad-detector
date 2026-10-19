---
id: B-07
title: run_radiation_damage default proton_energy_MeV=5.6 is outside the NIEL table (30-150 MeV) and silently clamps to kappa(30)
severity: medium
category: bug
status: open 2026-10-18 (final-verify: physics default, not changed; 91194a5 only logs the edge-kappa use once and documents it in run_radiation_damage; owner Giada Petringa / INFN-LNS, see B-INDEX 'Physics items left open')
location: etna/api/simulation.py:451,478 (default, docstring); etna/core/radiation_damage.py:61-66 (table), 407-429 (get_hardness_factor), ~663 (kappa = get_hardness_factor(energy_MeV)); etna/core/generation_profiles.py:235-243 (raises outside range)
labels: [radiation-damage, api, niel, extrapolation, raw-F7]
---

## Summary

`get_hardness_factor` returns `np.interp(energy_MeV, energies, factors)`. Outside the [30, 150] MeV table, `np.interp` silently clamps to the endpoint values. The public API `run_radiation_damage(..., proton_energy_MeV=5.6)` therefore uses kappa(30 MeV) = 0.5 while reporting `energy_MeV=5.6` in its metadata. Proton NIEL in SiC rises steeply below about 30 MeV, so the default API run understates damage per proton, and the label suggests a result specific to 5.6 MeV. `generation_profiles._interpolate_proton_range_water` raises `ValueError` outside 30-150 MeV, so the two modules behave inconsistently. The rest of the codebase defaults to 62 MeV (`cce_vs_fluence`, `compute_damaged_params`).

## Reproduction / Evidence

- ran: `uv run python -c "from etna.core.radiation_damage import get_hardness_factor as g; print(g(5.6), g(30), g(300), g(150))"` prints `0.5 0.5 0.22 0.22`.
- ran: `grep -rn "5\.6" etna/api app` finds only `simulation.py:451` (default) and `:478` (docstring "Default 5.6."). Neither gives a rationale.
- call chain: simulation.py forwards `energy_MeV=proton_energy_MeV`. `cce_vs_fluence` calls `compute_damaged_params(energy_MeV=...)`, which sets `kappa = get_hardness_factor(energy_MeV)`.
- The Streamlit page restricts the choice to [30, 62, 70, 150] MeV (from the raw audit, `app/workflows/radiation_damage.py:84`), so only API and notebook callers get the default.
- The kappa table values themselves are already-known placeholders and are not re-reported here.

## Impact

API users who keep the default get the damage of a 30 MeV proton labelled as 5.6 MeV. Any extrapolation is silent.

## Proposed fix

```diff
 def get_hardness_factor(energy_MeV, niel_table=None):
     ...
     energies = sorted(niel_table.keys())
+    if not (energies[0] <= energy_MeV <= energies[-1]):
+        raise ValueError(
+            f"energy_MeV={energy_MeV} outside NIEL table range "
+            f"[{energies[0]}, {energies[-1]}] MeV; no extrapolation"
+        )
```

```diff
-    proton_energy_MeV: float = 5.6,
+    proton_energy_MeV: float = 62.0,
 ...
-        Default 5.6.
+        Default 62.0 (must lie within the NIEL table, 30-150 MeV).
```

## Test plan

- `pytest.raises(ValueError)` for `get_hardness_factor(5.6)` and `get_hardness_factor(200)`.
- Test that the default `run_radiation_damage` energy lies within `NIEL_HARDNESS_PROTON_SIC` (inspect the signature; no solve needed).

## Effort

S
