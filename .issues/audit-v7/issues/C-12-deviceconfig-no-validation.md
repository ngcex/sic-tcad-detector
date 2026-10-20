---
id: C-12
source: raw-C F12
title: DeviceConfig has no validation; an unknown doping_profile silently becomes uniform, and uniform with N_D=None silently uses a hidden default
severity: low
category: bug
status: fixed (d7280f7; numeric validation cf64bef)
location:
  - etna/api/device.py:25-43 (DeviceConfig, no __post_init__)
  - etna/core/device.py:92-93 (`if N_D is None: N_D = DEFAULT_N_D`, DEFAULT_N_D = 1.07e15 at :33)
  - etna/core/device.py:228-247 (`if doping_profile == "graded": ... else: uniform`)
labels: [api, input-validation, public-contract]
---

## Summary

- `DeviceConfig(doping_profile="Graded")`, `"step"` or any typo builds a uniform device with no error.
- `DeviceConfig(doping_profile="uniform")` with N_D left at None builds with the hidden `DEFAULT_N_D = 1.07e15`.
- Geometry, doping, T and area can be zero, negative or non-finite.

DeviceConfig is the public contract and is serialized into every CSV header.

Line-number correction: the raw file cites core/device.py:51-52 and :187-206. The current
locations are :92-93 and :228-247.

Scope note: the app always supplies N_D when "uniform" is selected (device_sidebar.py:100-103)
and restricts the profile to two strings. The N_D=None and typo cases therefore affect library
and notebook users only. The zero/negative cases affect both (C-03, C-09).

## Reproduction / Evidence

```
uv run python -c "from etna.core.device import create_sic_device; \
  d = create_sic_device(device_name='x1', doping_profile='Graded'); print('built ok', d['doping_profile'])"
-> built ok with Graded: Graded          # no error; took the uniform else-branch
grep -n "__post_init__\|raise" etna/api/device.py  -> (no output)
```

## Impact

Typos produce a physically different device that is labeled with the user's intended name.
Missing N_D silently substitutes a value the user never chose.

## Proposed fix

```diff
 @dataclass
 class DeviceConfig:
     ...
+    def __post_init__(self):
+        if self.doping_profile not in ("graded", "uniform"):
+            raise ValueError(f"doping_profile must be 'graded' or 'uniform', got {self.doping_profile!r}")
+        if self.doping_profile == "uniform" and self.N_D is None:
+            raise ValueError("doping_profile='uniform' requires N_D")
+        for name in ("epi_thickness_um", "substrate_thickness_um", "N_A", "N_D_junction",
+                     "N_D_bulk", "T", "area_cm2"):
+            v = getattr(self, name)
+            if not (math.isfinite(v) and v > 0):
+                raise ValueError(f"{name} must be finite and > 0, got {v!r}")
+        if not (math.isfinite(self.L_transition_um) and self.L_transition_um >= 0):
+            raise ValueError("L_transition_um must be finite and >= 0")
+        if self.half_width_um is not None and not self.half_width_um > 0:
+            raise ValueError("half_width_um must be > 0 or None")
```

`dataclasses.replace` (used by ParametricSweep) re-runs `__post_init__`, so sweeps are covered
automatically.

## Test plan

A parametrized `pytest.raises(ValueError)` case for each rule. Also check that `DeviceConfig()`
and the sidebar's `assemble_config` outputs still construct, and that the existing
test_api_* / test_app_* suites pass.

## Effort

S

## Status 2026-10-18 (final-c): open (partial)

- Fixed in cf64bef: DeviceConfig.__post_init__ rejects non-positive / NaN geometry, doping, area, N_t, out-of-gap E_t, T outside 100-600 K, N_A < 1e17, donors >= N_A, inconsistent buffer / donor_table settings; dataclasses.replace re-validates.
- Still open (repro on final-c): `DeviceConfig(doping_profile="Graded")` and `DeviceConfig(doping_profile="uniform")` (N_D None) both construct; device_kwargs passes them on, and the uniform path substitutes DEFAULT_N_D = 1.07e15 (etna/api/simulation.py validated_domain_reason / compensation_fluence, etna/core/device.py).
- Why not fixed here: the __post_init__ docstring says doping_profile is "deliberately not validated here (downstream code decides)", and two tests pin that acceptance (tests/test_app_robustness.py::test_device_config_accepts_defaults_and_edges, tests/test_api_engine_switch.py::test_outside_validated_domain_uses_devsim with doping_profile="linear"). Rejecting them changes the public contract.
- Owner: Cesare (API owner). Decision needed: reject unknown doping_profile and require N_D for "uniform" in __post_init__ (then update the two pinned tests), or keep the current contract and document the DEFAULT_N_D fallback.

## Status 2026-10-21: fixed (d7280f7)

Contract decision taken with the API owner (Cesare): strict, no case folding.

- `doping_profile` is typed `DopingProfile = Literal["graded", "uniform"]`
  (`etna.api.device.DOPING_PROFILES`). Any other value, `"Graded"` included, raises
  ValueError listing the allowed values. There is no case folding, because the config is
  serialised into every CSV / figure header and its provenance must be unique per physics.
- `doping_profile="uniform"` with `N_D=None` raises ValueError naming N_D, so the hidden
  DEFAULT_N_D 1.07e15 is no longer reachable from DeviceConfig. Core `create_sic_device`
  keeps its own fallback (not the public contract).
- `doping_profile="graded"` with an explicit `N_D` raises. Checked first: the graded
  builders set the epi from N_D_junction / N_D_bulk / L_transition (core/device.py graded
  branch), `mobility_doping` uses N_D_bulk, the batched engine uses `nd` only where
  `uniform` is true (simulation.py `np.where(uniform, nd, ...)`). The only readers of
  `device_info["N_D"]` are the analytic `extract_depletion_width` (no facade calls it) and
  the 2D-engine gate, which treats any N_D as out of box (it would only switch solver).
  So the value was silently ignored, and it is now rejected.
- Presets: `legacy_config()` is graded with an explicit junction / bulk / L triplet and no
  N_D, so it does not rely on a hidden default and needed no change. The microdosimeter
  (the only uniform preset) already sets N_D = 1e14.
- Positive checks now also reject inf (and N_A must be finite). The checks are factored
  into `DeviceConfig.validate()` (pure), which `__post_init__` and `device_kwargs` both call
  (see C-03).
- Call sites: no notebook uses DeviceConfig. The app restricts the profile to the two
  strings and forces N_D=None for graded. Scripts and tests all pass N_D with uniform.
  Two pinned tests were updated: tests/test_app_robustness.py ("linear" now raises) and
  tests/test_api_engine_switch.py (out-of-box case is now a uniform N_D 5e13).
- Bit-identity: a sha256 fingerprint of x / y for every `run_all(DeviceConfig())` section
  and of `run_cv` for the legacy / microdosimeter / n+-calibrated / tabulated presets is
  identical before (main f3f00ee) and after. Pinned in tests/test_device_config_contract.py.
- Follow-up (not in this change, simulation.py is owned elsewhere): the
  `DEFAULT_N_D if config.N_D is None else config.N_D` fallbacks in
  etna/api/simulation.py (Debye length, engine batch, compensation fluence) and
  etna/export/labels.py are now dead code for the uniform branch.
