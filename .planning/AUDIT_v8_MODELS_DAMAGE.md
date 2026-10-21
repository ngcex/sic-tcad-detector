# Audit v8: physical models and radiation damage (2026-10-21)

Independent adversarial review of code written in the previous day and never second-read:
`physics_models.py`, `avalanche.py`, `batched/models_dd.py`, `devsim_models.py`,
`sic_material.py` (c-axis mobility default), `radiation_traps.py`, `trap_kinetics.py`,
`trap_presets.py`, `trap_calibration.py`, `batched/poisson.py` (trap kernels),
`radiation_damage.py`, `charge_collection.py` (damage paths), `api/simulation.py`
(damage dispatch, `"auto"`, NIEL conventions), `niel_tables.py`, `neutron_niel.py`, and
the USAGE.md NIEL-convention algebra. Base: main at 7746e8d. No code under `etna/` was
changed.

Every major finding and the numerical minor ones (N1-N4, N8) were reproduced
numerically. N5-N7 and N9-N11 are text, label and source checks. Each candidate was
also attacked, and the refuted ones are listed at the end. Repro scripts for the major findings and some minor
ones are in `lab/audit_v8/` (run with `uv run python lab/audit_v8/<script>.py`).

**Severity rubric.**

- Critical: wrong physics or numbers on the default path (`"auto"`, `sic-damage`, c-axis mobility).
- Major: wrong results on an opt-in path, a silent physics-changing fallback, or a
  default-path error that is real but small.
- Minor: docs, tests, labels, provenance.

**Totals: 0 critical, 4 major, 12 minor.** One of the minor items, N12, is an open
provenance question demoted from major after review (issue #14).

## Verdict

The core physics is sound. Each of the following was re-derived or tested:

- SRH occupancy (single and two-level negative-U) is exactly Fermi-Dirac at equilibrium (2e-16).
- The charge assignment, backward-Euler steps, 2x2 kinetic system and trap Jacobians
  all check out (FD to 1e-4 / 1.5e-3).
- The Hatakeyama 2004/2005 impact-ionisation coefficients and the anisotropic formula
  match the source.
- The McIntyre integrals, the DD signs and Jacobian, and the mobility sets agree with
  Ayalew and the 2410.06798 review tables.
- All Burin Table I, Capan Table I and Radmanovac Table 2 numbers used in the code match
  the opened PDFs.
- All documented NIEL / compensation numbers reproduce.

The NIEL convention / calibrated `d` algebra is consistent, and on the trap path
nothing is double counted (verdict below). The defects are in the plumbing:

- the convention guard has holes (M1, M4);
- the default lifetime average depends on the discretisation (M2);
- the provenance of the legacy scalar `eta_removal` is ambiguous: the code anchors it per
  neq, a review proposes re-citing it per proton (N12, open question);
- the tilted-track plasma slices ignore `DeviceConfig.mobility` (M5).

## NIEL convention vs calibrated removal d: independent re-derivation

Definitions, with N the NIEL:

    kappa_si-neq(E)     = N_p,Si(E)  / N_n,Si(1 MeV, ASTM 95 mb)
    kappa_sic-damage(E) = N_p,SiC(E) / (r_n N_n,Si(1 MeV)),  r_n = <N_SiC>_fission / <N_Si>_fission = 1.0687
    kappa_sic-1mev(E)   = N_p,SiC(E) / N_n,SiC(1 MeV)

With these, a fission neutron has the same kappa in si-neq and sic-damage. Rates anchored
on reactor neq (Burin Table I) therefore carry over unchanged. Per proton they change by
r_p/r_n, which is the intended correction.

Proton-calibrated data are per proton: g_p = kappa_c [G_lit(kappa_c Phi) + d_c]. G_lit
is not linear in kappa, so d_c has to be re-fitted per convention, not rescaled.

Reproduced (scripts of sub-audit C):

- kappa at 62 MeV: 1.942 / 1.552 / 1.470.
- kappa at 252.7 MeV: 0.8514 / 0.6390 / 0.6051.
- r_n/r_p: 1.133 / 1.251 / 1.303 at 5.6 / 62 / 150 MeV.
- Phi_comp at 62 MeV: 1.802e12 / 4.36e12 / 3.48e12 / 1.696e12.
- Rescaled d: 3.83.
- The per-proton slopes agree between conventions to 0.1 /p.

Single conversion points:

- kappa is applied once (`compute_damaged_params`, `from_config`, `calibrate_cv`).
- The trap path sets `carrier_removal=False`, so the scalar removal is not stacked on the traps.
- `donor_removal` is applied once per engine.
- `resolve_traps("auto", conv)` picks d = 4.69 / 2.87 / 5.10 correctly.
- The 12.09 /neq is derived from the preset's own levels:
  - Z1+Z2 contribute 2 x 2.5 = 5.0, EH4 2.4 and d 4.69.
  - The full two-level neutral solution gives 12.090 /neq up to 0.95 Phi_c.

**Verdict:** consistent; nothing double counted on the trap path. Exceptions:

- the guard holes M1 and M4 let a preset run with the wrong d;
- whether the legacy `eta_removal` should carry over unchanged depends on its anchor. The
  code labels it Burin, 1-MeV-neutron-equivalent, so it carries over by item 1. If the
  review's proposed per-proton re-citation is adopted, it must be re-fitted per
  convention (N12).

Open, unverified: Burin Table I gives "g_int [cm^-1]" and "1 MeV equivalent neutron
fluences" with no Si/SiC qualifier (PDF text read). The "Si-labelled" in the
sic-damage definition is therefore the community convention, not something the source
states.

## Major

### M1. `KineticTrapSolver.from_config` runs a preset given by name in the wrong NIEL convention (silent) (issue #12)

- **Where:** `etna/core/trap_kinetics.py:445-454` and `:364`.
  - `from_config` checks `resolve_traps(trap_model, conv)` but passes the raw string to `cls(...)`.
  - `__init__` then calls `resolve_traps(trap_model)` with no convention, which rebuilds the
    preset in the default `sic-damage` convention.
  - kappa still comes from `conv`.
- **Evidence:** `lab/audit_v8/m1_from_config_convention.py`. `from_config(cfg, [1e12], 252.7,
trap_model="radmanovac_calibrated", niel_convention="si-neq")` gives
  `model.niel_convention == "sic-damage"` and d = 4.69, against the expected 2.87
  (+63 %). Passing the model object raises, as it should.
- **Related:** `trap_calibration.py:148` (`calibrate_cv`) has the same missing convention.
  There it fails loudly: a spurious ValueError for `"auto"` under si-neq.
- **Fix:** `return cls(ref, F*kappa, resolve_traps(trap_model, conv))`. Use
  `resolve_traps(trap_model, conv)` in `calibrate_cv`.
- **Confidence:** high. Found independently by sub-audits B and C.

### M2. Default-path CCE lifetimes are a node-count mean over a mask that includes the n+ buffer and substrate (issue #13)

- **Where:** `etna/core/radiation_traps.py:514-517` (`occupancy_lifetimes`: `(...)[m].mean()`),
  called from `charge_collection.py:1040-1042` with `ref["epi_mask"]` (x >= junction).
  The default preset has `lifetimes="occupancy"`, so this is the default damage path.
- **What is wrong:**
  - The mean is unweighted, while the mesh spacing runs from 1.5e-4 to 0.099 um. The dense
    junction region and the N_D > 1e16 region dominate: 26 % of the nodes for 11 % of the
    length, where the acceptors are filled and cannot capture.
  - The docstring says "(the epi)".
- **Evidence:** `lab/audit_v8/m2_occupancy_lifetime_average.py` (default config, -40 V, 62 MeV).
  Trap part of 1/tau_n. The "same mask" column isolates the node-count weighting from
  the mask choice; averaging over the whole epi, neutral region included, is disclosed
  in USAGE, but the node weighting is not:

| Phi_p (p/cm²) | code, node mean (1/s) | length-weighted, same mask (1/s) | weighting alone | length-weighted, epi only (1/s) | epi / code |
| --- | --- | --- | --- | --- | --- |
| 1e13 | 5.61e7 | 3.24e7 | x0.58 | 3.56e7 | 0.64 |
| 1e14 | 1.30e9 | 1.90e9 | x1.47 | 2.09e9 | 1.61 |
| 1e15 | 1.68e10 | 2.16e10 | x1.29 | 2.38e10 | 1.42 |

CCE moves from 0.998181 / 0.992217 / 0.940703 to 0.998249 / 0.993286 / 0.947765, so 1-CCE
changes by a factor of 0.96 / 0.86 / 0.88. The result also depends on mesh refinement.

- **Severity:** major, not critical. It is on the default path but at most 7e-3 in CCE, and which
  average is physically correct (volume mean vs a tau(x) field) is a modelling choice.
  The node mean is not a defensible choice.
- **Fix:** take a control-volume-weighted mean over junction..epi_end, or better, pass tau(x) as
  a node field. Correct the docstring.
- **Confidence:** high that it depends on the discretisation; medium on which average is right.
  Found by B and D.

### M4. `trap_compensation_fluence` bypasses `check_trap_convention` (silent wrong number) (issue #15)

- **Where:** `etna/api/simulation.py:1828-1846`.
- **What is wrong:** an explicit `RadiationTrapModel` built in one convention is evaluated
  with kappa of another.
- **Evidence:** `lab/audit_v8/m4_trap_comp_convention.py`. A sic-damage preset under
  `niel_convention="si-neq"` gives 1.441e12 p/cm². The correct si-neq preset gives 1.696e12
  (-15 %). `compensation_fluence`, `default_fluence_grid` and `run_radiation_damage` all
  raise on the same input.
- **Same hole:** `trap_calibration.removal_observable`, flat-fluence branch
  (`trap_calibration.py:383-386`). A si-neq model under a sic-damage observation returns
  12.819 instead of raising; the matching model gives 12.742
  (`lab/audit_v8/n_flat_fluence_convention.py`).
- **Consequence:** USAGE item 3 ("every fluence sweep ... refuses") overstates the guard.
- **Fix:** call `check_trap_convention(model, resolve_niel_convention(niel_convention))` in
  both places.
- **Confidence:** high.

### M5. Tilted-alpha plasma slices ignore `DeviceConfig.mobility` (issue #16)

- **Where:**
  - `etna/core/track/cross_section.py:54` sets the module-level `PARAMS = SiC4H_Parameters()`,
    which `node_mobility` (:134) uses;
  - `CrossSectionDevice` has no mobility field;
  - `track/tilted.py::tilted_alpha_cce` never forwards `config.mobility`;
  - the column part honours it (`track/grid.py:164`).
- **Effect:** `run_cce_tilted_alpha(DeviceConfig(mobility="legacy_perp"))` mixes a legacy_perp
  column with c-axis slices. USAGE section 18 promises `legacy_perp` restores old numbers
  "bit for bit".
- **Evidence:**
  - `lab/audit_v8/m5_slice_mobility.py`: no mobility field, slice PARAMS is c_axis;
    mu_n(1e14) is 1162.4 vs 941.1.
  - `m5_slice_loss_mobility.py`: one slice (2e4 V/cm, 150 keV/um) has a recombined fraction
    of 2.193e-4 (c_axis) vs 2.363e-4 (legacy_perp), +7.8 %.
- **Severity:** major, as a silent drop of a physics option; the absolute CCE impact is small.
- **Fix:** add a `mobility` field to `CrossSectionDevice`, forward it, use `material_params(mobility)`,
  and include it in the PlasmaTable cache key.
- **Confidence:** high on the path; medium on the magnitude.

## Minor

- **N1. Fixed electron mobility of 950 in the resistivity and flat-C-V estimators.**
  `trap_calibration.py:356` (`semi_insulating_fluence(mu_n=950.0)`) and
  `radiation_traps.py:476` (resistivity proxy) hard-code the legacy_perp low-field
  electron mobility. The c-axis default is about 1160-1190 at low doping, so the proxy
  is about 20 % off and the flat-C-V estimate shifts slightly. Fix: take mu from
  `material_params(config.mobility)`.
- **N2. `calibrate_removal` turns NaN residuals into 1e3 with no warning.**
  Location: `trap_calibration.py:474-476`. With one observable forced to NaN
  (`lab/audit_v8/n_calibrate_removal_nan.py`): model_values [4.75, 504.0, 4.75],
  chi2 1.0e6, dof 2, warnings []. The fabricated value is reported as a model value.
  Fix: warn and drop the row, or raise.
- **N3. Out-of-range energy reported as a solver failure.**
  Location: `charge_collection.py:1208-1225`. In the core `cce_vs_fluence`, a broad
  `except Exception` turns an input error (energy outside the NIEL table) into a per-point
  NaN labelled "solver failure: ValueError ..." (`lab/audit_v8/n_energy_error_swallowed.py`,
  2000 MeV gives [0.99985, nan]). USAGE says this raises; the facade still does. Fix:
  validate kappa before the loop.
- **N4. `"auto"` evaluated at 0 V in the helpers and on the page.**
  `damage_model_of` / `compensation_fluence` default to V_bias = 0, and
  `app/workflows/radiation_damage.py:206-208` never passes the bias. For legacy geometry
  at -60 V the page announces the charged-trap model and a range up to 1e16 neq, but the
  run falls back to legacy. The `damage_model_note` reaches only the CSV header, and
  `run_all._compensation_note` ignores it (`lab/audit_v8/n_auto_model_vs_bias.py`:
  phi_c 4.70e12 at 0 / -20 V vs 1.14e13 at -60 V). Not silent at API level.
- **N5. `run_radiation_damage` docstring states the old default kappa.**
  `simulation.py:1624-1625` gives "kappa(5.6 MeV) = 6.67, Si-normalised"; the default
  (sic-damage) gives 5.89.
- **N6. Stale comments and docstrings.**
  - `radiation_damage.py:842`: "None = the default Si-neq proton table".
  - `radiation_damage.py:66`: kappa(252.7) "0.855"; it computes to 0.8514.
  - `tests/test_neutron_niel.py:130`: calls si-neq "the default".
- **N7. Stale module docstrings in the trap code.**
  - `trap_kinetics.py:44` lists the two-level Z1/2 under "does NOT do"; it is implemented
    and is in the default preset.
  - `radiation_traps.py:64-65`: "Lifetimes: unchanged (compute_K_tau)", while the default
    preset uses occupancy lifetimes.
- **N8. Z1/Z2 capture cross-sections taken from Capan's emission prefactor.**
  `trap_presets.py:103-118`. Capan et al. Table I (arXiv:1812.06462, opened) gives, for
  Z2(=/0), sigma_inf = 2e-16 cm² with a 0.03 eV barrier. The code's 1.02e-14 is 50x that.
  Emission matches Capan to 0.2-0.35 % over 250-350 K. Capture stays far faster than lab
  times either way, so there is no observable effect at present. Document it.
- **N9. A refit test that cannot catch a wrong d.**
  `test_refit_is_not_a_kappa_rescaling` (`tests/test_niel_convention_calibration.py:41-50`)
  only compares the shipped constants with each other. A wrong d would pass. Also, the
  `"auto"` fallback test checks `damage_model_of`, not the `run_radiation_damage`
  metadata end to end.
- **N10. Wrong source pointers.**
  `physics_models.py:396` and `:649` point to the wrong Ayalew thesis nodes. The BGN
  table 3.4 is on node63 and the anisotropy eqs. 3.80-3.81 on node65. The values are
  correct.
- **N11. Status labels disagree.**
  - USAGE section 18 labels `legacy_perp` "primary (Ayalew Table 3.5)"; the code says
    "verified-secondary".
  - The `c_axis` set is labelled "verified-primary" as a whole, although its T exponent and
    hole parameters are secondary.
- **N12 (open question; demoted from major after review; issue #14). Anchor of the
  legacy `eta_removal = 5 /neq` under the new default convention.**
  `radiation_damage.py:280`. It is used by the `"auto"` fallback, by
  `cce_vs_bias_at_fluence`, by C-V and dark current at fluence, and by
  `legacy_compensation_fluence`.
  - **Why not a bug:** the code labels the whole `RadiationDamageParams` set
    (`radiation_damage.py:283-284`) as `source="Burin ..."` and
    `reference_particle="1 MeV neutron equivalent"`. Under that anchor, carrying 5 /neq
    unchanged into sic-damage is exactly USAGE item 1. The per-proton drop
    (`lab/audit_v8/m3_legacy_eta_per_proton.py`) is the intended r_p/r_n correction:
    4.26 to 3.20 /p at 252.7 MeV, and 9.71 to 7.76 /p at 62 MeV.
  - **What is open:** `.planning/CARRIER_REMOVAL_REVIEW.md` section 4 calls the value
    mis-sourced and says Radmanovac's per-proton data (4.5-5.6 /p at 252.7 MeV) "happen to
    back it" in si-neq. Its section 5 recommends re-citing it to Radmanovac, which has not
    been done.
  - **If the re-citation is adopted:** the value becomes proton-calibrated and must be
    re-fitted per convention, about 7.0-8.8 /neq in sic-damage. Otherwise, drop the
    Radmanovac comparison in USAGE line 752.
  - Either way, the decision should be explicit.

Unverified (not counted): Ishikawa thesis eq. 2.16, the c-axis electron mobility. The
repository page could not be read; only the internal arithmetic was checked (1174.4 at
3.38e13), and the label should say so until someone opens it.

Out of scope, noted, not verified: `run_transit_time` (`etna/api/models.py:459-475`) does
not check `ramp.converged`, while `run_avalanche` does.

## Refuted candidates (checked, no defect)

- **Trap occupancy:**
  - SRH steady occupancy, single and two-level, equals Fermi-Dirac at E_F 0.2-1.8 eV
    (2.2e-16).
  - Charge per state is correct: acceptor -(P1+2P2), donor 2-P1-2P2.
  - Backward Euler and the 2x2 system are correct; `_pair_occ` derivative vs FD 7.8e-5.
  - Steady Poisson trap Jacobian vs FD 1.5e-3, which is truncation error.
  - dt = 0 and dt = inf limits are correct.
  - `capacitance_ac` series formula, `neff_from_cv` and `generation_current` units are correct.
- **Bulk trap numbers:**
  - Neutral-bulk removal: 7.400 /neq (Burin) and 12.090 /neq (default preset), constant to
    0.95 Phi_c.
  - USAGE bulk numbers reproduced: E_F 1.077 / 1.200 eV, N_eff = N_D + 1.432 Phi,
    EH6/7 10.5 % filled, emission times 0.040 s / 1774 s.
- **Sources match the code:**
  - Burin Table I (EH6/7 sigma_n deviation already disclosed, C-4).
  - Capan Table I energies and prefactors, [Z2]:[Z1] = 4.4.
  - Radmanovac Table 2, 4.48 +- 0.26 / 4.60 +- 0.33 / 5.60 +- 0.81, with the depth
    windows and fluences.
  - Hatakeyama APL 85 1380 Table I, and the SISPAD 2005 anisotropic formula (it reduces
    exactly to the axis values; the sqrt clamp is never active).
  - Review 2410.06798 Table XIII / eq. 42 (Niwa, Zhao, Loh, Stefanakis).
  - Ayalew Tables 3.5 / 3.6.
- **Device-model numerics:**
  - The McIntyre integrals and the 1 - I_n = e^{-phi}(1 - I_p) identity hold.
  - DD impact-ionisation sign and Jacobian are correct; devsim parity 3e-12;
    DD gain vs McIntyre 1e-4.
- **NIEL:** no double application of kappa or r_n; no unit errors (ASTM 95 MeV mb gives
  2.037e-3 MeV cm²/g); out-of-range interpolation raises; the neutron CSV column is
  correct; `arrhenius_sigma` matches Capan within 0.35 % over 250-350 K.
- **Damage paths:**
  - CCE clip to [0, 1] hides nothing: pre-clip values are <= 1 everywhere tested.
  - `trap_diagnostics` subtracts `donor_removal`.
  - `"auto"` fallbacks are logged and noted in metadata, not silent at API level.
  - No per-proton vs per-neq mixing on the trap path.
  - The kinetic dark state raises or NaNs on non-convergence.
  - Zero fluence is bit-identical to legacy, asserted with `==`.
- **Docs:** every compensation-fluence and kappa number in USAGE reproduces.
- **Tests:** the pinned CCE / fluence tests use tight tolerances (atol 2e-6 to 2e-5) and can fail.
- **Calibration:** the Gauss-Newton covariance is standard.

## Method

Four parallel sub-audits:

- A: mobility, avalanche and DD models;
- B: trap kernels and kinetics;
- C: NIEL, presets and calibration;
- D: damage dispatch and CCE paths.

Each one tried to refute its own candidates. The lead auditor then re-ran every major
candidate (all five `lab/audit_v8/m*` scripts) and the N2 and M4-flat-fluence scripts on
7746e8d; the quoted outputs come from those re-runs. The other minor numbers are the
sub-auditors' and were not re-run. After review, M3 was demoted to N12, because the code
anchors `eta_removal` per neq. Sources opened: Burin arXiv:2407.16710, Capan arXiv:1812.06462, Radmanovac
arXiv:2510.11304, Hatakeyama APL 85 (2004) 1380 and SISPAD 2005, Ayalew thesis (TU Wien),
and the review arXiv:2410.06798. The Ishikawa thesis was not reachable.
