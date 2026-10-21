# Data request checklist (INFN-LNS and foundry)

What ETNA needs from the measurement owners (INFN-LNS) and from the foundry to turn the
current calibrated-by-digitised-figure model into one calibrated on real data. Each
item names the ETNA input it feeds, the file template, and the checks the loader
applies. Formats: `data/measurements/README.md` (C-V, SIMS/SRP, I(T)) and
`data/foundry/layout_params.toml` (mask dimensions). A worked SYNTHETIC example of every
measurement file, with a round-trip script, is in `examples/measurements/`
(`uv run python examples/measurements/roundtrip.py`).

Priority: **A** = blocks a calibration in use today, **B** = removes a documented
assumption, **C** = useful cross-check.

Before sending, a file can be checked locally:
`uv run python -c "from etna import measurements as M; M.load_cv_csv('file.csv')"`
(or `load_sims_csv` / `load_it_csv`). The loader refuses malformed files and obvious unit
mix-ups (F in a pF column, nm in `depth_um`, m^-3 in `N_cm3`, pA in `I_A`, impossible
areas) with one message per problem naming the line and column, and warns
(`MeasurementPlausibilityWarning`) on unusual values (depletion width outside
0.05-500 um, |V| > 1.5 kV, T outside 77-700 K, C rising with reverse bias, ...).

## 1. Foundry (process and layout)

| #   | item                                                                                                                                                                                          | pri | ETNA input                                                                                                        | notes / what to state                                                                                                                                                                                                        |
| --- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --- | ----------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| F1  | **Mask dimensions** of both devices: active p+ side, gap to the first guard ring, guard-ring count / width / spacing, bond-pad size, die side, dicing street, oxide and passivation thickness | A   | `data/foundry/layout_params.toml` (set `class = "DOCUMENTED"`, `source = "mask set <name> rev <n>"`, uncertainty) | every PLACEHOLDER there keeps the GDSII top cell `REFERENCE-NOT-FAB`; the active area also sets the absolute C and dark current (`DeviceConfig.area_cm2`)                                                                    |
| F2  | **SIMS** of N and Al vs depth on a witness / test structure of the same wafer                                                                                                                 | A   | `templates/sims_template.csv` (`technique: SIMS`)                                                                 | depth from the surface, crater-depth calibration method, **the N and Al detection limits** (columns or `# detection_limit_*`), relative uncertainty; depth to at least the epi end + 1 um so the n-/n+ transition is visible |
| F3  | **SRP** (spreading resistance) or ECV of the same stack                                                                                                                                       | B   | same template, `technique: SRP`                                                                                   | net carrier density; gives the epi doping where SIMS N is below its detection limit (see the note below)                                                                                                                     |
| F4  | Epi specification: thickness and doping (nominal and wafer-map / certificate values), buffer layer, substrate doping                                                                          | A   | `DeviceConfig.epi_thickness_um`, `N_D_*`, `buffer_*`                                                              | the current epi doping is inferred from a digitised C-V figure                                                                                                                                                               |
| F5  | p+ anode: Al implant dose / energy / anneal or p+ epi thickness and doping                                                                                                                    | B   | `substrate_thickness_um`, `N_A`                                                                                   | sets the junction depth used to place SIMS on the model axis                                                                                                                                                                 |
| F6  | Contact / metal stack and the passivation material                                                                                                                                            | C   | export cross-section                                                                                              |                                                                                                                                                                                                                              |

**SIMS floor**: if the N detection limit of the SIMS lab is above the epi doping
(~1e13-1e14 cm^-3 here), the epi knots are flagged _censored_ and
`ProfileTable.to_device_config` refuses them (an upper bound is not a doping value).
Then SIMS still gives the p+ profile, the junction depth and the n-/n+ transition, and
the epi doping must come from C-V (INFN-LNS item L1) or SRP (F3). Ask the lab for its
actual limit; ETNA does not assume one.

## 2. INFN-LNS (device measurements)

| #   | item                                                                                                                       | pri | ETNA input                                       | notes / what to state                                                                                                                                                                                                                                                                                                                       |
| --- | -------------------------------------------------------------------------------------------------------------------------- | --- | ------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| L1  | **C-V from 0 to -200 V** on pristine devices (both designs), 0.5-1 V steps up to -20 V and <= 5 V beyond, several dies     | A   | `templates/cv_template.csv`                      | area, test frequency, temperature, open/short correction, guard ring biased or floating, bias convention (reverse negative, or `# bias_convention: reverse_positive`), 1-sigma of C if known. The model solves to -200 V on the batched engine (example file). Feeds `compare_cv`, `invert_cv` (donor-table inversion), `depletion_profile` |
| L2  | **Dark current I(T)** at fixed reverse biases (e.g. -10, -30, -60 V) for T from 290 to 380 K or wider                      | A   | `templates/it_template.csv`                      | one row per (bias, T); guard ring connected or not; settling time; temperature sensor and its uncertainty. Feeds `fit_arrhenius` / `analyze_it` (activation energy vs the model's) and the N_t / E_t dark-current fit                                                                                                                       |
| L3  | I-V at room temperature to -200 V (same dies as L1)                                                                        | B   | I(T) template with one T, several biases         | the current N_t fit rests on one bound (< 18 pA at -60 V)                                                                                                                                                                                                                                                                                 |
| L4  | C-V and I-V **after irradiation** at known fluences (and after each annealing step)                                        | B   | same templates, with `# notes:` fluence / anneal | for the radiation-damage chain (carrier removal, trap model)                                                                                                                                                                                                                                                                                |
| L5  | Charge-collection efficiency vs bias with the alpha source (energy, source-to-device distance, collimation, vacuum or air) | B   | CCE comparison (`run_cce`)                       | the measured microdosimeter CCE shape is not reproduced by any DD variant yet                                                                                                                                                                                                                                                               |

## 3. INFN-LNS (irradiation metadata) -- needed to use the NIEL scaling

| #   | item                                                                                                                                                  | pri | ETNA input                                             | notes                                                                                                              |
| --- | ----------------------------------------------------------------------------------------------------------------------------------------------------- | --- | ------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------ |
| R1  | Particle and energy at the device (after degraders / windows)                                                                                         | A   | `run_radiation_damage(proton_energy_MeV=, particle=)`  | kappa(E) is steep below ~20 MeV                                                                                    |
| R2  | Fluence and how it was measured (Faraday cup, ionisation chamber, activation foils) with its uncertainty                                              | A   | `fluences`                                             |                                                                                                                    |
| R3  | **Which NIEL convention any quoted "neq" fluence uses** (Si NIEL / ASTM E722 hardness, or a SiC-specific one, and from which table)                   | A   | `niel_convention` (`si-neq`, `sic-damage`, `sic-1mev`) | the default Si-neq factor and the SiC-damage factor differ by r_n/r_p (x1.25 at 62 MeV), see `data/niel/README.md` |
| R4  | For any **reactor / neutron** irradiation: the neutron spectrum (or at least the fast / epithermal / thermal flux split and the hardness factor used) | B   | r_n of `data/niel/README.md`                           | the SiC/Si ratio r_n depends on the spectrum (1.02-1.20 for plausible shapes); a measured spectrum fixes it        |
| R5  | Temperature during and after irradiation, time to measurement, annealing steps (T, t)                                                                 | B   | `compute_annealed_params`                              |                                                                                                                    |

## 4. What happens to each file in ETNA

- C-V -> `load_cv_csv` -> `compare_cv` (residual W and C vs the model) -> `invert_cv`
  (donor table) -> `DeviceConfig(donor_table=...)`.
- SIMS/SRP -> `load_sims_csv` -> `resample_profile` -> `ProfileTable.to_device_config`.
- I(T) -> `load_it_csv` -> `analyze_it` (data E_a vs model E_a at the same T and bias).
- Every output (JSON, Markdown report, figures) carries the file SHA-256 and metadata.

Values that cannot be verified from these data stay flagged (PLACEHOLDER dimensions,
censored SIMS knots, "calibrated" vs "predicted" labels in `docs/USAGE.md`).
