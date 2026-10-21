# Audit v8: MC transport, I/O, export, app (2026-10-21)

Independent adversarial review of code written in the last day and never read by a
second reader: `etna/core/{mc_transport,straggling,mc_response,mc_exchange,electron_range,microdosimetry,stopping_power}.py`,
`etna/measurements/*`, `etna/export/*`, `etna/viz/*`, `etna/_paths.py` and packaging,
`app/**`. Base: `main` at 7746e8d. No code under `etna/` or `app/` was changed.

Method: every candidate was reproduced numerically or traced in code and then attacked
(refutation attempt). Only demonstrable defects are listed as findings. Geant4 /
Allpix Squared statements cite only source files that were opened (Geant4 v11.3.0
`G4IonFluctuations.cc`; Allpix Squared master `ElectricFieldReaderModule.cpp`,
`DetectorField.tpp`, `field_parser.h`). Repro snippets are inline because the scratch
scripts are not kept.

Summary: **0 critical, 4 major, 13 minor.** GitHub issues (label bug): #17 (finding 1), #18 (2), #19 (3), #20 (4).

| #    | Sev   | Where                                      | One line                                                                                                                                       |
| ---- | ----- | ------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| 1    | major | `etna/core/mc_transport.py:294-300`        | per-event RNG streams overlap: event k+1 = event k shifted by one draw; SEs 1.6-2.3x too small                                                 |
| 2    | major | `etna/core/mc_response.py:66-75`           | `default_l_bar_um` for an oblique beam is t/cos(theta), not the mean chord V/A_proj; y up to 2x low                                            |
| 3    | major | `etna/core/mc_exchange.py:369-399,440-453` | Allpix field maps span the 11.3 um stack but the snippet sets no `field_depth`/`doping_depth`; Allpix stretches them over the sensor thickness |
| 4    | major | `app/workflows/measurements.py:225-230`    | "same geometry as the sidebar device" snippet rebuilds the inverted table on `DeviceConfig()` defaults                                         |
| 5-17 | minor | see below                                  |                                                                                                                                                |

---

## Major

### 1. Overlapping per-event random streams (`etna/core/mc_transport.py:294-300`)

`_seed_state(seed, event)` sets the SplitMix64 counter to `seed*M1 + event*GOLD + M2`,
and every `_next_u64` adds `GOLD`. So event k's draw j uses counter
`S0 + (k + 4 + j) GOLD`, which is the same counter as event k+1's draw j-1. Each event's
stream is the previous event's stream shifted by one draw. The module docstring and
USAGE ("A SplitMix64 stream per (seed, event)") describe independent streams. They are
not independent.

Evidence (seed 0):

```python
s0, s1 = mt._seed_state(0, 0), mt._seed_state(0, 1)
[mt._uniform(s0) for _ in range(6)][1:]  # 0.403089 0.233619 0.611955 0.592489 0.569001
[mt._uniform(s1) for _ in range(5)]      # 0.403089 0.233619 0.611955 0.592489 0.569001
ev = mt.sample_track_events("proton", 100.0, (100, 100, 10), n_events=20000)
np.max(np.abs(ev.entry_um[1:, 0] - ev.entry_um[:-1, 1]))   # 0.0 exactly
```

- Normal beam: event i+1's entry x equals event i's entry y, bit for bit.
- Lag-1 autocorrelation of `edep_keV`: 0.72 for 100 MeV p in 10 um, and 0.22 for the
  USAGE example (62 MeV p, 100x100x10 um, defaults).
- Error bars:
  - 100 MeV p, 10 um, 20 seeds x 20000 events: the seed-to-seed sd of the mean deposit
    is 0.301 keV. The iid SE is 0.133 keV, so the SE is 2.3x too small.
  - 62 MeV USAGE setup, 30 seeds x 50000 events: seed-to-seed sd / reported
    `lineal_moments` SE = 1.89 for y_F and 1.58 for y_D.

Impact:

- Each event's marginal distribution is valid, so means and spectra are unbiased in
  expectation.
- Every iid standard error is too small, and the effective sample size is several times
  lower than n.
- `adaptive_lineal_energy` (`mc_response.py:119-123`) stops when the understated SE is
  below `rel_tol`, so a run reported as "converged at 1 %" carries about 1.6-1.9 % real
  error.
- Hit lists exported to Allpix share random numbers between consecutive events.
- The benchmark error bars in `lab/mc_transport` are understated.

Fix: derive each event's start state by hashing, e.g.
`st = mix64(mix64(seed) ^ mix64(event))` (SplitMix64 finaliser as `mix64`). Another
option is to give each event a disjoint block (counter += event * 2^40 draws).
Acceptance test: lag-1 autocorrelation is about 0 within noise, and the seed-to-seed sd
of y_F is about the reported SE (ratio ~1.0 +- 0.15 over 30 seeds). Note: tests that pin
exact numbers for a given seed (e.g. the USAGE `2.53 2.85` example) will move.

Confidence: high (exact equality shown).

### 2. Oblique-beam "mean chord" is t/cos(theta) (`etna/core/mc_response.py:66-75`)

The docstring says "Mean chord: t / cos(theta) for a parallel beam". The sampler spreads a
tilted beam uniformly over the projected area, entering through the top face and one side
face (`_sample_entry`, mode 0). Its mean chord is therefore V / A_proj, with
A_proj = wx wy cos(theta) + wy t sin(theta), and not the infinite-slab t / cos(theta).
`adaptive_lineal_energy` uses this value as l_bar by default, and
`tests/test_mc_response.py:36` pins the wrong value (20.0 at 60 deg).

| SV (um), angle        | sampled mean chord | V/A_proj      | `default_l_bar_um` | y_F(default) / y_F(true chord) |
| --------------------- | ------------------ | ------------- | ------------------ | ------------------------------ |
| 100x100x10, 60 deg    | 17.061             | 17.047        | 20.000             | 0.853                          |
| 10x10x10, 45 deg      | 7.084              | 7.071         | 14.142             | 0.501                          |
| 100x100x10, isotropic | 16.694             | 16.667 (4V/S) | 16.667             | ok                             |
| 100x100x10, normal    | 10.000             | 10            | 10.000             | ok                             |

(200000 events, `fluctuations="none", delta_rays=False`; y_F(true chord) equals the
LET 1.948 keV/um.)

Fix: return `wx*wy*t / (wx*wy*cos th + wy*t*sin th)` for the x-z tilt (the sampler's
convention), and correct the test. Alternatively, document t/cos(theta) as a slab
convention, not the mean chord.

Confidence: high.

### 3. Allpix field maps silently rescaled to the Allpix sensor thickness (`etna/core/mc_exchange.py:369-399, 440-453`)

`export_allpix_fields` writes the E-field and doping over the whole modelled stack. For
the default device that is 11.3 um: 0.3 um p+, 10 um epi and 1.0 um n+ (the label says
the n+ is not a fabrication thickness).

What Allpix does (`ElectricFieldReaderModule.cpp`, master, opened):

- `depletion_depth` (alias `field_depth`) defaults to the sensor thickness.
- `thickness_domain = (sensor_max_z - depletion_depth, sensor_max_z)`.
- `DetectorField` bins z over that domain. The file's own thickness is never compared
  with the sensor's.
- The DopingProfileReader uses `doping_depth` the same way.

The returned config snippet sets neither key, and USAGE does not mention sensor
thickness. Consequences:

- In a 10 um Allpix sensor (the sampler's default SV, also used by
  `write_deposition_csv`), the map is compressed by 11.5 %. The junction moves, and the
  potential integral is about 47 V instead of 53.2 V (default device, -50 V).
- In a 350 um sensor the map is stretched about 31x.
- The hit file (SV-centred, thickness `sv_um[2]`) and the field file (top = p+ surface)
  do not share a z frame.

Checked and correct: INIT layout, V/cm units, sign (Ez > 0 in the depleted epi,
integral Ez dz = potential drop to 1e-12), 2D mirror sign.

Fix: emit `field_depth = <th>um` / `doping_depth = <th>um` in the snippet. State that the
Allpix sensor must be at least that thick. Offer `depth_range_um` defaulting to the epi,
and warn when the exported thickness differs from the SV used for the hits.

Confidence: high on the mechanism (Allpix source), medium on how often users hit it.

### 4. Measurements page: C-V inversion snippet builds on the wrong device (`app/workflows/measurements.py:225-230`)

`invert_cv` fits the table on the sidebar device (`etna/measurements/cv.py:641`,
`table_config`). The code snippet shown to the user is
`cfg = etna.DeviceConfig(donor_table=table)` with the comment "# same geometry as the
sidebar device". It uses the `DeviceConfig()` defaults instead. Reproduced:

- Sidebar epi 8 um: the snippet raises `ValueError: donor_table must cover the epi: last
depth 8 um < epi thickness 10 um`.
- Sidebar legacy preset: the snippet silently changes
  - substrate 1.0 -> 0.3 um
  - n+ None -> 1.0 um
  - mobility legacy_perp -> c_axis
  - N_t 2.8e12 -> 1.36e12
  - E_t 0.65 -> 0.67

  The user's re-simulation of the inverted profile is then a different device.

Fix: emit `dataclasses.replace(<sidebar cfg as code>, donor_table=table, ...)`, or ship
the fit config in the JSON, and correct the comment.

Confidence: high.

---

## Minor

5. **Zero-deposit events bias the binned y_F** (`etna/core/microdosimetry.py:162-182`).
   - The f(y) normalisation divides by all events, including epsilon = 0 and
     out-of-range events, so y_F scales by the non-zero fraction while y_D is unchanged.
     Only a log warning is emitted.
   - Repro: gamma(4,5)/10 sample gives y_F = 2.013. Appending an equal number of zeros
     gives y_F = 1.007, and y_D = 2.512 in both cases.
   - `run_microdosimetry` rejects only all-zero files.
   - `mc_response.lineal_moments` correctly drops zeros, so the two modules disagree.
   - Fix: normalise over events with y in range (or epsilon > 0), and report the dropped
     count in metadata.
6. **NaN `ekin_MeV` with a surround** (`etna/core/mc_exchange.py:231-251`).
   - Events that miss the SV but score delta energy (12 of 5000 for 20 MeV p, 20x20x10
     SV, 50 um surround) get `ekin_MeV = NaN` and position 0.
   - `process_mc_ensemble` keeps them, and `tissue_equivalence_correction(...,
particle_energy_MeV=...)` then raises "must be finite". This is a loud failure in
     the documented chain.
   - Fix: use the primary's energy at the delta's production point, or drop or flag
     these events.
7. **`iv_t` weights dropped on one empty `sigma_I`** (`etna/measurements/iv_t.py:310-313,
446-451`).
   - One missing sigma makes the activation-energy fit unweighted, with the note "no
     sigma in the file".
   - Repro: E_a = 0.981 +- 0.031 eV unweighted against 1.0000 +- 0.0016 eV with all
     sigmas, on data with a 50 % outlier carrying a 50 % sigma.
   - Fix: weight the rows that have a sigma, or refuse the file.
8. **C-V inversion discards all sigmas on one empty `sigma_C`** (`etna/measurements/cv.py:176,
867-868, 992-996`).
   - The missing value becomes sigma_W = 0, so every point gets the default 0.02 um
     (2-4x the given sigmas), with a false "no sigma in the file" note.
   - The weights and the reported uncertainty change. On noisy data the weighted
     estimate also moves (on the noise-free repro the fitted values did not).
   - Fix as in 8.
9. **Non-UTF-8 measurement files crash the page** (`etna/measurements/_csv.py:70,75,81`;
    `app/workflows/measurements.py:99,118,249,351` catch only `MeasurementFileError`).
    - A Latin-1 "°C"/"µm" comment raises `UnicodeDecodeError` and shows a traceback.
    - Fix: re-raise as `MeasurementFileError` or fall back to cp1252.
10. **Run All "About" text describes the old damage model** (`app/workflows/run_all.py:44-47`).
    - It says 8 fluences from 0.01-3x the compensation fluence, beyond-compensation
      points refused, and charged traps not used.
    - The default device now runs calibrated kinetic traps on 10 points up to about
      3574x, all solved. The report table itself is correct.
11. **Dark-current CSV records the sidebar N_t when the Advanced override is used**
    (`app/workflows/dark_current.py:183-200`, header from `app/components/results.py:658-660`).
    - Leakage is linear in N_t, so the file mislabels its own calibration.
12. **Radiation Damage page: switching neutron -> proton is not flagged as stale**
    (`app/workflows/radiation_damage.py:258-265`, `app/components/page_utils.py:107-110`).
    - Particle and convention are recorded only when non-default, and the stale check
      compares only current keys.
13. **1D "Net doping vs depth" drops the p+ anode** (`app/workflows/field_map.py:126-136`).
    - Signed values on a log axis: 20 negative nodes (down to -1.32e18) are not drawn.
14. **Other stale-state gaps:**
    - Figures page (`figures.py:620-627`): only the device is checked; V_stop, single
      bias and the LET table are not.
    - Microdosimetry page (`microdosimetry.py:110-112`): the upload key is name + size,
      so a different file with the same name and size is not flagged.
15. **Mislabelled axes and captions:**
    - `geometry_viewer.py:38-47` titles 1D bar charts "2D ...".
    - `radiation_damage.py:69` says "proton fluence" in neutron mode.
    - The neutron help says "to 18.9 MeV" while the input stops at 18.8.
    - `results.py:540` names the column `fluence_p_per_cm2` for neutrons.
16. **Spec sheet W(-10 V) interpolated on a 1.28 V grid** (`etna/export/spec_sheet.py:299-309`).
    - For `nplus_calibrated_config` it shows 9.4723 um against 9.5000 um with -10 V on
      the grid. That 0.3 % apparent miss is on a value documented as "solved exactly".
    - Fix: include 0/-10/-30 V in the spec's C-V grid.
17. **VTU comment sanitiser** (`etna/export/mesh_export.py:219`).
    - `replace("--", "- -")` turns `"---"` into `"- --"`, so a free-text note containing
      `--->` would end the XML comment early. No current string triggers it.
    - Fix: loop until no `--` remains.

Unverified (not counted): the `mc_exchange.py:36` statement "DopingProfileReader:
positive = n-type". The opened Allpix module source does not state the sign convention.

---

## Checked and refuted / not retained

- **Fluctuation sampler mean.** Single-seed runs looked biased by up to -3 %. That was
  an artefact of finding 1 plus Landau tails. The isolated `_sample_fluct` with 2e6 draws
  reproduces its input mean within 0.02 % (urban, universal, bohr; protons 10/100 MeV,
  1-10 um, tcut = tmax and 5 keV). A 40-seed transport check gave -0.6 % at -2.6
  seed-SE (skewed per-seed means). Inconclusive, not retained.
- **Geant4 port** against v11.3 source and Urban's formulas, line by line:
  - `G4UniversalFluctuation` Gaussian/Glandz branches;
  - Urban two-level a1/a2/a3, alfa;
  - the `ComputeFluctModel` Zeff/f1/f2/E2;
  - `G4IonFluctuations` threshold `parameter*charge*mass`, `sqrt` of the dispersion,
    the 0.2/0.2 slowing-down correction and the gauss/gamma/flat branches;
  - Ziegler He effective charge.

  Only difference: the port's 10 eV minLoss shortcut on the ion branch (Geant4 uses
  0.001 eV), which is negligible.

- **Delta rays:**
  - restricted-loss subtraction `pref [ln(Tmax/Tcut) - b2 (1 - Tcut/Tmax)]` and the
    cross section `pref [(1/Tc - 1/Tm) - b2/Tm ln(Tm/Tc)]` equal the Bethe difference;
  - 1/T^2 inversion and (1 - b2 T/Tmax) acceptance are correct;
  - emission angle cos^2 = T(Tmax + 2m)/(Tmax(T + 2m)) reduces to the exact
    T/((T + 2m) b2) for M >> m;
  - straight-line CSDA bookkeeping of entry/exit energies is consistent with
    `delta_range_factor`.
- **Energy conservation:**
  - 5.486 MeV alpha stopping in 60 um: mean 5473.65 keV against the CSDA electronic
    integral 5473.70 keV (the rest is nuclear loss, not deposited by design);
  - never above E0;
  - sum of hit energies equals edep to 3e-14 keV.
- **Geometry and analytic references:**
  - isotropic chord mean 16.69 um = 4V/S 16.67 um;
  - normal chord = t;
  - Vavilov density (`straggling.py`): normalisation, mean and variance reproduce
    xi Tmax (1 - b2/2) to < 1e-7;
  - Landau mode -0.2228 and FWHM 4.018 xi (scipy 1.18 convention mapped correctly).
- **Ranges and kappa:**
  - ESTAR power-law start R = T/(S(1 + p)) is consistent;
  - Bohr range straggling integral is correct;
  - kappa linear/mass conventions match the docstrings (memory: 0.43 -> 0.375).
- **Hit track ids:** delta-ray hits are labelled 2, 3, ... as documented (`mc_transport.py:955`, `nd` is incremented first).
- **Allpix INIT format, E sign, 2D mirror:** correct (see 3).
- **C-V inversion formulas:**
  - N = 2/(q eps d(1/C'^2)/dV) recovers 1e14 to 6.5e-6;
  - W = eps/C';
  - pF/mm^2 and reverse-sign handling are correct.
- **Measurement loaders:** semicolon / decimal-comma files are refused loudly.
- **GDS:** units 1e-6/1e-9, polygon areas and pad/ring/die extents match the TOML.
  Placeholders are labelled PLACEHOLDER, and the top cell is
  `...REFERENCE-NOT-FAB`. Spec-sheet numbers (V_fd, W, I_dark) are reproduced.
- **Overclaims:** no "validated" overclaim found in app/ or export labels.
- **Packaging:**
  - The wheel contains all 38 data files.
  - In a clean venv with `python -I`, the data root resolves to `etna/_data` and ESTAR,
    NIEL, mc_transport and the layout all load.
  - `scripts/wheel_smoke.py` does not exercise mc_transport, electron_range or
    neutron_niel (coverage gap only).
- **App caching:** no `st.cache_*` in use. Session caches keep the source object next to
  id keys, and caches are capped except on the figures page (memory only).

---

## Resolution (2026-10-21, branch worktree-agent-ad02bbb6dea5fe590)

Each fix has a regression test that fails on edf4483 and passes after the fix.

| # | commit | test |
| --- | --- | --- |
| 1 (#17) | c397508 | `tests/test_mc_rng_streams.py`: lag-1 0.72 -> 0.006; seed scatter / SE 1.99 / 1.53 -> 0.99 / 0.94 |
| 2 (#18) | 6cf3d35 | `tests/test_mc_response.py::test_default_l_bar_is_sampled_mean_chord` |
| 3 (#19) | 3a4dd96 | `tests/test_mc_exchange_allpix_depth.py`: the Allpix placement rule keeps 53.2 V |
| 4 (#20) | 72d0ace | `tests/test_measurements_rebuild_code.py` |
| 5-17 | 39dbe7f, d38387f | `tests/test_audit_v8_minors.py` |

Notes on the minor findings:

- **12.** On main, a neutron -> proton switch was already flagged, but only through the
  fluence and energy keys, which change with the particle. The particle key itself is
  now always compared.
- **13, 14, 15 and 10.** The page-level changes (1D doping traces, figures stale key,
  microdosimetry content hash, captions, Run All text) have no new dedicated AppTest.
  The existing page tests pass.
- **Unverified item (DopingProfileReader sign, positive = n-type).** Still unverified.
  The `mc_exchange` docstring now says so.
