# Measurement input files

ETNA reads owner-supplied measurements (INFN-LNS, foundry) as CSV files and uses them
through `etna.measurements` (Python) or the **Measurements** page of the app. This
directory holds **header-only templates** (`templates/`): copy one, fill the metadata
lines and add data rows. **No measurement data is committed here.** Put real files
wherever you like (they are read by path or uploaded) and record where they come from in
the metadata lines.

## Common format

```
# free comment
# key: value          <- metadata (known keys below; an empty value means "not given")
COLUMN_unit,COLUMN_unit,...    <- header: the UNIT IS PART OF THE COLUMN NAME
value,value,...
```

- Comma separated, one header line, `#` lines anywhere are comments / metadata.
- Every column name must be one of the accepted names below (unknown names are an error
  that lists the accepted ones); two names for the same quantity (e.g. `C_pF` and `C_F`)
  are an error.
- Required cells must be finite numbers; optional cells may be left empty.
- Errors name the **file line** (1-based, counting comment and header lines) and the
  **column**; all problems of a file are reported at once.
- Bias sign: ETNA convention, **reverse bias negative**. A file with reverse bias
  positive adds `# bias_convention: reverse_positive` (V is then negated on load).
  Mixed-sign files are allowed (forward points are kept but not compared with the
  model); an all-forward file is refused as a probable sign error.
- Free-text metadata keys (kept in the provenance of every output): `sample`, `device`,
  `operator`, `date`, `instrument`, `source`, `notes` (SIMS: `lab` instead of
  `operator`/`device`).

## C-V — `templates/cv_template.csv`

| column                                                         | unit               | required               | meaning                                                       |
| -------------------------------------------------------------- | ------------------ | ---------------------- | ------------------------------------------------------------- |
| `V_V`                                                          | V                  | yes                    | bias (reverse negative, or see `bias_convention`); no repeats |
| `C_F` / `C_nF` / `C_pF`                                        | F / nF / pF        | one capacitance column | absolute small-signal C (needs an area)                       |
| `C_per_area_F_cm2` / `C_per_area_nF_cm2` / `C_per_area_pF_mm2` | per area           | (alternative)          | C/A                                                           |
| `sigma_C_F` / `sigma_C_nF` / `sigma_C_pF`                      | same family as C   | no                     | 1-sigma of C                                                  |
| `sigma_C_per_area_F_cm2` / `_nF_cm2` / `_pF_mm2`               | same family as C/A | no                     | 1-sigma of C/A                                                |
| `sigma_V_V`                                                    | V                  | no                     | 1-sigma of the bias (enters W through dW/dV)                  |
| `area_cm2` / `area_mm2`                                        | cm^2 / mm^2        | no                     | device area (one value per file)                              |
| `frequency_Hz` / `frequency_kHz`                               | Hz / kHz           | no                     | test-signal frequency (one value per file)                    |
| `T_K`                                                          | K                  | no                     | temperature (one value per file)                              |

Metadata: `area_cm2`, `area_mm2`, `frequency_hz`, `temperature_k`, `bias_convention`.
Area precedence: column, then metadata, then the `area_cm2` argument of `load_cv_csv`
(the app passes the sidebar device area and says so); a disagreement is an error, an
absolute C without any area is an error. Internally C is per area (F/cm^2) and
`W = eps / (C/A)` with the model's eps (4H-SiC, eps_r 9.7).

## SIMS / SRP — `templates/sims_template.csv`

| column                                             | unit     | required      | meaning                                                        |
| -------------------------------------------------- | -------- | ------------- | -------------------------------------------------------------- |
| `depth_um` / `depth_nm`                            | um / nm  | yes           | depth **from the surface** (p+ top), strictly increasing, >= 0 |
| `N_cm3`                                            | cm^-3    | one of N / Al | nitrogen (donor) concentration; 0 = below detection            |
| `Al_cm3`                                           | cm^-3    | one of N / Al | aluminium (acceptor) concentration                             |
| `sigma_N_cm3` / `sigma_Al_cm3`                     | cm^-3    | no            | absolute 1-sigma                                               |
| `rel_sigma_N` / `rel_sigma_Al`                     | fraction | no            | relative 1-sigma (instead of the absolute one)                 |
| `detection_limit_N_cm3` / `detection_limit_Al_cm3` | cm^-3    | no            | per-row detection limit                                        |

Metadata: `technique` (`SIMS` = total atoms, default; `SRP` = net carrier density),
`quantity` (free text), `depth_reference` (must be `surface`), `junction_depth_um`
(p+/n- junction below the surface), `detection_limit_N_cm3`, `detection_limit_Al_cm3`,
`rel_sigma_N`, `rel_sigma_Al` (global values when there is no column).

Junction depth used to convert to the `donor_table` depth (below the junction): the
`junction_depth_um` argument, else the metadata, else the Al/N crossing, else the model
p+ thickness (`config.substrate_thickness_um`); the choice is recorded. The table stops
at the epi end; the n+ substrate stays `N_D_nplus`. Knots at or below the detection
limit are flagged and refused by `to_device_config` unless `allow_censored=True`.

## I(T) — `templates/it_template.csv`

| column                              | unit        | required | meaning                                              |
| ----------------------------------- | ----------- | -------- | ---------------------------------------------------- |
| `T_K` / `T_C`                       | K / degC    | yes      | temperature                                          |
| `I_A` / `I_nA` / `I_pA` / `I_fA`    | A ... fA    | yes      | dark current (sign kept, \|I\| fitted); non-zero     |
| `sigma_I_A` / `_nA` / `_pA` / `_fA` | same choice | no       | 1-sigma of I                                         |
| `V_V`                               | V           | no       | bias (several biases per file allowed: one fit each) |
| `sigma_T_K`                         | K           | no       | recorded only                                        |

Metadata: `bias_V` (when there is no `V_V` column), `bias_convention`, `area_cm2`
(informative: E_a does not depend on area). No repeated (bias, T) pair.

## Outputs

- donor-table JSON (`etna.donor_table/1`): `donor_table` [[depth_um, N_D], ...], per-knot
  relative uncertainty, method details, notes, limits, data hash and model provenance;
  `etna.measurements.load_donor_table_json(path)` returns the tuple for
  `DeviceConfig(donor_table=...)`.
- E_a JSON (`etna.measurements.activation_energy`), C-V comparison JSON, a Markdown
  report (`report_markdown`) and PNG/SVG/PDF figures (`etna.measurements.plots`).

Methods, limits and examples: `docs/USAGE.md`, section "Measurement inputs".
