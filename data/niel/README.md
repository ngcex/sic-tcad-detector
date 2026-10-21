# Proton NIEL data (SR-NIEL) for the radiation-damage hardness factors

Fetched 2026-10-20 with `fetch_sr_niel.py` (this directory) from the public SR-NIEL web
calculator, "Protons & Ions NIEL Calculator"
(https://www.sr-niel.org/index.php/sr-niel-web-calculators/niel-calculator-for-electrons-protons-and-ions/protons-ions-niel-calculator).
Code version printed on every result page: **SR-NIEL-7 ver. 11.1, 12 May 2026**.
Form settings: incident particle 1:H, model **Hadron + Coulomb**, hadronic contribution
**No Scaling**, energies 0.1-1000 MeV plus the extra energies listed in the script
(5.6, 62, 70, ... MeV), fluence 1 cm^-2. Robinson partition function (the free
calculator). Each `.html` is the raw result page as served; each `.csv` is its table
parsed by the script (columns: proton kinetic energy in MeV, NIEL in MeV cm^2 g^-1).
The result-page URLs in the CSV headers are temporary on the SR-NIEL server.

| file                           | target                        | displacement thresholds E_d                           |
| ------------------------------ | ----------------------------- | ----------------------------------------------------- |
| `sr_niel_proton_si_Ed21.*`     | Si (user defined, Z = 14)     | 21 eV (SR-NIEL default for Si)                        |
| `sr_niel_proton_sic_Ed21-21.*` | compound (IV) SiC, Si:C = 1:1 | Si 21, C 21 eV (calculator default)                   |
| `sr_niel_proton_sic_Ed35-21.*` | compound (IV) SiC             | Si 35, C 21 eV (sublattice-specific, sensitivity run) |

E_d sensitivity of the SiC proton NIEL (35/21 vs 21/21 eV): -5 % at 1 MeV, -4 % at
5.6 MeV, -2 % at 30 MeV, -1.6 % at 62 MeV, -1.4 % at 150 MeV.

## The 1 MeV neutron denominator

- **Silicon** (used by the default table): the ASTM E722 normalising value
  D_n(1 MeV) = 95 MeV mb, i.e. NIEL_n,Si(1 MeV) = 95e-27 cm^2 MeV x N_A / 28.0855 g =
  **2.0370e-3 MeV cm^2 g^-1** (Fretwurst, Lindstroem, Pintilie, Stahl,
  arXiv:physics/0211118: "normalizing value: Dn(1MeV) = 95 MeVmb", "100 MeVmb =
  2.144 keVcm2/g"). The point-wise ASTM E722-19 damage function around 1 MeV is
  resonant (1.1e-3 to 3.6e-3 MeV cm^2/g between 0.5 and 2 MeV, SR-NIEL ASTM neutron
  calculator), which is why the defined 95 MeV mb value is the normalisation, not a
  point value.
- **SiC: computed here (2026-10-21), see "Neutron NIEL (NJOY HEATR)" below.** SR-NIEL's
  NJOY-based neutron NIEL for SiC is behind a registered login, and no open table of
  NIEL_n,SiC(1 MeV) was found in the literature searched (SR-NIEL pages, De Napoli 2022,
  Gsponer arXiv:2310.02047, Radmanovac arXiv:2510.11304, Lazanu physics/0007011 (protons
  and pions only), IRDFF-II arXiv:1909.03336 (Si and GaAs 1-MeV damage, no SiC), IAEA
  INDC(NDS)-0648 (SiC only planned); the OSU thesis osu1173103438 states "The values for
  F_D,SiC(E) are not directly available in the literature" and only plots one; Heinisch
  et al., J. Nucl. Mater. 327 (2004) 175 tabulate SiC dpa cross sections (MD-based, not
  Lindhard NIEL; full text not opened)). So the SiC numbers below have no published
  point to compare with; their checks are the Si ones and the internal ones.

## Why the default hardness factor is Si-normalised

The defect introduction rates of `RadiationDamageParams` (Burin et al.,
arXiv:2407.16710, Table I, "g_int [cm^-1]") are per **1 MeV equivalent neutron fluence**
as labelled by the reactor-neutron data they were fitted to (Gsponer et al.,
arXiv:2310.02047, TRIGA Vienna, "1 MeV equivalent neutron fluences between 10^14 neq/cm2
and 10^16 neq/cm2"; Hazdra et al. NIMB 327 (2014) for g_int). The same group (Burin
co-author) converts protons to that fluence with **silicon** NIEL: Radmanovac et al.,
arXiv:2510.11304: "According to non-ionizing energy loss (NIEL) data in silicon, ...
252.7 MeV proton beam (k ~ 0.9)" and "equivalent to 9x10^13 neq./cm2 in silicon".
The default table therefore is

    kappa_Si(E) = NIEL_p,Si(E) / NIEL_n,Si(1 MeV)          (sr_niel_proton_si_Ed21.csv)

which reproduces that convention (SR-NIEL gives 0.855 at 250 MeV vs the quoted ~0.9).

What it assumes: the factor that converts a proton fluence into the damage of the
reactor-neutron calibration in **SiC** is

    kappa(E) = NIEL_p,SiC(E) / NIEL_n,SiC(reactor) x NIEL_n,Si(reactor) / NIEL_n,Si(1 MeV)

and kappa_Si equals it only if the SiC/Si NIEL ratio of protons, r_p(E), equals that of
the reactor neutrons, r_n. r_p is measured here: 0.97 (1 MeV), 0.94 (5.6), 0.88 (30),
0.85 (62), 0.82 (150), 0.80 (250 MeV). r_n was the blocked quantity until 2026-10-21;
it is now computed with NJOY (section "Neutron NIEL" below): r_n = 1.069 for a fission
spectrum (plausible range 1.02-1.20 for a reactor spectrum), > r_p at every energy, so
kappa_Si OVERSTATES the SiC damage per proton by r_n / r_p (x1.25 at 62 MeV with the
fission r_n). The default stays kappa_Si; the corrected factor is the opt-in
`sic-damage` convention.

## Known systematic of NIEL scaling at low energy

Measured silicon proton hardness factors are below the NIEL prediction at low energy:
3.07 +- 0.37 at 10.5 MeV and 2.19 +- 0.22 at 24.3 MeV (arXiv:2408.12410), vs SR-NIEL
kappa_Si = 4.62 at 10 MeV and 2.99 at 25 MeV here (35-50 % above). The 5.6 MeV API
default sits in this regime.

## Neutron NIEL (NJOY HEATR), 2026-10-21

Computed, not downloaded: `scripts/build_neutron_niel.py` runs **NJOY2016.78**
(conda-forge `njoy2016` 2016.78; RECONR 0.1 % then HEATR, partial kerma MT444 = damage
energy production, 0 K) on **ENDF/B-VIII.0** Si-28/29/30 and C-12/13 (NNDC
`ENDF-B-VIII.0_neutrons.zip`; the SHA-256 of every ENDF file is in
`etna/core/neutron_niel.py:ENDF_FILES` and in each output header). HEATR's damage energy
is the Lindhard partition in Robinson's fit with a sharp displacement threshold E_d,
recoil in its own element -- the same partition family as the SR-NIEL proton tables
("Robinson partition"). Output is byte-reproducible (run twice, identical SHA-256).

| file                                                            | content                                                                                                                                                                                                 |
| --------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `neutron_damage_<iso>_Ed<ed>.csv.gz`                            | pointwise MT444 (energy eV, damage eV barn) per isotope, NJOY grid (1.4k-9.2k points)                                                                                                                   |
| `neutron_niel_si_<tag>.csv`, `neutron_niel_sic_<tag>.csv`       | natural-abundance NIEL (MeV cm^2/g), lethargy average over 20 bins/decade, 10 eV-20 MeV, bin edges + centre; the header lists sha256[:16] of the pointwise inputs                                       |
| `si_neutron_damage_griffin_rd48.csv` (+ `convert_si_dfuncs.py`) | Si damage function D(E)/(95 MeV mb) of Griffin et al. SAND92-0094 as compiled by Vasilescu & Lindstroem (RD48, `neutrons.xls`, archived sesam.desy.de page), verbatim -- the reference for the Si check |

Tags: `Ed21` (Si 21 / C 21 eV, as the SR-NIEL proton tables) and `Ed35-21` (Si 35 / C
21 eV, the sublattice values of Devanathan & Weber, J. Nucl. Mater. 278 (2000) 258, also
used by Lazanu physics/0007011 "35 eV on the Si sublattice, and 20 eV on the C"). For
fast neutrons E_d changes nothing at the 5th digit (keV recoils; the threshold only
drops the few recoils below E_d) -- unlike protons (Coulomb recoils, -1.6 % at 62 MeV).
Natural abundances: Si 92.223 / 4.685 / 3.092 %, C 98.93 / 1.07 % (IUPAC representative
compositions); per-gram normalisation with 28.0855 and 12.011 g/mol; SiC = 1:1 molecule.

### Results (embedded in `etna/core/niel_tables.py`; a test recomputes them)

| quantity                                                                 | value                                                           |
| ------------------------------------------------------------------------ | --------------------------------------------------------------- |
| NIEL_n,Si(1.000 MeV), pointwise                                          | 2.0322e-3 MeV cm^2/g = **0.998 x ASTM 95 MeV mb**               |
| NIEL_n,SiC(1.000 MeV), pointwise                                         | **2.2989e-3 MeV cm^2/g** (76.5 MeV mb per atom); SiC/Si = 1.131 |
| Watt fission spectrum (a 0.988 MeV, b 2.249 /MeV, 1 keV-20 MeV) averages | Si 2.1806e-3 (kappa_Si,fission 1.0705), SiC 2.3303e-3           |
| **r_n = SiC/Si, fission spectrum**                                       | **1.0687** (`R_N_SIC_SI_FISSION`)                               |

### Checks

- **ASTM E722**: the pointwise Si value at 1 MeV is 0.998 of the 95 MeV mb normalisation
  (which is a defined value, not a point value: the damage function is resonant there;
  lethargy bins 0.9-1.0 / 1.0-1.1 MeV give 1.17 / 0.81 x 95 MeV mb).
- **Griffin (SAND92-0094) Si damage function**: Watt-fission average 1.0705 vs 1.0628
  (+0.7 %); per lethargy bin +0.1 to +2 % below 10 MeV, +3.6 % (10-15 MeV), +4.6 %
  (15-20 MeV; Griffin 1.79-1.82 at 14 MeV, HEATR 1.89).
- **Independent integration** (`lab/neutron_niel/elastic_crosscheck.py`): elastic damage
  of Si-28 (0.2-1.5 MeV) and C-12 (0.2-4 MeV) recomputed in Python from MF3 (RECONR) and
  MF4 (ENDF) with the same partition agrees with HEATR to <= 2e-4.
- **Library**: ENDF/B-VII.1 (natural C, older Si -- the VIII.0 Si files are the VII.1
  evaluations, so this is a carbon check only) changes SiC(1 MeV) by -0.09 % and r_n by
  -0.16 % (`lab/neutron_niel/library_spread.py`). JEFF-3.3 / JENDL-5 downloads were
  blocked (HTTP 403 / bot wall) and were not run.

### Systematics of r_n (not in the stored value)

- **Spectrum** (dominant): the reactor spectrum behind the introduction rates (TRIGA
  Vienna) is not available here. Watt only: 1.069 (1 keV cut), 1.067 (0.1 MeV), 1.019
  (1 MeV cut); mono-energetic 1 MeV: 1.131. Adding a 1/E slowing-down tail (1 keV-1 MeV,
  illustrative shape, not a measured spectrum) carrying 30 / 50 / 70 % of the >1 keV
  flux gives 1.10 / 1.14 / 1.20, because C recoils clear the threshold at energies where
  Si barely displaces (1-100 keV: SiC/Si = 3.4). **Plausible range 1.02-1.20.**
- **Compound partition**: HEATR stops each recoil in its own element. Stopping both Si and
  C recoils in a mean SiC lattice (Z 10, A 19.8) raises their fission-averaged elastic
  damage by +1.9 % each (+2-3 % with the other element's lattice), i.e. r_n +1.9 %.
  Whether the SR-NIEL compound proton table (r_p) uses a compound partition is not
  documented on the calculator page, so the proton/neutron consistency is uncertain at
  this ~2 % level.
- E_d and library: < 0.2 %.

### How it is used (opt-in; defaults unchanged)

`etna.core.radiation_damage.hardness_table(particle, convention)`, threaded as
`niel_table=` through `compute_damaged_params`, `cce_vs_fluence`,
`cce_vs_bias_at_fluence`, and as `particle=`, `niel_convention=` through
`run_radiation_damage`, `compensation_fluence`, `trap_compensation_fluence`,
`default_fluence_grid` and the Radiation Damage page:

- `si-neq` (default): NIEL_x,Si(E) / NIEL_n,Si(1 MeV, ASTM) -- the label of the rates.
- `sic-damage`: NIEL_x,SiC(E) / (r_n NIEL_n,Si(1 MeV)) -- the Si-labelled reactor-neq
  fluence giving the same SiC damage; for protons = kappa_Si r_p / r_n: 0.88x at 5.6 MeV,
  0.80x at 62 MeV, 0.77x at 150 MeV (compensation fluence 1.13x / 1.25x / 1.30x).
- `sic-1mev`: NIEL_x,SiC(E) / NIEL_n,SiC(1 MeV) -- SiC-native 1 MeV neq (for rates
  quoted that way; none of the shipped rates are).
- `particle="neutron"`: the bin-averaged tables above (10 eV-18.9 MeV centres),
  log-log interpolated, no extrapolation.

The default stays `si-neq` because the r_n spectrum range (1.02-1.20) is wider than the
choice it would fix is certain; `sic-damage` is the physically consistent option when
the calibration spectrum is fission-like. What stays unvalidated: the introduction rates
themselves, and NIEL scaling for low-energy protons (see above).
