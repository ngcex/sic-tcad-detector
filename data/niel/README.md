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
- **SiC: NOT AVAILABLE (blocked).** SR-NIEL's NJOY-based neutron NIEL for SiC is behind a
  registered login (the page redirects to `/index.php/component/users/login`), and no
  open table of NIEL_n,SiC(1 MeV) or of the reactor-spectrum-averaged SiC damage function
  was found (searched: SR-NIEL pages, De Napoli 2022 review, Gsponer et al.
  arXiv:2310.02047, Radmanovac et al. arXiv:2510.11304, the OSU thesis osu1173103438,
  which only plots the SiC damage function). Not fabricated.

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
0.85 (62), 0.82 (150), 0.80 (250 MeV). r_n is the blocked quantity. If r_n > r_p
(lighter C recoils take more energy per neutron collision, the usual expectation),
kappa_Si OVERSTATES the SiC damage per proton by r_n / r_p; no number is claimed for it.

## Known systematic of NIEL scaling at low energy

Measured silicon proton hardness factors are below the NIEL prediction at low energy:
3.07 +- 0.37 at 10.5 MeV and 2.19 +- 0.22 at 24.3 MeV (arXiv:2408.12410), vs SR-NIEL
kappa_Si = 4.62 at 10 MeV and 2.99 at 25 MeV here (35-50 % above). The 5.6 MeV API
default sits in this regime.
