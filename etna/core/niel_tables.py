"""Proton NIEL tables (SR-NIEL) behind the radiation-damage hardness factors.

Values are transcribed verbatim from ``data/niel/sr_niel_proton_*.csv`` (SR-NIEL-7
ver. 11.1, 12 May 2026; Hadron + Coulomb, no hadronic scaling, Robinson partition);
``tests/test_niel_tables.py`` checks every row against those files. Provenance, the
form settings and the convention discussion are in ``data/niel/README.md``.

Range 1-1000 MeV (rows below 1 MeV are in the CSVs but not used: such protons stop
within a few um and do not give a uniform displacement fluence).
"""

# ASTM E722 1 MeV neutron normalisation in silicon, D_n(1 MeV) = 95 MeV mb
# (Fretwurst et al., arXiv:physics/0211118), as NIEL: 95e-27 * N_A / 28.0855 g.
NIEL_NEUTRON_1MEV_SI_MEV_CM2_G = 95e-27 * 6.02214076e23 / 28.0855

# Proton NIEL in Si, E_d = 21 eV: (energy MeV, NIEL MeV cm^2/g).
SR_NIEL_PROTON_SI = (
    (1.0, 0.063442),
    (1.5, 0.043453),
    (2.0, 0.033142),
    (2.5, 0.026827),
    (3.0, 0.022565),
    (3.5, 0.01958),
    (4.0, 0.017491),
    (4.5, 0.015958),
    (5.0, 0.014754),
    (5.5, 0.013767),
    (5.6, 0.01359),
    (6.0, 0.012939),
    (6.5, 0.012239),
    (7.0, 0.011642),
    (7.5, 0.01113),
    (8.0, 0.010687),
    (8.5, 0.010303),
    (9.0, 0.0099662),
    (9.5, 0.0096687),
    (10.0, 0.0094037),
    (12.0, 0.0085741),
    (15.0, 0.0077203),
    (20.0, 0.0067586),
    (25.0, 0.0060892),
    (30.0, 0.0055935),
    (35.0, 0.0052095),
    (40.0, 0.0048988),
    (45.0, 0.0046373),
    (50.0, 0.0044102),
    (55.0, 0.0042079),
    (60.0, 0.0040246),
    (62.0, 0.0039558),
    (65.0, 0.0038567),
    (70.0, 0.0037017),
    (75.0, 0.0035579),
    (80.0, 0.0034243),
    (85.0, 0.0032999),
    (90.0, 0.0031841),
    (95.0, 0.0030762),
    (100.0, 0.0029757),
    (120.0, 0.0026379),
    (150.0, 0.0022803),
    (200.0, 0.0019291),
    (250.0, 0.0017408),
    (300.0, 0.0016326),
    (350.0, 0.0015658),
    (400.0, 0.0015217),
    (450.0, 0.0014906),
    (500.0, 0.0014677),
    (550.0, 0.0014505),
    (600.0, 0.0014375),
    (650.0, 0.0014277),
    (700.0, 0.0014203),
    (750.0, 0.0014142),
    (800.0, 0.0014085),
    (850.0, 0.0014017),
    (900.0, 0.0013926),
    (950.0, 0.0013794),
    (1000.0, 0.0013607),
)

# Proton NIEL in SiC (compound (IV) SiC, Si:C = 1:1), E_d = 21 / 21 eV.
SR_NIEL_PROTON_SIC = (
    (1.0, 0.061349),
    (1.5, 0.041795),
    (2.0, 0.031769),
    (2.5, 0.025655),
    (3.0, 0.021539),
    (3.5, 0.018683),
    (4.0, 0.01669),
    (4.5, 0.015181),
    (5.0, 0.013978),
    (5.5, 0.012996),
    (5.6, 0.012821),
    (6.0, 0.012181),
    (6.5, 0.011498),
    (7.0, 0.010918),
    (7.5, 0.01042),
    (8.0, 0.009987),
    (8.5, 0.0096069),
    (9.0, 0.0092693),
    (9.5, 0.0089668),
    (10.0, 0.0086937),
    (12.0, 0.0078188),
    (15.0, 0.006918),
    (20.0, 0.0059718),
    (25.0, 0.0053688),
    (30.0, 0.0049368),
    (35.0, 0.0045951),
    (40.0, 0.0043062),
    (45.0, 0.0040533),
    (50.0, 0.0038286),
    (55.0, 0.0036275),
    (60.0, 0.0034466),
    (62.0, 0.0033792),
    (65.0, 0.0032831),
    (70.0, 0.0031349),
    (75.0, 0.0029998),
    (80.0, 0.0028763),
    (85.0, 0.0027631),
    (90.0, 0.0026589),
    (95.0, 0.0025629),
    (100.0, 0.0024743),
    (120.0, 0.00218),
    (150.0, 0.0018702),
    (200.0, 0.001563),
    (250.0, 0.0013968),
    (300.0, 0.0013027),
    (350.0, 0.0012469),
    (400.0, 0.0012121),
    (450.0, 0.0011891),
    (500.0, 0.0011729),
    (550.0, 0.0011606),
    (600.0, 0.001151),
    (650.0, 0.0011434),
    (700.0, 0.0011375),
    (750.0, 0.001133),
    (800.0, 0.0011292),
    (850.0, 0.0011251),
    (900.0, 0.0011191),
    (950.0, 0.0011092),
    (1000.0, 0.0010927),
)


# ---------------------------------------------------------------------------
# Neutron NIEL (issue: SiC neutron 1-MeV normalisation), NJOY2016.78 HEATR on
# ENDF/B-VIII.0, Lindhard/Robinson partition, E_d 21 eV (Si) / 21 eV (C); computed by
# scripts/build_neutron_niel.py, stored in data/niel/neutron_*.csv(.gz).
# tests/test_neutron_niel.py recomputes every constant below from those files.
# Checks (data/niel/README.md): Si pointwise 1 MeV = 0.998 x the ASTM E722 95 MeV mb;
# Si Watt-fission average within 0.7 % of the Griffin (SAND92-0094) damage function;
# an independent Python elastic integration reproduces HEATR to 1e-4; E_d 35/21 eV
# changes nothing at the 5th digit; ENDF/B-VII.1 changes r_n by -0.16 %.
#
# Spectrum average: Watt fission spectrum (a 0.988 MeV, b 2.249 /MeV), 1 keV-20 MeV.

# Si, pointwise at 1.000 MeV (cf. ASTM 95 MeV mb = NIEL_NEUTRON_1MEV_SI_MEV_CM2_G)
NIEL_NEUTRON_1MEV_SI_NJOY_MEV_CM2_G = 2.03224e-3
# 4H-SiC (1:1, natural abundance), pointwise at 1.000 MeV
NIEL_NEUTRON_1MEV_SIC_MEV_CM2_G = 2.29889e-3
# Watt-fission-spectrum averages
NIEL_NEUTRON_FISSION_SI_MEV_CM2_G = 2.18058e-3
NIEL_NEUTRON_FISSION_SIC_MEV_CM2_G = 2.33030e-3
# SiC/Si neutron NIEL ratio of the fission spectrum (the "r_n" of data/niel/README.md)
R_N_SIC_SI_FISSION = NIEL_NEUTRON_FISSION_SIC_MEV_CM2_G / NIEL_NEUTRON_FISSION_SI_MEV_CM2_G
