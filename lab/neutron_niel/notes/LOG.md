# lab/neutron_niel LOG

2026-10-21. SiC neutron NIEL normalisation (the open item of data/niel/README.md).

## Tools (scratch, not ETNA dependencies)

    conda create -p <scratch>/njoy/env -c conda-forge --override-channels njoy2016   # 2016.78
    curl -O https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_neutrons.zip       # 310 MB
    curl -O https://www.nndc.bnl.gov/endf-b7.1/zips/ENDF-B-VII.1-neutrons.zip        # 226 MB
    # JEFF-3.3 (oecd-nea.org) 403, JENDL-5 (jaea) bot wall, IAEA download-endf 403: not run

## Build (29 s for 8 NJOY runs; byte-reproducible, checked twice)

    uv run python scripts/build_neutron_niel.py --njoy <njoy> --endf-dir <b80> --work <w>

## Checks

    uv run python lab/neutron_niel/bench.py                         # 0.05 s, stored tables
    uv run python lab/neutron_niel/elastic_crosscheck.py <b80> <w>  # 7 s, needs ENDF + PENDF
    uv run python lab/neutron_niel/library_spread.py <njoy> <b71> <w>

- Si HEATR pointwise 1 MeV / ASTM 95 MeV mb: 0.9977.
- Watt-fission kappa_Si: HEATR 1.0705, Griffin (RD48 neutrons.xls) 1.0628 (+0.72 %).
  Per bin HEATR/Griffin: 0.1-0.5 1.013, 0.5-1 1.017, 1-2 1.001, 2-5 1.005, 5-10 1.022,
  10-15 1.036, 15-19.9 1.046.
- Python elastic (MF3 from RECONR PENDF, MF4 LTT=3 from ENDF, same Robinson partition)
  vs HEATR MT444: Si-28 0.2-1.5 MeV ratio 0.9998-1.0002; C-12 0.2-4 MeV 1.0000.
- Partition lattice (fission-averaged elastic damage): Si-28 in a mean SiC lattice
  (Z 10, A 19.8) x1.0189, in C x1.0305; C-12 in the mean lattice x1.0191, in Si x1.0209.
- E_d 35/21 vs 21/21: identical to 5 digits at 0.1, 1, 14 MeV and in the fission average.
- ENDF/B-VII.1 (C-nat): SiC(1 MeV) x0.9991, SiC fission x0.9984, r_n x0.9984. Si
  identical (VIII.0 Si = VII.1 evaluations).

## Results

- NIEL_n,Si(1 MeV) 2.03224e-3, NIEL_n,SiC(1 MeV) 2.29889e-3 MeV cm2/g (76.5 MeV mb/atom).
- Watt fission averages: Si 2.18058e-3, SiC 2.33030e-3; r_n = 1.0687.
- r_n spectrum variants: 0.1-20 MeV 1.0670; 1-20 MeV 1.0194; mono 1 MeV 1.1312;
  1-100 keV 1/E part alone 3.381; Watt + 1/E tail (30/50/70 % of flux) 1.101/1.137/1.200.
- sic-damage / si-neq proton kappa: 0.883 (5.6 MeV), 0.826 (30), 0.799 (62), 0.767 (150),
  0.751 (250).

## Performance

Everything is table look-up (no solver): bench 0.05 s; no numba/GPU work is warranted.
