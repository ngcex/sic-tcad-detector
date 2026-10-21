# NIST ESTAR electron stopping powers (real data)

Downloaded **2026-10-21** from the NIST STAR text interface (M.J. Berger, J.S. Coursey,
M.A. Zucker, J. Chang, _ESTAR, PSTAR and ASTAR_, NIST,
https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html). Nothing here was typed by
hand. Used by `etna/core/electron_range.py` (delta-ray ranges for
`etna/core/mc_transport.py`, issue #7).

## Raw responses (`raw/`, unmodified HTTP bodies)

Request: `POST https://physics.nist.gov/cgi-bin/Star/e_table-t.pl`, form fields
`matno=<NNN>`, `ShowDefault=on`, and for the `lowE` files also
`Energies=0.001 0.0015 0.002 0.003 0.004 0.005 0.006 0.008` (MeV, one per line).

| file                         | content                                                                                                              |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| `estar_lowE_014.html`        | Si: T, collision, radiative, total stopping power (MeV cm2/g), density-effect delta; 1 keV - 1 GeV (no range column) |
| `estar_lowE_006.html`        | amorphous carbon (rho 2.0), same columns                                                                             |
| `estar_lowE_276.html`        | liquid water, same columns                                                                                           |
| `estar_default_<matno>.html` | default grid 10 keV - 1 GeV with NIST's CSDA range column (used only to test the range integration)                  |

SHA-256:

```
06f71b25cb3d2974f7a33a2ad25be47f1b8c017246799fc4915caee8f110f412  estar_default_014.html
19325a71ec9a783745ab33eb08f4914e8092392366ff7065027509b12f197c42  estar_default_006.html
429b1ae686498176c06508fc7c793cb74d7e1b6339443a1e26b868d15bce6dbe  estar_default_276.html
1a90decb293ebe7598f55b059ba1c282e9d0b130d51abd9ede5d7b6ef35c5370  estar_lowE_014.html
c802e47a5699f0d9809cf8577bf6e8e7ed996177965b91be8c7773da52029ab6  estar_lowE_006.html
e3643d7d3183e6a2cb58c895eab65ff04460bdbdc92ae1fbac963e38b4d3a7a9  estar_lowE_276.html
```

## Derived quantities

- CSDA range = integral of 1/S_total from 1 keV, started with R(1 keV) = T/(S(1+p))
  (power law S ~ T^-p below 1 keV, assumption). It reproduces NIST's own CSDA column
  within 0.2 % above 20 keV and 0.7 % at 10 keV for Si, C and water (test).
- SiC: Bragg additivity of the Si and amorphous-C total stopping powers (mass
  fractions 0.7004 / 0.2996). Electron range in SiC (rho 3.21): 0.024 um at 1 keV,
  1.02 um at 10 keV, 54.6 um at 100 keV.
