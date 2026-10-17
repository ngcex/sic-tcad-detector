# PERF-04: iv_sweep anode current is noise at 300 K

Found by an independent review of the balanced current extraction (2026-10-17).
`iv_sweep` (default `contact="anode"`, etna/core/drift_diffusion.py) uses the direct
edge current. At 300 K, -10 V it returns +5.5e-13 A/cm2 where the balanced/interior
value is -1.49e-10 (wrong sign, ~270x too small): the hole quasi-Fermi gradient across
p+ is below double precision. notebooks/02_electrical_characterization.ipynb uses it for
the reverse I-V. Fix: use `extract_steady_state_current` (or the cathode) in iv_sweep,
then re-run notebook 02 and check tests. Transient / single-particle paths also still
use the direct form.

## Status 2026-10-17: iv_sweep part FIXED
`iv_sweep` now uses `extract_steady_state_current` (tests/test_iv_sweep_balanced.py;
notebook 02 re-executes). Dark device, -10 V, 300 K: anode before +5.5e-13, after -4.67e-11
A/cm2 (= -cathode, 1e-6). Forward 2.0 V: direct anode 1.157e-9 vs balanced 1.2005e-9
(direct cathode -1.2001e-9 agrees with balanced, so direct anode was off by 3.6 %); 2.5 V
agrees to 1e-5. Transient / single-particle paths stay on the direct form on purpose
(they have a dQ/dt displacement term; the balanced identity is steady-state only).
