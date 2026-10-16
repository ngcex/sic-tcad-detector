# PERF-04: iv_sweep anode current is noise at 300 K

Found by an independent review of the balanced current extraction (2026-10-17).
`iv_sweep` (default `contact="anode"`, etna/core/drift_diffusion.py) uses the direct
edge current. At 300 K, -10 V it returns +5.5e-13 A/cm2 where the balanced/interior
value is -1.49e-10 (wrong sign, ~270x too small): the hole quasi-Fermi gradient across
p+ is below double precision. notebooks/02_electrical_characterization.ipynb uses it for
the reverse I-V. Fix: use `extract_steady_state_current` (or the cathode) in iv_sweep,
then re-run notebook 02 and check tests. Transient / single-particle paths also still
use the direct form.
