# PERF-06: notebook speed (2026-10-17)

nb12 1957 -> 41 s, nb10 89 -> 14 s, nb13 150 -> 16 s (fluence sweeps via the batched engine,
deviation ~3e-15). nb19 1373 -> 737 s (failed devices were left loaded and re-solved by every
later cell; fixed). nb05 (146 s) stays: dark_current_vs_fluence needs TAT/devsim (9.3 s/call).

Open notebook BUGS (not fixed: they change results, need a physics decision):

- nb19: FIXED 2026-10-18 (737 -> ~200 s, all cells run). Every structure ramps the cathode to
  +50 V; CCE via cce_with_generation / contact_current_sharing (inject -> solve -> balanced
  N-contact extraction -> restore); guard ring rebuilt (the p+ ran under the whole width: a
  short), delta-E/E rebuilt as two independent diodes, 3D anode contact got its wall nodes;
  lateral p+ edges need the aggressive solver tier (limit cycle, accepted via solve(info=True),
  verified by the contact-current residual ~1e-15); the no-ring reference diverges above
  ~11.5 V so the with/without comparison is at 10 V. Contact-node normalisation FIXED
  globally 2026-10-18: compute_cce_2d / contact_current_sharing normalise by the injected
  (non-contact) generation, alpha_stripe_generation / cce_lateral_scan default to
  zero_on_contacts=True, engine2d zeroes its contact rows (parity kept, ~1e-14); nb19 no
  longer passes the flag. Legacy stack only.
- nb13: reads df['epi_thickness_cm'] but radiation_hardness_sweep returns 'epi_um' (cells 4-7 error).
