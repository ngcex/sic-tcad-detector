# PERF-06: notebook speed (2026-10-17)
nb12 1957 -> 41 s, nb10 89 -> 14 s, nb13 150 -> 16 s (fluence sweeps via the batched engine,
deviation ~3e-15). nb19 1373 -> 737 s (failed devices were left loaded and re-solved by every
later cell; fixed). nb05 (146 s) stays: dark_current_vs_fluence needs TAT/devsim (9.3 s/call).

Open notebook BUGS (not fixed: they change results, need a physics decision):
- nb19: mesa/3D/guard-ring cells ramp the ANODE (p+) to +50 V = forward bias (planar ramps the
  cathode = reverse). Flipping the sign makes mesa converge (CCE 0.9987) but guard ring fails
  ('balanced extraction needs exactly anode+cathode contacts', third contact) and 3D fails at -25 V.
  delta-E cell calls ramp_voltage(V_target=...) (signature is now V_start): TypeError.
  3D radial / guard-ring edge / delta-E call compute_cce_2d without injecting generation: values
  meaningless. These structures need devsim 2D on non-uniform meshes (engine2d is uniform only).
- nb13: reads df['epi_thickness_cm'] but radiation_hardness_sweep returns 'epi_um' (cells 4-7 error).
