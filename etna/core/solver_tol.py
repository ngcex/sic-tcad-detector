"""Shared devsim solver tolerances (one source for every steady-state / transient site).

Newton is quadratic, so on 1D devices rel 1e-8 returns bit-identical outputs to 1e-10
(measured on cv_sweep, cce_vs_bias, dark current). On 2D devices at high bias the
per-device update stagnates near 1e-9, so a 1e-10 request hit the iteration cap
and fell back to the 1e-8 retry anyway; 2D outputs change by <= ~1e-10 relative
(CCE) and 1e-14 (potential). One constant keeps every steady-state DC site in step.
"""

DC_REL_ERROR = 1e-8
DC_ABS_ERROR = 1e10

# Fallback tiers of the DC solves (tried in this order after the strict pair above).
# relaxed: abs 1e12 at the same relative error; aggressive: converged to 1e-6 relative
# (zero-flux surfaces next to p+ implant edges limit-cycle at ~5e-8, see
# drift_diffusion._RAMP_SOLVE_AGGRESSIVE).
RELAXED_ABS_ERROR = 1e12
RELAXED_REL_ERROR = 1e-8
AGGRESSIVE_ABS_ERROR = 1e14
AGGRESSIVE_REL_ERROR = 1e-6

# Transient (BDF1) steps of etna.core.single_particle: strict, then relaxed retry.
TRANSIENT_ABS_ERROR = 1e10
TRANSIENT_REL_ERROR = 1e-8
TRANSIENT_RELAXED_ABS_ERROR = 1e12
TRANSIENT_RELAXED_REL_ERROR = 1e-6

# etna.core.flash_recombination keeps the pre-2026 rel 1e-10 on its high-injection DC
# solves: not moved to DC_REL_ERROR because no parity run covers the flash continuation
# (bit-identity at 1e-8 is measured only on cv_sweep, cce_vs_bias and the dark current).
FLASH_REL_ERROR = 1e-10
