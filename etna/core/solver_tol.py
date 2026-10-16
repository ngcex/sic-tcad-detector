"""Shared devsim DC tolerance.

Newton is quadratic, so on 1D devices rel 1e-8 returns bit-identical outputs to 1e-10
(measured on cv_sweep, cce_vs_bias, dark current). On 2D devices at high bias the
per-device update stagnates near 1e-9, so a 1e-10 request hit the iteration cap
and fell back to the 1e-8 retry anyway; 2D outputs change by <= ~1e-10 relative
(CCE) and 1e-14 (potential). One constant keeps every steady-state DC site in step.
"""

DC_REL_ERROR = 1e-8
