"""numba float64 2D drift-diffusion engine for laterally uniform 4H-SiC p+/n- diodes.

Promoted from lab/solver2d (PERF-03). Serves ``microdosimetric_sweep(solver="auto")``
only; see ``scan.lateral_cce_scan_2d``.
"""

from .scan import (
    FIELD_ENGINE_MAX,
    FIELD_NAN_MIN,
    ENGINE2D_BOX,
    Engine2DFallback,
    Engine2DPunchThrough,
    engine2d_unsupported_reason,
    lateral_cce_scan_2d,
)

__all__ = [
    "FIELD_ENGINE_MAX",
    "FIELD_NAN_MIN",
    "ENGINE2D_BOX",
    "Engine2DFallback",
    "Engine2DPunchThrough",
    "engine2d_unsupported_reason",
    "lateral_cce_scan_2d",
]
