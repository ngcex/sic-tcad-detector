"""One transverse plasma slice: default (c_axis) vs legacy_perp slice mobility.
The slice ignores DeviceConfig.mobility; we patch the module PARAMS to show the effect."""

import time
from etna.core.track import cross_section as cs
from etna.core.track.radial import thermalized, whk  # noqa
from etna.core.sic_material import SiC4H_Parameters

E, LET, beta, s_nm = 2e4, 150.0, 0.05, 10.0
rad = thermalized(whk(beta), s_nm * 1e-7)
out = {}
for name in ("c_axis", "legacy_perp"):
    cs.PARAMS = SiC4H_Parameters(mobility=name)
    t = time.time()
    r = cs.slice_loss(E, LET, rad)
    out[name] = r.recombined_fraction
    print(
        name,
        "loss",
        r.recombined_fraction,
        "conv",
        r.converged,
        f"{time.time()-t:.1f}s",
    )
print("rel diff", out["legacy_perp"] / out["c_axis"] - 1)
