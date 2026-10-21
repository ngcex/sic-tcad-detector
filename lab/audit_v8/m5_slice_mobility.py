import inspect
import numpy as np
from etna.core.track import cross_section as cs
from etna.core.track.tilted import tilted_alpha_cce
from etna.core.device import local_mobility
from etna.core.sic_material import SiC4H_Parameters

print(
    "CrossSectionDevice has mobility field:",
    "mobility" in cs.CrossSectionDevice.__dataclass_fields__,
)
print(
    "tilted_alpha_cce source mentions mobility:",
    "mobility" in inspect.getsource(tilted_alpha_cce),
)
for name in ("c_axis", "legacy_perp"):
    p = SiC4H_Parameters(mobility=name)
    print(
        name,
        [
            float(local_mobility(1e14, 0.0, 1.0, 300.0, c, p))
            for c in ("electron", "hole")
        ],
    )
print("slice PARAMS.mobility:", cs.PARAMS.mobility)
