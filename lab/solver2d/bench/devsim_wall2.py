"""devsim 2D ramp wall vs epi: last converged cathode bias in 0.5 V steps, plus the engine's
dark cathode-cell field at the last converged and first failing bias.

usage: uv run python lab/solver2d/bench/devsim_wall2.py HW EPI_UM V_START V_END ND [ND ...]
"""

import sys
import time

from etna.core.charge_collection_2d import create_2d_dd_device
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.drift_diffusion import ramp_bias
from etna.core.engine2d.device import Device2D
from etna.core.engine2d.grid import grid_from_devsim

hw, epi = float(sys.argv[1]), float(sys.argv[2]) * 1e-4
v0, v1 = float(sys.argv[3]), float(sys.argv[4])
for nd in map(float, sys.argv[5:]):
    t = time.perf_counter()
    V = v0
    try:
        info = create_2d_dd_device(
            half_width_um=hw, V_bias=v0, epi_thickness_cm=epi, N_D_bulk=nd
        )
        while V < v1 - 1e-9:
            ramp_bias(info, V_target=V + 0.5, contact="cathode", V_step=0.5)
            V += 0.5
        wall = None
    except Exception:
        wall = V + 0.5 if V > v0 or "info" in dir() else v0
    reset_devsim_fully(preserve_solver=True)
    g, ginfo = grid_from_devsim(half_width_um=hw, epi_thickness_cm=epi, N_D_bulk=nd)
    f = lambda v: Device2D(g, ginfo).dark(v)
    msg = f"WALL hw={hw:g} epi={epi * 1e4:g}um nd={nd:.2e}: last converged {V:.1f} V field {f(V):.4g}"
    if wall is not None:
        msg += f"; failed at {V + 0.5:.1f} V field {f(V + 0.5):.4g}"
    else:
        msg += f" (reached V_END)"
    print(msg + f" ({time.perf_counter() - t:.0f}s)", flush=True)
    info = None
