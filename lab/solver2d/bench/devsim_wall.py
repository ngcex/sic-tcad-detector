"""Where does devsim's 2D ramp fail (hw, N_D_bulk)? Ramp from 50 V in 0.5 V fixed steps.
usage: uv run python lab/solver2d/bench/devsim_wall.py HW ND [ND ...]"""
import sys, time
import devsim
from etna.core.charge_collection_2d import create_2d_dd_device
from etna.core.drift_diffusion import ramp_bias
from etna.core.devsim_reset import reset_devsim_fully
hw = float(sys.argv[1])
for nd in map(float, sys.argv[2:]):
    t = time.perf_counter()
    info = create_2d_dd_device(half_width_um=hw, V_bias=50.0, N_D_bulk=nd)
    V = 50.0
    try:
        while V < 80:
            ramp_bias(info, V_target=V + 0.5, contact="cathode", V_step=0.5)
            V += 0.5
    except Exception as e:
        pass
    print(f"WALL hw={hw:g} nd={nd:.2e}: last converged {V:.1f} V ({time.perf_counter()-t:.1f}s)", flush=True)
    reset_devsim_fully(preserve_solver=True)
