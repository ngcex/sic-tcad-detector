"""devsim DD build time (create_2d_dd_device hw 150, -50 V, default stack) per mesh candidate.
usage: uv run python lab/solver2d/bench/nplus_mesh_devsim.py CAND [CAND ...]"""
import sys
import time
import uuid
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import devsim
from nplus_mesh import CANDIDATES, mesh_variant

from etna.core.charge_collection_2d import create_2d_dd_device
from etna.core.device import DEFAULT_STACK

epi_end = DEFAULT_STACK["substrate_thickness_cm"] + 10e-4
for cand in sys.argv[1:]:
    name = f"nm_{uuid.uuid4().hex[:6]}"
    t0 = time.perf_counter()
    with mesh_variant(CANDIDATES[cand], epi_end):
        info = create_2d_dd_device(half_width_um=150.0, V_bias=50.0, device_name=name)
    t = time.perf_counter() - t0
    print(f"RESULT {cand}: num_nodes {info['num_nodes']} build {t:.1f} s", flush=True)
    devsim.delete_device(device=name)
