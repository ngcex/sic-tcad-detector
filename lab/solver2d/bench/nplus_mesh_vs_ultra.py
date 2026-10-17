"""fine / k3 / e05n02 vs the 'ultra' mesh on the 21-config gate (engine2d): is k3 within the
fine mesh's own discretisation error?  usage: uv run python lab/solver2d/bench/nplus_mesh_vs_ultra.py"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
from nplus_mesh import GATE, build, solve

print("hw epi np V | relW fine k3 e05n02 | relCCE fine k3 e05n02")
for hw, epi, npu, V in GATE:
    if hw == 150:
        continue  # laterally identical to hw 25/50 in y; ultra at hw150 is slow
    kw = dict(half_width_um=hw, epi_thickness_cm=epi * 1e-4, nplus_thickness_cm=npu * 1e-4)
    gu, iu, _ = build("ultra", **kw)
    su = solve(gu, iu, V)
    rw, rc = [], []
    for c in ("fine", "k3", "e05n02"):
        g, i, _ = build(c, **kw)
        s = solve(g, i, V)
        rw.append(f"{(s['W'] - su['W']) / su['W']:+.1e}")
        rc.append(f"{np.max(np.abs(s['cce'] / su['cce'] - 1)):.1e}")
    print(hw, epi, npu, V, "|", *rw, "|", *rc, flush=True)
