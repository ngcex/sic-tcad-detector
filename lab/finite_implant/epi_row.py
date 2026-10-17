"""1D (plane, no fringe) W(V) of a 9.75 um epi refitted on 0 / -10 V: the knee-table row.

usage: uv run python lab/finite_implant/epi_row.py > lab/finite_implant/notes/epi975_plane.out
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import recal  # noqa: E402

V = (0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0)
nj, nb = recal.solve_plane(
    {0.0: 0.0, 10.0: 0.0}, (2.958237e15, 8.252732e13), epi_um=9.75
)
c = recal.plane_C(nj, nb, volts=V, epi_um=9.75)
print(
    f"RESULT epi 9.75 um, no buffer, no fringe: N_D_junction {nj:.4e} N_D_bulk {nb:.4e}"
)
print("RESULT " + "  ".join(f"{-v:+.0f} V: {recal.EPS / c[v] * 1e4:.4f}" for v in V))
