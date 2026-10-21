"""Re-derive the tests/test_contact_mesh.py references (mesa centre CCE, uniform-bisection
limit of the graded strip; legacy-mesh values) for a low-field mobility set.

    uv run python lab/mobility_caxis/bench/contact_mesh_refs.py [c_axis|legacy_perp]
"""

import sys

import numpy as np

from etna.core.alternative_structures import create_mesa_device
from etna.core.fv2d import cce_lateral_scan_fv2d, mesh_from_builder, prepared_solver
from etna.core.fv2d.refine import refine_lines
from etna.core.fv2d.structures import build_mesh, mesa_structure
from etna.core.sic_material import SiC4H_Parameters

prm = SiC4H_Parameters(mobility=sys.argv[1] if len(sys.argv) > 1 else "c_axis")


def centre(xs, ys, V):
    m = build_mesh(mesa_structure(), xs[:3], ys, params=prm)
    s = prepared_solver(m, {"cathode": V}, dark_check=False)
    return cce_lateral_scan_fv2d(s, x_positions_cm=[0.0])["cce_values"][0]


meshes = {
    dm: mesh_from_builder(create_mesa_device, depth_mesh=dm)[0]
    for dm in ("legacy", "graded")
}
L = meshes["legacy"]
print(
    "legacy mesh", prm.mobility, [f"{centre(L.xs, L.ys, V):.8f}" for V in (10.0, 50.0)]
)
G = meshes["graded"]
for V in (10.0, 50.0):
    y, vals = G.ys, []
    for _ in range(6):
        vals.append(centre(G.xs, y, V))
        y = refine_lines(y, np.ones(len(y) - 1, bool))
    d = np.diff(vals)
    print(
        f"V {V}: graded {vals[0]:.9f}, limit {vals[-1] + d[-1]:.9f}, ratio {d[-1] / d[-2]:.3f}"
    )
