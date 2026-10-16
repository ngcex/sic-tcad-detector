"""Is edge-vs-centre CCE a lateral-mesh artefact? Engine lateral scans on the production
graded x mesh vs uniform x meshes (same y lines, doping, BCs, stripe sigma 2 um).

usage: uv run python lab/solver2d/bench/mesh_artefact.py
The uniform grids are built directly as Grid2D arrays (tensor grid + devsim's diagonal
orientation (i,j+1)-(i+1,j) for the I_gen triangle quadrature); rebuilding the production
grid the same way reproduces devsim's element list (checked below).
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from solver2d.engine import Device2D  # noqa: E402
from solver2d.grid import Grid2D, grid_from_devsim  # noqa: E402


def synth_grid(xs, ys):
    nx, ny = len(xs), len(ys)
    k = lambda i, j: j * nx + i
    tri = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a, b, c, d = k(i, j), k(i + 1, j), k(i, j + 1), k(i + 1, j + 1)
            tri += [(a, b, c), (b, d, c)]
    return Grid2D(
        xs=np.asarray(xs, float), ys=ys, perm=np.arange(nx * ny), tri=np.array(tri)
    )


def scan(g, xpos):
    dev = Device2D(g)
    dev.dark(50.0)
    dev.factor_continuity()
    out = []
    for x in xpos:
        G = dev.generation(x)
        I = dev.solve_illuminated(G)[3]
        out.append(abs(I) / g.integrate_tri(G))
    return np.array(out)


for hw in (25.0, 50.0):
    gp = grid_from_devsim(hw)
    gs = synth_grid(gp.xs, gp.ys)
    xp = np.linspace(0, hw * 1e-4, 10)
    c_prod = scan(gp, xp)
    c_syn = scan(gs, xp)
    print(
        f"\nhw={hw:g}: production nx={gp.nx}; synth-tri reproduces production: max diff {np.abs(c_prod - c_syn).max():.1e}"
    )
    print(
        f"  (a) production : centre {c_prod[0]:.8f} edge {c_prod[-1]:.8f} edge/centre-1 {c_prod[-1] / c_prod[0] - 1:+.2e}  spread(all 10) {np.ptp(c_prod):.2e}"
    )
    for h_um in (0.5, 0.25, 0.125):
        n = int(round(hw / h_um)) + 1
        g = synth_grid(np.linspace(0, hw * 1e-4, n), gp.ys)
        xm = np.array(
            [0.0, hw * 0.5e-4, hw * 1e-4]
        )  # centre, mid (far from both sides), edge
        c = scan(g, np.concatenate([xp, xm[1:2]]))
        c10, cmid = c[:10], c[10]
        print(
            f"  uniform {h_um} um (nx={n}): centre {c10[0]:.8f} mid {cmid:.8f} edge {c10[-1]:.8f}  "
            f"edge/centre-1 {c10[-1] / c10[0] - 1:+.2e}  centre/mid-1 {c10[0] / cmid - 1:+.2e}  spread(10) {np.ptp(c10):.2e}"
        )

# Check of the O(h) remainder: same scans, I_gen by the box rule sum(NodeVolume * G)
# instead of the triangle quadrature (one diagonal orientation).
print("\nbox-rule I_gen (sum V*G) instead of triangle quadrature:")
for hw, h_um in ((50.0, None), (50.0, 0.5), (50.0, 0.25)):
    gp = grid_from_devsim(hw)
    g = gp if h_um is None else synth_grid(np.linspace(0, hw * 1e-4, int(round(hw / h_um)) + 1), gp.ys)
    dev = Device2D(g)
    dev.dark(50.0)
    dev.factor_continuity()
    vol = np.outer(g.dyd, g.dxd).ravel()
    row = []
    for x in (0.0, hw * 0.5e-4, hw * 1e-4):
        G = dev.generation(x)
        I = dev.solve_illuminated(G)[3]
        row.append(abs(I) / np.sum(vol * G))
    c0, cm, ce = row
    print(f"  hw={hw:g} mesh={'production' if h_um is None else f'uniform {h_um} um'}: centre {c0:.9f} mid {cm:.9f} edge {ce:.9f}  edge/centre-1 {ce / c0 - 1:+.2e}  centre/mid-1 {c0 / cm - 1:+.2e}")
