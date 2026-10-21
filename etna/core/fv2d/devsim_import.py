"""Import a devsim 2D tensor-mesh device into an :class:`~etna.core.fv2d.mesh.FVMesh`
VERBATIM (parity mode): devsim's node volumes (NodeVolume, or CylindricalNodeVolume on
an axisymmetric device), edge couplings (EdgeCouple / CylindricalEdgeCouple) and lengths,
the evaluated Donors / Acceptors / NetDoping, the edge mobilities ``mu_n`` / ``mu_p``
(add_local_mobility_models) and the node set of every contact of the region.

Asserted, not assumed: the region's nodes form a full tensor grid, every edge is a
lateral or a depth edge or a diagonal whose coupling is negligible (devsim splits each
rectangle into two right triangles; the hypotenuse carries ~1e-33 (Cartesian) or
~1e-13 relative (cylindrical) coupling), and node volumes are positive. Violations
raise :class:`NotTensorMesh`.

Works for every builder of etna.core.alternative_structures except the two-region
delta-E/E telescope (import one layer region at a time is not wired) and for the planar
create_sic_2d_device.
"""

from __future__ import annotations

import copy
import uuid

import numpy as np

from etna.core.devsim_lock import serialized

from .mesh import FVMesh

DIAG_COUPLE_MAX = 1e-9  # max |diagonal coupling| / median axis coupling accepted


class NotTensorMesh(ValueError):
    """The devsim region is not representable by the tensor-grid 5-point engine."""


def _node(device, region, name):
    import devsim

    return np.array(
        devsim.get_node_model_values(device=device, region=region, name=name)
    )


def _edge(device, region, name):
    import devsim

    return np.array(
        devsim.get_edge_model_values(device=device, region=region, name=name)
    )


def mesh_from_devsim(device_info):
    """FVMesh of the live devsim device ``device_info`` (mesh-only or solved).

    Returns (mesh, perm) with ``perm[devsim_node] = j * nx + i`` (grid index, row-major
    over (ny, nx)), so a devsim node array v maps to the grid as
    ``g = np.empty(N); g[perm] = v; g.reshape(ny, nx)``.
    """
    import devsim

    from etna.core.drift_diffusion import _contact_node_indices, contacts_on_region

    d, r = device_info["device_name"], device_info["region_name"]
    cyl = device_info.get("coordinate_system") == "cylindrical"
    x, y = _node(d, r, "x"), _node(d, r, "y")
    xs, ys = np.unique(x), np.unique(y)
    nx, ny = len(xs), len(ys)
    i = np.searchsorted(xs, x)
    j = np.searchsorted(ys, y)
    perm = j * nx + i
    if len(x) != nx * ny or len(np.unique(perm)) != len(x):
        raise NotTensorMesh("region nodes are not a full tensor-product grid")
    if nx < 3 or ny < 3:
        raise NotTensorMesh(f"degenerate grid {nx} x {ny}")

    def to_grid(v):
        g = np.empty(nx * ny)
        g[perm] = v
        return g.reshape(ny, nx)

    vol = _node(d, r, "CylindricalNodeVolume" if cyl else "NodeVolume")
    if np.any(vol <= 0):
        raise NotTensorMesh("non-positive node volume")
    names = devsim.get_edge_model_list(device=d, region=r)
    if "node_index@n0" not in names:
        devsim.edge_from_node_model(device=d, region=r, node_model="node_index")
    n0 = _edge(d, r, "node_index@n0").astype(int)
    n1 = _edge(d, r, "node_index@n1").astype(int)
    couple = _edge(d, r, "CylindricalEdgeCouple" if cyl else "EdgeCouple")
    length = _edge(d, r, "EdgeLength")
    mun = _edge(d, r, "mu_n")
    mup = _edge(d, r, "mu_p")
    a, b = perm[n0], perm[n1]
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    jl, il = lo // nx, lo % nx
    jh, ih = hi // nx, hi % nx
    lateral = (jl == jh) & (ih == il + 1)
    depth = (il == ih) & (jh == jl + 1)
    other = ~(lateral | depth)
    axis = np.abs(couple[~other])
    if other.any():
        rel = np.max(np.abs(couple[other])) / np.median(axis)
        if rel > DIAG_COUPLE_MAX:
            raise NotTensorMesh(f"diagonal edge coupling {rel:.2e} of the axis value")
    gh = np.zeros((ny, nx - 1))
    gv = np.zeros((ny - 1, nx))
    mh = [np.zeros((ny, nx - 1)) for _ in range(2)]
    mv = [np.zeros((ny - 1, nx)) for _ in range(2)]
    g = couple / length
    gh[jl[lateral], il[lateral]] = g[lateral]
    gv[jl[depth], il[depth]] = g[depth]
    for k, mu in enumerate((mun, mup)):
        mh[k][jl[lateral], il[lateral]] = mu[lateral]
        mv[k][jl[depth], il[depth]] = mu[depth]
    if lateral.sum() != ny * (nx - 1) or depth.sum() != (ny - 1) * nx:
        raise NotTensorMesh("missing lateral / depth edges")
    contacts = {}
    for c in contacts_on_region(d, r):
        mask = np.zeros(nx * ny, dtype=bool)
        mask[perm[sorted(_contact_node_indices(d, r, c))]] = True
        contacts[c] = mask.reshape(ny, nx)
    info = {
        k: copy.deepcopy(v)
        for k, v in device_info.items()
        if k not in ("params",) and not callable(v)
    }
    info["params"] = device_info.get("params")
    mesh = FVMesh(
        xs=xs,
        ys=ys,
        vol=to_grid(vol),
        gh=gh,
        gv=gv,
        donors=to_grid(_node(d, r, "Donors")),
        acceptors=to_grid(_node(d, r, "Acceptors")),
        net=to_grid(_node(d, r, "NetDoping")),
        mun_h=mh[0],
        mun_v=mv[0],
        mup_h=mh[1],
        mup_v=mv[1],
        contacts=contacts,
        coordinate_system="cylindrical" if cyl else "cartesian",
        info=info,
    )
    return mesh, perm


@serialized
def mesh_from_builder(builder, **kwargs):
    """Mesh-only devsim build ``builder(device_name=..., **kwargs)`` -> (FVMesh, perm).

    The device (and its mesh) are deleted afterwards and, for an axisymmetric builder,
    devsim's Cartesian geometry models are restored (``restore_cartesian_coords``).
    """
    import devsim

    from etna.core.alternative_structures import restore_cartesian_coords

    name = f"fv2d_{uuid.uuid4().hex[:8]}"
    info = None
    try:
        info = builder(device_name=name, **kwargs)
        return mesh_from_devsim(info)
    finally:
        if name in devsim.get_device_list():
            devsim.delete_device(device=name)
        if f"{name}_mesh" in devsim.get_mesh_list():
            devsim.delete_mesh(mesh=f"{name}_mesh")
        if info is not None and info.get("coordinate_system") == "cylindrical":
            restore_cartesian_coords()
