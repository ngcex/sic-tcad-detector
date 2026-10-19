"""Mesh node positions of the 1D diode, taken from devsim itself.

devsim grades the mesh internally, so the exact node positions are obtained
by building the mesh-only device once per geometry (no equations, no solve)
and caching its ``x`` array. This keeps the batched engine on the identical
discretisation as the devsim path.
"""

from __future__ import annotations

import uuid
from functools import lru_cache

import devsim
import numpy as np

from etna.core.device import build_diode_mesh
from etna.core.devsim_lock import serialized


@lru_cache(maxsize=32)
@serialized
def node_positions(epi_thickness_cm, substrate_thickness_cm, nplus_thickness_cm=None):
    """Node coordinates x (cm) of the diode mesh, read-only (N,).

    ``nplus_thickness_cm`` (None = no n+ layer) is part of the cache key. Unlike
    create_sic_device (whose default is the stack's 1 um n+), the default here is None:
    a mesh probe for a default-stack device must pass the n+ thickness. The default is
    kept because etna.api.simulation's legacy reach-through branch calls
    node_positions(epi, sub) for devices that genuinely have no n+; make_batch, which
    turns the mesh into a device, requires nplus/buffer explicitly.
    """
    name = f"mesh_probe_{uuid.uuid4().hex[:8]}"
    region = "sic"
    build_diode_mesh(
        name,
        region,
        substrate_thickness_cm,
        substrate_thickness_cm + epi_thickness_cm,
        nplus_thickness_cm=nplus_thickness_cm,
    )
    devsim.create_device(mesh=f"{name}_mesh", device=name)
    try:
        x = np.array(devsim.get_node_model_values(device=name, region=region, name="x"))
    finally:
        devsim.delete_device(device=name)
        devsim.delete_mesh(mesh=f"{name}_mesh")
    x.setflags(write=False)
    return x
