"""General (laterally non-uniform) 2D drift-diffusion engine, Cartesian or axisymmetric.

numba float64, tensor-grid box method with Scharfetter-Gummel currents, banded LU with the
shorter dimension fastest. See docs/USAGE.md "General 2D engine (fv2d)" and
lab/geometry2d/notes/LOG.md. Opt-in: nothing in the default paths calls it.
"""

from .devsim_import import NotTensorMesh, mesh_from_builder, mesh_from_devsim
from .mesh import FVMesh, SolverOrder, box_geometry
from .solver import FV2DNotConverged, Solver2D, State
from .workflows import (
    alpha_stripe_generation_fv2d,
    cce_lateral_scan_fv2d,
    contact_sharing_fv2d,
    cut_node_sets,
    prepared_solver,
)

__all__ = [
    "FVMesh",
    "FV2DNotConverged",
    "NotTensorMesh",
    "Solver2D",
    "SolverOrder",
    "State",
    "alpha_stripe_generation_fv2d",
    "box_geometry",
    "cce_lateral_scan_fv2d",
    "contact_sharing_fv2d",
    "cut_node_sets",
    "mesh_from_builder",
    "mesh_from_devsim",
    "prepared_solver",
]
