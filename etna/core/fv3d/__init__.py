"""fv3d: 3D Poisson (frozen flat quasi-Fermi levels, Boltzmann carriers) on tensor grids.

Matrix-free box discretisation (devsim / fv2d), Newton with an energy line search,
geometric-multigrid-preconditioned CG (numba float64), Ramo weighting potentials, and
analyses for 3D-column lattices and square-cell microdosimeters (issue #4: the
axisymmetric proxy is not valid for column lattices). See docs/USAGE.md section 11.2.
"""

from .analysis import depletion_scan, depletion_voltage, ramo_signal
from .mg import Multigrid, pcg
from .solver import FV3DNotConverged, Poisson3D, Solution3D
from .structures import (
    Electrode,
    Structure3D,
    coax_cell,
    column_cell,
    extrude_2d,
    slab,
    trench_cell,
)

__all__ = [
    "Electrode",
    "FV3DNotConverged",
    "Multigrid",
    "Poisson3D",
    "Solution3D",
    "Structure3D",
    "coax_cell",
    "column_cell",
    "depletion_scan",
    "depletion_voltage",
    "extrude_2d",
    "pcg",
    "ramo_signal",
    "slab",
    "trench_cell",
]
