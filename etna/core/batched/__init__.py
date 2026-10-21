"""Batched 1D diode engine (numba, float64): thousands of reverse-bias C-V / field /
CCE points per second, cross-checked against the devsim path (solver parity, not a
validation against measurement).

Scope: 1D graded or uniform p+/n- diode, reverse bias, Poisson-with-quasi-Fermi
electrostatics (W, C, psi) and Gummel transport for CCE. Devices in one batch
share a mesh (epi/substrate thickness); doping, N_A and T may differ per device.
devsim remains the reference implementation for everything else.
"""

from etna.core.batched.mesh import node_positions
from etna.core.batched.sweeps import (
    alpha_generation,
    cce_sweep,
    make_batch,
    reverse_sweep,
)

__all__ = ["node_positions", "make_batch", "reverse_sweep", "cce_sweep", "alpha_generation"]
