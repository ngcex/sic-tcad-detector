"""Constants and doping for the batched solver, taken from the devsim path.

Material functions (n_i, mobility, lifetimes, acceptor ionization) are the
``etna.core`` ones, so there is a single source of truth. Only the graded
doping profile is restated in numpy, as the same expression that
``etna.core.device.set_graded_doping_profile`` hands to devsim.
"""

from __future__ import annotations

import numpy as np

from etna.core.device import EPS_0, K_B_J, Q
from etna.core.sic_material import SiC4H_Parameters

PARAMS = SiC4H_Parameters()
EPS = PARAMS.eps_r * EPS_0  # F/cm


def thermal_voltage(T):
    return K_B_J * T / Q


def graded_doping(x, junction_pos, N_A_ion, N_D_junction, N_D_bulk, L_transition):
    """Donors and Acceptors, each (B, N), for ETNA's exponentially graded epi profile.

    Scalar or (B,) parameters broadcast over the batch; ``x`` is shared (N,).
    Uniform doping is the special case N_D_junction == N_D_bulk. As in devsim,
    ``step(0) == 1``: the junction node carries both donors and acceptors.
    """
    x = np.asarray(x, dtype=float)[None, :]
    if np.any(np.asarray(L_transition) <= 0):
        raise ValueError("L_transition must be positive")
    col = lambda a: np.atleast_1d(np.asarray(a, dtype=float))[:, None]
    NA, NJ, NB, Lt = col(N_A_ion), col(N_D_junction), col(N_D_bulk), col(L_transition)
    n_side = (x >= junction_pos).astype(float)
    p_side = (junction_pos - x >= 0).astype(float)
    donors = (NB + (NJ - NB) * np.exp(-np.maximum(x - junction_pos, 0.0) / Lt)) * n_side
    acceptors = NA * p_side
    B = max(donors.shape[0], acceptors.shape[0])
    shape = (B, x.shape[1])
    return np.broadcast_to(donors, shape).copy(), np.broadcast_to(acceptors, shape).copy()
