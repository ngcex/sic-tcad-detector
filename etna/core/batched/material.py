"""Constants and doping for the batched solver, taken from the devsim path.

Material functions (n_i, mobility, lifetimes, acceptor ionization) are the
``etna.core`` ones, so there is a single source of truth. Only the graded
doping profile is restated in numpy, as the same expression that
``etna.core.device.set_graded_doping_profile`` hands to devsim.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from etna.core.device import (
    EPS_0,
    K_B_J,
    Q,
    acceptor_total_scale,
    donor_profile,
    edge_mobility,
    local_mobility,
)
from etna.core.incomplete_ionization import ionized_acceptor_concentration
from etna.core.sic_material import SiC4H_Parameters

PARAMS = SiC4H_Parameters()
EPS = PARAMS.eps_r * EPS_0  # F/cm


def thermal_voltage(T):
    return K_B_J * T / Q


@lru_cache(maxsize=4096)
def _acceptor_scale(N_A, T):
    """N_A / N_A_ionized (cached: the ionization root solve is per call)."""
    return acceptor_total_scale(N_A, ionized_acceptor_concentration(N_A, T))


def edge_mobilities(donors, acceptors, N_A, T):
    """Edge mobilities (mu_n, mu_p), each (N-1, B), of the devsim builders' local model.

    ``donors`` / ``acceptors`` (B, N) are the doping the solver sees (acceptors
    IONIZED, as devsim's Acceptors model); ``N_A`` and ``T`` scalar or (B,). Node
    mobility = Caughey-Thomas at Donors + (N_A / N_A_ionized) * Acceptors
    (etna.core.device.local_mobility), arithmetic edge average (edge_mobility): the
    numpy twin of etna.core.device.add_local_mobility_models.
    """
    donors = np.atleast_2d(np.asarray(donors, dtype=float))
    acceptors = np.broadcast_to(np.asarray(acceptors, dtype=float), donors.shape)
    B = donors.shape[0]
    N_A = np.broadcast_to(np.asarray(N_A, dtype=float), (B,))
    T = np.broadcast_to(np.asarray(T, dtype=float), (B,))
    scales = [_acceptor_scale(float(N_A[b]), float(T[b])) for b in range(B)]
    out = []
    for carrier in ("electron", "hole"):
        mu = np.empty((B, donors.shape[1] - 1))
        for b in range(B):
            scale = scales[b]
            node = local_mobility(
                donors[b], acceptors[b], scale, float(T[b]), carrier, PARAMS
            )
            mu[b] = edge_mobility(node)
        out.append(np.ascontiguousarray(mu.T))
    return out[0], out[1]


def graded_doping(
    x,
    junction_pos,
    N_A_ion,
    N_D_junction,
    N_D_bulk,
    L_transition,
    nplus=None,
    buffer=None,
):
    """Donors and Acceptors, each (B, N), for ETNA's exponentially graded epi profile.

    Scalar or (B,) parameters broadcast over the batch; ``x`` is shared (N,).
    Uniform doping is the special case N_D_junction == N_D_bulk. As in devsim,
    ``step(0) == 1``: the junction node carries both donors and acceptors.
    ``nplus`` = (x_epi_end, N_D_nplus) adds the n+ layer N_D_nplus*step(x - x_epi_end)
    (N_D_nplus scalar or (B,)), as ``etna.core.device._nplus_term``; None = no n+.
    ``buffer`` = (x_epi_end, Lb, N_end) (Lb, N_end scalar or (B,)) adds the graded
    n-/n+ transition of ``etna.core.device._buffer_term`` (``donor_profile`` is the
    shared numpy twin of the devsim expression); None = no buffer.
    """
    x = np.asarray(x, dtype=float)[None, :]
    if np.any(np.asarray(L_transition) <= 0):
        raise ValueError("L_transition must be positive")
    col = lambda a: np.atleast_1d(np.asarray(a, dtype=float))[:, None]
    NA, NJ, NB, Lt = col(N_A_ion), col(N_D_junction), col(N_D_bulk), col(L_transition)
    p_side = (junction_pos - x >= 0).astype(float)
    if buffer is None:
        donors = donor_profile(x, junction_pos, NJ, NB, Lt)
    else:
        x_end, Lb, N_end = buffer
        if np.any(np.asarray(Lb) <= 0) or np.any(np.asarray(N_end) <= 0):
            raise ValueError("buffer Lb and N_end must be positive")
        donors = donor_profile(
            x, junction_pos, NJ, NB, Lt, (x_end, col(Lb), col(N_end))
        )
    if nplus is not None:
        x_end, N_np = nplus
        donors = donors + col(N_np) * (x - x_end >= 0).astype(float)
    acceptors = NA * p_side
    B = max(donors.shape[0], acceptors.shape[0])
    shape = (B, x.shape[1])
    return (
        np.broadcast_to(donors, shape).copy(),
        np.broadcast_to(acceptors, shape).copy(),
    )
