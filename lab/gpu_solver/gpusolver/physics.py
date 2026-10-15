"""4H-SiC material numbers, 1D mesh and graded doping for the batched solver.

Standalone re-statement of the handful of ETNA formulas the electrostatics
depend on (etna/core/sic_material.py, incomplete_ionization.py, device.py),
kept numerically identical so results can be compared against devsim. CGS
units throughout (cm, cm^-3, V, F/cm).
"""

from __future__ import annotations

import numpy as np

EPS_0 = 8.854e-14  # F/cm
Q = 1.602e-19  # C
K_B_EV = 8.617e-5  # eV/K
K_B_J = 1.3806503e-23  # J/K
EPS_R = 9.7
N_I_300 = 5.0e-9  # cm^-3, calibrated anchor used by ETNA
EPS = EPS_R * EPS_0


def _ni_raw(T):
    """First-principles n_i(T) (uncalibrated), as etna.core.sic_material.compute_ni."""
    h, m0 = 6.626e-34, 9.109e-31
    E_g = 3.2965625 - 6.5e-4 * T**2 / (T + 1300.0)
    kT_J = K_B_EV * T * 1.602e-19
    NC = 2 * (2 * np.pi * 0.77 * m0 * kT_J / h**2) ** 1.5 * 1e-6
    NV = 2 * (2 * np.pi * 1.0 * m0 * kT_J / h**2) ** 1.5 * 1e-6
    return np.sqrt(NC * NV) * np.exp(-E_g / (2 * K_B_EV * T))


def intrinsic_density(T):
    """Calibrated n_i(T): n_i(300 K) == 5e-9 cm^-3 exactly."""
    T = np.asarray(T, dtype=float)
    return N_I_300 * _ni_raw(T) / _ni_raw(300.0)


def thermal_voltage(T):
    return K_B_J * np.asarray(T, dtype=float) / Q


def ionized_acceptor_fraction(N_A, T=300.0):
    """ETNA's hybrid Al-acceptor ionization model (Gibbs below 1e18, empirical above)."""
    N_A = np.asarray(N_A, dtype=float)
    T = np.asarray(T, dtype=float)
    E_A_eff = np.maximum(0.220 - 3.4e-8 * N_A ** (1.0 / 3.0), 0.0)
    f_gibbs = 1.0 / (1.0 + 4.0 * np.exp(E_A_eff / (K_B_EV * T)))
    log_NA = np.log10(np.maximum(N_A, 1e18))
    f_ref = np.clip(0.05 + 0.10 * (log_NA - 18.0), 0.01, 0.95)
    T_factor = np.maximum(1.0 + 0.5 * (T - 300.0) / 300.0, 0.5)
    f_emp = np.minimum(f_ref * T_factor, 0.95)
    w = 1.0 / (1.0 + np.exp(-(np.log10(N_A) - 18.0) / 0.5))
    return np.clip((1.0 - w) * f_gibbs + w * f_emp, 0.0, 1.0)


def make_mesh(epi_thickness_cm=10e-4, substrate_thickness_cm=1e-4):
    """Non-uniform 1D mesh with the same anchors/spacings as etna.core.device.

    Geometric grading between anchor spacings (ratio <= 1.25). Node count is
    close to devsim's (~330 for the default geometry) but positions differ;
    validation tests use devsim's own node positions instead.
    """
    xj = substrate_thickness_cm
    L = substrate_thickness_cm + epi_thickness_cm
    anchors = [(0.0, 1e-5), (xj - 5e-6, 1e-6), (xj, 1e-7)]
    for pos in (xj + 2e-4, xj + 5e-4):
        if pos < L - 1e-6:
            anchors.append((pos, 5e-6))
    anchors.append((L, 1e-5))
    x = [0.0]
    for (xa, ha), (xb, hb) in zip(anchors[:-1], anchors[1:]):
        x.extend(_grade(xa, xb, ha, hb)[1:])
    return np.array(x)


def _grade(xa, xb, ha, hb, r=1.25):
    """Points from xa to xb, spacing growing/shrinking geometrically from ha to hb."""
    seg = xb - xa
    pts = [xa]
    # grow from ha toward the middle, shrink toward hb; simple two-sided march
    left, right = [xa], [xb]
    hl, hr = min(ha, seg), min(hb, seg)
    while (right[-1] - left[-1]) > max(hl, hr) * 1.5:
        if hl <= hr:
            left.append(left[-1] + hl)
            hl = min(hl * r, 1e-5)
        else:
            right.append(right[-1] - hr)
            hr = min(hr * r, 1e-5)
    gap = right[-1] - left[-1]
    n_mid = max(1, int(round(gap / max(hl, hr))))
    mid = np.linspace(left[-1], right[-1], n_mid + 1)[1:-1]
    pts = np.concatenate([left, mid, right[::-1]])
    return pts


def graded_doping(x, junction_pos, N_A_ion, N_D_junction, N_D_bulk, L_transition):
    """Donors/Acceptors (B, N) for ETNA's exponentially graded epi profile.

    Scalar or (B,) inputs broadcast over the batch; ``x`` is shared (N,).
    """
    x = np.asarray(x, dtype=float)[None, :]
    NA = np.atleast_1d(np.asarray(N_A_ion, dtype=float))[:, None]
    NJ = np.atleast_1d(np.asarray(N_D_junction, dtype=float))[:, None]
    NB = np.atleast_1d(np.asarray(N_D_bulk, dtype=float))[:, None]
    Lt = np.atleast_1d(np.asarray(L_transition, dtype=float))[:, None]
    # devsim's step(0) == 1: the junction node carries BOTH donors and acceptors
    n_side = (x >= junction_pos).astype(float)
    p_side = (junction_pos - x >= 0).astype(float)
    donors = (NB + (NJ - NB) * np.exp(-np.maximum(x - junction_pos, 0.0) / Lt)) * n_side
    acceptors = NA * p_side
    B = max(donors.shape[0], acceptors.shape[0])
    return (
        np.broadcast_to(donors, (B, x.shape[1])).copy(),
        np.broadcast_to(acceptors, (B, x.shape[1])).copy(),
    )
