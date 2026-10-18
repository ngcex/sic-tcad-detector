"""Batched dark current (the "TAT" effective-generation model) and the steady-state
FLASH CCE (SRH + Auger) on the numba 1D engine.

Same equations as the devsim path, term by term (no reformulation):

* dark current (etna.core.dark_current.setup_tat_model): SRH with the Z1/2 trap level,
  U_SRH = (n p - n_i^2) / (tau_p (n + n1_tat Gamma) + tau_n (p + p1_tat Gamma)), and the
  effective generation G = G0 Gamma min(E_node / E_ref, 1) (times step(epi_end - x) with
  an n+ layer), G0 = N_t n_i(T)/n_i(300). E_node is the 1D node field of
  ``_compute_node_efield`` (end nodes: their edge; interior: mean of the two adjacent
  |E|; floored at 1e3 V/cm). The devsim path iterates the frozen field to its fixed point
  E_node = E(state) (``_resolve_self_consistent``); here G is re-evaluated from the
  present potential at every Gummel iteration, which converges to that same fixed point.
  Extraction as ``extract_dark_current_components``: J_total in the balanced form,
  J_SRH / J_TAT as trapezoid integrals over all nodes, J_SRV at the cathode node.
  The Hurkx Gamma is 1 below ~4.5 MV/cm (above the breakdown field); a device whose
  field would make Gamma != 1 is reported as not converged (the caller falls back to
  devsim), so the engine never runs a regime it was not cross-checked in.
* FLASH (etna.core.flash_recombination.cce_vs_dose_rate): midgap SRH plus Auger
  (C_n n + C_p p)(n p - n_i^2) with a proton generation profile; CCE = |J| / (q trapz G).

Parity with devsim: tests/test_dark_current_engine.py, tests/test_flash_engine.py.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from etna.core.batched.material import Q, edge_mobilities, thermal_voltage
from etna.core.batched.poisson import PoissonBatch
from etna.core.batched.sweeps import MAX_STEP_V, _Ramp
from etna.core.batched.transport import ContinuityBatch
from etna.core.devsim_lock import serialized

# Gummel stop for the dark current: relative change of the cathode current between two
# iterations. The dark current is a generation integral (well conditioned), so this is
# reached in a few iterations; 1e-11 puts the engine-vs-devsim difference at the
# devsim DC tolerance (tests/test_dark_current_engine.py).
DARK_GUMMEL_TOL = 1e-11
DARK_GUMMEL_MAX_ITER = 200

_KEYS = ("J_total", "J_SRH", "J_TAT", "J_SRV")


@dataclass
class DarkResult:
    voltages: np.ndarray  # (K,) anode bias, reverse negative
    J_total: np.ndarray  # (K, B) A/cm^2, positive under reverse bias
    J_SRH: np.ndarray  # (K, B)
    J_TAT: np.ndarray  # (K, B)
    J_SRV: np.ndarray  # (K, B)
    converged: np.ndarray  # (K, B) bool, cumulative along the sweep
    gamma_active: np.ndarray  # (K, B) bool: Gamma != 1 somewhere (outside the engine)


def _col(a, B):
    return np.broadcast_to(np.asarray(a, dtype=float), (B,)).copy()


class TATGeneration:
    """G(psi) of the TAT model for a batch on mesh x: ``__call__(psi_T)`` -> (N, B)."""

    def __init__(self, x, G0, Kt_numerator, Kt_denom_const, E_ref, gen_mask):
        from etna.core.dark_current import gamma_from_kt, node_field_1d

        self._node_field = node_field_1d
        self._gamma = gamma_from_kt
        self.h = np.diff(np.asarray(x, dtype=float))[:, None]  # (N-1, 1)
        self.G0 = np.asarray(G0, dtype=float)[None, :]  # (1, B)
        self.Kt_num = np.asarray(Kt_numerator, dtype=float)[None, :]
        self.Kt_den = np.asarray(Kt_denom_const, dtype=float)[None, :]
        self.E_ref = float(E_ref)
        self.mask = np.asarray(gen_mask, dtype=float)[:, None]  # (N, 1)
        self.gamma_active = None

    def field(self, psi_T):
        """Node field (N, B), V/cm, as _compute_node_efield (1D)."""
        abs_E = np.abs((psi_T[:-1] - psi_T[1:]) / self.h)  # (N-1, B)
        return self._node_field(abs_E.T).T

    def __call__(self, psi_T):
        E = self.field(psi_T)
        Kt = self.Kt_num / (self.Kt_den * E)
        Gamma = self._gamma(Kt)
        self.gamma_active = np.any(Gamma != 1.0, axis=0)
        G = self.G0 * Gamma * np.minimum(E / self.E_ref, 1.0)
        return G * self.mask


@serialized
def dark_sweep(
    x,
    donors,
    acceptors,
    T,
    N_A,
    n_i,
    voltages,
    *,
    tau_n,
    tau_p,
    n1,
    p1,
    G0,
    Kt_numerator,
    Kt_denom_const,
    E_ref,
    epi_end=None,
    S_n=0.0,
    S_p=0.0,
    tol=1e-8,
    gummel_tol=DARK_GUMMEL_TOL,
    max_iter=DARK_GUMMEL_MAX_ITER,
    max_step_V=MAX_STEP_V,
):
    """Dark-current components of B devices on mesh ``x`` at each anode bias of
    ``voltages`` (<= 0, in the order given).

    ``donors`` / ``acceptors`` (B, N): the doping the devsim device carries (acceptors
    ionized, as devsim's Acceptors model); ``T``, ``N_A`` (total, for the local mobility),
    ``n_i``, ``tau_n``, ``tau_p``, ``n1``, ``p1``, ``G0``, ``Kt_*``: scalar or (B,)
    (``etna.core.dark_current.tat_constants``). ``epi_end`` (n+ layer present): the
    effective generation is confined to x <= epi_end (step(epi_end - x), step(0) = 1).
    ``S_n``/``S_p`` (cm/s): the post-solve SRV diagnostic at the cathode (0 = off).
    Bias convention: the cathode is raised to -V (equivalent to the devsim path's anode
    at V; n and p depend only on psi - phi).
    """
    x = np.asarray(x, dtype=float)
    donors = np.atleast_2d(np.asarray(donors, dtype=float))
    B, N = donors.shape
    acceptors = np.broadcast_to(np.asarray(acceptors, dtype=float), (B, N))
    T, N_A, n_i = _col(T, B), _col(N_A, B), _col(n_i, B)
    voltages = np.asarray(voltages, dtype=float)
    if not np.all(np.isfinite(voltages)) or np.any(voltages > 0):
        raise ValueError("voltages must be finite reverse biases (<= 0)")
    pb = PoissonBatch(x, donors, acceptors, n_i, thermal_voltage(T))
    mu_n, mu_p = edge_mobilities(donors, acceptors, N_A, T)
    n1, p1 = _col(n1, B), _col(p1, B)
    cb = ContinuityBatch(pb, mu_n, mu_p, _col(tau_n, B), _col(tau_p, B), n1=n1, p1=p1)
    mask = np.ones(N) if epi_end is None else (epi_end - x >= 0).astype(float)
    gen = TATGeneration(
        x, _col(G0, B), _col(Kt_numerator, B), _col(Kt_denom_const, B), E_ref, mask
    )
    S_n, S_p = _col(S_n, B), _col(S_p, B)
    ramp = _Ramp(pb, tol, max_step_V)
    out = {k: [] for k in _KEYS}
    conv, gam = [], []
    ok = ramp.converged.copy()
    ni2 = n_i**2
    tau_n_c, tau_p_c = cb.tau_n[:, None], cb.tau_p[:, None]
    for V in voltages:
        if abs(V) > 1e-12 or ramp.Vc != 0.0:
            ramp.to(-V)
        tr = cb.solve(
            ramp.psi,
            V_cathode=-V,
            tol=gummel_tol,
            max_iter=max_iter,
            poisson_tol=tol,
            G_fn=gen,
        )
        n, p, G = tr.n, tr.p, tr.G
        # extract_dark_current_components: U_SRH_only (Z1/2 level, no Gamma) and the
        # effective generation, trapezoid over all nodes, contact-current sign
        U_srh = (n * p - ni2[:, None]) / (
            tau_p_c * (n + n1[:, None]) + tau_n_c * (p + p1[:, None])
        )
        out["J_SRH"].append(-Q * np.trapezoid(U_srh, x, axis=1))
        out["J_TAT"].append(-Q * np.trapezoid(-G, x, axis=1))
        n_c, p_c = n[:, -1], p[:, -1]
        with np.errstate(divide="ignore", invalid="ignore"):
            srv = -Q * (n_c * p_c - ni2) / ((n_c + n_i) / S_p + (p_c + n_i) / S_n)
        out["J_SRV"].append(np.where((S_n > 0) & (S_p > 0), srv, 0.0))
        # the engine's balanced cathode current is positive INTO the device for electrons
        # leaving at the cathode; extract_steady_state_current's sign is the opposite
        # (positive under reverse bias): measured, not assumed (identical magnitudes)
        out["J_total"].append(-tr.J_cathode)
        g_act = np.asarray(gen.gamma_active, dtype=bool)
        ok = ok & ramp.converged & tr.converged & ~g_act & np.isfinite(tr.J_cathode)
        conv.append(ok.copy())
        gam.append(g_act)
    return DarkResult(
        voltages,
        *(np.array(out[k]) for k in _KEYS),
        converged=np.array(conv),
        gamma_active=np.array(gam),
    )


@serialized
def generation_cce(
    x,
    donors,
    acceptors,
    T,
    N_A,
    n_i,
    V_bias,
    G,
    *,
    tau_n,
    tau_p,
    C_n=0.0,
    C_p=0.0,
    tol=1e-8,
    gummel_tol=1e-11,
    max_iter=DARK_GUMMEL_MAX_ITER,
    max_step_V=MAX_STEP_V,
):
    """CCE = |J_cathode| / (q trapz(G, x)) of B devices at one anode bias ``V_bias``,
    each with its own generation ``G`` (B, N), midgap SRH and (optionally) Auger
    C_n/C_p: the steady-state FLASH sweep of cce_vs_dose_rate. Returns (cce, converged),
    both (B,); cce is clipped to [0, 1] as compute_cce_from_dd."""
    x = np.asarray(x, dtype=float)
    G = np.atleast_2d(np.asarray(G, dtype=float))
    B, N = G.shape
    donors = np.broadcast_to(np.asarray(donors, dtype=float), (B, N)).copy()
    acceptors = np.broadcast_to(np.asarray(acceptors, dtype=float), (B, N))
    T, N_A, n_i = _col(T, B), _col(N_A, B), _col(n_i, B)
    pb = PoissonBatch(x, donors, acceptors, n_i, thermal_voltage(T))
    mu_n, mu_p = edge_mobilities(donors, acceptors, N_A, T)
    cb = ContinuityBatch(
        pb, mu_n, mu_p, _col(tau_n, B), _col(tau_p, B), C_n=C_n, C_p=C_p
    )
    ramp = _Ramp(pb, tol, max_step_V)
    V_c = -float(V_bias)
    if abs(V_c) > 1e-12:
        ramp.to(V_c)
    tr = cb.solve(
        ramp.psi, V_cathode=V_c, G=G, tol=gummel_tol, max_iter=max_iter, poisson_tol=tol
    )
    J_gen = Q * np.trapezoid(G, x, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cce = np.where(J_gen > 0, np.abs(tr.J_cathode) / J_gen, 0.0)
    cce = np.clip(cce, 0.0, 1.0)
    return cce, ramp.converged & tr.converged & np.isfinite(cce)
