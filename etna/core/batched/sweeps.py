"""Reverse-bias sweeps on a batch of diodes: W(V), C(V), psi(x; V) and CCE(V).

All B devices share one mesh (geometry) and are swept over the same voltages;
doping, N_A and T may differ per device. Bias points are solved exactly where
requested (direct jumps with a secant predictor from the last two converged
biases), so unlike devsim's ``cv_sweep`` there is no 0.5 V lattice.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from etna.core.batched.extract import depletion_width
from etna.core.batched.material import (
    PARAMS,
    Q,
    edge_mobilities,
    graded_doping,
    material_params,
    thermal_voltage,
)
from etna.core.batched.poisson import PoissonBatch
from etna.core.batched.transport import ContinuityBatch
from etna.core.device import (
    DEFAULT_N_A,
    acceptor_total_scale,
    device_ionization_fields,
    ionized_acceptor_concentration_for,
)
from etna.core.devsim_lock import serialized
from etna.core.generation_profiles import alpha_generation_profile
from etna.core.incomplete_ionization import (
    ionized_acceptor_concentration,
    resolve_ionization,
)
from etna.core.sic_material import intrinsic_concentration, srh_lifetime


@lru_cache(maxsize=4096)
def _n_i(T):
    return intrinsic_concentration(T, PARAMS)[0]


@lru_cache(maxsize=4096)
def _n_a_ionized(N_A, T):
    return ionized_acceptor_concentration(N_A, T)


@lru_cache(maxsize=4096)
def _n_a_ionized_model(N_A, T, model):
    return ionized_acceptor_concentration_for(N_A, T, model)


def _per_device(fn, *args):
    """Apply scalar ``fn`` over scalar-or-(B,) arguments -> (B,) float array."""
    arrays = np.broadcast_arrays(
        *(np.atleast_1d(np.asarray(a, dtype=float)) for a in args)
    )
    return np.array([fn(*vals) for vals in zip(*arrays)], dtype=float)


def _check_voltages(voltages):
    voltages = np.asarray(voltages, dtype=float)
    if not np.all(np.isfinite(voltages)):
        raise ValueError("voltages must be finite")
    return voltages


@serialized
def make_batch(
    x,
    junction_pos,
    N_A,
    N_D_junction,
    N_D_bulk,
    L_transition,
    T=300.0,
    *,
    nplus,
    buffer,
    donor_table=None,
    incomplete_ionization=None,
    mobility=None,
):
    """PoissonBatch for graded p+/n- diodes on mesh ``x`` (uniform: N_D_junction == N_D_bulk).

    ``mobility`` (the create_sic_device keyword; one value for the batch, None = the
    default "c_axis"): the low-field mobility set, stored as ``pb.mobility`` and used by
    ``srh_continuity`` / ``material.edge_mobilities``.

    Parameters are scalars or (B,) arrays; lengths in cm, densities in cm^-3, T in K.
    ``nplus`` = (x_epi_end, N_D_nplus) for an n+ layer under the epi (the mesh must
    come from ``node_positions(..., nplus_thickness_cm)``); None = no n+.
    ``buffer`` = (x_epi_end, Lb, N_end) for the graded n-/n+ transition over the last
    Lb of the epi (``etna.core.device.create_sic_device`` buffer keywords); None = none.
    Both are REQUIRED keywords (since 2026-10-18; they defaulted to None, so a caller
    that forgot them silently got the legacy no-n+ / abrupt-end geometry while the
    devsim builders default to the stack's n+ and transition).
    ``donor_table`` (optional, one ``etna.core.device.donor_table_spec`` for the whole
    batch): the tabulated epi profile instead of the exponential one (the
    ``create_sic_device`` donor_table keyword); then buffer must be None and the
    N_D_junction / N_D_bulk / L_transition values are not used.
    ``incomplete_ionization`` (opt-in, the create_sic_device keyword; one value for the
    batch): None / False = legacy static doping; True or an IonizationModel = the
    carrier-dependent ionization (audit-v7 A-02), with the same per-node fields and
    static NetDoping as the devsim builder (etna.core.device.device_ionization_fields).
    """
    model = resolve_ionization(incomplete_ionization)
    params = material_params(mobility)  # validates the name before any work
    T = np.atleast_1d(np.asarray(T, dtype=float))
    if model is None:
        N_A_ion = _per_device(_n_a_ionized, N_A, T)
    else:
        N_A_ion = _per_device(lambda a, t: _n_a_ionized_model(a, t, model), N_A, T)
    donors, acceptors = graded_doping(
        x,
        junction_pos,
        N_A_ion,
        N_D_junction,
        N_D_bulk,
        L_transition,
        nplus=nplus,
        buffer=buffer,
        donor_table=donor_table,
    )
    B = donors.shape[0]
    T = np.broadcast_to(T, (B,))
    n_i = _per_device(_n_i, T)
    ionization = None
    if model is not None:
        N_A_b = np.broadcast_to(np.asarray(N_A, dtype=float), (B,))
        N_A_ion_b = np.broadcast_to(N_A_ion, (B,))
        acc = np.broadcast_to(acceptors, donors.shape)
        fields, netdop = [], np.empty(donors.shape)
        for b in range(B):
            f, netdop[b] = device_ionization_fields(
                donors[b], acc[b], N_A_b[b], N_A_ion_b[b], T[b], n_i[b], model
            )
            fields.append(f)
        ionization = (fields, netdop)
    pb = PoissonBatch(
        x, donors, acceptors, n_i, thermal_voltage(T), ionization=ionization
    )
    if model is not None:
        # mobility total doping: Acceptors carry the model's ionized p+ value (A-02)
        pb.acceptor_scales = np.array(
            [acceptor_total_scale(a, i) for a, i in zip(N_A_b, N_A_ion_b)]
        )
    pb.mobility = params.mobility
    return pb


# Largest Poisson bias step of the sweeps (V): longer jumps are walked in sub-steps.
MAX_STEP_V = 10.0


@dataclass
class _Ramp:
    """Walks a batch along increasing cathode bias, remembering the secant history."""

    pb: PoissonBatch
    tol: float = 1e-8
    # Largest cathode step of one Newton solve (V); a longer jump is walked in equal
    # sub-steps (None = jump directly). A cold 0 -> -50 V jump does not converge on
    # low-doped devices (the microdosimeter preset, uniform 3.06e14).
    max_step: "float | None" = None

    def __post_init__(self):
        res = self.pb.solve(self.pb.neutral_guess(), tol=self.tol)
        self.psi, self.Vc, self._psi_prev, self._Vc_prev = res.psi, 0.0, None, None
        self.converged = res.converged.copy()

    def to(self, V_cathode):
        """Solve all devices at ``V_cathode`` (V), via sub-steps of at most
        ``max_step``; returns the PoissonResult of the last solve."""
        if self.max_step is not None:
            n = int(np.ceil(abs(V_cathode - self.Vc) / self.max_step - 1e-9))
            for V in np.linspace(self.Vc, V_cathode, max(n, 1) + 1)[1:-1]:
                self._to(V)
        return self._to(V_cathode)

    def _to(self, V_cathode):
        guess = self.psi
        if self._psi_prev is not None and abs(self.Vc - self._Vc_prev) > 1e-12:
            guess = self.psi + (self.psi - self._psi_prev) * (
                (V_cathode - self.Vc) / (self.Vc - self._Vc_prev)
            )
        res = self.pb.solve(guess, V_cathode=V_cathode, tol=self.tol)
        self._psi_prev, self._Vc_prev = self.psi, self.Vc
        self.psi, self.Vc = res.psi, V_cathode
        self.converged &= res.converged
        return res


@dataclass
class SweepResult:
    voltages: np.ndarray  # (K,) conventional sign (reverse bias negative)
    W: "np.ndarray | None"  # (K, B) cm, depletion EDGE (n/N_D = 0.5); None if not asked
    C: "np.ndarray | None"  # (K, B) F/cm^2 (area = 1), small-signal dQ/dV; None if not asked
    psi: "np.ndarray | None"  # (K, B, N) V; None unless keep_psi
    converged: np.ndarray  # (K, B) bool, cumulative along the sweep


@serialized
def reverse_sweep(
    pb,
    voltages,
    junction_pos,
    epi_thickness_cm,
    tol=1e-8,
    keep_psi=False,
    capacitance=True,
    max_step_V=MAX_STEP_V,
    width=True,
):
    """Solve all devices at each of ``voltages`` (<= 0, in the order given).

    ``C`` is the small-signal capacitance dQ/dV (``PoissonBatch.small_signal_capacitance``),
    the observable of a measured C-V; ``W`` is the n/N_D = 0.5 depletion edge, which is
    NOT eps/C (they differ near full depletion: 9.500 vs 9.764 um at -10 V on the
    pre-2026-10-17 default stack). ``capacitance=False`` skips C (``C`` is None), for
    callers that only need W (reach-through refinement); ``width=False`` skips W (``W``
    is None), for callers that only need C (the C-V inversion). Jumps longer than
    ``max_step_V`` (default 10 V) are walked in Poisson-only sub-steps.

    Jumps of ~100 V (less at low T: the Newton clamp is in V_t units) between
    consecutive points may not converge in one solve; check ``converged``
    before trusting a result. ``keep_psi`` stores the (K, B, N) potentials (large
    for big batches).
    """
    voltages = _check_voltages(voltages)
    ramp = _Ramp(pb, tol, max_step_V)
    W, C, psis, conv = [], [], [], []
    for V in voltages:
        res = ramp.to(-V) if abs(V) > 1e-12 or ramp.Vc != 0.0 else None
        if res is None:
            res = pb.solve(ramp.psi, tol=tol)
        if width:
            W.append(
                depletion_width(pb.x, res.n, pb.donors, junction_pos, epi_thickness_cm)
            )
        if capacitance:
            C.append(pb.small_signal_capacitance(res))
        if keep_psi:
            psis.append(ramp.psi.copy())
        conv.append(ramp.converged.copy())
    return SweepResult(
        voltages,
        np.array(W) if width else None,
        np.array(C) if capacitance else None,
        np.array(psis) if keep_psi else None,
        np.array(conv),
    )


@dataclass
class CCEResult:
    voltages: np.ndarray  # (K,)
    cce: np.ndarray  # (K, B)
    J_collected: np.ndarray  # (K, B) A/cm^2, |Jn + Jp| at the cathode edge
    J_generated: float  # A/cm^2
    converged: np.ndarray  # (K, B) bool


def alpha_generation(x, junction_pos, alpha_range_cm, generation_rate, epi_end=None):
    """Generation profile G(x) (cm^-3 s^-1) exactly as ``cce_vs_bias`` builds it.

    ``epi_end`` (n+ layer present, generation confined to the epi as
    ``confine_generation_to_epi`` with nplus_generation=False): G = 0 beyond it.
    """
    x_epi = np.asarray(x, dtype=float) - junction_pos
    profile = alpha_generation_profile(x_epi, alpha_range_cm=alpha_range_cm)
    peak = np.max(profile)
    G = profile * (generation_rate / peak) if peak > 0 else np.zeros_like(x_epi)
    G[x_epi < 0] = 0.0
    if epi_end is not None:
        G[np.asarray(x, dtype=float) > epi_end + 1e-12] = 0.0
    return G


def srh_continuity(pb, T, N_A):
    """ContinuityBatch of the devsim DD builders (create_dd_device): local Caughey-Thomas
    edge mobility at the batch doping (``material.edge_mobilities``), midgap SRH with
    ``srh_lifetime(T)``. Shared by ``cce_sweep`` and ``dark_carriers``."""
    T = np.broadcast_to(np.asarray(T, dtype=float), (pb.B,))
    N_A = np.broadcast_to(np.asarray(N_A, dtype=float), (pb.B,))
    mu_n, mu_p = edge_mobilities(
        pb.donors,
        pb.acceptors,
        N_A,
        T,
        getattr(pb, "acceptor_scales", None),
        mobility=getattr(pb, "mobility", None),
    )
    return ContinuityBatch(
        pb,
        mu_n,
        mu_p,
        _per_device(lambda t: srh_lifetime(t, "electron", PARAMS), T),
        _per_device(lambda t: srh_lifetime(t, "hole", PARAMS), T),
    )


# Gummel stop of dark_carriers (relative change of the cathode current, and of psi / n /
# p floored at 2 * poisson_tol): 2-3 iterations; n, p agree with devsim's DD Electrons /
# Holes to <= ~1.5e-10 relative at every node (tests/test_field_engine_1d.py).
DARK_CARRIER_GUMMEL_TOL = 1e-10


def dark_carriers(pb, psi, V_cathode, T, N_A, tol=1e-8):
    """Dark drift-diffusion carriers (n, p) at a converged reverse-bias Poisson state.

    ``psi`` (B, N) is a flat-quasi-Fermi solution at ``V_cathode`` (``reverse_sweep``);
    a dark Gummel solve (``srh_continuity``) gives the minority tails that the flat
    quasi-Fermi levels miss (e.g. electrons at the p+ anode, n_i^2 / N_A there instead of
    ~0). Returns the TransportResult; its ``psi`` moves by ~1e-14 V and is not meant to
    replace the Poisson one.
    """
    cb = srh_continuity(pb, T, N_A)
    return cb.solve(
        psi,
        V_cathode=V_cathode,
        tol=DARK_CARRIER_GUMMEL_TOL,
        max_iter=200,
        poisson_tol=tol,
    )


@serialized
def cce_sweep(
    pb,
    voltages,
    junction_pos,
    G,
    T=300.0,
    N_A=DEFAULT_N_A,
    tol=1e-8,
    gummel_tol=1e-6,
    max_step_V=MAX_STEP_V,
):
    """CCE vs reverse bias, as ``etna.core.charge_collection.cce_vs_bias``.

    Per bias: dark Poisson, then Gummel with generation ``G`` (N,) and
    CCE = |J_cathode| / (q * trapz(G)). Like the devsim path, the mobility is LOCAL:
    Caughey-Thomas per node at the batch's own doping (donors + total acceptors, the
    total from ``N_A`` per device), arithmetic edge average (material.edge_mobilities
    = etna.core.device.add_local_mobility_models); lifetimes from ``srh_lifetime(T)``.
    Before 2026-10-18 both mobilities were one scalar per device at the epi doping
    ``N_D`` (and before 2026-10-17 holes at ``N_A``, .issues/audit-v7 A-03); the ignored
    ``N_D`` keyword was removed on 2026-10-18 (the doping is in ``pb``).
    The dark Poisson state is walked from 0 V in steps <= ``max_step_V`` (as
    run_cce_tilted_alpha): a cold jump to e.g. -50 V does not converge on low-doped
    devices; the Gummel solve runs only at the requested biases.
    """
    voltages = _check_voltages(voltages)
    cb = srh_continuity(pb, T, N_A)
    J_gen = Q * np.trapezoid(G, pb.x)
    ramp = _Ramp(pb, tol, max_step_V)
    cce, J, conv = [], [], []
    for V in voltages:
        if abs(V) > 1e-12 or ramp.Vc != 0.0:
            ramp.to(-V)
        tr = cb.solve(ramp.psi, V_cathode=-V, G=G, tol=gummel_tol, poisson_tol=tol)
        J.append(np.abs(tr.J_cathode))
        cce.append(np.clip(np.abs(tr.J_cathode) / J_gen, 0.0, 1.0))
        conv.append(ramp.converged & tr.converged)
    return CCEResult(voltages, np.array(cce), np.array(J), float(J_gen), np.array(conv))
