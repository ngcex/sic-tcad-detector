"""OPT-IN charged radiation-induced deep traps in Poisson (beyond-compensation model).

The default damage model (``etna.core.radiation_damage``) removes donors,
N_D(x) -> max(N_D - eta_removal * Phi_neq, 0), and refuses the points at or beyond the
compensation fluence Phi_comp (``charge_collection.full_compensation_reason``). Real
irradiated 4H-SiC does not stop working there: it turns semi-insulating, with the Fermi
level pinned by the deep levels, and keeps detecting up to ~1e15-1e16 neq (Gsponer et al.,
arXiv:2310.02047; De Napoli, Front. Phys. 10 (2022) 898833; review in
.planning/CARRIER_REMOVAL_REVIEW.md). This module replaces the donor-removal scalar by
the deep levels themselves as CHARGED species in the Poisson equation, the TCAD
approach of Burin et al. (arXiv:2407.16710; their Table I is the trap table used here).

Formulation (opt-in, ``RadiationTrapModel``; nothing here runs unless asked):

1. Trap densities. N_t,k(x) = eta_k * Phi_neq on the epi nodes (x >= junction, the mask
   of the legacy removal), Phi_neq = kappa(E) * Phi_p. eta_k, E_k, sigma_n,k,
   sigma_p,k and the type are ``RadiationDamageParams`` (Burin Table I: Z1/2 acceptor
   Ec-0.67, EH4 acceptor Ec-1.03, EH6/7 donor Ec-1.60; EH6/7 sigma_n is the audit-C-4
   value 1.4e-14, NOT Burin's 9e-12, see radiation_damage.py). The donors N_D(x) are
   NOT reduced: removal emerges from the electrons the acceptor traps bind (stacking
   both would count the removal twice).

2. Charge. Acceptor-like trap: -N_t f; donor-like: +N_t (1 - f), with the steady-state
   SRH (Shockley-Read-Hall) electron occupancy
       f = (c_n n + e_p) / (c_n (n + n1) + c_p (p + p1)),
       c = sigma v_th, e_n = c_n n1, e_p = c_p p1,
       n1 = N_C exp(-(Ec - Et)/kT),  p1 = n_i^2 / n1,
   so that at equilibrium (n p = n_i^2 with the solver's own n_i) f is exactly the
   Fermi-Dirac occupancy n / (n + n1). (The solver's n_i is the calibrated 5e-9 at
   300 K, below sqrt(N_C N_V) exp(-Eg/2kT) ~ 7.6e-9: p1 is tied to it for detailed
   balance, which makes e_p ~0.42x the N_V-based value; EH6/7 is then ~90 % ionized
   in depletion instead of ~74 %.) v_th = sqrt(3 kT / m*) with m_e 0.77 and m_h
   1.0 m0, the convention of ``compute_K_tau`` (only c_p / c_n enters f). Degeneracy
   factor g = 1 (ASSUMED: Burin Table I gives none). Single-level traps: the Z1/2
   double acceptor (-/0 and =/-, negative-U) is one acceptor level as in Burin Table I.

3. Which carrier state sets f. The DARK state at the applied bias, self-consistent
   with Poisson (Newton on psi with the trap charge and its derivative), the TCAD DC
   choice (Burin et al. simulate DC with Sentaurus). For a CCE point the dark trap
   charge is then FROZEN (fixed space charge) in the illuminated Gummel solve. Reason:
   capture at the injected density of the steady-state alpha proxy (G = 1e18
   cm^-3 s^-1: n ~ 1e8 in the depleted epi, c_n n ~ 25 /s) is as fast as Z1/2
   emission, so letting f follow it would make CCE depend on the proxy amplitude,
   while a real alpha pulse (ns) does not repopulate traps with ms..1e14 s time
   constants: between pulses the occupancy is the dark one.

KNOWN CONSEQUENCES (reported, not tuned; see docs/USAGE.md):
- Emergent removal in the NEUTRAL n- region is (eta_Z12 + eta_EH4) = 7.4 /neq (both
  acceptors filled), vs eta_removal = 5 of the legacy scalar and 5-7.5 measured at
  252.7 MeV (Radmanovac et al., arXiv:2510.11304). Compensation (bulk Fermi level at the
  Z1/2 level) at ~N_D / 7.4 neq.
- In the DEPLETED region the steady-state occupancy empties both acceptors (Z1/2
  emission time ~40 ms, EH4 ~30 min at 300 K) and leaves EH6/7 ~90 % ionized: the
  depleted space charge is N_D + ~1.43 Phi_neq, POSITIVE and growing with fluence. The
  model therefore never type-inverts (consistent with SiC; the neutral-bulk Fermi
  level stays above midgap, ~Ec-1.2 eV at 1e16 neq), but it does NOT reproduce the
  measured C-V removal (falling 1/C^2 slope) in the depleted region: that needs
  defects charged in depletion that this three-level table lacks, or a history-frozen
  occupancy (EH6/7 emission ~1e14 s never reaches steady state; frozen at the
  neutral-bulk value the compensated epi is charge-neutral, a dielectric with a
  uniform field, see ``trap_diagnostics``). The depleted N_eff, and with it the CCE
  beyond compensation at a given bias, is therefore an ASSUMPTION of the occupancy
  choice, stated, not a prediction.
- Lifetimes: unchanged (``compute_K_tau``, midgap SRH), so CCE trapping is the legacy
  one; only the electrostatics changes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

_K_B_J = 1.3806e-23  # J/K, as radiation_damage.compute_K_tau
_M0 = 9.109e-31  # kg
_Q = 1.602e-19  # C, etna.core.device.Q (the solvers' constant)


def _solver_vt(T):
    """kT/q with the solvers' constants (etna.core.device K_B_J, Q)."""
    from etna.core.device import K_B_J, Q

    return K_B_J * float(T) / Q


@dataclass(frozen=True)
class DeepTrap:
    """One radiation-induced deep level: introduction rate ``eta`` (cm^-1 per neq),
    level ``E_below_Ec`` (eV), capture cross-sections (cm^2), ``kind`` "acceptor"
    (charge -N_t f) or "donor" (charge +N_t (1 - f)).

    Optional SECOND level (a two-electron defect, e.g. the negative-U carbon vacancy
    Z1/2 = V_C(=/-/0)): ``E2_below_Ec``, ``sigma_n2``, ``sigma_p2`` of the transition
    that adds the SECOND electron (- -> = for an acceptor, + -> 0 for a donor); the
    first level adds the first electron (0 -> -, ++ -> +). Negative U means E2 is
    DEEPER than E1. Charge states hold 0, 1, 2 electrons: "acceptor" 0, -1, -2;
    "donor" +2, +1, 0. The occupancy is carried per LEVEL row (``level_table``) as
    g1 = P(>= 1 electron), g2 = P(2 electrons), so the charge is N (z - g1 - g2), z = 0
    (acceptor) or 2 (donor), the single-level formula N (don - f) summed over rows.
    Degeneracy factors are 1 (assumed, as for single levels)."""

    name: str
    kind: str
    eta: float
    E_below_Ec: float
    sigma_n: float
    sigma_p: float
    E2_below_Ec: Optional[float] = None
    sigma_n2: Optional[float] = None
    sigma_p2: Optional[float] = None

    def __post_init__(self):
        if self.kind not in ("acceptor", "donor"):
            raise ValueError(
                f"DeepTrap.kind must be 'acceptor' or 'donor' (got {self.kind!r})"
            )
        for name in ("eta", "sigma_n", "sigma_p"):
            v = getattr(self, name)
            if not (np.isfinite(v) and v > 0):
                raise ValueError(f"DeepTrap.{name} must be > 0 and finite (got {v!r})")
        if not (0.0 < self.E_below_Ec < 3.26):
            raise ValueError(
                f"DeepTrap.E_below_Ec must be inside the gap (got {self.E_below_Ec!r})"
            )
        second = (self.E2_below_Ec, self.sigma_n2, self.sigma_p2)
        if any(v is not None for v in second):
            if any(v is None for v in second):
                raise ValueError(
                    "DeepTrap second level needs E2_below_Ec, sigma_n2 and sigma_p2"
                )
            for name in ("sigma_n2", "sigma_p2"):
                v = getattr(self, name)
                if not (np.isfinite(v) and v > 0):
                    raise ValueError(f"DeepTrap.{name} must be > 0 (got {v!r})")
            if not (0.0 < self.E2_below_Ec < 3.26):
                raise ValueError("DeepTrap.E2_below_Ec must be inside the gap")

    @property
    def n_levels(self) -> int:
        return 1 if self.E2_below_Ec is None else 2


def burin_table_traps(params=None) -> Tuple[DeepTrap, ...]:
    """The three Burin Table I levels from ``RadiationDamageParams`` (single source)."""
    from etna.core.radiation_damage import RadiationDamageParams

    r = RadiationDamageParams() if params is None else params
    return (
        DeepTrap("Z1/2", r.type_Z12, r.eta_Z12, r.E_Z12, r.sigma_n_Z12, r.sigma_p_Z12),
        DeepTrap("EH4", r.type_EH4, r.eta_EH4, r.E_EH4, r.sigma_n_EH4, r.sigma_p_EH4),
        DeepTrap(
            "EH6/7", r.type_EH67, r.eta_EH67, r.E_EH67, r.sigma_n_EH67, r.sigma_p_EH67
        ),
    )


@dataclass(frozen=True)
class RadiationTrapModel:
    """Parameters of the opt-in charged-trap damage model (hashable, strict).

    ``traps``: the deep levels (default: Burin Table I via ``RadiationDamageParams``).
    The occupancy is the steady-state SRH one of the dark state (module docstring).

    ``donor_removal`` (cm^-1 per neq, default 0): DIRECT donor deactivation, N_D(x) ->
    max(N_D(x) - donor_removal * Phi_neq, 0) on the trap nodes: an occupancy-independent
    loss of shallow donors (e.g. nitrogen bound into a complex), on top of the traps.
    ``occupancy`` (default "steady", the documented behaviour): the trap occupancy that
    sets the dark charge of a CCE point. "kinetic": the ``etna.core.trap_kinetics``
    history instead: the ``initial`` state at 0 V, a ramp to the bias at
    ``ramp_V_per_s``, a hold of ``hold_s`` s (a measurement after biasing), then the
    trap charge frozen in the illuminated solve. Two-level traps need "kinetic" in the
    CCE path (the steady Poisson kernel is single-level; ``trap_kinetics`` handles
    both, its steady limit included). ``label``: provenance of a preset
    (``etna.core.trap_presets``).
    """

    traps: Tuple[DeepTrap, ...] = burin_table_traps()
    donor_removal: float = 0.0
    occupancy: str = "steady"
    ramp_V_per_s: float = 1.0
    hold_s: float = 60.0
    initial: str = "equilibrium"
    label: str = ""

    def __post_init__(self):
        traps = tuple(self.traps)
        if not traps:
            raise ValueError("RadiationTrapModel.traps must not be empty")
        if not all(isinstance(t, DeepTrap) for t in traps):
            raise TypeError("RadiationTrapModel.traps must be DeepTrap entries")
        object.__setattr__(self, "traps", traps)
        if not (np.isfinite(self.donor_removal) and self.donor_removal >= 0):
            raise ValueError("RadiationTrapModel.donor_removal must be >= 0")
        if self.occupancy not in ("steady", "kinetic"):
            raise ValueError(
                "RadiationTrapModel.occupancy must be 'steady' or 'kinetic'"
            )
        if not (self.ramp_V_per_s > 0 and self.hold_s >= 0):
            raise ValueError("ramp_V_per_s must be > 0 and hold_s >= 0")
        if self.initial not in ("equilibrium", "in_beam", "neutral"):
            raise ValueError("initial must be 'equilibrium', 'in_beam' or 'neutral'")

    @property
    def neutral_removal_rate(self) -> float:
        """Removal rate (cm^-1 per neq) of a neutral n-type region with E_F above every
        level: every acceptor level holds its electrons (two for a two-level acceptor),
        donors are neutral; plus the direct ``donor_removal``."""
        return float(
            self.donor_removal
            + sum(t.eta * t.n_levels for t in self.traps if t.kind == "acceptor")
        )

    @property
    def multilevel(self) -> bool:
        return any(t.n_levels > 1 for t in self.traps)


def level_table(model):
    """Per-LEVEL rows of ``model`` (a two-level trap gives two rows, first level
    first): dicts with name, trap (index), eta, E (eV below Ec), sigma_n, sigma_p, don
    (1 donor / 0 acceptor) and pair (row index of the second level for a first level,
    -2 for a second level, -1 for a single level)."""
    rows = []
    for k, t in enumerate(model.traps):
        don = 1.0 if t.kind == "donor" else 0.0
        base = dict(trap=k, eta=t.eta, don=don)
        if t.n_levels == 1:
            rows.append(dict(base, name=t.name, E=t.E_below_Ec, sigma_n=t.sigma_n,
                             sigma_p=t.sigma_p, pair=-1))  # fmt: skip
        else:
            i = len(rows)
            rows.append(dict(base, name=t.name, E=t.E_below_Ec, sigma_n=t.sigma_n,
                             sigma_p=t.sigma_p, pair=i + 1))  # fmt: skip
            rows.append(dict(base, name=t.name + "(2)", E=t.E2_below_Ec,
                             sigma_n=t.sigma_n2, sigma_p=t.sigma_p2, pair=-2))  # fmt: skip
    return rows


DEFAULT_TRAP_MODEL = RadiationTrapModel()


def resolve_traps(spec):
    """``radiation_traps`` keyword -> RadiationTrapModel or None (legacy removal).

    False / None -> None; True -> DEFAULT_TRAP_MODEL; a RadiationTrapModel as is.
    Anything else (also 1 / "yes") is a TypeError: strict contract."""
    if spec is None or spec is False:
        return None
    if spec is True:
        return DEFAULT_TRAP_MODEL
    if isinstance(spec, RadiationTrapModel):
        return spec
    raise TypeError(
        "radiation_traps must be a bool, None or a RadiationTrapModel "
        f"(got {type(spec).__name__})"
    )


@dataclass(frozen=True)
class TrapFields:
    """Per-node trap inputs, arrays (L, N) over the LEVEL rows of ``level_table`` (L =
    K for single-level traps): density N_t (cm^-3), n1, p1 (cm^-3), r = c_p / c_n, and
    don = 1 for donor-like / 0 for acceptor-like; ``pair`` (L,) int64 as
    ``level_table`` (None: every level single)."""

    N_t: np.ndarray
    n1: np.ndarray
    p1: np.ndarray
    r: np.ndarray
    don: np.ndarray
    pair: Optional[np.ndarray] = None

    @property
    def multilevel(self) -> bool:
        return self.pair is not None and bool(np.any(self.pair >= 0))


def _vth(T, m_rel):
    return np.sqrt(3 * _K_B_J * T / (m_rel * _M0)) * 100.0


def trap_fields(model, fluence_neq, mask, T, n_i) -> TrapFields:
    """TrapFields on the nodes of ``mask`` (bool (N,), traps there, zero elsewhere),
    one row per level (``level_table``)."""
    from etna.core.sic_material import SiC4H_Parameters, effective_dos

    p = SiC4H_Parameters()
    NC, _ = effective_dos(float(T))
    V_t = _solver_vt(T)
    mask = np.asarray(mask, dtype=bool)
    N = mask.size
    rows = level_table(model)
    L = len(rows)
    ratio = _vth(T, p.m_h_dos) / _vth(T, p.m_e_dos)
    out = {k: np.zeros((L, N)) for k in ("N_t", "n1", "p1", "r", "don")}
    for k, t in enumerate(rows):
        n1 = NC * np.exp(-t["E"] / V_t)
        out["N_t"][k, mask] = t["eta"] * float(fluence_neq)
        out["n1"][k] = n1
        out["p1"][k] = float(n_i) ** 2 / n1
        out["r"][k] = t["sigma_p"] / t["sigma_n"] * ratio
        out["don"][k] = t["don"]
    pair = np.array([t["pair"] for t in rows], dtype=np.int64)
    return TrapFields(**out, pair=pair if np.any(pair != -1) else None)


def occupancy(fields, n, p):
    """Steady-state SRH electron occupancy (L, N) of every level row at carrier
    densities n, p (N,). Single level: f = (n + r p1) / (n + n1 + r (p + p1)).
    Two-level trap (rows a, b): the chain 0 <-> 1 <-> 2 electrons has zero net flux on
    each link, x = R+_a / R-_a, y = R+_b / R-_b with R+ = n + r p1, R- = n1 + r p (in
    units of c_n), and the rows carry g1 = (x + x y) / (1 + x + x y) and
    g2 = x y / (1 + x + x y)."""
    n = np.asarray(n, dtype=float)[None, :]
    p = np.asarray(p, dtype=float)[None, :]
    r = fields.r
    f = (n + r * fields.p1) / (n + fields.n1 + r * (p + fields.p1))
    if fields.multilevel:
        for a in np.nonzero(fields.pair >= 0)[0]:
            b = int(fields.pair[a])
            Rpa = n[0] + r[a] * fields.p1[a]
            Rma = fields.n1[a] + r[a] * p[0]
            Rpb = n[0] + r[b] * fields.p1[b]
            Rmb = fields.n1[b] + r[b] * p[0]
            den = Rpa * Rpb + Rpa * Rmb + Rma * Rmb  # (1 + x + xy) Rma Rmb
            f[a] = Rpa * (Rpb + Rmb) / den
            f[b] = Rpa * Rpb / den
    return f


def trap_charge(fields, n, p):
    """Net trap charge density (N,), cm^-3 (positive = positive charge)."""
    f = occupancy(fields, n, p)
    return (fields.N_t * (fields.don - f)).sum(axis=0)


def neutral_state(static_net, fields, n_i, iterations=200):
    """Local neutrality p - n + static_net + rho_trap(n, p) = 0 at every node, n p = n_i^2.

    Bisection in u = ln(n / n_i) (the residual decreases strictly in u: dn/du > 0,
    dp/du < 0, df/du >= 0), as ``incomplete_ionization.neutral_state``. Returns
    (n, p, rho_trap), each (N,)."""
    static_net = np.asarray(static_net, dtype=float)
    n_i = float(n_i)
    N_max = max(
        float(np.max(np.abs(static_net))), float(np.max(fields.N_t.sum(axis=0))), 1.0
    )
    half = np.log(N_max / n_i) + 5.0
    lo = np.full(static_net.shape, -half)
    hi = np.full(static_net.shape, half)

    def residual(u):
        n = n_i * np.exp(u)
        pp = n_i * np.exp(-u)
        return pp - n + static_net + trap_charge(fields, n, pp)

    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        pos = residual(mid) > 0.0
        lo = np.where(pos, mid, lo)
        hi = np.where(pos, hi, mid)
    u = 0.5 * (lo + hi)
    n = n_i * np.exp(u)
    pp = n_i * np.exp(-u)
    return n, pp, trap_charge(fields, n, pp)


def trap_diagnostics(N_D, fluence_neq, T=300.0, model=None):
    """Bulk (uniform N_D) state of the trap model at ``fluence_neq``: no solver, pure
    local physics, for validation and docs.

    Returns a dict: n_neutral (cm^-3), Ec_minus_EF (eV) of the neutral bulk,
    resistivity proxy (Ohm cm, electrons only, mu_n = 950: an upper-mobility estimate),
    N_eff_depleted_ss = net charge with n = p = 0 under steady-state occupancy (the
    space charge this model puts in a depleted region), N_eff_depleted_frozen = the same
    with the occupancy frozen at the neutral-bulk value (the history-frozen
    alternative, not used by the solver; ~0 by construction: a compensated epi would
    then be a dielectric with a uniform field), ``neutral_removal_rate`` and the per-trap
    emission times 1 / (e_n + e_p) (s)."""
    from etna.core.sic_material import intrinsic_concentration, effective_dos

    model = DEFAULT_TRAP_MODEL if model is None else model
    n_i = intrinsic_concentration(float(T))[0]
    NC, _ = effective_dos(float(T))
    f = trap_fields(model, fluence_neq, np.array([True]), T, n_i)
    N_D_eff = max(float(N_D) - model.donor_removal * float(fluence_neq), 0.0)
    n, p, rho = neutral_state(np.array([N_D_eff]), f, n_i)
    V_t = _solver_vt(T)
    f_neu = occupancy(f, n, p)
    f_dep = occupancy(f, np.zeros(1), np.zeros(1))
    Ndep_ss = float(N_D + (f.N_t * (f.don - f_dep)).sum())
    Ndep_frozen = float(N_D + (f.N_t * (f.don - f_neu)).sum())
    from etna.core.sic_material import SiC4H_Parameters

    prm = SiC4H_Parameters()
    vn = _vth(T, prm.m_e_dos)
    tau_e = {}
    for k, t in enumerate(level_table(model)):
        e_n = t["sigma_n"] * vn * f.n1[k, 0]
        e_p = f.r[k, 0] * t["sigma_n"] * vn * f.p1[k, 0]  # c_p p1, c_p = r c_n
        tau_e[t["name"]] = 1.0 / (e_n + e_p)
    rm = float(N_D) - N_D_eff  # direct removal: both cases
    Ndep_ss -= rm
    Ndep_frozen -= rm
    return {
        "n_neutral": float(n[0]),
        "Ec_minus_EF": float(-V_t * np.log(n[0] / NC)),
        "resistivity_ohm_cm": float(1.0 / (_Q * 950.0 * n[0])),
        "N_eff_depleted_ss": Ndep_ss,
        "N_eff_depleted_frozen": Ndep_frozen,
        "neutral_removal_rate": model.neutral_removal_rate,
        "emission_time_s": tau_e,
    }
