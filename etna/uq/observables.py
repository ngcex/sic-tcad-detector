"""Batched evaluators of the UQ observables, on the numba engines only.

Each evaluator takes the base ``DeviceConfig``, the parameter matrix of n samples and
the ``ParameterSpace`` and returns an ``EvalOut``: the output axis ``x`` (bias,
fluence, ...), the outputs ``Y`` (n, m) with NaN where a sample did not converge or
lies outside the engine's scope, and bookkeeping. NO devsim fallback: an unconverged
engine sample is NaN and counted (``n_failed``), never silently re-solved on devsim
(which would cost seconds per sample).

Samples are grouped by mesh geometry (epi thickness) and solved as one batch per group
(``etna.api.simulation._engine_batch``, the facades' own batch builder, so the doping,
mesh and material constants are those of ``run_cv`` / ``run_cce``). A continuous epi
thickness would give one group per sample; it is snapped to a lattice of
``epi_quantum_um`` (default 0.02 um; 0 = off) and the largest snap is reported
(``epi_snap_max_um``): on the default 10 um device it moves C by at most ~0.1 %.

Model inputs that are not ``DeviceConfig`` fields are applied here, per device:

- ``mu_n_scale`` / ``mu_p_scale``: factors on the local (Caughey-Thomas, edge-averaged)
  low-field mobilities of every device (C-V and V_fd do not depend on them: an index of
  ~0 there is the correct answer);
- ``alpha_range_um``: the Am-241 alpha range of the CCE generation profile (per device);
- ``E_pair_eV``: the pair energy of the absolute collected charge (``alpha_charge``);
- ``d_removal`` and ``r_n``: the calibrated direct donor removal and the SiC/Si fission
  NIEL ratio of the damage observables, tied by ``etna.uq.params.tied_removal``.
"""

from __future__ import annotations

import dataclasses
import logging
from dataclasses import dataclass, field

import numpy as np

from etna.uq.params import ParameterSpace, tied_removal

logger = logging.getLogger(__name__)

DEFAULT_EPI_QUANTUM_UM = 0.02


@dataclass
class EvalOut:
    x: np.ndarray  # (m,) output axis
    Y: np.ndarray  # (n, m)
    n_failed: int = 0  # samples with any non-finite output
    info: dict = field(default_factory=dict)


@dataclass(frozen=True)
class Observable:
    name: str
    x_label: str
    y_label: str
    evaluate: object  # fn(config, values (n, d), space, **kw) -> EvalOut
    resolution: float  # absolute output spread below which indices are noise
    description: str
    output_names: tuple = ()  # names of the columns for multi-scalar observables
    log_y: bool = False


# ------------------------------------------------------------------ sample -> configs


def _model_values(space: ParameterSpace, values, name, default):
    """Column of model parameter ``name`` (n,), or ``default`` broadcast."""
    for j, p in enumerate(space.params):
        if p.target == f"model:{name}":
            return np.asarray(values[:, j], dtype=float)
    return np.full(values.shape[0], float(default))


def sample_configs(config, values, space: ParameterSpace, epi_quantum_um=None):
    """One DeviceConfig per sample (``config:`` targets replaced), plus the largest epi
    snap (um)."""
    from etna.uq.params import check_targets

    check_targets(space, config)
    q = DEFAULT_EPI_QUANTUM_UM if epi_quantum_um is None else float(epi_quantum_um)
    cols = [(j, p.target.split(":", 1)[1]) for j, p in enumerate(space.params)
            if p.target.startswith("config:")]  # fmt: skip
    out, snap = [], 0.0
    for row in np.atleast_2d(values):
        kw = {f: float(row[j]) for j, f in cols}
        if "epi_thickness_um" in kw and q > 0:
            e = kw["epi_thickness_um"]
            kw["epi_thickness_um"] = round(round(e / q) * q, 9)
            snap = max(snap, abs(kw["epi_thickness_um"] - e))
        out.append(dataclasses.replace(config, **kw) if kw else config)
    return out, snap


def _groups(configs):
    """Engine groups {mesh key: [indices]} and the out-of-box indices (NaN)."""
    from etna.api.device import device_kwargs
    from etna.api.simulation import validated_domain_reason

    groups, outside = {}, []
    for i, c in enumerate(configs):
        if c.half_width_um is not None:
            raise NotImplementedError("etna.uq: 1D devices only (half_width_um=None)")
        if validated_domain_reason(c) is not None:
            outside.append(i)
            continue
        kw = device_kwargs(c)
        key = (
            kw["epi_thickness_cm"], kw["substrate_thickness_cm"],
            kw.get("nplus_thickness_cm"), kw.get("buffer_thickness_cm") is not None,
            kw.get("donor_table"), kw.get("incomplete_ionization"), kw.get("mobility"),
        )  # fmt: skip
        groups.setdefault(key, []).append(i)
    return groups, outside


def _batch(configs, idxs, refine=0):
    from etna.api.device import device_kwargs
    from etna.api.models import refine_mesh
    from etna.api.simulation import _engine_batch
    from etna.core.batched.mesh import node_positions

    kws = [device_kwargs(configs[i]) for i in idxs]
    x = None
    if refine:
        k = kws[0]
        x = refine_mesh(node_positions(k["epi_thickness_cm"], k["substrate_thickness_cm"],
                                       k.get("nplus_thickness_cm")), refine)  # fmt: skip
    pb, x = _engine_batch(kws, x)
    return pb, x, kws


def _run_groups(configs, m, solve, refine=0):
    """Y (n, m) filled group by group by ``solve(pb, x, kws, idxs) -> (vals (B, m),
    ok (B,))``; NaN outside the box / unconverged."""
    n = len(configs)
    Y = np.full((n, m), np.nan)
    groups, outside = _groups(configs)
    unconv = 0
    for idxs in groups.values():
        pb, x, kws = _batch(configs, idxs, refine)
        vals, ok = solve(pb, x, kws, np.asarray(idxs))
        vals = np.asarray(vals, dtype=float).reshape(len(idxs), m)
        vals[~np.asarray(ok, dtype=bool)] = np.nan
        unconv += int(np.sum(~np.asarray(ok, dtype=bool)))
        Y[idxs] = vals
    return Y, {"n_outside_box": len(outside), "n_unconverged": unconv,
               "n_groups": len(groups)}  # fmt: skip


def _scaled_continuity(pb, kws, s_n, s_p):
    """ContinuityBatch of ``cce_sweep`` (srh_continuity) with per-device mobility
    factors."""
    from etna.core.batched.material import PARAMS, edge_mobilities
    from etna.core.batched.transport import ContinuityBatch
    from etna.core.sic_material import srh_lifetime

    T = np.array([k["T"] for k in kws], dtype=float)
    N_A = np.array([k["N_A"] for k in kws], dtype=float)
    mu_n, mu_p = edge_mobilities(pb.donors, pb.acceptors, N_A, T,
                                 getattr(pb, "acceptor_scales", None),
                                 mobility=getattr(pb, "mobility", None))  # fmt: skip
    tau_n = np.array([srh_lifetime(t, "electron", PARAMS) for t in T])
    tau_p = np.array([srh_lifetime(t, "hole", PARAMS) for t in T])
    return ContinuityBatch(pb, mu_n * np.asarray(s_n)[None, :],
                           mu_p * np.asarray(s_p)[None, :], tau_n, tau_p)  # fmt: skip


# ------------------------------------------------------------------ electrostatics


def fd_threshold_cm(config) -> float:
    """Depletion-edge depth (cm from the junction) at which ``run_cv`` calls the epi
    fully depleted: the rule of ``etna.api.simulation._reach_through_V`` (buffer start
    minus one Debye length on the default stack)."""
    import etna.api.simulation as sim

    epi_cm = config.epi_thickness_um * 1e-4
    if config.donor_table is not None:
        return sim.table_sensitive_end_cm(
            sim.donor_table_spec(epi_cm, config.donor_table)
        ) - sim._debye_length_epi_end(config)
    if config.buffer_thickness_um is not None:
        return (
            epi_cm
            - config.buffer_thickness_um * 1e-4
            - sim._debye_length_epi_end(config)
        )
    if config.nplus_thickness_um is not None:
        return epi_cm - sim._debye_length_epi_end(config)
    x = sim.node_positions(epi_cm, config.substrate_thickness_um * 1e-4)
    return epi_cm - (x[-1] - x[-2])


def _first_crossing(V, W, thr):
    """|V| where W first reaches thr, linear interpolation inside the bracketing grid
    interval (NaN if never)."""
    hit = np.flatnonzero(W >= thr)
    if not hit.size:
        return np.nan
    k = hit[0]
    if k == 0:
        return float(abs(V[0]))
    w0, w1 = W[k - 1], W[k]
    t = (thr - w0) / (w1 - w0) if w1 > w0 else 1.0
    return float(abs(V[k - 1]) + t * (abs(V[k]) - abs(V[k - 1])))


def eval_v_fd(config, values, space, *, v_max=15.0, step=0.25, epi_quantum_um=None,
              refine=0, **_):  # fmt: skip
    """Full-depletion voltage |V| (V): first crossing of the ``run_cv`` threshold by the
    depletion edge on a ``step`` grid to ``v_max``, interpolated linearly (continuous in
    the inputs; ``run_cv``'s reach_through_V is the same rule on a 0.05 V lattice)."""
    from etna.core.batched.sweeps import reverse_sweep

    configs, snap = sample_configs(config, values, space, epi_quantum_um)
    V = -np.arange(0.0, v_max + 0.5 * step, step)

    def solve(pb, x, kws, idxs):
        res = reverse_sweep(pb, V, kws[0]["substrate_thickness_cm"],
                            kws[0]["epi_thickness_cm"], capacitance=False)  # fmt: skip
        out = np.array([
            _first_crossing(V, res.W[:, j], fd_threshold_cm(configs[i]))
            for j, i in enumerate(idxs)
        ])  # fmt: skip
        return out, res.converged.all(axis=0)

    Y, info = _run_groups(configs, 1, solve, refine)
    info["epi_snap_max_um"] = snap
    return EvalOut(np.array([0.0]), Y, int(np.sum(~np.isfinite(Y).all(axis=1))), info)


def eval_cv(
    config, values, space, *, voltages=None, epi_quantum_um=None, refine=0, **_
):
    """Small-signal C (F, for the config area) at each bias of ``voltages``."""
    from etna.core.batched.sweeps import reverse_sweep

    V = np.linspace(0.0, -30.0, 31) if voltages is None else np.asarray(voltages, float)
    configs, snap = sample_configs(config, values, space, epi_quantum_um)

    def solve(pb, x, kws, idxs):
        res = reverse_sweep(pb, V, kws[0]["substrate_thickness_cm"],
                            kws[0]["epi_thickness_cm"], width=False)  # fmt: skip
        area = np.array([configs[i].area_cm2 for i in idxs])
        return (res.C * area[None, :]).T, res.converged.all(axis=0)

    Y, info = _run_groups(configs, V.size, solve, refine)
    info["epi_snap_max_um"] = snap
    return EvalOut(V, Y, int(np.sum(~np.isfinite(Y).all(axis=1))), info)


# ------------------------------------------------------------------ transport


def _alpha_G(x, junction, range_cm, epi_end):
    from etna.core.batched.sweeps import alpha_generation

    return np.array([alpha_generation(x, junction, r, 1e18, epi_end=epi_end)
                     for r in np.atleast_1d(range_cm)])  # fmt: skip


def _cce_solve(V, s_n, s_p, rng_cm, gummel_tol, tol=1e-8):
    from etna.core.batched.material import Q
    from etna.core.batched.sweeps import MAX_STEP_V, _Ramp

    def solve(pb, x, kws, idxs):
        sub = kws[0]["substrate_thickness_cm"]
        epi_end = (
            sub + kws[0]["epi_thickness_cm"]
            if kws[0].get("nplus_thickness_cm")
            else None
        )
        cb = _scaled_continuity(pb, kws, s_n[idxs], s_p[idxs])
        G = _alpha_G(x, sub, rng_cm[idxs], epi_end)  # (B, N)
        J_gen = Q * np.trapezoid(G, x, axis=1)
        ramp = _Ramp(pb, tol, MAX_STEP_V)
        out, conv = [], np.ones(len(idxs), dtype=bool)
        for v in V:
            if abs(v) > 1e-12 or ramp.Vc != 0.0:
                ramp.to(-v)
            tr = cb.solve(ramp.psi, V_cathode=-v, G=G, tol=gummel_tol, poisson_tol=tol)
            out.append(np.clip(np.abs(tr.J_cathode) / J_gen, 0.0, 1.0))
            conv &= ramp.converged & tr.converged
        return np.array(out).T, conv

    return solve


def eval_cce(config, values, space, *, voltages=None, gummel_tol=1e-8,
             epi_quantum_um=None, refine=0, **_):  # fmt: skip
    """CCE of the alpha-like generation in the n- epi (``run_cce``'s model) at each bias,
    with per-device mobility factors and alpha range. ``gummel_tol`` 1e-8 (the facade
    uses 1e-6): the spreads near saturation are ~1e-5."""
    from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM

    V = (
        np.linspace(-1.0, -30.0, 10)
        if voltages is None
        else np.asarray(voltages, float)
    )
    configs, snap = sample_configs(config, values, space, epi_quantum_um)
    s_n = _model_values(space, values, "mu_n_scale", 1.0)
    s_p = _model_values(space, values, "mu_p_scale", 1.0)
    rng_cm = (
        _model_values(space, values, "alpha_range_um", ALPHA_RANGE_SIC_CM * 1e4) * 1e-4
    )
    Y, info = _run_groups(
        configs, V.size, _cce_solve(V, s_n, s_p, rng_cm, gummel_tol), refine
    )
    info["epi_snap_max_um"] = snap
    return EvalOut(V, Y, int(np.sum(~np.isfinite(Y).all(axis=1))), info)


def epi_energy_fraction(range_um, epi_um, n=4001):
    """Fraction of the alpha's deposited energy inside the first ``epi_um`` of its path
    (``alpha_depth_shape``, the CCE generation shape, integrated on a fine grid)."""
    from etna.core.generation_profiles import alpha_depth_shape

    out = []
    for r, e in zip(
        np.atleast_1d(range_um),
        np.broadcast_to(epi_um, np.shape(np.atleast_1d(range_um))),
    ):
        s = np.linspace(0.0, 1.3 * r, n) * 1e-4
        f = alpha_depth_shape(s, alpha_range_cm=r * 1e-4)
        tot = np.trapezoid(f, s)
        inside = s <= e * 1e-4
        out.append(np.trapezoid(f[inside], s[inside]) / tot if tot > 0 else np.nan)
    return np.array(out)


def eval_alpha_charge(config, values, space, *, voltages=None, E0_MeV=5.486, **kw):
    """Collected charge (fC) of one Am-241 alpha entering at the junction: CCE x q x
    (E0 x f_epi / E_pair), f_epi the share of its energy deposited in the epi."""
    from etna.core.batched.material import Q
    from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM
    from etna.core.sic_material import SiC4H_Parameters

    out = eval_cce(config, values, space, voltages=voltages, **kw)
    E_pair = _model_values(space, values, "E_pair_eV", SiC4H_Parameters().E_pair_eV)
    rng_um = _model_values(space, values, "alpha_range_um", ALPHA_RANGE_SIC_CM * 1e4)
    configs, _ = sample_configs(config, values, space, kw.get("epi_quantum_um"))
    epi = np.array([c.epi_thickness_um for c in configs])
    f = epi_energy_fraction(rng_um, epi)
    charge_fC = Q * (E0_MeV * 1e6 * f / E_pair) * 1e15
    out.Y = out.Y * charge_fC[:, None]
    out.info["epi_energy_fraction_nominal"] = float(epi_energy_fraction(
        [ALPHA_RANGE_SIC_CM * 1e4], [config.epi_thickness_um])[0])  # fmt: skip
    return out


def eval_dark_current(config, values, space, *, voltages=None, epi_quantum_um=None,
                      refine=0, **_):  # fmt: skip
    """Total dark current (A, for the config area) at each bias: the TAT/SRH model of
    ``run_dark_current`` (``etna.core.batched.dark.dark_sweep``), per-device N_t and
    mobility factors."""
    from etna.core.batched.dark import dark_sweep
    from etna.core.batched.material import PARAMS
    from etna.core.dark_current import E_REF_TAT, tat_constants
    from etna.core.sic_material import intrinsic_concentration, srh_lifetime

    V = np.linspace(0.0, -60.0, 7) if voltages is None else np.asarray(voltages, float)
    configs, snap = sample_configs(config, values, space, epi_quantum_um)
    s_n = _model_values(space, values, "mu_n_scale", 1.0)
    s_p = _model_values(space, values, "mu_p_scale", 1.0)

    def solve(pb, x, kws, idxs):
        cs = [configs[i] for i in idxs]
        T = np.array([c.T for c in cs])
        n_i, _, _, E_g = zip(*[intrinsic_concentration(t, PARAMS) for t in T])
        consts = [tat_constants(t, ni, eg, PARAMS, E_t=c.E_t, N_t=c.N_t)
                  for t, ni, eg, c in zip(T, n_i, E_g, cs)]  # fmt: skip
        col = lambda k: np.array([d[k] for d in consts])
        sub = kws[0]["substrate_thickness_cm"]
        epi_end = (
            sub + kws[0]["epi_thickness_cm"]
            if kws[0].get("nplus_thickness_cm")
            else None
        )
        res = dark_sweep(
            x, pb.donors, pb.acceptors, T, [c.N_A for c in cs], np.array(n_i), V,
            tau_n=[srh_lifetime(t, "electron", PARAMS) for t in T],
            tau_p=[srh_lifetime(t, "hole", PARAMS) for t in T],
            n1=col("n1_tat"), p1=col("p1_tat"), G0=col("G0"), F_gamma=col("F_gamma"),
            E_hurkx_max=col("E_hurkx_max"), E_ref=E_REF_TAT, epi_end=epi_end,
            S_n=PARAMS.S_n, S_p=PARAMS.S_p, mobility=getattr(pb, "mobility", None),
            mu_n_scale=s_n[idxs], mu_p_scale=s_p[idxs],
        )  # fmt: skip
        area = np.array([c.area_cm2 for c in cs])
        return (res.J_total * area[None, :]).T, res.converged.all(axis=0)

    Y, info = _run_groups(configs, V.size, solve, refine)
    info["epi_snap_max_um"] = snap
    return EvalOut(V, Y, int(np.sum(~np.isfinite(Y).all(axis=1))), info)


def eval_transit_time(config, values, space, *, voltages=None, carrier="electron",
                      epi_quantum_um=None, refine=0, **_):  # fmt: skip
    """Low-field drift transit time (s) across the depleted epi (``run_transit_time``'s
    integral t = sum h / (mu |E|) over the epi edges with |E| > 1e3 V/cm), per device
    mobility factors."""
    from etna.core.batched.sweeps import _Ramp

    if carrier not in ("electron", "hole"):
        raise ValueError("carrier must be 'electron' or 'hole'")
    V = (
        np.linspace(-5.0, -100.0, 20)
        if voltages is None
        else np.asarray(voltages, float)
    )
    configs, snap = sample_configs(config, values, space, epi_quantum_um)
    s = _model_values(
        space, values, "mu_n_scale" if carrier == "electron" else "mu_p_scale", 1.0
    )
    order = np.argsort(np.abs(V))

    def solve(pb, x, kws, idxs):
        cb = _scaled_continuity(pb, kws, np.ones(len(idxs)), np.ones(len(idxs)))
        mu = (cb.mu_n if carrier == "electron" else cb.mu_p) * s[idxs][
            None, :
        ]  # (N-1, B)
        sub = kws[0]["substrate_thickness_cm"]
        epi_end = sub + kws[0]["epi_thickness_cm"]
        xe = 0.5 * (x[1:] + x[:-1])
        in_epi = ((xe > sub) & (xe < epi_end))[:, None]
        h = np.diff(x)[:, None]
        T = min(float(kws[0]["T"]), 300.0)
        ramp = _Ramp(pb, 1e-10, 10.0 * T / 300.0)
        out = np.full((len(idxs), V.size), np.nan)
        for k in order:
            mag = abs(float(V[k]))
            if mag > 0 or ramp.Vc != 0.0:
                ramp.to(mag)
            E = np.abs(np.diff(ramp.psi, axis=1)).T / h  # (N-1, B)
            sel = in_epi & (E > 1e3)
            with np.errstate(divide="ignore", invalid="ignore"):
                out[:, k] = np.sum(np.where(sel, h / (mu * E), 0.0), axis=0)
        return out, ramp.converged

    Y, info = _run_groups(configs, V.size, solve, refine)
    info["epi_snap_max_um"] = snap
    return EvalOut(V, Y, int(np.sum(~np.isfinite(Y).all(axis=1))), info)


def eval_breakdown(config, values, space, *, impact_ionization="hatakeyama2004",
                   epi_quantum_um=None, refine=0, **_):  # fmt: skip
    """Avalanche breakdown voltage |V| (V), McIntyre integral = 1, one impact-ionization
    set (``etna.core.avalanche.breakdown_voltage``, batched)."""
    from etna.core.avalanche import breakdown_voltage
    from etna.core.physics_models import PhysicsModels

    configs, snap = sample_configs(config, values, space, epi_quantum_um)
    models = PhysicsModels(impact_ionization=impact_ionization)

    def solve(pb, x, kws, idxs):
        Ts = {float(k["T"]) for k in kws}
        if len(Ts) != 1:  # one coefficient set per batch
            raise NotImplementedError(
                "breakdown UQ: T must not vary within a mesh group"
            )
        bd = breakdown_voltage(pb, models.impact_coefficients(Ts.pop()))
        return bd.V_BD[:, None], bd.converged

    Y, info = _run_groups(configs, 1, solve, refine)
    info["epi_snap_max_um"] = snap
    info["impact_ionization"] = impact_ionization
    return EvalOut(np.array([0.0]), Y, int(np.sum(~np.isfinite(Y).all(axis=1))), info)


# ------------------------------------------------------------------ damage


def _damage_constants(config, proton_energy_MeV):
    from etna.api.simulation import _kappa, _min_epi_donors
    from etna.core.trap_presets import radmanovac_calibrated

    m = radmanovac_calibrated()
    return dict(kappa=_kappa(proton_energy_MeV), R_lit=m.neutral_removal_rate - m.donor_removal,
                d0=m.donor_removal, min_nd=_min_epi_donors(config))  # fmt: skip


def eval_compensation_fluence(config, values, space, *, proton_energy_MeV=62.0, **_):
    """Proton fluence (p/cm^2) at which the calibrated trap model compensates the
    lowest-doped epi: min N_D / (kappa (R_lit + d)) (``compensation_fluence``), closed
    form, with r_n moving kappa by s = 1.069 / r_n and d re-fitted (tied_removal)."""
    from etna.uq.params import R_N_NOMINAL

    c = _damage_constants(config, proton_energy_MeV)
    configs, _ = sample_configs(config, values, space, 0)
    from etna.api.simulation import _min_epi_donors

    nd = np.array([_min_epi_donors(cf) for cf in configs])
    d = _model_values(space, values, "d_removal", c["d0"])
    r_n = _model_values(space, values, "r_n", R_N_NOMINAL)
    s = R_N_NOMINAL / r_n
    d_tied = tied_removal(d, r_n)
    phi = nd / (s * c["kappa"] * (c["R_lit"] + d_tied))
    return EvalOut(np.array([proton_energy_MeV]), phi[:, None],
                   int(np.sum(~np.isfinite(phi))),
                   {"kappa_nominal": c["kappa"], "R_lit_per_neq": c["R_lit"]})  # fmt: skip


def eval_cce_fluence(config, values, space, *, fluences=None, V_bias=-40.0,
                     proton_energy_MeV=62.0, epi_quantum_um=None, **_):  # fmt: skip
    """CCE at ``V_bias`` vs proton fluence (``run_radiation_damage`` with the calibrated
    kinetic trap preset, batched engine only, unconverged points NaN). One engine sweep
    per sample (~0.3 s): the slowest observable. r_n scales the neq fluence by
    s = 1.069 / r_n with d re-fitted (tied_removal)."""
    from etna.api.simulation import default_fluence_grid, run_radiation_damage
    from etna.core.trap_presets import radmanovac_calibrated
    from etna.uq.params import R_N_NOMINAL

    if fluences is None:
        fluences = default_fluence_grid(config, proton_energy_MeV, V_bias=V_bias)
    fluences = np.asarray(fluences, dtype=float)
    configs, snap = sample_configs(config, values, space, epi_quantum_um)
    d0 = radmanovac_calibrated().donor_removal
    d = _model_values(space, values, "d_removal", d0)
    r_n = _model_values(space, values, "r_n", R_N_NOMINAL)
    s = R_N_NOMINAL / r_n
    d_tied = tied_removal(d, r_n)
    from etna.api.simulation import validated_domain_reason

    Y = np.full((len(configs), fluences.size), np.nan)
    outside = unconv = 0
    for i, c in enumerate(configs):
        if validated_domain_reason(c) is not None:
            outside += 1
            continue
        model = radmanovac_calibrated(slow_removal=float(d_tied[i]))
        r = run_radiation_damage(dataclasses.replace(c, radiation_traps=model),
                                 fluences=fluences * s[i], V_bias=V_bias,
                                 proton_energy_MeV=proton_energy_MeV)  # fmt: skip
        per = [str(v) for v in r.metadata.get("solver_per_point", [])]
        y = np.asarray(r.y, dtype=float)
        bad = np.array([("devsim" in p) for p in per] or [False] * y.size)
        y[bad] = np.nan
        unconv += int(np.sum(~np.isfinite(y)))
        Y[i] = y
    return EvalOut(fluences, Y, int(np.sum(~np.isfinite(Y).all(axis=1))),
                   {"n_outside_box": outside, "n_unconverged_points": unconv,
                    "epi_snap_max_um": snap, "V_bias": V_bias})  # fmt: skip


# ------------------------------------------------------------------ microdosimetry


def _event_energies(mc_csv_path):
    from etna._paths import data_path
    from etna.core.mc_coupling import load_mc_events_csv

    path = mc_csv_path
    if path is None:  # the repository's synthetic events (sdist / checkout only)
        path = data_path("synthetic_mc_events.csv")
        if not path.is_file():
            raise ValueError(
                "microdosimetry UQ: pass mc_csv_path= (the default "
                "data/synthetic_mc_events.csv is not shipped in the wheel)"
            )
    ev = load_mc_events_csv(path)
    e = ev.groupby("event_id")["edep_keV"].sum().to_numpy(dtype=float)
    return e[e > 0]


def eval_microdosimetry(config, values, space, *, mc_csv_path=None, sv_width_um=100.0,
                        sv_depth_um=None, chord_model="cauchy", **_):  # fmt: skip
    """(y_F, y_D) in keV/um from the MC events (``run_microdosimetry``'s pipeline),
    exact event moments y_F = <eps>/l, y_D = <eps^2>/(<eps> l) (the facade's binned
    values differ by the 50-bins-per-decade histogram); SV thickness = the epi thickness.
    The deposits are held fixed (events not re-transported for each thickness)."""
    from etna.core.microdosimetry import mean_chord_length

    e = _event_energies(mc_csv_path)
    m1, m2 = e.mean(), (e**2).mean()
    configs, _ = sample_configs(config, values, space, 0)
    depth = sv_width_um if sv_depth_um is None else sv_depth_um
    lbar = np.array([
        mean_chord_length(c.epi_thickness_um, sv_width_um=sv_width_um, sv_depth_um=depth)
        if chord_model == "cauchy" else mean_chord_length(c.epi_thickness_um)
        for c in configs
    ])  # fmt: skip
    Y = np.column_stack([m1 / lbar, (m2 / m1) / lbar])
    return EvalOut(np.array([0.0, 1.0]), Y, 0, {"n_events": int(e.size)})


OBSERVABLES = {
    "v_fd": Observable(
        "v_fd",
        "",
        "full-depletion |V| (V)",
        eval_v_fd,
        1e-3,
        "full-depletion voltage (run_cv rule, interpolated)",
    ),
    "cv": Observable(
        "cv", "reverse bias (V)", "C (F)", eval_cv, 1e-16, "small-signal C-V"
    ),
    "cce": Observable(
        "cce",
        "reverse bias (V)",
        "CCE",
        eval_cce,
        1e-6,
        "alpha CCE vs bias (run_cce model)",
    ),
    "alpha_charge": Observable(
        "alpha_charge",
        "reverse bias (V)",
        "collected charge (fC)",
        eval_alpha_charge,
        1e-6,
        "collected charge of one Am-241 alpha",
    ),
    "dark_current": Observable(
        "dark_current",
        "reverse bias (V)",
        "I (A)",
        eval_dark_current,
        1e-17,
        "dark current (TAT/SRH)",
        log_y=True,
    ),
    "compensation_fluence": Observable(
        "compensation_fluence",
        "",
        "Phi_c (p/cm^2)",
        eval_compensation_fluence,
        0.0,
        "full-compensation proton fluence (closed form)",
    ),
    "cce_fluence": Observable(
        "cce_fluence",
        "proton fluence (p/cm^2)",
        "CCE",
        eval_cce_fluence,
        1e-6,
        "CCE vs proton fluence",
    ),
    "transit_time": Observable(
        "transit_time",
        "reverse bias (V)",
        "transit time (s)",
        eval_transit_time,
        0.0,
        "drift transit time",
        log_y=True,
    ),
    "breakdown": Observable(
        "breakdown", "", "V_BD (V)", eval_breakdown, 1e-2, "avalanche breakdown voltage"
    ),
    "microdosimetry": Observable(
        "microdosimetry",
        "",
        "lineal energy (keV/um)",
        eval_microdosimetry,
        0.0,
        "y_F and y_D",
        output_names=("y_F", "y_D"),
    ),
}
