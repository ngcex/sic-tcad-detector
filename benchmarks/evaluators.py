"""Evaluators of the benchmark registry (benchmarks/registry.toml).

Each function takes no argument and returns ``{"S": code value, "D": reference value
(only when the reference is computed, e.g. an analytic formula or a devsim run),
"details": {...}}``. Fixed references (literature, digitised data) live in the registry,
not here. Every evaluator reuses the public or core API exactly as the tests and lab
benches do; the provenance of each comparison is in the registry entry.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

EPS_0 = 8.854187817e-14  # F/cm
EPS_SIC = 9.7 * EPS_0
Q = 1.602176634e-19
K_B_EV = 8.617333262e-5


# ---------------------------------------------------------------- helpers
def _quiet():
    from etna.core.devsim_lock import quiet_devsim

    return quiet_devsim()


def _reset_devsim():
    from etna.core.devsim_reset import reset_devsim_fully

    with _quiet():
        reset_devsim_fully()


@lru_cache(maxsize=1)
def _uniform_cv():
    """Uniform p+/n- diode, N_D 1e15, 30 um epi, 300 K: run_cv (batched) 0..-40 V."""
    from etna import DeviceConfig, run_cv

    cfg = DeviceConfig(doping_profile="uniform", N_D=1e15, epi_thickness_um=30.0)
    r = run_cv(cfg, 0.0, -40.0, 41)
    assert r.metadata["solver"] == "batched"
    return cfg, r


def _v_bi(cfg):
    from etna.core.analytical import built_in_potential
    from etna.core.incomplete_ionization import ionized_acceptor_concentration
    from etna.core.sic_material import SiC4H_Parameters, intrinsic_concentration

    n_i = intrinsic_concentration(cfg.T, SiC4H_Parameters())[0]
    na_i = float(ionized_acceptor_concentration(cfg.N_A, cfg.T))
    return built_in_potential(na_i, cfg.N_D, n_i, T=cfg.T), na_i


# ------------------------------------------------- VERIFIED: analytic limits
def an_depletion_width_uniform():
    """eps*A/C at -40 V vs the one-sided abrupt junction with the 2kT/q correction
    W = sqrt(2 eps (V_bi - V - 2kT/q) (1 + N_D/N_A^-) / (q N_D)) (Sze & Ng Ch. 2)."""
    cfg, r = _uniform_cv()
    v_bi, na_i = _v_bi(cfg)
    i = int(np.argmin(np.abs(r.x + 40.0)))
    W_c = EPS_SIC * cfg.area_cm2 / r.y[i]
    kT = K_B_EV * cfg.T
    W_a = np.sqrt(
        2 * EPS_SIC * (v_bi - r.x[i] - 2 * kT) * (1 + cfg.N_D / na_i) / (Q * cfg.N_D)
    )
    return {"S": W_c * 1e4, "D": W_a * 1e4, "details": {"V_bi": v_bi, "V": r.x[i]}}


def an_mott_schottky_nd():
    """N_D from the 1/C^2 slope (-5..-40 V) of the uniform diode vs its input N_D."""
    cfg, r = _uniform_cv()
    m = r.x <= -5.0
    slope, _ = np.polyfit(r.x[m], 1.0 / r.y[m] ** 2, 1)
    nd = -2.0 / (Q * EPS_SIC * cfg.area_cm2**2 * slope)
    return {"S": nd, "D": cfg.N_D}


def _const_field(alpha_n, alpha_p, W=10e-4, N=201):
    x = np.linspace(0.0, W, N)
    psi = (-3e6 * x)[None, :]
    return x, psi, (alpha_n, 1e-9, 0.0, 1.0, alpha_p, 1e-9, 0.0, 1.0)


def an_mcintyre_alpha_eq_beta():
    """Constant field, alpha = beta: ionization integral I = alpha W and M(x) = 1/(1 - aW)
    (McIntyre 1966). alpha W = 0.5 -> M = 2 everywhere; S = max |M(x) - 2| + 2."""
    from etna.core.avalanche import ionization_integrals

    a = 0.5 / 10e-4
    x, psi, coef = _const_field(a, a)
    I = ionization_integrals(x, psi, coef)
    dev = float(np.max(np.abs(I.M_x[0] - 2.0)))
    return {"S": 2.0 + dev, "D": 2.0, "details": {"I_n": float(I.I_n[0])}}


def an_mcintyre_electrons_only():
    """Constant field, only electrons ionize: M_n = exp(alpha W) (closed form)."""
    from etna.core.avalanche import ionization_integrals

    a = 0.5 / 10e-4
    x, psi, coef = _const_field(a, 0.0)
    I = ionization_integrals(x, psi, coef)
    return {"S": float(I.M_n[0]), "D": float(np.exp(0.5))}


def an_breakdown_uniform_field():
    """p+ / 1e12 / n+ (nearly constant field), alpha = beta = a exp(-b/E): breakdown at
    a W exp(-b W / V) = 1, V_BD = b W / ln(a W). The field reaches ~0.1 um into p+ and
    n+, so W_eff > W: an O(1 %) analytic-model remainder."""
    from etna.core.avalanche import breakdown_voltage
    from etna.core.batched.mesh import node_positions
    from etna.core.batched.sweeps import make_batch

    SUB, NPL, W = 0.3e-4, 1e-4, 10e-4
    x = node_positions(W, SUB, NPL)
    pb = make_batch(
        x, SUB, 1e19, 1e12, 1e12, 1e-4, 300.0, nplus=(SUB + W, 1e19), buffer=None
    )
    a, b = 3.0e6, 1.5e7
    bd = breakdown_voltage(pb, (a, b, 0.0, 1.0, a, b, 0.0, 1.0))
    return {"S": float(bd.V_BD[0]), "D": float(b * W / np.log(a * W))}


# ------------------------------------------------ CROSS-CHECKED: code vs code
def cc_batched_vs_devsim_cv():
    """run_cv (batched) vs devsim cv_sweep, default stack, 0..-20 V on the 2.5 V
    lattice: S = 1 + max |C_engine / C_devsim - 1| (D = 1)."""
    from etna import DeviceConfig
    from etna.api import simulation as sim
    from etna.api.device import build_device
    from etna.core.cv_analysis import cv_sweep

    cfg = DeviceConfig()
    res = sim.run_cv(cfg, v_start=0.0, v_stop=-20.0, n_points=9)
    assert res.metadata["solver"] == "batched"
    _reset_devsim()
    try:
        with _quiet():
            ref = cv_sweep(
                build_device(cfg), np.linspace(0.0, -20.0, 9), area=cfg.area_cm2
            )
    finally:
        _reset_devsim()
    dc = float(np.max(np.abs(res.y / ref["capacitance"] - 1)))
    dw = float(np.max(np.abs(
        np.asarray(res.metadata["depletion_widths"]) / ref["depletion_widths"] - 1)))  # fmt: skip
    return {"S": 1.0 + dc, "D": 1.0, "details": {"max_rel_dC": dc, "max_rel_dW": dw}}


def cc_batched_vs_devsim_cce():
    """run_cce (batched) vs devsim cce_vs_bias, default stack, -10..-40 V:
    S = max |dCCE| (abs, D = 0)."""
    from etna import DeviceConfig
    from etna.api import simulation as sim
    from etna.api.device import device_kwargs
    from etna.core.charge_collection import cce_vs_bias

    cfg = DeviceConfig()
    res = sim.run_cce(cfg, v_start=-10.0, v_stop=-40.0, n_points=4)
    assert res.metadata["solver"] == "batched"
    kw = device_kwargs(cfg)
    _reset_devsim()
    try:
        with _quiet():
            ref = cce_vs_bias(
                np.linspace(-10.0, -40.0, 4),
                epi_thickness_cm=kw.pop("epi_thickness_cm"),
                device_kwargs=kw,
            )
    finally:
        _reset_devsim()
    return {"S": float(np.max(np.abs(res.y - ref["cce_values"]))), "D": 0.0,
            "details": {"cce_engine": res.y}}  # fmt: skip


def cc_batched_vs_devsim_dark():
    """run_dark_current engine vs solver='devsim', default stack, 0..-100 V (6 points):
    S = 1 + max |I_engine / I_devsim - 1|."""
    from etna import DeviceConfig
    from etna.api import simulation as sim

    cfg = DeviceConfig()
    e = sim.run_dark_current(cfg, 0.0, -100.0, 6)
    with _quiet():
        d = sim.run_dark_current(cfg, 0.0, -100.0, 6, solver="devsim")
    _reset_devsim()
    assert e.metadata["solver"] == "batched" and d.metadata["solver"] == "devsim"
    dev = float(np.max(np.abs(e.y / d.y - 1)))
    return {"S": 1.0 + dev, "D": 1.0, "details": {"I_engine_A": e.y}}


def cc_engine2d_vs_devsim_potential():
    """2D run_field, default stack, half width 25 um, -20 V: engine2d vs devsim on the
    same mesh, S = max |dpsi| (V)."""
    from etna import DeviceConfig
    from etna.api import simulation as sim

    cfg = DeviceConfig(half_width_um=25.0)
    e = sim.run_field(cfg, -20.0)
    with _quiet():
        d = sim.run_field(cfg, -20.0, solver="devsim")
    _reset_devsim()
    assert e.metadata["solver"] == "engine2d" and d.metadata["solver"] == "devsim"
    dpsi = np.abs(e.mesh.node_values["Potential"] - d.mesh.node_values["Potential"])
    Ee, Ed = e.mesh.node_values["ElectricField"], d.mesh.node_values["ElectricField"]
    return {"S": float(dpsi.max()), "D": 0.0,
            "details": {"max_rel_dE": float(np.max(np.abs(Ee - Ed)) / Ed.max())}}  # fmt: skip


def _geometry2d_refs():
    import json
    from pathlib import Path

    notes = Path(__file__).resolve().parents[1] / "lab" / "geometry2d" / "notes"
    return json.loads((notes / "devsim_baselines.json").read_text())


def cc_fv2d_vs_devsim_mesa():
    """fv2d on devsim's mesa mesh (legacy depth mesh), 50 V, 9-point cathode scan vs the
    stored devsim CCE (lab/geometry2d/notes/devsim_baselines.json): S = max |dCCE|."""
    from etna.core.alternative_structures import create_mesa_device
    from etna.core.fv2d import cce_lateral_scan_fv2d, mesh_from_builder, prepared_solver

    ref = _geometry2d_refs()["mesa"]
    with _quiet():
        mesh, _ = mesh_from_builder(create_mesa_device, depth_mesh="legacy")
    _reset_devsim()
    s = prepared_solver(mesh, {"cathode": ref["V"]})
    sc = cce_lateral_scan_fv2d(s, n_points=len(ref["cce"]))
    return {"S": float(np.max(np.abs(sc["cce_values"] - np.asarray(ref["cce"])))),
            "D": 0.0}  # fmt: skip


def cc_fv2d_vs_engine2d():
    """Default (n+) laterally uniform diode, hw 25 um, 20 V: fv2d on devsim's mesh vs
    engine2d, 4-point lateral CCE scan, S = max |dCCE|."""
    from etna.core.device2d import create_sic_2d_device
    from etna.core.engine2d import lateral_cce_scan_2d
    from etna.core.fv2d import cce_lateral_scan_fv2d, mesh_from_builder, prepared_solver

    hw, V, n = 25.0, 20.0, 4
    e2 = lateral_cce_scan_2d(hw, V, n_lateral_points=n)
    with _quiet():
        mesh, _ = mesh_from_builder(create_sic_2d_device, half_width_um=hw)
    _reset_devsim()
    s = prepared_solver(mesh, {"cathode": V})
    sc = cce_lateral_scan_fv2d(s, n_points=n)
    return {"S": float(np.max(np.abs(sc["cce_values"] - e2["cce_values"]))), "D": 0.0}


def cc_avalanche_dd_vs_mcintyre():
    """Coupled DD with Hatakeyama impact ionization vs the generation-weighted McIntyre
    multiplication on the same potential, default dosimeter at 0.99 V_BD:
    S = (1 - 1/gain) / (1 - 1/M) (D = 1)."""
    from etna import DeviceConfig
    from etna.api.models import _batch
    from etna.core.avalanche import breakdown_voltage, ionization_integrals
    from etna.core.batched.models_dd import ModelContinuityBatch
    from etna.core.batched.sweeps import _Ramp, alpha_generation, srh_continuity
    from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM
    from etna.core.physics_models import PhysicsModels

    SUB = 0.3e-4
    cfg = DeviceConfig()
    pb, x, kw = _batch(cfg)
    G = alpha_generation(x, SUB, ALPHA_RANGE_SIC_CM, 1e18, epi_end=SUB + 10e-4)
    cb = srh_continuity(pb, 300.0, kw["N_A"])
    coef = PhysicsModels(impact_ionization="hatakeyama2004").impact_coefficients(300.0)
    V = 0.99 * breakdown_voltage(pb, coef).V_BD[0]
    ramp = _Ramp(pb, 1e-10, 10.0)
    ramp.to(V)
    base = ModelContinuityBatch(pb, cb.mu_n, cb.mu_p, cb.tau_n, cb.tau_p)
    mc = ModelContinuityBatch(pb, cb.mu_n, cb.mu_p, cb.tau_n, cb.tau_p, ii=coef)
    r0 = base.solve_models(ramp.psi, V_cathode=V, G=G)
    r = mc.solve_models(ramp.psi, V_cathode=V, G=G)
    gain = r.J[0] / r0.J[0]
    M = ionization_integrals(x, ramp.psi, coef).generation_weighted(x, G)[0]
    return {"S": float((1 - 1 / gain) / (1 - 1 / M)), "D": 1.0,
            "details": {"gain": float(gain), "M": float(M), "V": float(V)}}  # fmt: skip


def cc_csda_integral_vs_nist_si():
    """ETNA's 1/S CSDA integration (the path used for SiC, which NIST lacks) vs NIST's
    own tabulated CSDA range in Si, 62 MeV protons: S/D - 1."""
    from etna.core import stopping_power as sp

    t = sp.load_star("proton", "si")
    E = 62.0
    nist = float(np.exp(np.interp(np.log(E), np.log(t["energy_MeV"]),
                                  np.log(t["csda_range_g_cm2"]))))  # fmt: skip
    e, r = sp._csda_table("proton", "si")
    return {"S": float(np.interp(E, e, r)), "D": nist}


def cc_csda_integral_vs_nist_si_alpha():
    """As above for 5.5 MeV alphas (ASTAR Si)."""
    from etna.core import stopping_power as sp

    t = sp.load_star("alpha", "si")
    E = 5.5
    nist = float(np.exp(np.interp(np.log(E), np.log(t["energy_MeV"]),
                                  np.log(t["csda_range_g_cm2"]))))  # fmt: skip
    e, r = sp._csda_table("alpha", "si")
    return {"S": float(np.interp(E, e, r)), "D": nist}


def cc_kappa_62mev_independent():
    """SiC -> water linear LET ratio kappa at 62 MeV from etna.core.stopping_power vs an
    independent recomputation from the NIST CSVs (log-log interpolation, Bragg weights
    0.7004 Si + 0.2996 C amorphous, electronic stopping, rho_SiC 3.21 g/cm3). Checks the table pipeline; the
    Bragg-additivity approximation itself is NOT validated (data/nist_star/README.md).
    """
    from etna.core import stopping_power as sp

    E = 62.0
    k = sp.kappa_table("proton")
    S = float(np.exp(np.interp(np.log(E), np.log(k["energy_MeV"]), np.log(k["kappa"]))))

    def s_el(material):
        t = sp.load_star("proton", material)
        return float(np.exp(np.interp(np.log(E), np.log(t["energy_MeV"]),
                                      np.log(t["S_electronic_MeV_cm2_g"]))))  # fmt: skip

    sic = 0.7004 * s_el("si") + 0.2996 * s_el("c_amorphous")
    D = s_el("water") * 1.0 / (sic * 3.21)
    return {"S": S, "D": D}


# ------------------------------------------- CALIBRATED: agreement with fit data
def ca_photons_plateau_rms():
    """Default stack eps*A/C vs the digitised Photons Fig. 4b plateau (-5..-59 V,
    lab/cv_gap/notes/photons_w.csv band points, weight >= 10): rms residual (um)."""
    plateau, _ = _photons_digitised()
    W_at = _default_w_at()
    res = np.array([w - W_at(v) for v, w in plateau])
    return {"S": float(np.sqrt(np.mean(res**2))), "D": 0.0,
            "details": {"n_points": len(plateau), "max_abs_um": float(np.max(np.abs(res)))}}  # fmt: skip


@lru_cache(maxsize=1)
def _default_w_at():
    from etna.api import simulation as sim
    from etna.api.device import measured_cv_config

    cfg = measured_cv_config(area_cm2=0.04)
    V = np.round(np.arange(0.0, -60.001, -0.25), 2)
    r = sim.run_cv(cfg, 0.0, -60.0, len(V))
    W = EPS_SIC * cfg.area_cm2 / r.y * 1e4
    return lambda v: W[np.argmin(np.abs(r.x - v))]


def _photons_digitised():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "lab" / "cv_gap" / "notes"
    band = np.genfromtxt(root / "photons_w.csv", delimiter=",", names=True)
    mk = np.genfromtxt(root / "photons_w_markers.csv", delimiter=",", names=True)
    skip = (-10.0, -30.0, -9.75, -10.25, -10.5, -9.5, -29.75, -30.25)
    plateau = [
        (float(v), float(w))
        for v, w, wt in zip(band["V"], band["W_um"], band["weight"])
        if np.isfinite(w) and wt >= 10 and v <= -5.0 and v not in skip
    ]
    knee = [(round(float(v) * 4) / 4, float(w)) for v, w in zip(mk["V"], mk["W_um"])
            if -3.3 < v < -0.1]  # fmt: skip
    return plateau, knee


@lru_cache(maxsize=1)
def _w_creep():
    from etna.viz.validation import w_creep

    return w_creep().metadata


def ca_w_label_m10():
    """eps*A/C at -10 V (w_creep, batched) vs the printed label 9.52 um."""
    md = _w_creep()
    return {"S": float(md["W_C_at_labels_um"][1]), "D": float(md["label_W_um"][1])}


def ca_w_label_0():
    """eps*A/C at 0 V vs the printed label 1.68 um."""
    md = _w_creep()
    return {"S": float(md["W_C_at_labels_um"][0]), "D": float(md["label_W_um"][0])}


def ca_w_label_m30():
    """eps*A/C at -30 V vs the printed label 9.73 um."""
    md = _w_creep()
    return {"S": float(md["W_C_at_labels_um"][2]), "D": float(md["label_W_um"][2])}


def ca_dark_current_m60():
    """run_dark_current (batched), default stack, I(-60 V, 4 mm2) vs the fitted
    "< 18 pA up to -60 V" bound: N_t is fitted to it, so E ~ 0 by construction."""
    from etna import DeviceConfig
    from etna.api import simulation as sim

    r = sim.run_dark_current(DeviceConfig(), 0.0, -60.0, 7)
    return {"S": float(r.y[-1]), "details": {"N_t": r.metadata["N_t"]}}


def ca_radmanovac_onsemi_slope():
    """Calibrated kinetic preset (default damage model), onsemi 50 um 4e13 device of
    Radmanovac et al. Table 2: C-V removal slope g at fixed width 10-40 um (per p/cm2),
    3 of the 5 published fluences (linear regime)."""
    from etna.core.trap_calibration import RemovalObservation, removal_observable
    from etna.core.trap_presets import trap_preset

    o = RemovalObservation(
        "onsemi", "slope", 5.60, 0.81, 4e13, 50.0, 252.7, (0.0, 1.4e12, 2.8e12),
        (10.0, 40.0), 80.0, 5.0, 1e6,
    )  # fmt: skip
    return {"S": float(removal_observable(o, trap_preset("radmanovac_calibrated")))}


def _hatakeyama(name):
    from etna.core.avalanche import breakdown_voltage
    from etna.core.batched.mesh import node_positions
    from etna.core.batched.sweeps import make_batch
    from etna.core.physics_models import PhysicsModels

    SUB, NPL = 0.3e-4, 1e-4
    x = node_positions(10e-4, SUB, NPL)
    ND = np.array([3e16, 1e17])
    pb = make_batch(
        x, SUB, 1e19, ND, ND, 1e-4, 300.0, nplus=(SUB + 10e-4, 1e19), buffer=None
    )
    coef = PhysicsModels(impact_ionization=name).impact_coefficients(300.0)
    bd = breakdown_voltage(pb, coef)
    assert bd.converged.all()
    return dict(zip(("3e16", "1e17"), bd.V_BD.tolist()))


def ca_hatakeyama_vbd_1e17():
    """hatakeyama2004 set, p+/10 um/n+ at N_D 1e17: V_BD vs Hatakeyama's measured fit
    1940 (1e16/N_D)^0.8 V. The set was fitted on these diodes: consistency."""
    return {"S": _hatakeyama("hatakeyama2004")["1e17"]}


def va_niwa_vbd_3e16():
    """niwa2015 set (independent of Hatakeyama's diodes), N_D 3e16: V_BD vs the
    measured fit 1940 (1e16/N_D)^0.8 V."""
    return {"S": _hatakeyama("niwa2015")["3e16"]}


def va_niwa_vbd_1e17():
    return {"S": _hatakeyama("niwa2015")["1e17"]}


# ------------------------- CALIBRATED: solver scale on the Bruzzi diode (circular)
def ca_bruzzi_w40_vbi_matched():
    """Uniform p-n configured to the Bruzzi & Verroi 2023 Schottky epi (N_D 2.46e15,
    30 um, r 1 mm, 293.15 K): eps*A/C read at the same total band bending as the
    paper's -40 V (V_R + V_bi, V_bi Schottky 1.2 V held in their fit) vs W = 4.2 um.
    Same procedure as scripts/validation/cv_vs_bruzzi2023.py. CIRCULAR: N_D 2.46e15 is
    the paper's own fit of the same C-V from which W comes, so this checks the solver
    scale, not a model prediction (.planning/VALIDATION_LITERATURE.md)."""
    from etna import DeviceConfig, run_cv

    cfg = DeviceConfig(doping_profile="uniform", N_D=2.46e15, epi_thickness_um=30.0,
                       area_cm2=np.pi * 0.1**2, T=293.15)  # fmt: skip
    r = run_cv(cfg, v_start=0.0, v_stop=-50.0, n_points=51)
    eps = 9.7 * 8.854e-14  # the script's constants
    W_um = eps * cfg.area_cm2 / r.y * 1e4
    inv = 1.0 / r.y**2
    m = r.x <= -5.0
    slope, icpt = np.polyfit(r.x[m], inv[m], 1)
    v_bi_pn = -icpt / slope
    v_r = -r.x
    order = np.argsort(v_r)
    W_matched = float(np.interp(40.0 + 1.2 - v_bi_pn, v_r[order], W_um[order]))
    W40 = float(W_um[np.argmin(np.abs(r.x + 40.0))])
    return {"S": W_matched, "details": {"W_at_minus40_um": W40, "V_bi_pn": v_bi_pn}}


def ca_bruzzi_w40_raw():
    """As above without the V_bi matching: eps*A/C at -40 V vs 4.2 um (the Schottky vs
    p-n built-in voltage difference dominates E)."""
    return {"S": ca_bruzzi_w40_vbi_matched()["details"]["W_at_minus40_um"]}


# -------------------------------------------- EXPECTED-FAIL: documented misses
def ef_photons_knee_rms():
    """Default stack eps*A/C vs the 13 digitised knee markers (-3.25..-0.25 V): rms (um).
    The exponential epi cannot make the 1.3e13 dip at 4-6 um (lab/buffer_profile LOG §3).
    """
    _, knee = _photons_digitised()
    W_at = _default_w_at()
    res = np.array([w - W_at(v) for v, w in knee])
    return {"S": float(np.sqrt(np.mean(res**2))), "D": 0.0,
            "details": {"n_points": len(knee)}}  # fmt: skip


@lru_cache(maxsize=1)
def _tilted_shape():
    from etna.viz.validation import tilted_alpha_shape

    return tilted_alpha_shape().metadata


def ef_tilted_alpha_shape_rms():
    """run_cce_tilted_alpha(microdosimeter_config()) normalised at -50 V vs the digitised
    Micro Fig. 4b CCE (7 points): rms shape residual (etna.viz.validation)."""
    md = _tilted_shape()
    return {"S": float(md["rms"]), "D": 0.0,
            "details": {"model_at_0V": float(md["model_at_data"][0])}}  # fmt: skip


def ef_tilted_alpha_0v():
    """Same, the 0 V point: model vs data 0.9083 (normalised at -50 V)."""
    return {"S": float(_tilted_shape()["model_at_data"][0])}


# ------------------------------------------- VERIFIED: more analytic limits
def an_full_depletion_knee_1e14():
    """Uniform 1e14, 10 um epi without n+ (batched engine): the Mott-Schottky knee
    ((d/eps)^2 - b)/s of a 1/C^2 fit below full depletion vs the analytic punch-through
    voltage q N_D d^2 / (2 eps) - V_bi + 2kT/q (Sze & Ng Ch. 2 with the Debye-tail
    correction)."""
    from etna.core.analytical import built_in_potential
    from etna.core.batched import make_batch, reverse_sweep
    from etna.core.batched.material import EPS, Q as QB
    from etna.core.batched.sweeps import _n_a_ionized, _n_i

    EPI, SUB, NA, T, ND = 10e-4, 1e-4, 1e19, 300.0, 1e14
    x = np.unique(np.r_[np.linspace(0, SUB, 401), np.linspace(SUB, SUB + EPI, 2001)])
    pb = make_batch(x, SUB, NA, ND, ND, 1e-4, T, nplus=None, buffer=None)
    v_bi = built_in_potential(_n_a_ionized(NA, T), ND, _n_i(T), T=T)
    v_fd = QB * ND * EPI**2 / (2 * EPS) - v_bi
    V = np.linspace(0.0, -0.6 * v_fd, 31)
    C = reverse_sweep(pb, V, SUB, EPI, width=False).C[:, 0]
    s, b = np.polyfit(-V[5:], 1.0 / C[5:] ** 2, 1)
    knee = ((EPI / EPS) ** 2 - b) / s
    return {"S": float(knee), "D": float(v_fd + 2 * K_B_EV * T),
            "details": {"V_fd_no_kT": float(v_fd)}}  # fmt: skip


def an_hecht_closed_form():
    """hecht_cce (two-carrier, depth-averaged Hecht 1932 closed form) vs an independent
    trapezoid average of q(x)/q0 = (le/d)(1 - e^-(d-x)/le) + (lh/d)(1 - e^-x/lh), in the
    low-CCE regime (tau 2 ps, 50 V, 10 um), far from the CCE -> 1 limit. SRH
    drift-diffusion does NOT reduce to Hecht (SRH needs both carriers): this verifies
    the closed-form implementation used by run_temperature_sweep(method="hecht")."""
    from etna.core.charge_collection import hecht_cce

    d, mu_e, mu_p, tau, V = 10e-4, 900.0, 120.0, 2e-12, 50.0
    le, lh = mu_e * tau * V / d, mu_p * tau * V / d
    xs = np.linspace(0, d, 20001)
    q = (le / d) * (1 - np.exp(-(d - xs) / le)) + (lh / d) * (1 - np.exp(-xs / lh))
    return {"S": float(hecht_cce(V, d, mu_e, tau, mu_p, tau)),
            "D": float(np.trapezoid(q, xs) / d)}  # fmt: skip


def an_srh_generation_current():
    """Batched dark-current kernel, pure midgap SRH (TAT off: G0 = 0; Hurkx off), n_i
    inflated to 1e8 so the current is far above round-off; p+/1e15/10 um at -40 V vs
    q n_i W_eff / (tau_n + tau_p), W_eff = W - dn - dp the width where n, p < n_i (Sze &
    Ng Ch. 2 generation current, edge-corrected)."""
    from etna.core.batched.dark import dark_sweep
    from etna.core.batched.material import EPS, Q as QB

    SUB, EPI, ND, NAion, ni, tau = 1e-4, 10e-4, 1e15, 1e18, 1e8, 1e-6
    x = np.unique(np.r_[np.linspace(0, SUB, 201), np.linspace(SUB, SUB + EPI, 2001)])
    donors = np.where(x >= SUB, ND, 0.0)[None, :]
    acc = np.where(x < SUB, NAion, 0.0)[None, :]
    V = np.array([-40.0])
    r = dark_sweep(x, donors, acc, 300.0, 1e19, ni, V, tau_n=tau, tau_p=tau, n1=ni,
                   p1=ni, G0=0.0, F_gamma=1e30, E_hurkx_max=1e30, E_ref=3.5e4)  # fmt: skip
    kT = K_B_EV * 300.0
    v_bi = kT * np.log(NAion * ND / ni**2)
    W = np.sqrt(2 * EPS * (v_bi - V[0] - 2 * kT) / (QB * ND))
    dn = np.sqrt(2 * EPS * kT * np.log(ND / ni) / (QB * ND))
    dp = kT * np.log(NAion / ni) / (QB * ND * W / EPS)
    J_a = QB * ni * (W - dn - dp) / (2 * tau)
    naive = r.J_SRH[0, 0] / (QB * ni * W / (2 * tau))
    return {"S": float(r.J_SRH[0, 0]), "D": float(J_a),
            "details": {"ratio_to_naive_qniW_2tau": float(naive)}}  # fmt: skip


@lru_cache(maxsize=1)
def _alpha_sample():
    from etna.core import mc_transport as mt

    return mt.sample_track_events("alpha", 5.486, sv_um=(100, 100, 2.0),
                                  n_events=50_000, delta_rays=False, seed=2)  # fmt: skip


def an_bohr_variance_alpha():
    """Energy-loss sampler (Geant4 Urban/Universal port, delta rays off), 5.486 MeV
    alphas through 2 um SiC (kappa ~ 22, Gaussian regime), 50k events, seed 2: sample
    variance vs the Bohr/Vavilov variance xi T_max (1 - beta^2/2). The +3.9 % is a
    model-form difference of the Urban/Universal fluctuation port (~6x the sqrt(2/N)
    sampling error at this fixed seed), not statistics; the tolerance is the existing
    gate of tests/test_mc_transport.py (0.8 < var/var_Bohr < 1.2)."""
    from etna.core import mc_transport as mt

    ev = _alpha_sample()
    return {"S": float(ev.edep_keV.var()),
            "D": float(mt.bohr_variance_keV2("alpha", 5.486, 2.0))}  # fmt: skip


def an_mean_deposit_alpha():
    """Same sample: mean deposit vs the CSDA integral of L_el/L_tot over the path."""
    from etna.core import mc_transport as mt

    ev = _alpha_sample()
    return {
        "S": float(ev.edep_keV.mean()),
        "D": float(mt.mean_csda_deposit_keV("alpha", 5.486, 2.0)),
    }


def an_vavilov_variance():
    """Vavilov density (CERNLIB G116 form, numerical Fourier integral), kappa 0.5,
    beta^2 0.3: second cumulant vs kappa (1 - beta^2/2) (lambda in T_max units)."""
    from etna.core import straggling as St

    kappa, beta2 = 0.5, 0.3
    lam = np.linspace(-kappa * (1 + beta2 - St.EULER) - 8, 25, 6001)
    p = St.vavilov_pdf_lambda(lam, kappa, beta2)
    m1 = np.trapezoid(p * lam, lam)
    var = np.trapezoid(p * (lam - m1) ** 2, lam)
    return {"S": float(var), "D": kappa * (1 - beta2 / 2)}


def an_fv2d_coaxial_order():
    """fv2d box scheme, (r, z) Laplace between coaxial electrodes a < r < b on UNIFORM r
    meshes of 21/41/81 lines: observed order of the max potential error against
    psi = V ln(b/r) / ln(b/a) (finest pair). D = 2 (the scheme's formal order)."""
    from etna.core.fv2d import FVMesh, Solver2D, box_geometry
    from etna.core.sic_material import SiC4H_Parameters, intrinsic_concentration

    info = dict(T=300.0, n_i=intrinsic_concentration(300.0, SiC4H_Parameters())[0])

    def err(xs, a, b, V=3.0):
        ys = np.linspace(0, 1e-4, 4)
        ny, nx = len(ys), len(xs)
        vol, gh, gv = box_geometry(xs, ys, "cylindrical")
        z = np.zeros((ny, nx))
        c1 = np.zeros((ny, nx), bool)
        c1[:, 0] = True
        c2 = np.zeros((ny, nx), bool)
        c2[:, -1] = True
        mh, mv = np.full((ny, nx - 1), 100.0), np.full((ny - 1, nx), 100.0)
        m = FVMesh(xs=xs, ys=ys, vol=vol, gh=gh, gv=gv, donors=z, acceptors=z, net=z,
                   mun_h=mh, mun_v=mv, mup_h=mh, mup_v=mv,
                   contacts={"in": c1, "out": c2}, coordinate_system="cylindrical",
                   info=info)  # fmt: skip
        s = Solver2D(m).equilibrium().set_bias({"in": V})
        d = s.potential()[1] - s.potential()[1, -1]
        return np.max(np.abs(d - V * np.log(b / xs) / np.log(b / a)))

    a, b = 5e-4, 50e-4
    errs = [err(np.linspace(a, b, n), a, b) for n in (21, 41, 81)]
    order = np.log2(errs[1] / errs[2])
    return {"S": float(order), "D": 2.0,
            "details": {"errors_V": errs,
                        "geometric_mesh_error_V": float(err(np.geomspace(a, b, 21), a, b))}}  # fmt: skip


# ------------------------------- radiation damage vs literature (different devices)
def va_gsponer_cce_1e15():
    """Default damage model (calibrated kinetic traps, sic-damage NIEL), Gsponer-like
    diode (50 um, uniform 1.5e14), alpha CCE at 1e15 neq and -300 V (entered as 62 MeV
    protons, 1e15 / kappa_p; the paper's bias is not stated in its abstract) vs the
    abstract's bound "CCE better than 50 % up to 1e15 n_eq/cm2" (arXiv:2310.02047).
    E = S - 0.5 > 0 means consistent with the bound; a bound is a weak comparison."""
    from etna import DeviceConfig, run_radiation_damage
    from etna.core.radiation_damage import get_hardness_factor, hardness_table

    g = DeviceConfig(epi_thickness_um=50.0, doping_profile="uniform", N_D=1.5e14)
    kappa = get_hardness_factor(62.0, hardness_table("proton", "sic-damage"))
    r = run_radiation_damage(g, fluences=[1e15 / kappa], V_bias=-300.0,
                             proton_energy_MeV=62.0)  # fmt: skip
    return {"S": float(r.y[0]), "details": {"kappa": kappa,
                                             "convention": r.metadata.get("niel_convention")}}  # fmt: skip


@lru_cache(maxsize=1)
def _denapoli():
    from etna import DeviceConfig, run_radiation_damage

    dn = DeviceConfig(epi_thickness_um=20.0, doping_profile="uniform", N_D=1e14)
    r = run_radiation_damage(dn, fluences=[0.0, 1.31e14, 7.29e14], V_bias=-300.0,
                             proton_energy_MeV=14.0, particle="neutron")  # fmt: skip
    loss = 1.0 - np.asarray(r.y[1:]) / r.y[0]
    return loss, float(r.metadata["kappa"])


def ef_denapoli_loss_729():
    """Default damage model, 14 MeV neutrons in n/cm2 (particle="neutron", NJOY SiC NIEL,
    sic-damage), 20 um uniform 1e14 epi at -300 V, alpha CCE: relative CCE loss at
    7.29e14 n/cm2 vs the measured ~22.5 % (De Napoli, Front. Phys. 10 (2022) 898833,
    review of a 239Pu-alpha measurement on a Schottky diode)."""
    loss, kappa = _denapoli()
    return {"S": float(loss[1]), "details": {"kappa_neq_per_n": kappa}}


def ef_denapoli_loss_131():
    """As above at 1.31e14 n/cm2 vs the measured ~7.0 %."""
    loss, kappa = _denapoli()
    return {"S": float(loss[0]), "details": {"kappa_neq_per_n": kappa}}


# ------------------------------------------------------------------- slow entries
def cc_fv2d_vs_devsim_guard_ring():
    """fv2d on devsim's guard-ring mesh (legacy depth mesh), 10 V, 11-point anode scan
    with the mid-gap node-set cut vs the stored devsim CCE: S = max |dCCE|."""
    from etna.core.alternative_structures import create_guard_ring_device
    from etna.core.fv2d import (
        cce_lateral_scan_fv2d,
        cut_node_sets,
        mesh_from_builder,
        prepared_solver,
    )

    ref = _geometry2d_refs()["gr"]
    with _quiet():
        mesh, _ = mesh_from_builder(create_guard_ring_device, depth_mesh="legacy")
    _reset_devsim()
    s = prepared_solver(mesh, {"cathode": ref["V"]})
    info = s.mesh.info
    ns = cut_node_sets(s.mesh, info["sv_half_width_cm"] + 0.5 * info["gap_cm"])
    sc = cce_lateral_scan_fv2d(s, contact="anode", node_sets=ns,
                               x_positions_cm=ref["x_cm"])  # fmt: skip
    dev = float(np.max(np.abs(sc["cce_values"] - np.asarray(ref["cce_anode"]))))
    return {"S": dev, "D": 0.0}


def cc_batched_vs_devsim_cv_deep():
    """As cc_batched_vs_devsim_cv to -100 V (the deepest cross-checked bias with n+),
    21 points on the 5 V lattice."""
    from etna import DeviceConfig
    from etna.api import simulation as sim
    from etna.api.device import build_device
    from etna.core.cv_analysis import cv_sweep

    cfg = DeviceConfig()
    V = np.linspace(0.0, -100.0, 21)
    res = sim.run_cv(cfg, v_start=0.0, v_stop=-100.0, n_points=21)
    _reset_devsim()
    try:
        with _quiet():
            ref = cv_sweep(build_device(cfg), V, area=cfg.area_cm2)
    finally:
        _reset_devsim()
    dc = float(np.max(np.abs(res.y / ref["capacitance"] - 1)))
    return {"S": 1.0 + dc, "D": 1.0}
