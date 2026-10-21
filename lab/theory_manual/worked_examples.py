"""Worked examples of docs/THEORY.md: one function per chapter, each a few seconds.

Every example evaluates the code path the chapter documents and, next to it, the same
quantity from the chapter's formula written out by hand, so the manual's equations are
checked against the code. Run all (about 10 s on an M-series laptop):

    uv run python lab/theory_manual/worked_examples.py
    uv run python lab/theory_manual/worked_examples.py materials niel   # a subset

The printed numbers are quoted in docs/THEORY.md (notes/LOG.md keeps the last run).
"""

from __future__ import annotations

import math
import sys
import warnings

import numpy as np

Q = 1.602176634e-19  # C
K_B_EV = 8.617333262e-5  # eV/K
EPS0 = 8.8541878128e-14  # F/cm


def materials():
    """Ch. 2: n_i(T), N_C, N_V, Varshni gap and the low-field mobility sets."""
    from etna.core.physics_models import LOW_FIELD_MOBILITY
    from etna.core.sic_material import (
        SiC4H_Parameters,
        bandgap,
        compute_ni,
        intrinsic_concentration,
        mobility_caughey_thomas_T,
    )

    p = SiC4H_Parameters()
    ni_raw, NC, NV, Eg = compute_ni(300.0)
    print(f"E_g(300 K) = {bandgap(300.0):.6f} eV (Varshni, E_g0 3.2965625)")
    print(
        f"N_C = {NC:.4e}, N_V = {NV:.4e} cm^-3; sqrt(NC NV) exp(-Eg/2kT) = {ni_raw:.3e}"
    )
    print(
        f"solver n_i(300) = {intrinsic_concentration(300.0)[0]:.3e} (anchored to n_i_300)"
    )
    print(f"n_i(400 K) = {intrinsic_concentration(400.0)[0]:.3e} cm^-3")
    for N in (3.382506e13, 1e17, 1e19):
        mu_c = mobility_caughey_thomas_T(N, 300.0, "electron")
        hand = 10.0 + (1190.0 - 10.0) / (1.0 + (N / 1e17) ** 0.54)
        mu_l = LOW_FIELD_MOBILITY["legacy_perp"].mobility(N, 300.0, "electron")
        mu_h = mobility_caughey_thomas_T(N, 300.0, "hole")
        print(
            f"N = {N:.3e}: mu_n c_axis {mu_c:.2f} (hand {hand:.2f}), legacy_perp "
            f"{float(mu_l):.2f}, mu_p {mu_h:.2f} cm2/Vs"
        )
    mu_350 = mobility_caughey_thomas_T(3.382506e13, 350.0, "electron")
    print(f"mu_n(3.38e13, 350 K) = {mu_350:.1f} cm2/Vs ((350/300)^-2.40 on mu_max)")
    print(f"eps = {p.eps_r} eps0, rho = {p.rho} g/cm3, E_pair = {p.E_pair_eV} eV")


def scharfetter_gummel():
    """Ch. 4: the Bernoulli function and the SG flux limits."""
    from etna.core.batched.transport import bernoulli

    x = np.array([-50.0, -1.0, -1e-9, 0.0, 1e-9, 1.0, 50.0, 800.0])
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        exact = np.where(x == 0, 1.0, x / np.expm1(x))
    print("x        B(x) code        x/(e^x-1)")
    for xi, b, e in zip(x, bernoulli(x), exact):
        print(f"{xi:8.1e} {b:.15e} {e:.15e}")
    d = np.array([0.3, -2.0, 7.0])
    print("B(-x) - B(x) - x =", bernoulli(-d) - bernoulli(d) - d, "(identity, 0)")
    # SG flux of a uniform density in a linear potential: pure drift, J = q mu n E
    mu, n, h, V_t, dpsi = 1000.0, 1e14, 1e-5, 0.025852, 0.01
    D = dpsi / V_t
    Jn = (
        Q
        * mu
        * V_t
        / h
        * (n * bernoulli(np.array([-D]))[0] - n * bernoulli(np.array([D]))[0])
    )
    print(
        f"SG Jn (uniform n, dpsi {dpsi} V over {h*1e4} um) = {Jn:.6e} A/cm2; "
        f"q mu n E = {Q * mu * n * dpsi / h:.6e}"
    )


def recombination():
    """Ch. 6: SRH lifetimes, Hurkx Gamma, the TAT constants, Auger and Hecht."""
    from etna.core.charge_collection import hecht_cce
    from etna.core.dark_current import hurkx_gamma, tat_constants
    from etna.core.sic_material import (
        SiC4H_Parameters,
        intrinsic_concentration,
        mobility_caughey_thomas_T,
        srh_lifetime,
    )

    p = SiC4H_Parameters()
    ni, _, _, Eg = intrinsic_concentration(300.0)
    c = tat_constants(300.0, ni, Eg, p)
    hb, m0 = 1.054571817e-34, 9.1093837015e-31
    kT = K_B_EV * 300.0 * Q
    FG = math.sqrt(24 * 0.25 * m0 * kT**3) / (Q * hb) / 100.0
    print(
        f"F_Gamma code {c['F_gamma']:.4e} V/cm, hand {FG:.4e}; E_hurkx_max "
        f"{c['E_hurkx_max']:.3e} V/cm"
    )
    for E in (5e4, 1e5, 1.7e5):
        f = E / c["F_gamma"]
        print(
            f"Gamma({E:.1e}) = {float(hurkx_gamma(E, c['F_gamma'])):.3f}  "
            f"(hand 2 sqrt(3 pi) f exp(f^2) = {2*math.sqrt(3*math.pi)*f*math.exp(f*f):.3f})"
        )
    print(
        f"n1_tat = {c['n1_tat']:.3e}, p1_tat = {c['p1_tat']:.3e} cm^-3 (E_t = Ec - 0.67)"
    )
    print(f"N_t = {c['N_t']:.6e} cm^-3 s^-1 (G0 at 300 K = {c['G0']:.6e})")
    for T in (300.0, 350.0):
        print(
            f"tau_n({T:.0f}) = {srh_lifetime(T, 'electron'):.3e} s, "
            f"tau_p = {srh_lifetime(T, 'hole'):.3e} s"
        )
    N = 1e19
    print(
        f"Auger minority-hole lifetime in 1e19 n+: 1/(C_n N^2) = {1/(p.C_n*N*N):.2e} s"
    )
    # Hecht (depth-averaged two-carrier form)
    V, d = 40.0, 10e-4
    cce = float(hecht_cce(V, d))
    nd = 3.382506e13  # bulk epi doping of the default stack (mobility_doping)
    mu_e = mobility_caughey_thomas_T(nd, 300.0, "electron")
    mu_h = mobility_caughey_thomas_T(nd, 300.0, "hole")
    le, lh = mu_e * 1e-9 * V / d, mu_h * 6e-7 * V / d
    h = sum((l / d) - (l / d) ** 2 * (1 - math.exp(-d / l)) for l in (le, lh))
    print(
        f"Hecht CCE(40 V, 10 um) code {cce:.6f}, hand {h:.6f} (lambda_e "
        f"{le*1e4:.0f} um, lambda_h {lh*1e4:.0f} um)"
    )


def avalanche():
    """Ch. 7: impact-ionisation coefficients and the parallel-plane breakdown."""
    from etna import DeviceConfig
    from etna.api.models import run_breakdown
    from etna.core.physics_models import IMPACT_IONIZATION

    m = IMPACT_IONIZATION["hatakeyama2004"]
    for E in (1e6, 2e6, 3e6):
        an = float(m.alpha(E, carrier="electron"))
        ap = float(m.alpha(E, carrier="hole"))
        print(
            f"E = {E:.0e} V/cm: alpha_n {an:.3e} (hand {1.76e8*math.exp(-3.30e7/E):.3e}),"
            f" alpha_p {ap:.3e} (hand {3.41e8*math.exp(-2.50e7/E):.3e}) 1/cm"
        )
    bd = run_breakdown(DeviceConfig())
    md = bd.metadata
    print(
        f"default dosimeter V_BD = {md['V_BD']:.2f} V, mesh_error_V = "
        f"{md['mesh_error_V']:.3f} V, E_max = {md['E_max_at_breakdown']:.4e} V/cm"
    )


def generation():
    """Ch. 8: dose rate -> G, NIST stopping in SiC, alpha range and Bragg curve."""
    from etna.core.generation_profiles import (
        ALPHA_RANGE_SIC_CM,
        dose_rate_to_generation,
        proton_generation_profile,
    )
    from etna.core.stopping_power import (
        alpha_bragg_summary,
        bohr_range_straggling_cm,
        linear_stopping,
        mass_stopping,
    )

    G1 = float(dose_rate_to_generation(1.0))
    print(
        f"G per 1 Gy/s in SiC: {G1:.4e} cm^-3 s^-1 (hand 3.21e4/(7.8*1.602e-12) = "
        f"{3.21e4/(7.8*1.602e-12):.4e})"
    )
    for E in (5.6, 62.0, 150.0):
        L = float(linear_stopping("proton", "sic", E))  # MeV/cm
        s_si = float(mass_stopping("proton", "si", E))
        s_c = float(mass_stopping("proton", "c_amorphous", E))
        hand = (0.7004 * s_si + 0.2996 * s_c) * 3.21
        print(
            f"p {E:6.1f} MeV: L_SiC = {L*0.1:.4f} keV/um (Bragg hand {hand*0.1:.4f}),"
            f" {L*1e6/7.8*1e-4:.1f} pairs/um"
        )
    G = proton_generation_profile(np.array([0.0, 5e-4]), 62.0, fluence_rate_cm2_s=1e8)
    print(f"62 MeV, 1e8 p/cm2/s: G = {G[0]:.4e} cm^-3 s^-1 (fluence mode)")
    s = alpha_bragg_summary()
    print(
        f"alpha 5.486 MeV: CSDA range {s['csda_range_um']:.3f} um "
        f"(ALPHA_RANGE_SIC_CM {ALPHA_RANGE_SIC_CM*1e4:.3f}), Bohr sigma_R "
        f"{s['bohr_straggling_um']:.3f} um, peak/entrance {s['peak_to_entrance']:.3f}"
        f" at {s['peak_depth_um']:.2f} um, pairs {5.486e6/7.8:.0f}"
    )
    print(
        f"proton 62 MeV Bohr range straggling {bohr_range_straggling_cm('proton', 62.0)*1e4:.1f} um"
    )


def cce():
    """Ch. 9-10: the steady-state DD CCE of the default device and the LET table."""
    import json

    from etna import DeviceConfig, run_cce
    from etna._paths import data_path

    r = run_cce(DeviceConfig(), v_start=-10.0, v_stop=-40.0, n_points=4)
    for v, y in zip(r.x, r.y):
        print(f"CCE({v:.0f} V) = {y:.6f}")
    print("solver:", r.metadata.get("solver"))
    t = json.loads(data_path("cce_let_table_100um.json").read_text())
    L, C = np.array(t["LET_keV_um"]), np.array(t["CCE"])
    qg = np.array(t["Q_generated_fC"])
    print(
        f"CCE(LET) table 100 um, {t['bias_V']} V, mobility {t['mobility']}: CCE "
        f"{C.min():.6f}..{C.max():.6f} over LET {L.min():.2g}..{L.max():.3g} keV/um"
    )
    print(
        f"Q_generated/LET = {qg[0]/L[0]:.5f} fC per keV/um (constant: {np.ptp(qg/L):.1e})"
    )


def flash():
    """Ch. 11: FLASH dose-rate CCE (sensitivity bound)."""
    from etna.core.flash_recombination import cce_vs_dose_rate

    r = cce_vs_dose_rate(np.array([20.0, 230.0]))
    c = r["cce_values"]
    print(
        f"CCE(20 Gy/s) = {c[0]:.9f}, CCE(230 Gy/s) = {c[1]:.9f}, no-Auger ref "
        f"{r['cce_no_auger_ref']:.9f} ({r['solver']})"
    )
    print(
        f"dose-rate change {c[1]-c[0]:.1e}; Auger offset {r['cce_no_auger_ref']-c[0]:.2e}"
    )


def niel():
    """Ch. 12: hardness factors in the three NIEL conventions, by hand from the tables."""
    from etna.core.niel_tables import (
        NIEL_NEUTRON_1MEV_SIC_MEV_CM2_G,
        NIEL_NEUTRON_1MEV_SI_MEV_CM2_G,
        R_N_SIC_SI_FISSION,
    )
    from etna.core.radiation_damage import (
        NIEL_PROTON_SI,
        NIEL_PROTON_SIC,
        get_hardness_factor,
        hardness_table,
    )

    print(
        f"NIEL_n,Si(1 MeV) = 95 MeV mb = {NIEL_NEUTRON_1MEV_SI_MEV_CM2_G:.5e} MeV cm2/g;"
        f" NIEL_n,SiC(1 MeV) = {NIEL_NEUTRON_1MEV_SIC_MEV_CM2_G:.5e}; r_n = "
        f"{R_N_SIC_SI_FISSION:.4f}"
    )
    E = 62.0
    si, sic = NIEL_PROTON_SI[E], NIEL_PROTON_SIC[E]
    hand = {
        "si-neq": si / NIEL_NEUTRON_1MEV_SI_MEV_CM2_G,
        "sic-damage": sic / (R_N_SIC_SI_FISSION * NIEL_NEUTRON_1MEV_SI_MEV_CM2_G),
        "sic-1mev": sic / NIEL_NEUTRON_1MEV_SIC_MEV_CM2_G,
    }
    print(f"62 MeV proton: NIEL Si {si:.5e}, SiC {sic:.5e} MeV cm2/g")
    for conv, h in hand.items():
        k = get_hardness_factor(E, hardness_table("proton", conv))
        print(
            f"  kappa[{conv}] = {k:.4f} (hand {h:.4f}); 1e12 p/cm2 -> {k*1e12:.3e} neq/cm2"
        )


def damage():
    """Ch. 12: trapping constants, trap emission, removal and compensation."""
    from etna import DeviceConfig
    from etna.api.simulation import compensation_fluence, legacy_compensation_fluence
    from etna.core.radiation_damage import RadiationDamageParams, compute_K_tau
    from etna.core.radiation_traps import _vth
    from etna.core.sic_material import effective_dos
    from etna.core.trap_presets import CALIBRATED_SLOW_REMOVAL_BY_CONVENTION

    prm = RadiationDamageParams()
    Kn, Kp = compute_K_tau(prm, "electron"), compute_K_tau(prm, "hole")
    vn = _vth(300.0, 0.77)
    print(f"v_th,n = {vn:.4e} cm/s; K_tau,n = {Kn:.4e}, K_tau,p = {Kp:.4e} cm2/s")
    print(f"hand K_tau,n = {(5.0*2e-14 + 1.6*1.4e-14 + 2.4*5e-13)*vn:.4e}")
    NC, _ = effective_dos(300.0)
    V_t = K_B_EV * 300.0
    for name, E, s in (
        ("Z1/2", 0.67, 2e-14),
        ("EH4", 1.03, 5e-13),
        ("EH6/7", 1.60, 1.4e-14),
    ):
        en = s * vn * NC * math.exp(-E / V_t)
        print(
            f"{name}: e_n = sigma v_th N_C exp(-E/kT) = {en:.3e} /s, 1/e_n = {1/en:.3e} s"
        )
    print(
        "calibrated direct removal d (per neq):", CALIBRATED_SLOW_REMOVAL_BY_CONVENTION
    )
    cfg = DeviceConfig()
    print(
        f"compensation fluence, 62 MeV p: default {compensation_fluence(cfg, 62.0):.3e}, "
        f"legacy scalar + si-neq {legacy_compensation_fluence(cfg, 62.0, niel_convention='si-neq'):.3e} p/cm2"
    )


def microdosimetry():
    """Ch. 14-15: mean chord, kappa, moments, straggling width."""
    from etna.core.mc_response import lineal_moments
    from etna.core.microdosimetry import mean_chord_length
    from etna.core.stopping_power import kappa_table, sic_electron_density_cm3

    lb = mean_chord_length(10.0, 100.0, 100.0)
    w, t, d = 100.0, 10.0, 100.0
    print(
        f"l_bar(100x100x10 um) = {lb:.4f} um (4V/S = {4*w*t*d/(2*(w*t+w*d+t*d)):.4f}),"
        f" slab 2t = {mean_chord_length(10.0):.1f} um"
    )
    k = kappa_table("proton", "water")
    for E in (1.0, 10.0, 62.0, 150.0):
        print(
            f"kappa_p(SiC->water, {E:g} MeV) = {np.interp(E, k['energy_MeV'], k['kappa']):.4f}"
        )
    rng = np.random.default_rng(1)
    eps = rng.gamma(4.0, 2.5, 20000)  # keV, synthetic deposits
    m = lineal_moments(eps, 10.0)
    yF, yD = eps.mean() / 10.0, (eps**2).mean() / eps.mean() / 10.0
    print(
        f"y_F = {m['y_F']:.4f} (hand {yF:.4f}) +- {m['se_y_F']:.4f}, y_D = "
        f"{m['y_D']:.4f} (hand {yD:.4f}) +- {m['se_y_D']:.4f} keV/um"
    )
    Ne = sic_electron_density_cm3()
    e2 = 1.439964e-13  # MeV cm
    om = math.sqrt(4 * math.pi * e2**2 * Ne * 10e-4) * 1e3
    from etna.core.optimization import estimate_noise_floor

    nf = estimate_noise_floor(18e-12, 1e-6, 10.0, 100.0, 100.0)
    hand = 3 * math.sqrt(2 * Q * 18e-12 * 1e-6) * 1e15
    print(
        f"shot-noise floor, 18 pA, 1 us: Q_min {nf['Q_min_fC']:.4f} fC (3 sigma, hand "
        f"{hand:.4f}), E_min {nf['E_min_keV']:.3f} keV, y_min {nf['y_min_keV_um']:.4f} keV/um"
    )
    print(
        f"N_e(SiC) = {Ne:.4e} cm^-3; Bohr sigma of a proton over 10 um SiC = {om:.3f} keV"
    )


def tracks():
    """Ch. 9.4: radial track models conserve the pairs."""
    from etna.core.track import radial

    r = np.array([1e-7, 1e-6, 1e-5, 1e-4])
    for name, F in (
        ("gaussian(50 nm)", radial.gaussian(5e-6)),
        ("core_penumbra(1 nm, 1 um)", radial.core_penumbra(1e-7, 1e-4)),
    ):
        print(name, "F(r) =", np.round(F.cdf(r), 5))
    from etna import DeviceConfig
    from etna.core.track.phd import MESH_FAST, simulate_track

    ev = simulate_track(
        DeviceConfig(), LET_keV_um=1000.0, V_bias=50.0, radial=radial.gaussian(50e-7),
        mesh=MESH_FAST,
    )
    print(f"track LET 1000 keV/um, 50 V, gaussian 50 nm: CCE {ev.cce:.6f}, recombined "
          f"{ev.recombined_fraction:.3e}, ledger residual {ev.ledger_residual:.1e}")


def device():
    """Ch. 1-3 and 13: stack, ionization, V_bi, W(V), C-V inversion, thermal constants
    (lab/theory_manual/device_example.py)."""
    import importlib.util
    from pathlib import Path

    path = Path(__file__).with_name("device_example.py")
    spec = importlib.util.spec_from_file_location("device_example", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.main()


def numerics():
    """Ch. 16: engine vs devsim parity of the small-signal C on the default device."""
    from etna import DeviceConfig, run_cv
    from etna.api.simulation import _run_cv_devsim

    cfg = DeviceConfig()
    eng = run_cv(cfg, v_start=0.0, v_stop=-10.0, n_points=3)
    dev = _run_cv_devsim(cfg, np.array([0.0, -5.0, -10.0]))
    print("solvers:", eng.metadata.get("solver"), "/", dev.metadata.get("solver"))
    for v, ce, cd in zip(eng.x, eng.y, dev.y):
        print(f"V = {v:6.1f}: C engine {ce:.10e} F, devsim {cd:.10e} F, rel diff {abs(ce-cd)/cd:.1e}")


CHAPTERS = {
    "device": device,
    "materials": materials,
    "sg": scharfetter_gummel,
    "recombination": recombination,
    "avalanche": avalanche,
    "generation": generation,
    "cce": cce,
    "flash": flash,
    "niel": niel,
    "damage": damage,
    "microdosimetry": microdosimetry,
    "tracks": tracks,
    "numerics": numerics,
}


def main(names=None):
    warnings.simplefilter("ignore")
    for name in names or CHAPTERS:
        print(f"\n## {name}")
        CHAPTERS[name]()


if __name__ == "__main__":
    main(sys.argv[1:] or None)
