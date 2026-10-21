"""Worked example, chapter "device": stack, 4H-SiC material constants, Poisson/contacts,
ionization, V_bi / depletion width / C-V, thermal constants.

Every number below is computed by ETNA's own functions, with a hand formula alongside.
Run: cd <worktree> && uv run python <this file>
"""

import numpy as np

from etna import DeviceConfig
from etna.api.simulation import run_cv
from etna.core.analytical import (
    EPS_0,
    Q,
    built_in_potential,
    depletion_width,
    full_depletion_voltage_profile,
)
from etna.core.batched.poisson import contact_values
from etna.core.device import DEFAULT_STACK, K_B_J, donor_profile
from etna.core.incomplete_ionization import (
    ALUMINUM_ACHATZ_2008,
    NITROGEN_KAGAMIHARA_2004,
    ionized_acceptor_fraction,
    neutral_ionized_acceptor_concentration,
)
from etna.core.incomplete_ionization import (
    ionization_fields,
    neutral_state,
)
from etna.core.sic_material import (
    SiC4H_Parameters,
    bandgap,
    compute_ni,
    effective_dos,
    intrinsic_concentration,
)
from etna.core.track.thermal import C_V, line_source_peak


def section(t):
    print("\n=== " + t + " ===")


def main():
    p = SiC4H_Parameters()
    T = 300.0
    Vt = K_B_J * T / Q  # the solvers' thermal voltage
    kT = 8.617e-5 * T  # eV, the k_B of compute_ni / analytical.py

    section("1. material constants at 300 K")
    ni_raw, NC, NV, Eg = compute_ni(T)
    ni_cal = intrinsic_concentration(T, p)[0]
    hand_ni = np.sqrt(NC * NV) * np.exp(-Eg / (2 * kT))
    print(
        f"E_g(300)  code {bandgap(T):.6f} eV   hand 3.2965625-6.5e-4*300^2/1600 = "
        f"{3.2965625 - 6.5e-4 * 300**2 / 1600:.6f}"
    )
    print(f"E_g(350)  code {bandgap(350.0):.6f} eV ; E_g(400) {bandgap(400.0):.6f} eV")
    print(f"N_C(300)  code {NC:.4e}  (tabulated NC_300 {p.NC_300:.3e})")
    print(f"N_V(300)  code {NV:.4e}  (tabulated NV_300 {p.NV_300:.3e})")
    print(f"n_i raw sqrt(NcNv)exp(-Eg/2kT): code {ni_raw:.4e}  hand {hand_ni:.4e}")
    print(
        f"n_i used by solvers (anchored to n_i_300={p.n_i_300:g}): {ni_cal:.4e}"
        f"  -> anchor factor n_i_300/raw = {p.n_i_300 / ni_raw:.4f}"
    )
    for TT in (250.0, 350.0, 400.0):
        print(f"  n_i({TT:.0f} K) used = {intrinsic_concentration(TT, p)[0]:.4e} cm^-3")
    print(
        f"eps = eps_r*eps0 = {p.eps_r}*{EPS_0:g} = {p.eps_r * EPS_0:.5e} F/cm ; "
        f"V_t = k_B T/q = {Vt:.6f} V"
    )

    section("2. p+ Al ionization (legacy scalar, the DEFAULT) at 1e19, 300 K")
    NA = 1e19
    f_code = ionized_acceptor_fraction(NA, T)
    EA = 0.220 - 3.4e-8 * NA ** (1 / 3)
    f_gibbs = 1 / (1 + 4 * np.exp(EA / kT))
    w = 1 / (1 + np.exp(-(np.log10(NA) - 18) / 0.5))
    f_emp = (0.05 + 0.10 * (np.log10(NA) - 18)) * (1 + 0.5 * (T - 300) / 300)
    f_hand = (1 - w) * f_gibbs + w * f_emp
    print(
        f"code f = {f_code:.6f} ; hand (1-w) f_Gibbs + w f_emp = (1-{w:.4f})*{f_gibbs:.3e}"
        f" + {w:.4f}*{f_emp:.3f} = {f_hand:.6f} ; N_A- = {NA * f_code:.4e}"
    )
    print(f"legacy fraction at 400 K: {ionized_acceptor_fraction(NA, 400.0):.4f}")

    section("3. carrier-dependent ionization (opt-in), neutral region, 300 K")
    NCc, NVc = effective_dos(T)
    lv = ALUMINUM_ACHATZ_2008[0]
    p1 = NVc * np.exp(-float(lv.energy(NA)) / Vt) / 4.0  # p1/g_A
    # neutral p+: p = N_A- = N_A / (1 + p/p1)  ->  p^2 + p1 p - p1 N_A = 0
    p_hand = 0.5 * (-p1 + np.sqrt(p1 * p1 + 4 * p1 * NA))
    print(f"Al 1e19: dE = 0.220 - 2.32e-8*N^(1/3) = {float(lv.energy(NA)):.5f} eV")
    print(
        f"  code N_A-/N_A = {neutral_ionized_acceptor_concentration(NA, T) / NA:.6f}"
        f" ; hand quadratic = {p_hand / NA:.6f}"
    )
    for ND in (DEFAULT_STACK["N_D_junction"], DEFAULT_STACK["buffer_end_doping"], 1e19):
        f = ionization_fields(np.array([ND]), np.array([0.0]), T)
        n_i = intrinsic_concentration(T, p)[0]
        Dp = neutral_state(f, n_i)[2][0]
        # hand: n = sum_l (N/2) / (1 + n / (Nc exp(-dE_l/Vt)/2)), bisection on n
        refs = [
            NCc * np.exp(-float(l.energy(ND)) / Vt) / 2.0
            for l in NITROGEN_KAGAMIHARA_2004
        ]
        lo, hi = 0.0, ND
        for _ in range(200):
            n = 0.5 * (lo + hi)
            r = sum(0.5 * ND / (1 + n / rf) for rf in refs) - n
            lo, hi = (n, hi) if r > 0 else (lo, n)
        print(
            f"N {ND:.3e}: code N_D+/N_D = {Dp / ND:.5f} ; hand two-site = {n / ND:.5f}"
        )

    section("4. built-in potential (default stack, legacy ionization)")
    NA_ion = NA * f_code
    NDj = DEFAULT_STACK["N_D_junction"]
    Vbi_an = built_in_potential(NA_ion, NDj, ni_cal, T=T)
    print(
        f"analytical V_bi = kT ln(N_A- N_Dj / n_i^2) = {Vbi_an:.4f} V "
        f"(N_A- {NA_ion:.3e}, N_Dj {NDj:.4e})"
    )
    psi_a, psi_c = contact_values(
        np.array([-NA_ion]), np.array([1e19]), ni_cal, Vt, 0.0, 0.0
    )
    psi_a, psi_c = float(np.ravel(psi_a)[0]), float(np.ravel(psi_c)[0])
    print(
        f"solver contacts: psi_anode {float(psi_a):.4f} V, psi_cathode(n+ 1e19) "
        f"{float(psi_c):.4f} V -> total built-in {float(psi_c - psi_a):.4f} V"
    )
    print(
        f"  hand V_t ln(N_A- N_D+ / n_i^2) with n+ 1e19: "
        f"{Vt * np.log(NA_ion * 1e19 / ni_cal**2):.4f} V "
        f"(= p+/n- junction + n-/n+ high-low steps)"
    )

    section("5. depletion width: analytic forms vs run_cv (batched engine) ")
    eps = p.eps_r * EPS_0
    xj = 0.0
    y = np.linspace(0.0, 10e-4, 200001)
    nd = donor_profile(
        y,
        xj,
        NDj,
        DEFAULT_STACK["N_D_bulk"],
        DEFAULT_STACK["L_transition"],
        buffer=(
            10e-4,
            DEFAULT_STACK["buffer_thickness_cm"],
            DEFAULT_STACK["buffer_end_doping"],
        ),
    )
    moment = np.concatenate(
        [[0.0], np.cumsum(0.5 * (y[1:] * nd[1:] + y[:-1] * nd[:-1]) * np.diff(y))]
    )
    V_of_W = Q / eps * moment  # V_bi - V needed to deplete to W (abrupt, one-sided)

    def W_graded(Vbias):
        return float(np.interp(Vbi_an - Vbias, V_of_W, y))

    res = run_cv(DeviceConfig(), v_start=0.0, v_stop=-10.0, n_points=41)
    V = res.x
    WC = res.metadata["W_from_C"]
    We = res.metadata["depletion_widths"]
    C = res.y
    A = res.metadata["area_cm2"]
    print(
        f"run_cv solver = {res.metadata.get('solver')}, area {A} cm^2, "
        f"reach_through_V = {res.metadata.get('reach_through_V')}"
    )
    print(
        " V(V)  eps*A/C(um) W_edge(um) | uniform N_Dj(um) uniform N_Db(um) graded-moment(um)"
    )
    for Vb in (0.0, -0.5, -1.0, -2.0, -5.0, -10.0):
        k = int(np.argmin(np.abs(V - Vb)))
        Wu_j = depletion_width(Vbi_an, Vb, NDj, eps_r=p.eps_r, epi_thickness=10e-4)
        Wu_b = depletion_width(
            Vbi_an, Vb, DEFAULT_STACK["N_D_bulk"], eps_r=p.eps_r, epi_thickness=10e-4
        )
        print(
            f"{V[k]:6.2f} {WC[k]*1e4:10.4f} {We[k]*1e4:10.4f} | {Wu_j*1e4:14.4f} "
            f"{Wu_b*1e4:16.4f} {W_graded(Vb)*1e4:16.4f}"
        )
    k0 = int(np.argmin(np.abs(V)))
    print(f"C(0 V) = {C[k0]:.4e} F ; eps*A/W(0)= {eps * A / WC[k0]:.4e} F (identity)")
    vfd = full_depletion_voltage_profile(
        y[y <= 10e-4 - DEFAULT_STACK["buffer_thickness_cm"]],
        nd[y <= 10e-4 - DEFAULT_STACK["buffer_thickness_cm"]],
        Vbi_an,
    )
    print(f"analytic V_fd (moment integral to transition start 8.863 um) = {vfd:.3f} V")

    section("6. Mott-Schottky / differential doping from the run_cv C(V)")
    from etna.measurements.cv import CVData, depletion_profile

    data = CVData(
        V=V,
        C_per_area=C / A,
        sigma_C_per_area=None,
        sigma_V=None,
        area_cm2=A,
        area_source="run_cv",
        frequency_Hz=None,
        T_K=300.0,
        source="run_cv(DeviceConfig())",
        sha256="",
    )
    dep = depletion_profile(data, window=3)  # etna.measurements.cv: 3-point local fit
    inv = 1.0 / (C / A) ** 2  # (F/cm^2)^-2
    dinv = np.gradient(inv, V)  # d(1/C'^2)/dV  (negative: V negative = reverse)
    N_app = 2.0 / (Q * eps * (-dinv))  # hand formula, central differences
    print(" V(V)  W=epsA/C(um)  N(W) code  N(W) hand  L_D(um)  N_D(x=W) of the profile")
    for Vb in (-0.25, -0.5, -1.0, -1.5, -2.0, -3.0, -5.0, -8.0):
        k = int(np.argmin(np.abs(V - Vb)))
        print(
            f"{V[k]:6.2f} {dep.depth_um[k]:10.4f} {dep.N_cm3[k]:10.3e} {N_app[k]:10.3e}"
            f" {dep.debye_um[k]:7.3f} {float(np.interp(WC[k], y, nd)):12.3e}"
        )

    section("7. thermal constants (track engine, opt-in)")
    print(f"C_v = c_p*rho = 0.69*{p.rho} = {C_V:.4f} J/cm^3/K")
    E_per_cm = 1000 * 1e3 * 1.602e-19 / 1e-4  # LET 1000 keV/um -> J/cm
    for t in (1e-12, 1e-11, 1e-10):
        print(
            f"line source LET 1000 keV/um, w 10 nm, t {t:.0e} s: dT axis = "
            f"{line_source_peak(E_per_cm, 10e-7, t):.2f} K "
            f"(hand E/(pi c_v (w^2+4 a t)) = "
            f"{E_per_cm / (np.pi * C_V * ((10e-7) ** 2 + 4 * 4.71 / C_V * t)):.2f})"
        )


if __name__ == "__main__":
    main()
