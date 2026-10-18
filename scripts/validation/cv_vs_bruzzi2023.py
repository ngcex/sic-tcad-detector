"""Generate the C-V solver-scale check: ETNA (uniform p-n, configured to the Bruzzi &
Verroi 2023 Schottky diode's epi) vs the paper's reported W(V) endpoints.

Framing (corrected 2026-10-18, see .planning/VALIDATION_LITERATURE.md): this is a
SOLVER-SCALE check on a DIFFERENT device (a Schottky diode on another wafer, 30x lower
doping than the ETNA default), NOT a validation of the ETNA dosimeter. The +3.1 % at
-40 V is the Schottky (V_bi = 1.2 V) vs p-n (V_bi ~ 2.9 V) built-in-voltage difference:
reading the simulation at the same total band bending (V_R + V_bi), i.e. with V_bi
matched, the absolute scale agrees to < 1 %. Panel 1 overlays the abrupt-junction
W = sqrt(2 eps (V_bi + V_R) / (q N_D)) for both V_bi values to make that explicit.
Panel 2 (1/C^2 slope) is a solver self-consistency round-trip -- N_D was a simulation
INPUT, not extracted from Bruzzi's data -- so it does not validate against the paper.
"""

from pathlib import Path

from etna import DeviceConfig, run_cv
import numpy as np
import matplotlib.pyplot as plt

cfg = DeviceConfig(
    doping_profile="uniform",
    N_D=2.46e15,
    epi_thickness_um=30.0,
    area_cm2=np.pi * (0.1) ** 2,
    T=293.15,
)
# 1 V grid so -40 V is a solved point: the default linspace(0, -50, 40) has no -40 V
# (nearest -39.74 V, which read +2.8 % instead of +3.1 %; audit-v7 A-01).
result = run_cv(cfg, v_start=0.0, v_stop=-50.0, n_points=51)

eps0 = 8.854e-14
eps_r = 9.7
W_um = eps_r * eps0 * cfg.area_cm2 / result.y * 1e4
V_reverse = -result.x  # plot as positive reverse bias, matching paper convention

q = 1.602e-19
eps = eps_r * eps0
V_BI_SCHOTTKY = 1.2  # V, held in Bruzzi & Verroi's fit
W_PAPER = {0.0: 0.7, 40.0: 4.2}  # um, reported endpoints

# p-n built-in voltage of the simulated diode from the Mott-Schottky intercept
inv_C2 = 1.0 / result.y**2
mask = result.x <= -5.0
slope, intercept = np.polyfit(result.x[mask], inv_C2[mask], 1)
N_D_recovered = -2.0 / (q * eps * cfg.area_cm2**2 * slope)
V_bi_pn = -intercept / slope  # 1/C^2 = 0 at V = +V_bi (forward)


def w_abrupt_um(v_r, v_bi):
    return np.sqrt(2.0 * eps * (v_bi + v_r) / (q * cfg.N_D)) * 1e4


# V_bi-matched read-out: same total band bending V_R + V_bi as the paper at 40 V
V_R_matched = 40.0 + V_BI_SCHOTTKY - V_bi_pn
order = np.argsort(V_reverse)
W_matched = float(np.interp(V_R_matched, V_reverse[order], W_um[order]))
W_at_40 = float(W_um[np.argmin(np.abs(result.x - (-40.0)))])

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Panel 1: W(V)
ax = axes[0]
ax.plot(
    V_reverse, W_um, "b-", lw=2, label="ETNA DD solver\n(uniform p-n, same epi $N_D$)"
)
v_grid = np.linspace(0.0, 45.0, 200)
ax.plot(
    v_grid,
    w_abrupt_um(v_grid, V_bi_pn),
    "b:",
    lw=1.2,
    label=f"abrupt analytic, $V_{{bi}}$={V_bi_pn:.2f} V (p-n)",
)
ax.plot(
    v_grid,
    w_abrupt_um(v_grid, V_BI_SCHOTTKY),
    "r:",
    lw=1.2,
    label=f"abrupt analytic, $V_{{bi}}$={V_BI_SCHOTTKY:.1f} V (Schottky)",
)
ax.plot(
    list(W_PAPER.keys()),
    list(W_PAPER.values()),
    "ro",
    ms=10,
    mfc="none",
    mew=2,
    label="Bruzzi & Verroi 2023\n(Schottky, other wafer; reported endpoints)",
)
ax.set_xlabel("Reverse bias |V| (V)")
ax.set_ylabel("Depletion width W (µm)")
ax.set_xlim(0, 45)
ax.set_title("Depletion width vs bias (solver-scale check)")
ax.legend(fontsize=7, loc="lower right")
ax.grid(alpha=0.3)
ax.annotate(
    "Offset = $V_{bi}$ (p-n vs Schottky)",
    xy=(0, W_um[np.argmin(np.abs(result.x))]),
    xytext=(14, 1.4),
    fontsize=8,
    color="darkred",
    arrowprops=dict(arrowstyle="->", color="darkred", alpha=0.6),
)
ax.annotate(
    f"-40 V: {W_at_40:.3f} vs 4.2 um ({100 * (W_at_40 / 4.2 - 1):+.1f} %)\n"
    f"$V_{{bi}}$ matched: {W_matched:.3f} um ({100 * (W_matched / 4.2 - 1):+.1f} %)",
    xy=(40, W_at_40),
    xytext=(12, 4.7),
    fontsize=8,
    color="darkgreen",
    arrowprops=dict(arrowstyle="->", color="darkgreen", alpha=0.6),
)

# Panel 2: 1/C^2 vs V (Mott-Schottky) -- the junction-type-independent check
ax = axes[1]
ax.plot(V_reverse, inv_C2, "b-", lw=2, label="ETNA simulation")
ax.plot(
    V_reverse[mask],
    np.polyval([slope, intercept], result.x[mask]),
    "k--",
    lw=1,
    alpha=0.7,
    label=f"Linear fit → $N_D$={N_D_recovered:.2e} cm$^{{-3}}$\n(sim input: {cfg.N_D:.2e} cm$^{{-3}}$)",
)
ax.set_xlabel("Reverse bias |V| (V)")
ax.set_ylabel(r"1/C² (cm⁴/F²)")
ax.set_title("Mott-Schottky slope: solver self-consistency (not vs paper)")
ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)
ax.grid(alpha=0.3)

plt.suptitle(
    "Solver-scale C-V check on a different device (Schottky diode, Bruzzi & Verroi 2023,"
    " MDPI Materials 16:3643) -- V_bi matched; not a validation of the ETNA dosimeter",
    fontsize=10,
    y=1.02,
)
plt.tight_layout()
OUTPUT_PATH = (
    Path(__file__).resolve().parents[2] / "figures" / "cv_validation_bruzzi2023.png"
)
plt.savefig(OUTPUT_PATH, dpi=150, bbox_inches="tight")
print("Saved plot.")
print(
    f"W(0V) sim={W_um[np.argmin(np.abs(result.x-0.0))]:.3f} um, paper endpoint=0.7 um"
)
print(f"W(-40V) sim={W_at_40:.3f} um, paper endpoint=4.2 um ({100 * (W_at_40 / 4.2 - 1):+.1f} %)")
print(
    f"V_bi(p-n, 1/C^2 intercept)={V_bi_pn:.2f} V; V_bi-matched W(V_R={V_R_matched:.2f} V)"
    f"={W_matched:.3f} um ({100 * (W_matched / 4.2 - 1):+.1f} %)"
)
print(
    f"abrupt analytic W(40 V): V_bi={V_BI_SCHOTTKY} V -> {w_abrupt_um(40.0, V_BI_SCHOTTKY):.3f} um,"
    f" V_bi={V_bi_pn:.2f} V -> {w_abrupt_um(40.0, V_bi_pn):.3f} um"
)
print(
    f"N_D recovered={N_D_recovered:.3e}, input={cfg.N_D:.3e}, agreement={100*(1-abs(N_D_recovered-cfg.N_D)/cfg.N_D):.1f}%"
)
