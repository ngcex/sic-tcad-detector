"""Generate the C-V comparison figure: ETNA simulation (configured to Bruzzi &
Verroi 2023 device parameters) vs the paper's reported W(V) endpoints and Vbi.

See .planning/VALIDATION_LITERATURE.md for the full writeup and caveats.
Panel 1 (W vs bias) is the genuine independent check: high-bias agreement
(2%) confirms correct absolute C-V scale for an independently-fabricated
device; the low-bias offset is an expected p-n vs Schottky Vbi difference,
not a discrepancy. Panel 2 (1/C^2 slope) is a solver self-consistency
round-trip -- N_D was a simulation INPUT (set to match the paper's device),
not something extracted from Bruzzi's raw data -- so it does NOT independently
validate against the paper; it confirms the DD solver reproduces the analytic
Mott-Schottky relation, which it should by construction.
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
result = run_cv(cfg)

eps0 = 8.854e-14
eps_r = 9.7
W_um = eps_r * eps0 * cfg.area_cm2 / result.y * 1e4
V_reverse = -result.x  # plot as positive reverse bias, matching paper convention

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Panel 1: W(V)
ax = axes[0]
ax.plot(
    V_reverse, W_um, "b-", lw=2, label="ETNA simulation\n(p-n junction, same $N_D$)"
)
ax.plot(
    [0, 40],
    [0.7, 4.2],
    "ro",
    ms=10,
    mfc="none",
    mew=2,
    label="Bruzzi & Verroi 2023\n(Schottky, reported endpoints)",
)
ax.set_xlabel("Reverse bias |V| (V)")
ax.set_ylabel("Depletion width W (µm)")
ax.set_xlim(0, 45)
ax.set_title("Depletion width vs bias")
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3)
ax.annotate(
    "Junction-type offset\n(V$_{bi}$ differs: p-n vs Schottky)",
    xy=(0, 1.117),
    xytext=(8, 2.3),
    fontsize=8,
    color="darkred",
    arrowprops=dict(arrowstyle="->", color="darkred", alpha=0.6),
)
ax.annotate(
    "Converges at high |V|\n(V$_R \\gg V_{bi}$ for both)",
    xy=(40, 4.284),
    xytext=(15, 4.6),
    fontsize=8,
    color="darkgreen",
    arrowprops=dict(arrowstyle="->", color="darkgreen", alpha=0.6),
)

# Panel 2: 1/C^2 vs V (Mott-Schottky) -- the junction-type-independent check
ax = axes[1]
inv_C2 = 1.0 / result.y**2
ax.plot(V_reverse, inv_C2, "b-", lw=2, label="ETNA simulation")
mask = result.x <= -5.0
slope, intercept = np.polyfit(result.x[mask], inv_C2[mask], 1)
q = 1.602e-19
N_D_recovered = -2.0 / (q * eps_r * eps0 * cfg.area_cm2**2 * slope)
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
    "ETNA C-V validation vs Bruzzi & Verroi 2023 (MDPI Materials 16:3643, open access)",
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
print(
    f"W(-40V) sim={W_um[np.argmin(np.abs(result.x-(-40.0)))]:.3f} um, paper endpoint=4.2 um"
)
print(
    f"N_D recovered={N_D_recovered:.3e}, input={cfg.N_D:.3e}, agreement={100*(1-abs(N_D_recovered-cfg.N_D)/cfg.N_D):.1f}%"
)
