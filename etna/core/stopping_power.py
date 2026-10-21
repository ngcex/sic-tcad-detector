"""Tabulated NIST PSTAR / ASTAR stopping powers and the quantities derived from them.

Data (downloaded 2026-10-20 from the NIST STAR text interface, see
``data/nist_star/README.md`` for URLs, request parameters and SHA-256 of the raw
responses):

- PSTAR (protons) and ASTAR (alpha particles) for silicon (014), amorphous carbon
  (006, rho 2.0 g/cm^3, I = 81 eV), graphite (906, rho 1.7, I = 78 eV), liquid
  water (276, I = 75 eV) and ICRU striated muscle (202, rho 1.04, I = 74.7 eV).
- 1 keV to 10 GeV (PSTAR) / 1 keV to 1 GeV (ASTAR), 7 columns per row
  (T, S_e, S_n, S_tot in MeV cm^2/g, CSDA and projected range in g/cm^2, detour).

NIST has no SiC table. 4H-SiC is composed here by **Bragg additivity** of the
elemental mass stopping powers, ``(S/rho)_SiC = w_Si (S/rho)_Si + w_C (S/rho)_C``
(w_Si = 0.7004, w_C = 0.2996). Caveats, stated rather than corrected:

- Bragg additivity ignores chemical-binding and phase effects (error sizes below are
  an UNVERIFIED literature estimate, not checked in this repo). For protons above
  ~1 MeV the error is usually a few percent or less; near the stopping-power maximum
  (protons < ~0.3 MeV, alphas < ~2 MeV) it can reach ~10%.
- The carbon component is NIST "amorphous carbon" by default. Graphite instead
  changes the composed SiC stopping by <= 1% above 1 MeV (alphas: ~3% at 0.5 MeV).
- The SiC CSDA range is computed by integrating 1/S_tot of the composed table (the
  same integration reproduces NIST's own Si CSDA range to < 0.1%, tested).
"""

from __future__ import annotations

import functools
import re
from pathlib import Path

import numpy as np

from etna._paths import data_path

DATA_DIR = data_path("nist_star")  # repo data/ in a checkout, etna/_data in a wheel

# 4H-SiC mass fractions (Si 28.0855, C 12.011 g/mol; 1:1 stoichiometry).
SIC_MASS_FRACTION_SI = 0.7004
SIC_MASS_FRACTION_C = 0.2996

# Densities (g/cm^3) of the NIST target materials (from the NIST composition pages,
# saved in data/nist_star/raw/compos_*.html). SiC from SiC4H_Parameters.
MATERIAL_DENSITY = {
    "water": 1.0,
    "muscle_icru": 1.04,
    "si": 2.33,
}

# NIST material numbers (matno) -> file stem used in data/nist_star/
NIST_MATERIALS = {
    "si": "014",
    "c_amorphous": "006",
    "c_graphite": "906",
    "water": "276",
    "muscle_icru": "202",
}
PROGRAMS = {"proton": "pstar", "alpha": "astar"}

COLUMNS = (
    "energy_MeV",
    "S_electronic_MeV_cm2_g",
    "S_nuclear_MeV_cm2_g",
    "S_total_MeV_cm2_g",
    "csda_range_g_cm2",
    "projected_range_g_cm2",
    "detour_factor",
)

_ROW = re.compile(r"^\s*(\d\.\d+E[+-]\d+)((?:\s+[-+0-9.E]+){6})\s*$")


def parse_nist_star_html(text):
    """Parse a NIST PSTAR/ASTAR text-interface response (HTML with a <pre> table).

    Returns ``(material_name, table)`` with ``table`` an (N, 7) float array in the
    column order of ``COLUMNS``.
    """
    start = text.index("<pre>")
    stop = text.index("</pre>", start)
    body = text[start + 5 : stop].replace("<br>", "\n")
    lines = body.splitlines()
    name = ""
    for i, line in enumerate(lines):
        if "STOPPING POWER AND RANGE" in line:
            # the material name is the first non-blank line after the title
            for nxt in lines[i + 1 :]:
                if nxt.strip():
                    name = nxt.strip()
                    break
            break
    rows = []
    for line in lines:
        m = _ROW.match(line)
        if m:
            rows.append([float(m.group(1))] + [float(v) for v in m.group(2).split()])
    if not rows:
        raise ValueError("no NIST STAR table rows found in the response")
    return name, np.asarray(rows, dtype=float)


def _csv_path(particle, material):
    return DATA_DIR / f"{PROGRAMS[particle]}_{material}.csv"


@functools.lru_cache(maxsize=None)
def load_star(particle, material):
    """Load a derived NIST table as a dict of read-only arrays keyed by ``COLUMNS``.

    ``particle`` in {"proton", "alpha"}; ``material`` in ``NIST_MATERIALS``.
    """
    if particle not in PROGRAMS:
        raise ValueError(
            f"particle must be one of {sorted(PROGRAMS)} (got {particle!r})"
        )
    if material not in NIST_MATERIALS:
        raise ValueError(
            f"material must be one of {sorted(NIST_MATERIALS)} (got {material!r})"
        )
    path = _csv_path(particle, material)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing. Regenerate it from data/nist_star/raw with "
            "`uv run python scripts/build_nist_star_tables.py`."
        )
    arr = _read_csv(path)
    out = {}
    for j, col in enumerate(COLUMNS):
        a = np.array(arr[:, j], dtype=float)
        a.setflags(write=False)
        out[col] = a
    return out


def _read_csv(path):
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("energy_MeV"):
            continue
        rows.append([float(v) for v in line.split(",")])
    return np.asarray(rows, dtype=float)


def _loglog(x, xp, fp):
    x = np.asarray(x, dtype=float)
    return np.exp(np.interp(np.log(x), np.log(xp), np.log(fp)))


def _column(component):
    cols = {
        "electronic": "S_electronic_MeV_cm2_g",
        "nuclear": "S_nuclear_MeV_cm2_g",
        "total": "S_total_MeV_cm2_g",
    }
    if component not in cols:
        raise ValueError(f"component must be one of {sorted(cols)} (got {component!r})")
    return cols[component]


def mass_stopping(particle, material, energy_MeV=None, component="electronic"):
    """Mass stopping power (MeV cm^2/g) of a NIST material, log-log interpolated.

    ``material`` may be any key of ``NIST_MATERIALS`` or ``"sic"`` (Bragg additivity,
    see module docstring). Without ``energy_MeV`` the NIST grid is returned as
    ``(energy_MeV, S)``.
    """
    if material == "sic":
        return sic_mass_stopping(particle, energy_MeV, component=component)
    t = load_star(particle, material)
    e = t["energy_MeV"]
    s = t[_column(component)]
    if energy_MeV is None:
        return e, s
    _check_energy(energy_MeV, e, particle)
    return _loglog(energy_MeV, e, s)


def _check_energy(energy_MeV, grid, particle):
    E = np.asarray(energy_MeV, dtype=float)
    if np.any(~np.isfinite(E)) or np.any(E < grid[0]) or np.any(E > grid[-1]):
        raise ValueError(
            f"{particle} energy outside the NIST table [{grid[0]:g}, {grid[-1]:g}] MeV"
        )


def sic_mass_stopping(
    particle, energy_MeV=None, component="electronic", carbon="c_amorphous"
):
    """4H-SiC mass stopping power (MeV cm^2/g) by Bragg additivity of NIST Si + C."""
    if carbon not in ("c_amorphous", "c_graphite"):
        raise ValueError("carbon must be 'c_amorphous' or 'c_graphite'")
    si = load_star(particle, "si")
    c = load_star(particle, carbon)
    col = _column(component)
    grid = si["energy_MeV"]
    s = SIC_MASS_FRACTION_SI * si[col] + SIC_MASS_FRACTION_C * _loglog(
        grid, c["energy_MeV"], c[col]
    )
    if energy_MeV is None:
        return grid.copy(), s
    _check_energy(energy_MeV, grid, particle)
    return _loglog(energy_MeV, grid, s)


def _rho_sic():
    from etna.core.sic_material import SiC4H_Parameters

    return float(SiC4H_Parameters().rho)


def density(material):
    """Density (g/cm^3) used for linear stopping power / ranges."""
    if material == "sic":
        return _rho_sic()
    if material in MATERIAL_DENSITY:
        return MATERIAL_DENSITY[material]
    raise ValueError(f"no density for material {material!r}")


def linear_stopping(particle, material, energy_MeV, component="electronic"):
    """Linear stopping power (MeV/cm) = mass stopping x density."""
    return mass_stopping(particle, material, energy_MeV, component) * density(material)


@functools.lru_cache(maxsize=None)
def _csda_table(particle, material, n=200001):
    """Fine log grid (E, CSDA range g/cm^2) by integrating 1/S_tot (trapezoid)."""
    if material == "sic":
        e_tab, s_tab = sic_mass_stopping(particle, component="total")
        si = load_star(particle, "si")
        c = load_star(particle, "c_amorphous")
        # range at the first grid point (1 keV): mass-fraction harmonic mix of the
        # elemental NIST CSDA ranges (it is ~1e-5 g/cm^2, i.e. ~0.03 um in SiC).
        r0 = 1.0 / (
            SIC_MASS_FRACTION_SI / si["csda_range_g_cm2"][0]
            + SIC_MASS_FRACTION_C / c["csda_range_g_cm2"][0]
        )
    else:
        t = load_star(particle, material)
        e_tab, s_tab = t["energy_MeV"], t["S_total_MeV_cm2_g"]
        r0 = t["csda_range_g_cm2"][0]
    e = np.geomspace(e_tab[0], e_tab[-1], n)
    inv = 1.0 / _loglog(e, e_tab, s_tab)
    r = r0 + np.concatenate([[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(e))])
    e.setflags(write=False)
    r.setflags(write=False)
    return e, r


def csda_range_cm(particle, material, energy_MeV):
    """CSDA range (cm): NIST's own column (log-log interpolated) for NIST materials,
    integration of 1/S_tot of the Bragg-composed table for ``"sic"``."""
    if material != "sic":
        t = load_star(particle, material)
        _check_energy(energy_MeV, t["energy_MeV"], particle)
        return _loglog(energy_MeV, t["energy_MeV"], t["csda_range_g_cm2"]) / density(
            material
        )
    e, r = _csda_table(particle, material)
    _check_energy(energy_MeV, e, particle)
    return np.interp(energy_MeV, e, r) / density(material)


def residual_energy_MeV(particle, material, E0_MeV, depth_cm):
    """Kinetic energy (MeV) left after ``depth_cm`` of CSDA path; 0 beyond the range."""
    e, r = _csda_table(particle, material)
    _check_energy(E0_MeV, e, particle)
    rho = density(material)
    r0 = np.interp(E0_MeV, e, r)
    left = r0 - np.asarray(depth_cm, dtype=float) * rho
    out = np.interp(np.clip(left, r[0], None), r, e)
    return np.where(left > r[0], out, 0.0)


def kappa_table(
    particle="proton",
    tissue="water",
    ratio="linear",
    carbon="c_amorphous",
    component="electronic",
):
    """SiC -> tissue stopping-power ratio vs KINETIC energy on the NIST grid.

    ratio="linear" (default): kappa = L_tissue / L_SiC = (S/rho)_t rho_t /
    ((S/rho)_SiC rho_SiC), the LET ratio. For a track crossing the SV, lineal energy
    y ~ LET, so y_tissue = kappa_linear * y_SiC (the convention of the Si
    microdosimetry literature, e.g. ~0.55-0.58 for Si -> water/muscle in proton and
    carbon beams; for SiC it is ~0.38-0.45).
    ratio="mass": (S/rho)_t / (S/rho)_SiC (~1.2-1.5); appropriate for converting
    absorbed DOSE between media (Bragg-Gray), not lineal energy.
    """
    if tissue not in ("water", "muscle_icru"):
        raise ValueError("tissue must be 'water' or 'muscle_icru'")
    if ratio not in ("linear", "mass"):
        raise ValueError("ratio must be 'linear' or 'mass'")
    e_t, s_t = mass_stopping(particle, tissue, component=component)
    s_sic = sic_mass_stopping(particle, e_t, component=component, carbon=carbon)
    k = s_t / s_sic
    if ratio == "linear":
        k = k * density(tissue) / _rho_sic()
    return {
        "energy_MeV": np.array(e_t, dtype=float),
        "kappa": k,
        "particle": particle,
        "tissue": tissue,
        "ratio": ratio,
        "source": "bragg",
    }


# ---------------------------------------------------------------------------
# Alpha Bragg curve in SiC
# ---------------------------------------------------------------------------

AM241_ALPHA_E0_MEV = 5.486

# Physical constants for Bohr straggling
_E2_MEV_CM = 1.439964e-13  # e^2 / (4 pi eps0), MeV cm
_N_A = 6.02214076e23


def sic_electron_density_cm3():
    """Electron density of 4H-SiC (cm^-3) from rho and the mass fractions."""
    return (
        _rho_sic()
        * _N_A
        * (SIC_MASS_FRACTION_SI * 14 / 28.0855 + SIC_MASS_FRACTION_C * 6 / 12.011)
    )


def bohr_range_straggling_cm(particle, E0_MeV, n=20001):
    """Bohr range straggling sigma_R (cm) in SiC: sigma_R^2 = int Omega_B^2' / L^3 dE.

    dOmega_B^2/dx = 4 pi z^2 e^4 N_e (Bohr, non-relativistic, electrons free). The
    Bohr formula overestimates straggling at low velocity, so this is an upper-ish
    estimate (computed, not fitted).
    """
    z = {"proton": 1, "alpha": 2}[particle]
    d_omega2 = 4.0 * np.pi * z**2 * _E2_MEV_CM**2 * sic_electron_density_cm3()
    e_tab, _ = sic_mass_stopping(particle, component="total")
    e = np.geomspace(e_tab[0], E0_MeV, n)
    L = sic_mass_stopping(particle, e, component="total") * _rho_sic()
    return float(np.sqrt(np.trapezoid(d_omega2 / L**3, e)))


@functools.lru_cache(maxsize=8)
def alpha_bragg_curve_sic(E0_MeV=AM241_ALPHA_E0_MEV, n=4001, straggle=True):
    """Ionization dE/dx (MeV/cm, electronic) vs CSDA path depth (cm) in SiC.

    Built from ASTAR Si + C (Bragg additivity): depth s(E) = R(E0) - R(E), local
    S_e(E(s)); optionally convolved with a Gaussian of the computed Bohr range
    straggling. Returns ``(s_cm, dEdx_MeV_cm, csda_range_cm, sigma_cm)``; the curve
    is on a uniform grid 0..1.3 R.
    """
    R = float(csda_range_cm("alpha", "sic", E0_MeV))
    sigma = bohr_range_straggling_cm("alpha", E0_MeV) if straggle else 0.0
    s = np.linspace(0.0, 1.3 * R, n)
    E = residual_energy_MeV("alpha", "sic", E0_MeV, s)
    e_min = load_star("alpha", "si")["energy_MeV"][0]
    f = np.where(
        E > e_min,
        linear_stopping("alpha", "sic", np.clip(E, e_min, None), "electronic"),
        0.0,
    )
    if sigma > 0:
        ds = s[1] - s[0]
        half = max(1, int(np.ceil(5 * sigma / ds)))
        k = np.arange(-half, half + 1) * ds
        ker = np.exp(-0.5 * (k / sigma) ** 2)
        ker /= ker.sum()
        padded = np.concatenate([np.full(half, f[0]), f, np.zeros(half)])
        f = np.convolve(padded, ker, mode="valid")
    s.setflags(write=False)
    f.setflags(write=False)
    return s, f, R, sigma


def alpha_bragg_summary(E0_MeV=AM241_ALPHA_E0_MEV):
    """Data-derived alpha numbers in SiC: range, straggling, peak/entrance ratio."""
    s, f, R, sigma = alpha_bragg_curve_sic(E0_MeV)
    s0, f0, _, _ = alpha_bragg_curve_sic(E0_MeV, straggle=False)
    i_pk = int(np.argmax(f))
    e_grid = np.geomspace(0.01, E0_MeV, 20001)
    s_e = sic_mass_stopping("alpha", e_grid, component="electronic")
    return {
        "E0_MeV": E0_MeV,
        "csda_range_um": R * 1e4,
        "bohr_straggling_um": sigma * 1e4,
        "peak_to_entrance": float(f[i_pk] / f[0]),
        "peak_to_entrance_unstraggled": float(np.max(f0) / f0[0]),
        "peak_depth_um": float(s[i_pk] * 1e4),
        "E_peak_MeV": float(e_grid[np.argmax(s_e)]),
    }
