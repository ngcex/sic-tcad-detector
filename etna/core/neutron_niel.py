"""Neutron NIEL (displacement damage energy) of Si and 4H-SiC from NJOY HEATR.

Data: ``data/niel/neutron_damage_<iso>_Ed<ed>.csv`` are the pointwise damage-energy
production cross sections (ENDF MF3/MT444, eV barn) that NJOY2016 HEATR computes from
ENDF/B-VIII.0 with the Lindhard/Robinson partition and a sharp displacement threshold
E_d; ``data/niel/neutron_niel_<material>_<tag>.csv`` are the natural-abundance NIEL
(MeV cm^2/g) lethargy-averaged on ``NIEL_GRID_MEV``. Both are written by
``scripts/build_neutron_niel.py`` (needs NJOY + the ENDF files; nothing here calls NJOY).
Provenance, the checks against ASTM E722 and the literature, and the caveats are in
``data/niel/README.md``.

This module only reads the stored files and averages them; the derived constants that
the radiation-damage chain uses are embedded in ``etna.core.niel_tables`` and a test
recomputes them from these files.
"""

from __future__ import annotations

import gzip
from functools import lru_cache
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "niel"

AVOGADRO = 6.02214076e23
# Standard atomic weights (g/mol) used for the per-gram normalisation; Si as in the
# ASTM E722 conversion of niel_tables.NIEL_NEUTRON_1MEV_SI_MEV_CM2_G.
M_SI = 28.0855
M_C = 12.011

NJOY_VERSION = "NJOY2016.78 (conda-forge njoy2016 2016.78)"

# isotope -> (ENDF/B-VIII.0 file, MAT, natural atom fraction, sha256 of the file)
# Atom fractions: IUPAC representative isotopic compositions (CIAAW), as listed by
# NIST "Atomic Weights and Isotopic Compositions".
ENDF_FILES = {
    "Si28": (
        "n-014_Si_028.endf",
        1425,
        0.92223,
        "34031d6a70ca4fda8a39b4cb037d966bd8c2cdc0b2df885382a8be1a5df304ab",
    ),
    "Si29": (
        "n-014_Si_029.endf",
        1428,
        0.04685,
        "f212f322a224d8cd2cc6d164ecba8f598e492d1f912b121043994b6ae096cfe1",
    ),
    "Si30": (
        "n-014_Si_030.endf",
        1431,
        0.03092,
        "56a633f775a581c0ac51031694853f2cff5a00b5ee196d56dbbd930af56a2a04",
    ),
    "C12": (
        "n-006_C_012.endf",
        625,
        0.9893,
        "daa7901e0067f9ac6d9b2392514ba00e71e37d00c8a8ecbbde0620237815f981",
    ),
    "C13": (
        "n-006_C_013.endf",
        628,
        0.0107,
        "abaa8cc34398e6ada0a8c451efd5eef1019a37c7a0a708a4280af71ab77621de",
    ),
}

# Displacement-threshold sets (eV). "Ed21": the SR-NIEL calculator default, the same
# thresholds as the stored proton tables (so the proton/neutron ratios are consistent);
# "Ed35-21": sublattice-specific sensitivity run (Si 35, C 21; as the SR-NIEL 35/21
# proton run). For fast neutrons E_d barely matters (recoils are keV-scale; the sharp
# threshold only removes the few recoils below E_d): see data/niel/README.md.
ED_SETS = {
    "Ed21": {"Si": 21.0, "C": 21.0},
    "Ed35-21": {"Si": 35.0, "C": 21.0},
}

# Lethargy grid of the averaged tables: 20 bins per decade, 10 eV .. 20 MeV.
NIEL_GRID_MEV = np.logspace(
    -5, np.log10(20.0), int(round(20 * (np.log10(20.0) + 5))) + 1
)


def damage_csv_name(iso: str, ed_eV: float) -> str:
    return f"neutron_damage_{iso}_Ed{ed_eV:g}.csv.gz"


def niel_csv_name(material: str, tag: str) -> str:
    return f"neutron_niel_{material}_{tag}.csv"


# ---------------------------------------------------------------- ENDF reading
def _endf_float(s: str) -> float:
    s = s.strip()
    if not s:
        return 0.0
    # ENDF floats drop the 'e': 1.234567+5, -1.2-3
    for i in range(len(s) - 1, 0, -1):
        if s[i] in "+-" and s[i - 1] not in "eE":
            return float(s[:i] + "e" + s[i:])
    return float(s)


def parse_mf3_section(text: str, mat: int, mt: int):
    """(E [eV], sigma) of the TAB1 record of MF3/``mt`` of ``mat`` in an ENDF tape."""
    tag = f"{mat:4d} 3{mt:3d}"
    lines = [ln for ln in text.splitlines() if ln[66:75] == tag]
    if not lines:
        raise ValueError(f"MAT {mat} MF3 MT{mt} not found")
    # line 0: HEAD; line 1: TAB1 control (C1 C2 L1 L2 NR NP); then NR interp pairs
    ctrl = lines[1]
    nr = int(ctrl[44:55])
    np_ = int(ctrl[55:66])
    ri = 2 + (2 * nr + 5) // 6
    vals: list[float] = []
    for ln in lines[ri:]:
        for k in range(6):
            vals.append(_endf_float(ln[11 * k : 11 * k + 11]))
            if len(vals) == 2 * np_:
                break
        if len(vals) == 2 * np_:
            break
    a = np.asarray(vals).reshape(np_, 2)
    return a[:, 0], a[:, 1]


# ---------------------------------------------------------------- NJOY driver
def njoy_input(mat: int, ed_eV: float, title: str) -> str:
    """RECONR (0.1 % reconstruction tolerance, 0 K) then HEATR with MT444 only."""
    return "\n".join(
        [
            "reconr",
            "20 22",
            f"'{title}'/",
            f"{mat} 0/",
            "0.001/",
            "0/",
            "heatr",
            "20 22 23/",
            # matd npk nqa ntemp local iprint ed
            f"{mat} 1 0 0 0 0 {ed_eV:.1f}/",
            "444/",
            "stop",
            "",
        ]
    )


def run_njoy(njoy: Path, endf: Path, mat: int, ed_eV: float, work: Path):
    """Run RECONR + HEATR for one material in ``work``; (E eV, MT444 eV b, input).

    Only used by scripts/build_neutron_niel.py and lab/neutron_niel (needs NJOY)."""
    import subprocess

    work.mkdir(parents=True, exist_ok=True)
    (work / "tape20").write_bytes(endf.read_bytes())
    inp = njoy_input(mat, ed_eV, f"{endf.stem} Ed={ed_eV:g}")
    (work / "input").write_text(inp)
    with open(work / "input") as fin, open(work / "output.log", "w") as fout:
        subprocess.run(
            [str(njoy)],
            stdin=fin,
            stdout=fout,
            stderr=subprocess.STDOUT,
            cwd=work,
            check=True,
        )
    text = (work / "tape23").read_text()
    E, sig = parse_mf3_section(text, mat, 444)
    return E, sig, inp



# ---------------------------------------------------------------- tables
def read_text(path: Path) -> str:
    """Text of a stored table (``.gz`` transparently)."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.decompress(path.read_bytes()).decode("utf-8")
    return path.read_text(encoding="utf-8")


def _load_numeric(path: Path) -> np.ndarray:
    rows = []
    for ln in read_text(path).splitlines():
        if not ln or ln.startswith("#") or ln[0].isalpha():
            continue
        rows.append([float(x) for x in ln.split(",")])
    return np.asarray(rows)


@lru_cache(maxsize=None)
def isotope_damage(iso: str, ed_eV: float, data_dir: str | None = None):
    """Pointwise (E [MeV], damage-energy cross section [MeV cm^2]) of one isotope."""
    d = Path(data_dir) if data_dir else DATA_DIR
    a = _load_numeric(d / damage_csv_name(iso, ed_eV))
    return a[:, 0] * 1e-6, a[:, 1] * 1e-6 * 1e-24


def material_niel_from_isotopes(material: str, eds: dict, data_dir=None):
    """Pointwise natural-abundance NIEL (E [MeV], MeV cm^2/g) on the union grid.

    Si: N_A/M_Si * sum_i f_i sigma_D,i. SiC (1:1): N_A/(M_Si+M_C) * (sigma_D,Si +
    sigma_D,C), each element's sigma_D the abundance-weighted isotope sum. NJOY's
    grids are lin-lin interpolable (RECONR/HEATR output), so the union grid with
    linear interpolation is exact.
    """
    dd = str(data_dir) if data_dir else None
    elements = {"si": ("Si",), "sic": ("Si", "C")}[material]
    isos = [i for i in ENDF_FILES if any(i.startswith(e) for e in elements)]
    curves = []
    for iso in isos:
        el = "Si" if iso.startswith("Si") else "C"
        E, s = isotope_damage(iso, float(eds[el]), dd)
        curves.append((E, s, ENDF_FILES[iso][2]))
    grid = np.unique(np.concatenate([c[0] for c in curves]))
    tot = np.zeros_like(grid)
    for E, s, f in curves:
        tot += f * np.interp(grid, E, s, left=0.0, right=0.0)
    mass = M_SI if material == "si" else M_SI + M_C
    return grid, tot * AVOGADRO / mass


def lethargy_average(E: np.ndarray, y: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Average of the lin-lin curve y(E) over each [edges[i], edges[i+1]] in lethargy.

    Integrates y(E)/E dE exactly on the lin-lin segments (sub-sampled at the curve's
    own points plus the bin edges) and divides by ln(E_hi/E_lo).
    """
    out = np.empty(len(edges) - 1)
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        m = (E > lo) & (E < hi)
        x = np.concatenate([[lo], E[m], [hi]])
        v = np.interp(x, E, y, left=0.0, right=0.0)
        # exact integral of a linear segment divided by E
        x0, x1, v0, v1 = x[:-1], x[1:], v[:-1], v[1:]
        slope = (v1 - v0) / (x1 - x0)
        seg = (v0 - slope * x0) * np.log(x1 / x0) + slope * (x1 - x0)
        out[i] = seg.sum() / np.log(hi / lo)
    return out


def spectrum_average(
    E: np.ndarray, y: np.ndarray, phi, e_min: float, e_max: float, n: int = 20001
) -> float:
    """<y>_phi = int y phi dE / int phi dE on [e_min, e_max] MeV (log-spaced trapezoid,
    plus every tabulated point of y so resonances are resolved)."""
    x = np.unique(
        np.concatenate([np.geomspace(e_min, e_max, n), E[(E > e_min) & (E < e_max)]])
    )
    w = phi(x)
    v = np.interp(x, E, y, left=0.0, right=0.0)
    return float(np.trapezoid(v * w, x) / np.trapezoid(w, x))


def watt_spectrum(a: float = 0.988, b: float = 2.249):
    """Watt fission spectrum chi(E) ~ exp(-E/a) sinh(sqrt(b E)), E in MeV.

    Defaults: the ENDF/B thermal U-235 parameters a = 0.988 MeV, b = 2.249 /MeV.
    """
    return lambda E: np.exp(-E / a) * np.sinh(np.sqrt(b * E))


@lru_cache(maxsize=None)
def neutron_niel_table(material: str, tag: str = "Ed21") -> tuple:
    """Lethargy-averaged neutron NIEL of ``material`` ("si"/"sic") as
    ((E_centre MeV, NIEL MeV cm^2/g), ...) from ``data/niel/neutron_niel_*.csv``."""
    a = _load_numeric(DATA_DIR / niel_csv_name(material, tag))
    return tuple((float(e), float(v)) for e, v in zip(a[:, 0], a[:, 3]))
