"""NIST ESTAR electron stopping powers and CSDA ranges (delta-ray transport input).

Data: ``data/nist_estar/raw/`` holds the unmodified NIST ESTAR text-interface
responses (POST https://physics.nist.gov/cgi-bin/Star/e_table-t.pl, matno=014 silicon,
006 amorphous carbon, 276 liquid water; retrieved 2026-10-21, see the README there):

- ``estar_lowE_<matno>.html``: ``Energies`` = 1, 1.5, ..., 8 keV plus ``ShowDefault=on``;
  NIST then prints T, collision, radiative, total stopping power and the density-effect
  parameter from 1 keV to 1 GeV (no range column).
- ``estar_default_<matno>.html``: the default grid (10 keV .. 1 GeV) WITH NIST's CSDA
  range column; used only to check the range integration below (test: < 1 % above
  20 keV for Si, C and water).

SiC is composed by Bragg additivity of NIST Si + amorphous C (mass fractions of
``etna.core.stopping_power``), as for the ion tables. The CSDA range is the integral of
1/S_tot from 1 keV, started with the power-law range R(1 keV) = T / (S (1 + p)), where
p = -dlnS/dlnT between 1 and 1.5 keV (S ~ T^-p below 1 keV is an assumption; R(1 keV)
is ~0.3-0.4 um in SiC, small against the ranges that matter for escape).
"""

from __future__ import annotations

import functools
import re
from pathlib import Path

import numpy as np

from etna.core.stopping_power import SIC_MASS_FRACTION_C, SIC_MASS_FRACTION_SI

ESTAR_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "nist_estar"
ESTAR_MATNO = {"si": "014", "c_amorphous": "006", "water": "276"}
_ROW = re.compile(r"^\s*(\d\.\d+E[+-]\d+)((?:\s+[-+0-9.E]+){4,6})\s*$")


def parse_estar_html(text):
    """Rows (N, k) of a NIST ESTAR text response (``<pre>`` table, ``<br>`` lines)."""
    start = text.index("<pre>")
    stop = text.index("</pre>", start)
    body = text[start + 5 : stop].replace("<br>", "\n")
    rows = []
    for line in body.splitlines():
        m = _ROW.match(line)
        if m:
            rows.append([float(m.group(1))] + [float(v) for v in m.group(2).split()])
    if not rows:
        raise ValueError("no NIST ESTAR rows found")
    return np.asarray(rows, dtype=float)


@functools.lru_cache(maxsize=None)
def estar_total_stopping(material):
    """(T_MeV, S_total MeV cm^2/g) from 1 keV, for "si", "c_amorphous", "water", "sic"."""
    if material == "sic":
        e_si, s_si = estar_total_stopping("si")
        e_c, s_c = estar_total_stopping("c_amorphous")
        s_c_on = np.exp(np.interp(np.log(e_si), np.log(e_c), np.log(s_c)))
        return e_si, SIC_MASS_FRACTION_SI * s_si + SIC_MASS_FRACTION_C * s_c_on
    if material not in ESTAR_MATNO:
        raise ValueError(f"no ESTAR table for {material!r}")
    path = ESTAR_DIR / "raw" / f"estar_lowE_{ESTAR_MATNO[material]}.html"
    t = parse_estar_html(path.read_text(encoding="latin-1"))
    # columns: T, collision, radiative, total, density-effect delta
    e, s = t[:, 0], t[:, 3]
    e.setflags(write=False)
    s.setflags(write=False)
    return e, s


def estar_default_csda(material):
    """NIST's own CSDA range column (T_MeV, R g/cm^2) on the default grid (tests)."""
    path = ESTAR_DIR / "raw" / f"estar_default_{ESTAR_MATNO[material]}.html"
    t = parse_estar_html(path.read_text(encoding="latin-1"))
    return t[:, 0], t[:, 4]


@functools.lru_cache(maxsize=None)
def electron_csda_table(material, n=4001):
    """Fine log grid (T_MeV, CSDA range g/cm^2) from 1 keV by integrating 1/S_tot."""
    e_tab, s_tab = estar_total_stopping(material)
    p = -np.log(s_tab[1] / s_tab[0]) / np.log(e_tab[1] / e_tab[0])
    r0 = e_tab[0] / (s_tab[0] * (1.0 + p))
    e = np.geomspace(e_tab[0], e_tab[-1], n)
    inv = 1.0 / np.exp(np.interp(np.log(e), np.log(e_tab), np.log(s_tab)))
    r = r0 + np.concatenate([[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(e))])
    e.setflags(write=False)
    r.setflags(write=False)
    return e, r


def electron_csda_range_cm(material, T_MeV, density_g_cm3):
    """Electron CSDA range (cm) at kinetic energy ``T_MeV`` (log-log interpolation)."""
    e, r = electron_csda_table(material)
    T = np.asarray(T_MeV, dtype=float)
    out = np.exp(np.interp(np.log(T), np.log(e), np.log(r))) / density_g_cm3
    return out
