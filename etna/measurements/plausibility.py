"""Physical-plausibility and unit checks for the measurement loaders.

The schema checks of ``_csv.parse_table`` and the loaders catch malformed files; these
catch files that are well formed but in the wrong unit or convention (F in a ``C_pF``
column, nm in ``depth_um``, m^-3 in ``N_cm3``, degC in ``T_K``, A vs pA). Two levels:

- **errors**: values no detector of this kind can produce whatever the device (W below
  1 nm or above 1 cm, N above the atomic density of SiC, a dark current above 0.1 A,
  ...). They are appended to the loader's issue list, so the file is refused with every
  other problem, each message naming the column and the most likely unit mix-up.
- **warnings**: unusual but possible values (W outside 0.05-500 um, |V| > 1.5 kV, T
  outside 77-700 K, C rising with reverse bias, ...). Emitted as
  ``MeasurementPlausibilityWarning`` and returned by ``plausibility_report``.

The bounds are deliberately wide (orders of magnitude, not device tolerances): they
are sanity limits for 4H-SiC pad diodes, not model validity ranges.
"""

from __future__ import annotations

import math
import warnings

import numpy as np

# 4H-SiC atomic density: 3.21 g/cm^3 * N_A / 40.097 g/mol * 2 atoms = 9.64e22 cm^-3
SIC_ATOM_DENSITY_CM3 = 9.64e22

W_ERROR_UM = (1e-3, 1e4)  # depletion width implied by C/A
W_WARN_UM = (0.05, 500.0)
AREA_ERROR_CM2 = (1e-6, 100.0)
V_ERROR_V = 5000.0
V_WARN_V = 1500.0
FREQ_WARN_HZ = (10.0, 1e8)
T_ERROR_K = (4.0, 2000.0)
T_WARN_K = (77.0, 700.0)
I_ERROR_A = 0.1
I_WARN_A = (1e-17, 1e-4)
DEPTH_ERROR_UM = 5000.0
N_WARN_CM3 = 1e12


class MeasurementPlausibilityWarning(UserWarning):
    """A measurement file is well formed but holds unusual values (unit? convention?)."""


def _decades_to(value: float, lo: float, hi: float) -> int:
    """Power of ten k such that value * 10**k is inside [lo, hi] (closest to the
    log-centre); 0 if already inside."""
    if lo <= value <= hi or not value > 0:
        return 0
    centre = 0.5 * (math.log10(lo) + math.log10(hi))
    return int(round(centre - math.log10(value)))


def _unit_hint(k: int, what: str) -> str:
    known = {
        12: "values look like F in a pF column (x1e12)",
        -12: "values look like pF in an F column (x1e-12)",
        9: "values look like F in an nF column (x1e9)",
        -9: "values look like nF in an F column",
        3: "values look like nF in a pF column, or um vs nm",
        -3: "values look like pF in an nF column, or nm in an um column",
        6: "a factor 1e6 (cm^-3 vs m^-3, or A vs uA)",
        -6: "a factor 1e-6 (m^-3 in a cm^-3 column?)",
        2: "a factor 100 (mm^2 vs cm^2?)",
        -2: "a factor 100 (mm^2 vs cm^2?)",
    }
    hint = known.get(k, f"off by about 10^{k}")
    return f"{what}: {hint}"


def check_cv(V, C_per_area, area_cm2, freq_hz, T_K, cname, eps_F_cm):
    """(errors, warnings) for a C-V file; C_per_area in F/cm^2, V reverse-negative."""
    err, warn = [], []
    CA = np.asarray(C_per_area, float)
    ok = CA > 0
    if ok.any():
        W_um = eps_F_cm / CA[ok] * 1e4
        med = float(np.median(W_um))
        if med < W_ERROR_UM[0] or med > W_ERROR_UM[1]:
            k = _decades_to(med, 0.5, 50.0)
            err.append(
                f"column '{cname}': the per-area capacitance gives a depletion width of "
                f"{med:.3g} um (median), impossible for a diode; "
                + _unit_hint(-k, "check the capacitance unit / the area")
            )
        elif med < W_WARN_UM[0] or med > W_WARN_UM[1]:
            warn.append(
                f"column '{cname}': median depletion width {med:.3g} um (eps/(C/A)) is "
                f"outside {W_WARN_UM[0]:g}-{W_WARN_UM[1]:g} um; check the unit and area"
            )
    if area_cm2 is not None and not (
        AREA_ERROR_CM2[0] <= area_cm2 <= AREA_ERROR_CM2[1]
    ):
        k = _decades_to(area_cm2, 1e-3, 1.0)
        err.append(
            f"area {area_cm2:g} cm^2 is outside {AREA_ERROR_CM2[0]:g}-"
            f"{AREA_ERROR_CM2[1]:g} cm^2; "
            + _unit_hint(k, "check area_cm2 vs area_mm2")
        )
    Vmax = float(np.max(np.abs(V))) if len(V) else 0.0
    if Vmax > V_ERROR_V:
        err.append(f"|V| up to {Vmax:g} V (> {V_ERROR_V:g} V): mV in a V_V column?")
    elif Vmax > V_WARN_V:
        warn.append(f"|V| up to {Vmax:g} V (> {V_WARN_V:g} V)")
    if freq_hz is not None and not (FREQ_WARN_HZ[0] <= freq_hz <= FREQ_WARN_HZ[1]):
        warn.append(f"test frequency {freq_hz:g} Hz is unusual (Hz vs kHz?)")
    if T_K is not None and not (T_WARN_K[0] <= T_K <= T_WARN_K[1]):
        warn.append(f"temperature {T_K:g} K is unusual (degC in a K field?)")
    # reverse branch ordered from 0 towards the deepest bias: C should not rise
    V = np.asarray(V, float)
    rev = V <= 0
    if rev.sum() >= 3:
        o = np.argsort(-V[rev])
        c = CA[rev][o]
        rise = np.where(c[1:] > c[:-1] * 1.05)[0]
        if len(rise):
            warn.append(
                f"C rises by > 5 % with reverse bias at {len(rise)} step(s) (first at "
                f"V = {V[rev][o][rise[0] + 1]:g} V): bias sign/convention, a forward "
                "branch, or series resistance / frequency effects?"
            )
    return err, warn


def check_sims(depth_um, columns, names):
    """(errors, warnings) for a SIMS/SRP file; ``columns``: {"N": arr, "Al": arr}."""
    err, warn = [], []
    d = np.asarray(depth_um, float)
    dname = names.get("depth", "depth")
    if len(d) and d.max() > DEPTH_ERROR_UM:
        err.append(
            f"column '{dname}': depth up to {d.max():g} um (> {DEPTH_ERROR_UM:g} um): "
            "nm or Angstrom in a depth_um column? (use depth_nm)"
        )
    elif len(d) and d.max() < 0.2:
        warn.append(
            f"column '{dname}': the profile ends at {d.max():g} um; um vs nm, or a "
            "surface-only profile?"
        )
    for q, v in columns.items():
        if v is None:
            continue
        v = np.asarray(v, float)
        pos = v[v > 0]
        if not len(pos):
            continue
        cname = names.get(q, q)
        if pos.max() > SIC_ATOM_DENSITY_CM3:
            err.append(
                f"column '{cname}': {pos.max():.3g} cm^-3 exceeds the atomic density of "
                f"SiC ({SIC_ATOM_DENSITY_CM3:.3g}); m^-3 in a cm^-3 column (/1e6)?"
            )
        elif pos.max() <= 30.0:
            err.append(
                f"column '{cname}': every value <= 30: these look like log10(N); give "
                "the concentration itself in cm^-3"
            )
        elif pos.min() < N_WARN_CM3:
            warn.append(
                f"column '{cname}': values down to {pos.min():.3g} cm^-3, below any "
                "SIMS/SRP detection limit in SiC; unit? (0 = below detection)"
            )
    return err, warn


def check_it(T_K, I_A, tname, iname, t_in_kelvin_column: bool):
    """(errors, warnings) for an I(T) file (T in K after conversion, I in A)."""
    err, warn = [], []
    T = np.asarray(T_K, float)
    I = np.abs(np.asarray(I_A, float))
    if len(T) and (T.min() < T_ERROR_K[0] or T.max() > T_ERROR_K[1]):
        err.append(
            f"column '{tname}': T from {T.min():g} to {T.max():g} K is outside "
            f"{T_ERROR_K[0]:g}-{T_ERROR_K[1]:g} K"
        )
    elif len(T) and (T.min() < T_WARN_K[0] or T.max() > T_WARN_K[1]):
        warn.append(
            f"column '{tname}': T from {T.min():g} to {T.max():g} K"
            + (
                "; values look like degC: use a T_C column"
                if t_in_kelvin_column and T.max() < 250
                else " is unusual for a SiC dark-current I(T)"
            )
        )
    if len(I) and I.max() > I_ERROR_A:
        err.append(
            f"column '{iname}': |I| up to {I.max():.3g} A (> {I_ERROR_A:g} A) is not a "
            "dark current; check the unit (pA in an I_A column?)"
        )
    elif len(I) and I.max() > I_WARN_A[1]:
        warn.append(
            f"column '{iname}': |I| up to {I.max():.3g} A (> {I_WARN_A[1]:g} A), large "
            "for a SiC pad-diode dark current; unit? breakdown?"
        )
    if len(I) and 0 < I.min() < I_WARN_A[0]:
        warn.append(
            f"column '{iname}': |I| down to {I.min():.3g} A, below a 10 aA floor; unit?"
        )
    return err, warn


def emit(source: str, msgs) -> list:
    """Warn once per message (MeasurementPlausibilityWarning); return the messages."""
    msgs = list(msgs)
    for m in msgs:
        warnings.warn(f"{source}: {m}", MeasurementPlausibilityWarning, stacklevel=3)
    return msgs


def plausibility_report(data) -> list:
    """Warnings (strings) of a loaded ``CVData`` / ``DepthProfile`` / ``ITData``, without
    emitting them (errors cannot occur: the loader refused such files)."""
    from etna.measurements.cv import CVData, _eps
    from etna.measurements.iv_t import ITData
    from etna.measurements.sims import DepthProfile

    if isinstance(data, CVData):
        return check_cv(
            data.V,
            data.C_per_area,
            data.area_cm2,
            data.frequency_Hz,
            data.T_K,
            data.input_columns.get("C", data.input_columns.get("CA", "C")),
            _eps(),
        )[1]
    if isinstance(data, DepthProfile):
        return check_sims(
            data.depth_um, {"N": data.N, "Al": data.Al}, data.input_columns
        )[1]
    if isinstance(data, ITData):
        cols = data.input_columns
        return check_it(
            data.T_K,
            data.I_A,
            cols.get("T", cols.get("TC", "T")),
            cols.get("I", "I"),
            "T" in cols,
        )[1]
    raise TypeError(f"not a measurement data object: {type(data).__name__}")
