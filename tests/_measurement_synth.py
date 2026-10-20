"""SYNTHETIC measurement files for the etna.measurements tests.

Every file written here is generated in-test from the ETNA model itself (or from an
analytic Arrhenius law) and lands in pytest's tmp_path: none of it is measurement data
and none of it is committed. Each file carries a '# notes: SYNTHETIC ...' line.
"""

from __future__ import annotations

import numpy as np

SYNTH_NOTE = (
    "# notes: SYNTHETIC test data generated from the ETNA model (not a measurement)"
)


def g(v) -> str:
    """Full-precision decimal of a (numpy) float."""
    return f"{float(v):.17g}"


def cv_csv(V, C_F, area_cm2=0.04, unit="pF", meta=(), sigma_rel=None, convention=None):
    """C-V CSV text: absolute C in ``unit`` with ``# area_cm2`` metadata."""
    f = {"F": 1.0, "nF": 1e9, "pF": 1e12}[unit]
    lines = [SYNTH_NOTE]
    if area_cm2 is not None:
        lines.append(f"# area_cm2: {g(area_cm2)}")
    if convention:
        lines.append(f"# bias_convention: {convention}")
    lines += list(meta)
    sign = -1.0 if convention == "reverse_positive" else 1.0
    head = f"V_V,C_{unit}" + (f",sigma_C_{unit}" if sigma_rel is not None else "")
    lines.append(head)
    for v, c in zip(V, C_F):
        row = f"{g(sign * v)},{g(c * f)}"
        if sigma_rel is not None:
            row += f",{g(c * f * sigma_rel)}"
        lines.append(row)
    return "\n".join(lines) + "\n"


def sims_csv(depth_um, N, Al=None, meta=(), unit="um"):
    f = {"um": 1.0, "nm": 1e3}[unit]
    lines = [SYNTH_NOTE, "# technique: SIMS"] + list(meta)
    lines.append(f"depth_{unit},N_cm3" + (",Al_cm3" if Al is not None else ""))
    for i, d in enumerate(depth_um):
        row = f"{g(d * f)},{g(N[i])}"
        if Al is not None:
            row += f",{g(Al[i])}"
        lines.append(row)
    return "\n".join(lines) + "\n"


def it_csv(T, I_A, V=None, sigma_A=None, meta=(), unit="pA"):
    f = {"A": 1.0, "nA": 1e9, "pA": 1e12, "fA": 1e15}[unit]
    lines = [SYNTH_NOTE] + list(meta)
    head = ([] if V is None else ["V_V"]) + ["T_K", f"I_{unit}"]
    if sigma_A is not None:
        head.append(f"sigma_I_{unit}")
    lines.append(",".join(head))
    for i in range(len(T)):
        row = ([] if V is None else [g(V[i])]) + [g(T[i]), g(I_A[i] * f)]
        if sigma_A is not None:
            row.append(g(sigma_A[i] * f))
        lines.append(",".join(row))
    return "\n".join(lines) + "\n"


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def rng(seed=12345):
    return np.random.default_rng(seed)
