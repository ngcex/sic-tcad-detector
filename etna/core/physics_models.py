"""Selectable physical-model library for 4H-SiC (opt-in; the defaults are unchanged).

One frozen, hashable ``PhysicsModels`` selects, by registry name or by instance:

* ``low_field_mobility``: doping/temperature Caughey-Thomas (``LOW_FIELD_MOBILITY``);
* ``high_field_mobility``: velocity saturation, Caughey-Thomas/Canali form
  mu = mu0 / (1 + (mu0 E / v_sat)^beta)^(1/beta) (``HIGH_FIELD_MOBILITY``);
* ``anisotropy``: mobility ratio parallel / perpendicular to the c-axis for the
  vertical (c-axis) transport of the 1D diode (``ANISOTROPY``);
* ``impact_ionization``: avalanche coefficients
  alpha = a (1 + c (T - 300)) exp(-[b (1 + d (T - 300)) / E]^m) (Chynoweth m = 1,
  Okuto-Crowell-type m != 1; ``IMPACT_IONIZATION``);
* ``bandgap_narrowing``: Slotboom-type Delta E_g(N) (``BANDGAP_NARROWING``;
  informational: not coupled into the solvers yet, see docs/USAGE.md).

``PhysicsModels()`` (all None / "etna_default") is the existing ETNA model, bit for bit
(tests/test_physics_models.py). Every parameter set records its provenance
(``Provenance``): citation, DOI/URL actually opened, and whether the numbers were read in
the primary source ("verified-primary"), in a secondary compilation that restates it
("verified-secondary"; mostly J. Burin et al., "TCAD parameters for 4H-SiC: a review",
arXiv:2410.06798, and the TU Wien thesis of T. Ayalew, 2004,
https://www.iue.tuwien.ac.at/phd/ayalew/), or are an ETNA assumption ("assumption").

Direction convention: the 1D diode is a vertical device on a (0001) wafer, so its
current and field are along the c-axis (the 4 deg off-cut is neglected). Since
2026-10-21 the solvers' default low-field set is ``"c_axis"`` (parallel to c: electrons
from the drift mobility of Ishikawa's 2025 Kyoto thesis; holes unchanged, the primaries
disagree on them; the opt-in ``"c_axis_holes_0p9"`` takes them x 0.90); ``"legacy_perp"`` is the previous default, the Caughey-Thomas set measured
PERPENDICULAR to c (Schaffer 1994 via Ayalew Table 3.5, Hall data in the basal plane),
which ETNA used for vertical transport until then (DeviceConfig(mobility=...) selects
the set of the device; ``low_field_mobility="etna_default"`` here means "the device's
own"). An ``anisotropy`` model rotates a perpendicular set only. Impact-ionization sets
carry their own direction; the 1D solvers use the parallel-to-c coefficients and refuse
a set that has none.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, fields, replace
from typing import Literal, Optional, Tuple, Union

import numpy as np

from etna.core.sic_material import DEFAULT_MOBILITY, SiC4H_Parameters

Status = Literal["verified-primary", "verified-secondary", "assumption", "legacy"]
STATUSES = ("verified-primary", "verified-secondary", "assumption", "legacy")

REVIEW_2024 = (
    "J. Burin et al., 'TCAD parameters for 4H-SiC: a review', arXiv:2410.06798 (v8, 2025)"
)
REVIEW_2024_URL = "https://arxiv.org/abs/2410.06798"
AYALEW_2004 = (
    "T. Ayalew, 'SiC semiconductor devices technology, modeling and simulation', "
    "PhD thesis, TU Wien (2004)"
)
AYALEW_URL = "https://www.iue.tuwien.ac.at/phd/ayalew/"


@dataclass(frozen=True)
class Provenance:
    """Where a parameter set comes from. ``source`` is what was actually opened."""

    citation: str  # the primary publication
    doi: Optional[str]
    status: Status
    source: str  # the document the numbers were read in (primary or secondary)
    note: str = ""

    def __post_init__(self):
        if self.status not in STATUSES:
            raise ValueError(f"Provenance.status must be one of {STATUSES}")

    def as_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in fields(self)}


# --------------------------------------------------------------------------- mobility


class MobilityExtrapolationWarning(UserWarning):
    """A low-field mobility set is evaluated outside the fitted range of its T law."""



@dataclass(frozen=True)
class LowFieldMobility:
    """Doping/temperature low-field mobility, N = total doping, t = T / 300 K.

    ``form="ct"`` (Caughey-Thomas, Ayalew eqs. 3.77-3.79):
        mu = mu_min t^beta + (mu_max t^gamma - mu_min t^beta) / (1 + (N / N_ref)^alpha)
    ``form="arora"`` (Ishikawa 2025 thesis Eqs. 2.10 / 4.5, an Arora-type form):
        mu = (mu_max - mu_min) t^gamma / (1 + (N / N_ref)^alpha t^delta) + mu_min t^beta
    The two forms coincide at 300 K. ``direction``: "c_axis" (parallel to c, the
    vertical (0001) diode) or "perp" (basal plane, Hall data on (0001) wafers).
    """

    name: str
    mu_max_n: float
    mu_min_n: float
    N_ref_n: float
    alpha_n: float
    gamma_n: float  # exponent of mu_max ("ct") / of (mu_max - mu_min) ("arora")
    beta_n: float  # exponent of mu_min (0 = T-independent)
    mu_max_p: float
    mu_min_p: float
    N_ref_p: float
    alpha_p: float
    gamma_p: float
    beta_p: float
    provenance: Provenance
    direction: str = "perp"
    form: str = "ct"
    delta_n: float = 0.0  # "arora" only: exponent on (N / N_ref)^alpha
    delta_p: float = 0.0
    T_range_K: Optional[Tuple[float, float]] = None  # fitted T range of the T law

    def __post_init__(self):
        if self.direction not in ("c_axis", "perp"):
            raise ValueError("LowFieldMobility.direction must be 'c_axis' or 'perp'")
        if self.form not in ("ct", "arora"):
            raise ValueError("LowFieldMobility.form must be 'ct' or 'arora'")
        if self.form == "ct" and (self.delta_n != 0.0 or self.delta_p != 0.0):
            raise ValueError("LowFieldMobility: delta is an 'arora'-form parameter")

    def _carrier(self, carrier):
        s = "n" if carrier == "electron" else "p" if carrier == "hole" else None
        if s is None:
            raise ValueError(f"carrier must be 'electron' or 'hole' (got {carrier!r})")
        return lambda k: getattr(self, f"{k}_{s}")

    def coefficients(self, T, carrier):
        """T-resolved constants (A, B, C, N_ref, alpha) of the evaluation
        mu = C + A / (1 + (N / N_ref)^alpha) ("ct", B = None) or
        mu = A / (1 + (N / N_ref)^alpha * B) + C ("arora"); python floats, shared by
        the numpy evaluation and the devsim expression (same order of operations)."""
        g = self._carrier(carrier)
        t = float(T) / 300.0
        mu_min = g("mu_min") * t ** g("beta") if g("beta") != 0.0 else g("mu_min")
        if self.form == "ct":
            mu_max = g("mu_max") * t ** g("gamma")
            return mu_max - mu_min, None, mu_min, g("N_ref"), g("alpha")
        A = (g("mu_max") - g("mu_min")) * t ** g("gamma")
        return A, t ** g("delta"), mu_min, g("N_ref"), g("alpha")

    def mobility(self, N, T, carrier):
        """mu (cm^2/Vs) at total doping N (cm^-3) and T (K); numpy, broadcasts."""
        g = self._carrier(carrier)
        N = np.asarray(N, dtype=float)
        if self.form == "ct":
            t = np.asarray(T, dtype=float) / 300.0
            mu_max = g("mu_max") * t ** g("gamma")
            mu_min = g("mu_min") * t ** g("beta") if g("beta") != 0.0 else g("mu_min")
            return mu_min + (mu_max - mu_min) / (1.0 + (N / g("N_ref")) ** g("alpha"))
        self.check_T(T)
        if np.ndim(T) != 0:  # array T (track engine): the same expression, broadcast
            t = np.asarray(T, dtype=float) / 300.0
            A = (g("mu_max") - g("mu_min")) * t ** g("gamma")
            C = g("mu_min") * t ** g("beta") if g("beta") != 0.0 else g("mu_min")
            return A / (1.0 + (N / g("N_ref")) ** g("alpha") * t ** g("delta")) + C
        A, B, C, N_ref, alpha = self.coefficients(T, carrier)
        return A / (1.0 + (N / N_ref) ** alpha * B) + C

    def check_T(self, T):
        """Warn when T leaves the fitted range of the T law
        (``T_range_K``; None = no stated range): the law is then an extrapolation."""
        if self.T_range_K is None:
            return
        lo, hi = self.T_range_K
        T = np.asarray(T, dtype=float)
        if T.size and (np.min(T) < lo - 1e-9 or np.max(T) > hi + 1e-9):
            import warnings

            warnings.warn(
                f"low-field mobility {self.name!r}: T outside its fitted range "
                f"{lo:g}-{hi:g} K (the T law is extrapolated; below 300 K the electron "
                "(T/300)^-3.4 grows fast)",
                MobilityExtrapolationWarning,
                stacklevel=3,
            )

    def devsim_expression(self, N_expr, T, carrier):
        """devsim string of ``mobility`` at the density expression ``N_expr`` (the T-
        dependent constants inlined at full precision, same order of operations)."""
        self.check_T(T)
        A, B, C, N_ref, alpha = self.coefficients(T, carrier)
        if self.form == "ct":
            # written mu_min + (mu_max(T) - mu_min) / (...) as the numpy evaluation
            g = self._carrier(carrier)
            mu_max = g("mu_max") * (float(T) / 300.0) ** g("gamma")
            return (
                f"({float(C)!r} + ({float(mu_max)!r} - {float(C)!r}) / "
                f"(1.0 + pow(({N_expr}) / {float(N_ref)!r}, {float(alpha)!r})))"
            )
        return (
            f"({float(A)!r} / (1.0 + pow(({N_expr}) / {float(N_ref)!r}, "
            f"{float(alpha)!r}) * {float(B)!r}) + {float(C)!r})"
        )


_P = SiC4H_Parameters(mobility="legacy_perp")
ISHIKAWA_2025 = (
    "R. Ishikawa, 'Low-Field and High-Field Anisotropic Carrier Transport in SiC', PhD "
    "thesis, Kyoto University (2025), doi:10.14989/doctor.k26017"
)
ISHIKAWA_2025_URL = "http://hdl.handle.net/2433/296110"
ISHIKAWA_SBD_2024 = (
    "R. Ishikawa, M. Kaneko, T. Kimoto, Solid State Phenom. 360, 205 (2024), "
    "doi:10.4028/p-MAe9pZ (CC BY; http://hdl.handle.net/2433/290505)"
)
_LEGACY_PERP_PROVENANCE = Provenance(
    "W. J. Schaffer et al., MRS Proc. 339, 595 (1994) via Ayalew Table 3.5",
    None, "verified-secondary", AYALEW_URL + " (node65, Table 3.5)",
    "PERPENDICULAR to c (Ayalew's Table 3.5 header is mu_perp,300; eq. 3.81 gives "
    "mu_n,par = 1.2 mu_n): ETNA's default until 2026-10-21, used for vertical transport. "
    "Differs from Ayalew in T: mu_min is held T-independent (Ayalew: (T/300)^-0.5). The "
    "electron gamma -2.40 is Schaffer's <0001> (parallel) value, -2.15 in the basal plane "
    "(review Table XXVIII footnote c); the hole -2.15 is borrowed from electrons "
    "(Schaffer stated it for electrons only, review text)",
)  # fmt: skip
LOW_FIELD_MOBILITY = {
    # ETNA's solver mobility until 2026-10-21 (perpendicular-to-c set, vertical use)
    "legacy_perp": LowFieldMobility(
        "legacy_perp", _P.mu_n_max, _P.mu_n_min, _P.N_ref_n, _P.alpha_n, _P.gamma_n,
        0.0, _P.mu_p_max, _P.mu_p_min, _P.N_ref_p, _P.alpha_p, _P.gamma_p, 0.0,
        _LEGACY_PERP_PROVENANCE, direction="perp",
    ),  # fmt: skip
    # ETNA's default since 2026-10-21: parallel to c (vertical (0001) diode)
    "c_axis": LowFieldMobility(
        "c_axis", 1190.0, 10.0, 1.0e17, 0.54, -2.40, 0.0,
        _P.mu_p_max, _P.mu_p_min, _P.N_ref_p, _P.alpha_p, _P.gamma_p, 0.0,
        Provenance(
            "electrons: " + ISHIKAWA_2025 + ", Eq. (2.16) (drift mobility parallel to c, "
            "Hall data / Hall scattering factor, (11-20) epilayers, N_D 5e14-3e18); holes: "
            "the legacy set (Schaffer 1994 via Ayalew Table 3.5), taken isotropic as Ayalew "
            "eq. 3.81",
            None, "verified-primary", ISHIKAWA_2025_URL + " (pp. 44, 51)",
            "electrons at 300 K: mu_max 1190, mu_min 10, N_ref 1e17, alpha 0.54 exactly as "
            "Eq. 2.16; cross-checks at low doping: vertical SBD drift 1070 +- 290 at 9.7e14 "
            "(Ishikawa 2024; Eq. 2.16 gives 1101); the thesis's own perpendicular drift fit "
            "(Eq. 2.17) equals the legacy set within 3.5 % up to 5e15; mu_par / mu_perp 1.16 "
            "(Kimoto JJAP 58, 018002 (2019), Hall) to 1.2 (Ayalew eq. 3.81; Eq. 2.16 / 2.17: "
            "1.20 at 1e15); uncertainty about +-8 % (Hall factor). The Hall fit (Eq. 2.11) "
            "and the pulsed-conductance values (Table 3.1, Hall density) are 8-30 % higher: "
            "the electron Hall factor. T law NOT the thesis's drift fit (-3.4 on mu_max - "
            "mu_min, 300-600 K only; Hall fits give -2.4 to -3.2): ETNA keeps Schaffer's "
            "<0001> gamma -2.40 on mu_max (the legacy T law, now direction-consistent). "
            "HOLES UNCHANGED (strict reading: the primaries disagree on the hole magnitude, "
            "Hall 83-124 cm2/Vs at low N_A, and on the anisotropy: Ishikawa 0.9, Hatakeyama "
            "2003 1/1.15 via the review, Ayalew isotropic): uncertainty -35 / +0 % along c. "
            "The 0.90 ratio is the opt-in set 'c_axis_holes_0p9'. Electrons in the p+ "
            "(Al-doped, minority) are an extrapolation of an N-doped fit. The status "
            "'verified-primary' covers the 300 K electron parameters only: the T exponent "
            "and the hole set are secondary (Schaffer 1994 via Ayalew Table 3.5)",
        ),
        direction="c_axis",
    ),  # fmt: skip
    # OPT-IN alternative: the c_axis electrons and holes x 0.90 (thesis drift ratio)
    "c_axis_holes_0p9": LowFieldMobility(
        "c_axis_holes_0p9", 1190.0, 10.0, 1.0e17, 0.54, -2.40, 0.0,
        112.5, 14.31, 1.76e19, 0.34, -2.15, 0.0,
        Provenance(
            "electrons as 'c_axis'; holes: the legacy set x 0.90, the drift-mobility ratio "
            "mu_par / mu_perp of " + ISHIKAWA_2025 + " (Sec. 4.3.1, resistivity ratio, about "
            "0.9 for N_A 1e16-7e18, 300-600 K)",
            None, "verified-primary", ISHIKAWA_2025_URL + " (pp. 44, 51, 95-97)",
            "OPT-IN, not the default: one opened primary for the hole ratio (Hatakeyama "
            "2003 1/1.15 secondary; Ayalew eq. 3.81 isotropic). The hole magnitude stays "
            "uncertain (Hall 83 parallel in the thesis vs 112.5 here)",
        ),
        direction="c_axis",
    ),  # fmt: skip
    "ishikawa2025": LowFieldMobility(
        "ishikawa2025", 1190.0, 10.0, 1.0e17, 0.54, -3.4, 0.8,
        83.0, 20.0, 6.4e18, 0.83, -2.3, -2.2,
        Provenance(
            ISHIKAWA_2025 + ", Eq. (2.16) (electron DRIFT mobility parallel to c) and Eq. "
            "(4.6) (hole HALL mobility parallel to c)",
            None, "verified-primary", ISHIKAWA_2025_URL + " (pp. 44, 51, 88)",
            "full thesis T laws, fitted 300-600 K (below 300 K an extrapolation; the "
            "electron (T/300)^-3.4 grows fast below 300 K). Holes are HALL mobilities (the "
            "hole Hall factor was not converted in the thesis; conductance values in Fig. "
            "5.8 are about 5-10 % lower)",
        ),
        direction="c_axis", form="arora", delta_n=-2.2, delta_p=-0.9,
        T_range_K=(300.0, 600.0),
    ),  # fmt: skip
    "ayalew2004": LowFieldMobility(
        "ayalew2004", 950.0, 40.0, 1.94e17, 0.61, -2.40, -0.5,
        125.0, 15.9, 1.76e19, 0.34, -2.15, -0.5,
        Provenance(
            AYALEW_2004 + ", Table 3.5, eqs. 3.77-3.79", None, "verified-primary",
            AYALEW_URL + " (node65)",
            "full Ayalew T law (beta -0.5 on mu_min); a mixture of Schaffer 1994 "
            "(N_ref, alpha) and Roschke & Schwierz 2001 (T exponents), perpendicular to c",
        ),
        direction="perp",
    ),  # fmt: skip
}
# "etna_default" = the device's own selection (DeviceConfig.mobility, default "c_axis"):
# the facades and the devsim twin leave the engine's mobility untouched for it
LOW_FIELD_MOBILITY["etna_default"] = replace(
    LOW_FIELD_MOBILITY[DEFAULT_MOBILITY], name="etna_default"
)


@dataclass(frozen=True)
class HighFieldMobility:
    """Caughey-Thomas/Canali high-field mobility. v_sat(T) = v_sat (T/300)^delta,
    beta(T) = beta (T/300)^gamma_beta. ``*_perp`` perpendicular, ``*_par`` parallel to
    c (the ``_par`` values are used only with an ``anisotropy`` model)."""

    name: str
    vsat_perp_n: float  # cm/s
    vsat_par_n: float
    beta_n: float
    delta_n: float
    gamma_beta_n: float
    vsat_perp_p: float
    vsat_par_p: float
    beta_p: float
    delta_p: float
    gamma_beta_p: float
    provenance: Provenance

    def vsat_beta(self, T, carrier, parallel=False):
        s = "n" if carrier == "electron" else "p"
        g = lambda k: getattr(self, f"{k}_{s}")
        t = float(T) / 300.0
        v = g("vsat_par" if parallel else "vsat_perp") * t ** g("delta")
        return v, g("beta") * t ** g("gamma_beta")

    def mobility(self, mu0, E, T, carrier, parallel=False):
        """mu (cm^2/Vs) at low-field mobility mu0 and |E| (V/cm); numpy reference of
        ``etna.core.batched.models_dd.edge_mobility_kernel``."""
        v, b = self.vsat_beta(T, carrier, parallel)
        mu0 = np.asarray(mu0, dtype=float)
        return mu0 / (1.0 + (mu0 * np.abs(E) / v) ** b) ** (1.0 / b)


HIGH_FIELD_MOBILITY = {
    "lades2000": HighFieldMobility(
        "lades2000",
        2.2e7, 1.64e7, 1.2, -0.44, 1.0,
        1.1e7, 6.5e6, 1.2, -0.44, 1.0,
        Provenance(
            "electrons: M. Lades, PhD thesis TU Munich (2000), fit to I. A. Khan and "
            "J. A. Cooper, IEEE TED 47, 269 (2000); holes and the parallel-to-c values: "
            "H.-E. Nilsson et al., Monte Carlo (2000)",
            None, "verified-secondary",
            AYALEW_URL + " (node67, Table 3.6) and " + REVIEW_2024_URL
            + " Tables XXXIII-XXXIV",
            "electron v_sat 2.2e7 / beta 1.2 / delta -0.44 / gamma_beta 1.0 are Ayalew "
            "Table 3.6 (perpendicular to c); v_sat parallel 1.64e7 (e) and hole v_sat "
            "1.1e7 / 6.5e6 are Nilsson 2000 MC (review Tables XXXIII/XXXIV). ASSUMPTION: "
            "the hole beta and T exponents are set equal to the electron ones (no hole "
            "fit is published, review Sec. VIII). Check against the measured parallel-to-c "
            "hole v(E) (Ishikawa 2025 thesis Fig. 5.7a, N_A 1e16, read off the plot): "
            "about 4-5e6 cm/s at 100 kV/cm and 5e6 at 200 kV/cm vs 4.1e6 / 5.2e6 from the "
            "parallel hole parameters with mu0 = 85",
        ),
    ),  # fmt: skip
    "ishikawa2025": HighFieldMobility(
        "ishikawa2025",
        1.7e7, 1.7e7, 0.76, 0.0, 0.0,
        1.1e7, 6.5e6, 1.2, -0.44, 1.0,
        Provenance(
            "electrons: " + ISHIKAWA_2025 + ", Ch. 3, Table 3.1 (pulsed conductance on "
            "(11-20) bridges, 1-200 kV/cm); holes: the lades2000 hole values",
            None, "verified-primary", ISHIKAWA_2025_URL + " (pp. 73-75)",
            "electron v_sat 1.7e7 cm/s in BOTH directions (mean of N_D <= 1e16, where "
            "self-heating is negligible; computed with the Hall electron density, so "
            "possibly over-estimated by the Hall factor > 1, thesis Sec. 3.4.3); beta 0.76 = "
            "mean of the parallel fits 0.62 / 0.78 / 0.87 at N_D 2e15-1e16 (perpendicular "
            "0.58-1.02); no T law measured (300 K values at every T). Holes: no saturation "
            "observed up to 200 kV/cm (thesis Ch. 5), so the lades2000 hole set is kept",
        ),
    ),  # fmt: skip
}


@dataclass(frozen=True)
class Anisotropy:
    """mu_parallel / mu_perpendicular to c (low field), per carrier."""

    name: str
    ratio_n: float
    ratio_p: float
    provenance: Provenance


ANISOTROPY = {
    "ayalew2004": Anisotropy(
        "ayalew2004", 1.2, 1.0,
        Provenance(
            AYALEW_2004 + ", eqs. 3.80-3.81", None, "verified-primary",
            AYALEW_URL + " (node65)",
            "mu_n,par = 1.2 mu_n, mu_p,par = mu_p (holes isotropic: outdated, see the "
            "review)",
        ),
    ),  # fmt: skip
    "burin2024": Anisotropy(
        "burin2024", 1.0 / 0.83, 1.0 / 1.15,
        Provenance(
            "review consensus: electrons mu_perp/mu_par = 0.83 +- 0.07 (Schaffer 1994, "
            "Hatakeyama 2003, Nilsson 1996), holes 1.15 +- 0.04 (Hatakeyama 2003)",
            None, "verified-secondary", REVIEW_2024_URL + " Table XXVII and Sec. VIII",
        ),
    ),  # fmt: skip
}


# ----------------------------------------------------------------- impact ionization


@dataclass(frozen=True)
class AlphaParams:
    """alpha = a (1 + c (T - 300)) exp(-[b (1 + d (T - 300)) / E]^m) for one carrier
    and one direction (a 1/cm, b V/cm, c and d 1/K)."""

    a: float
    b: float
    m: float = 1.0
    c: float = 0.0
    d: float = 0.0

    def resolved(self, T):
        """(a(T), b(T), m) of the kernel form a exp(-(b/E)^m)."""
        dT = float(T) - 300.0
        return self.a * (1.0 + self.c * dT), self.b * (1.0 + self.d * dT), self.m

    def alpha(self, E, T=300.0):
        a, b, m = self.resolved(T)
        E = np.abs(np.asarray(E, dtype=float))
        with np.errstate(divide="ignore", over="ignore"):
            return np.where(E > 1.0, a * np.exp(-((b / np.maximum(E, 1.0)) ** m)), 0.0)


@dataclass(frozen=True)
class ImpactIonizationModel:
    """Electron/hole coefficients parallel (``*_par``) and, if published,
    perpendicular (``*_perp``) to c. ``T_range_K`` is the range the coefficients were
    measured over (None = 300 K only); ``field_range_V_cm`` the fitted field range."""

    name: str
    n_par: Optional[AlphaParams]
    p_par: Optional[AlphaParams]
    n_perp: Optional[AlphaParams]
    p_perp: Optional[AlphaParams]
    T_range_K: Optional[Tuple[float, float]]
    field_range_V_cm: Tuple[float, float]
    provenance: Provenance

    def coefficients(self, T, direction="c_axis", allow_T_extrapolation=False):
        """(a_n, b_n, k_n, m_n, a_p, b_p, k_p, m_p) at T for the kernels
        (``etna.core.avalanche``); k = 0 (no E prefactor) for every set here."""
        if direction not in ("c_axis", "basal"):
            raise ValueError(f"direction must be 'c_axis' or 'basal' (got {direction!r})")
        n, p = (self.n_par, self.p_par) if direction == "c_axis" else (
            self.n_perp, self.p_perp
        )
        if n is None or p is None:
            raise ValueError(
                f"impact ionization {self.name!r} has no {direction} coefficients"
            )
        self.check_T(T, allow_T_extrapolation)
        an, bn, mn = n.resolved(T)
        ap, bp, mp = p.resolved(T)
        if not (an > 0 and bn > 0 and ap > 0 and bp > 0):
            raise ValueError(
                f"impact ionization {self.name!r}: non-positive coefficient at T={T} K "
                "(the linear temperature law has left its range)"
            )
        return (an, bn, 0.0, mn, ap, bp, 0.0, mp)

    def check_T(self, T, allow_T_extrapolation=False):
        lo, hi = self.T_range_K if self.T_range_K is not None else (300.0, 300.0)
        if not (lo - 1e-9 <= T <= hi + 1e-9) and not allow_T_extrapolation:
            raise ValueError(
                f"impact ionization {self.name!r} was measured at {lo:g}-{hi:g} K; "
                f"T={T:g} K needs allow_T_extrapolation=True (the coefficients are then "
                "used outside their range: an extrapolation, not a model prediction)"
            )

    def alpha(self, E, T=300.0, carrier="electron", direction="c_axis"):
        """numpy alpha (1/cm) at |E| (V/cm): reference of the numba kernel."""
        par = direction == "c_axis"
        prm = {
            ("electron", True): self.n_par, ("hole", True): self.p_par,
            ("electron", False): self.n_perp, ("hole", False): self.p_perp,
        }[(carrier, par)]  # fmt: skip
        if prm is None:
            raise ValueError(f"{self.name!r} has no {direction} {carrier} coefficients")
        return prm.alpha(E, T)


def hatakeyama_anisotropic_alpha(model, F_c, F_basal, carrier="electron"):
    """Hatakeyama's anisotropic interpolation (SISPAD 2005, Eq. 1) of a model with both
    directions (``hatakeyama2004``): alpha for the field components F_c (along c) and
    F_basal (in the basal plane), V/cm.

        F^2 = F_b^2 + F_c^2,   F^2 / B^2 = F_b^2 / b_b^2 + F_c^2 / b_c^2,
        a = a_b^(B^2 F_b^2 / (b_b^2 F^2)) a_c^(B^2 F_c^2 / (b_c^2 F^2)),
        A = ln(a_c / a_b),
        alpha = a exp(-(B/F) sqrt(1 - A^2 (B F_c F_b / (F b_c b_b))^2))

    For 2D/3D use (not yet coupled into the 2D solvers); reduces to the parallel /
    perpendicular Chynoweth law on the axes (tests/test_physics_models.py).
    """
    if carrier == "electron":
        c, b_ = model.n_par, model.n_perp
    else:
        c, b_ = model.p_par, model.p_perp
    if c is None or b_ is None:
        raise ValueError(f"{model.name!r} needs both directions")
    Fc = np.abs(np.asarray(F_c, dtype=float))
    Fb = np.abs(np.asarray(F_basal, dtype=float))
    F2 = Fc**2 + Fb**2
    with np.errstate(divide="ignore", invalid="ignore"):
        B2 = F2 / (Fb**2 / b_.b**2 + Fc**2 / c.b**2)
        wb = B2 * Fb**2 / (b_.b**2 * F2)
        wc = B2 * Fc**2 / (c.b**2 * F2)
        a = np.exp(wb * np.log(b_.a) + wc * np.log(c.a))
        A = math.log(c.a / b_.a)
        F = np.sqrt(F2)
        B = np.sqrt(B2)
        root = np.sqrt(np.maximum(1.0 - A**2 * (B * Fc * Fb / (F * c.b * b_.b)) ** 2, 0.0))
        out = a * np.exp(-(B / F) * root)
    return np.where(F2 > 0, out, 0.0)


_HATA = Provenance(
    "T. Hatakeyama, T. Watanabe, T. Shinohe, K. Kojima, K. Arai, N. Sano, 'Impact "
    "ionization coefficients of 4H silicon carbide', Appl. Phys. Lett. 85, 1380 (2004), "
    "Table I",
    "10.1063/1.1784520", "verified-primary",
    "https://tsukuba.repo.nii.ac.jp/record/18553/files/APL_85-8.pdf",
    "300 K only; photomultiplication of p+n mesas on (0001) and (11-20) faces, "
    "N_D 3e16-2e17. The anisotropic interpolation is Hatakeyama, SISPAD 2005, pp. "
    "171-174 (http://in4.iue.tuwien.ac.at/pdfs/sispad2005/1562052.pdf)",
)  # fmt: skip

IMPACT_IONIZATION = {
    "hatakeyama2004": ImpactIonizationModel(
        "hatakeyama2004",
        AlphaParams(1.76e8, 3.30e7), AlphaParams(3.41e8, 2.50e7),
        AlphaParams(2.10e7, 1.70e7), AlphaParams(2.96e7, 1.60e7),
        None, (2.0e6, 5.0e6), _HATA,
    ),  # fmt: skip
    "niwa2015": ImpactIonizationModel(
        "niwa2015",
        AlphaParams(0.143e6, 4.93e6, m=2.37),
        AlphaParams(3.14e6, 11.8e6, m=1.02, c=6.3e-3, d=1.23e-3),
        None, None, (300.0, 448.0), (1.0e6, 3.3e6),
        Provenance(
            "H. Niwa, J. Suda, T. Kimoto, 'Impact ionization coefficients in 4H-SiC "
            "toward ultrahigh-voltage power devices', IEEE TED 62, 3326 (2015)",
            None, "verified-secondary",
            REVIEW_2024_URL + " Table XIII (Eq. 42)",
            "parallel to c; electrons 2.2-3.3 MV/cm (no T dependence detected), holes "
            "1-2.35 MV/cm. The review notes a_p is printed as both 3.12e6 and 3.14e6 and "
            "Steinmann (JAP 133, 235705, 2023) refits c_p 2.05e-3, d_p 0.65e-3. "
            "ASSUMPTION: T range 300-448 K (the primary's measured range was not "
            "opened; outside it pass allow_T_extrapolation)",
        ),
    ),  # fmt: skip
    "zhao2019": ImpactIonizationModel(
        "zhao2019",
        AlphaParams(0.339e6, 5.15e6, m=2.37),
        AlphaParams(3.56e6, 11.7e6, m=1.02, c=6.19e-3, d=1.15e-3),
        None, None, (150.0, 673.0), (1.0e6, 3.2e6),
        Provenance(
            "D. Zhao, H. Niwa, T. Kimoto, Jpn. J. Appl. Phys. 58, 018001 (2019)",
            None, "verified-secondary",
            REVIEW_2024_URL + " Table XIII",
            "parallel to c, 150-673 K. The review warns that the hole c, d make beta rise "
            "from 150 to 300 K, which the paper's own plots do not show. The primary was not "
            "opened (citation as listed by the review)",
        ),
    ),  # fmt: skip
    "loh2008": ImpactIonizationModel(
        "loh2008",
        AlphaParams(2.78e6, 10.5e6, m=1.37), AlphaParams(3.51e6, 10.3e6, m=1.09),
        None, None, None, (1.0e6, 5.0e6),
        Provenance(
            "W. S. Loh et al., 'Impact ionization coefficients in 4H-SiC', IEEE TED 55, "
            "1984 (2008)",
            None, "verified-secondary",
            REVIEW_2024_URL + " Table XIII; A. S. Kyuregyan, arXiv:1506.00980",
            "direction not stated in the review (single values); ASSUMPTION: used as "
            "parallel to c (OBIC on (0001) epi). The primary was not opened",
        ),
    ),  # fmt: skip
    "stefanakis2021": ImpactIonizationModel(
        "stefanakis2021",
        AlphaParams(2.8e6, 20.7e6, c=-1.0e-3, d=-0.29e-3),
        AlphaParams(2.5e6, 12.1e6, c=1.74e-3, d=0.59e-3),
        None, None, (300.0, 573.0), (1.25e6, 10.0e6),
        Provenance(
            "D. Stefanakis et al., global fit of published 4H-SiC data, IEEE TED (2021); "
            "temperature coefficients supplied by Steinmann et al. (review footnote m)",
            None, "verified-secondary", REVIEW_2024_URL + " Table XIII",
            "parallel to c, 1.25-10 MV/cm. ASSUMPTION: T range 300-573 K (not opened)",
        ),
    ),  # fmt: skip
    "ayalew2004_legacy": ImpactIonizationModel(
        "ayalew2004_legacy",
        AlphaParams(3.44e6, 2.58e7), AlphaParams(3.5e6, 1.7e7),
        None, None, None, (1.0e6, 5.0e6),
        Provenance(
            AYALEW_2004 + ", Sec. 3.5.3 (Chynoweth, 4H)", None, "legacy",
            AYALEW_URL + " (node74)",
            "LEGACY, not recommended: the review (App. A4) finds b_n changed from the 6H "
            "2.58e6 to 2.58e7 V/cm and the anisotropy 'most likely derived for 6H'; kept "
            "because TCAD tools ship it as a default. The optical-phonon T factor is not "
            "applied (300 K)",
        ),
    ),  # fmt: skip
}


# ---------------------------------------------------------------- bandgap narrowing


@dataclass(frozen=True)
class BandgapNarrowing:
    """Slotboom-type Delta E_g = C (F + sqrt(F^2 + G)), F = ln(N / N_ref) (eV), for
    n-type (donors) and p-type (acceptors) material."""

    name: str
    C_n: float
    C_p: float
    N_ref: float
    G: float
    provenance: Provenance

    def delta_Eg(self, N, carrier_type):
        C = self.C_n if carrier_type == "n" else self.C_p
        N = np.maximum(np.asarray(N, dtype=float), 1e-300)
        F = np.log(N / self.N_ref)
        return C * (F + np.sqrt(F * F + self.G))


BANDGAP_NARROWING = {
    "lades2000": BandgapNarrowing(
        "lades2000", 2.0e-2, 9.0e-3, 1.0e17, 0.5,
        Provenance(
            "M. Lades, PhD thesis TU Munich (2000), Slotboom fit to U. Lindefelt, J. "
            "Appl. Phys. 84, 2628 (1998)",
            None, "verified-primary", AYALEW_URL + " (node63, Table 3.4)",
            "Ayalew eq. 3.72 / Table 3.4 (the thesis is the opened source; G = 0.5 per "
            "the review Table XII). The two secondary restatements of Lindefelt's own "
            "four-term formula disagree on exponents and signs, so that form is not "
            "offered",
        ),
    ),  # fmt: skip
}


# ------------------------------------------------------------------- the selection

_REGISTRIES = {
    "low_field_mobility": (LOW_FIELD_MOBILITY, LowFieldMobility),
    "high_field_mobility": (HIGH_FIELD_MOBILITY, HighFieldMobility),
    "anisotropy": (ANISOTROPY, Anisotropy),
    "impact_ionization": (IMPACT_IONIZATION, ImpactIonizationModel),
    "bandgap_narrowing": (BANDGAP_NARROWING, BandgapNarrowing),
}


@dataclass(frozen=True)
class PhysicsModels:
    """Model selection (registry names or parameter-set instances); frozen, hashable.

    ``PhysicsModels()`` = ETNA's existing physics (the device's own low-field set,
    "c_axis" by default, no velocity saturation, no avalanche, no bandgap narrowing). ``allow_T_extrapolation`` lets an impact-
    ionization set run outside its measured temperature range (flagged in metadata).
    """

    low_field_mobility: Union[str, LowFieldMobility] = "etna_default"
    high_field_mobility: Union[None, str, HighFieldMobility] = None
    anisotropy: Union[None, str, Anisotropy] = None
    impact_ionization: Union[None, str, ImpactIonizationModel] = None
    bandgap_narrowing: Union[None, str, BandgapNarrowing] = None
    allow_T_extrapolation: bool = False
    _resolved: dict = field(default=None, init=False, repr=False, compare=False, hash=False)

    def __post_init__(self):
        problems = []
        out = {}
        for key, (reg, cls) in _REGISTRIES.items():
            v = getattr(self, key)
            if v is None:
                if key == "low_field_mobility":
                    problems.append("low_field_mobility cannot be None")
                out[key] = None
            elif isinstance(v, str):
                if v not in reg:
                    problems.append(
                        f"{key}={v!r} is not a registered model (exact names: "
                        f"{sorted(reg)})"
                    )
                else:
                    out[key] = reg[v]
            elif isinstance(v, cls):
                out[key] = v
            else:
                problems.append(
                    f"{key} must be None, a name in {sorted(reg)} or a {cls.__name__} "
                    f"(got {type(v).__name__})"
                )
        if not isinstance(self.allow_T_extrapolation, bool):
            problems.append("allow_T_extrapolation must be a bool")
        lf, an = out.get("low_field_mobility"), out.get("anisotropy")
        if lf is not None and an is not None and lf.direction != "perp":
            problems.append(
                f"anisotropy={an.name!r} rotates a PERPENDICULAR-to-c low-field set to "
                f"the c-axis, but low_field_mobility={lf.name!r} is already parallel to c "
                "(the default since 2026-10-21): pass low_field_mobility='legacy_perp' "
                "(or 'ayalew2004') with an anisotropy, or no anisotropy"
            )
        ii = out.get("impact_ionization")
        if ii is not None and (ii.n_par is None or ii.p_par is None):
            problems.append(
                f"impact_ionization {ii.name!r} has no parallel-to-c coefficients (the "
                "1D diode's field is along c)"
            )
        if problems:
            raise ValueError("Invalid PhysicsModels: " + "; ".join(problems))
        object.__setattr__(self, "_resolved", out)

    def get(self, key):
        """The resolved parameter-set instance (or None) of ``key``."""
        return self._resolved[key]

    @property
    def is_default(self) -> bool:
        return (
            self.get("low_field_mobility") is LOW_FIELD_MOBILITY["etna_default"]
            and all(
                self.get(k) is None
                for k in ("high_field_mobility", "anisotropy", "impact_ionization",
                          "bandgap_narrowing")
            )
        )  # fmt: skip

    def check_temperature(self, T):
        ii = self.get("impact_ionization")
        if ii is not None:
            ii.check_T(T, self.allow_T_extrapolation)

    def T_extrapolated(self, T) -> bool:
        ii = self.get("impact_ionization")
        if ii is None:
            return False
        lo, hi = ii.T_range_K if ii.T_range_K is not None else (300.0, 300.0)
        return not (lo - 1e-9 <= T <= hi + 1e-9)

    # --- solver inputs (1D vertical = c-axis) ---------------------------------------

    def low_field_set(self, device_mobility=None):
        """The low-field set in effect: the selected one, or for "etna_default" the
        device's own (``device_mobility``, DeviceConfig.mobility; None = default)."""
        lf = self.get("low_field_mobility")
        if lf.name == "etna_default":
            return LOW_FIELD_MOBILITY[device_mobility or DEFAULT_MOBILITY]
        return lf

    def vertical_is_parallel(self, device_mobility=None) -> bool:
        """True when the vertical transport is parallel to c: a c-axis low-field set, or
        a perpendicular one rotated by an ``anisotropy`` model. Selects the parallel
        high-field v_sat."""
        if self.get("anisotropy") is not None:
            return True
        return self.low_field_set(device_mobility).direction == "c_axis"

    def low_field_edge_mobility(self, donors, acceptors_total, T, carrier,
                                device_mobility=None):  # fmt: skip
        """Low-field node mobility of the selected set at total doping, times the
        anisotropy ratio (vertical = parallel to c) when one is selected."""
        mu = self.low_field_set(device_mobility).mobility(
            np.asarray(donors) + np.asarray(acceptors_total), T, carrier
        )
        an = self.get("anisotropy")
        if an is not None:
            mu = mu * (an.ratio_n if carrier == "electron" else an.ratio_p)
        return mu

    def high_field_params(self, T, device_mobility=None):
        """(vsat_n, beta_n, vsat_p, beta_p) for vertical transport, or None: the
        parallel-to-c v_sat when the vertical transport is parallel to c
        (``vertical_is_parallel``), the perpendicular one with a "legacy_perp" device."""
        hf = self.get("high_field_mobility")
        if hf is None:
            return None
        par = self.vertical_is_parallel(device_mobility)
        vn, bn = hf.vsat_beta(T, "electron", parallel=par)
        vp, bp = hf.vsat_beta(T, "hole", parallel=par)
        return vn, bn, vp, bp

    def impact_coefficients(self, T):
        ii = self.get("impact_ionization")
        if ii is None:
            return None
        return ii.coefficients(T, "c_axis", self.allow_T_extrapolation)

    def describe(self, device_mobility=None) -> dict:
        """JSON-able provenance record of the selection (for result metadata)."""
        d = {}
        for key in _REGISTRIES:
            m = self.get(key)
            d[key] = None if m is None else {
                "name": m.name, **m.provenance.as_dict()
            }  # fmt: skip
        lf = self.low_field_set(device_mobility)
        if self.get("low_field_mobility").name == "etna_default":
            d["low_field_mobility"]["device_set"] = lf.name
        d["allow_T_extrapolation"] = self.allow_T_extrapolation
        if self.get("anisotropy") is not None:
            d["transport_direction"] = "parallel to c (perpendicular set x anisotropy)"
        elif lf.direction == "c_axis":
            d["transport_direction"] = f"parallel to c ({lf.name!r} set)"
        else:
            d["transport_direction"] = (
                f"vertical, with the perpendicular-to-c set {lf.name!r} (legacy)"
            )
        return d


def resolve_models(models) -> PhysicsModels:
    """None -> PhysicsModels(); a PhysicsModels as is; a dict -> PhysicsModels(**d)."""
    if models is None:
        return PhysicsModels()
    if isinstance(models, PhysicsModels):
        return models
    if isinstance(models, dict):
        return PhysicsModels(**models)
    raise TypeError(
        f"models must be None, a PhysicsModels or a dict (got {type(models).__name__})"
    )


__all__ = [
    "Provenance", "LowFieldMobility", "MobilityExtrapolationWarning",
    "HighFieldMobility", "Anisotropy", "AlphaParams",
    "ImpactIonizationModel", "BandgapNarrowing", "PhysicsModels", "resolve_models",
    "LOW_FIELD_MOBILITY", "HIGH_FIELD_MOBILITY", "ANISOTROPY", "IMPACT_IONIZATION",
    "BANDGAP_NARROWING", "hatakeyama_anisotropic_alpha", "replace",
]  # fmt: skip
