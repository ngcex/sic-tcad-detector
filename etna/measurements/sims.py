"""SIMS / SRP depth profiles: loader, validation, resampling into a ``donor_table``.

File (template ``data/measurements/templates/sims_template.csv``)::

    # technique: SIMS            (SIMS = total N / Al atoms; SRP = net carrier density)
    # junction_depth_um: <xj>    (depth of the p+/n- junction below the surface; optional)
    # detection_limit_N_cm3: <DL>   (optional; also as a per-row column)
    # rel_sigma_N: <fraction>    (optional global relative 1-sigma of N)
    depth_um,N_cm3,sigma_N_cm3,Al_cm3,sigma_Al_cm3

Depth is measured from the SURFACE (p+ top). ``DeviceConfig.donor_table`` depth is below
the p+/n- junction, so ``resample_profile`` needs the junction depth: the argument, else
the file's ``junction_depth_um``, else the depth where the Al and N curves cross (both
columns present), else ``config.substrate_thickness_um`` (the model's p+ thickness); the
choice is recorded. The table is truncated at the epi end (``config.epi_thickness_um``
below the junction): the n+ substrate belongs to ``N_D_nplus`` / ``nplus_thickness_um``,
not to the table.

Resampling (``method``): "bin_mean" (default, noise-robust) = at each knot the mean of
log10 N over the bin half-way to the neighbouring knots (log-linear interpolation of the
data when a bin is empty); it smooths structure narrower than the knot spacing (a steep
n-/n+ rise is rounded by up to a bin). "interp" = log-linear interpolation of the data at
the knot (exact on a noise-free profile, noise passed through). Both give the log-linear
form ``donor_table_spec`` interpolates between knots. Uncertainty per knot (relative 1-sigma of N): the
mean relative measurement sigma in the bin (SIMS calibration errors are correlated along
depth, so it is NOT divided by sqrt(n)) combined with the scatter of log10 N in the bin /
sqrt(n). Knots at or below the detection limit are flagged ``censored`` (the value is an
upper bound): ``to_device_config`` refuses them unless ``allow_censored=True`` (a SIMS
floor is never silently used as the epi doping).

Acceptors: ``DeviceConfig`` has a scalar ``N_A`` and a p+ thickness only (no acceptor
table). From an Al column the resampling reports the p+ plateau (median Al above the
junction), the junction depth and the Al profile itself as data; ``to_device_config``
applies the junction depth as ``substrate_thickness_um`` and, with ``use_acceptor``, the
plateau as ``N_A``. The Al tail into the epi is not modelled.

SIMS measures total N (incl. any inactive N); SRP measures the net free-carrier density
(|N_D - N_A| with the ionisation at the measurement T). The file's ``technique`` is
recorded; the table is used as N_D in both cases.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np

from etna.measurements._csv import (
    Column,
    MeasurementFileError,
    parse_table,
    read_text,
)
from etna.measurements._common import (
    ENGINE_TABLE_N_D,
    data_provenance,
    dumps,
    model_provenance,
)

SIMS_COLUMNS: "dict[str, Column]" = {
    "depth_um": Column("depth_um", "depth", 1.0, "um", "um"),
    "depth_nm": Column("depth_nm", "depth", 1e-3, "nm", "um"),
}
for _el in ("N", "Al"):
    SIMS_COLUMNS[f"{_el}_cm3"] = Column(f"{_el}_cm3", _el, 1.0, "cm^-3", "cm^-3")
    SIMS_COLUMNS[f"sigma_{_el}_cm3"] = Column(
        f"sigma_{_el}_cm3", f"s{_el}", 1.0, "cm^-3", "cm^-3"
    )
    SIMS_COLUMNS[f"rel_sigma_{_el}"] = Column(
        f"rel_sigma_{_el}", f"r{_el}", 1.0, "fraction", "fraction"
    )
    SIMS_COLUMNS[f"detection_limit_{_el}_cm3"] = Column(
        f"detection_limit_{_el}_cm3", f"dl{_el}", 1.0, "cm^-3", "cm^-3"
    )

SIMS_META = {
    "technique": str,
    "quantity": str,
    "depth_reference": str,
    "junction_depth_um": float,
    "detection_limit_n_cm3": float,
    "detection_limit_al_cm3": float,
    "rel_sigma_n": float,
    "rel_sigma_al": float,
    "sample": str,
    "lab": str,
    "date": str,
    "instrument": str,
    "source": str,
    "notes": str,
}
TECHNIQUES = ("SIMS", "SRP")


@dataclass
class DepthProfile:
    """A validated SIMS / SRP profile (depth from the surface, um; cm^-3)."""

    depth_um: np.ndarray
    N: "np.ndarray | None"
    Al: "np.ndarray | None"
    rel_sigma_N: "np.ndarray | None"  # relative 1-sigma per row (NaN = not given)
    rel_sigma_Al: "np.ndarray | None"
    detection_limit_N: "np.ndarray | None"  # per row (NaN = not given)
    detection_limit_Al: "np.ndarray | None"
    technique: str
    source: str
    sha256: str
    meta: dict = field(default_factory=dict)
    comments: list = field(default_factory=list)
    input_columns: dict = field(default_factory=dict)

    def provenance(self) -> dict:
        md = dict(self.meta)
        md.update(technique=self.technique, columns=self.input_columns)
        return data_provenance(
            self.source, self.sha256, f"{self.technique} depth profile", md
        )


def load_sims_csv(src) -> DepthProfile:
    """Load and validate a SIMS / SRP depth-profile CSV (path, bytes, file-like or text).

    Depth strictly increasing and >= 0; concentrations >= 0 (0 = below detection);
    at least one of N_cm3 / Al_cm3. Raises ``MeasurementFileError`` listing every
    problem with its file line and column.
    """
    text, source = read_text(src)
    raw = parse_table(
        text,
        source,
        SIMS_COLUMNS,
        required=[("depth",), ("N", "Al")],
        exclusive=[("sN", "rN"), ("sAl", "rAl")],
        meta_keys=SIMS_META,
        kind="SIMS/SRP",
    )
    cols, names, issues = raw.columns, raw.column_names, []
    lines = np.asarray(raw.lines)
    d = np.asarray(cols["depth"], float)
    for ln, v in zip(lines, d):
        if v < 0:
            issues.append(
                f"line {ln}, column '{names['depth']}': depth must be >= 0 (got {v:g})"
            )
    for k in np.where(np.diff(d) <= 0)[0]:
        issues.append(
            f"line {lines[k + 1]}, column '{names['depth']}': depth {d[k + 1]:g} is not "
            f"above {d[k]:g} on line {lines[k]} (depth must increase strictly)"
        )
    for q in ("N", "Al", "sN", "sAl", "rN", "rAl", "dlN", "dlAl"):
        if q not in cols:
            continue
        for ln, v in zip(lines, cols[q]):
            if not math.isnan(v) and v < 0:
                issues.append(
                    f"line {ln}, column '{names[q]}': must be >= 0 (got {v:g})"
                )
    for q in ("N", "Al"):
        if q in cols and not np.any(np.asarray(cols[q]) > 0):
            issues.append(f"column '{names[q]}': no value above zero")
    tech = raw.meta.get("technique", "SIMS").strip().upper()
    if tech not in TECHNIQUES:
        issues.append(f"metadata technique must be one of {TECHNIQUES} (got {tech!r})")
    ref = raw.meta.get("depth_reference", "surface").strip().lower()
    if ref != "surface":
        issues.append(
            f"metadata depth_reference must be 'surface' (got {ref!r}); shift the depth "
            "column so 0 is the device surface"
        )
    xj = raw.meta.get("junction_depth_um")
    if xj is not None and not xj > 0:
        issues.append(f"metadata junction_depth_um must be > 0 (got {xj:g})")
    for k in (
        "rel_sigma_n",
        "rel_sigma_al",
        "detection_limit_n_cm3",
        "detection_limit_al_cm3",
    ):
        v = raw.meta.get(k)
        if v is not None and v < 0:
            issues.append(f"metadata {k} must be >= 0 (got {v:g})")
    if issues:
        raise MeasurementFileError(source, issues)

    def arr(q):
        return np.asarray(cols[q], float) if q in cols else None

    def rel(el):
        n = arr(el)
        if n is None:
            return None
        if f"s{el}" in cols:
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.where(n > 0, arr(f"s{el}") / n, np.nan)
        if f"r{el}" in cols:
            return arr(f"r{el}")
        g = raw.meta.get(f"rel_sigma_{el.lower()}")
        return None if g is None else np.full(len(n), g)

    def dl(el):
        n = arr(el)
        if n is None:
            return None
        if f"dl{el}" in cols:
            return arr(f"dl{el}")
        g = raw.meta.get(f"detection_limit_{el.lower()}_cm3")
        return None if g is None else np.full(len(n), g)

    return DepthProfile(
        depth_um=d,
        N=arr("N"),
        Al=arr("Al"),
        rel_sigma_N=rel("N"),
        rel_sigma_Al=rel("Al"),
        detection_limit_N=dl("N"),
        detection_limit_Al=dl("Al"),
        technique=tech,
        source=source,
        sha256=raw.sha256,
        meta=dict(raw.meta),
        comments=raw.comments,
        input_columns=dict(names),
    )


def default_sims_knots_um(epi_um: float, spacing_um: float = 0.25) -> np.ndarray:
    """Uniform knots 0 .. epi every ``spacing_um`` (the epi end is always a knot)."""
    n = max(2, int(math.ceil(epi_um / spacing_um - 1e-9)) + 1)
    return np.linspace(0.0, float(epi_um), n)


def _crossing(depth, N, Al):
    """First depth where N rises above Al (log-linear between rows); None if absent."""
    with np.errstate(divide="ignore"):
        g = np.log10(np.maximum(N, 1e-300)) - np.log10(np.maximum(Al, 1e-300))
    for i in range(1, len(g)):
        if g[i - 1] < 0 <= g[i]:
            t = -g[i - 1] / (g[i] - g[i - 1])
            return float(depth[i - 1] + t * (depth[i] - depth[i - 1]))
    return None


@dataclass
class ProfileTable:
    """A depth profile resampled onto donor_table knots (depth below the junction)."""

    knots_um: np.ndarray
    N_D: np.ndarray  # cm^-3 at the knots
    rel_uncertainty: np.ndarray  # relative 1-sigma of N_D (NaN = unknown)
    n_rows: np.ndarray  # data rows in each knot bin (0 = interpolated)
    censored: np.ndarray  # knot at / below the detection limit
    extrapolated: np.ndarray  # knot outside the measured depth (edge value held)
    junction_depth_um: float
    junction_source: str
    epi_thickness_um: float
    N_A_plateau: "float | None"  # median Al above the junction (p+ plateau)
    acceptor_depth_um: "np.ndarray | None"  # Al profile as data (depth from surface)
    acceptor_cm3: "np.ndarray | None"
    data_provenance: dict
    method: str = "bin_mean"
    notes: list = field(default_factory=list)

    @property
    def donor_table(self):
        return tuple((float(d), float(n)) for d, n in zip(self.knots_um, self.N_D))

    def to_device_config(
        self,
        base,
        use_acceptor: bool = True,
        allow_censored: bool = False,
        allow_extrapolated: bool = True,
    ):
        """``base`` with this donor table (parametric n-/n+ transition off), the junction
        depth as the p+ thickness and, with ``use_acceptor`` and an Al plateau, N_A.

        Raises ValueError on censored knots (unless ``allow_censored``) or extrapolated
        knots (unless ``allow_extrapolated``, default True: edge values held).
        """
        if self.censored.any() and not allow_censored:
            raise ValueError(
                f"{int(self.censored.sum())} knot(s) at or below the detection limit "
                f"(depths {self.knots_um[self.censored].round(3).tolist()} um): the SIMS "
                "floor is an upper bound, not the epi doping; pass allow_censored=True to "
                "use it anyway"
            )
        if self.extrapolated.any() and not allow_extrapolated:
            raise ValueError(
                f"{int(self.extrapolated.sum())} knot(s) outside the measured depth"
            )
        if abs(base.epi_thickness_um - self.epi_thickness_um) > 1e-9:
            raise ValueError(
                f"table built for a {self.epi_thickness_um:g} um epi, base config has "
                f"{base.epi_thickness_um:g} um"
            )
        kw = dict(
            donor_table=self.donor_table,
            doping_profile="graded",
            N_D=None,
            buffer_thickness_um=None,
            buffer_end_doping=None,
            N_D_junction=float(self.N_D[0]),
            N_D_bulk=float(np.min(self.N_D)),
        )
        if self.junction_source != "config.substrate_thickness_um":
            kw["substrate_thickness_um"] = float(self.junction_depth_um)
        if use_acceptor and self.N_A_plateau is not None:
            kw["N_A"] = float(self.N_A_plateau)
        return replace(base, **kw)

    def to_dict(self, config=None) -> dict:
        return {
            "schema": "etna.donor_table/1",
            "method": f"sims_resample ({self.method})",
            "donor_table": [list(p) for p in self.donor_table],
            "rel_uncertainty": self.rel_uncertainty,
            "uncertainty_note": (
                "relative 1-sigma of N_D per knot: mean measurement sigma in the bin "
                "(not divided by sqrt(n), correlated calibration) combined with the bin "
                "scatter / sqrt(n); NaN where the file gives no sigma and the bin has < 2 rows"
            ),
            "censored": self.censored,
            "extrapolated": self.extrapolated,
            "rows_per_knot": self.n_rows,
            "junction_depth_um": self.junction_depth_um,
            "junction_source": self.junction_source,
            "epi_thickness_um": self.epi_thickness_um,
            "acceptor": {
                "N_A_plateau_cm3": self.N_A_plateau,
                "depth_from_surface_um": self.acceptor_depth_um,
                "Al_cm3": self.acceptor_cm3,
                "note": (
                    "DeviceConfig has a scalar N_A and a p+ thickness only: the Al tail "
                    "into the epi is reported, not modelled"
                ),
            },
            "notes": self.notes,
            "data": self.data_provenance,
            "model": model_provenance(config, "") if config is not None else {},
        }

    def to_json(self, config=None) -> str:
        return dumps(self.to_dict(config))


def resample_profile(
    profile: DepthProfile,
    config,
    knots_um=None,
    junction_depth_um: "float | None" = None,
    spacing_um: float = 0.25,
    method: str = "bin_mean",
) -> ProfileTable:
    """Resample the N profile onto donor_table knots below the junction of ``config``.

    ``knots_um`` (depth below the junction, from 0 to >= the epi) defaults to uniform
    knots every ``spacing_um``. See the module docstring for the junction choice, the
    resampling ``method`` ("bin_mean" | "interp") and the uncertainty.
    """
    if method not in ("bin_mean", "interp"):
        raise ValueError(f"method must be 'bin_mean' or 'interp' (got {method!r})")
    if profile.N is None:
        raise ValueError("the profile has no N_cm3 column: no donor table can be built")
    epi = float(config.epi_thickness_um)
    knots = (
        default_sims_knots_um(epi, spacing_um)
        if knots_um is None
        else np.asarray(knots_um, float)
    )
    from etna.core.device import donor_table_spec

    donor_table_spec(epi * 1e-4, [(k, 1e14) for k in knots])
    notes = []
    d = profile.depth_um
    if junction_depth_um is not None:
        xj, src = float(junction_depth_um), "argument"
    elif profile.meta.get("junction_depth_um") is not None:
        xj, src = float(profile.meta["junction_depth_um"]), "file metadata"
    elif (
        profile.Al is not None
        and (x := _crossing(d, profile.N, profile.Al)) is not None
    ):
        xj, src = x, "Al/N crossing"
    else:
        xj, src = float(config.substrate_thickness_um), "config.substrate_thickness_um"
    if not xj > 0:
        raise ValueError(f"junction depth must be > 0 um (got {xj!r})")
    notes.append(f"junction depth {xj:.4g} um from {src}")
    z = d - xj  # depth below the junction
    N = profile.N
    good = N > 0
    if good.sum() < 2:
        raise ValueError("fewer than two N rows above zero")
    zg, lg = z[good], np.log10(N[good])
    rel_s = profile.rel_sigma_N
    dl = profile.detection_limit_N
    edges = np.concatenate([[-np.inf], 0.5 * (knots[1:] + knots[:-1]), [np.inf]])
    K = len(knots)
    logN = np.empty(K)
    unc = np.full(K, np.nan)
    nrow = np.zeros(K, int)
    cens = np.zeros(K, bool)
    extra = (knots < z.min() - 1e-9) | (knots > z.max() + 1e-9)
    for k in range(K):
        lo_e = max(edges[k], 0.0 if k == 0 else edges[k])
        hi_e = min(edges[k + 1], epi if k == K - 1 else edges[k + 1])
        m = (
            (z >= lo_e) & (z < hi_e) & good
            if k < K - 1
            else (z >= lo_e) & (z <= hi_e) & good
        )
        nrow[k] = int(m.sum())
        if nrow[k] > 0 and method == "bin_mean":
            logN[k] = float(np.mean(np.log10(N[m])))
        else:
            logN[k] = float(np.interp(knots[k], zg, lg))
        parts = []
        if rel_s is not None:
            if nrow[k] > 0 and np.any(np.isfinite(rel_s[m])):
                parts.append(float(np.nanmean(rel_s[m])))
            else:
                fin = np.isfinite(rel_s) & good
                if fin.any():
                    parts.append(float(np.interp(knots[k], z[fin], rel_s[fin])))
        if nrow[k] >= 2:
            parts.append(
                float(
                    np.std(np.log10(N[m]), ddof=1) * np.log(10.0) / math.sqrt(nrow[k])
                )
            )
        if parts:
            unc[k] = math.sqrt(sum(p * p for p in parts))
        mk = (z >= lo_e) & (z <= hi_e)
        if dl is not None:
            dlk = dl[mk] if mk.any() else dl[[int(np.argmin(np.abs(z - knots[k])))]]
            dlk = dlk[np.isfinite(dlk)]
            if dlk.size and 10.0 ** logN[k] <= float(np.max(dlk)) * (1 + 1e-12):
                cens[k] = True
        if mk.any() and np.any(N[mk] <= 0):
            cens[k] = True
    ND = 10.0**logN
    if extra.any():
        notes.append(
            f"{int(extra.sum())} knot(s) outside the measured depth "
            f"({z.min():.3g}..{z.max():.3g} um below the junction): edge value held"
        )
    if cens.any():
        notes.append(
            f"{int(cens.sum())} knot(s) at or below the detection limit (upper bounds)"
        )
    if np.isnan(unc).any():
        notes.append(
            f"{int(np.isnan(unc).sum())} knot(s) without an uncertainty estimate"
        )
    lo, hi = ENGINE_TABLE_N_D
    if ND.min() < lo or ND.max() > hi:
        notes.append(
            f"N_D range [{ND.min():.3g}, {ND.max():.3g}] leaves the batched engine's "
            f"validated table range [{lo:g}, {hi:g}]: runs fall back to devsim"
        )
    if z.max() > epi:
        notes.append(
            f"data beyond the epi end ({epi:g} um below the junction) not used: the n+ is "
            "config.N_D_nplus"
        )
    NA = None
    adepth = acm = None
    if profile.Al is not None:
        adepth, acm = d.copy(), profile.Al.copy()
        above = (d < xj) & (profile.Al > 0)
        if profile.N is not None:
            above &= profile.Al > 10.0 * np.maximum(profile.N, 0)
        if above.any():
            NA = float(np.median(profile.Al[above]))
            notes.append(
                f"p+ plateau N_A ~ {NA:.3g} cm^-3 (median Al above the junction)"
            )
            if NA < 1e17:
                notes.append(
                    "Al plateau below 1e17 cm^-3: DeviceConfig requires N_A >= 1e17"
                )
    if profile.technique == "SRP":
        notes.append(
            "SRP: net free-carrier density, used as N_D (ionisation / compensation not separated)"
        )
    return ProfileTable(
        knots_um=knots,
        N_D=ND,
        rel_uncertainty=unc,
        n_rows=nrow,
        censored=cens,
        extrapolated=extra,
        junction_depth_um=xj,
        junction_source=src,
        epi_thickness_um=epi,
        N_A_plateau=NA,
        acceptor_depth_um=adepth,
        acceptor_cm3=acm,
        data_provenance=profile.provenance(),
        method=method,
        notes=notes,
    )


def model_donor_profile(config, depth_below_junction_um) -> np.ndarray:
    """N_D (cm^-3) of ``config``'s model at depths below the junction (um): the donor
    table, the uniform N_D or the graded profile with its n-/n+ transition (the numpy
    twins the solvers use, ``etna.core.device``), N_D_nplus beyond the epi end (0 there
    without an n+ layer)."""
    from etna.core.device import (
        DEFAULT_N_D,
        buffer_spec,
        donor_profile,
        donor_table_spec,
        tabulated_donor_profile,
    )

    z = np.asarray(depth_below_junction_um, float) * 1e-4  # cm
    epi = config.epi_thickness_um * 1e-4
    if config.donor_table is not None:
        spec = donor_table_spec(epi, config.donor_table)
        nd = tabulated_donor_profile(z, 0.0, spec, n_side=False)
    elif config.doping_profile == "uniform":
        nd = np.full(z.shape, DEFAULT_N_D if config.N_D is None else config.N_D)
    else:
        buf = None
        if config.buffer_thickness_um is not None:
            buf = buffer_spec(
                epi, epi, config.buffer_thickness_um * 1e-4, config.buffer_end_doping
            )
        nd = donor_profile(
            z,
            0.0,
            config.N_D_junction,
            config.N_D_bulk,
            config.L_transition_um * 1e-4,
            buffer=buf,
            n_side=False,
        )
    nplus = config.N_D_nplus if config.nplus_thickness_um is not None else 0.0
    return np.where(z > epi, nplus, nd)


def profile_vs_model(table: ProfileTable, config) -> dict:
    """rms and max of log10(N_table / N_model) at the knots (epi only)."""
    m = model_donor_profile(config, table.knots_um)
    r = np.log10(table.N_D / m)
    return {
        "rms_log10": float(np.sqrt(np.mean(r**2))),
        "max_abs_log10": float(np.max(np.abs(r))),
        "knots_um": table.knots_um,
        "log10_ratio": r,
        "N_model": m,
    }
