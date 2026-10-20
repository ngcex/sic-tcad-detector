"""Donor-table JSON (the file the C-V inversion / SIMS resampling writes) and reports.

JSON schema ``etna.donor_table/1``::

    {"schema": "etna.donor_table/1",
     "method": "...",                         # cv_inversion | sims_resample
     "donor_table": [[depth_um, N_D_cm3], ...],   # below the junction, from 0, covers the epi
     "rel_uncertainty": [...],                # relative 1-sigma per knot (null = unknown)
     ... method details, notes, limits ...,
     "data": {source_file, source_sha256, metadata},   # provenance of the measurement
     "model": {etna_version, config_hash, solver, date_utc}}

``load_donor_table_json`` returns the ``donor_table`` tuple ready for
``DeviceConfig(donor_table=...)``; the full record is in ``read_donor_table_json``.
"""

from __future__ import annotations

import json
from pathlib import Path

from etna.measurements._common import MEASUREMENT_LABEL, dumps, now_utc

SCHEMA = "etna.donor_table/1"


def read_donor_table_json(src) -> dict:
    """The full donor-table record from a path, JSON text or dict (schema-checked)."""
    if isinstance(src, dict):
        rec = src
    else:
        s = str(src)
        if isinstance(src, Path) or not s.lstrip().startswith("{"):
            s = Path(src).read_text()
        rec = json.loads(s)
    if rec.get("schema") != SCHEMA:
        raise ValueError(f"not an {SCHEMA} record (schema={rec.get('schema')!r})")
    tab = rec.get("donor_table")
    if not isinstance(tab, list) or len(tab) < 2:
        raise ValueError(
            "donor_table must be a list of at least two [depth_um, N_D] pairs"
        )
    for i, p in enumerate(tab):
        if not (isinstance(p, (list, tuple)) and len(p) == 2):
            raise ValueError(
                f"donor_table[{i}] must be a [depth_um, N_D] pair (got {p!r})"
            )
    return rec


def load_donor_table_json(src, epi_thickness_um: "float | None" = None):
    """``donor_table`` tuple ((depth_um, N_D), ...) from a donor-table JSON record.

    With ``epi_thickness_um`` the table is also validated against the epi
    (``etna.core.device.donor_table_spec``: from 0, increasing, covering the epi).
    """
    rec = read_donor_table_json(src)
    table = tuple((float(d), float(n)) for d, n in rec["donor_table"])
    if epi_thickness_um is not None:
        from etna.core.device import donor_table_spec

        donor_table_spec(float(epi_thickness_um) * 1e-4, table)
    return table


def write_json(result, path, **kw) -> Path:
    """``result.to_json()`` to ``path``."""
    p = Path(path)
    p.write_text(result.to_json(**kw))
    return p


def _fmt(v, f=".4g"):
    return "n/a" if v is None else format(v, f)


def report_markdown(
    cv_comparison=None, cv_inversion=None, profile_table=None, it_analysis=None
) -> str:
    """One Markdown report of whatever results are given (numbers, notes, limits,
    provenance)."""
    out = [
        "# ETNA measurement comparison report",
        "",
        f"Generated {now_utc()}.",
        "",
        f"_{MEASUREMENT_LABEL}._",
        "",
    ]

    def prov(d, m=None):
        lines = [
            f"- data: `{d['source_file']}` (sha256 `{d['source_sha256'][:16]}`), {d['kind']}"
        ]
        if m:
            lines.append(
                f"- model: ETNA {m.get('etna_version')}, solver {m.get('solver') or 'n/a'}, "
                f"config `{m.get('config_hash')}`"
            )
        return lines

    def bullets(title, items):
        if not items:
            return []
        return [f"**{title}**", ""] + [f"- {x}" for x in items] + [""]

    if cv_comparison is not None:
        c = cv_comparison
        out += (
            ["## C-V vs model", ""] + prov(c.data_provenance, c.model_provenance) + [""]
        )
        out += [
            "| quantity | value |",
            "| --- | --- |",
            f"| points compared | {c.n_points} |",
            f"| rms W (um) | {_fmt(c.rms_W_um)} |",
            f"| max abs W (um) | {_fmt(c.max_abs_W_um)} |",
            f"| rms C (relative) | {_fmt(c.rms_C_rel)} |",
            f"| chi2 | {_fmt(c.chi2)} |",
            "",
        ]
        out += bullets("Notes", c.notes) + bullets("Limits", c.limits)
    if cv_inversion is not None:
        v = cv_inversion
        out += (
            ["## C-V inversion to a donor_table", ""]
            + prov(v.data_provenance, v.model_provenance)
            + [""]
        )
        out += [
            f"rms W {_fmt(v.rms_W_um)} um, chi2 {_fmt(v.chi2)} / {v.dof} dof, lam {v.lam:g}, "
            f"probed depth {v.probed_depth_um[0]:.3g}..{v.probed_depth_um[1]:.3g} um.",
            "",
            "| depth (um) | N_D (cm^-3) | rel. 1-sigma | constrained | Debye (um) |",
            "| --- | --- | --- | --- | --- |",
        ]
        import numpy as np

        for d, l, s, c_, ld in zip(
            v.knots_um, v.log10_N, v.sigma_log10_N, v.constrained, v.debye_um
        ):
            out.append(
                f"| {d:.4g} | {10**l:.4g} | {np.log(10) * s:.3g} | {'yes' if c_ else 'no'} | {ld:.3g} |"
            )
        out += [""] + bullets("Notes", v.notes) + bullets("Limits", v.limits)
    if profile_table is not None:
        t = profile_table
        out += (
            ["## SIMS/SRP profile resampled to a donor_table", ""]
            + prov(t.data_provenance)
            + [""]
        )
        out += [
            f"Junction depth {t.junction_depth_um:.4g} um ({t.junction_source}); epi "
            f"{t.epi_thickness_um:g} um; p+ plateau N_A {_fmt(t.N_A_plateau, '.3g')} cm^-3.",
            "",
            "| depth (um) | N_D (cm^-3) | rel. 1-sigma | rows | censored | extrapolated |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
        for d, n, u, r, c_, e in zip(
            t.knots_um, t.N_D, t.rel_uncertainty, t.n_rows, t.censored, t.extrapolated
        ):
            out.append(
                f"| {d:.4g} | {n:.4g} | {_fmt(None if u != u else u, '.3g')} | {r} | {'yes' if c_ else ''} | {'yes' if e else ''} |"
            )
        out += [""] + bullets("Notes", t.notes)
    if it_analysis is not None:
        a = it_analysis
        out += (
            ["## I(T): activation energy", ""]
            + prov(a.data_provenance, a.model_provenance or None)
            + [""]
        )
        out += [
            f"Fit: ln(|I|/T^{a.T_power:g}) vs 1/kT.",
            "",
            "| bias (V) | E_a data (eV) | sigma (eV) | E_a model (eV) | data - model (eV) | z | E_g/2 (eV) |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
        comp = {c.bias_V: c for c in a.comparisons}
        for f in a.fits:
            c = comp.get(f.bias_V)
            out.append(
                f"| {_fmt(f.bias_V, 'g')} | {f.E_a_eV:.4f} | {f.sigma_E_a_eV:.4f} | "
                f"{_fmt(c.model_fit.E_a_eV if c else None, '.4f')} | {_fmt(c.delta_eV if c else None, '+.4f')} | "
                f"{_fmt(c.z if c else None, '.2f')} | {_fmt(c.E_g_half_eV if c else None, '.4f')} |"
            )
        out += [""] + bullets("Notes", a.notes) + bullets("Limits", a.limits)
    return "\n".join(out) + "\n"


def results_json(**results) -> str:
    """All given results (name -> result with ``to_dict``) in one JSON document."""
    return dumps({k: v.to_dict() for k, v in results.items() if v is not None})
