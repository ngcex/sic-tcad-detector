"""Device spec sheet from a DeviceConfig and an ``etna.run_all`` Report.

``build_spec(config, report)`` assembles a plain dict (layers incl. the graded profile,
area, bias, full-depletion numbers with their definitions, CCE, dark current, solver
provenance and validity per section, known systematics). The serialisers turn it into
JSON, CSV (key/value rows), Markdown and a PDF (matplotlib; already a dependency).
Pure: never runs a simulation; sections missing from the report become labelled rows.
"""

from __future__ import annotations

import csv
import io
import json
from dataclasses import asdict

import numpy as np

from etna._version import __version__
from etna.api.results import peak_field_1d
from etna.export.cross_section import stack_layers
from etna.export.labels import (
    ALPHA_NOTE,
    CALIBRATION_NOTE,
    CV_SYSTEMATIC_NOTE,
    CV_SYSTEMATIC_OTHER,
    DARK_CALIBRATION_NOTE,
    MEASURED_W_30V_UM,
    NPLUS_CALIBRATION_NOTE,
    CV_SYSTEMATIC_NPLUS_NOTE,
    CV_SYSTEMATIC_TABULATED_NOTE,
    TABULATED_CALIBRATION_NOTE,
    VALIDATION_LABEL,
    VFD_DEFINITION,
    analytic_vfd,
    doping_formula,
    epi_donor_profile,
    has_buffer,
    has_nplus,
    has_table,
    is_measured_cv_device,
    is_measured_cv_tabulated_device,
    is_nplus_calibrated_device,
    is_reference_device,
    reach_through_definition,
    table_sensitive_end_um,
)

# The run_all sections a spec sheet reads (the cheap ones; see etna.run_all).
SPEC_SECTIONS = ("cv", "cce", "field", "dark")
_DOPING_DEPTHS_UM = (0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0)
_DARK_BIAS_V = -60.0


def _section(report, name):
    if report is None:
        return None
    try:
        return report.get(name)
    except KeyError:
        return None


def _ok_result(report, name):
    sec = _section(report, name)
    if sec is None or sec.status != "ok" or sec.result is None:
        return None
    return sec.result


def _interp_at(x, y, v):
    """y at bias v (V, negative = reverse) if v lies inside the swept range, else None."""
    x = np.abs(np.asarray(x, dtype=float))
    y = np.asarray(y, dtype=float)
    if not x.size or not (x.min() - 1e-9 <= abs(v) <= x.max() + 1e-9):
        return None
    o = np.argsort(x)
    return float(np.interp(abs(v), x[o], y[o]))


def _f(v):
    return None if v is None else float(v)


def _provenance_rows(report) -> "list[dict]":
    rows = []
    for name in SPEC_SECTIONS:
        sec = _section(report, name)
        if sec is None:
            rows.append(
                dict(
                    section=name,
                    status="not run",
                    solver="",
                    validated_domain=None,
                    domain_note="",
                    reason="section not in the report",
                )
            )
            continue
        md = sec.result.metadata if sec.result is not None else {}
        rows.append(
            dict(
                section=name,
                status=sec.status,
                solver=sec.solver,
                validated_domain=md.get("validated_domain"),
                domain_note=md.get("domain_note") or "",
                reason=sec.reason or md.get("fallback_reason", "") or "",
            )
        )
    return rows


def build_spec(config, report=None) -> dict:
    """Spec-sheet dict for ``config`` from ``report`` (an ``etna.run_all`` Report).

    ``report`` may be None or partial (2D device, skipped / failed sections): the
    corresponding entries are None with a reason; nothing is recomputed here.
    """
    if report is not None and report.config != config:
        raise ValueError("report was computed for a different DeviceConfig")
    nplus = has_nplus(config)
    ref = is_reference_device(config)
    measured_cv = is_measured_cv_device(config)  # == ref since the 2026-10-18 flip
    prev_default = is_nplus_calibrated_device(config)
    table = has_table(config)
    tab_preset = is_measured_cv_tabulated_device(config)
    layers = stack_layers(config)
    t_e = config.epi_thickness_um
    depths = [d for d in _DOPING_DEPTHS_UM if d < t_e]
    if has_buffer(config):  # sample the transition too (start, middle)
        t0 = t_e - config.buffer_thickness_um
        depths = sorted(set(depths) | {round(t0, 6), round((t0 + t_e) / 2, 6)})
    if table:  # the knots themselves, plus the end of the low-doped epi
        t0 = table_sensitive_end_um(config)
        depths = sorted(
            set(depths) | {round(t0, 6)} | {round(d, 6) for d, _ in config.donor_table}
        )
    depths = [d for d in depths if d < t_e] + [t_e]
    profile = [
        dict(
            depth_below_junction_um=float(d),
            N_D_cm3=float(epi_donor_profile(config, d)),
        )
        for d in depths
    ]
    if table:
        model = (
            "opt-in tabulated measured-C-V profile (19 knots: inversion of the digitised "
            "Fig. 4b C-V, not a SIMS profile)"
            if tab_preset
            else f"tabulated epi profile ({len(config.donor_table)} knots, user table)"
        ) + (" + n+ layer" if nplus else ", no n+ layer")
    elif has_buffer(config):
        model = (
            "default stack: measured-C-V profile (graded n-/n+ transition, fit to the "
            "digitised Fig. 4b)"
            if ref
            else "stack with a graded n-/n+ transition (buffer)"
        ) + (" + n+ layer" if nplus else ", no n+ layer")
    elif nplus:
        model = (
            "previous default (3-point n+ calibration, abrupt epi/n+ step)"
            if prev_default
            else "n+ stack with an abrupt epi/n+ step"
        )
    else:
        model = "legacy model (no n+, cathode on the n- epi)"

    spec = {
        "label": VALIDATION_LABEL,
        "generator": f"etna {__version__}",
        "model": model,
        "reference_device": ref,
        "measured_cv_preset": measured_cv,
        "measured_cv_tabulated_preset": tab_preset,
        "nplus_calibrated_preset": prev_default,
        "dimension": "1D" if config.half_width_um is None else "2D",
        "config": asdict(config),
        "area_cm2": config.area_cm2,
        "area_mm2": config.area_cm2 * 100.0,
        "temperature_K": config.T,
        "layers": [
            {
                k: L[k]
                for k in (
                    "name",
                    "top_um",
                    "bottom_um",
                    "thickness_um",
                    "doping",
                    "note",
                )
            }
            for L in layers
        ],
        "epi_doping": {
            "profile": config.doping_profile,
            "formula": doping_formula(config),
            "N_D_junction_cm3": config.N_D_junction,
            "N_D_bulk_cm3": config.N_D_bulk,
            "L_transition_um": config.L_transition_um,
            "buffer_thickness_um": config.buffer_thickness_um,
            "buffer_end_doping_cm3": config.buffer_end_doping,
            # tabulated profile: the knots as given (depth_um, N_D) and where the low-doped
            # epi ends by the decade rule (None without a table)
            "donor_table": (
                None
                if not table
                else [
                    dict(depth_below_junction_um=float(d), N_D_cm3=float(n))
                    for d, n in config.donor_table
                ]
            ),
            "table_sensitive_epi_end_um": (
                None if not table else float(table_sensitive_end_um(config))
            ),
            "table_note": (
                TABULATED_CALIBRATION_NOTE
                if tab_preset
                else (
                    "user-supplied table: provenance not known to etna"
                    if table
                    else None
                )
            ),
            "samples": profile,
        },
    }

    # full depletion: both numbers, each with its definition
    cv = _ok_result(report, "cv")
    cce = _ok_result(report, "cce")
    rt = None
    rt_method = ""
    for r in (cv, cce):
        if r is not None and r.metadata.get(
            "reach_through_known", "reach_through_V" in r.metadata
        ):
            rt = r.metadata.get("reach_through_V")
            rt_method = r.metadata.get("reach_through_method") or ""
            break
    spec["full_depletion"] = {
        "reach_through_V": _f(rt),
        "reach_through_definition": reach_through_definition(config),
        "reach_through_method": rt_method,
        "reach_through_note": (
            "" if (cv or cce) else "C-V / CCE not available in the report"
        ),
        "analytic_V_fd_V": analytic_vfd(config),
        "analytic_V_fd_definition": VFD_DEFINITION,
    }

    # C-V: W at the measurement biases. W_um is eps*A/C of the small-signal C (y), the
    # quantity a measured W = eps*A/C is; W_edge_um is the n/N_D = 0.5 depletion edge.
    cvd = {"available": cv is not None}
    if cv is not None:
        from etna.core.batched.material import EPS

        W = EPS * config.area_cm2 / np.asarray(cv.y, dtype=float) * 1e4
        We = np.asarray(cv.metadata.get("depletion_widths", []), dtype=float) * 1e4
        at = lambda arr: {
            f"{v:g} V": _f(_interp_at(cv.x, arr, v)) for v in (0.0, -10.0, -30.0)
        }
        cvd.update(
            bias_range_V=[float(np.max(cv.x)), float(np.min(cv.x))],
            W_um=at(W),
            W_definition="eps*A/C, C = small-signal dQ/dV (as the measured W)",
            W_edge_um=at(We) if We.size == np.size(cv.x) else {},
            C_at_minus30V_F=_f(_interp_at(cv.x, cv.y, -30.0)),
        )
    w30 = cvd.get("W_um", {}).get("-30 V")
    d30 = None if w30 is None else 100.0 * (w30 - MEASURED_W_30V_UM) / MEASURED_W_30V_UM
    if ref:
        cvd["systematic"] = (
            CV_SYSTEMATIC_NOTE.format(w=w30, d=d30)
            if w30 is not None
            else "No out-of-fit C-V comparison: the whole digitised curve is a fit input "
            "(see cv.calibration)."
        )
        cvd["measured_W_minus30V_um"] = MEASURED_W_30V_UM
    elif tab_preset:
        cvd["systematic"] = (
            CV_SYSTEMATIC_TABULATED_NOTE.format(w=w30, d=d30)
            if w30 is not None
            else "No out-of-fit C-V comparison: the whole digitised curve is the "
            "inversion's input (see cv.calibration)."
        )
        cvd["measured_W_minus30V_um"] = MEASURED_W_30V_UM
    elif prev_default and w30 is not None:
        cvd["systematic"] = CV_SYSTEMATIC_NPLUS_NOTE.format(w=w30, d=d30)
        cvd["measured_W_minus30V_um"] = MEASURED_W_30V_UM
    else:
        cvd["systematic"] = CV_SYSTEMATIC_OTHER
    cvd["calibration"] = (
        CALIBRATION_NOTE
        if ref
        else (
            TABULATED_CALIBRATION_NOTE
            if tab_preset
            else (
                NPLUS_CALIBRATION_NOTE
                if prev_default
                else "The C-V fit applies to the reference device (default stack) only."
            )
        )
    )
    spec["cv"] = cvd

    ccd = {"available": cce is not None}
    if cce is not None:
        y = np.asarray(cce.y, dtype=float)
        ccd.update(
            bias_range_V=[float(np.max(cce.x)), float(np.min(cce.x))],
            CCE_min=float(np.nanmin(y)),
            CCE_max=float(np.nanmax(y)),
            CCE_at_deepest_bias=float(y[np.argmax(np.abs(cce.x))]),
        )
    ccd["generation"] = (
        ALPHA_NOTE
        if nplus
        else ("Alpha generation confined to the n- epi (legacy model, no n+ layer).")
    )
    spec["cce"] = ccd

    dark = _ok_result(report, "dark")
    from etna.core.device import (
        DEFAULT_STACK,
        MEASURED_CV_TABULATED_STACK,
        NPLUS_CALIBRATION,
    )

    n_t = (
        float(dark.metadata.get("N_t", config.N_t)) if dark is not None else config.N_t
    )
    calibrated = (
        ref and n_t == DEFAULT_STACK["N_t"] and config.E_t == DEFAULT_STACK["E_t"]
    )
    calibrated_prev = (
        prev_default
        and n_t == NPLUS_CALIBRATION["N_t"]
        and config.E_t == NPLUS_CALIBRATION["E_t"]
    )
    calibrated_tab = (
        tab_preset
        and n_t == MEASURED_CV_TABULATED_STACK["N_t"]
        and config.E_t == MEASURED_CV_TABULATED_STACK["E_t"]
    )
    dd = {
        "available": dark is not None,
        "N_t_used": n_t,
        "E_t_eV": config.E_t,
        "note": (
            DARK_CALIBRATION_NOTE
            if calibrated
            else (
                (
                    DARK_CALIBRATION_NOTE
                    + f" (N_t = {n_t:.7g} fitted on the previous default profile.)"
                )
                if calibrated_prev
                else (
                    (
                        DARK_CALIBRATION_NOTE
                        + f" (N_t = {n_t:.7g} refitted on the tabulated profile; the "
                        "parametric default's 2.849024e12 gives 17.999 pA there: the proxy is "
                        "profile-independent to 1e-4.)"
                    )
                    if calibrated_tab
                    else (
                        f"N_t = {n_t:.4g}, E_t = {config.E_t:g} eV on this device: the dark-current "
                        "calibration (default stack, N_t = "
                        f"{DEFAULT_STACK['N_t']:.7g}, E_t = {DEFAULT_STACK['E_t']:g} eV, fitted to the "
                        "measured < 18 pA bound at -60 V on 4 mm2) does not apply; uncalibrated "
                        "estimate."
                    )
                )
            )
        ),
    }
    if dark is not None:
        I60 = _interp_at(dark.x, np.abs(dark.y), _DARK_BIAS_V)
        dd.update(
            bias_range_V=[float(np.max(dark.x)), float(np.min(dark.x))],
            area_cm2=config.area_cm2,
            I_at_minus60V_A=_f(I60),
            J_at_minus60V_A_per_cm2=None if I60 is None else I60 / config.area_cm2,
            I_at_deepest_bias_A=float(np.abs(dark.y)[np.argmax(np.abs(dark.x))]),
            truncated=bool(dark.metadata.get("truncated", False)),
        )
    spec["dark_current"] = dd

    field = _ok_result(report, "field")
    fd = {"available": field is not None}
    if field is not None:
        fd["bias_V"] = _f(field.metadata.get("bias_V"))
        pk = peak_field_1d(field.x, field.y, config) if np.size(field.x) else None
        if pk is not None:
            fd.update(
                peak_E_V_per_cm=pk[0],
                peak_depth_um=pk[1],
                peak_note=(
                    "excl. the nm-scale spikes at the abrupt doping steps (mesh/estimator dependent)"
                    if pk[2]
                    else ""
                ),
            )
        elif field.mesh is not None and "ElectricField" in field.mesh.node_values:
            E = np.asarray(field.mesh.node_values["ElectricField"], dtype=float)
            fd.update(
                peak_E_V_per_cm=float(np.nanmax(E)),
                peak_note="2D map, max over all nodes",
            )
    spec["field"] = fd

    spec["provenance"] = _provenance_rows(report)
    spec["statements"] = [
        VALIDATION_LABEL,
        "Measurement comparisons (the 'stated' exceptions): "
        + (
            "the full digitised Photons Fig. 4b C-V is a FIT INPUT of the default stack "
            "(see cv.calibration: plateau rms 0.013 um, knee rms 0.24 um), not a "
            "validation; dark current is a calibration to the measured bound."
            if ref
            else (
                (
                    "the full digitised Photons Fig. 4b C-V is the INPUT of the tabulated "
                    "profile's inversion (see cv.calibration: knee rms 0.014 um, plateau rms "
                    "0.007 um), not a validation, and the table is not a SIMS profile; dark "
                    "current is a calibration to the measured bound."
                )
                if tab_preset
                else (
                    "C-V eps*A/C at -30 V vs Photons Fig. 4b (see cv.systematic; 0 / -10 V are "
                    "calibration inputs of this previous default); dark current is a "
                    "calibration to the measured bound."
                    if prev_default
                    else "none for this config (they apply to the reference device only)."
                )
            )
        ),
        "No fabrication layout: the GDSII export is a labelled REFERENCE-NOT-FAB schematic; "
        "guard ring / pad dimensions are not available and are not invented.",
    ]
    return spec


# --------------------------------------------------------------------------- formats


def spec_to_json(spec: dict) -> str:
    return json.dumps(
        spec, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)
    )


def _flatten(prefix, obj, out):
    if isinstance(obj, dict):
        for k, v in obj.items():
            _flatten(f"{prefix}.{k}" if prefix else str(k), v, out)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            _flatten(f"{prefix}[{i}]", v, out)
    else:
        out.append((prefix, "" if obj is None else obj))


def spec_to_csv(spec: dict) -> str:
    """Key/value CSV; the first row is the validation label."""
    rows = []
    _flatten("", spec, rows)
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["# " + VALIDATION_LABEL])
    w.writerow(["key", "value"])
    for k, v in rows:
        w.writerow([k, v])
    return buf.getvalue()


def _g(v, fmt="{:.4g}", unit=""):
    if v is None:
        return "n/a"
    return fmt.format(v) + (f" {unit}" if unit else "")


def spec_to_markdown(spec: dict) -> str:
    L = []
    L.append("# Device spec sheet")
    L.append("")
    L.append(f"> **{spec['label']}**")
    L.append("")
    L.append(
        f"- Generator: {spec['generator']}; model: {spec['model']}; {spec['dimension']}"
    )
    L.append(
        f"- Active area: {spec['area_cm2']:g} cm^2 ({spec['area_mm2']:g} mm^2); T = {spec['temperature_K']:g} K"
    )
    L.append("")
    L.append("## Layers (anode to cathode)")
    L.append("")
    L.append("| layer | depth (um) | thickness (um) | doping | note |")
    L.append("|---|---|---|---|---|")
    for x in spec["layers"]:
        L.append(
            f"| {x['name']} | {x['top_um']:g} - {x['bottom_um']:g} | {x['thickness_um']:g} | {x['doping']} | {x['note']} |"
        )
    ed = spec["epi_doping"]
    L.append("")
    L.append(f"Epi doping ({ed['profile']}): `{ed['formula']}`")
    if ed.get("buffer_thickness_um") is not None:
        L.append("")
        L.append(
            f"Graded n-/n+ transition: last {ed['buffer_thickness_um']:g} um of the epi, "
            f"N_D -> {ed['buffer_end_doping_cm3']:.3g} cm^-3 at the epi end (fitted form, "
            "not a measured profile)."
        )
    if ed.get("donor_table"):
        L.append("")
        L.append(
            f"Tabulated profile: {len(ed['donor_table'])} knots (log10 N_D linear between "
            f"them); low-doped epi ends at {ed['table_sensitive_epi_end_um']:.2f} um below "
            "the junction (one decade above the table minimum). "
            + (ed.get("table_note") or "")
        )
        L.append("")
        L.append("| knot | depth below junction (um) | N_D (cm^-3) | log10 N_D |")
        L.append("|---|---|---|---|")
        for k, s in enumerate(ed["donor_table"]):
            L.append(
                f"| {k} | {s['depth_below_junction_um']:g} | {s['N_D_cm3']:.4g} | "
                f"{np.log10(s['N_D_cm3']):.3f} |"
            )
    L.append("")
    L.append("| depth below junction (um) | N_D (cm^-3) |")
    L.append("|---|---|")
    for s in ed["samples"]:
        L.append(f"| {s['depth_below_junction_um']:g} | {s['N_D_cm3']:.4g} |")
    fd = spec["full_depletion"]
    L.append("")
    L.append("## Full depletion")
    L.append("")
    L.append(
        f"- **Reach-through: {_g(fd['reach_through_V'], '{:.2f}', 'V')}** "
        f"({fd['reach_through_method'] or fd['reach_through_note'] or 'n/a'}). "
        f"Definition: {fd['reach_through_definition']}."
    )
    L.append(
        f"- **Analytic graded V_fd: {_g(fd['analytic_V_fd_V'], '{:.2f}', 'V')}**. "
        f"Definition: {fd['analytic_V_fd_definition']}."
    )
    cv = spec["cv"]
    L.append("")
    L.append("## C-V")
    L.append("")
    if cv["available"]:
        L.append(
            f"- Bias range: {cv['bias_range_V'][0]:g} to {cv['bias_range_V'][1]:g} V"
        )
        L.append(
            "- W = eps*A/C (small-signal C, as measured): "
            + ", ".join(f"{k}: {_g(v, '{:.3f}', 'um')}" for k, v in cv["W_um"].items())
        )
        if cv.get("W_edge_um"):
            L.append(
                "- depletion edge (n/N_D = 0.5): "
                + ", ".join(
                    f"{k}: {_g(v, '{:.3f}', 'um')}" for k, v in cv["W_edge_um"].items()
                )
            )
        L.append(f"- C(-30 V): {_g(cv['C_at_minus30V_F'], '{:.4g}', 'F')}")
    else:
        L.append("- not available in the report")
    L.append(f"- {cv['calibration']}")
    L.append(f"- {cv['systematic']}")
    cc = spec["cce"]
    L.append("")
    L.append("## Charge collection (CCE)")
    L.append("")
    if cc["available"]:
        L.append(
            f"- Bias range: {cc['bias_range_V'][0]:g} to {cc['bias_range_V'][1]:g} V; "
            f"CCE {cc['CCE_min']:.4f} - {cc['CCE_max']:.4f}; at deepest bias {cc['CCE_at_deepest_bias']:.4f}"
        )
    else:
        L.append("- not available in the report")
    L.append(f"- {cc['generation']}")
    dk = spec["dark_current"]
    L.append("")
    L.append("## Dark current")
    L.append("")
    if dk["available"]:
        L.append(
            f"- I(-60 V) = {_g(dk['I_at_minus60V_A'], '{:.3g}', 'A')} on {dk['area_cm2']:g} cm^2 "
            f"(J = {_g(dk['J_at_minus60V_A_per_cm2'], '{:.3g}', 'A/cm^2')}); "
            f"bias range {dk['bias_range_V'][0]:g} to {dk['bias_range_V'][1]:g} V"
        )
    else:
        L.append("- not available in the report")
    L.append(f"- {dk['note']}")
    fl = spec["field"]
    L.append("")
    L.append("## Electric field")
    L.append("")
    if fl["available"] and fl.get("peak_E_V_per_cm") is not None:
        where = f" at {fl['peak_depth_um']:.2f} um" if "peak_depth_um" in fl else ""
        L.append(
            f"- Peak |E| at {fl['bias_V']:g} V: {fl['peak_E_V_per_cm']:.4g} V/cm{where} {fl.get('peak_note', '')}"
        )
    else:
        L.append("- not available in the report")
    L.append("")
    L.append("## Solver provenance and validity")
    L.append("")
    L.append(
        "| section | status | solver | inside engine/devsim cross-check box "
        "| domain note | reason |"
    )
    L.append("|---|---|---|---|---|---|")
    for p in spec["provenance"]:
        L.append(
            f"| {p['section']} | {p['status']} | {p['solver']} | {p['validated_domain']} | {p['domain_note']} | {p['reason']} |"
        )
    L.append("")
    L.append(
        "Inside the cross-check box (JSON key `validated_domain`) the batched / 2D "
        "engines and devsim agree; it is a solver cross-check, not a validation "
        "against measurement."
    )
    L.append("")
    L.append("## Statements")
    L.append("")
    for s in spec["statements"]:
        L.append(f"- {s}")
    L.append("")
    return "\n".join(L)


def spec_to_pdf(spec: dict, config=None) -> bytes:
    """PDF bytes: the Markdown text on text pages, plus the cross-section if ``config``."""
    import textwrap

    from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
    from matplotlib.figure import Figure

    lines = []
    for raw in spec_to_markdown(spec).splitlines():
        lines.extend(textwrap.wrap(raw, 110) or [""])
    per_page = 62
    buf = io.BytesIO()
    # pages are laid out by hand: layout="none" opts out of the rcParams constrained
    # layout that etna.core.plotting sets
    with PdfPages(buf) as pdf:
        for start in range(0, len(lines), per_page):
            fig = Figure(figsize=(8.27, 11.69), layout="none")
            FigureCanvasPdf(fig)
            fig.text(
                0.05,
                0.97,
                VALIDATION_LABEL,
                fontsize=7,
                color="#b2182b",
                style="italic",
            )
            fig.text(
                0.05,
                0.95,
                "\n".join(lines[start : start + per_page]),
                va="top",
                family="monospace",
                fontsize=6.6,
            )
            pdf.savefig(fig)
        if config is not None:
            from etna.export.cross_section import _draw_cross_section

            fig = Figure(figsize=(11.69, 8.27), layout="none")
            FigureCanvasPdf(fig)
            _draw_cross_section(fig, config)
            pdf.savefig(fig)
    return buf.getvalue()
