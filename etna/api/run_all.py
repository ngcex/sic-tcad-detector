"""One-shot run of every simulation for a single ``DeviceConfig``.

``run_all(config)`` calls the public facades (``run_cv``, ``run_cce``, ...) one
after the other and returns a ``Report`` with one ``SectionResult`` per section:
status, solver, failure reason, wall time and the underlying ``SimResult``.

A section that raises, is outside the model's validity or does not apply to the
configuration (1D-only facades on a 2D device) becomes a labelled row; it never
aborts the run. Sections run in one process: the facades leave no devsim device
alive, which ``run_all`` asserts between sections (a leaked device would make a
later build depend on an earlier one, so it is reported as an error).
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from typing import Callable

import numpy as np

from etna.api.device import DeviceConfig
from etna.api.results import SimResult, peak_field_1d

_ONE_D_ONLY = "needs a 1D device (set half_width_um to None)"


def _sections() -> "dict[str, tuple[str, Callable, bool]]":
    """name -> (title, runner, one_d_only). Imported lazily to keep import light."""
    from etna.api import simulation as s

    return {
        "cv": ("C-V", lambda c: s.run_cv(c), True),
        "cce": ("Charge collection", lambda c: s.run_cce(c), True),
        "field": ("Field map", lambda c: s.run_field(c), False),
        "damage": ("Radiation damage", lambda c: s.run_radiation_damage(c), True),
        "dark": ("Dark current", lambda c: s.run_dark_current(c), True),
        "temperature": (
            "Temperature sweep",
            lambda c: s.run_temperature_sweep(c),
            True,
        ),
        "flash": ("Flash recombination", lambda c: s.run_flash_recombination(c), True),
        "sweep": (SWEEP_TITLE, _run_sweep, False),
    }


# The default graded n-/n+ transition (Lb, N_end) was fitted on the C-V of the 10 um epi;
# the sweep applies the config's transition unchanged at any epi thickness.
_TRANSITION_FIT_EPI_UM = 10.0
TRANSITION_EXTRAPOLATION_NOTE = (
    "transition Lb/N_end fitted on the 10 um epi, applied unchanged at this epi "
    "thickness: model extrapolation"
)
_SWEEP_HALF_WIDTHS = (25.0, 50.0)
_SWEEP_BIASES = (20.0, 40.0, 60.0)
SWEEP_TITLE = "Microdosimetry sweep (2D engine)"


def _run_sweep(config: DeviceConfig) -> SimResult:
    """Small operating-point grid around the sidebar epi/doping, engine only.

    Never falls back to the 2D devsim path (tens of seconds per config): if any
    grid point is outside the 2D engine box the section raises, and is reported
    as skipped by ``run_all``.
    """
    from etna.core.engine2d import engine2d_unsupported_reason
    from etna.core.optimization import microdosimetric_sweep

    epi_cm = config.epi_thickness_um * 1e-4
    # the p+ anode, graded doping and n+ layer (if any) go to the sweep; the engine serves
    # them inside its parity box
    # nplus_thickness_cm always explicit: None (legacy config) must not fall back to the
    # sweep's default n+ layer
    nplus = dict(nplus_thickness_cm=None)
    if config.nplus_thickness_um is not None:
        nplus = dict(
            nplus_thickness_cm=config.nplus_thickness_um * 1e-4,
            N_D_nplus=config.N_D_nplus,
        )
    nplus.update(
        substrate_thickness_cm=config.substrate_thickness_um * 1e-4,
        N_D_junction=config.N_D_junction,
        L_transition=config.L_transition_um * 1e-4,
    )
    # graded n-/n+ transition always explicit too: None (abrupt epi end, e.g. legacy or
    # nplus_calibrated configs) must not fall back to the builders' default transition
    nplus.update(
        buffer_thickness_cm=(
            None
            if config.buffer_thickness_um is None
            else config.buffer_thickness_um * 1e-4
        ),
        buffer_end_doping=config.buffer_end_doping,
    )
    # tabulated epi profile (opt-in): the table replaces the graded doping of every
    # config; the sweep's N_D_bulk column is then only the config's informative value
    if config.donor_table is not None:
        nplus.update(donor_table=config.donor_table)
    for hw in _SWEEP_HALF_WIDTHS:
        for vb in _SWEEP_BIASES:
            # the sweep builds graded, N_A 1e19, 300 K devices: the gate sees the config's
            # values so a config the sweep would not represent is reported, not re-labelled
            why = engine2d_unsupported_reason(
                hw,
                vb,
                epi_cm,
                config.N_D_bulk,
                10,
                doping_profile=config.doping_profile,
                N_A=config.N_A,
                T=config.T,
                **nplus,
            )
            if why:
                raise _SweepUnsupported(why)
    df = microdosimetric_sweep(
        half_widths_um=_SWEEP_HALF_WIDTHS,
        epi_thicknesses_cm=(epi_cm,),
        N_D_bulks=(config.N_D_bulk,),
        V_biases=_SWEEP_BIASES,
        solver="auto",
        **nplus,
    )
    n_valid = int(df["is_valid"].sum())
    best = df.iloc[0]
    summary = f"{n_valid}/{len(df)} valid"
    if config.nplus_thickness_um is not None:
        summary += (
            f" (n+ {config.nplus_thickness_um:g} um, {config.N_D_nplus:.1e} cm^-3)"
        )
    else:
        summary += " (legacy model: no n+ layer)"
    if config.buffer_thickness_um is not None:
        summary += (
            f" (graded n-/n+ transition {config.buffer_thickness_um:g} um to "
            f"{config.buffer_end_doping:.1e} cm^-3)"
        )
        if config.epi_thickness_um != _TRANSITION_FIT_EPI_UM:
            summary += f" [{TRANSITION_EXTRAPOLATION_NOTE}]"
    if config.donor_table is not None:
        summary += (
            f" (tabulated epi profile, {len(config.donor_table)} knots: the table defines "
            "the doping of every row; the N_D_bulk column is the config's informative value)"
        )
    if best["is_valid"]:
        summary += (
            f"; best: half-width {best['half_width_um']:g} um at "
            f"{best['V_bias']:g} V, CCE {best['center_cce']:.4f}"
        )
    per_v = df.groupby("V_bias")["center_cce"].max()
    return SimResult(
        config=config,
        sim_type="sweep",
        x=per_v.index.to_numpy(dtype=float),
        y=per_v.to_numpy(dtype=float),
        metadata={
            "solver": "+".join(sorted(set(df["solver"]))),
            "summary": summary,
            "table": df,
        },
    )


class _SweepUnsupported(Exception):
    pass


SECTION_NAMES = tuple(_sections())
DEFAULT_SECTIONS = tuple(n for n in SECTION_NAMES if n != "sweep")


@dataclass
class SectionResult:
    name: str
    title: str
    status: str  # "ok" | "warning" (ran, no valid point) | "error" | "skipped"
    solver: str = ""
    reason: str = ""  # failure / skip reason, empty when ok
    wall_s: float = 0.0
    summary: str = ""  # one-line key figures
    result: "SimResult | None" = field(default=None, repr=False)


@dataclass
class Report:
    config: DeviceConfig
    sections: "list[SectionResult]"
    wall_s: float

    @property
    def ok(self) -> bool:
        """True when no section failed or warned (skipped sections do not count)."""
        return all(s.status in ("ok", "skipped") for s in self.sections)

    def get(self, name: str) -> SectionResult:
        for s in self.sections:
            if s.name == name:
                return s
        raise KeyError(name)

    def rows(self) -> "list[dict]":
        return [
            {
                "section": s.title,
                "status": s.status,
                "solver": s.solver,
                "key figures": s.summary,
                "reason": s.reason,
                "wall time (s)": round(s.wall_s, 2),
            }
            for s in self.sections
        ]

    def to_json(self) -> str:
        """Summary table plus the primary curve of every section that ran."""
        out = {"config": asdict(self.config), "wall_s": self.wall_s, "sections": []}
        for s in self.sections:
            d = {k: v for k, v in asdict_light(s).items()}
            if s.result is not None and len(np.atleast_1d(s.result.x)):
                d["x"] = np.asarray(s.result.x, dtype=float).tolist()
                d["y"] = np.asarray(s.result.y, dtype=float).tolist()
            out["sections"].append(d)
        return json.dumps(out, indent=2, default=_json_default)


def asdict_light(s: SectionResult) -> dict:
    return {
        "name": s.name,
        "title": s.title,
        "status": s.status,
        "solver": s.solver,
        "reason": s.reason,
        "wall_s": s.wall_s,
        "summary": s.summary,
    }


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def _truncation_note(res: SimResult) -> str:
    """' (TRUNCATED: reached X of requested Y V)' for a sweep cut short, else ''."""
    md = res.metadata
    if not md.get("truncated"):
        return ""
    req = md.get("requested_v_stop")
    x = np.asarray(res.x, dtype=float)
    reached = f"{x[np.argmax(np.abs(x))]:.3g} V" if x.size else "no point"
    req_txt = f" of requested {req:.3g} V" if req is not None else ""
    return f" (TRUNCATED: reached {reached}{req_txt})"


def _valid_counts(res: SimResult) -> "tuple[int, int]":
    """(valid, total) points of a curve result: metadata["valid"] if present, else finite y."""
    y = np.asarray(res.y, dtype=float)
    if not y.size:
        return 0, 0
    valid = res.metadata.get("valid")
    valid = None if valid is None else np.asarray(valid, dtype=bool)
    if valid is None or valid.shape != y.shape:  # absent or not per point: use finite y
        valid = np.isfinite(y)
    return int(np.count_nonzero(valid)), int(y.size)


def _compensation_note(res: SimResult) -> str:
    """'; full donor compensation at X p/cm^2' for a damage sweep (else '')."""
    phi_c = res.metadata.get("compensation_fluence")
    if res.sim_type != "damage" or phi_c is None:
        return ""
    return f"; full donor compensation at {phi_c:.3g} p/cm^2"


def _summarize(res: SimResult) -> str:
    if res.metadata.get("summary"):
        return res.metadata["summary"] + _truncation_note(res)
    x, y = np.asarray(res.x, dtype=float), np.asarray(res.y, dtype=float)
    cfg = res.config
    if res.sim_type == "field" and x.size:
        # exclude the nm-scale mesh/estimator-dependent spikes at the abrupt doping steps
        pk = peak_field_1d(x, y, cfg)
        if pk is not None and pk[2]:
            return (
                f"{x.size} pts, peak |E| {pk[0]:.3g} V/cm at x={pk[1]:.2f} um "
                "(excl. spikes at the abrupt doping steps)"
            )
    if not x.size and res.mesh is not None and "ElectricField" in res.mesh.node_values:
        E = np.asarray(res.mesh.node_values["ElectricField"], dtype=float)
        return f"2D map, {E.size} nodes, max |E| {np.nanmax(E):.3g} V/cm"
    if not x.size or not y.size:
        return "no data" + _truncation_note(res)
    finite = y[np.isfinite(y)]
    n_valid, n_pts = _valid_counts(res)
    partial = f"{n_valid} of {n_pts} points valid; " if n_valid < n_pts else ""
    if not finite.size:
        return (
            f"{x.size} points, all NaN"
            + _compensation_note(res)
            + _truncation_note(res)
        )
    return (
        (
            f"{partial}{x.size} pts, x {np.nanmin(x):.3g}..{np.nanmax(x):.3g}, "
            f"y {finite.min():.3g}..{finite.max():.3g}"
        )
        + _compensation_note(res)
        + _truncation_note(res)
    )


def _solver_of(res: SimResult) -> str:
    return str(res.metadata.get("solver") or "")


def _delete_device(name: str) -> None:
    from etna.core.devsim_reset import delete_device_and_mesh

    delete_device_and_mesh(name)


def _devices_alive() -> "list[str]":
    import devsim

    return list(devsim.get_device_list())


def run_all(config: DeviceConfig, sections: "tuple[str, ...] | None" = None) -> Report:
    """Run the selected sections (default: all except the opt-in ``sweep``) for ``config``.

    Never raises for a section failure; see ``Report.rows()``.
    """
    table = _sections()
    names = tuple(sections) if sections is not None else DEFAULT_SECTIONS
    unknown = [n for n in names if n not in table]
    if unknown:
        raise ValueError(f"unknown section(s) {unknown}; choose from {list(table)}")

    t_all = time.perf_counter()
    preexisting = set(_devices_alive())  # not ours; never flagged or deleted
    done: "list[SectionResult]" = []
    for name in names:
        title, runner, one_d_only = table[name]
        if one_d_only and config.half_width_um is not None:
            done.append(SectionResult(name, title, "skipped", reason=_ONE_D_ONLY))
            continue
        t0 = time.perf_counter()
        try:
            res = runner(config)
        except _SweepUnsupported as exc:
            done.append(
                SectionResult(
                    name, title, "skipped", reason=f"not available here: {exc}"
                )
            )
        except Exception as exc:  # a section never aborts the run
            done.append(
                SectionResult(
                    name,
                    title,
                    "error",
                    reason=f"{type(exc).__name__}: {exc}",
                    wall_s=time.perf_counter() - t0,
                )
            )
        else:
            err = res.metadata.get("error")
            n_valid, n_pts = _valid_counts(res)
            status, why = ("error", str(err)) if err else ("ok", "")
            if not err and n_pts and n_valid == 0:
                status, why = "warning", f"no valid point (0 of {n_pts})" + (
                    f": {res.metadata['reason'][-1]}"
                    if isinstance(res.metadata.get("reason"), (list, tuple))
                    and res.metadata["reason"]
                    else ""
                )
            sec = SectionResult(
                name,
                title,
                status,
                solver=_solver_of(res),
                reason=why,
                wall_s=time.perf_counter() - t0,
                summary="" if err else _summarize(res),
                result=res,
            )
            done.append(sec)
        leaked = [d for d in _devices_alive() if d not in preexisting]
        if leaked:  # would couple this section to the next one
            last = done[-1]
            last.status = "error"
            last.reason = (last.reason + "; " if last.reason else "") + (
                f"leaked devsim device(s) {leaked}"
            )
            for d in leaked:
                _delete_device(d)
    return Report(config, done, time.perf_counter() - t_all)
