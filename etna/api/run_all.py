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

import hashlib
import json
import time
from dataclasses import asdict, dataclass, field
from typing import Callable

import numpy as np

from etna.api.device import DeviceConfig
from etna.api.results import SimResult

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
    }


SECTION_NAMES = tuple(_sections())


@dataclass
class SectionResult:
    name: str
    title: str
    status: str  # "ok" | "error" | "skipped"
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
        return all(s.status == "ok" for s in self.sections)

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


def config_key(config: DeviceConfig) -> str:
    """Stable hash of a DeviceConfig, for caching a Report."""
    return hashlib.sha256(
        json.dumps(asdict(config), sort_keys=True, default=str).encode()
    ).hexdigest()[:16]


def _summarize(res: SimResult) -> str:
    x, y = np.asarray(res.x, dtype=float), np.asarray(res.y, dtype=float)
    if not x.size or not y.size:
        return "no data"
    finite = y[np.isfinite(y)]
    if not finite.size:
        return f"{x.size} points, all NaN"
    return (
        f"{x.size} pts, x {np.nanmin(x):.3g}..{np.nanmax(x):.3g}, "
        f"y {finite.min():.3g}..{finite.max():.3g}"
    )


def _solver_of(res: SimResult) -> str:
    return str(res.metadata.get("solver") or "")


def _devices_alive() -> "list[str]":
    import devsim

    return list(devsim.get_device_list())


def run_all(config: DeviceConfig, sections: "tuple[str, ...] | None" = None) -> Report:
    """Run the selected sections (default: all) for ``config``.

    Never raises for a section failure; see ``Report.rows()``.
    """
    table = _sections()
    names = tuple(sections) if sections is not None else tuple(table)
    unknown = [n for n in names if n not in table]
    if unknown:
        raise ValueError(f"unknown section(s) {unknown}; choose from {list(table)}")

    t_all = time.perf_counter()
    done: "list[SectionResult]" = []
    for name in names:
        title, runner, one_d_only = table[name]
        if one_d_only and config.half_width_um is not None:
            done.append(SectionResult(name, title, "skipped", reason=_ONE_D_ONLY))
            continue
        t0 = time.perf_counter()
        try:
            res = runner(config)
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
            sec = SectionResult(
                name,
                title,
                "error" if err else "ok",
                solver=_solver_of(res),
                reason=str(err) if err else "",
                wall_s=time.perf_counter() - t0,
                summary="" if err else _summarize(res),
                result=res,
            )
            done.append(sec)
        leaked = _devices_alive()
        if leaked:  # would couple this section to the next one
            last = done[-1]
            last.status = "error"
            last.reason = (last.reason + "; " if last.reason else "") + (
                f"leaked devsim device(s) {leaked}"
            )
            import devsim

            for d in leaked:
                devsim.delete_device(device=d)
    return Report(config, done, time.perf_counter() - t_all)
