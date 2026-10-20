"""User-editable layout / mask dimensions (``data/foundry/layout_params.toml``).

The GDSII / SVG layout (``etna.export.gds``), the spec sheet and the cross-section all
read their lateral and passivation dimensions from one TOML file, so mask data can be
dropped in without touching code. Every dimension is a table::

    [dosimeter.guard_ring_gap]
    value = 150.0          # > 0 (counts: integer >= 0)
    unit = "um"            # "mm" | "um" | "nm"; "count" for counts
    class = "PLACEHOLDER"  # DOCUMENTED | MEASURED_FROM_FIGURE | PLACEHOLDER
    source = "..."         # non-empty
    uncertainty = "unknown"  # >= 0 in `unit`, or "unknown"
    doc_item = "..."       # optional: row of deliverables/DIMENSIONS-provenance.md
    note = "..."           # optional

File choice: ``load_layout_params(path)`` > env var ``ETNA_LAYOUT_PARAMS`` > the
repository default. ``parse_layout_params(mapping)`` validates an in-memory mapping
(e.g. an uploaded file). Every problem raises ``LayoutParamsError`` naming the field
(``dosimeter.guard_ring_width.value: ...``). Values still classed PLACEHOLDER are listed
by ``LayoutParams.placeholders(device)`` and announced by ``warn_placeholders`` with a
``PlaceholderDimensionWarning``.
"""

from __future__ import annotations

import math
import os
import tomllib
import warnings
from dataclasses import dataclass, field
from pathlib import Path

ENV_VAR = "ETNA_LAYOUT_PARAMS"
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PARAMS_PATH = REPO_ROOT / "data" / "foundry" / "layout_params.toml"
SCHEMA_VERSION = 1
CLASSES = ("DOCUMENTED", "MEASURED_FROM_FIGURE", "PLACEHOLDER")
DEVICES = ("dosimeter", "microdosimeter")
LENGTH_UNITS_UM = {"mm": 1e3, "um": 1.0, "nm": 1e-3}

# name -> kind ("length" | "count"); all required except the OPTIONAL ones
DIMENSION_KINDS = {
    "active_side": "length",
    "die_side": "length",
    "guard_ring_count": "count",
    "guard_ring_gap": "length",
    "guard_ring_width": "length",
    "guard_ring_spacing": "length",
    "pad_size": "length",
    "dicing_street": "length",
    "oxide_thickness": "length",
    "passivation_thickness": "length",
}
OPTIONAL = frozenset({"passivation_thickness"})
_KEYS_REQUIRED = ("value", "unit", "class", "source", "uncertainty")
_KEYS_OPTIONAL = ("doc_item", "note")


class LayoutParamsError(ValueError):
    """Invalid layout parameter file; the message names the offending field."""


class PlaceholderDimensionWarning(UserWarning):
    """An output was produced from dimensions still classed PLACEHOLDER."""


@dataclass(frozen=True)
class Dimension:
    device: str
    name: str
    value: float
    unit: str
    cls: str
    source: str
    uncertainty: "float | None"  # None = "unknown"
    doc_item: str = ""
    note: str = ""

    @property
    def is_placeholder(self) -> bool:
        return self.cls == "PLACEHOLDER"

    @property
    def um(self) -> float:
        """Value in um (lengths only)."""
        if self.unit == "count":
            raise TypeError(f"{self.device}.{self.name} is a count, not a length")
        return float(self.value) * LENGTH_UNITS_UM[self.unit]

    def text(self) -> str:
        """'2.7 +- 0.15 mm' / '1 count'."""
        v = f"{self.value:g}" if self.unit != "count" else f"{int(self.value)}"
        unc = "" if self.uncertainty is None else f" +- {self.uncertainty:g}"
        return f"{v}{unc} {self.unit}"

    def label(self, what: str = "") -> str:
        """'CLASS | what: value | src: source' (the layout label convention)."""
        head = f"{what or self.name}: {self.text()}"
        return f"{self.cls} | {head} | src: {self.source}"

    def as_dict(self) -> dict:
        return dict(
            name=self.name,
            value=self.value,
            unit=self.unit,
            value_um=None if self.unit == "count" else self.um,
            cls=self.cls,
            source=self.source,
            uncertainty=self.uncertainty,
            doc_item=self.doc_item,
            note=self.note,
        )


@dataclass(frozen=True)
class LayoutParams:
    path: str  # file it was read from, or "<mapping>"
    devices: dict = field(default_factory=dict)  # device -> {name: Dimension}

    @property
    def display_path(self) -> str:
        """``path`` relative to the repository when inside it (no local paths in exports)."""
        try:
            return str(Path(self.path).resolve().relative_to(REPO_ROOT))
        except (ValueError, OSError):
            return self.path

    def dim(self, device: str, name: str) -> Dimension:
        return self.devices[device][name]

    def get(self, device: str, name: str) -> "Dimension | None":
        return self.devices.get(device, {}).get(name)

    def used(self, device: str) -> "list[str]":
        """Dimension names that enter the outputs: the ring gap / width are unused
        with no ring, the ring spacing with fewer than two rings."""
        d = self.devices[device]
        n = int(d["guard_ring_count"].value)
        skip = set()
        if n < 2:
            skip.add("guard_ring_spacing")
        if n == 0:
            skip |= {"guard_ring_gap", "guard_ring_width"}
        return [k for k in d if k not in skip]

    def placeholders(self, device: "str | None" = None) -> "list[str]":
        """'device.name' of every PLACEHOLDER dimension that is used (see ``used``),
        of one device or of all."""
        devs = [device] if device else list(self.devices)
        return [
            f"{d}.{n}"
            for d in devs
            for n in self.used(d)
            if self.devices[d][n].is_placeholder
        ]

    def placeholder_note(self, device: str) -> str:
        ph = self.placeholders(device)
        if not ph:
            return f"no PLACEHOLDER dimension left for {device} ({self.display_path})"
        names = ", ".join(p.split(".", 1)[1] for p in ph)
        return (
            f"PLACEHOLDER dimensions ({device}), drawing values NOT mask data: "
            f"{names}; overwrite them in {self.display_path}"
        )


# --------------------------------------------------------------------------- loading


def resolve_params_path(path=None) -> Path:
    if path is not None:
        return Path(path)
    env = os.environ.get(ENV_VAR)
    return Path(env) if env else DEFAULT_PARAMS_PATH


def load_layout_params(path=None) -> LayoutParams:
    """Read and validate the layout parameter file (argument > $ETNA_LAYOUT_PARAMS >
    the repository default ``data/foundry/layout_params.toml``)."""
    p = resolve_params_path(path)
    try:
        raw = tomllib.loads(p.read_text(encoding="utf-8"))
    except FileNotFoundError as e:
        raise LayoutParamsError(f"layout parameter file not found: {p}") from e
    except tomllib.TOMLDecodeError as e:
        raise LayoutParamsError(f"{p}: not valid TOML ({e})") from e
    return parse_layout_params(raw, source=str(p))


def as_layout_params(params=None) -> LayoutParams:
    """LayoutParams from None (default file), a path, a mapping or a LayoutParams."""
    if isinstance(params, LayoutParams):
        return params
    if isinstance(params, dict):
        return parse_layout_params(params)
    return load_layout_params(params)


def _err(where: str, msg: str) -> LayoutParamsError:
    return LayoutParamsError(f"{where}: {msg}")


def _number(where, v):
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise _err(where, f"must be a number, got {v!r}")
    if not math.isfinite(v):
        raise _err(where, f"must be finite, got {v!r}")
    return float(v)


def _parse_dim(device: str, name: str, t) -> Dimension:
    where = f"{device}.{name}"
    if not isinstance(t, dict):
        raise _err(where, "must be a table with value/unit/class/source/uncertainty")
    for k in _KEYS_REQUIRED:
        if k not in t:
            raise _err(where, f"missing key '{k}'")
    extra = set(t) - set(_KEYS_REQUIRED) - set(_KEYS_OPTIONAL)
    if extra:
        raise _err(where, f"unknown key(s) {sorted(extra)}")
    kind = DIMENSION_KINDS[name]
    unit = t["unit"]
    if kind == "count":
        if unit != "count":
            raise _err(f"{where}.unit", f"must be 'count', got {unit!r}")
        v = t["value"]
        if isinstance(v, bool) or not isinstance(v, int) or v < 0:
            raise _err(f"{where}.value", f"must be an integer >= 0, got {v!r}")
        value = float(v)
    else:
        if unit not in LENGTH_UNITS_UM:
            raise _err(
                f"{where}.unit",
                f"must be one of {sorted(LENGTH_UNITS_UM)}, got {unit!r}",
            )
        value = _number(f"{where}.value", t["value"])
        if value <= 0:
            raise _err(f"{where}.value", f"must be > 0, got {t['value']!r}")
    cls = t["class"]
    if cls not in CLASSES:
        raise _err(f"{where}.class", f"must be one of {CLASSES}, got {cls!r}")
    src = t["source"]
    if not isinstance(src, str) or not src.strip():
        raise _err(f"{where}.source", "must be a non-empty string")
    unc = t["uncertainty"]
    if unc == "unknown":
        unc_v = None
    else:
        unc_v = _number(f"{where}.uncertainty", unc)
        if unc_v < 0:
            raise _err(
                f"{where}.uncertainty", f"must be >= 0 or 'unknown', got {unc!r}"
            )
    for k in _KEYS_OPTIONAL:
        if k in t and not isinstance(t[k], str):
            raise _err(f"{where}.{k}", "must be a string")
    return Dimension(
        device,
        name,
        value,
        unit,
        cls,
        src.strip(),
        unc_v,
        t.get("doc_item", ""),
        t.get("note", ""),
    )


def _check_geometry(device: str, d: dict) -> None:
    a = d["active_side"].um / 2.0
    b = d["die_side"].um / 2.0
    if b <= a:
        raise _err(
            f"{device}.die_side",
            f"die ({2 * b:g} um) must be larger than the active area ({2 * a:g} um)",
        )
    n = int(d["guard_ring_count"].value)
    r_out = a
    if n:
        r_out = (
            a
            + d["guard_ring_gap"].um
            + n * d["guard_ring_width"].um
            + (n - 1) * d["guard_ring_spacing"].um
        )
        if r_out >= b:
            raise _err(
                f"{device}.guard_ring_*",
                f"guard ring(s) do not fit inside the die: outer ring edge at "
                f"{r_out - a:g} um from the p+ edge, die edge at {b - a:g} um",
            )
    pad = d["pad_size"].um
    if pad > b - r_out:
        raise _err(
            f"{device}.pad_size",
            f"pad ({pad:g} um) does not fit between the outermost ring and the die edge "
            f"({b - r_out:g} um available)",
        )


def parse_layout_params(raw: dict, source: str = "<mapping>") -> LayoutParams:
    """Validate a parsed TOML mapping (see the module docstring for the schema)."""
    if not isinstance(raw, dict):
        raise LayoutParamsError(f"{source}: top level must be a table")
    ver = raw.get("schema_version")
    if ver != SCHEMA_VERSION:
        raise _err(
            f"{source}: schema_version", f"must be {SCHEMA_VERSION}, got {ver!r}"
        )
    unknown = set(raw) - {"schema_version", *DEVICES}
    if unknown:
        raise LayoutParamsError(
            f"{source}: unknown top-level key(s) {sorted(unknown)}; devices are {DEVICES}"
        )
    devices = {}
    for device in DEVICES:
        if device not in raw:
            raise LayoutParamsError(f"{source}: missing device table [{device}.*]")
        tbl = raw[device]
        if not isinstance(tbl, dict):
            raise LayoutParamsError(f"{source}: {device} must be a table")
        bad = set(tbl) - set(DIMENSION_KINDS)
        if bad:
            raise _err(
                f"{device}.{sorted(bad)[0]}",
                f"unknown dimension; known: {sorted(DIMENSION_KINDS)}",
            )
        for name in DIMENSION_KINDS:
            if name not in tbl and name not in OPTIONAL:
                raise _err(f"{device}.{name}", "missing dimension")
        dims = {n: _parse_dim(device, n, t) for n, t in tbl.items()}
        _check_geometry(device, dims)
        devices[device] = dims
    return LayoutParams(source, devices)


# --------------------------------------------------------------------------- helpers


def warn_placeholders(params: LayoutParams, device: str, where: str) -> "list[str]":
    """Warn (PlaceholderDimensionWarning) if ``device`` still has PLACEHOLDER values."""
    ph = params.placeholders(device)
    if ph:
        warnings.warn(
            f"{where}: {params.placeholder_note(device)}",
            PlaceholderDimensionWarning,
            stacklevel=3,
        )
    return ph


def device_for_config(config, params: LayoutParams) -> "str | None":
    """Layout device whose active square matches ``config.area_cm2`` (rel 1e-6)."""
    for device in DEVICES:
        side_cm = params.dim(device, "active_side").um * 1e-4
        if math.isclose(side_cm * side_cm, config.area_cm2, rel_tol=1e-6):
            return device
    return None
