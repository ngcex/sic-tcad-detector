"""Shared reader for the measurement CSV files (C-V, SIMS/SRP, I(T)).

File layout (all three kinds)::

    # free comment lines start with '#'
    # key: value        <- metadata line (known keys are parsed, unknown keys kept as text)
    V_V,C_pF,sigma_C_pF  <- header: one column per quantity, the UNIT IS IN THE NAME
    ...data rows...

Units live in the column names (``C_pF``, ``depth_nm``, ``I_pA``, ``T_K``): each loader
has a fixed table of accepted names with their conversion factor to the internal unit, and
an unknown name is an error that lists the accepted ones. Every problem found is
collected and raised together as one ``MeasurementFileError`` whose lines name the file
line (1-based, counting comment and header lines) and the column.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field
from pathlib import Path

MAX_REPORTED_ISSUES = 25


class MeasurementFileError(ValueError):
    """A measurement file does not follow its documented schema.

    ``issues`` holds one message per problem, each naming the file line and/or column.
    """

    def __init__(self, source: str, issues):
        self.source = source
        self.issues = list(issues)
        shown = self.issues[:MAX_REPORTED_ISSUES]
        more = len(self.issues) - len(shown)
        msg = f"{source}: {len(self.issues)} problem(s)\n  - " + "\n  - ".join(shown)
        if more > 0:
            msg += f"\n  - ... and {more} more"
        super().__init__(msg)


@dataclass(frozen=True)
class Column:
    """Accepted column: header name -> (quantity, factor to the internal unit)."""

    name: str
    quantity: str
    factor: float
    unit: str  # human-readable unit of the header name
    internal_unit: str


@dataclass
class RawTable:
    """Parsed file: numeric columns by quantity (internal units), metadata, provenance."""

    source: str
    sha256: str
    columns: dict  # quantity -> list[float] (internal unit; NaN = empty optional cell)
    column_names: dict  # quantity -> header name used in the file
    lines: list  # file line number (1-based) of each data row
    meta: dict = field(default_factory=dict)  # parsed metadata (known keys typed)
    comments: list = field(default_factory=list)


def _decode(data: bytes, source: str) -> str:
    """UTF-8 (BOM allowed); anything else is a ``MeasurementFileError`` naming the line
    and byte, not a bare ``UnicodeDecodeError`` (e.g. a Latin-1 'µm' / '°C' comment)."""
    try:
        return data.decode("utf-8-sig")
    except UnicodeDecodeError as e:
        line = data[: e.start].count(b"\n") + 1
        raise MeasurementFileError(
            source,
            [
                f"line {line}: byte 0x{data[e.start]:02x} is not UTF-8 (the file is "
                "probably Latin-1 / Windows-1252, e.g. a 'µ' or '°' in a comment); "
                "save it as UTF-8"
            ],
        ) from None


def read_text(src) -> "tuple[str, str]":
    """``src`` (path, bytes, str content or file-like) -> (text, source name)."""
    if isinstance(src, (bytes, bytearray)):
        return _decode(bytes(src), "<uploaded bytes>"), "<uploaded bytes>"
    if hasattr(src, "read"):
        data = src.read()
        name = str(getattr(src, "name", "<stream>"))
        if isinstance(data, bytes):
            data = _decode(data, name)
        return data, name
    if isinstance(src, Path) or (
        isinstance(src, str) and "\n" not in src and len(src) < 4096
    ):
        p = Path(src)
        if not p.exists():
            raise FileNotFoundError(f"measurement file not found: {p}")
        text = _decode(p.read_bytes(), p.name)
        # universal newlines, as Path.read_text gave (keeps the recorded sha256)
        return text.replace("\r\n", "\n").replace("\r", "\n"), p.name
    return str(src), "<text>"


def _number(cell: str):
    """Float of ``cell`` or None when it is not a finite number."""
    try:
        v = float(cell)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def parse_table(
    text: str,
    source: str,
    accepted: "dict[str, Column]",
    required: "list[tuple[str, ...]]",
    exclusive: "list[tuple[str, ...]]" = (),
    meta_keys: "dict[str, type]" = None,
    kind: str = "measurement",
    complete_if_present: "tuple[str, ...]" = (),
) -> RawTable:
    """Parse and validate one CSV against an accepted-column table.

    ``required``: groups of quantities, at least one of each group must be present.
    ``exclusive``: groups of quantities of which at most one may be present.
    ``meta_keys``: metadata key -> type (float or str); an empty value means "not given".
    Required-quantity cells must be finite numbers on every row; optional cells may be
    empty (NaN) but, when filled, must be finite numbers. ``complete_if_present``:
    optional quantities (uncertainties) that must be filled on every row or on none.
    """
    meta_keys = meta_keys or {}
    issues: list[str] = []
    meta: dict = {}
    comments: list[str] = []
    header = None
    header_line = None
    rows: list[tuple[int, list[str]]] = []
    for ln, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            body = line.lstrip("#").strip()
            key, sep, value = body.partition(":")
            key_n = key.strip().lower().replace(" ", "_")
            if sep and key_n in meta_keys:
                value = value.strip()
                if value == "":
                    continue  # template placeholder: not given
                typ = meta_keys[key_n]
                if typ is float:
                    v = _number(value)
                    if v is None:
                        issues.append(
                            f"line {ln}: metadata '{key_n}' must be a finite number "
                            f"(got {value!r})"
                        )
                        continue
                    meta[key_n] = v
                else:
                    meta[key_n] = value
            elif body:
                comments.append(body)
            continue
        cells = [c.strip() for c in line.split(",")]
        if header is None:
            header, header_line = cells, ln
        else:
            rows.append((ln, cells))
    if header is None:
        raise MeasurementFileError(
            source, [f"no header row found (expected a {kind} CSV header line)"]
        )

    # ---- header
    names_ok = ", ".join(sorted(accepted))
    by_quantity: dict[str, tuple[int, Column]] = {}
    seen = set()
    for j, name in enumerate(header):
        if name == "":
            issues.append(f"line {header_line}: column {j + 1} has an empty name")
            continue
        if name in seen:
            issues.append(f"line {header_line}: column '{name}' appears twice")
            continue
        seen.add(name)
        col = accepted.get(name)
        if col is None:
            issues.append(
                f"line {header_line}: unknown column '{name}' (accepted: {names_ok})"
            )
            continue
        if col.quantity in by_quantity:
            other = by_quantity[col.quantity][1].name
            issues.append(
                f"line {header_line}: columns '{other}' and '{name}' give the same "
                f"quantity ({col.quantity}); keep one"
            )
            continue
        by_quantity[col.quantity] = (j, col)
    for group in required:
        if not any(q in by_quantity for q in group):
            opts = sorted(n for n, c in accepted.items() if c.quantity in group)
            issues.append(
                f"line {header_line}: missing required column (one of: {', '.join(opts)})"
            )
    for group in exclusive:
        present = [by_quantity[q][1].name for q in group if q in by_quantity]
        if len(present) > 1:
            issues.append(
                f"line {header_line}: columns {', '.join(present)} are exclusive; keep one"
            )
    if issues:
        raise MeasurementFileError(source, issues)
    if not rows:
        raise MeasurementFileError(
            source,
            [f"no data rows after the header (line {header_line}): the file is empty"],
        )

    required_q = {q for g in required for q in g if q in by_quantity}
    cols: dict[str, list[float]] = {q: [] for q in by_quantity}
    lines = []
    for ln, cells in rows:
        if len(cells) != len(header):
            issues.append(
                f"line {ln}: {len(cells)} cells, the header has {len(header)} columns"
            )
            continue
        lines.append(ln)
        for q, (j, col) in by_quantity.items():
            cell = cells[j]
            if cell == "":
                if q in required_q:
                    issues.append(f"line {ln}, column '{col.name}': empty cell")
                cols[q].append(float("nan"))
                continue
            v = _number(cell)
            if v is None:
                issues.append(
                    f"line {ln}, column '{col.name}': {cell!r} is not a finite number"
                )
                cols[q].append(float("nan"))
                continue
            cols[q].append(v * col.factor)
    for q in complete_if_present:
        if q not in by_quantity or q in required_q:
            continue
        empty = [ln for ln, v in zip(lines, cols[q]) if math.isnan(v)]
        if empty and len(empty) < len(lines):
            shown = ", ".join(str(x) for x in empty[:10])
            more = f" and {len(empty) - 10} more" if len(empty) > 10 else ""
            issues.append(
                f"column '{by_quantity[q][1].name}': empty cell on line(s) {shown}{more} "
                "while other rows give one; fill every row or leave the column empty "
                "(a partial uncertainty column would silently drop all weights)"
            )
    if issues:
        raise MeasurementFileError(source, issues)
    return RawTable(
        source=source,
        sha256=hashlib.sha256(text.encode("utf-8")).hexdigest(),
        columns=cols,
        column_names={q: c.name for q, (_, c) in by_quantity.items()},
        lines=lines,
        meta=meta,
        comments=comments,
    )


def constant_column(raw: RawTable, quantity: str, issues: list, rel_tol=1e-6):
    """The single value of an optional per-row column (None if absent / all empty).

    Every filled cell must agree (relative ``rel_tol``); otherwise an issue naming the
    first differing line is appended.
    """
    vals = raw.columns.get(quantity)
    if vals is None:
        return None
    filled = [(ln, v) for ln, v in zip(raw.lines, vals) if not math.isnan(v)]
    if not filled:
        return None
    v0 = filled[0][1]
    for ln, v in filled[1:]:
        if abs(v - v0) > rel_tol * max(abs(v0), 1e-300):
            issues.append(
                f"line {ln}, column '{raw.column_names[quantity]}': {v:g} differs from "
                f"{v0:g} on line {filled[0][0]} (one file = one value; split the file)"
            )
            return v0
    return v0
