"""Licence inventory of ETNA's RUNTIME dependency closure (issue #3 P3.2).

    uv run python scripts/license_inventory.py            # print the markdown table
    uv run python scripts/license_inventory.py --write    # refresh THIRD_PARTY_LICENSES.md
    uv run python scripts/license_inventory.py --check    # exit 1 if the file is stale or
                                                          # a package is copyleft/unknown
                                                          # without a reviewed note

The package list is the locked runtime closure, with no dev, viz or other extras:
`uv export --frozen --no-dev --no-emit-project`. Each licence is read from the INSTALLED
distribution's metadata (`importlib.metadata`), so run the script in the project
environment (`uv sync --extra dev`). The fields are tried in this order:

1. `License-Expression` (PEP 639);
2. a short `License` field;
3. the `License ::` trove classifiers.

The licence is reported as declared by each package and not legally reviewed. Packages
whose closure entry is not installed on this platform (environment markers) are listed
as such.

Flags:

- `copyleft`: a GPL-family licence (GPL, LGPL, AGPL, MPL, EPL);
- `unknown`: no licence metadata.

A flagged package needs a reviewed note in `REVIEWED`. Only the table between the
GENERATED markers of THIRD_PARTY_LICENSES.md is rewritten; the prose around it is kept.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "THIRD_PARTY_LICENSES.md"
BEGIN = "<!-- BEGIN GENERATED: scripts/license_inventory.py -->"
END = "<!-- END GENERATED -->"

COPYLEFT = re.compile(
    r"\b(A?GPL|LGPL|GNU (Lesser |Library |Affero )?General Public|MPL|Mozilla Public|EPL|Eclipse Public)",
    re.I,
)

# Reviewed notes for flagged packages (keep short; details in THIRD_PARTY_LICENSES.md).
REVIEWED: dict[str, str] = {
    "devsim": (
        "Apache-2.0 for devsim itself; its wheel bundles libumfpack_lgpl (UMFPACK 5.1, "
        "LGPL-2.1+ per the wheel's NOTICE). See the devsim section below."
    ),
    "certifi": "MPL-2.0 (CA bundle data, unmodified; weak copyleft, file-level).",
}

# Platform-specific closure entries (environment markers): their licence as declared on
# PyPI (checked 2026-10-21), used on EVERY platform so the table does not depend on the
# machine that generates it. A new marker-restricted package missing here is flagged.
PLATFORM_SPECIFIC: dict[str, str] = {
    "colorama": "BSD License",
    "tzdata": "Apache-2.0",
    "watchdog": "Apache-2.0",
}


def runtime_closure() -> list[tuple[str, str, str]]:
    """(name, locked version, environment marker or '') of the locked runtime closure."""
    out = subprocess.run(
        [
            "uv",
            "export",
            "--frozen",
            "--no-dev",
            "--no-emit-project",
            "--no-hashes",
            "--no-header",
            "--no-annotate",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    pkgs = []
    for line in out.splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "-")):
            continue
        req, _, marker = line.partition(";")
        name, _, ver = req.strip().partition("==")
        pkgs.append((name.strip(), ver.strip(), marker.strip()))
    return sorted(pkgs, key=lambda p: p[0].lower())


def declared_licence(dist: metadata.Distribution) -> str:
    md = dist.metadata
    expr = md.get("License-Expression")
    if expr:
        return expr.strip()
    lic = (md.get("License") or "").strip()
    if lic and len(lic) <= 60 and "\n" not in lic and lic.upper() != "UNKNOWN":
        return lic
    classifiers = [
        c.split("::")[-1].strip()
        for c in md.get_all("Classifier") or []
        if c.startswith("License ::")
    ]
    classifiers = [c for c in classifiers if c not in ("OSI Approved",)]
    if classifiers:
        return " / ".join(dict.fromkeys(classifiers))
    if lic and lic.upper() != "UNKNOWN":
        first = lic.splitlines()[0].strip()
        return f"{first[:60]} (full text in metadata)"
    return ""


def inventory() -> list[dict]:
    rows = []
    for name, locked, marker in runtime_closure():
        if marker and name.lower() in PLATFORM_SPECIFIC:
            lic = PLATFORM_SPECIFIC[name.lower()]
            rows.append(dict(name=name, version=locked, licence=lic, flag="", marker=marker))
            continue
        try:
            dist = metadata.distribution(name)
        except metadata.PackageNotFoundError:
            rows.append(
                dict(
                    name=name,
                    version=locked,
                    licence="",
                    flag="unknown (not installed here)",
                    marker=marker,
                )
            )
            continue
        lic = declared_licence(dist)
        flag = "unknown" if not lic else ("copyleft" if COPYLEFT.search(lic) else "")
        if dist.version != locked:
            flag = (flag + "; " if flag else "") + f"installed {dist.version} != locked"
        rows.append(
            dict(name=name, version=locked, licence=lic, flag=flag, marker=marker)
        )
    return rows


def table(rows: list[dict]) -> str:
    lines = [
        "| Package | Locked version | Declared licence (package metadata) | Flag | Note |",
        "| --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        note = REVIEWED.get(r["name"].lower(), "")
        if r["marker"]:
            note = (note + " " if note else "") + f"Only when `{r['marker']}`; licence as declared on PyPI."
        lic = r["licence"].replace("|", "/") or "-"
        lines.append(
            f"| {r['name']} | {r['version']} | {lic} | {r['flag'] or '-'} | {note or '-'} |"
        )
    return "\n".join(lines)


def problems(rows: list[dict]) -> list[str]:
    bad = []
    for r in rows:
        if ("copyleft" in r["flag"] or "unknown" in r["flag"]) and r[
            "name"
        ].lower() not in REVIEWED:
            bad.append(
                f"{r['name']}: {r['flag']} ({r['licence'] or 'no licence metadata'})"
            )
    return bad


def render_file(current: str, body: str) -> str:
    if BEGIN not in current or END not in current:
        raise SystemExit(f"{OUT.name}: GENERATED markers not found")
    head, _, rest = current.partition(BEGIN)
    _, _, tail = rest.partition(END)
    return f"{head}{BEGIN}\n{body}\n{END}{tail}"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    g = ap.add_mutually_exclusive_group()
    g.add_argument(
        "--write", action="store_true", help=f"rewrite the table in {OUT.name}"
    )
    g.add_argument(
        "--check",
        action="store_true",
        help=f"fail if {OUT.name} is stale or a flag is unreviewed",
    )
    args = ap.parse_args(argv)

    rows = inventory()
    body = table(rows)
    bad = problems(rows)
    if args.write:
        OUT.write_text(
            render_file(OUT.read_text(encoding="utf-8"), body), encoding="utf-8"
        )
        print(f"wrote {OUT.relative_to(ROOT)} ({len(rows)} packages)")
    elif args.check:
        current = OUT.read_text(encoding="utf-8")
        if render_file(current, body) != current:
            print(
                f"{OUT.name} is stale: run scripts/license_inventory.py --write",
                file=sys.stderr,
            )
            return 1
    else:
        print(body)
    for b in bad:
        print(f"UNREVIEWED: {b}", file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
