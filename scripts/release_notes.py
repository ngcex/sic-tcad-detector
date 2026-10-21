"""Print the CHANGELOG.md section of one release (used by .github/workflows/release.yml).

    uv run python scripts/release_notes.py 5.0.0            # section body to stdout
    uv run python scripts/release_notes.py --check-version  # CHANGELOG has a section for
                                                            # etna/_version.py; exit 0/1

Exits 1 when the section is missing or empty, so a release cannot go out without notes.
No third-party imports: runs on a bare Python.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def package_version(root: Path = ROOT) -> str:
    text = (root / "etna" / "_version.py").read_text(encoding="utf-8")
    m = re.search(r'__version__\s*=\s*"([^"]+)"', text)
    if not m:
        raise SystemExit("etna/_version.py: no __version__ found")
    return m.group(1)


def changelog_section(version: str, changelog: Path = ROOT / "CHANGELOG.md") -> str:
    """Body of the ``## [version]`` section (heading excluded), stripped; '' if absent.

    The section ends at the next level-2 heading or at the link-reference block.
    """
    lines = changelog.read_text(encoding="utf-8").splitlines()
    head = re.compile(rf"^## \[{re.escape(version)}\](\s|$)")
    out: list[str] = []
    inside = False
    for line in lines:
        if inside:
            if line.startswith("## ") or re.match(r"^\[[^\]]+\]:\s", line):
                break
            out.append(line)
        elif head.match(line):
            inside = True
    return "\n".join(out).strip()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("version", nargs="?", help="version without the leading v")
    ap.add_argument(
        "--check-version",
        action="store_true",
        help="check that CHANGELOG.md has a non-empty section for etna/_version.py",
    )
    args = ap.parse_args(argv)
    version = args.version or package_version()
    if args.version and args.version.startswith("v"):
        version = args.version[1:]
    body = changelog_section(version)
    if not body:
        print(
            f"CHANGELOG.md has no non-empty section '## [{version}]'", file=sys.stderr
        )
        return 1
    if not args.check_version:
        print(body)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
