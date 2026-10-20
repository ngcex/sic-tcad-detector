"""Regenerate data/nist_star/*.csv from the raw NIST responses in data/nist_star/raw/.

The raw files are the unmodified HTTP response bodies of the NIST STAR text interface
(POST https://physics.nist.gov/cgi-bin/Star/ap_table-t.pl, prog=PSTAR|ASTAR,
matno=<NIST material number>, ShowDefault=on), downloaded 2026-10-20. See
data/nist_star/README.md. This script only re-formats them; it adds no numbers.

    uv run python scripts/build_nist_star_tables.py
"""

import hashlib
from pathlib import Path

from etna.core.stopping_power import (
    COLUMNS,
    DATA_DIR,
    NIST_MATERIALS,
    PROGRAMS,
    parse_nist_star_html,
)

URL = "https://physics.nist.gov/cgi-bin/Star/ap_table-t.pl"
RETRIEVED = "2026-10-20"


def main():
    raw = DATA_DIR / "raw"
    for particle, prog in PROGRAMS.items():
        for material, matno in NIST_MATERIALS.items():
            src = raw / f"{prog}_{matno}.html"
            text = src.read_text(encoding="latin-1")
            name, table = parse_nist_star_html(text)
            sha = hashlib.sha256(src.read_bytes()).hexdigest()
            out = DATA_DIR / f"{prog}_{material}.csv"
            lines = [
                f"# NIST {prog.upper()} ({particle}) -- {name}",
                f"# source: POST {URL} prog={prog.upper()} matno={matno} ShowDefault=on",
                f"# retrieved: {RETRIEVED}; raw response: raw/{src.name} sha256={sha}",
                "# units: energy MeV; stopping powers MeV cm2/g; ranges g/cm2",
                ",".join(COLUMNS),
            ]
            for row in table:
                lines.append(",".join(f"{v:.6g}" for v in row))
            out.write_text("\n".join(lines) + "\n", encoding="utf-8")
            print(
                f"wrote {out.relative_to(Path.cwd()) if out.is_relative_to(Path.cwd()) else out} ({len(table)} rows, {name})"
            )


if __name__ == "__main__":
    main()
