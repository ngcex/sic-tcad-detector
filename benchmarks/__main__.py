"""``uv run python -m benchmarks [--slow] [--only ID ...] [--out DIR] [--list]``

Evaluates the fast subset of ``benchmarks/registry.toml`` (7-15 s on an M4 with a warm
numba cache; a cold cache adds the JIT compilation of every engine: 162 s measured with
an empty NUMBA_CACHE_DIR at load average ~50) and writes ``<out>/fast.json`` and
``<out>/fast.md`` (``all.*`` with ``--slow``). Exit status 1 if any gate fails.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from benchmarks.runner import (
    RESULTS,
    load_registry,
    run,
    select,
    update_report,
    write,
)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m benchmarks", description=__doc__)
    ap.add_argument("--slow", action="store_true", help="also run the slow entries")
    ap.add_argument("--only", nargs="+", metavar="ID", help="run only these ids")
    ap.add_argument("--out", type=Path, default=RESULTS, help="output directory")
    ap.add_argument("--list", action="store_true", help="list the registry and exit")
    ap.add_argument("--update-report", action="store_true",
                    help="regenerate the registry table in docs/VV_REPORT.md and exit")  # fmt: skip
    a = ap.parse_args(argv)
    entries = load_registry()
    if a.update_report:
        changed = update_report(entries)
        print("docs/VV_REPORT.md", "updated" if changed else "already in sync")
        return 0
    if a.list:
        for e in entries:
            print(f"{e.id:34s} {e.status:14s} {e.speed:7s} {e.quantity}")
        return 0
    chosen = select(entries, slow=a.slow, ids=set(a.only) if a.only else None)
    report = run(chosen)
    stem = "selection" if a.only else ("all" if a.slow else "fast")
    j, m = write(report, a.out, stem)
    for r in report["rows"]:
        e = r.get("E")
        e_txt = "-" if e is None else f"{e:+.3e}"
        print(f"{r['id']:34s} {r['status']:14s} E={e_txt:>11s}  {r['verdict']}  "
              f"({r['runtime_s']} s)")  # fmt: skip
    print(f"\n{len(report['rows'])} entries, {report['n_failed']} failed gates, "
          f"{report['wall_s']} s -> {j}, {m}")  # fmt: skip
    return 1 if report["n_failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
