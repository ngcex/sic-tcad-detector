"""Run a list of run_one.py configs, one process each, N in parallel.

usage: uv run python lab/finite_implant/sweep.py SWEEPFILE [N]
  SWEEPFILE lines: "<tag> key=value key=value ..." -> notes/<tag>.json, notes/<tag>.log
"""

import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent


def run(line):
    tag, *args = line.split()
    with open(HERE / "notes" / f"{tag}.log", "w") as log:
        rc = subprocess.call(
            [
                sys.executable,
                str(HERE / "run_one.py"),
                str(HERE / "notes" / f"{tag}.json"),
                *args,
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    print(f"{tag} rc={rc}", flush=True)


if __name__ == "__main__":
    lines = [
        l
        for l in Path(sys.argv[1]).read_text().splitlines()
        if l.strip() and not l.startswith("#")
    ]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    with ThreadPoolExecutor(n) as ex:
        list(ex.map(run, lines))
