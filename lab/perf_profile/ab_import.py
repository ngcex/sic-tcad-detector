"""A/B `import etna` time: OLD tree (PYTHONPATH=<old>) vs the worktree, interleaved.

usage: uv run python lab/perf_profile/ab_import.py <old_tree_dir> [repeats]
  old_tree_dir: a directory holding an `etna/` package (e.g. `git archive HEAD etna | tar -x -C DIR`)
Prints min / median of wall and cpu (user+sys of the child) over the repeats.
"""

from __future__ import annotations

import os
import statistics
import subprocess
import sys

CODE = (
    "import time,resource,sys;t=time.perf_counter();import etna;w=time.perf_counter()-t;"
    "r=resource.getrusage(resource.RUSAGE_SELF);"
    "assert etna.__file__.startswith(sys.argv[1]), etna.__file__;"
    "print(w, r.ru_utime+r.ru_stime, int('matplotlib.pyplot' in sys.modules), int('scipy.stats' in sys.modules))"
)


def once(root: str) -> tuple[float, float, int, int]:
    env = {**os.environ, "PYTHONPATH": root}
    out = subprocess.run(
        [sys.executable, "-c", CODE, root],
        capture_output=True,
        text=True,
        env=env,
        cwd=root,
    )
    line = [ln for ln in out.stdout.splitlines() if ln and ln[0].isdigit()][-1]
    w, c, plt, st = line.split()
    return float(w), float(c), int(plt), int(st)


if __name__ == "__main__":
    old = os.path.abspath(sys.argv[1])
    new = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    res = {"old": [], "new": []}
    for _ in range(reps):
        for k, root in (("old", old), ("new", new)):
            res[k].append(once(root))
    for k, rows in res.items():
        w = [r[0] for r in rows]
        c = [r[1] for r in rows]
        print(
            f"{k}: wall min {min(w):.3f} med {statistics.median(w):.3f} s | cpu min {min(c):.3f} "
            f"med {statistics.median(c):.3f} s | pyplot loaded {rows[0][2]} scipy.stats loaded {rows[0][3]}"
        )
