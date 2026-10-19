"""Check that every scripts/create_notebook_*.py reproduces its committed notebook.

Each generator runs with ``nbformat.write`` intercepted: nothing is written to the
repository and nothing is executed (generators 17/18 execute the notebook only
after writing it, so they stop at the intercepted write). The cell types and
sources of the generated notebook are compared with notebooks/<name>.ipynb.

Run from anywhere:
    uv run python scripts/check_notebook_generators.py          # all generators
    uv run python scripts/check_notebook_generators.py 16 19    # a subset
Exit code 1 when any generator drifts from its notebook (diff printed).
"""

import difflib
import os
import pathlib
import runpy
import sys

import nbformat

REPO = pathlib.Path(__file__).resolve().parents[1]


class _Captured(Exception):
    pass


def _cells(nb):
    lines = []
    for i, c in enumerate(nb.cells):
        lines.append(f"### cell {i} [{c.cell_type}]")
        lines.extend(c.source.splitlines())
    return lines


def main(only):
    captured = {}
    real_write = nbformat.write

    def fake_write(nb, path, *args, **kwargs):
        captured["nb"] = nb
        captured["path"] = path if isinstance(path, (str, os.PathLike)) else path.name
        raise _Captured

    nbformat.write = fake_write
    drift = 0
    try:
        for gen in sorted((REPO / "scripts").glob("create_notebook_*.py")):
            if only and not any(o in gen.name for o in only):
                continue
            captured.clear()
            cwd = os.getcwd()
            # Generator 03 writes ../notebooks (run from scripts/); 15-19 write
            # notebooks/ (repo root); 08 uses an absolute path. Only the basename
            # of the target matters here.
            os.chdir(REPO / "scripts")
            try:
                runpy.run_path(str(gen), run_name="__main__")
            except _Captured:
                pass
            finally:
                os.chdir(cwd)
            if "nb" not in captured:
                print(f"{gen.name}: no notebook written")
                drift += 1
                continue
            target = REPO / "notebooks" / os.path.basename(str(captured["path"]))
            if not target.exists():
                print(f"{gen.name}: target {target.name} is not committed")
                drift += 1
                continue
            committed = nbformat.read(target, as_version=4)
            diff = list(
                difflib.unified_diff(
                    _cells(committed),
                    _cells(captured["nb"]),
                    f"notebooks/{target.name}",
                    gen.name,
                    lineterm="",
                    n=1,
                )
            )
            if diff:
                drift += 1
                print(f"{gen.name}: DRIFT vs {target.name}")
                print("\n".join(diff))
            else:
                print(f"{gen.name}: in sync with {target.name}")
    finally:
        nbformat.write = real_write
    return 1 if drift else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
