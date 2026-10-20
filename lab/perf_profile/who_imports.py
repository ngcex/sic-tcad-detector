"""Which module first imports a heavy dependency during `import etna`?

usage: uv run python lab/perf_profile/who_imports.py matplotlib.pyplot scipy.stats pandas
"""

from __future__ import annotations

import builtins
import sys
import traceback

targets = set(sys.argv[1:] or ["matplotlib.pyplot", "scipy.stats", "pandas"])
_orig = builtins.__import__
seen = set()


def hook(name, globals=None, locals=None, fromlist=(), level=0):
    if name in targets and name not in sys.modules and name not in seen:
        seen.add(name)
        stack = [f for f in traceback.extract_stack()[:-1] if "/etna/" in f.filename]
        where = stack[-1] if stack else None
        print(
            f"{name} <- {where.filename.split('/etna/')[-1]}:{where.lineno}"
            if where
            else f"{name} <- (non-etna)"
        )
    return _orig(name, globals, locals, fromlist, level)


builtins.__import__ = hook
import etna  # noqa: E402,F401

builtins.__import__ = _orig
print("loaded:", {t: t in sys.modules for t in targets})
