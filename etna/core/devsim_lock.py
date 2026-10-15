"""Process-wide lock for devsim, whose state is global and not thread-safe.

Concurrent solves (e.g. two Streamlit sessions) delete each other's devices or
segfault the process, so every public entry point that touches devsim runs
under this re-entrant lock.
"""

import functools
import threading

DEVSIM_LOCK = threading.RLock()


def serialized(fn):
    """Run ``fn`` while holding the devsim lock."""

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        with DEVSIM_LOCK:
            return fn(*args, **kwargs)

    return wrapper
