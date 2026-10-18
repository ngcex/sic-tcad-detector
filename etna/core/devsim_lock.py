"""Process-wide lock for devsim, whose state is global and not thread-safe.

Concurrent solves (e.g. two Streamlit sessions) delete each other's devices or
segfault the process, so every public entry point that touches devsim runs
under this re-entrant lock.

devsim also prints ~65 lines per Newton solve through ``sys.stdout`` and has no
quiet ``debug_level`` ("info" still prints, any other value is verbose), so a
facade call used to flood the console (``run_all``: ~45k lines). ``serialized``
therefore also runs ``fn`` inside :func:`quiet_devsim`, which drops what the
calling thread writes to ``sys.stdout`` (etna itself reports through
``logging``/``warnings``, never ``print``). Set ``ETNA_DEVSIM_VERBOSE=1`` to see
devsim's output again. Use ``with quiet_devsim():`` around direct core calls in
notebooks or scripts.

Only Python-level ``sys.stdout`` is redirected: devsim writes its solver log
through the Python stream (the 5 math-library lines printed once at import are
the only C-level output). Redirecting file descriptor 1 would also silence the
other threads of the process (Streamlit server, pytest capture), so it is not done.
"""

from __future__ import annotations

import contextlib
import functools
import os
import sys
import threading

DEVSIM_LOCK = threading.RLock()

VERBOSE_ENV = "ETNA_DEVSIM_VERBOSE"


class _ThreadQuietStdout:
    """``sys.stdout`` proxy that drops writes from threads inside ``quiet_devsim``.

    Writes from every other thread go to the wrapped stream unchanged, so a
    devsim solve in one thread does not silence the rest of the process.
    """

    def __init__(self, wrapped):
        self._wrapped = wrapped
        self._quiet = {}  # thread id -> nesting depth

    def _is_quiet(self) -> bool:
        return threading.get_ident() in self._quiet

    def write(self, text):
        if self._is_quiet():
            return len(text)
        return self._wrapped.write(text)

    def writelines(self, lines):
        if not self._is_quiet():
            self._wrapped.writelines(lines)

    def flush(self):
        if not self._is_quiet():
            self._wrapped.flush()

    def __getattr__(self, name):
        return getattr(self._wrapped, name)


_PROXY_LOCK = threading.Lock()


def _verbose() -> bool:
    return os.environ.get(VERBOSE_ENV, "").strip() not in ("", "0")


@contextlib.contextmanager
def quiet_devsim():
    """Silence devsim's solver log (``sys.stdout`` of the calling thread).

    Re-entrant and thread-aware; a no-op when ``ETNA_DEVSIM_VERBOSE=1``. No
    effect on results: it only discards text.
    """
    if _verbose():
        yield
        return
    tid = threading.get_ident()
    with _PROXY_LOCK:
        proxy = sys.stdout
        if not isinstance(proxy, _ThreadQuietStdout):
            proxy = _ThreadQuietStdout(sys.stdout)
            sys.stdout = proxy
        proxy._quiet[tid] = proxy._quiet.get(tid, 0) + 1
    try:
        yield
    finally:
        with _PROXY_LOCK:
            depth = proxy._quiet.get(tid, 1) - 1
            if depth:
                proxy._quiet[tid] = depth
            else:
                proxy._quiet.pop(tid, None)
            # uninstall once nobody is quiet, unless someone replaced sys.stdout since
            if not proxy._quiet and sys.stdout is proxy:
                sys.stdout = proxy._wrapped


def serialized(fn):
    """Run ``fn`` while holding the devsim lock, with devsim's log silenced."""

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        with DEVSIM_LOCK, quiet_devsim():
            return fn(*args, **kwargs)

    return wrapper
