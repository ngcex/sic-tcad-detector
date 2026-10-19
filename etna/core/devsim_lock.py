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
notebooks or scripts, or call :func:`set_devsim_quiet` once (first notebook
cell, script ``main``) to silence every devsim call of the process, from any
thread: each public ``devsim`` function is then wrapped so that only the time
spent inside devsim is quiet; ``print`` from the notebook itself still shows.

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


_GLOBAL_LOCK = threading.Lock()
_patched: list[tuple[object, str, object]] = []  # (module, name, original)
_enabled = False  # read by the wrappers at call time


def _quiet_call(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        if (
            not _enabled
        ):  # a copy taken while on (e.g. a later ``from devsim import *``)
            return fn(*args, **kwargs)
        with quiet_devsim():
            return fn(*args, **kwargs)

    wrapper.__etna_quiet_original__ = fn
    return wrapper


def set_devsim_quiet(enabled: bool = True) -> None:
    """Process-wide switch: silence devsim's log for every devsim call, any thread.

    ``set_devsim_quiet(True)`` wraps each public function of the ``devsim``
    module (and the aliases bound by ``from devsim import *`` in the already
    imported ``devsim.*`` packages) in :func:`quiet_devsim`, so direct core
    calls from notebooks, scripts and labs stop streaming the Newton log, while
    the caller's own ``print`` output between devsim calls is unaffected.
    ``set_devsim_quiet(False)`` restores the original functions (the default
    state); wrappers copied while it was on (a module imported later with
    ``from devsim import *``) then pass through. ``ETNA_DEVSIM_VERBOSE=1`` still forces devsim's output (checked at
    call time). Idempotent; no effect on results: it only discards text.
    """
    global _enabled
    import devsim

    with _GLOBAL_LOCK:
        _enabled = bool(enabled)
        if enabled and not _patched:
            originals = {
                name: obj
                for name, obj in vars(devsim).items()
                if not name.startswith("_")
                and callable(obj)
                and not isinstance(obj, type)
            }
            by_id = {id(obj): name for name, obj in originals.items()}
            wrappers = {name: _quiet_call(obj) for name, obj in originals.items()}
            modules = [devsim] + [
                m
                for k, m in list(sys.modules.items())
                if k.startswith("devsim.") and m is not None
            ]
            for mod in modules:
                for name, obj in list(vars(mod).items()):
                    src = by_id.get(id(obj))
                    if src is not None and originals[src] is obj:
                        _patched.append((mod, name, obj))
                        setattr(mod, name, wrappers[src])
        elif not enabled and _patched:
            for mod, name, obj in reversed(_patched):
                setattr(mod, name, obj)
            _patched.clear()


def devsim_quiet_enabled() -> bool:
    """True while :func:`set_devsim_quiet` (process-wide) is on."""
    return _enabled


def serialized(fn):
    """Run ``fn`` while holding the devsim lock, with devsim's log silenced."""

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        with DEVSIM_LOCK, quiet_devsim():
            return fn(*args, **kwargs)

    return wrapper
