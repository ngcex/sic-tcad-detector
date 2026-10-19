"""Concurrent facade calls must not crash the process (devsim is process-global)."""

import subprocess
import sys
import textwrap

import pytest

pytest.importorskip("devsim")

_SCRIPT = textwrap.dedent(
    """
    import threading
    from etna.api.device import DeviceConfig
    from etna.api.simulation import run_cce, run_cv

    errors = []
    solvers = set()

    def work(fn, cfg=None, **kw):
        try:
            res = fn(cfg or DeviceConfig(), **kw)
            solvers.add(res.metadata["solver"])
        except BaseException as exc:  # noqa: BLE001
            errors.append(repr(exc))

    def work_hot(fn, **kw):
        work(fn, cfg=DeviceConfig(T=500.0), **kw)

    # Default configs take the numba engine (aborts the process if several threads
    # enter it at once); DeviceConfig(T=500) is outside its domain and takes devsim.
    threads = [
        threading.Thread(target=work, args=(run_cv,), kwargs=dict(n_points=3)),
        threading.Thread(target=work, args=(run_cv,), kwargs=dict(n_points=3)),
        threading.Thread(target=work, args=(run_cv,), kwargs=dict(n_points=3)),
        threading.Thread(target=work, args=(run_cv,), kwargs=dict(n_points=3)),
        threading.Thread(target=work, args=(run_cce,), kwargs=dict(n_points=3)),
        threading.Thread(target=work_hot, args=(run_cv,), kwargs=dict(n_points=3)),
        threading.Thread(target=work_hot, args=(run_cv,), kwargs=dict(n_points=3)),
        threading.Thread(target=work_hot, args=(run_cce,), kwargs=dict(n_points=3)),
    ]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not errors, errors
    assert solvers == {"batched", "devsim"}, solvers
    print("OK")
    """
)


@pytest.mark.slow
def test_concurrent_facades_do_not_crash():
    # Subprocess: a regression is a SIGSEGV, which would kill pytest itself.
    proc = subprocess.run(
        [sys.executable, "-c", _SCRIPT], capture_output=True, text=True, timeout=300
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "OK" in proc.stdout
