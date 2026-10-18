"""quiet_devsim: facades do not flood stdout; ETNA_DEVSIM_VERBOSE=1 restores it; no numeric effect."""

import io
import sys
import threading

import numpy as np

from etna.api.device import DeviceConfig
from etna.api.simulation import run_field
from etna.core.devsim_lock import VERBOSE_ENV, quiet_devsim


def _field(capsys):
    res = run_field(DeviceConfig(), bias_V=-20.0)
    return res, capsys.readouterr().out


def test_facade_is_quiet_and_verbose_env_restores_output(capsys, monkeypatch):
    monkeypatch.delenv(VERBOSE_ENV, raising=False)
    quiet, out_quiet = _field(capsys)
    assert out_quiet.count("\n") == 0, out_quiet[:500]

    monkeypatch.setenv(VERBOSE_ENV, "1")
    loud, out_loud = _field(capsys)
    assert out_loud.count("\n") > 50  # devsim's Newton log

    # silencing only discards text
    np.testing.assert_array_equal(quiet.x, loud.x)
    np.testing.assert_array_equal(quiet.y, loud.y)


def test_quiet_devsim_is_thread_local_and_restores_stdout(monkeypatch):
    monkeypatch.delenv(VERBOSE_ENV, raising=False)
    sink = io.StringIO()
    monkeypatch.setattr(sys, "stdout", sink)
    other_ready, other_done = threading.Event(), threading.Event()

    def other():
        other_ready.wait()
        print("from-other-thread")
        other_done.set()

    t = threading.Thread(target=other)
    t.start()
    with quiet_devsim():
        with quiet_devsim():  # re-entrant
            print("dropped")
        print("dropped-too")
        other_ready.set()
        other_done.wait()
    t.join()
    print("after")
    assert sys.stdout is sink
    assert sink.getvalue() == "from-other-thread\nafter\n"
