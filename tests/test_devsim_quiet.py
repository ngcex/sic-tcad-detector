"""quiet_devsim: facades do not flood stdout; ETNA_DEVSIM_VERBOSE=1 restores it; no numeric effect."""

import io
import sys
import threading

import numpy as np

from etna.api.device import DeviceConfig
from etna.api.simulation import run_field
from etna.core.devsim_lock import VERBOSE_ENV, quiet_devsim


def _field(capsys):
    res = run_field(DeviceConfig(), bias_V=-20.0, solver="devsim")
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


_probe = iter(range(10**6))


def _direct_core_field():
    """A direct core call (no facade, no ``serialized``): builds and solves devsim."""
    import devsim

    from etna.core.charge_collection_2d import create_2d_dd_device

    info = create_2d_dd_device(
        half_width_um=20.0, V_bias=5.0, device_name=f"quiet_probe{next(_probe)}"
    )
    try:
        return np.array(
            devsim.get_node_model_values(
                device=info["device_name"], region=info["region_name"], name="Potential"
            )
        )
    finally:
        devsim.delete_device(device=info["device_name"])


def test_set_devsim_quiet_silences_direct_core_calls(capsys, monkeypatch):
    import devsim

    from etna.core.devsim_lock import devsim_quiet_enabled, set_devsim_quiet

    monkeypatch.delenv(VERBOSE_ENV, raising=False)
    original_solve = devsim.solve
    loud = _direct_core_field()
    out_loud = capsys.readouterr().out
    assert out_loud.count("\n") > 20  # direct calls are not covered by the facades

    set_devsim_quiet(True)
    try:
        set_devsim_quiet(True)  # idempotent
        assert devsim_quiet_enabled()
        quiet = _direct_core_field()
        print("caller-print-survives")
        out_quiet = capsys.readouterr().out
        assert out_quiet == "caller-print-survives\n", out_quiet[:500]

        # from other threads too
        box = {}
        t = threading.Thread(target=lambda: box.setdefault("v", _direct_core_field()))
        t.start()
        t.join()
        assert capsys.readouterr().out == ""

        # the verbose override still wins
        monkeypatch.setenv(VERBOSE_ENV, "1")
        forced = _direct_core_field()
        assert capsys.readouterr().out.count("\n") > 20
        monkeypatch.delenv(VERBOSE_ENV)
        leaked = devsim.solve  # a copy taken while on (``from devsim import *``)
    finally:
        set_devsim_quiet(False)
    assert not devsim_quiet_enabled()
    assert devsim.solve is original_solve
    # off means off everywhere: a wrapper copied while on passes devsim's log through
    from etna.core import devsim_lock

    assert not devsim_lock._enabled
    assert leaked.__etna_quiet_original__ is original_solve
    monkeypatch.setattr(devsim, "solve", leaked)
    capsys.readouterr()
    after = _direct_core_field()
    assert capsys.readouterr().out.count("\n") > 20
    np.testing.assert_array_equal(after, loud)

    # zero numeric impact
    np.testing.assert_array_equal(quiet, loud)
    np.testing.assert_array_equal(box["v"], loud)
    np.testing.assert_array_equal(forced, loud)
