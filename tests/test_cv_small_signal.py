"""C-V observable: run_cv returns the small-signal C = dQ/dV on both paths.

Batched engine: exact linearisation of the flat-quasi-Fermi Poisson state
(PoissonBatch.small_signal_capacitance). devsim: central difference (+-1 mV) of the
terminal hole charge on the DD state (cv_analysis.small_signal_capacitance). They are
independent and must agree to <= 1e-6 (measured <= 3e-8). A devsim small-signal AC
solve pins both (own process); gate: lab/microdosimeter/bench/cv_gate.py
(lab/microdosimeter/notes/LOG.md section 6).
"""

import devsim
import numpy as np
import pytest

from etna.api import simulation as sim
from etna.api.device import (
    DeviceConfig,
    build_device,
    legacy_config,
    microdosimeter_config,
)
from etna.core.batched import make_batch, node_positions
from etna.core.batched.material import EPS, Q
from etna.core.batched.sweeps import _Ramp
from etna.core.cv_analysis import small_signal_capacitance
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.drift_diffusion import ramp_bias, extract_steady_state_current
from etna.core.poisson import extract_depletion_width_numerical

CONFIGS = {
    "default": DeviceConfig(),
    "legacy_no_nplus": legacy_config(),
    "uniform": DeviceConfig(doping_profile="uniform", N_D=1e15),
    "T250": DeviceConfig(T=250.0),
    "T400": DeviceConfig(T=400.0),
    "microdosimeter": microdosimeter_config(),
}


@pytest.mark.parametrize("name", CONFIGS)
def test_engine_and_devsim_ac_capacitance_agree(name):
    cfg = CONFIGS[name]
    eng = sim.run_cv(cfg, 0.0, -20.0, 9)
    assert eng.metadata["solver"] == "batched"
    dev = sim._run_cv_devsim(cfg, np.linspace(0.0, -20.0, 9))
    assert dev.metadata["solver"] == "devsim"
    np.testing.assert_allclose(eng.y, dev.y, rtol=1e-6)
    np.testing.assert_allclose(
        eng.metadata["depletion_widths"], dev.metadata["depletion_widths"], rtol=1e-8
    )
    np.testing.assert_allclose(
        eng.metadata["W_from_C"], EPS * cfg.area_cm2 / eng.y, rtol=1e-12
    )
    # C is not eps*A/W_edge: the small-signal width is never below the edge
    assert np.all(eng.metadata["W_from_C"] >= eng.metadata["depletion_widths"] * 0.999)


def test_engine_capacitance_is_the_terminal_charge_derivative():
    """Analytic dQ/dV == central difference of the hole charge q*sum(v p) (plane-free)."""
    cfg = DeviceConfig()
    sub, epi, t_np = 0.3e-4, 10e-4, 1e-4
    x = node_positions(epi, sub, t_np)
    pb = make_batch(
        x,
        sub,
        cfg.N_A,
        cfg.N_D_junction,
        cfg.N_D_bulk,
        cfg.L_transition_um * 1e-4,
        nplus=(sub + epi, cfg.N_D_nplus),
    )
    v = pb.v.copy()
    v[0] = v[-1] = 0.0
    ramp = _Ramp(pb, 1e-12)
    dv = 0.05
    for Vc in (0.0, 10.0, 30.0):
        r_lo = ramp.to(Vc - dv)
        r_mid = ramp.to(Vc)
        C = pb.small_signal_capacitance(r_mid)[0]
        r_hi = ramp.to(Vc + dv)
        dQ = Q * np.sum(v * (r_lo.p[0] - r_hi.p[0]))
        assert C == pytest.approx(dQ / (2 * dv), rel=2e-5)


def test_small_signal_solve_leaves_devsim_state_intact():
    reset_devsim_fully()
    info = build_device(DeviceConfig())
    try:
        ramp_bias(info, 10.0, contact="cathode", V_step=0.5)
        I0 = extract_steady_state_current(info)
        W0 = extract_depletion_width_numerical(info)
        small_signal_capacitance(info)
        assert devsim.get_parameter(
            device=info["device_name"], name="cathode_bias"
        ) == pytest.approx(10.0, abs=0)
        assert extract_depletion_width_numerical(info) == W0
        ramp_bias(info, 12.0, contact="cathode", V_step=0.5)
        ramp_bias(info, 10.0, contact="cathode", V_step=0.5)
        assert extract_steady_state_current(info) == pytest.approx(I0, rel=1e-6)
        assert extract_depletion_width_numerical(info) == pytest.approx(W0, rel=1e-10)
    finally:
        reset_devsim_fully()


def test_failed_side_solve_restores_the_state(monkeypatch):
    reset_devsim_fully()
    info = build_device(DeviceConfig())
    real_solve = devsim.solve
    calls = {"n": 0}

    def solve(**kw):
        calls["n"] += 1
        if calls["n"] >= 2:  # the -step side solve and its relaxed retry
            raise devsim.error("forced failure")
        return real_solve(**kw)

    try:
        ramp_bias(info, 5.0, contact="cathode", V_step=0.5)
        W0 = extract_depletion_width_numerical(info)
        monkeypatch.setattr(devsim, "solve", solve)
        with pytest.raises(devsim.error):
            small_signal_capacitance(info)
        monkeypatch.setattr(devsim, "solve", real_solve)
        assert devsim.get_parameter(
            device=info["device_name"], name="cathode_bias"
        ) == pytest.approx(5.0, abs=0)
        assert extract_depletion_width_numerical(info) == W0
        assert np.isfinite(small_signal_capacitance(info))
    finally:
        reset_devsim_fully()


def test_devsim_cv_does_not_perturb_later_devsim_runs():
    cfg = DeviceConfig()
    before = sim.run_dark_current(cfg, 0.0, -20.0, 3)
    sim._run_cv_devsim(cfg, np.array([0.0, -5.0]))
    after = sim.run_dark_current(cfg, 0.0, -20.0, 3)
    np.testing.assert_allclose(after.y, before.y, rtol=1e-9)


_AC_SNIPPET = """
import json, os, sys
sys.path.insert(0, os.path.join("lab", "microdosimeter", "bench"))
out = os.dup(1)
os.dup2(os.open(os.devnull, os.O_WRONLY), 1)  # devsim's C-level chatter
import cv_gate
C = cv_gate.devsim_ac(volts=(0.0, -10.0, -30.0), freqs=(1.0,))
os.write(out, ("AC " + json.dumps({str(k): v for k, v in C.items()})).encode())
"""


def test_devsim_ac_solve_agrees_with_both_paths():
    """Independent pin: devsim's small-signal AC solve (in its own process: devsim 2.10
    cannot remove a circuit again) == engine and devsim-path dQ/dV."""
    import json
    import pathlib
    import subprocess
    import sys

    root = pathlib.Path(__file__).resolve().parents[1]
    r = subprocess.run(
        [sys.executable, "-c", _AC_SNIPPET],
        cwd=root,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert r.returncode == 0, r.stderr[-2000:]
    ac = json.loads(r.stdout.split("AC ", 1)[1])
    cfg = DeviceConfig()
    eng = sim.run_cv(cfg, 0.0, -30.0, 4)  # 0 / -10 / -20 / -30 V
    for i, V in ((0, "0.0"), (1, "-10.0"), (3, "-30.0")):
        assert eng.y[i] / cfg.area_cm2 == pytest.approx(ac[V], rel=1e-6)
