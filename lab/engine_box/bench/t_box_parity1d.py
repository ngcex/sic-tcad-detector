"""Parity engine vs devsim for run_cv / run_cce outside the 250-400 K box."""

import sys
import time
from dataclasses import replace
import numpy as np
from etna.api import simulation as sim
from etna.api.device import (
    DeviceConfig,
    legacy_config,
    nplus_calibrated_config,
    microdosimeter_config,
)

sim._ENGINE_BOX["T"] = (50.0, 1000.0)
temps = [float(t) for t in sys.argv[1].split(",")]
cfgs = [
    ("default", DeviceConfig()),
    ("nplus_cal", nplus_calibrated_config()),
    ("legacy", legacy_config()),
    ("micro", microdosimeter_config()),
]
from etna.core.devsim_lock import quiet_devsim
_q = quiet_devsim(); _q.__enter__()
import builtins, sys as _s
_out = _s.stdout._wrapped
def pr(*a, **k):
    k.pop("flush", None); _out.write(" ".join(map(str, a)) + "\n"); _out.flush()
for name, base in cfgs:
    for T in temps:
        cfg = replace(base, T=T)
        vstop = -40.0 if name == "legacy" else -100.0
        bias = np.linspace(0.0, vstop, 21)
        t = time.perf_counter()
        e = sim.run_cv(cfg, 0.0, vstop, 21)
        te = time.perf_counter() - t
        t = time.perf_counter()
        try:
            d = sim._run_cv_devsim(cfg, bias)
        except Exception as ex:
            pr(f"{name} T={T} cv devsim failed {ex}")
            continue
        td = time.perf_counter() - t
        n = len(d.x)
        dC = np.max(np.abs(e.y[:n] / d.y - 1))
        dW = np.max(
            np.abs(e.metadata["depletion_widths"][:n] - d.metadata["depletion_widths"])
        )
        cb = np.linspace(-10.0, vstop, 7)
        ec = sim.run_cce(cfg, -10.0, vstop, 7)
        dc = sim._run_cce_devsim(cfg, cb)
        m = len(dc.x)
        dcce = np.max(np.abs(ec.y[:m] - dc.y)) if m else float("nan")
        pr(
            f"{name:9s} T={T:5.0f} solver={e.metadata['solver']:8s}/{ec.metadata['solver']:8s} "
            f"devsim pts {n}/{len(bias)} cce {m}/7 |dC/C|={dC:.1e} |dW|={dW:.1e} cm |dCCE|={dcce:.1e} "
            f"t_eng {te:.2f}s t_dev {td:.2f}s",
            flush=True,
        )
