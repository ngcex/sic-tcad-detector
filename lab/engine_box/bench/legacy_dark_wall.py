import importlib.util

import devsim
import numpy as np

from etna.api.device import device_kwargs, legacy_config
from etna.core import dark_current as new
from etna.core.devsim_lock import quiet_devsim
from etna.core.devsim_reset import reset_devsim_fully

spec = importlib.util.spec_from_file_location("dc_old", "/tmp/wpb/dc_old.py")
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
cfg = legacy_config()
bias = np.linspace(0, -100, 20)
out = {}
for label, mod, step in (
    ("old0.5", old, 0.5),
    ("new0.5", new, 0.5),
    ("new2", new, 2.0),
    ("old2", old, 2.0),
):
    with quiet_devsim():
        reset_devsim_fully()
        info = mod.create_dark_current_device(
            N_t=cfg.N_t, E_t=cfg.E_t, **device_kwargs(cfg)
        )
        r = mod.dark_current_sweep(info, bias, area=cfg.area_cm2, V_step=step)
        devsim.delete_device(device=info["device_name"])
    out[label] = r
    print("L", label, "reached", r["reached_v_stop"], "n", len(r["voltages"]))
n = min(len(out["old0.5"]["I_total"]), len(out["new2"]["I_total"]))
print(
    "L rel new2 vs old0.5 on common pts",
    np.max(np.abs(out["new2"]["I_total"][1:n] / out["old0.5"]["I_total"][1:n] - 1)),
)
print("L new2 I", " ".join(f"{v:.4e}" for v in out["new2"]["I_total"]))
