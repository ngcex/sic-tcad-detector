import importlib.util

import numpy as np

from etna.core import dark_current as new
from etna.core.devsim_lock import quiet_devsim
from etna.core.devsim_reset import reset_devsim_fully

# pre-change module: git show 2f86879:etna/core/dark_current.py > /tmp/wpb/dc_old.py
spec = importlib.util.spec_from_file_location("dc_old", "/tmp/wpb/dc_old.py")
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)

lines = []
fl = np.array([0.0, 1e12, 1e13, 2e13])
with quiet_devsim():
    for V in (-30.0, -80.0):
        a = old.dark_current_vs_fluence(fl, V_bias=V, area=0.04)
        b = new.dark_current_vs_fluence(fl, V_bias=V, area=0.04)
        lines.append(f"FL V={V} valid old {a['valid']} new {b['valid']}")
        lines.append(f"FL   old I {a['I_total']}")
        lines.append(f"FL   rel   {b['I_total'] / a['I_total'] - 1}")
    for args in ((1e13, 600.0, 1800.0), (2e13, 400.0, 600.0)):
        a = old.dark_current_post_anneal(*args, area=0.04)
        b = new.dark_current_post_anneal(*args, area=0.04)
        lines.append(
            f"AN {args} valid {a['valid']} {b['valid']} I {a['I_total']:.4e} rel {b['I_total'] / a['I_total'] - 1:.2e}"
        )
    for pn, vals in (("N_t", [1e12, 5e12]), ("epi_thickness_cm", [8e-4, 12e-4])):
        for T in (300.0, 400.0):
            reset_devsim_fully()
            a = old.sensitivity_sweep(pn, vals, V_eval=-30.0, T=T, area=0.04)
            reset_devsim_fully()
            b = new.sensitivity_sweep(pn, vals, V_eval=-30.0, T=T, area=0.04)
            lines.append(
                f"SE {pn} T={T} I {a['I_total'].to_numpy()} rel {(b['I_total'] / a['I_total'] - 1).to_numpy()}"
            )
print("\n".join(lines))
