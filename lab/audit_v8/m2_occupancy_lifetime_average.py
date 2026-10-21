"""Default damage path (auto -> radmanovac_calibrated, lifetimes='occupancy'): occupancy
lifetimes are a NODE-COUNT mean over x >= junction (epi + n+ buffer + substrate).
Compare with a length-weighted mean over the epi only, then the CCE impact."""

import numpy as np

import etna.core.radiation_traps as rt
from etna import DeviceConfig, run_radiation_damage
from etna.api.device import device_kwargs
from etna.core.charge_collection import _fluence_reference, _trap_dark_states_kinetic
from etna.core.radiation_damage import get_hardness_factor, hardness_table
from etna.core.sic_material import SiC4H_Parameters
from etna.core.radiation_traps import (
    resolve_traps,
    trap_fields,
    level_table,
    capture_availability,
)

cfg = DeviceConfig()
model = resolve_traps("auto")
print("model:", model.label[:60], "lifetimes", model.lifetimes)
ref = _fluence_reference(device_kwargs(cfg))
x = ref["x"]
jp = ref["info"]["junction_pos"]
xe = ref["info"]["epi_end"]
kappa = get_hardness_factor(62.0, hardness_table("proton", "sic-damage"))
Fp = np.array([1e13, 1e14, 1e15])
dark = _trap_dark_states_kinetic(ref, Fp * kappa, -40.0, model)
st = dark["state"]
T = float(ref["info"]["T"])
prm = SiC4H_Parameters()
vn = rt._vth(T, prm.m_e_dos)
rows = level_table(model)
sn = np.array([r["sigma_n"] for r in rows])[:, None]
mask = ref["epi_mask"]
epi_only = (x >= jp) & (x <= xe)
# control-volume weights
h = np.diff(x)
w = np.zeros_like(x)
w[1:] += 0.5 * h
w[:-1] += 0.5 * h
for b, F in enumerate(Fp * kappa):
    fl = trap_fields(model, F, mask, T, float(ref["info"]["n_i"]))
    lo, up = capture_availability(fl, st.f[b])
    rate = (fl.N_t * sn * vn * lo).sum(axis=0)  # 1/s per node
    node_mean_mask = rate[mask].mean()
    lw_epi = (rate * w)[epi_only].sum() / w[epi_only].sum()
    lw_mask = (rate * w)[mask].sum() / w[mask].sum()
    print(
        f"Phi_p {Fp[b]:.0e}: 1/tau_n trap part: node-mean(mask, code) {node_mean_mask:.3e}  "
        f"length-weighted(mask) {lw_mask:.3e}  length-weighted(epi only) {lw_epi:.3e}  "
        f"ratio epi/code {lw_epi / node_mean_mask:.3f}"
    )

# CCE impact: patch occupancy_lifetimes to length-weight over the epi only
orig = rt.occupancy_lifetimes
fl_ = Fp
base = run_radiation_damage(cfg, fluences=Fp, V_bias=-40.0, proton_energy_MeV=62.0)


def patched(model, fields, f, mask_, T, tau0_n, tau0_p, x=None):  # x: post-#13 call
    prm = SiC4H_Parameters()
    vn_, vp_ = rt._vth(T, prm.m_e_dos), rt._vth(T, prm.m_h_dos)
    rows_ = level_table(model)
    sn_ = np.array([r["sigma_n"] for r in rows_])[:, None]
    sp_ = np.array([r["sigma_p"] for r in rows_])[:, None]
    lo_, up_ = capture_availability(fields, f)
    rn = ((fields.N_t * sn_ * vn_ * lo_).sum(axis=0) * w)[epi_only].sum() / w[
        epi_only
    ].sum()
    rp = ((fields.N_t * sp_ * vp_ * up_).sum(axis=0) * w)[epi_only].sum() / w[
        epi_only
    ].sum()
    return 1.0 / (1.0 / tau0_n + rn), 1.0 / (1.0 / tau0_p + rp)


rt.occupancy_lifetimes = patched
alt = run_radiation_damage(cfg, fluences=Fp, V_bias=-40.0, proton_energy_MeV=62.0)
rt.occupancy_lifetimes = orig
print("CCE code   :", np.round(base.y, 6), base.metadata["damage_model"])
print("CCE epi-lw :", np.round(alt.y, 6))
print("1-CCE ratio:", np.round((1 - alt.y) / (1 - base.y), 3))
