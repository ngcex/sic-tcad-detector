import numpy as np
from etna.core import trap_calibration as tc
from etna.core.trap_presets import literature_kinetic

obs = [
    tc.RemovalObservation(
        "a", "slope", 5.0, 0.5, 4e13, fluences=(0, 1e12), energy_MeV=None
    ),
    tc.RemovalObservation(
        "b", "slope", 4.0, 0.5, 4e13, fluences=(0, 1e12), energy_MeV=None
    ),
    tc.RemovalObservation(
        "c", "slope", 4.5, 0.5, 4e13, fluences=(0, 1e12), energy_MeV=None
    ),
]


def fake(o, m):  # linear in donor_removal; observation "b" unmeasurable (NaN)
    if o.name == "b":
        return np.nan
    return 2.0 + m.donor_removal


tc.removal_observable = fake
cal = tc.calibrate_removal(obs, ("donor_removal",), base_model=literature_kinetic())
print("d", cal.scales, "+-", cal.sigma, "chi2", cal.chi2, "dof", cal.dof)
print("model_values", cal.settings["model_values"])
print("warnings", cal.warnings)
