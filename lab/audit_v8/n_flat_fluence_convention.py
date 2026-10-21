from etna.core.trap_calibration import RemovalObservation, removal_observable
from etna.core.trap_presets import trap_preset

o = RemovalObservation(
    "flat", "flat_fluence", 12.7, 0.3, 4e13, energy_MeV=252.7, freq_Hz=1e6
)  # default convention sic-damage
m_si = trap_preset("radmanovac_calibrated", "si-neq")
m_c = trap_preset("radmanovac_calibrated")
print(
    "si-neq-calibrated model under sic-damage obs (should raise):",
    removal_observable(o, m_si),
)
print("matching model:", removal_observable(o, m_c))
