from etna.api.device import legacy_config
from etna.api.simulation import damage_model_of, compensation_fluence

cfg = legacy_config()
print("radiation_traps", cfg.radiation_traps, "nplus", cfg.nplus_thickness_um)
for V in (0.0, -20.0, -60.0):
    m, note = damage_model_of(cfg, 62.0, V_bias=V)
    print(
        V,
        "model:",
        None if m is None else m.label[:30],
        "| note:",
        note[:90],
        "| phi_c",
        compensation_fluence(cfg, 62.0, V_bias=V),
    )
