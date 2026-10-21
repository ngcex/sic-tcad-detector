from etna.core.radiation_damage import (
    RadiationDamageParams,
    get_hardness_factor,
    hardness_table,
)

eta = RadiationDamageParams().eta_removal
for c in ("si-neq", "sic-damage", "sic-1mev"):
    k = get_hardness_factor(252.7, hardness_table("proton", c))
    k62 = get_hardness_factor(62.0, hardness_table("proton", c))
    print(
        f"{c:10s} legacy removal per proton at 252.7 MeV: {eta * k:.2f} /p (Radmanovac 4.48-5.60);"
        f" eta needed for 4.48-5.60 /p: {4.48 / k:.2f}-{5.60 / k:.2f} /neq; 62 MeV: {eta * k62:.2f} /p"
    )
