from etna import DeviceConfig
from etna.core.trap_kinetics import KineticTrapSolver
from etna.core.trap_presets import trap_preset

cfg = DeviceConfig(epi_thickness_um=50.0, doping_profile="uniform", N_D=4e13)
for spec in ("auto", "radmanovac_calibrated"):
    s = KineticTrapSolver.from_config(
        cfg, [1e12], 252.7, trap_model=spec, niel_convention="si-neq"
    )
    print(
        spec,
        "si-neq run: model convention",
        s.model.niel_convention,
        "d",
        s.model.donor_removal,
    )
# what it should be
print("expected", trap_preset("radmanovac_calibrated", "si-neq").donor_removal)
# explicit object in wrong convention raises?
try:
    KineticTrapSolver.from_config(
        cfg,
        [1e12],
        252.7,
        trap_model=trap_preset("radmanovac_calibrated"),
        niel_convention="si-neq",
    )
    print("object: no raise")
except ValueError as e:
    print("object raises OK")
# calibrate_cv with string trap_model
import inspect
from etna.core import trap_calibration as tc

print(
    "calibrate_cv resolve without conv:",
    "resolve_traps(trap_model)" in inspect.getsource(tc.calibrate_cv),
)
from etna.core.radiation_traps import resolve_traps
from etna.core.radiation_traps import check_trap_convention

m = resolve_traps("auto")
print("calibrate_cv('auto', si-neq) check would raise:", end=" ")
try:
    check_trap_convention(m, "si-neq")
    print("no")
except ValueError:
    print("yes")
