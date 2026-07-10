"""Unit tests for the pure app.components.device_sidebar.assemble_config seam.

UI-02 coverage: config assembly (all 11 DeviceConfig fields), doping-mode
mapping (graded <-> uniform / N_D consistency), and dimensionality mapping
(1D <-> 2D / half_width_um consistency). No Streamlit runtime is required —
assemble_config is a pure function, imported directly.
"""

from petringa import DeviceConfig
from app.components.device_sidebar import assemble_config

# Non-default values for all 11 DeviceConfig fields, uniform-profile mode
# (2D). Chosen to differ from every DeviceConfig field default so the
# Task-2 stub (which returns DeviceConfig()) fails these assertions.
ALL_FIELDS_VALUES = {
    "epi_thickness_um": 5.0,
    "substrate_thickness_um": 2.0,
    "half_width_um": 50.0,
    "N_A": 2e19,
    "doping_profile": "uniform",
    "N_D": 1e15,
    "N_D_junction": 3.5e15,
    "N_D_bulk": 9.0e13,
    "L_transition_um": 1.5,
    "T": 310.0,
    "area_cm2": 2e-4,
}


def test_assemble_config_all_fields():
    result = assemble_config(ALL_FIELDS_VALUES)

    assert isinstance(result, DeviceConfig)
    for field_name in DeviceConfig.__dataclass_fields__:
        assert getattr(result, field_name) == ALL_FIELDS_VALUES[field_name], (
            f"field {field_name!r} did not round-trip: "
            f"expected {ALL_FIELDS_VALUES[field_name]!r}, "
            f"got {getattr(result, field_name)!r}"
        )


def test_doping_mode_mapping():
    graded_values = dict(ALL_FIELDS_VALUES)
    graded_values["doping_profile"] = "graded"
    graded_values["N_D"] = None
    graded_values["N_D_junction"] = 3.5e15
    graded_values["N_D_bulk"] = 9.0e13
    graded_values["L_transition_um"] = 1.5

    graded_result = assemble_config(graded_values)
    assert graded_result.doping_profile == "graded"
    assert graded_result.N_D is None
    assert graded_result.N_D_junction == 3.5e15
    assert graded_result.N_D_bulk == 9.0e13
    assert graded_result.L_transition_um == 1.5

    uniform_values = dict(ALL_FIELDS_VALUES)
    uniform_values["doping_profile"] = "uniform"
    uniform_values["N_D"] = 1e15

    uniform_result = assemble_config(uniform_values)
    assert uniform_result.doping_profile == "uniform"
    assert isinstance(uniform_result.N_D, float)
    assert uniform_result.N_D == 1e15


def test_dimensionality_mapping():
    values_1d = dict(ALL_FIELDS_VALUES)
    values_1d["half_width_um"] = None

    result_1d = assemble_config(values_1d)
    assert result_1d.half_width_um is None

    values_2d = dict(ALL_FIELDS_VALUES)
    values_2d["half_width_um"] = 50.0

    result_2d = assemble_config(values_2d)
    assert result_2d.half_width_um == 50.0
