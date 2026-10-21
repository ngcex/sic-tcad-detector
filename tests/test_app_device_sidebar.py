"""Unit tests for the pure app.components.device_sidebar.assemble_config seam.

UI-02 coverage: config assembly (all 11 DeviceConfig fields), doping-mode
mapping (graded <-> uniform / N_D consistency), and dimensionality mapping
(1D <-> 2D / half_width_um consistency). No Streamlit runtime is required —
assemble_config is a pure function, imported directly.
"""

from etna import DeviceConfig
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
    "nplus_thickness_um": 2.0,  # added with the n+ layer DeviceConfig fields
    "N_D_nplus": 5e18,
    "N_t": 3e12,  # dark-current fit fields, added with the n+ default flip
    "E_t": 0.66,
    "buffer_thickness_um": 1.0,  # graded n-/n+ transition (opt-in measured-C-V preset)
    "buffer_end_doping": 1e17,
    # opt-in tabulated profile: None here (a table needs the graded profile and excludes
    # the transition above; its round trip is tests/test_measured_cv_tabulated.py)
    "donor_table": None,
    # opt-in A-02 model: the sidebar has no toggle yet (app/ not touched), so the
    # round-trip is the default (legacy) value
    "incomplete_ionization": False,
    # opt-in charged radiation traps: set by the Radiation Damage page's own checkbox,
    # not the sidebar, so the round-trip is the default
    "radiation_traps": False,
    # low-field mobility set (sidebar selectbox; non-default value here)
    "mobility": "legacy_perp",
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


def test_historical_model_checkbox_applies_the_legacy_preset():
    """The sidebar's 'historical model' checkbox merges LEGACY_PRESET into the values."""
    from etna import LEGACY_PRESET, legacy_config

    values = dict(ALL_FIELDS_VALUES, doping_profile="graded", N_D=None)
    values.update(LEGACY_PRESET)
    cfg = assemble_config(values)
    ref = legacy_config(
        epi_thickness_um=5.0, half_width_um=50.0, N_A=2e19, T=310.0, area_cm2=2e-4
    )
    assert cfg == ref
    assert cfg.nplus_thickness_um is None and (cfg.N_t, cfg.E_t) == (2.8e12, 0.65)


def test_sidebar_default_is_the_nplus_stack():
    from app.components.device_sidebar import render_device_sidebar  # noqa: F401
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file("app/main.py", default_timeout=120)
    at.run()
    # the sidebar models the 4 mm2 dosimeter == DeviceConfig's default area (2026-10-18;
    # before, the sidebar overrode DeviceConfig's 1e-4 cm2 with 0.04)
    import dataclasses

    assert at.session_state["device_config"] == DeviceConfig()
    assert not at.sidebar.checkbox(key="cfg_legacy_on").value
    texts = " ".join(c.value for c in at.sidebar.caption)
    assert "Inside the 1D engine/devsim cross-check domain" in texts
    at.sidebar.checkbox(key="cfg_legacy_on").check()
    at.run()
    from etna import legacy_config, measured_cv_config, nplus_calibrated_config

    assert at.session_state["device_config"] == legacy_config()
    # the default IS the measured-C-V profile (2026-10-18); the previous calibrated model
    # (abrupt n+, no transition) is the opt-in checkbox now
    at.sidebar.checkbox(key="cfg_legacy_on").uncheck()
    at.run()
    assert at.session_state["device_config"] == measured_cv_config()
    assert not at.sidebar.checkbox(key="cfg_nplus_cal_on").value
    at.sidebar.checkbox(key="cfg_nplus_cal_on").check()
    at.run()
    assert at.session_state["device_config"] == nplus_calibrated_config()
    infos = " ".join(c.value for c in at.sidebar.info)
    assert "Previous calibrated model" in infos and "no graded transition" in infos
    texts = " ".join(c.value for c in at.sidebar.caption)
    assert "Inside the 1D engine/devsim cross-check domain" in texts
    # the transition checkbox off = abrupt epi end with the default doping
    at.sidebar.checkbox(key="cfg_nplus_cal_on").uncheck()
    at.sidebar.checkbox(key="cfg_buffer_on").uncheck()
    at.run()
    assert at.session_state["device_config"] == dataclasses.replace(
        DeviceConfig(), buffer_thickness_um=None, buffer_end_doping=None
    )
