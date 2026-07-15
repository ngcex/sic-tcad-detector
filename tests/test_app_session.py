"""UI-07 persistence-key contract test.

Confirms that a DeviceConfig assembled by assemble_config, placed in a plain
dict under the key "device_config" and retrieved, is field-for-field
identical to the config independently built from the same input values. A
plain dict round-trip simulates st.session_state cross-page persistence
without needing a running Streamlit server.
"""

from etna import DeviceConfig
from app.components.device_sidebar import assemble_config

# Non-default values (mirrors tests/test_app_device_sidebar.py) so this test
# genuinely exercises assemble_config's field mapping rather than passing
# vacuously against a default-returning stub.
NON_DEFAULT_VALUES = {
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

# Independently constructed expected DeviceConfig — built directly from the
# same non-default values, NOT derived from assemble_config's output. This
# is the anchor that keeps the test from being a vacuous self-round-trip.
EXPECTED_CONFIG = DeviceConfig(**NON_DEFAULT_VALUES)


def test_config_persistence_key():
    cfg = assemble_config(NON_DEFAULT_VALUES)

    session_state = {}
    session_state["device_config"] = cfg

    retrieved = session_state["device_config"]

    for field_name in DeviceConfig.__dataclass_fields__:
        assert getattr(retrieved, field_name) == getattr(EXPECTED_CONFIG, field_name), (
            f"field {field_name!r} lost persistence-key round-trip fidelity: "
            f"expected {getattr(EXPECTED_CONFIG, field_name)!r}, "
            f"got {getattr(retrieved, field_name)!r}"
        )
