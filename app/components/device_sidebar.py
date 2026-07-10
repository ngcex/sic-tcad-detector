"""Device configuration sidebar: pure config-assembly seam + Streamlit renderer.

`assemble_config(values: dict) -> DeviceConfig` is a pure function with ZERO
Streamlit dependency — it is imported directly by unit tests
(tests/test_app_device_sidebar.py, tests/test_app_session.py) so the
config-assembly and mode-mapping logic (UI-02, UI-07) is testable without a
running Streamlit server.

`render_device_sidebar()` is the only function in this module that touches
`st.*` — it renders the reactive (non-form) sidebar widgets, collects their
values into a dict, calls `assemble_config()`, and stores the result under
the single non-widget key `st.session_state["device_config"]`.

STUB (Task 2 / RED): `assemble_config` currently returns the DeviceConfig
default, ignoring `values`, so the Wave-0 unit tests fail on assertions
(not on import). The real implementation lands in Task 3 (GREEN).
"""

from __future__ import annotations

from petringa import DeviceConfig


def assemble_config(values: dict) -> DeviceConfig:
    """Turn sidebar form values into a DeviceConfig.

    STUB: returns the DeviceConfig default, ignoring `values`. Real mapping
    logic (mode-consistency rules for doping_profile/N_D) lands in Task 3.
    """
    return DeviceConfig()
