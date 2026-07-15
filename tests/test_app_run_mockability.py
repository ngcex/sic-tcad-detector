"""Wave 0 spike (Phase 39): prove etna.run_cv is mockable via monkeypatch
under AppTest.from_function.

This is the empirical confirmation of RESEARCH A6: pages must reference
facades as a module attribute (`import etna; etna.run_cv(cfg)`)
rather than `from etna import run_cv`, so that `monkeypatch.setattr`
on the `etna` module is visible to the page code under test. Page
plans 39-03/39-04 depend on this seam to avoid real, expensive devsim
solves in their AppTest suites.
"""

from __future__ import annotations

import numpy as np
from streamlit.testing.v1 import AppTest

import etna
from etna import DeviceConfig


def test_run_cv_mockable_via_module_attribute(monkeypatch):
    """etna.run_cv, referenced as a module attribute, is interceptable
    by monkeypatch under AppTest.from_function."""

    def fake_run_cv(cfg, **kwargs):
        return etna.SimResult(
            config=cfg,
            sim_type="cv",
            x=np.array([0.0, -10.0]),
            y=np.array([1e-9, 2e-9]),
            metadata={
                "one_over_C_squared": np.array([1e18, 2.5e17]),
                "depletion_widths": np.array([1e-4, 2e-4]),
                "area_cm2": 1e-4,
            },
            mesh=None,
        )

    monkeypatch.setattr(etna, "run_cv", fake_run_cv)

    def _run_cv_wrapper():
        import streamlit as st

        import etna
        from etna import DeviceConfig

        result = etna.run_cv(DeviceConfig())
        st.session_state["cv_result"] = result

    at = AppTest.from_function(_run_cv_wrapper)
    at.run()

    assert at.exception == [], f"wrapper raised: {at.exception}"
    cv_result = at.session_state["cv_result"]
    assert cv_result.sim_type == "cv"
    assert len(cv_result.x) == 2
